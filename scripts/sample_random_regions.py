#!/usr/bin/env python3

import os as os
import sys as sys
import logging as log
import io as io
import traceback as trb
import argparse as argp
import collections as col
import multiprocessing as mp

import numpy.random as rng
import intervaltree as ivt


FASTQ_QUAL_ENCODING = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""

Read = col.namedtuple('Read', 'line seqid sequence separator qualities record')


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser(prog="sample_random_regions.py", description=__doc__)
    parser.add_argument(
        "--debug",
        "-d",
        action="store_true",
        default=False,
        dest="debug",
        help="Print status and error messages to STDOUT. Otherwise, only"
        " errors/warnings will be reported to STDERR.",
    )
    parser.add_argument(
        "--input",
        "-i",
        type=str,
        required=True,
        dest="input",
        help="Full path to BED-like file with input regions (chrom, start, and end)",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        dest="output",
        help="Full path to output file (BED format). Directories will"
             " be created if they do not exist.",
    )
    parser.add_argument(
        "--chrom-sizes",
        '-cs',
        type=str,
        required=True,
        dest='chromsizes',
        help="Full path to 2-column file with chromosome sizes"
    )
    parser.add_argument(
        "--num-sets",
        "-s",
        type=int,
        default=1,
        dest='numsets',
    )
    parser.add_argument(
        "--num-cpu",
        "-n",
        type=int,
        default=1,
        dest="numcpu",
        help="Specify number of CPU cores to use for parallel computation. Note"
        " that there is no sanity checking for actually free/available cores."
        " Default: 1",
    )
    return parser.parse_args()


def read_chromosome_sizes(table):
    """
    :param table:
    :return:
    """
    chrom_sizes = dict()
    with open(table, 'r') as tsv:
        for line in tsv:
            if line.startswith('#'):
                continue
            else:
                chrom, size = line.strip().split()[:2]
                chrom_sizes[chrom] = int(size)
    return chrom_sizes


def read_region_lengths(bedfile):
    """
    :param bedfile:
    :return:
    """
    region_lengths = col.defaultdict(list)
    with open(bedfile, 'r') as tsv:
        for line in tsv:
            if line.startswith('#'):
                continue
            chrom, start, end = line.strip().split()[:3]
            length = int(end) - int(start)
            assert length > 0, 'Invalid region in file {}: {}'.format(bedfile, line.strip())
            region_lengths[chrom].append(length)
    return region_lengths


def sample_chromosome_regions(params):
    """
    :param params:
    :return:
    """
    chrom, size, lengths = params

    sample_lengths = col.deque(sorted(lengths, reverse=True))
    sample_lengths.append(None)
    coords = []
    ivtree = ivt.IntervalTree()
    matched = 0
    cycles = 0

    while matched < len(lengths):
        query_length = sample_lengths.popleft()
        if query_length is None:
            # encountered sentinel
            cycles += 1
            sample_lengths.append(None)
            if cycles > 5 and matched >= round(len(lengths) * 0.95, 0):
                break
            elif cycles > 10:
                break
            else:
                continue
        upper_bound = size - query_length - 1
        assert upper_bound > 0, 'Malformed something: {} - {} - {}'.format(size, upper_bound, query_length)
        rand_start = rng.randint(0, upper_bound)
        rand_end = rand_start + query_length
        if ivtree.overlaps(rand_start, rand_end):
            sample_lengths.append(query_length)
        else:
            ivtree.addi(rand_start, rand_end)
            matched += 1
            coords.append((rand_start, rand_end))

    return chrom, sorted(coords)


def main(logger, cargs):
    """
    :param logger:
    :param cargs:
    :return:
    """
    logger.debug("Starting computations")
    chrom_sizes = read_chromosome_sizes(cargs.chromsizes)
    logger.debug('Read chromosome sizes')
    region_lengths = read_region_lengths(cargs.input)
    logger.debug('Read region lengths')

    params = []
    for chrom, lengths in region_lengths.items():
        params.append((chrom, chrom_sizes[chrom], tuple(lengths)))

    logger.debug('Created parameter list of size {} for processing'.format(len(params)))

    for n in range(cargs.numsets):
        output_file = cargs.output
        if cargs.numsets > 1 and n > 0:
            output_file = output_file.rsplit('.', 1)[0] + '.{}.bed'.format(n)

        os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
        with open(output_file, 'w'):
            pass

        with mp.Pool(min(len(params), cargs.numcpu)) as pool:
            resit = pool.imap_unordered(sample_chromosome_regions, params)
            for chrom, coords in resit:
                logger.debug('Received results for chromosome {}'.format(chrom))

                with open(output_file, 'a') as dump:
                    for s, e in coords:
                        dump.write('{}\t{}\t{}\n'.format(chrom, s, e))

    return


if __name__ == "__main__":
    logger = None
    rc = 0
    try:
        rng.seed(0)  # fix seed for replicability
        log_msg_format = "%(asctime)s | %(levelname)s | %(message)s"
        cargs = parse_command_line()
        if cargs.debug:
            log.basicConfig(stream=sys.stdout, level=log.DEBUG, format=log_msg_format)
        else:
            log.basicConfig(stream=sys.stderr, level=log.WARNING, format=log_msg_format)
        logger = log.getLogger()
        logger.debug("Logger initiated")
        main(logger, cargs)
        logger.debug("Run completed - exit")
        log.shutdown()
    except Exception as exc:
        rc = 1
        if logger is not None:
            logger.error("Unrecoverable error: {}".format(str(exc)))
            logger.debug("=== TRACEBACK ===\n\n")
            buf = io.StringIO()
            trb.print_exc(file=buf)
            logger.error(buf.getvalue())
            logger.debug("Exit\n")
            log.shutdown()
        else:
            trb.print_exc()
    finally:
        sys.exit(rc)
