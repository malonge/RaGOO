import time
import subprocess

from utilities.SeqReader import SeqReader

""" A collection of various helper functions"""

complements = str.maketrans("ACGTNURYSWKMBVDH", "TGCANAYRSWMKVBHD")


def reverse_complement(seq):
    """
    Reverse complement a nucleotide sequence.
    :param seq: Sequence to be reverse complemented
    :return: A reverse complemented sequence
    """
    return seq.translate(complements)[::-1]


def run(cmnd):
    """ Run command and report status. """
    log(' ---- Running : %s' % cmnd)
    if subprocess.call(cmnd, shell=True) != 0:
        raise RuntimeError('Failed : %s ' % cmnd)


def log(message):
    """ Log messages to standard output. """
    print (time.ctime() + ' ' + message)


def read_contigs(in_file):
    d = dict()
    x = SeqReader(in_file)
    for header, seq in x.parse_fasta():
        d[header.replace('>', '')] = seq
    return d