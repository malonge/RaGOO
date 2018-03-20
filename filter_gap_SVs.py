import os
import re

from ragoo_utilities.SeqReader import SeqReader
from ragoo_utilities.utilities import run


class BaseSequence(object):

    def __init__(self, in_sequence):
        if not isinstance(in_sequence, str):
            raise AttributeError('Only a string can be used to instantiate this class.')
        self.sequence = in_sequence.upper()


class GapSequence(BaseSequence):
    def get_gap_coords(self):
        """ Find all of the gap string indices for this sequence. """
        return re.finditer(r'N+', self.sequence)


def make_gaps_bed():
    x = SeqReader('../ragoo.fasta')
    with open('ragoo.gaps.bed', 'w') as f:
        for header, sequence in x.parse_fasta():
            gap_sequence = GapSequence(sequence)
            all_coordinates = [(m.start(0), m.end(0)) for m in gap_sequence.get_gap_coords()]
            for i in all_coordinates:
                f.write(header[1:] + '\t' + str(i[0]) + '\t' + str(i[1]) + '\n')


def make_svs_bed():
    f_in = open('assemblytics_out.Assemblytics_structural_variants.bed', 'r')
    f_in.readline()
    with open('svs_wrt_query.bed', 'w') as f_out:
        for line in f_in:
            L1 = line.split('\t')
            query_coords = L1[9]
            L2 = query_coords.split(':')
            header = L2[0]
            L3 = L2[1].split('-')
            start, stop = int(L3[0]), int(L3[1])
            # Disregard strandedness
            if stop < start:
                start, stop = stop, start

            f_out.write('%s\t%r\t%r\t' %(header, start, stop))
            f_out.write(line)


def filter():
    # Bedtools intersect to filter out SVs mostly caused by gaps
    cmd = 'bedtools intersect -v -f 0.25 -a svs_wrt_query.bed -b ragoo.gaps.bed -nonamecheck > svs_wrt_query.filtered.bed'
    if not os.path.isfile('svs_wrt_query.filtered.bed'):
        run(cmd)

    if not os.path.isfile('ragoo.SVs.bed'):
        run('cut -f 4,5,6,7,8,9,10,11,12,13,14 svs_wrt_query.filtered.bed > ragoo.SVs.bed')




def main():
    # Make a new bed file w.r.t the query
    make_gaps_bed()
    make_svs_bed()
    filter()


if __name__ == "__main__":
    main()
