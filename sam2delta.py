#!/usr/bin/python3
import argparse
from collections import defaultdict

"""
This utility converts a SAM file to a nucmer delta file.
"""


class SAMAlignment:

    def __init__(self, in_ref_header, in_query_header, in_ref_start, in_cigar, in_flag, in_seq_len):
        self.seq_len = in_seq_len
        self.flag = int(in_flag)
        self.is_rc = None

        # Check if the query was reverse complemented
        if self.flag & 16:
            self.is_rc = True
        else:
            self.is_rc = False

        self.cigar = in_cigar
        self.parsed_cigar = []
        self._parse_cigar()

        if self.is_rc:
            self.parsed_cigar = self.parsed_cigar[::-1]

        self.ref_header = in_ref_header
        self.query_header = in_query_header
        self.ref_start = int(in_ref_start)
        self.ref_end = None
        self.query_start = None
        self.query_end = None
        self.query_aln_len = 0
        self.ref_aln_len = 0

        num_hard_clipped = 0
        for i in self.parsed_cigar:
            if i[-1] == 'H':
                num_hard_clipped += int(i[:-1])

        if self.seq_len == 1:
            self.seq_len = 0
            for i in self.parsed_cigar:
                if i[-1] == 'M' or i[-1] == 'I' or i[-1] == 'S' or i[-1] == 'X' or i[-1] == '=':
                    self.seq_len += int(i[:-1])

        self.query_len = self.seq_len + num_hard_clipped

        self._get_query_start()
        self._get_query_aln_len()
        self._get_ref_aln_len()
        self._get_query_end()
        self._get_ref_end()

    def __str__(self):
        return ' '.join([self.ref_header, self.query_header])

    def __repr__(self):
        return '<' +  ' '.join([self.ref_header, self.query_header]) + '>'

    def _get_query_start(self):
        """
        An alignment will either start with soft clipping, hard clipping, or with one of the alignment
        codes (e.g. Match/Mismatch (M) or InDel (I/D)). If hard or soft clipped, the start of the alignment is
        the very next base after clipping. If the alignment has commenced, then the start of the alignment is 1.

        Notes:
        At least for now, I am only specifying this for minimap2 output, which I believe only has the following
        characters in cigar strings:
        {'I', 'S', 'M', 'D', 'H'}
        so I will only handle these cases.
        """
        if self.parsed_cigar[0][-1] == 'I' or self.parsed_cigar[0][-1] == 'D' or self.parsed_cigar[0][-1] == 'M':
            self.query_start = 1
        else:
            # The alignment has started with either soft or hard clipping
            # Need to double check if a 1 needs to be added here
            self.query_start = int(self.parsed_cigar[0][:-1]) + 1

    def _get_query_aln_len(self):
        """ Just the addition of all cigar codes which "consume" the query (M, I)"""
        for i in self.parsed_cigar:
            if i[-1] == 'M' or i[-1] == 'I':
                self.query_aln_len += int(i[:-1])

    def _get_ref_aln_len(self):
        for i in self.parsed_cigar:
            if i[-1] == 'M' or i[-1] == 'D':
                self.ref_aln_len += int(i[:-1])

    def _get_query_end(self):
        """ This is the length of the alignment + query start """
        self.query_end = self.query_start + self.query_aln_len

    def _get_ref_end(self):
        """ This is the length of the alignment + query start """
        self.ref_end = self.ref_start + self.ref_aln_len

    def _parse_cigar(self):
        cigar_chars = {
            'M',
            'I',
            'D',
            'N',
            'S',
            'H',
            'P',
            '=',
            'X'
        }

        this_field = ''
        for char in self.cigar:
            this_field += char
            if char in cigar_chars:
                self.parsed_cigar.append(this_field)
                this_field = ''


def write_delta(in_alns, in_file_name):
    with open(in_file_name, 'w') as f:
        f.write('Fake first line\n')
        f.write('NUCMER\n')
        for aln in in_alns.keys():
            query_len = in_alns[aln][0].query_len
            #print (query_len)
            f.write('>%s\t%s\t%r\t%r\n' % (aln[0], aln[1], ref_chr_lens[aln[0]], query_len))
            for i in in_alns[aln]:
                f.write('%r\t%r\t%r\t%r\t%r\t%r\t%r\n' % (
                    i.ref_start,
                    i.ref_end,
                    i.query_start,
                    i.query_end,
                    0,
                    0,
                    0
                ))
                # Continue with the cigar string
                offsets = []
                cigar = i.parsed_cigar
                if cigar[0][-1] == 'S' or cigar[0][-1] == 'H':
                    cigar = cigar[1:-1]
                else:
                    cigar = cigar[:-1]

                counter = 1
                for code in cigar:
                    if code[-1] == 'M':
                        counter += int(code[:-1])
                    elif code[-1] == 'I':
                        offsets.append(counter)
                        num_I = int(code[:-1])
                        for i in range(1, num_I):
                            offsets.append(1)
                        counter = 1
                    elif code[-1] == 'D':
                        offsets.append(-1*counter)
                        num_I = int(code[:-1])
                        for i in range(1, num_I):
                            offsets.append(-1)
                        counter = 1
                    else:
                        raise ValueError('Unexpected CIGAR code')
                offsets.append(0)
                offsets = [str(i) for i in offsets]
                f.write('\n'.join(offsets) + '\n')


parser = argparse.ArgumentParser(description='Convert a SAM file to a nucmer delta file.')
parser.add_argument("sam_file", metavar="<alns.sam>", type=str, help="SAM file to convert")

args = parser.parse_args()
sam_file = args.sam_file

# Make a dictionary storing all alignments between and query and a reference sequence.
# key = (reference_header, query_header)
# value = list of alignments between those sequences. only store info needed for the delta file

alns = defaultdict(list)

# Read through the sam file
ref_chr_lens = dict()
with open(sam_file) as f:
    for line in f:
        # Skip SAM headers
        if line.startswith('@SQ'):
            header_list = line.split('\t')
            chrom = header_list[1].replace('SN:', '')
            chr_len = int(header_list[2].replace('LN:', ''))
            ref_chr_lens[chrom] = chr_len
            continue

        if line.startswith('@'):
            continue

        fields = line.split('\t')
        ref_header = fields[2]
        query_header = fields[0]
        ref_start = fields[3]
        cigar = fields[5]
        flag = fields[1]
        seq_len = len(fields[9])
        if not cigar == '*':
            x = SAMAlignment(
                ref_header,
                query_header,
                ref_start,
                cigar,
                flag,
                seq_len
            )
            alns[(ref_header, query_header)].append(x)

write_delta(alns, sam_file + '.delta')