#!/usr/bin/env python


class GFFLine:

    def __init__(self, in_fields):
        self.seqname = in_fields[0]
        self.source = in_fields[1]
        self.feature = in_fields[2]
        self.start = int(in_fields[3])
        self.end = int(in_fields[4])
        self.score = in_fields[5]
        self.strand = in_fields[6]
        self.frame = in_fields[7]
        self.attribute = in_fields[8]

    def __str__(self):
        all_atts = [
            self.seqname,
            self.source,
            self.feature,
            self.start,
            self.end,
            self.score,
            self.strand,
            self.frame,
            self.attribute
        ]
        all_atts = [str(i) for i in all_atts]
        return '\t'.join(all_atts)


class GFFReader:

    def __init__(self, in_gff_file):
        self.gff_file = in_gff_file

    def parse_gff(self):
        with open(self.gff_file) as f:
            for line in f:
                if not line.startswith('#'):
                    L1 = line.rstrip().split('\t')
                    yield GFFLine(L1)