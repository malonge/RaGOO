class PAFLine:
    """ Object to represent a single alignment in a minimap PAF file. """

    def __init__(self, in_line):
        """
        start positions should be before end positions for both query and target
        """
        self.line = in_line.rstrip().split('\t')
        self.contig = self.line[0]
        self.query_len = int(self.line[1])
        self.query_start = int(self.line[2])
        self.query_end = int(self.line[3])
        self.strand = self.line[4]
        self.ref_header = self.line[5]
        self.ref_len = int(self.line[6])
        self.ref_start = int(self.line[7])
        self.ref_end = int(self.line[8])
        self.num_match = int(self.line[9])
        self.aln_len = int(self.line[10])
        self.mapq = int(self.line[11])

        assert self.query_start <= self.query_end
        assert self.ref_start <= self.ref_end

    def __str__(self):
        return '\t'.join(self.line)

    def __eq__(self, other):
        return self.line == other.line


class PAFReader:

    def __init__(self, paf_file):
        self.paf_file = paf_file

    def parse_paf(self):
        with open(self.paf_file) as f:
            for line in f:
                yield PAFLine(line)