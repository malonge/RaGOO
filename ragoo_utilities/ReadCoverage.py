
class ReadCoverage:

    def __init__(self, in_paf):
        self.paf = in_paf
        self.coverage_map = dict()
        self._make_coverage_map()

    @staticmethod
    def _tabulate_coverage(cov_list):
        current_coverage = 0
        coverage_list = []
        seen = set()

        for header, pos in cov_list:
            if header in seen:
                current_coverage -= 1
                coverage_list.append((pos, current_coverage))
                seen.remove(header)
            else:
                current_coverage += 1
                coverage_list.append((pos, current_coverage))
                seen.add(header)
        return coverage_list

    def _make_coverage_map(self):
        """
        Populate self.coverage_map. This is a dictionary that associates each contig header with a list of alignment
        positions and their coverage levels.
        """
        # Associate with each contig header, a list of (query header, start), (query header, end)
        alns_pos = dict()
        with open(self.paf, 'r') as f:
            for line in f:
                L1 = line.rstrip().split('\t')

                # Only consider an alignment if at least 75% of the read aligned.
                if abs((int(L1[3]) - int(L1[2])))/int(L1[1]) >= 0.75:
                    if L1[5] in alns_pos:
                        alns_pos[L1[5]].append((L1[0], int(L1[7])))
                        alns_pos[L1[5]].append((L1[0], int(L1[8])))
                    else:
                        alns_pos[L1[5]] = [(L1[0], int(L1[7])), (L1[0], int(L1[8]))]

        # Sort these coverage positions and get the coverage map
        for i in alns_pos:
            self.coverage_map[i] = self._tabulate_coverage(sorted(alns_pos[i], key=lambda x: x[1]))

