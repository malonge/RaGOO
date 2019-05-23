import bisect

import numpy as np


class ReadCoverage:

    def __init__(self, in_paf):
        self.paf = in_paf
        self.glob_mean = None
        self.glob_std = None
        self.coverage_map = dict()
        self.ctg_lens = dict()
        self._make_coverage_map()
        self._get_glob_mean()

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

    @staticmethod
    def _smooth_supported_breaks(sup_breaks, break_types):
        """
        If there are multiple low coverage breaks in close proximity,
        merge it into one break at the lowest coverage point.
        """
        i = 0
        j = 1
        while j < len(sup_breaks):
            if break_types[i] == 'l' and break_types[j] == 'l':
                if abs(sup_breaks[i][0] - sup_breaks[j][0]) < 100000:
                    # Merge these two break points
                    sup_breaks[i] = min([sup_breaks[i], sup_breaks[j]], key=lambda x: x[1])
                    sup_breaks.pop(j)
                else:
                    i += 1
                    j += 1
            else:
                i += 1
                j += 1

        return [z[0] for z in sup_breaks]

    def _trim_ends(self, dist=25000):
        """ Remove the ends of the contigs from the coverage map. """
        for i in self.coverage_map:
            # Start with the beginning of the contig
            start_idx = 0
            end_idx = 0
            for j in range(len(self.coverage_map[i])):
                if self.coverage_map[i][j][0] < dist:
                    start_idx = j

                if self.coverage_map[i][j][0] > self.ctg_lens[i] - dist:
                    end_idx = j-1
                    break
            self.coverage_map[i] = self.coverage_map[i][start_idx:end_idx]

        # Remove contigs which don't have coverage info.
        header_keys = list(self.coverage_map.keys())
        for i in header_keys:
            if not self.coverage_map[i]:
                self.coverage_map.pop(i)

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
                        self.ctg_lens[L1[5]] = int(L1[6])

        # Sort these coverage positions and get the coverage map
        for i in alns_pos:
            self.coverage_map[i] = self._tabulate_coverage(sorted(alns_pos[i], key=lambda x: x[1]))

        self._trim_ends()

    def _get_glob_mean(self):
        L1 = []
        for i in self.coverage_map:

            # In the case where we have multiple coverage values for one position, take the last one.
            last_pos = 0
            curr_val = 0
            for j in self.coverage_map[i]:
                if j[0] == last_pos:
                    curr_val = j[1]
                else:
                    last_pos = j[0]
                    L1.append(curr_val)
                    cur_val = j[1]

        L1 = np.asarray(L1, dtype=np.int32)
        self.glob_mean = np.median(L1)
        self.glob_std = np.sqrt(self.glob_mean)

    def _get_index_range(self, header, start_ind, distance):
        """
        Get the list of indices that are contained within a distance around the start index
        """
        all_inds = []
        low_counter = 1

        # Check if the start point is at the end of the contig
        if start_ind == len(self.coverage_map[header]):
            start_ind -= 1

        start_pos = self.coverage_map[header][start_ind][0]
        is_low = False

        # Get all coverage map indices representing regions 50kbp upstream of the start
        while not is_low:
            next_ind = start_ind-low_counter
            if next_ind < 0:
                is_low = True
            else:
                next_pos = self.coverage_map[header][next_ind][0]
                if start_pos - next_pos > distance:
                    is_low = True
                else:
                    all_inds.append(next_ind)
                    low_counter += 1

        # Repeat for 50kbp downstream
        high_counter = 1
        is_high = False
        while not is_high:
            next_ind = start_ind + high_counter
            if next_ind >= len(self.coverage_map[header]):
                is_high = True
            else:
                next_pos = self.coverage_map[header][next_ind][0]
                if next_pos - start_pos > distance:
                    is_high = True
                else:
                    all_inds.append(next_ind)
                    high_counter += 1

        return sorted(all_inds)

    def check_break_cov(self, header, in_breaks, min_cov=None, max_cov=None):
        """
        Given a list of potential break points, verify if those break points occur around low or high coverage
        areas. If so, replace the candidate break point with the low/high coverage break point.
        :param header: contig header for these breaks
        :param in_breaks: list of candidate break points
        :param min_cov: break at coverage levels below this value
        :param max_cov: break at coverage levels above this value
        :return: list of real break points, or empty list if not breaking is recommended
        """
        # Check that we have coverage info for this contig.
        if header not in self.coverage_map:
            return []

        if min_cov is None or max_cov is None:
            # Automatically calculate min and max coverage
            min_cov = max(self.glob_mean - (self.glob_std*3), 0)
            max_cov = self.glob_mean + (self.glob_std*3)

        supported_breaks = []
        break_types = []

        for i in in_breaks:
            # Get the coverage for the position closest to this potential break point
            ins_ind = bisect.bisect_left(self.coverage_map[header], (i, 0))
            # Get the coverage for positions within 50 kbp of the candidate coverage position
            # Exclude positions near the ends of the sequence.
            ind_range = self._get_index_range(header, ins_ind, 50000)

            if len(set(ind_range)) > 1:
                # Check for low coverage
                lowest_cov = min(self.coverage_map[header][ind_range[0]:ind_range[-1]], key=lambda x: x[1])
                if lowest_cov[1] < min_cov:
                    supported_breaks.append(lowest_cov)
                    break_types.append('l')
                    continue

                # Check for high coverage
                highest_cov = max(self.coverage_map[header][ind_range[0]:ind_range[-1]], key=lambda x: x[1])
                if highest_cov[1] > max_cov:
                    supported_breaks.append(highest_cov)
                    break_types.append('h')

        return self._smooth_supported_breaks(supported_breaks, break_types)


