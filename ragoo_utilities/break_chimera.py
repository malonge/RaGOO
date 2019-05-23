from collections import defaultdict
import copy

from intervaltree import IntervalTree

from ragoo_utilities.ContigAlignment import UniqueContigAlignment


def get_ref_parts(alns, l, p, r):
    """

    :param alns: A ContigAlignment object for the contig in question
    :param l: Minimum alignment length to consider
    :param p: Minimum percentage of a reference chromosome that must be covered to be considered significant
    :param r: Minimum number of reference chromosome nucleotides needed to be covered to be considered significant
    :return: A list of reference chromosomes to which the input contig significantly aligns to
    """
    all_intervals = defaultdict(list)

    # Iterate over all alignments
    for i in range(len(alns.ref_headers)):
        if alns.aln_lens[i] < l:
            continue

        all_intervals[alns.ref_headers[i]].append((alns.ref_starts[i], alns.ref_ends[i]))

    ranges = dict()
    for i in all_intervals.keys():
        ranges[i] = 0

    # For each chromosome, sort the range and get the union interval length.
    # Algorithm and pseudocode given by Michael Kirsche
    for i in all_intervals.keys():
        sorted_intervals = sorted(all_intervals[i], key=lambda tup: tup[0])
        max_end = -1
        for j in sorted_intervals:
            start_new_terr = max(j[0], max_end)
            ranges[i] += max(0, j[1] - start_new_terr)
            max_end = max(max_end, j[1])

    total_sum = sum(ranges.values())
    chimeric_refs = []
    for i in ranges.keys():
        if (ranges[i]/total_sum)*100 > float(p) and ranges[i] > r:
            chimeric_refs.append(i)

    return chimeric_refs


def cluster_contig_alns(contig, alns, chroms, l):
    ctg_alns = alns[contig]

    # Filter contig alignments so as to only include alignments to desired chromosomes
    ctg_alns.filter_ref_chroms(chroms)

    # Intialize the list of start and end positions w.r.t the query
    query_pos = []

    for i in range(len(ctg_alns.ref_headers)):
        query_pos.append((ctg_alns.query_starts[i], ctg_alns.query_ends[i], i))

    final_order = [i[2] for i in sorted(query_pos)]
    final_refs = []
    for i in final_order:
        final_refs.append(ctg_alns.ref_headers[i])

    return get_borders(final_refs, ctg_alns, final_order)


def get_borders(ordered_refs, alns, order):
    borders = [0]
    current_ref = ordered_refs[0]
    for i in range(1, len(ordered_refs)-1):
        if ordered_refs[i] != current_ref and ordered_refs[i+1] != current_ref:
            current_ref = ordered_refs[i]
            borders.append(alns.query_ends[order[i-1]])

    borders.append(alns.query_lens[0])
    intervals = [(borders[i], borders[i+1]) for i in range(len(borders)-1)]
    return intervals


def avoid_gff_intervals(borders, gff_ins):
    """
    :param borders:
    :param gff_ins: all features
    :return:
    """

    # Make an interval tree from the intervals of features associated with this contig
    t = IntervalTree()
    for line in gff_ins:
        # If the interval is one bp long, skip
        if line.start == line.end:
            continue
        t[line.start:line.end] = (line.start, line.end)

    new_borders = []
    for i in borders:
        # Returns an empty set if a border does not fall in any interval.
        # Otherwise, returns a list of Interval objects
        interval_query = t[i[0]]
        if interval_query:
            while interval_query:
                # Make new border larger than the max end position of all intervals
                i = (i[0] + 1, i[1])
                interval_query = t[i[0]]

        interval_query = t[i[1]]
        if interval_query:
            while interval_query:
                # Make new border larger than the max end position of all intervals
                i = (i[0], i[1] + 1)
                interval_query = t[i[1]]

        new_borders.append(i)

    new_borders = new_borders[:-1] + [(new_borders[-1][0], borders[-1][1])]
    return new_borders


def update_gff(features, borders, contig):
    # Pop the features to be updated
    contig_feats = features.pop(contig)

    # Initialize lists for new broken contig headers
    for i in range(len(borders)):
        features[contig + '_chimera_broken:' + str(borders[i][0]) + '-' + str(borders[i][1])] = []

    # Could just use binary search instead of this tree since intervals are non-overlapping.
    # but input so small both are trivial
    t = IntervalTree()
    for i in borders:
        t[i[0]:i[1]] = i

    for i in contig_feats:
        query = t[i.start]
        assert len(query) == 1
        break_start = list(query)[0].begin
        break_end = list(query)[0].end
        query_border = (break_start, break_end)
        break_number = borders.index(query_border)
        i.seqname = contig + '_chimera_broken:' + str(borders[break_number][0]) + '-' + str(borders[break_number][1])
        i.start = i.start - break_start
        i.end = i.end - break_start
        features[contig + '_chimera_broken:' + str(borders[break_number][0]) + '-' + str(borders[break_number][1])].append(i)

    return features


def break_contig(contigs_dict, header, break_points):
    seq = contigs_dict.pop(header)
    test_seq = ''
    for i in range(len(break_points)):
        contigs_dict[header + '_chimera_broken:' + str(break_points[i][0]) + '-' + str(break_points[i][1])] = seq[break_points[i][0]: break_points[i][1]]
        test_seq += seq[break_points[i][0]: break_points[i][1]]

    assert test_seq == seq
    return contigs_dict


def get_intra_contigs(alns, l, d, c):
    """
    Flag contigs as being intrachromosomal chimeras
    :param alns:
    :param l: Minimum alignment length to consider
    :param d: Distance between consecutive adjacent alignments with respect to the reference. If larger than this, flag
    :param c: Distance between consecutive adjacent alignments with respect to the query. If larger than this, flag
    :return: dict of contigs and break points.
    """

    # Get only the header to which this contig mostly aligns to and filter out smaller alignments.
    uniq_aln = UniqueContigAlignment(alns)
    best_header = uniq_aln.ref_chrom
    ctg_alns = copy.deepcopy(alns)
    ctg_alns.filter_ref_chroms([best_header])
    ctg_alns.filter_lengths(l)

    # If there are no longer any alignments after length filtering, give up
    if not len(ctg_alns.ref_headers):
        return

    # Sort the alignments with respect to the reference start and end positions.
    ctg_alns.sort_by_ref()

    # Make a list of distance between alignments
    # first with respect to (wrt) the reference.
    distances_wrt_ref = []
    for i in range(len(ctg_alns.ref_headers)-1):
        distances_wrt_ref.append(ctg_alns.ref_starts[i+1] - ctg_alns.ref_starts[i])

    # next, with respect to (wrt) the contig.
    distances_wrt_ctg = []
    for i in range(len(ctg_alns.ref_headers) - 1):
        distances_wrt_ctg.append(abs(ctg_alns.query_starts[i + 1] - ctg_alns.query_starts[i]))

    # Next, assign the following two identities.
    #  1. When ordered by the reference, the alignments start at the beginning or the end of the query
    #  2. For the alignment which will be broken on, is it on the forward or reverse strand.

    is_query_start = True

    if ctg_alns.query_starts[0] >= ctg_alns.query_lens[0]/5:
        is_query_start = False

    # This conditional essentially checks if there are any break points for this contig.
    # Returns None otherwise (no return statement)
    if distances_wrt_ref:
        if max(distances_wrt_ref) > d:
            gap_index = distances_wrt_ref.index(max(distances_wrt_ref))
            a_alns_strands = ctg_alns.strands[:gap_index]
            if is_query_start:
                if a_alns_strands.count('-') > a_alns_strands.count('+'):
                    # The first subcontig is on the reverse strand
                    return (ctg_alns.contig, [(0, ctg_alns.query_ends[0]), (ctg_alns.query_ends[0], ctg_alns.query_lens[0])])
                else:
                    # The first subcontig is on the forward strand.
                    return (ctg_alns.contig, [(0, ctg_alns.query_ends[gap_index]), (ctg_alns.query_ends[gap_index], ctg_alns.query_lens[0])])
            else:
                # The first subcontig starts at the end of the contig
                if a_alns_strands.count('-') > a_alns_strands.count('+'):
                    # The first subcontig is on the reverse strand
                    return (ctg_alns.contig, [(0, ctg_alns.query_starts[gap_index]), (ctg_alns.query_starts[gap_index], ctg_alns.query_lens[0])])
                else:
                    # The first subcontig is on the forward strand.
                    return (ctg_alns.contig, [(0, ctg_alns.query_starts[0]), (ctg_alns.query_starts[0], ctg_alns.query_lens[0])])

        if max(distances_wrt_ctg) > c:
            gap_index = distances_wrt_ctg.index(max(distances_wrt_ctg)) + 1
            a_alns_strands = ctg_alns.strands[:gap_index]
            if is_query_start:
                if a_alns_strands.count('-') > a_alns_strands.count('+'):
                    # The first subcontig is on the reverse strand
                    return (ctg_alns.contig, [(0, ctg_alns.query_ends[0]), (ctg_alns.query_ends[0], ctg_alns.query_lens[0])])
                else:
                    # The first subcontig is on the forward strand.
                    return (ctg_alns.contig, [(0, ctg_alns.query_ends[gap_index]), (ctg_alns.query_ends[gap_index], ctg_alns.query_lens[0])])
            else:
                # The first subcontig starts at the end of the contig
                if a_alns_strands.count('-') > a_alns_strands.count('+'):
                    # The first subcontig is on the reverse strand
                    return (ctg_alns.contig, [(0, ctg_alns.query_starts[gap_index-1]), (ctg_alns.query_starts[gap_index-1], ctg_alns.query_lens[0])])
                else:
                    # The first subcontig is on the forward strand.
                    return (ctg_alns.contig, [(0, ctg_alns.query_starts[0]), (ctg_alns.query_starts[0], ctg_alns.query_lens[0])])

