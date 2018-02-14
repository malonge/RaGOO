#! /usr/bin/env python


# Author: Maria Nattestad
# Email: mnattest@cshl.edu
# This script is part of Assemblytics, a program to detect and analyze structural variants from an assembly aligned to a reference genome using MUMmer. 


import argparse
import gzip
# from intervaltree import *
import time

import numpy as np
import operator



def run(args):
    filename = args.delta
    unique_length = args.unique_length
    output_filename = args.out
    keep_small_uniques = args.keep_small_uniques
    
    f = open(filename)
    header1 = f.readline()
    if header1[0:2]=="\x1f\x8b":
        f.close()
        f = gzip.open(filename)


    linecounter = 0

    current_query_name = ""
    current_header = ""

    lines_by_query = {}
    header_lines_by_query = {}

    before = time.time()
    last = before

    existing_query_names = set()

    for line in f:
        if line[0]==">":

            fields = line.strip().split()
            current_query_name = fields[1]
            current_header = line.strip()
            if current_query_name not in existing_query_names:
                lines_by_query[current_query_name] = []
                header_lines_by_query[current_query_name] = []
                existing_query_names.add(current_query_name)
        else:
            fields = line.strip().split()
            if len(fields) > 4:
                # sometimes start and end are the other way around, but for this they need to be in order
                query_min = min([int(fields[2]),int(fields[3])])
                query_max = max([int(fields[2]),int(fields[3])])

                ##########  TESTING ONLY  ###########
                # lines_by_query[current_query_name] = (query_min,query_max)
                # test_list = test_list + [(query_min,query_max)]
                #####################################

                lines_by_query[current_query_name].append((query_min,query_max))
                header_lines_by_query[current_query_name].append(current_header)
        # linecounter += 1
        # if linecounter % 10000000 == 0:
        #     print "%d,%f" % (linecounter, time.time()-last)
        #     last = time.time()
        

    f.close()
    

    before = time.time()
    alignments_to_keep = {}
    num_queries = len(lines_by_query)
    
    num_query_step_to_report = num_queries/100
    if num_queries < 100:
        num_query_step_to_report = num_queries/10
    if num_queries < 10:
        num_query_step_to_report = 1

    query_counter = 0

    for query in lines_by_query:

        ################   TESTING    ####################   

        # results_intervaltree = summarize_intervaltree(lines_by_query[query], unique_length_required = unique_length)
        # intervaltree_filtered_out = set(range(0,len(lines_by_query[query]))) - set(results_intervaltree)
    
        # results_planesweep = summarize_planesweep(lines_by_query[query], unique_length_required = unique_length) 
        # planesweep_filtered_out = set(range(0,len(lines_by_query[query]))) - set(results_planesweep)
        # if intervaltree_filtered_out == planesweep_filtered_out :
        #     num_matches += 1
        # else:
        #     num_mismatches += 1
        #     print "MISMATCH:"
        #     print "number of alignments:", len(lines_by_query[query])
        #     print "results_intervaltree:"
        #     print results_intervaltree
        #     for i in results_intervaltree:
        #         print lines_by_query[query][i]
        #     print "results_planesweep:"
        #     print results_planesweep
        #     for i in results_planesweep:
        #         print lines_by_query[query][i]
        ################   TESTING    ####################

        alignments_to_keep[query] = summarize_planesweep(lines_by_query[query], unique_length_required = unique_length,keep_small_uniques=keep_small_uniques)

        query_counter += 1

    before = time.time()


    fout = open(output_filename + ".Assemblytics.unique_length_filtered_l%d.delta" % (unique_length),'w')
    

    f = open(filename)
    header1 = f.readline()
    if header1[0:2]=="\x1f\x8b":
        f.close()
        f = gzip.open(filename)
        header1 = f.readline()

    fout.write(header1) # write the first line that we read already
    fout.write(f.readline())
    
    linecounter = 0

    # For filtered delta file:
    list_of_alignments_to_keep = []
    alignment_counter = {}
    keep_printing = False

    # For coords:
    current_query_name = ""
    current_query_position = 0
    fcoords_out_tab = open(output_filename + ".coords.tab",'w')
    fcoords_out_csv = open(output_filename + ".coords.csv",'w')
    fcoords_out_csv.write("ref_start,ref_end,query_start,query_end,ref_length,query_length,ref,query,tag\n")


    # For basic assembly stats:
    ref_sequences = set()
    query_sequences = set()
    ref_lengths = []
    query_lengths = []

    f_stats_out = open(output_filename + ".Assemblytics_assembly_stats.txt","w")

    for line in f:
        linecounter += 1
        if line[0]==">":
            fields = line.strip().split()
            
            # For delta file output:
            query = fields[1]
            list_of_alignments_to_keep = alignments_to_keep[query]

            header_needed = False
            for index in list_of_alignments_to_keep:
                if line.strip() == header_lines_by_query[query][index]:
                    header_needed = True
            if header_needed == True:
                fout.write(line) # if we have any alignments under this header, print the header
            alignment_counter[query] = alignment_counter.get(query,0)

            # For coords:
            current_reference_name = fields[0][1:]
            current_query_name = fields[1]

            current_reference_size = int(fields[2])
            current_query_size = int(fields[3])

            # For basic assembly stats:
            if not current_reference_name in ref_sequences:
                ref_lengths.append(current_reference_size)
                ref_sequences.add(current_reference_name)
            if not current_query_name in query_sequences:
                query_lengths.append(current_query_size)
                query_sequences.add(current_query_name)

        else:
            fields = line.strip().split()
            if len(fields) > 4:
                # For coords:
                ref_start = int(fields[0])
                ref_end = int(fields[1])
                query_start = int(fields[2])
                query_end = int(fields[3])
                csv_tag = "repetitive"
                if alignment_counter[query] in list_of_alignments_to_keep:
                    fout.write(line)
                    fcoords_out_tab.write("\t".join(map(str,[ref_start,ref_end,query_start, query_end,current_reference_size,current_query_size,current_reference_name,current_query_name])) + "\n")
                    csv_tag = "unique"
                    keep_printing = True
                else:
                    keep_printing = False
                fcoords_out_csv.write(",".join(map(str,[ref_start,ref_end,query_start, query_end,current_reference_size,current_query_size,current_reference_name.replace(",","_"),current_query_name.replace(",","_"),csv_tag])) + "\n")
                alignment_counter[query] = alignment_counter[query] + 1

            elif keep_printing == True:
                fout.write(line)

    fcoords_out_tab.close()
    fcoords_out_csv.close()


    ref_lengths.sort()
    query_lengths.sort()

    # Assembly statistics
    ref_lengths = np.array(ref_lengths)
    query_lengths = np.array(query_lengths)

    f_stats_out.write("Reference: %s\n" % (header1.split()[0].split("/")[-1]))
    f_stats_out.write( "Number of sequences: %s\n" % intWithCommas(len(ref_lengths)))
    f_stats_out.write( "Total sequence length: %s\n" %  gig_meg(sum(ref_lengths)))
    f_stats_out.write( "Mean: %s\n" % gig_meg(np.mean(ref_lengths)))
    f_stats_out.write( "Min: %s\n" % gig_meg(np.min(ref_lengths)))
    f_stats_out.write( "Max: %s\n" % gig_meg(np.max(ref_lengths)))
    f_stats_out.write( "N50: %s\n" % gig_meg(N50(ref_lengths)))
    f_stats_out.write( "\n\n")
    f_stats_out.write( "Query: %s\n" % header1.split()[1].split("/")[-1])
    f_stats_out.write( "Number of sequences: %s\n" % intWithCommas(len(query_lengths)))
    f_stats_out.write( "Total sequence length: %s\n" % gig_meg(sum(query_lengths)))
    f_stats_out.write( "Mean: %s\n" % gig_meg(np.mean(query_lengths)))
    f_stats_out.write( "Min: %s\n" % gig_meg(np.min(query_lengths)))
    f_stats_out.write( "Max: %s\n" % gig_meg(np.max(query_lengths)))
    f_stats_out.write( "N50: %s\n" % gig_meg(N50(query_lengths)))


    f.close()
    fout.close()
    f_stats_out.close()

def N50(sorted_list):
    # List should be sorted as increasing

    # We flip the list around here so we start with the largest element
    cumsum = 0
    for length in sorted_list[::-1]:
        cumsum += length
        if cumsum >= sum(sorted_list)/2:
            return length


def gig_meg(number,digits = 2):
    gig = 1000000000.
    meg = 1000000.
    kil = 1000.

    if number > gig:
        return str(round(number/gig,digits)) + " Gbp"
    elif number > meg:
        return str(round(number/meg,digits)) + " Mbp"
    elif number > kil:
        return str(round(number/kil,digits)) + " Kbp"
    else:
        return str(number) + " bp"

def intWithCommas(x):
    if type(x) not in [type(0)]:
        raise TypeError("Parameter must be an integer.")
    if x < 0:
        return '-' + intWithCommas(-x)
    result = ''
    while x >= 1000:
        x, r = divmod(x, 1000)
        result = ",%03d%s" % (r, result)
    return "%d%s" % (x, result)


def summarize_planesweep(lines,unique_length_required, keep_small_uniques=False):

    alignments_to_keep = []
    # print len(lines)

    # If no alignments:
    if len(lines)==0:
        return []

    # If only one alignment:
    if len(lines) == 1:
        if keep_small_uniques == True or abs(lines[0][1] - lines[0][0]) >= unique_length_required:
            return [0]
        else:
            return []

    starts_and_stops = []
    for query_min,query_max in lines:
        # print query_min, query_max
        starts_and_stops.append((query_min,"start"))
        starts_and_stops.append((query_max,"stop"))


    sorted_starts_and_stops = sorted(starts_and_stops,key=operator.itemgetter(0))
    # print sorted_starts_and_stops

    current_coverage = 0
    last_position = -1
    # sorted_unique_intervals = []
    sorted_unique_intervals_left = []
    sorted_unique_intervals_right = []
    for pos,change in sorted_starts_and_stops:
        # print sorted_starts_and_stops[i]
        # pos = sorted_starts_and_stops[i][0]
        # change = sorted_starts_and_stops[i][1]
        
        # print pos,change
        # First alignment only:
        # if last_position == -1:
        #     last_position = pos
        #     continue

        # print last_position,pos,current_coverage

        if current_coverage == 1:
            # sorted_unique_intervals.append((last_position,pos))
            sorted_unique_intervals_left.append(last_position)
            sorted_unique_intervals_right.append(pos)

        if change == "start":
            current_coverage += 1
        else:
            current_coverage -= 1
        last_position = pos


    linecounter = 0
    for query_min,query_max in lines:

        i = binary_search(query_min,sorted_unique_intervals_left,0,len(sorted_unique_intervals_left))

        exact_match = False
        if sorted_unique_intervals_left[i] == query_min and sorted_unique_intervals_right[i] == query_max:
            exact_match = True
        sum_uniq = 0
        while i < len(sorted_unique_intervals_left) and sorted_unique_intervals_left[i] >= query_min and sorted_unique_intervals_right[i] <= query_max:
            sum_uniq += sorted_unique_intervals_right[i] - sorted_unique_intervals_left[i]
            i += 1

        # print query_min,query_max,sum_uniq
        if sum_uniq >= unique_length_required:
            alignments_to_keep.append(linecounter)
        elif keep_small_uniques == True and exact_match == True:
            alignments_to_keep.append(linecounter)
            # print "Keeping small alignment:", query_min, query_max
            # print sorted_unique_intervals_left[i-1],sorted_unique_intervals_right[i-1]

        linecounter += 1

    return alignments_to_keep



def binary_search(query, numbers, left, right):
    #  Returns index of the matching element or the first element to the right
    
    if left >= right:
        return right
    mid = (right+left)//2
    

    if query == numbers[mid]:
        return mid
    elif query < numbers[mid]:
        return binary_search(query,numbers,left,mid)
    else: # if query > numbers[mid]:
        return binary_search(query,numbers,mid+1,right)


# def summarize_intervaltree(lines, unique_length_required):

#     alignments_to_keep = []
#     # print len(lines)

#     if len(lines)==0:
#         return alignments_to_keep

#     if len(lines) == 1:
#         if abs(lines[0][1] - lines[0][0]) >= unique_length_required:
#             return [0]


#     starts_and_stops = []
#     for query_min,query_max in lines:
#         starts_and_stops.append((query_min,query_max))

#     # build full tree
#     tree = IntervalTree.from_tuples(starts_and_stops) 
    

#     # for each interval (keeping the same order as the lines in the input file)
#     line_counter = 0
#     for query_min,query_max in lines:
        
#         # create a tree object from the current interval
#         this_interval = IntervalTree.from_tuples([(query_min,query_max)])

#         # create a copy of the tree without this one interval
#         rest_of_tree = tree - this_interval

#         # find difference between this interval and the rest of the tree by subtracting out the other intervals one by one
#         for other_interval in rest_of_tree:
#             this_interval.chop(other_interval.begin, other_interval.end)
        
#         # loop through to count the total number of unique basepairs
#         total_unique_length = 0
#         for sub_interval in this_interval:
#             total_unique_length += sub_interval.end - sub_interval.begin

#         # if the total unique length is above our threshold, add the index to the list we are reporting       
#         if total_unique_length >= unique_length_required:
#             alignments_to_keep.append(line_counter)
#         line_counter += 1


#     return alignments_to_keep


def main():
    parser=argparse.ArgumentParser(description="Filters alignments in delta file based whether each alignment has a unique sequence anchoring it")
    parser.add_argument("--delta",help="delta file" ,dest="delta", type=str, required=True)
    parser.add_argument("--out",help="output file" ,dest="out", type=str, required=True)
    parser.add_argument("--unique-length",help="The total length of unique sequence an alignment must have on the query side to be retained. Default: 10000" ,dest="unique_length",type=int, default=10000)
    parser.add_argument("--keep-small-uniques",help="Keep small aligments (below the unique anchor length) if they are completely unique without any part of the alignment mapping multiple places" ,dest="keep_small_uniques",action="store_true")
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main()
