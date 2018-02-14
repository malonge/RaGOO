#!/usr/bin/perl -w

# Authors: Maria Nattestad and Mike Schatz
# Email: mnattest@cshl.edu

use strict;
my @chromosome_filter_choices =  ("all-chromosomes","primary-chromosomes");
my @longrange_filter_choices = ("include-longrange","exclude-longrange","longrange-only");
my @output_file_choices = ("bed","bedpe");

my $USAGE = "Usage:\nAssemblytics_between_alignments.pl coords.tab minimum_event_size maximum_event_size [@chromosome_filter_choices] [@longrange_filter_choices] [@output_file_choices] > fusions.svs.bedpe ";

my $coordsfile = shift @ARGV or die $USAGE;
my $minimum_event_size = int(shift @ARGV);
my $maximum_event_size = int(shift @ARGV);
my $chromosome_filter = shift @ARGV or die $USAGE;
my $longrange_filter = shift @ARGV or die $USAGE;
my $output_file = shift @ARGV or die $USAGE;



# How close do alignments have to be in order to call deletions and insertions? (as opposed to contractions and expansions)
my $narrow_threshold = 50;



# Number of basepairs of distance in either the reference or the query before we call an SV long-range
my $longrange = $maximum_event_size;



# What is the longest two alignments can map apart in the query before we throw the variant between them away?
my $max_query_dist = 100000;



my %chromosome_filter_choices_hash = map { $_, 1 } @chromosome_filter_choices;
my %longrange_filter_choices_hash = map { $_, 1 } @longrange_filter_choices;
my %output_file_choices_hash = map { $_, 1 } @output_file_choices;

if ( $chromosome_filter_choices_hash{ $chromosome_filter } &&  $longrange_filter_choices_hash{ $longrange_filter } && $output_file_choices_hash { $output_file }) {
  # All is well with the world
} else {
  die $USAGE;
}

if ($longrange_filter ne "exclude-longrange" && $output_file eq "bed"){
  die "Cannot output bed while allowing long-range variants\n$USAGE";
}

# open COORDS, "./bin/show-coords -rclHT $deltafile |"
#   or die "Can't process $deltafile ($!)\n";

open COORDS, "$coordsfile"  or die "Can't process $coordsfile ($!)\n";

##open COORDS, "show-coords -rclHT $deltafile |"
##  or die "Can't process $deltafile ($!)\n";



## Require the flanking alignments are at least this long to call an SV
## Note there is no minimum length for fusions, this is determined by how 
## the delta file was filtered

my $MIN_SV_ALIGN = 100;


#my $minimum_event_size = 50;
my $approximately_zero = $narrow_threshold;


my %alignments;
my $numalignments = 0;


while (<COORDS>)
{
  chomp;
  my @vals = split /\s+/, $_;

  my $rid = $vals[6];
  my $qid = $vals[7];

  my $a;
  $a->{"rstart"} = $vals[0];
  $a->{"rend"}   = $vals[1];
  $a->{"qstart"} = $vals[2];
  $a->{"qend"}   = $vals[3];
  $a->{"rlen"}   = $vals[4];
  $a->{"qlen"}   = $vals[5];
  $a->{"rid"}    = $vals[6];
  $a->{"qid"}    = $vals[7];
  $a->{"str"}    = $_;

  $a->{"qidx"}   = 0;

  $a->{"qrc"} = ($a->{"qend"} > $a->{"qstart"}) ? 0 : 1;

  push @{$alignments{$qid}->{$rid}}, $a; # a is a hash with all the info for one alignment

  $numalignments++;
}



my $candidatefusions = 0;
my $candidatesvs     = 0;

my $sv_id_counter = 0;

my %svstats;

foreach my $qid (sort keys %alignments) # query name is the key for the alignments hash
{
  my @refs = sort keys %{$alignments{$qid}}; # grab all alignments of that query
  my $numref = scalar @refs;

  ## scan for fusions
  # if ($numref > 1) # if query aligns to multiple chromosomes
  # {
  #   my $allrefs = join " ", @refs; # join the names together for output

  #   print "== $qid [$numref] $allrefs\n"; # output the names of the chromosomes 
  #   $candidatefusions++;

  #   my $rcnt = 0;
  #   foreach my $rid (@refs)
  #   {
  #     print "--\n" if ($rcnt > 0);
  #     $rcnt++;

  #     foreach my $a (@{$alignments{$qid}->{$rid}})
  #     {
  #       my $str = $a->{"str"};
  #       print "$str\n";
  #     }
  #   }

  #   print "\n";
  # }


  ## Resort the alignments by query sort position
  my @qaligns;
  foreach my $rid (@refs)
  {
    foreach my $a (@{$alignments{$qid}->{$rid}})
    {
      push @qaligns, $a;
    }
  }

  ## Now record the index of the sorted query indices
  @qaligns = sort { $a->{"qstart"} <=> $b->{"qstart"}} @qaligns;
  for (my $i=0; $i < scalar @qaligns; $i++)
  {
    $qaligns[$i]->{"qidx"} = $i;
  }

  ## scan for SVs
  my $numalign = scalar @qaligns;

  if ($numalign > 1) # if the query has more than 1 alignment
  {
    ## note skip first one
    for (my $j = 1; $j < $numalign; $j++)
    {
      my $ai = $qaligns[$j-1];
      my $aj = $qaligns[$j];

      my $istr = $ai->{"str"};
      my $jstr = $aj->{"str"};

      # if ($ai->{"rid"} ne $aj->{"rid"})
      # {
      #   ## skip the fusions for now #############################################################################################
      #   next;
      # }

      my $rid = $ai->{"rid"};

      if (($ai->{"rlen"} >= $MIN_SV_ALIGN) && 
          ($aj->{"rlen"} >= $MIN_SV_ALIGN))
      {
        ## r alignments are always forward, q alignments may be flipped

        my $rpos;
        my $qpos;
        my $rdist = 0;
        my $qdist = 0;
        my $svtype = 0;

        my $chromi = $ai->{"rid"};
        my $chromj = $aj->{"rid"};
        my $posi;
        my $posj;
        my $strandi;
        my $strandj;

        $sv_id_counter++;

        if (($ai->{"qrc"} == 0) && ($aj->{"qrc"} == 0))
        {
          ## ri: [1 - 1000] | j: [2000 - 3000] => 1000 
          ## qi: [1 - 1000] | j: [2000 - 3000] => 1000 

          $svtype = "FF";

          $qdist = $aj->{"qstart"} - $ai->{"qend"};
          $rdist = $aj->{"rstart"} - $ai->{"rend"};

          if ($rdist >= 0) { $rpos = sprintf("%s:%d-%d:+", $rid, $ai->{"rend"}, $aj->{"rstart"}); }
          else             { $rpos = sprintf("%s:%d-%d:-", $rid, $aj->{"rstart"}, $ai->{"rend"}); }

          if ($qdist >= 0) { $qpos = sprintf("%s:%d-%d:+", $qid, $ai->{"qend"}, $aj->{"qstart"}); }
          else             { $qpos = sprintf("%s:%d-%d:-", $qid, $aj->{"qstart"}, $ai->{"qend"}); }

          # When the alignments are forward-forward, the connection point is at the end of the first (i: rend) and at the beginning of the second (j: rstart)
          #    i  + -   j
          # ------> -------->
          $posi = $ai->{"rend"}; 
          $posj = $aj->{"rstart"};
          $strandi = "+";
          $strandj = "-";

        }
        elsif (($ai->{"qrc"} == 1) && ($aj->{"qrc"} == 1))
        {
          ## ri: [2000 - 3000] | j: [1 - 1000] => 1000 
          ## qi: [1000 - 1]    | j: [3000 - 2000] => 1000 

          $svtype = "RR";

          $rdist = $ai->{"rstart"} - $aj->{"rend"}; 
          $qdist = $aj->{"qend"} - $ai->{"qstart"};

          if ($rdist >= 0) { $rpos = sprintf("%s:%d-%d:+", $rid, $aj->{"rend"}, $ai->{"rstart"}); }
          else             { $rpos = sprintf("%s:%d-%d:-", $rid, $ai->{"rstart"}, $aj->{"rend"}); }

          if ($qdist >= 0) { $qpos = sprintf("%s:%d-%d:+", $qid, $ai->{"qstart"}, $aj->{"qend"}); }
          else             { $qpos = sprintf("%s:%d-%d:-", $qid, $aj->{"qend"}, $ai->{"qstart"}); }

          # When the alignments are reverse-reverse, the connection point is at the beginning of the first (i: rstart) and at the end of the second (j: rend)
          #     j  + -    i
          # <------- <--------
          $posi = $ai->{"rstart"};  # rstart means first reference coordinate, not with respect to the contig
          $posj = $aj->{"rend"}; # rend means last reference coordinate, not with respect to the contig
          $strandi = "-";
          $strandj = "+";
        }
        elsif (($ai->{"qrc"} == 0) && ($aj->{"qrc"} == 1))
        {
          ## ri: [1 - 1000] | j: [2000 - 3000] => 1000 
          ## qi: [1 - 1000] | j: [3000 - 2000] => 1000 

          $svtype = "FR";

          $qdist = $aj->{"qend"} - $ai->{"qend"};
          $rdist = $aj->{"rstart"} - $ai->{"rend"};

          if ($rdist >= 0) { $rpos = sprintf("%s:%d-%d:+", $rid, $ai->{"rend"}, $aj->{"rstart"}); }
          else             { $rpos = sprintf("%s:%d-%d:-", $rid, $aj->{"rstart"}, $ai->{"rend"}); }

          if ($qdist >= 0) { $qpos = sprintf("%s:%d-%d:+", $qid, $ai->{"qend"}, $aj->{"qend"}); }
          else             { $qpos = sprintf("%s:%d-%d:-", $qid, $aj->{"qend"}, $ai->{"qend"}); }

          # When the alignments are forward-reverse, the connection point is at the beginning of the first (i: rstart) and at the end of the second (j: rend)
          #    i   +     j   +
          # -------> <--------
          $posi = $ai->{"rend"}; 
          $posj = $aj->{"rend"};
          $strandi = "+";
          $strandj = "+";

        }
        elsif (($ai->{"qrc"} == 1) && ($aj->{"qrc"} == 0))
        {
          ## ri: [1 - 1000] | j: [2000 - 3000] => 1000 
          ## qi: [1000 - 1] | j: [2000 - 3000] => 1000 

          $svtype = "RF";

          $qdist = $ai->{"qend"} - $aj->{"qend"};
          $rdist = $aj->{"rstart"} - $ai->{"rend"};

          if ($rdist >= 0) { $rpos = sprintf("%s:%d-%d:+", $rid, $ai->{"rend"}, $aj->{"rstart"}); }
          else             { $rpos = sprintf("%s:%d-%d:-", $rid, $aj->{"rstart"}, $ai->{"rend"}); }

          if ($qdist >= 0) { $qpos = sprintf("%s:%d-%d:+", $qid, $aj->{"qend"}, $ai->{"qend"}); }
          else             { $qpos = sprintf("%s:%d-%d:-", $qid, $ai->{"qend"}, $aj->{"qend"}); }

          # When the alignments are reverse-forward:
          # -   i    -    j   
          # <------- -------->
          $posi = $ai->{"rstart"}; 
          $posj = $aj->{"rstart"};
          $strandi = "-";
          $strandj = "-";
        }
        else
        {
          my $irc = $ai->{"qrc"};
          my $jrc = $aj->{"qrc"};

          print "ERROR: Unknown SV: $irc $jrc\n";
          print "$istr\n";
          print "$jstr\n";
          die "ERROR: Unknown SV: $irc $jrc\n";
        }
        
        my $totaldist = $rdist + $qdist;
        my $typeguess = "";

        my $abs_event_size = abs($rdist-$qdist);

        if ($chromi ne $chromj) { # interchromosomal
          $typeguess = "Interchromosomal";
          $rdist = 0;
        } else { # same chromosome
          if ($strandi eq $strandj) {
            $typeguess = "Inversion";
            $abs_event_size = $rdist;
          }
          elsif ($qdist > $rdist) {
            # both are significantly negative: (means the size of an overlapping region got larger, so tandem element expansion)
            if ($rdist > -1*$approximately_zero && $rdist < $approximately_zero  && $qdist > -1*$approximately_zero) {
              $typeguess = "Insertion";
                # split into out of nowhere (rdist ~ 0) vs. rdist is > 0: insertion_in_unmapped_region
            }
            else {
              if ($rdist < 0 || $qdist < 0) {
                $typeguess = "Tandem_expansion";  
              } else {
                $typeguess = "Repeat_expansion";  
              }
            }
          }
          elsif ($qdist < $rdist) {
            # both are significantly negative: (means the size of an overlapping region got smaller, so tandem element contraction)
            if ($rdist > -1*$approximately_zero && $qdist > -1*$approximately_zero && $qdist < $approximately_zero) {
              $typeguess = "Deletion";
              # split into out of nowhere (rdist ~ 0) vs. rdist is > 0: deletion_in_unmapped_region
            }
            else {
              if ($rdist < 0 || $qdist < 0) {
                $typeguess = "Tandem_contraction";  
              } else {
                $typeguess = "Repeat_contraction";  
              }
            }
          }
          else {
            $typeguess = "None";
          }
          
          if ($abs_event_size > $longrange) {   #  || abs($rdist) > $longrange || abs($qdist) > $longrange 
            $typeguess = "Longrange";
            if (abs($qdist) > $max_query_dist) {
              $typeguess = "None";
            }
          }
          # my $ratio; 
          # if ($qdist != 0){
          # #   $ratio = abs(($rdist/$qdist)-1);
          # #   if ($ratio < 0.1) {
          # #     $typeguess = "Equilibrium";
          # #   }
            #   if ($rdist==$qdist || abs($qdist) > $longrange) {
            #     $typeguess = "None";
            #   }
          # }
        }

# my @chromosome_filter_choices = ("all-chromosomes","primary-chromosomes");
# my @longrange_filter_choices = ("include-longrange","exclude-longrange");

        my $chromi_length = length $chromi; # length of the chromosome names: a way to filter to primary chromosomes and cut out alts and patches from the assembly
        my $chromj_length = length $chromj;
        if ($typeguess ne "Inversion" && $typeguess ne "None" && $abs_event_size >= $minimum_event_size) { # always required
          if ($chromosome_filter eq "all-chromosomes" || ($chromi_length < 6 && $chromj_length < 6)) { # test for primary chromosomes unless "all-chromosomes" is chosen
            if ($longrange_filter ne "exclude-longrange" || ($typeguess ne "Interchromosomal" && $typeguess ne "Longrange")) {
              if ($longrange_filter ne "longrange-only" || ($typeguess eq "Interchromosomal" || $typeguess eq "Longrange")) {
                if ($output_file eq "bedpe") {
                  print "$chromi\t$posi\t@{[$posi + 1]}\t$chromj\t$posj\t@{[$posj + 1]}\tAssemblytics_b_$sv_id_counter\t$abs_event_size\t$strandi\t$strandj\t$typeguess\t$rdist\t$qdist\t$qpos\t$abs_event_size\t$svtype\tbetween_alignments\n";  
                }
                else {
                  use List::Util qw(min max);
                  my $ref_start = min(($posi, $posj));
                  my $ref_stop = max(($posi, $posj));
                  if ($ref_stop eq $ref_start) {
                    $ref_stop = $ref_start + 1;
                  }
                  # "chrom","start","stop","name","event.size","strand","event.type","ref.dist","query.dist","contig.name"
                  print "$chromi\t$ref_start\t$ref_stop\tAssemblytics_b_$sv_id_counter\t$abs_event_size\t+\t$typeguess\t$rdist\t$qdist\t$qpos\tbetween_alignments\n";  
                }
              }
            }
          }
          #if ($filter_type ~~ ("primary-allsizes","primary-shortrange") {
           # && $typeguess ne "Interchromosomal" && $typeguess ne "Inversion" && $chromi_length < 6 && $chromj_length < 6 && $abs_event_size >= $minimum_event_size) {
        }
        $candidatesvs++;
        #push @{$svstats{$svtype}}, $totaldist;
      }
    }
  }
}



# print "Processed $numalignments alignments found $candidatefusions fusions and $candidatesvs SVs\n";
# print STDERR "Processed $numalignments alignments found $candidatefusions fusions and $candidatesvs SVs\n";

# foreach my $svtype (keys %svstats)
# {
#   my @events = @{$svstats{$svtype}};
#   my $cnt = scalar @events;

#   my $sum = 0.0;
#   foreach my $e (@events)
#   {
#     $sum += $e;
#   }

#   my $mean = sprintf ("%0.02f", $sum/$cnt);

#   print "svtype[$svtype]: $cnt $mean\n";
#   print STDERR "svtype[$svtype]: $cnt $mean\n";
# }









