#!/usr/bin/python3

## usage: ./make_scaf_group.py chr_size.txt first_N_scaf_groups
##  first_N_scaf_groups means that the first N scafs will be 
##  turned into their own individual scaffold group
## ex. ~/scripts/make_scaf_group.py ~/refs/Omykiss_Arlee_genome/Arlee-genome.chr-lengths.txt 32 

import sys

# set global variables
# the first N scafs will get their own groups
first_N_scaf_groups = int(sys.argv[2]) 
# set max size for a given scaffold group
scaf_size_cutoff = 50000000  
# list of lines to output
out_lines_list = []
# header line for output file
header_line = "id\tchrom\tstart\tstop\tangsd_chrom\tmh_label"

# initialize counter variables 
scaf_counter = 0
size_sum = 0

# add header line to ouput list of lines
out_lines_list.append(header_line)

with open(sys.argv[1]) as scaf_size_file:
  for line in scaf_size_file:
    values = line.strip().split() 
    scaf_name = values[0]
    scaf_size = values[1]
    # if one of the first N scaffolds, just update the counter
    if scaf_counter <= first_N_scaf_groups:
      scaf_counter+=1
    # if we've passed the first N scaffolds, then start creating
    # multi-scaffold groups
    else:
      # if adding this scaffold will push the group size over
      # the size limit, start a new scaffold 
      if size_sum + int(scaf_size) > scaf_size_cutoff:
        # reset size sum and add +1 to the scaffold counter
        size_sum = 0
        scaf_counter+=1
      # otherwise add the scaffold length to the size sum
      else:
        size_sum = size_sum + int(scaf_size)

    # create output line and add to list 
    out_line = "scaff_grp_"+str(scaf_counter)+"\t"+scaf_name+"\t1\t"+scaf_size+"\t"+scaf_name.replace("_","-")+"\t"+str(scaf_counter)
    out_lines_list.append(out_line)

# write outfile
with open('scaffold_groups.tsv', mode='wt', encoding='utf-8') as outfile:
  outfile.write('\n'.join(out_lines_list))
