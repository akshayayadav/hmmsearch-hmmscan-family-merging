#!/usr/bin/env python

from __future__ import division
import re
import sys
import os

#fam1 - family id that is processed - large fam
#fam2 - family id from which the sequences were predicted as missing by fam1 - small fam

##########################################################################################################################################

def fam_fasta_dicts(fam_fasta_dirName):
	seqid_famid_dict={}
	famid_famsize_dict={}
	for fam_fasta_fileName in os.listdir(fam_fasta_dirName):
		fam_fasta_file = open(fam_fasta_dirName+"/"+fam_fasta_fileName,"r")
		famid_famsize_dict[fam_fasta_fileName]=0
		for line in fam_fasta_file:
			line = line.rstrip()
			if not (re.match(r'^>',line)):
				continue
			seqid_famid_dict[line[1:]]=fam_fasta_fileName
			famid_famsize_dict[fam_fasta_fileName]+=1
		fam_fasta_file.close()

	return([seqid_famid_dict, famid_famsize_dict])		

def get_large_fam_small_fam_seqcount_dict(predicted_missing_seqs_fileName, seqid_famid_dict):
	large_fam_small_fam_seqcount_dict={}
	predicted_missing_seqs_file = open(predicted_missing_seqs_fileName,"r")
	for line in predicted_missing_seqs_file:
		line = line.rstrip()
		linearr = re.split(r'\s+',line)
		large_famid = linearr[0]
		seqid = linearr[1]
		if not(seqid_famid_dict.has_key(seqid)):
			continue
		small_famid = seqid_famid_dict[seqid]
		if(large_fam_small_fam_seqcount_dict.has_key(large_famid)):
			if(large_fam_small_fam_seqcount_dict[large_famid].has_key(small_famid)):
				large_fam_small_fam_seqcount_dict[large_famid][small_famid]+=1
			else:
				large_fam_small_fam_seqcount_dict[large_famid][small_famid]=1
		else:
			large_fam_small_fam_seqcount_dict[large_famid]={}
			large_fam_small_fam_seqcount_dict[large_famid][small_famid]=1
		
	return(large_fam_small_fam_seqcount_dict)

def get_family_mergings(large_fam_small_fam_seqcount_dict, famid_famsize_dict, small_fam_large_fam_overlap_cutoff):
	print "#fam1 fam2 fam1-fam2_overlap fam1_size fam2_size"
	for large_fam in large_fam_small_fam_seqcount_dict:
		for small_fam in large_fam_small_fam_seqcount_dict[large_fam]:
			small_fam_seqcount=large_fam_small_fam_seqcount_dict[large_fam][small_fam]
			large_fam_size = famid_famsize_dict[large_fam]
			small_fam_size = famid_famsize_dict[small_fam]
			small_fam_large_fam_overlap = small_fam_seqcount/small_fam_size
			if (small_fam_large_fam_overlap>small_fam_large_fam_overlap_cutoff):
				print '{0} {1} {2} {3} {4}'.format(large_fam, small_fam, small_fam_large_fam_overlap, large_fam_size, small_fam_size)
##########################################################################################################################################

fam_fasta_dirName="/data/family_fasta/"
predicted_missing_seqs_fileName="/data/combined-family.predicted_missing_sequences"
small_fam_large_fam_overlap_cutoff=0.5

seqid_famid_dict, famid_famsize_dict = fam_fasta_dicts(fam_fasta_dirName)
large_fam_small_fam_seqcount_dict = get_large_fam_small_fam_seqcount_dict(predicted_missing_seqs_fileName, seqid_famid_dict)
get_family_mergings(large_fam_small_fam_seqcount_dict, famid_famsize_dict, small_fam_large_fam_overlap_cutoff)

