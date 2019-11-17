#!/usr/bin/env python

import re
import sys
import os
import subprocess


def build_family_msa(fam_fasta_dir, fam_msa_dir, fam_id_arr):
	for fam_fasta_filename in fam_id_arr:
		print 'building msa for family {0}'.format(fam_fasta_filename)
		fam_msa_outfile = open(fam_msa_dir+"/"+fam_fasta_filename, "w")
		run_mafft = run_mafft = subprocess.Popen(["mafft", "--auto", "--amino", fam_fasta_dir+"/"+fam_fasta_filename], stdout=fam_msa_outfile, stderr=subprocess.PIPE)
		run_mafft.communicate()
		fam_msa_outfile.close

def build_family_hmm(fam_msa_dir, fam_hmm_model_dir, fam_id_arr):
	for fam_msa_filename in fam_id_arr:
		print 'building hmm for family {0}'.format(fam_msa_filename)
		run_hmmbuild = subprocess.Popen(["hmmbuild", "-n", fam_msa_filename, fam_hmm_model_dir+"/"+fam_msa_filename, fam_msa_dir+"/"+fam_msa_filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		run_hmmbuild.communicate()

def build_hmm_database(outpath, fam_hmm_model_dir, fam_hmm_database_filename, fam_id_arr):
	concatenate_family_hmm_model_files(outpath, fam_hmm_model_dir, fam_hmm_database_filename, fam_id_arr)
	run_hmmpress =  subprocess.Popen(["hmmpress", "-f", fam_hmm_database_filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	run_hmmpress.communicate()

def concatenate_family_hmm_model_files(outpath, fam_hmm_model_dir, fam_hmm_database_filename, fam_id_arr):
	hmm_database_file = open(fam_hmm_database_filename, "w")
	for hmm_model_filename in fam_id_arr:
		hmm_model_file = open(fam_hmm_model_dir+"/"+hmm_model_filename)
		for line in hmm_model_file:
			hmm_database_file.write(line)
		hmm_model_file.close()
	hmm_database_file.close()
	
def get_sequence_dictionaries(fam_fasta_dir, hmmsearch_fasta_database_filename):
	seqid_famid_dict = {}
	seqid_sequence_dict = {}
	for fam_fasta_filename in os.listdir(fam_fasta_dir):
		fam_fasta_file = open(fam_fasta_dir+"/"+fam_fasta_filename, "r")
		for line in fam_fasta_file:
			line = line.rstrip()
			if(re.match(r'^\>', line)):
				linearr = re.split(r'\s+', line)
				seqid = linearr[0]
				seqid = seqid[1:]
				seqid_famid_dict[seqid] = fam_fasta_filename
				seqid_sequence_dict[seqid] = ''
			else:
				seqid_sequence_dict[seqid] = seqid_sequence_dict[seqid]+line
		
		fam_fasta_file.close()
	
	seqid_sequence_dict = update_seqid_sequence_dict(hmmsearch_fasta_database_filename, seqid_sequence_dict)
	return([seqid_famid_dict, seqid_sequence_dict])

def update_seqid_sequence_dict(hmmsearch_fasta_database_filename, seqid_sequence_dict):
	hmmsearch_fasta_database_file = open(hmmsearch_fasta_database_filename, "r")
	for line in hmmsearch_fasta_database_file:
		line = line.rstrip()
		if(re.match(r'^\>', line)):
			linearr = re.split(r'\s+', line)
			seqid = linearr[0]
			seqid = seqid[1:]
			seqid_sequence_dict[seqid] = ''
		else:
			seqid_sequence_dict[seqid] = seqid_sequence_dict[seqid]+line
	
	hmmsearch_fasta_database_file.close()
	return(seqid_sequence_dict)
		

def create_family_results_dir(outpath, fam_id):
    if os.path.exists(outpath+"/"+fam_id):
        return 1
    else:
        os.makedirs(outpath+"/"+fam_id)

def execute_hmmsearch(outpath, fam_hmm_model_dir, hmmsearch_fasta_database_filename, fam_id, hmmsearch_num_threads):
	run_hmmsearch = subprocess.Popen(["hmmsearch", "--tblout", outpath+"/"+fam_id+"/"+fam_id + ".hmmsearch", "--noali", "-E", "1e-5", "--cpu", hmmsearch_num_threads, fam_hmm_model_dir+"/"+fam_id, hmmsearch_fasta_database_filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	run_hmmsearch.communicate()


def get_closer_ingroup_sequences(fam_id, outpath, outgroup_id_arr, seqid_famid_dict):
	hmmsearch_file = open(outpath+"/"+fam_id+"/"+fam_id + ".hmmsearch", "r")
	close_seqid_arr=list()
	outgroup_seq_flag = 0
	for line in hmmsearch_file:
		line = line.rstrip()
		if (re.match(r'^\#', line)):
			continue
		linearr = re.split(r'\s+', line)
		subject_seqid = linearr[0]
		if not (check_outgroup_sequence(subject_seqid, outgroup_id_arr)):
			if(seqid_famid_dict.has_key(subject_seqid)):
				if not (seqid_famid_dict[subject_seqid] == fam_id):
					close_seqid_arr.append(subject_seqid)
			else:
				close_seqid_arr.append(subject_seqid)
		else:
			outgroup_seq_flag = 1
			break
	hmmsearch_file.close()

	if(outgroup_seq_flag == 1):
		return (close_seqid_arr)
	else:
		close_seqid_arr = list()
		return (close_seqid_arr)

def check_outgroup_sequence(seqid, outgroup_id_arr):
	detection_flag = 0
	for outgroup_id in outgroup_id_arr:
		if(re.match(outgroup_id, seqid)):
			detection_flag = 1
			break

	return(detection_flag)

def write_closer_sequences_fasta(fam_id, outpath, close_seqid_arr, seqid_sequence_dict):
	closer_sequences_fasta_file = open(outpath+"/"+fam_id+"/"+fam_id+".closer_sequences", "w")
	for seqid in close_seqid_arr:
		closer_sequences_fasta_file.write(">"+seqid+"\n")
		closer_sequences_fasta_file.write(seqid_sequence_dict[seqid]+"\n")
	
	closer_sequences_fasta_file.close()

def execute_hmmscan(fam_id, outpath, fam_hmm_database_filename, hmmscan_num_threads):
	run_hmmscan = subprocess.Popen(["hmmscan", "--tblout", outpath+"/"+fam_id+"/"+fam_id + ".hmmscan", "--noali", "-E", "1e-5", "--cpu", hmmscan_num_threads, fam_hmm_database_filename, outpath+"/"+fam_id+"/"+fam_id+".closer_sequences"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	run_hmmscan.communicate()

def read_hmmscan_output(fam_id, outpath):
	hmmscan_dict={}
	hmmscan_output_file = open(outpath+"/"+fam_id+"/"+fam_id + ".hmmscan", "r")
	for line in hmmscan_output_file:
		line = line.rstrip()
		if(re.match(r'^\#', line)):
			continue
		linearr = re.split(r'\s+', line)
		query_id = linearr[2]
		target_id = linearr[0]
		if(hmmscan_dict.has_key(query_id)):
			hmmscan_dict[query_id]["targets"].append(target_id)
		else:
			hmmscan_dict[query_id]={}
			hmmscan_dict[query_id]["targets"]=list()
			hmmscan_dict[query_id]["targets"].append(target_id)
	
	hmmscan_output_file.close()
	return(hmmscan_dict)

def predict_missing_sequences(fam_id, hmmscan_dict, seqid_famid_dict):
	missing_seqs_arr = list()
	for seqid in hmmscan_dict:
		target_arr = hmmscan_dict[seqid]["targets"]
		if(seqid_famid_dict.has_key(seqid)):
			if(seqid_famid_dict[seqid]==target_arr[0]):
				if(fam_id == target_arr[1]):
					missing_seqs_arr.append(seqid)
			
			elif(fam_id == target_arr[0]):
				missing_seqs_arr.append(seqid)

		elif(fam_id == target_arr[0]):
			missing_seqs_arr.append(seqid)
	
	return(missing_seqs_arr)

def print_missing_sequences_to_file(fam_id, outpath, missing_seqs_arr):
	missing_sequences_file = open(outpath+"/"+fam_id+"/"+fam_id+".predicted_missing_sequences", "w")
	for seqid in missing_seqs_arr:
		missing_sequences_file.write(seqid+"\n")
	missing_sequences_file.close()

def execute_worflow_for_family(fam_id, outpath, fam_fasta_dir, fam_msa_dir, fam_hmm_model_dir, hmmsearch_fasta_database_filename, fam_hmm_database_filename, outgroup_id_arr, hmmsearch_num_threads, hmmscan_num_threads):
	print 'processing family {0}'.format(fam_id)
	seqid_famid_dict, seqid_sequence_dict = get_sequence_dictionaries(fam_fasta_dir, hmmsearch_fasta_database_filename)
	create_family_results_dir(outpath, fam_id)
	execute_hmmsearch(outpath, fam_hmm_model_dir, hmmsearch_fasta_database_filename, fam_id, hmmsearch_num_threads)
	close_seqid_arr = get_closer_ingroup_sequences(fam_id, outpath, outgroup_id_arr, seqid_famid_dict)
	if(len(close_seqid_arr)==0):
		print "No ingroup sequences closer than outgroup sequences found"
		return 0
	write_closer_sequences_fasta(fam_id, outpath, close_seqid_arr, seqid_sequence_dict)
	execute_hmmscan(fam_id, outpath, fam_hmm_database_filename, hmmscan_num_threads)
	hmmscan_dict = read_hmmscan_output(fam_id, outpath)
	missing_seqs_arr = predict_missing_sequences(fam_id, hmmscan_dict, seqid_famid_dict)
	if(len(missing_seqs_arr)==0):
		print "No missing sequences found"
		return 0
	print_missing_sequences_to_file(fam_id, outpath, missing_seqs_arr)
	
def execute_worflow_for_family_arr(fam_id_arr, outpath, fam_fasta_dir, fam_msa_dir, fam_hmm_model_dir, hmmsearch_fasta_database_filename, fam_hmm_database_filename, outgroup_id_arr, hmmsearch_num_threads, hmmscan_num_threads):
	for fam_id in fam_id_arr:
		execute_worflow_for_family(fam_id, outpath, fam_fasta_dir, fam_msa_dir, fam_hmm_model_dir, hmmsearch_fasta_database_filename, fam_hmm_database_filename, outgroup_id_arr, hmmsearch_num_threads, hmmscan_num_threads)

def get_fam_id_arr(fam_fasta_dir, min_fam_size):
	fam_id_arr = list()
	for fam_fasta_filename in os.listdir(fam_fasta_dir):
		fam_size = get_famsize(fam_fasta_dir+"/"+fam_fasta_filename)
		if(fam_size>=min_fam_size):
			fam_id_arr.append(fam_fasta_filename)

	return(fam_id_arr)
	

def get_famsize(fam_fasta_filename):
	fam_size = 0
	fam_fasta_file = open(fam_fasta_filename, "r")
	for line in fam_fasta_file:
		line = line.rstrip()
		if(re.match(r'^\>', line)):
			fam_size+=1
	fam_fasta_file.close()
	return(fam_size)
########################################################################################################################################
def get_outgroup_ids(outgroup_ids_filename):
	outgroup_id_arr = list()
	outgroup_ids_file = open(outgroup_ids_filename, "r")
	for line in outgroup_ids_file:
		line = line.rstrip()
		outgroup_id_arr.append(line)
	outgroup_ids_file.close()
	return(outgroup_id_arr) 
#########################################################################################################################################
outpath = "/data/results/"
fam_fasta_dir = "/data/family_fasta/"
fam_msa_dir = "/data/family_msa/"
fam_hmm_model_dir = "/data/family_hmm/"
hmmsearch_fasta_database_filename = "/data/family-outgroup-sequences.fa"
fam_hmm_database_filename = "/data/combined-family.hmm"
outgroup_ids_filename = "/data/outgroups.list"
min_fam_size = 2
hmmsearch_num_threads = "10"
hmmscan_num_threads = "10"

outgroup_id_arr = get_outgroup_ids(outgroup_ids_filename)

fam_id_arr = get_fam_id_arr(fam_fasta_dir, min_fam_size)
print 'Total number of families {0}'.format(len(fam_id_arr))

build_family_msa(fam_fasta_dir, fam_msa_dir, fam_id_arr)
build_family_hmm(fam_msa_dir, fam_hmm_model_dir, fam_id_arr)
build_hmm_database(outpath, fam_hmm_model_dir, fam_hmm_database_filename, fam_id_arr)


execute_worflow_for_family_arr(fam_id_arr, outpath, fam_fasta_dir, fam_msa_dir, fam_hmm_model_dir, hmmsearch_fasta_database_filename, fam_hmm_database_filename, outgroup_id_arr, hmmsearch_num_threads, hmmscan_num_threads)
