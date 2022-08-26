#!/usr/bin/env python3
import sys
import math
import numpy as np
import argparse

#  
# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
# ===========================================================================
#
#

def get_args():
	parser = argparse.ArgumentParser(description="Post processing of BLAST output. Returns coverage of query input.")
	parser.add_argument("-s", type=argparse.FileType('r', encoding='UTF-8'), help="Please provide the BLAST output for a given accession", required=True)
	parser.add_argument('-p', type=int, help="perc identity threshold", default=90)
	parser.add_argument('-a', type=int, help="perc alignment length threshold", default=50)
	parser.add_argument('-q', type=str, help="input query", required=True)
	return parser.parse_args()


def GetCoverage(querylen_arr, last_qlen):
	qbp=0
	for y in range(len(querylen_arr)):
		if querylen_arr[y]>0:
			qbp=qbp+1
	percq_cov = (qbp/last_qlen)*100
	assert(qbp <= last_qlen)
	#print("qbp:", qbp) 
	#print("{0} {1} {2} {3} {4} {5:4.2f}".format(last_query, srr_id, good_matches, last_qlen, last_slen, percq_cov))
	#if last_query == myquery:
	#	outfilename=srr_id+"_qcov.txt"
	#	fout = open(outfilename, 'w')
	#	for y in range(len(querylen_arr)):
	#		fout.write("{0} {1}\n".format(y, querylen_arr[y]))
	#	fout.close()
	return qbp,percq_cov


# parse tabular file using -outfmt '6 std qlen slen qcovs'
# Read in 15 fields.  Look at query and subject length
def GetMatches(f, ident_cutoff, alignlen_cutoff, srr_id, myquery):
	last_subj = "NONE"
	last_query = "NONE"
	good_matches = 0
	my_qlen = 0
	percq_cov=0.0
	start = 0
	stop = 0
	myquery_exists = "NONE" 
	for x in f:
		line = x.split()
		subj = line[0]
		query = line[1]
		ident = float(line[2])
		alignlen = int(line[3])
		qstart = int(line[8])
		qend = int(line[9])
		slen = int(line[12])
		qlen = int(line[13])
		scov = float(line[14])
# only first HSP per subject
		if last_subj == subj:
		    continue
		#align_len = math.floor(slen*(alignlen_cutoff/100))
		start = qstart if qstart < qend else qend
		stop = qend if qstart < qend else qstart
		if query == myquery:
			if myquery_exists == "NONE":
				my_qlen = qlen
				myquery_exists = myquery
				querylen_arr = np.zeros(qlen+1,dtype=int)
			if ident >= ident_cutoff and scov >= alignlen_cutoff:
				#print("{0} {1} {2:4.2f} {3} {4} {5} {6} {7}".format(query, subj, ident, alignlen, qlen, slen, start, stop))
				good_matches = good_matches + 1
				for x in range(start,stop+1):
					querylen_arr[x] = querylen_arr[x] + 1
		last_subj = line[0]
		last_query = line[1]

	#printing the summary for myquery
	if myquery_exists == myquery:
		qbp, percq_cov = GetCoverage(querylen_arr, my_qlen)
		print("query, accession, number-of-good-matches, query-length, query-bp-covered, perc-query-coverage")
		print("{0} {1} {2} {3} {4} {5:4.2f}\n".format(myquery, srr_id, good_matches, my_qlen, qbp, percq_cov))
	if myquery_exists == "NONE":
		print("BLAST output for accession {0} did not find any matches for {1}".format(srr_id, myquery))


def main():
	args = get_args()
	f1 = open(args.s.name, "r")
	first=f1.readline()
	for x in f1:
		line = x.split()
		subj = line[0]
		srr_line=subj.split(".")
		srr_id=srr_line[0]
		break
	f1.close()
	f2 = open(args.s.name, "r")
	perc_ident = args.p
	perc_alignment_length = args.a
	user_query=(args.q).strip()
	GetMatches(f2, perc_ident, perc_alignment_length, srr_id, user_query)
	f2.close()
	

if __name__ == "__main__":
	main()
