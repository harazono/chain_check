#! /usr/bin/env python3

import sys
import argparse
from Bio import SeqIO
from Bio import pairwise2
import pprint
pp = pprint.PrettyPrinter(indent = 2)


class chain_data():
	def __init__(self, source_chrom, source_start, source_end, target_chrom, target_start, target_end, cigar, rev, dsc = None, color = None):
		self.source_chrom = source_chrom
		self.source_start = source_start
		self.source_end = source_end
		self.target_chrom = target_chrom
		self.target_start = target_start
		self.target_end = target_end
		self.cigar = cigar
		self.rev = rev
		self.dsc = dsc

	def __str__(self):
		return f"source_chrom: {self.source_chrom}, source_start: {self.source_start}, source_end: {self.source_end}, target_chrom: {self.target_chrom}, target_start: {self.target_start}, target_end: {self.target_end}, cigar: {self.cigar}, rev: {self.rev}, dsc: {self.dsc}"


def chain_parser(filename):
	ret_array = []
	score = 0
	tName = ""
	tSize = 0
	tStart = 0
	tEnd = 0
	qName = ""
	qSize = 0
	qStrand = ""
	qStart = 0
	qEnd = 0
	id_ = 0
	source_current_start_pos = 0
	target_current_start_pos = 0
	with open(filename) as f:
		for line in f:
			if line.startswith("chain"):
				chain, score, tName, tSize, tStrand, tStart, tEnd, qName, qSize, qStrand, qStart, qEnd, id_ = line.strip().split()
				source_current_start_pos = int(tStart)
				if qStrand == "+":#forward strand
					target_current_start_pos = int(qStart)
				else:#reverse strand
					target_current_start_pos = int(qSize) - int(qStart)#右端。表示用ということで、画面の左端側が開始、右端側が終了という風に扱う。

			elif line != "\n":
				block = line.strip().split()
				source_start = source_current_start_pos #sourceは常にforward strandと想定。assume that sStrand is always "+"
				source_end   = source_start + int(block[0])
				if qStrand == "+":
					target_start = target_current_start_pos
					target_end   = target_start + int(block[0])
				else:
					#assert(int(qEnd) > int(qStart), "unexpected")
					target_end   = target_current_start_pos
					target_start = target_end - int(block[0])

				tmp_chain_data = chain_data(
					target_chrom = tName,
					target_start = target_start,
					target_end = target_end,
					source_chrom = qName,
					source_start = source_start,
					source_end = source_end,
					cigar = target_end - target_start,
					rev = qStrand,
					dsc = filename
				)
				assert(target_end - target_start == source_end - source_start), "chain file reports different length between reference and target"
				ret_array.append(tmp_chain_data)
				if len(block) == 3 and qStrand == "+":
					source_current_start_pos = source_end + int(block[1])
					target_current_start_pos = target_end + int(block[2])
				if len(block) == 3 and qStrand == "-":
					source_current_start_pos = source_end + int(block[1])
					target_current_start_pos = target_end - int(block[0]) - int(block[2])


	return ret_array


def main():
	parser = argparse.ArgumentParser(description = "compare two regions mapped by chain file each other")
	parser.add_argument('genome1', metavar='genome1',   type=str, help='genome1 file name')
	parser.add_argument('genome2', metavar='genome2',   type=str, help='genome2 file name')
	parser.add_argument('chain', metavar='chain', type=str, help='chain file name(genome1Togenome2.chain)')
	parser.add_argument('report_prefix', metavar='report_prefix', type=str, help='file prefix for report')
	args = parser.parse_args()

	print(f"start  reading files", file = sys.stderr)
	genome1_seq = SeqIO.to_dict(SeqIO.parse(args.genome1, "fasta"))
	print(f"finish reading {args.genome1}", file = sys.stderr)
	genome2_seq = SeqIO.to_dict(SeqIO.parse(args.genome2, "fasta"))
	print(f"finish reading {args.genome2}", file = sys.stderr)
	chain    = chain_parser(args.chain)
	print(f"finish reading {args.chain}", file = sys.stderr)

	chain_result = []
	same_cnt = 0
	different_cnt = 0
	pre = args.report_prefix
	failure_report_file               = open(f"{pre}_fail.txt", "w")
	mismatch_report_output_file       = open(f"{pre}_mismatch_rate.tsv", "w")
	exact_mismatch_report_output_file = open(f"{pre}_completely_mismatch", "w")
	final_report_output_file          = open(f"{pre}_result", "w")
	total_n_of_chain = len(chain)
	print("len\tmismatch\tsrc_loci\ttgt_loci", file = mismatch_report_output_file)

	for i, c in enumerate(chain):
		target_chrom = c.target_chrom
		target_start = c.target_start
		target_end   = c.target_end
		source_chrom = c.source_chrom
		source_start = c.source_start
		source_end   = c.source_end
		rev          = c.rev

		print(f"\r{i}/{total_n_of_chain}", end = "", file = sys.stderr)
		try:
			target_seq = str(genome2_seq[target_chrom].seq[target_start:target_end]).lower()
			if rev == "+":
				source_seq = str(genome1_seq[source_chrom].seq[source_start:source_end]).lower()
			else:
				source_seq = str(genome1_seq[source_chrom].seq[source_start:source_end].reverse_complement()).lower()
			if target_seq != source_seq:
				print(f"genome2:{target_chrom}:{target_start}-{target_end}\tgenome1:{source_chrom}:{source_start}-{source_end}\tare different({source_end - source_start})", file = final_report_output_file)
				different_cnt      = different_cnt + 1
				different_base_cnt = 0
				f = False
				for i in range(len(target_seq)):
					f = target_seq[i] == source_seq[i]
					if f == False:
						different_base_cnt = different_base_cnt + 1
						print(f">genome2 {target_chrom}{i}th base: {target_seq[i]}", file = final_report_output_file)
						print(f" genome1 {source_chrom}{i}th base: {source_seq[i]}", file = final_report_output_file)
				print(f"{different_base_cnt}/{target_end - target_start}({different_base_cnt/(target_end - target_start)})", file = final_report_output_file)
				print(f"{target_end - target_start}\t{different_base_cnt}\t{source_chrom}:{source_start}-{source_end}\t{target_chrom}:{target_start}-{target_end}", file = mismatch_report_output_file)
				if target_end - target_start == different_base_cnt:
					print(f">genome2, {target_chrom}:{target_start}-{target_end}\t{target_seq}",   file = exact_mismatch_report_output_file)
					print(f">genome1, {source_chrom}:{source_start}-{source_end}\t{source_seq}\n", file = exact_mismatch_report_output_file)
			else:
				print(f"genome2:{target_chrom}:{target_start}-{target_end}\tgenome1:{source_chrom}:{source_start}-{source_end}\tare same({source_end - source_start})", file = final_report_output_file)
				same_cnt = same_cnt + 1
		except Exception as e:
			#print(e, file = sys.stderr)
			print(f"{target_chrom}, {source_chrom}, {e}", file = failure_report_file)
			continue

	total = same_cnt + different_cnt
	print(f"", file = sys.stderr)
	print(f"same...{same_cnt}, different...{different_cnt}",             file = final_report_output_file)
	print(f"same...{same_cnt/total}, different...{different_cnt/total}", file = final_report_output_file)
	print(f"same...{same_cnt}, different...{different_cnt}",             file = sys.stderr)
	print(f"same...{same_cnt/total}, different...{different_cnt/total}", file = sys.stderr)


if __name__ == '__main__':
	main()