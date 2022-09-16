#! /usr/bin/env python3

import sys
import argparse
from Bio import SeqIO
import pprint
from chain_parser import chain_data
from chain_parser import chain_parser


def main():
	parser = argparse.ArgumentParser(description = "compare two regions mapped by chain file each other")
	parser.add_argument('genome1', metavar='genome1',   type=str, help='genome1 file name')
	parser.add_argument('genome2', metavar='genome2',   type=str, help='genome2 file name')
	parser.add_argument('chain',   metavar='chain', type=str, help='chain file name(genome1Togenome2.chain)')

	print(f"start  reading files", file = sys.stderr)
	genome1_seq = SeqIO.to_dict(SeqIO.parse(args.genome1, "fasta"))
	genome2_seq = SeqIO.to_dict(SeqIO.parse(args.genome2, "fasta"))
	print(f"finish reading {args.genome2}", file = sys.stderr)
	print(f"finish reading {args.chain}", file = sys.stderr)

	different_cnt = 0
	pre = args.report_prefix
	failure_report_file               = open(f"{pre}_fail.txt", "w")
	mismatch_report_output_file       = open(f"{pre}_mismatch_rate.tsv", "w")
	exact_mismatch_report_output_file = open(f"{pre}_completely_mismatch", "w")
	final_report_output_file          = open(f"{pre}_result", "w")
	total_n_of_chain = len(chain)
	print("len\tmismatch\tsrc_loci\ttgt_loci", file = mismatch_report_output_file)

	for i, c in enumerate(chain):
		target_chrom = c.tName
		target_start = c.tStart
		target_end   = c.tEnd
		source_chrom = c.qName
		source_start = c.qStart
		source_end   = c.qEnd
		rev          = c.qStrand

		print(f"\r{i}/{total_n_of_chain}", end = "", file = sys.stderr)
		try:
			target_seq = str(genome2_seq[target_chrom].seq[target_start:target_end]).lower()
			if rev == "+":
				source_seq = str(genome1_seq[source_chrom].seq[source_start:source_end]).lower()
			else:
				source_seq = str(genome1_seq[source_chrom].seq[source_start:source_end].reverse_complement()).lower()
			if target_seq != source_seq:
				different_cnt      = different_cnt + 1
				different_base_cnt = 0
				f = False
				for i in range(len(target_seq)):
					if target_seq[i] == source_seq[i]:
						different_base_cnt = different_base_cnt + 1
				print(f"{different_base_cnt}\t{len(target_seq)}\t{target_chrom}:{target_start}-{target_end}\t{source_chrom}:{source_start}-{source_end}", file = final_report_output_file)
				if target_end - target_start == different_base_cnt:
					print(f">genome2, {target_chrom}:{target_start}-{target_end}\t{target_seq}",   file = exact_mismatch_report_output_file)
					print(f">genome1, {source_chrom}:{source_start}-{source_end}\t{source_seq}\n", file = exact_mismatch_report_output_file)
			else:
				print(f"0\t{len(target_seq)}\t{target_chrom}:{target_start}-{target_end}\t{source_chrom}:{source_start}-{source_end}", file = final_report_output_file)
				same_cnt = same_cnt + 1
		except Exception as e:
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