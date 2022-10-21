#! /usr/bin/env python3

import sys
import argparse
from Bio import SeqIO
import pprint
pp = pprint.PrettyPrinter(indent = 2)


class chain_data():
	def __init__(self, tName, tSize, tStrand, tStart, tEnd, qName, qSize, qStrand, qStart, qEnd, id):
		self.tName   = tName
		self.tSize   = tSize
		self.tStrand = tStrand
		self.tStart  = tStart
		self.tEnd    = tEnd
		self.qName   = qName
		self.qSize   = qSize
		self.qStrand = qStrand
		self.qStart  = qStart
		self.qEnd    = qEnd
		self.id      = id
	def __str__(self):
		return "\t".join([str(x) for x in [self.tName, self.tSize, self.tStrand, self.tStart, self.tEnd, self.qName, self.qSize, self.qStrand, self.qStart, self.qEnd, self.id]])


def chain_parser(filename):
	ret_array      = []
	score          = 0
	tName          = ""
	tSize          = 0
	tStrand        = ""
	tStart         = 0
	tEnd           = 0
	qName          = ""
	qSize          = 0
	qStrand        = ""
	qStart         = 0
	qEnd           = 0
	id_            = 0
	current_chains = []
	with open(filename) as f:
		each_line = f.readlines()
		for index, line in enumerate(each_line):
			if line.startswith("chain") or index == len(each_line) - 1:#最後のレコードだけはどうしようもない。EOFを検知して処理にすすむ
				if current_chains != []:#chainsがすでにあるとき
					#chainsで漸化式を回す処理
					if qStrand == "+":
						ts = tStart
						qs = qStart
						for i in range(len(current_chains)):
							tmp_chain_data = chain_data(
								tName   = tName,
								tSize   = tSize,
								tStrand = tStrand,
								tStart  = ts,
								tEnd    = ts + current_chains[i][0],
								qName   = qName,
								qSize   = qSize,
								qStrand = qStrand,
								qStart  = qs,
								qEnd    = qs + current_chains[i][0],
								id      = id_
							)
							ret_array.append(tmp_chain_data)
							if i != len(current_chains) - 1:
								ts = ts + current_chains[i][0] + current_chains[i][1]
								qs = qs + current_chains[i][0] + current_chains[i][2]
							else:
								pass
					if qStrand == "-":
						ts = tStart
						qe = qSize - qStart
						for i in range(len(current_chains)):
							tmp_chain_data = chain_data(
								tName   = tName,
								tSize   = tSize,
								tStrand = tStrand,
								tStart  = ts,
								tEnd    = ts + current_chains[i][0],
								qName   = qName,
								qSize   = qSize,
								qStrand = qStrand,
								qStart  = qe - current_chains[i][0],
								qEnd    = qe,
								id      = id_
							)
							ret_array.append(tmp_chain_data)
							if i != len(current_chains) - 1:
								ts = ts + current_chains[i][0] + current_chains[i][1]
								qe = qe - current_chains[i][0] - current_chains[i][2]
							else:
								pass
				else:#chainsが入ってないとき（ファイルの先頭）
					pass
				#chainをつないだら、初項にあたる情報を初期化
				if index == len(each_line) - 1:
					break
				chain, score, tName, tSize, tStrand, tStart, tEnd, qName, qSize, qStrand, qStart, qEnd, id_ = line.strip().split()
				score   = int(score)
				tSize   = int(tSize)
				tStart  = int(tStart)
				tEnd    = int(tEnd)
				qSize   = int(qSize)
				qStart  = int(qStart)
				qEnd    = int(qEnd)
				current_chains = []
			elif line != "\n":
				current_chains.append([int(x) for x in line.strip().split("\t")])
	return ret_array

if __name__ == "__main__":
	for i in chain_parser("./small_chain2"):
		print(i)