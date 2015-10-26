#!/usr/bin/python
# coding: UTF-8

import os,sys,errno
from optparse import Option, OptionParser

__version__ = "1.0.0"

def start():
	parser = OptionParser(usage="usage: %prog [options] input_network_file", version="controllability version %s" % __version__,description="Determine the Driver and the Critical nodes of a network with the maximum matching. This script is an implementation of the Nature paper by Liu et al. in 2011; Controllability of complex networks.")
	
	(options, args) = parser.parse_args()

	if len(args) != 1:
		parser.print_help()
		sys.exit()

	input_network_file, = args

	create_result_dir()

	graph, Tnum, RUck, RUnode, Vnode, Unode, numnode, interaction, nodes = create_network(input_network_file)
	bip = bipartiteMatch(graph)
	output_results(input_network_file, bip, Tnum, RUck, RUnode, Vnode, Unode, numnode, interaction, nodes)


def create_result_dir():
	if os.path.exists('results/'):
		msg = 'Directory results already exists.'
		print msg
	else:
		os.mkdir('results/')

def get_progressbar_str(progress):
	MAX_LEN = 30
	BAR_LEN = int(MAX_LEN * progress)
	return ('[' + '=' * BAR_LEN +
			('>' if BAR_LEN < MAX_LEN else '') +
			' ' * (MAX_LEN - BAR_LEN) +
			'] %.1f%%' % (progress * 100.))


def create_network(input_network_file):
	nodeslist = {}
	interaction = {}

	with open(input_network_file, 'r') as fi:
		line = fi.readline()
		while line:
			itemList = line[:-1].split('\t')
			nodeslist[itemList[0]] = 1
			nodeslist[itemList[1]] = 1
			if itemList[0] not in interaction:
				interaction[itemList[0]] = {}
			interaction[itemList[0]][itemList[1]] = 1
			line = fi.readline()

	nodes = []
	for item in nodeslist:
		nodes.append(item)

	numnode = len(nodes)

	RUck = {}
	RUnode = {}
	Tnum = numnode
	graph = {}
	Vnode = {}
	Unode = {}

	for i in range(numnode):
		Vnode[nodes[i]] = i
		Unode[nodes[i]] = i + numnode
		RUck[i + numnode] = 0
		RUnode[i + numnode] = nodes[i]

	for temp1 in interaction:
		for temp2 in interaction[temp1]:
			tempNode1 = Vnode[temp1]
			tempNode2 = Unode[temp2]
			if tempNode2 not in graph:
				graph[tempNode2] = {}
			graph[tempNode2][tempNode1] = 1
			print(str(tempNode1) + "\t" + str(tempNode2) + "\t" + temp1 + "\t" + temp2)

	return graph, Tnum, RUck, RUnode, Vnode, Unode, numnode, interaction, nodes

# Hopcroft-Karp bipartite max-cardinality matching and max independent set
# David Eppstein, UC Irvine, 27 Apr 2002

def bipartiteMatch(graph):
	'''Find maximum cardinality matching of a bipartite graph (U,V,E).
	The input format is a dictionary mapping members of U to a list
	of their neighbors in V.  The output is a triple (M,A,B) where M is a
	dictionary mapping members of V to their matches in U, A is the part
	of the maximum independent set in U, and B is the part of the MIS in V.
	The same object may occur in both U and V, and is treated as two
	distinct vertices if this happens.'''
	
	# initialize greedy matching (redundant, but faster than full search)
	matching = {}
	for u in graph:
		for v in graph[u]:
			if v not in matching:
				matching[v] = u
				break
	
	while 1:
		# structure residual graph into layers
		# pred[u] gives the neighbor in the previous layer for u in U
		# preds[v] gives a list of neighbors in the previous layer for v in V
		# unmatched gives a list of unmatched vertices in final layer of V,
		# and is also used as a flag value for pred[u] when u is in the first layer
		preds = {}
		unmatched = []
		pred = dict([(u,unmatched) for u in graph])
		for v in matching:
			del pred[matching[v]]
		layer = list(pred)
		
		# repeatedly extend layering structure by another pair of layers
		while layer and not unmatched:
			newLayer = {}
			for u in layer:
				for v in graph[u]:
					if v not in preds:
						newLayer.setdefault(v,[]).append(u)
			layer = []
			for v in newLayer:
				preds[v] = newLayer[v]
				if v in matching:
					layer.append(matching[v])
					pred[matching[v]] = v
				else:
					unmatched.append(v)
		
		# did we finish layering without finding any alternating paths?
		if not unmatched:
			unlayered = {}
			for u in graph:
				for v in graph[u]:
					if v not in preds:
						unlayered[v] = None
			return (matching,list(pred),list(unlayered))

		# recursively search backward through layers to find alternating paths
		# recursion returns true if found path, false otherwise
		def recurse(v):
			if v in preds:
				L = preds[v]
				del preds[v]
				for u in L:
					if u in pred:
						pu = pred[u]
						del pred[u]
						if pu is unmatched or recurse(pu):
							matching[v] = u
							return 1
			return 0

		for v in unmatched: recurse(v)

def output_results(input_network_file, bip, Tnum, RUck, RUnode, Vnode, Unode, numnode, interaction, nodes):

	value = bip[0].values()
	key = bip[0].keys()

	input_file = input_network_file.split('/')[-1]

	outtest = os.path.join('results',input_file.replace('.txt', '_test.txt'))
	outNdrivers = os.path.join('results',input_file.replace('.txt', '_NonDriver.txt'))
	outDrivers = os.path.join('results',input_file.replace('.txt', '_Driver.txt'))
	outreport = os.path.join('results',input_file.replace('.txt', '_Analysis_report.txt'))
	outnodes = os.path.join('results',input_file.replace('.txt', '_List_Nodes.txt'))

	numS1 = len(value)
	numS2 = len(key)

	with open(outtest, 'w') as fot:
		for i in range(numS1):
			fot.write(str(key[i]) + "\n")

	with open(outNdrivers, 'w') as fond:
		for i in range(numS1):
			fond.write(RUnode[value[i]] + "\n")
			RUck[value[i]] = 1

	with open(outDrivers, 'w') as fod:
		for i in range(numnode):
			if RUck[i + numnode] != 1:
				fod.write(RUnode[i + numnode] + "\n")

	numD = len(key)
	ND = Tnum - numD
	FD = float(ND) / Tnum

	with open(outreport, 'w') as fore:
		fore.write("Number of driver nodes:\t" + str(ND) + "\n")
		fore.write("Fraction of driver nodes:\t" + str(FD) + "\n")
		fore.write("Total number of nodes:\t" + str(Tnum) + "\n")

	NCN = 0
	NRN = 0
	NON = 0
	criticality = {}

	for j in range(numnode):

		graph = {}

		for temp1 in interaction:
			if temp1 != nodes[j]:
				for temp2 in interaction[temp1]:
					if temp2 != nodes[j]:
						tempNode1 = Vnode[temp1]
						tempNode2 = Unode[temp2]
						if tempNode2 not in graph:
							graph[tempNode2] = {}
						graph[tempNode2][tempNode1] = 1

		bip = bipartiteMatch(graph)

		value = bip[0].values()
		key = bip[0].keys()

		numD = len(key)
		TND = Tnum - numD - 1
		criticality[nodes[j]] = TND - ND
		if TND > ND:
			NCN += 1
		elif TND == ND:
			NON += 1
		else:
			NRN += 1

	FCR = float(NCN) / Tnum
	FON = float(NON) / Tnum
	FRN = float(NRN) / Tnum

	with open(outreport, 'a') as fore:
		fore.write("Number of critical nodes:\t" + str(NCN) + "\n")
		fore.write("Number of redundant nodes:\t" + str(NRN) + "\n")
		fore.write("Number of ordinary nodes:\t" + str(NON) + "\n")
		fore.write("Fraction of critical nodes:\t" + str(FCR) + "\n")
		fore.write("Fraction of redundant nodes:\t" + str(FRN) + "\n")
		fore.write("Fraction of ordinary nodes:\t" + str(FON) + "\n")

	with open(outnodes, 'w') as fon:
		for i in range(numnode):
			if criticality[nodes[i]] > 0:
				fon.write(nodes[i] + "\tC\t" + str(criticality[nodes[i]]) + "\n")
			elif criticality[nodes[i]] == 0:
				fon.write(nodes[i] + "\tO\t" + str(criticality[nodes[i]]) + "\n")
			else:
				fon.write(nodes[i] + "\tR\t" + str(criticality[nodes[i]]) + "\n")

if __name__ == "__main__":
	try:
		start()
	except KeyboardInterrupt:
		pass
	except IOError as e:
		if e.errno == errno.EPIPE:
			pass
		else:
			raise

