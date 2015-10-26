#!/usr/bin/python
# coding: UTF-8

import os,sys,errno
import networkx as nx
from optparse import Option, OptionParser

__version__ = "1.0.0"

def start():
	parser = OptionParser(usage="usage: %prog [options] edgelist_file ST_file species_file", version="bowtie_simple30 version %s" % __version__,description="Determine bowtie score of a network considering connecting paths between Sources (S) and Targets (T).")
	parser.add_option("--shortest",default=False,action="store_true",help="use shortest paths between Sources and Targets to calculate bowtie score")
	parser.add_option("-l","--length", default=30,help="Specify the length of simple paths to be considered. If you set the length more than 50, time for the calculation will be very long.")
	parser.add_option("-d","--DE-file",default=None,help="Specify mRNA differential expression data file. This parameter is valid only when --weighted is True")
	parser.add_option("-s", "--species-file",help="specify CellDesigner species file")
	parser.add_option("-o", "--output", help="sets the prefix of output file. If species file is given, it is not necessary, because the part of species file name will be used.")

	(options, args) = parser.parse_args()

	if len(args) != 2:
		parser.print_help()
		sys.exit()

	a, b = args

	create_result_dir()

	if options.DE_file != None:
		bowtie_w, bowtie, nonzero_connection, nodes, Sources, Targets = bowtie_weighted(options, a, b)
		if options.species_file != None:
			print_weighted_results_species(options, bowtie_w, bowtie, nonzero_connection, Sources, Targets)
		else:
			print_weighted_results(options, bowtie_w, bowtie, nonzero_connection, nodes, Sources, Targets)
	else:
		bowtie, nonzero_connection, nodes, Sources, Targets = bowtie_unweighted(options, a, b)
		if options.species_file != None:
			print_results_species(options, bowtie, nonzero_connection, Sources, Targets)
		else:
			print_results(options, bowtie, nonzero_connection, nodes, Sources, Targets)

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

def read_STfile(fname):
	Sources = []
	Targets = []
	try:
		with open(fname,'r') as fi:
			lines = fi.readlines()
			if len(lines) != 2:
				print "specify Sources and Targets in each line"
				sys.exit()
			itemList = lines[0][:-1].split('\t')
			if len(itemList) < 1:
				print "number of Sources should be > 0"
				sys.exit()
			for i in range(1,len(itemList)):
				Sources.append(itemList[i])
			if len(itemList) < 1:
				print "number of Targets should be > 0"
				sys.exit()
			itemList = lines[1][:-1].split('\t')
			for i in range(1,len(itemList)):
				Targets.append(itemList[i])
			return Sources, Targets
	except IOError as e:
		if e.errno == errno.EPIPE:
			pass
		else:
			raise

def read_DEfile(fname):
	logFC = {}
	with open(fname,'r') as fi:
		line = fi.readline()
		line = fi.readline()
		while line:
			itemList = line[:-1].split('\t')
			#print ref2eg_dic[itemList[0]]
			if itemList[1] != "NA" and itemList[1] != '':
				logFC[itemList[0]] = float(itemList[1])
			else:
				logFC[itemList[0]] = 0
			line = fi.readline()
	return logFC

def result_file_name(options):
	if options.DE_file != None:
		weight = '_weighted'
	else:
		weight = ''
	if options.species_file != None:
		species_file = options.species_file
		itemList = species_file.split('/')
		if options.shortest:
			filename = weight+'_shortest.txt'
		else:
			filename = weight+'_simple'+str(options.length)+'.txt'
		bowtieResult = os.path.join('results',itemList[-1].replace('_Species.txt',filename))
	else:
		if options.shortest:
			filename = weight+'_shortest.txt'
		else:
			filename = weight+'_simple'+str(options.length)+'.txt'
		bowtieResult = os.path.join('results',options.output+filename)

	return bowtieResult

def print_results_species(options, bowtie, nonzero_connection, Sources, Targets):
	bowtieResult = result_file_name(options)

	with open(bowtieResult,'w') as fo:
		with open(options.species_file,'r') as fi:
			line = fi.readline()
			line = fi.readline()
			fo.write("Name\tID\tType\tCompartment\tbowtie\n")
			while line:
				itemList = line[:-1].split('\t')
				if itemList[0] != 'GENE' and itemList[0] != 'RNA' and itemList[0] != 'PHENOTYPE':
					spID = itemList[1]
					if spID in Sources:
						bt = "Source"
					elif spID in Targets:
						bt = "Target"
					elif spID not in bowtie:
						bt = "0"
					else:
						bt = str(float(bowtie[spID])/nonzero_connection)
					fo.write(itemList[2] + "\t" + spID + "\t" + itemList[0] +"\t" + itemList[4] + "\t" + bt + "\n")
				line = fi.readline()

def print_weighted_results_species(options, bowtie_w, bowtie, nonzero_connection, Sources, Targets):
	bowtieResult = result_file_name(options)

	with open(bowtieResult,'w') as fo:
		with open(options.species_file,'r') as fi:
			line = fi.readline()
			line = fi.readline()
			fo.write("\t".join(["Name","ID","Type","Compartment","bowtie_weighted","bowtie","ratio\n"]))
			while line:
				itemList = line[:-1].split('\t')
				if itemList[0] != 'GENE' and itemList[0] != 'RNA':
					spID = itemList[1]
					if spID in Sources:
						bt = "Source"
						bt_w = "Source"
						ratio = "Source"
					elif spID not in bowtie:
						bt_w = "0"
						bt = "0"
						ratio = "0"
					else:
						bt_w = str(float(bowtie_w[spID])/nonzero_connection)
						bt = str(float(bowtie[spID])/nonzero_connection)
						if bowtie[spID] != 0:
							ratio = str(bowtie_w[spID]/bowtie[spID])
						else:
							ratio = "0"
					fo.write("\t".join([itemList[2],spID,itemList[0],itemList[4],bt_w,bt,ratio])+"\n")
				line = fi.readline()

def print_results(options, bowtie, nonzero_connection, nodes, Sources, Targets):
	bowtieResult = result_file_name(options)

	with open(bowtieResult,'w') as fo:
		fo.write("\t".join(["ID","bowtie\n"]))
		for node in nodes:
			if node in Sources:
				bt = "Source"
			elif node in Targets:
				bt = "Target"
			elif node not in bowtie:
				bt = "0"
			else:
				bt = str(float(bowtie[node])/nonzero_connection)
			fo.write("\t".join([node,bt])+"\n")

def print_weighted_results(options, bowtie_w, bowtie, nonzero_connection, nodes, Sources, Targets):
	bowtieResult = result_file_name(options)

	with open(bowtieResult,'w') as fo:
		fo.write("\t".join(["ID","bowtie_weighted","bowtie","ratio\n"]))
		for node in nodes:
			if node in Sources:
				bt = "Source"
				bt_w = "Source"
				ratio = "Source"
			elif node in Targets:
				bt = "Target"
				bt_w = "Target"
				ratio = "Target"
			elif node not in bowtie:
				bt = "0"
				bt_w = "0"
				ratio = "0"
			else:
				bt_w = str(float(bowtie_w[node])/nonzero_connection)
				bt = str(float(bowtie[node])/nonzero_connection)
				if bowtie[node] != 0:
					ratio = str(bowtie_w[node]/bowtie[node])
				else:
					ratio = "0"
			fo.write("\t".join([node,bt_w,bt,ratio])+"\n")

def bowtie_weighted(options, a, b):
	edgelistFile = a
	STFile = b
	DEfile = options.DE_file

	logFC = read_DEfile(DEfile)

	with open(edgelistFile,'r') as fi:
		G=nx.read_edgelist(fi, create_using=nx.DiGraph())

	nds = nx.nodes(G)

	Sources, Targets = read_STfile(STFile)

	total_sp = 0
	nonzero_connection = 0
	bowtie_w = {}
	bowtie = {}

	for tg in Targets:
		if tg not in logFC:
			logFC[tg] = 0

	x = 1.0
	total_num = len(Sources) * len(Targets)

	msg = "calculating weighted bowtie score..."
	print msg

	for sc in Sources:
		#print sc
		for tg in Targets:
			path_num = {}
			#print tg
			progress = x / total_num
			sys.stderr.write('\r\033[K' + get_progressbar_str(progress))
			sys.stderr.flush()
			x += 1
			try:
				if options.shortest:
					simple_paths = list(nx.all_shortest_paths(G, source=sc, target=tg))
				else:
					simple_paths = list(nx.all_simple_paths(G, source=sc, target=tg, cutoff=30))
			except nx.NetworkXNoPath:
				simple_paths = []
			total_path_num = len(simple_paths)
			if total_path_num != 0:
				for path in simple_paths:
					for i in range(1,len(path)-1):
						if path[i] not in path_num:
							path_num[path[i]] = 0
						path_num[path[i]] += 1
				for p in path_num:
					if p not in bowtie:
						bowtie[p] = 0
						bowtie_w[p] = 0						
					bowtie[p] += path_num[p]/total_path_num
					bowtie_w[p] += path_num[p]/total_path_num*logFC[tg]
				nonzero_connection += 1
	sys.stderr.write('\n')
	sys.stderr.flush()

	return bowtie_w, bowtie, nonzero_connection, nds, Sources, Targets

def bowtie_unweighted(options, a, b):
	edgelistFile = a
	STFile = b

	with open(edgelistFile,'r') as fi:
		G=nx.read_edgelist(fi, create_using=nx.DiGraph())

	nds = nx.nodes(G)

	Sources, Targets = read_STfile(STFile)

	total_sp = 0
	nonzero_connection = 0
	bowtie = {}

	x = 1.0
	total_num = len(Sources) * len(Targets)

	msg = "calculating bowtie score..."
	print msg

	for sc in Sources:
		for tg in Targets:
			path_num = {}
			progress = x / total_num
			sys.stderr.write('\r\033[K' + get_progressbar_str(progress))
			sys.stderr.flush()
			x += 1
			try:
				if options.shortest:
					simple_paths = list(nx.all_shortest_paths(G, source=sc, target=tg))
				else:
					simple_paths = list(nx.all_simple_paths(G, source=sc, target=tg, cutoff=30))
			except nx.NetworkXNoPath:
				simple_paths = []
			total_path_num = len(simple_paths)
			if total_path_num != 0:
				for path in simple_paths:
					for i in range(1,len(path)-1):
						if path[i] not in path_num:
							path_num[path[i]] = 0
						path_num[path[i]] += 1
				for p in path_num:
					if p not in bowtie:
						bowtie[p] = 0
					bowtie[p] += float(path_num[p])/total_path_num
				nonzero_connection += 1

	sys.stderr.write('\n')
	sys.stderr.flush()

	return bowtie, nonzero_connection, nds, Sources, Targets

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
