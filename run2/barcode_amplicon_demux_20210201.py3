#!/usr/bin/python3

licenc = '''
Copyright (C) 2020  Balint Mark Vasarhelyi
Contact: balint.mark.vasarhelyi@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>
'''

#1. demux
#PE merge, csak az atfedo readek kellenek
#q30
#
#mincov 10?
#mintankent hany read
#min var freq: 4
#
#
#2. refmap
#3. transzlacio
#4. variant call - az egyes aminosavvariansok

import csv
import time
import os
import copy
import pickle
import argparse
import glob
import itertools
from Bio import SeqIO, Seq
import multiprocessing as mp
#import more_itertools as mit
import subprocess
import pandas as pd
import sys
import random
from collections import *
import yaml
import logging
from pathlib import Path
import numpy as np
import re
import math
import gzip

parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-c','--configfile',help="Configfile",default="config.yaml")
parser.add_argument("-f","--force",help="Force execution of all tasks",action='store_true')
parser.add_argument("-r","--revise",help="Revise done jobs. If it's used, nothing is executed.",action='store_true')
parser.add_argument("-F","--do",help="Force execution of special tasks (separated by a comma)")
parser.add_argument("-x","--only",help="Only executing special tasks (separated by a comma) (priority is higher than that of -F and -f)")
parser.add_argument("-d","--clear",help="Clear lock file and exit",action="store_true")

parser.add_argument("-n", help="CPU number", type=int)
pargs = parser.parse_args()

BLACK = '\033[30m'
RED = '\033[31m'
GREEN = '\033[32m'
YELLOW = '\033[33m'
BLUE = '\033[34m'
MAGENTA = '\033[35m'
CYAN = '\033[36m'
WHITE = '\033[37m'
UNDERLINE = '\033[4m'
RESET = '\033[0m'
clearline = '\x1b[1A\x1b[2K'

lock = mp.Lock()

def colorprint(txt,color):
	print(color + txt + RESET)


def append_to_lockfile(txt):
	with open(config['lockfile'],"a") as f:
		f.write(txt + "\n")

def read_lockfile():
	if pargs.force:
		os.remove(config['lockfile'])
	if not os.path.isfile(config['lockfile']):
		os.system("touch " + config['lockfile'])
		return []
	print("Reading " + config['lockfile'])
	with open(config['lockfile']) as f:
		flines = set(f.read().strip().split("\n"))
	try:
		flines.remove("")
	except:
		pass
	return sorted(flines)

def unconditional_execution(job,*parameters):
	job(*parameters)
	append_to_lockfile(job.__name__)


def conditional_execution(job,*parameters):
	executed = False
	if len(only_tasks) > 0:
		if job.__name__ in only_tasks:
			unconditional_execution(job,*parameters)
			executed = True
		else:
			msg("Skipping " + job.__name__ + " as it is not required to perform",CYAN)
	else:
		if job.__name__ not in donejobs or job.__name__ in forced_tasks:
			unconditional_execution(job,*parameters)
			executed = True
			#except:
			#	msg("Error when doing task " + condition)
		else:
			msg("Skipping "+ job.__name__+" as it is already done",CYAN)
	print()
	return executed
	
def remove_sure(f):
	'''Removes file without yielding error message in case of inexisting file'''
	if os.path.isfile(f):
		os.remove(f)




def msg(cmd,color=""):
	'''Prints to screen'''
	message = time.strftime("%H:%M:%S ") + cmd
	if color == "":
		print(message)
	else:
		colorprint(message,color)
	logging.info(cmd)


def sysexec(cmd): 
	msg(cmd)
	os.system(cmd)

def read_cmd(cmd):
	return subprocess.check_output(cmd,shell=True,universal_newlines=True).strip()


def execute(procs):
	'''Execute processes in list'''
	for p in procs:
		p.start()
	for p in procs:
		p.join()





def merge_dic(diclist):
	'''Merges dictionaries of a list'''
	dic = {}
	for d in diclist:
		for key in d:
			dic[key] = d[key]
	return dic

def emptyfile(f,header=""):
	if header == "":
		if os.path.isfile(f):
			os.remove(f)
		cmd = "touch " + f
		msg(cmd)
		os.system(cmd)
	else:
		with open(f,"w") as out:
			out.write("\t".join(list(header))+"\n")

def writerows(rows,name,header = False):
	'''Kiiratas csv tablazatba'''
	msg("Printing to " + name,MAGENTA)
	if header:
		rows = [header] + rows
	with open(name,"w") as c:
		wtr = csv.writer(c,delimiter="\t",lineterminator="\n")
		wtr.writerows(rows)

def default_configs(config):
	if "project" not in config:
		config['project'] = 'myproject'
	if 'resdir' not in config:
		config['resdir'] = config['project'] + '_results'
	if 'tmpdir' not in config:
		config['tmpdir'] = config['project'] + '_tmp'
	if 'log' not in config:
		config['log'] = config['project'] + '_log'
	if 'lockfile' not in config:
		config['lockfile'] = config['tmpdir'] + "/" + config['project'] + '.donejobs'
	return config

def surely_make_directory(mydir):
	if not os.path.isdir(mydir):
		Path(mydir).mkdir(parents=True, exist_ok=True)
		msg("Directory " + mydir + " is created",GREEN)


def make_directories(szotar):
	for k,v in szotar.items():
		if type(v) == dict:
			make_directories(v)
		elif type(v) == str and not "http:" in v and not "https:" in v and "/" in v:
			surely_make_directory(v.rsplit("/",1)[0])

def replace_wildcard_in_dic(dic,old,new):
	outdict = {}
	for k,v in dic.items():
		if type(v) == dict:
			outdict[k] = replace_wildcard_in_dic(v,old,new) #Itt a rekurzio
		elif type(v) == str and old in v:
			outdict[k] = v.replace(old,new)
		else:
			outdict[k] = v
	return outdict




def parse_config():
	colorprint("\nParsing " + pargs.configfile,YELLOW+UNDERLINE)
	config = yaml.safe_load(open(pargs.configfile))
	config = default_configs(config)
	config = replace_wildcard_in_dic(config,"{project}",config["project"])
	config = replace_wildcard_in_dic(config,"{tmpdir}",config["tmpdir"])
	config = replace_wildcard_in_dic(config,"{resdir}",config["resdir"])
	config = replace_wildcard_in_dic(config,"{datadir}",config["datadir"])
	if 'resdir' in config:
		surely_make_directory(config['resdir'])
	if 'tmpdir' in config:
		surely_make_directory(config['tmpdir'])
	if 'datadir' in config:
		surely_make_directory(config['datadir'])
	make_directories(config)
	return config

def perform_cmd_if_missing_outfile(cmd,outlist):
	do = False
	#if pargs.overwrite:
	#	do = True
	#else:
	if type(outlist) == list:
		for out in outlist:
			if not os.path.isfile(out):
				do = True
				break
	elif not os.path.isfile(outlist):
		do = True
	if do:
		os.system(cmd)
	else:
		msg(out + " exists, skipping this step.")# For overwriting existing files, please use -overwrite when running the script")

def get_tasks_to_do():
	only_tasks = []
	forced_tasks = []
	if pargs.only:
		only_tasks = set(map(lambda x:x.strip(),pargs.only.split(",")))
	elif pargs.do:
		forced_tasks = set(map(lambda x:x.strip(),pargs.do.split(",")))
		print (forced_tasks)
	if pargs.revise:
		message = "\nReady jobs:\n" + "\n".join(donejobs) +"\n\n"
		if pargs.only:
			message += "Performing only tasks:\n" + "\n".join(only_tasks) + "\n\n"
		if pargs.do:
			message += "Forced jobs:\n" + "\n".join(forced_tasks) + "\n\n"
		print(message)
		os.abort()
	elif pargs.clear:
		os.remove(config['lockfile'])
		print("Removing " + config['lockfile'])
		os.abort()
	return[forced_tasks,only_tasks]
	
def initlog():
	logging.basicConfig(format='%(asctime)s. %(levelname)s:%(message)s', filename=config['log'] + "_" + time.strftime("%Y%m%d%H%M%S") + ".log" , level=logging.INFO) #Ez ele nem kerulhet logging (msg/sysexec) parancs!!!
	logging.Formatter(fmt='%(asctime)s',datefmt="%Y")
	print("Logging to: " + config['log'])
	msg("Config file: " + pargs.configfile)
	msg("Project name: " + config['project'])
	msg('Using {} threads'.format(threads))
	emptyfile(config['versions'])


def check_barcodes():
	msg("Correcting barcodes with revcomps")
	bcdict = defaultdict(lambda:"")
	with open(config['data']['barcodes'],"r") as f:
		rdr = csv.reader(f,delimiter="\t")
		next(rdr)
		for row in rdr:
			name = row[0]
			bc1 = row[1]
			bc2 = row[2]
			rcbc1 = str(Seq.Seq(bc1).reverse_complement())
			rcbc2 = str(Seq.Seq(bc2).reverse_complement())
			for i in [bc1,rcbc1]:
				for j in [bc2,rcbc2]:
					bcdict[i+'+'+j] = name
					bcdict[j+'+'+i] = name
	#with open(config['project'] + "_allpairs","w") as g:	
	#	for elem in bcdict:
	#		g.write(elem+"\t"+ bcdict[elem] +"\n")
	msg("Checking barcodes")
	barcodes = Counter(read_cmd("bioawk -c fastx '{print $comment}' " + config['data']['read1'] + " | cut -d':' -f 4").split())
	outrows = []
	inforows = []
	for bc in sorted(barcodes,key=lambda x:-barcodes[x]):
		if bc in bcdict:
			[bc1,bc2] = bc.split("+")
			outrows.append([bcdict[bc],bc1,bc2])
			inforows.append([bc,bcdict[bc],barcodes[bc]])
	outrows.sort(key = lambda x:x[0])
	inforows.sort(key = lambda x:x[-1],reverse=True)
	writerows(outrows,config['data']['corrbarcodes'].rsplit("/",1)[-1], ["Sample","Bc1","Bc2"])
	writerows(inforows,config['tmpdir'] + "/" + config['project']+"_bc_stats.tsv",['Barcode','Sample','BC_count'])

def parse_fastq_and_gz(fq,gz="?"):
	if gz=="?":
		if fq.endswith("gz"):
			gz = True
		else:
			gz = False
	if gz:
		return SeqIO.parse(gzip.open(fq,"rt"),"fastq")
	else:
		return SeqIO.parse(open(fq,"r"),"fastq")

def demux_fq(fq,r1,bcdict):
		if r1:
			suffix = "_R1.fastq"
		else:
			suffix = "_R2.fastq"
		demuxed_files = {}
		for k,v in bcdict.items():
			outf = config['datadir'] + "/" + v + suffix
			emptyfile(outf)
			demuxed_files[k] = outf
		readszam = 0
		msg("Demuxing " + fq)
		demuxed = {}
		filestreams = defaultdict(list)
		for rcd in parse_fastq_and_gz(fq):
			readszam += 1
			bc = rcd.description.split(':')[-1]
			if bc in bcdict:
				rcd.id = rcd.description.split("_")[0].split(" ")[0]
				rcd.name = rcd.id
				rcd.description = rcd.id
				filestreams[bc].append(rcd.format("fastq"))
			if readszam%1000000 == 0: #Lenullazzuk, hogy ne fogyasszon olyan sok memoriat
				msg("Processed " + str(readszam) + " reads from " + fq)
				for key in filestreams.keys():
					with open(demuxed_files[key],"a") as g:
						g.writelines(filestreams[key])
				filestreams = defaultdict(list)
		for key in filestreams.keys():
			with open(demuxed_files[key],"a") as g:
				g.writelines(filestreams[key])
		msg("Finished with demuxing " + fq)

def demux():
	msg("Reading barcodes")
	bcdict = {}
	with open(config['project'] + "_" + config['data']['barcodes'].rsplit("/",1)[-1],"r") as f:
		rdr = csv.reader(f,delimiter="\t")
		next(rdr)
		for row in rdr:
			bcdict["{}+{}".format(row[1],row[2])] = row[0]
	demux_fq(config['data']['read1'],True,bcdict)
	demux_fq(config['data']['read2'],False,bcdict)

def refmap():
	#with bwasema:
	pass

def get_samples():
	samples = read_cmd("cut -f 1 {} | tail -n +2".format(config['data']['corrbarcodes'])).split()
	return sorted(set(samples))

def get_versions():
	#os.system('bwa 2>&1 | grep Version | cut -f2 -d " " >> {}'.format(config['versions']))
	os.system('cutadapt --version >> {}'.format(config['versions']))

def get_readnum(sample):
	return int(read_cmd("wc -l {}/{}_R1.fastq".format(config['datadir'],sample)).split()[0])//4




class Amplicon:
	def __init__(self,name,r1start,r1end,r2start,r2end,bclength,*args):
		self.name = name
		self.r1start = r1start
		self.r1end = r1end
		self.r2start = r2start
		self.r2end = r2end
		self.bclength = bclength

def read_adapters():
	with open(config['data']['amp']) as f:
		rdr = csv.reader(f,delimiter="\t")
		next(rdr)
		return [Amplicon(*row) for row in rdr]

def cutadapt(sample,start,end,name,bclength,read):
	if start == "":
		adapter = end
	else:
		adapter = "^" + start + "..." + end
	sysexec("cutadapt --quiet -a {a} -o {tmp}/{s}.{amp}.{r}.fastq {datadir}/{s}_{r}.fastq --trimmed-only -e 0 -m {bc} -M {bc}".format(a=adapter,tmp=config['tmpdir'],s=sample,amp=name,datadir=config['datadir'],bc=bclength,r=read))
	sysexec("seqtk comp {tmp}/{s}.{amp}.{r}.fastq | cut -f 1 | sort > {tmp}/{s}.{amp}.{r}.list".format(tmp=config['tmpdir'],s=sample,amp=name,r=read))


def check_amplicon(sample,amplicon):
	msg("Cutting adapters for {} {}".format(sample, amplicon.name))
	cutadapt(sample,amplicon.r1start,amplicon.r1end,amplicon.name,amplicon.bclength,"R1")
	cutadapt(sample,amplicon.r2start,amplicon.r2end,amplicon.name,amplicon.bclength,"R2")
	sysexec("comm -12 {tmp}/{s}.{amp}.R1.list {tmp}/{s}.{amp}.R2.list > {tmp}/{s}.{amp}.paired.list".format(tmp=config['tmpdir'],s=sample,amp=amplicon.name))
	if read_cmd("wc -l {tmp}/{s}.{amp}.paired.list | cut -d ' ' -f 1".format(tmp=config['tmpdir'],s=sample,amp=amplicon.name)) == "0":
		totalreadnum = 0
		goodreadnum = 0
	else:
		msg("Selecting matching read pairs for {} {}".format(sample,amplicon.name))
		r1 = {row.split()[0]:row.split()[1] for row in read_cmd("seqtk subseq -t {tmp}/{s}.{amp}.R1.fastq {tmp}/{s}.{amp}.paired.list | cut -f 1,3".format(tmp=config['tmpdir'],s=sample,amp=amplicon.name)).split("\n")}
		r2 = {row.split()[0]:row.split()[1] for row in read_cmd("seqtk subseq -t {tmp}/{s}.{amp}.R2.fastq {tmp}/{s}.{amp}.paired.list | cut -f 1,3".format(tmp=config['tmpdir'],s=sample,amp=amplicon.name)).split("\n")}
		accepted_reads = defaultdict(list)
		goodreadnum = 0
		for read in r1.keys():
			if r1[read] == str(Seq.Seq(r2[read]).reverse_complement()):
				accepted_reads[r1[read]].append(read)
				goodreadnum += 1
		barcodes = sorted(list(accepted_reads.keys()))
		bcout = ""
		with open("{}/{}.{}.goodreads.csv".format(config['tmpdir'],sample,amplicon.name),'w') as f:
			for bc in barcodes:
				bcout += ("{}\t{}\t{}\t{}\n".format(sample,amplicon.name,bc,len(accepted_reads[bc])))
				for read in accepted_reads[bc]:
					f.write("{}\t{}\n".format(read,bc))
		with lock:
			with open(config['result']['bcnumbers'],"a") as g:
				g.write(bcout)
		totalreadnum = int(read_cmd("wc -l {}/{}_R1.fastq".format(config['datadir'],sample)).split()[0])//4
	with lock:
		with open(config['result']['readnums'],'a') as g:
			g.write("{}\t{}\t{}\t{}\n".format(sample,amplicon.name,totalreadnum,goodreadnum))
		
		
	#&writerows(list(zip(samples,readnums)),config['data']['readnums'],['Sample','Read number'])
		
if __name__=="__main__":
	config = parse_config()
	donejobs = read_lockfile()
	if pargs.n:
		threads = pargs.n
	elif 'threads' in config:
		threads = config['threads']
	else:
		threads = 1
	forced_tasks,only_tasks = get_tasks_to_do()
	sema = mp.Semaphore(threads)
	bwasema = mp.Semaphore(config['bwathreads'])
	lock = mp.Lock()
	multilock = [mp.Lock() for i in range(threads)]
	initlog()
	#Print versions to logfile
	##get_versions()
	#Demultiplex
	conditional_execution(check_barcodes)
	conditional_execution(demux)
	samples = get_samples()
	amplicons = read_adapters()
	emptyfile(config['result']['bcnumbers'],['Sample','Up/Down','Barcode','Readnumber'])
	emptyfile(config['result']['readnums'],['Sample','Up/Down','Total reads in sample','Matching reads'])
	#for sample in samples:
	#	for amplicon in amplicons:
	#		check_amplicon(sample,amplicon)
	pool = mp.Pool(threads)
	pool.starmap(check_amplicon,[(sample,amplicon,) for sample in samples for amplicon in amplicons])
	#Reference mapping
	#First indexing reference
	##perform_cmd_if_missing_outfile("bwa index {}".format(config['data']['amp1']),[config['data']['amp1'] + '.' + x for x in ["amb","ann","bwt","pac","sa"]])
	##perform_cmd_if_missing_outfile("bwa index {}".format(config['data']['amp2']),[config['data']['amp2'] + '.' + x for x in ["amb","ann","bwt","pac","sa"]])
	##hatch_process(refmap,samples)
