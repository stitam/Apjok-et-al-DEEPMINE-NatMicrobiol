#!/usr/bin/python3

###from multiprocessing import set_start_method
###set_start_method("spawn")
###from multiprocessing import get_context


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

lock = mp.Lock()

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


def colorprint(txt,color):
	print(color + txt + RESET)



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

def initlog():
	logging.basicConfig(format='%(asctime)s. %(levelname)s:%(message)s', filename=config['log'] + "_" + time.strftime("%Y%m%d%H%M%S") + ".log" , level=logging.INFO) #Ez ele nem kerulhet logging (msg/sysexec) parancs!!!
	logging.Formatter(fmt='%(asctime)s',datefmt="%Y")
	print("Logging to: " + config['log'])
	msg("Config file: " + pargs.configfile)
	msg("Project name: " + config['project'])
	msg('Using {} threads'.format(threads))
	emptyfile(config['versions'])

def get_versions():
	os.system('cutadapt --version > {}'.format(config['versions']))
	os.system(' seqtk 2>&1 | head -n 3 | tail -n 1 >> {}'.format(config['versions']))
	

class Amplicon:
	def __init__(self,name,r1start,r1end,r2start,r2end,bclength,*args):
		self.name = name
		self.r1start = r1start
		self.r1end = r1end
		self.r2start = r2start
		self.r2end = r2end
		self.bclength = bclength

def read_adapters():
	amplicons = {}
	with open(config['amp']) as f:
		rdr = csv.reader(f,delimiter="\t")
		next(rdr)
		for row in rdr:
			amplicons[row[0]] = Amplicon(*row)
	return amplicons

def cutadapt(sample,start,end,name,bclength,read):
	adapter = "X" + start + "..." + end + "X"
	sysexec("cutadapt --quiet -g {a} -o {tmp}/{s}.{amp}.{r}.fastq {datadir}/{s}.{r}.fastq --trimmed-only -e 0".format(a=adapter,tmp=config['tmpdir'],s=sample,amp=name,datadir=config['datadir'],bc=bclength,r=read))
	sysexec("seqtk comp {tmp}/{s}.{amp}.{r}.fastq | cut -f 1 | sort > {tmp}/{s}.{amp}.{r}.list".format(tmp=config['tmpdir'],s=sample,amp=name,r=read))


def get_barcodes_from_fq(fq_name,amp):
	msg("Reading " + fq_name)
	fq = SeqIO.parse(fq_name,"fastq")
	bc_dict = {}
	for rcd in fq:
		bc,mid = ["",""]
		if amp not in rcd.seq:
			rcd.seq = rcd.seq.reverse_complement()
		if amp in rcd.seq:
			ix = rcd.seq.index(amp)
			bc = rcd.seq[:ix][-8:]
			if len(bc) == 8:
				mid = rcd.seq[ix + len(amp):ix + len(amp)+10] #a kozepso barcode 10 hosszu
			else:
				bc = ""
		bc_dict[rcd.id] = [str(bc),str(mid)]
	return bc_dict



###	mynames = list(prevalence.keys())
###	mynames.sort()
###	msg("Writing out results for " +  sample)
###	with open(config['result']['bcnumbers'],"a") as g:
###		for key in mynames:
###			##g.write("\t".join([sample,amplicon.name,key,str(prevalence[key]),"{:.02f}".format(prevalence[key]/totalreadnum*100)]) + "\n")
###			g.write("\t".join([sample,amplicon.name,key,str(prevalence[key])]) + "\n")
###			good_readnum += prevalence[key]
		
def check_amplicon(sample,amplicon,names_for_bcpairs):
	####if (sample.startswith("UP") and amplicon.name == "Downstream") or (sample.startswith("DOWN") and amplicon.name == "Upstream"):
	####	return
	msg("Processing " + sample)
	##sysexec("cutadapt -a {0} -A {1} --action=none --trimmed-only --discard-untrimmed -o {2}/{3}.amplicon.R1.fastq -p {2}/{3}.amplicon.R2.fastq {2}/{3}.trim.R1.fastq {2}/{3}.trim.R2.fastq ".format(amplicon.r1start,amplicon.r2start,config['datadir'],sample))
	bc_dict_1 = get_barcodes_from_fq("{}/{}.trim.R1.fastq".format(config['datadir'],sample),amplicon.r1start)
	bc_dict_2 = get_barcodes_from_fq("{}/{}.trim.R2.fastq".format(config['datadir'],sample),amplicon.r2start)
	bc_mismatch = 0
	good_reads = 0
	no_name = 0
	prevalence = defaultdict(lambda:0)
	msg("Checking barcodes for " + sample)
	good_read_num = 0
	for key in bc_dict_1:
		bc1,mid1 = bc_dict_1[key]
		bc2,mid2 = bc_dict_2[key]
		if bc1 and bc2 and mid1 and mid2 and (mid1 == mid2 or mid1 == revcomp(mid2)) and (sample,bc1,bc2,) in names_for_bcpairs:
			good_read_num += 1
			prevalence[(names_for_bcpairs[(sample,bc1,bc2,)],sample,bc1,bc2,mid1,)] += 1
	totalreadnum = len(bc_dict_1.keys())
	msg("Collecting results for " +  sample)
	bcout = ""
	prevalence_keys = list(prevalence.keys())
	prevalence_keys.sort(key = lambda x: (x[1],x[0],-prevalence[x],))
	for key in prevalence_keys:
		bcout += "\t".join(key) + "\t" +amplicon.name + "\t" + str(prevalence[key]) + "\n"
		##if key[1] == key[4] and key[2] == key[5] and key[3] == key[6] and len(key[1]) == 8 and len(key[2]) == 8:
		##	accepted_read_num += prevalence[key]
		##	if (sample,key[1],key[2],) in names_for_bcpairs:
		##		bcout += names_for_bcpairs[(sample,key[1],key[2],)]
		##		good_reads += "\t".join([sample,names_for_bcpairs[(sample,key[1],key[2],)],key[1],key[2],key[3],str(prevalence[key])])+"\n"
		##		good_read_num += prevalence[key]
		##	else:
		##		bcout += "-"
		##else:
		##	bcout += "0\t-"
		##bcout += "\n"
	with lock:
		msg("Writing out results for " +  sample)
		with open(config['result']['bcnumbers'],"a") as g:
			g.write(bcout)
		#with open(config['result']['good_reads'],'a') as g:
		#	g.write(good_reads)
		with open(config['result']['readnums'],'a') as g:
			g.write("\t".join(map(str,[sample,amplicon.name,totalreadnum,good_read_num]))+"\n")


def concat_lanes(sample):
	ls1 = sorted(glob.glob(config['datadir'] + "/" + sample + "*R1_001.fastq.gz"))
	ls2 = sorted(glob.glob(config['datadir'] + "/" + sample + "*R2_001.fastq.gz"))
	emptyfile("{}/{}.R1.fastq.gz".format(config['datadir'],sample))
	emptyfile("{}/{}.R2.fastq.gz".format(config['datadir'],sample))
	for item in ls1:
		os.system("cat {} >> {}/{}.R1.fastq.gz".format(item,config['datadir'],sample))
	for item in ls2:
		os.system("cat {} >> {}/{}.R2.fastq.gz".format(item,config['datadir'],sample))
	os.system("gunzip {}/{}.R1.fastq.gz -f -q".format(config['datadir'],sample))
	os.system("gunzip {}/{}.R2.fastq.gz -f -q".format(config['datadir'],sample))

def sort_and_add_header(f,header):
	with open(f + ".tmp","w") as g:
		g.write("\t".join(map(str,header)) + "\n")
	cmd = "sort {0} >> {0}.tmp".format(f)
	sysexec(cmd)
	os.remove(f)
	os.rename(f + ".tmp",f)

def revcomp(s):
	return str(Seq.Seq(s).reverse_complement())

def read_bc_details():
	with open(config['barcodes']) as f:
		rdr = csv.reader(f,delimiter="\t")
		next(rdr)
		names_for_bcpairs = {}
		for row in rdr:
			if len(row) >= 4:
				names_for_bcpairs[(row[0],row[1],row[2],)] = row[3]
				names_for_bcpairs[(row[0],revcomp(row[1]),row[2],)] = row[3]
				names_for_bcpairs[(row[0],row[1],revcomp(row[2]),)] = row[3]
				names_for_bcpairs[(row[0],revcomp(row[1]),revcomp(row[2]),)] = row[3]
	return names_for_bcpairs

def trim_reads(sample):
	sysexec("cutadapt --quiet -q 30 -o {0}/{1}.trim.R1.fastq -p {0}/{1}.trim.R2.fastq {0}/{1}.R1.fastq {0}/{1}.R2.fastq".format(config['datadir'],sample))




if __name__=="__main__":
	config = parse_config() # -->  config = {'threads': 60, 'samples': 'samples.csv',...}
	if pargs.n:
		threads = pargs.n
	elif 'threads' in config:
		threads = config['threads']
	else:
		threads = 1
	initlog()
	#Print versions to logfile
	get_versions()
	#Demultiplex
	pool = mp.Pool(config['threads'])
	samples = open(config['samples'],"r").read().strip().split()
	#!pool.map(concat_lanes,samples)
	pool.map(trim_reads,samples)
	amplicons = read_adapters()
	emptyfile(config['result']['bcnumbers'],["Sample","Batch","BC1","BC2","Mid BC","UP/DOWN","Readnumber"])
	emptyfile(config['result']['readnums'],["Batch name","UP/DOWN","Total number of reads","Good reads with given barcodes"])
	names_for_bcpairs = read_bc_details()
	check_amplicon("Uptag",amplicons["Upstream"],names_for_bcpairs)
	check_amplicon("Downtag",amplicons["Downstream"],names_for_bcpairs)
	##check_amplicon("Uptag",amplicons["Downstream"],names_for_bcpairs)
	##check_amplicon("Downtag",amplicons["Upstream"],names_for_bcpairs)

