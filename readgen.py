#!/usr/bin/env python
import sys
import random 
import os
import argparse
from collections import defaultdict


def genReads (iter1, length, dna, readstring, stposval, endposval,totlength):

	L1=length-int(0.05*length)
	L2=length+int(0.05*length)
	#print L1,L2
	for i in range(totiter,totiter+iter1):

		length=random.randint(L1,L2)
		lengthlist.append(length)
		pos=random.randint(stposval,(endposval-length))
		for j in range(length):
			hashposlength[totlength+j]=length
			hashstpos[totlength+j]=totlength
			hashgenpos[totlength+j]=pos
			hashrounditer[totlength+j]=i
			#print i,totlength+j,totlength,pos+j
		hashstpos1[i]=totlength
		read=dna[pos:pos+length]
		position.append(pos)		
		readstring+=read
		totlength+=length
	return readstring, totlength

def errmutate(dna, delidcount):
	i=0
	nucleotides=['A','C','G','T']
	while i<numerror:
		i+=1
		pos=random.randint(0,len(dna)-1)
		readstpos=hashstpos[pos]
		genpos=hashgenpos[pos]
		ai=hashrounditer[pos]
		rem=pos-readstpos
		if genpos in varpos1:
			i-=1
			continue		
		newbase=random.choice(nucleotides)
		a=dna[pos:pos+1]
		if a==newbase:
			nucs = [x for x in nucleotides if x != a]
			newbase=random.choice(nucs)
		dna=dna[:pos]+newbase+dna[pos+1:]
	return dna		

def errdeletion(errdna, delidcount):
	errdelreadid=[]
	bases=['A','C','G','T']
	i=0
	while i<numerrdel:		
		pos = random.randint(0,len(errdna)-1)
		length=hashposlength[pos]
		readstpos=hashstpos[pos]
		genpos=hashgenpos[pos]
		ia=hashrounditer[pos]
		rem=pos-readstpos
		#print pos,length,ia,genpos,readstpos,rem
		endcoord=genpos+length
		readid=str(ia)+":"+str(genpos)+":"+str(endcoord)
		if genpos in varpos1:
			i-=1
			continue
		a=errdna[0:pos]
		b=errdna[readstpos+rem+1:readstpos+1+length]
		c=errdna[pos+length-rem:]
		oriread=errdna[readstpos:readstpos+length]
		parta=errdna[readstpos:pos]
		#print oriread,pos,length,ia,genpos,readstpos,rem 
		if len(parta)<2 or len(b)<2:
			i-=1
			continue
		try:
			delidcount[readid]
		except:
			delidcount[readid]=0

		delidcount[readid]+=1

		if delidcount[readid]>1:
			add=delidcount[readid]-1
			addnuc=dna[genpos+length+add:genpos+length+add+1]
		elif delidcount[readid]==1:
			add=1
			addnuc=dna[genpos+length:genpos+length+1]	
		else:
			add=delidcount[readid]-1
			addnuc=dna[genpos+length+add:genpos+length+add+1]	

		if addnuc is "":  ##deletion in near end of genome
			i-=1
			continue	

		newread=errdna[readstpos:readstpos+rem]+errdna[readstpos+rem+1:readstpos+length]+addnuc
		#print readid,pos,length,ia,genpos,readstpos,rem,len (parta),len (b), "\n",oriread,"\n",newread 
		errdna=errdna[:readstpos]+newread+errdna[readstpos+length:]
		ori1read=dna[genpos:genpos+length]

		delsite=genpos+rem
		genposend=genpos+length
		readid=readstpos / length
		errdelreadid.append(readid)
		i+=1
	return errdna

def errinsertion(errdna, delidcount):
	bases=['A','C','G','T']
	i=0
	while i<numerrins:
		i+=1
		nbases=[]
		n1bases=[]
		nucleotides=['A','C','G','T']
		pos = random.randint(0,len(errdna)-1)
		length=hashposlength[pos]
		readstpos=hashstpos[pos]
		genpos=hashgenpos[pos]
		ai=hashrounditer[pos]
		rem=pos-readstpos
		if genpos in varpos1:
			i-=1
			continue
		a=errdna[:pos]
		b=errdna[readstpos+rem:readstpos+length-1]
		c=errdna[pos+length-rem:]
		
		prvnuc=a[:1]
		a=a+prvnuc
		nextnuc=b[:1]
		oriread=errdna[readstpos:readstpos+length]
		parta=errdna[readstpos:pos]
		readid=str(ai)+":"+str(genpos)+":"+str(genpos+length)
		if nextnuc in nucleotides:
			nbases = [x for x in nucleotides if x != nextnuc]
		else:
			nbases=nucleotides
		if prvnuc in nbases:
			n1bases=[x for x in nucleotides if x != prvnuc]
		else:
			n1bases=nbases

		insertnuc=random.choice(n1bases)
		if (len (parta)<2 or len(b)<2):
			i-=1
			continue
		
		newread=errdna[readstpos:readstpos+rem]+insertnuc+errdna[readstpos+rem:readstpos+length-1]
		#print readid,pos,length,ai,genpos,readstpos,rem,len (parta),len (b), "\n",oriread,"\n",newread
		errdna=errdna[:readstpos]+newread+errdna[readstpos+length:]
		insite=genpos+rem
	return errdna

def mutate(varpos,readstring,varsize,varnewallele):
	allvarpos,snpreadpos=[],[]
	count=0
	bases=['A','C','G','T']

	stringmutnuc=""
	iter1=len(position)
	for i in range (iter1):
		length=lengthlist[i]
		if position[i]>=varpos-length+1 and position[i]<=varpos:
			endcoord=position[i]+length
			readid=str(i)+":"+str(position[i])+":"+str(endcoord)
			for jj in range (varsize):
				hashofarray[varpos+1+jj].append(readid)
			count+=1
		if position[i]>=varpos-length+1 and position[i]<=varpos-1:
			allvarpos.append(position[i])

	numvarread=int(0.01*varrate*count)
	if len(allvarpos)<numvarread:
		numvarread=len(allvarpos)
		#print "adjustment:","alt:",numvarread, "Total:", count,"len(varpos):",len(allvarpos)
	for i in range(numvarread):
		tarpos=random.choice(allvarpos)
		allvarpos.remove(tarpos)
		snpreadpos.append(tarpos)
	for i in range(iter1):
		if position[i] in snpreadpos:
			snpreadpos.remove(position[i])
			npos=hashstpos1[i]
			stpos=npos+(varpos-position[i])
			length=lengthlist[i]
			rem=varpos-position[i]
			endpos=npos+length
			jj=0
			if rem<varsize:
				inc=rem
			else:
				inc=varsize
			for jj in range(inc):
				currnucl=readstring[stpos+jj:stpos+jj+1]
				snpnucs=[x for x in bases if x != currnucl]
				if len(varnewallele)<1:
					snp=random.choice(snpnucs)
				else:
					snp=varnewallele[jj]
				readstring=readstring[:stpos+jj]+snp+readstring[stpos+jj+1:]
				varpos1=varpos+1+jj
				hashofarrayalt[varpos1].append(readid)	
				a=readstring[npos:stpos+jj-2]
				b=readstring[stpos+1+jj:npos+length+2]
				newread=a+snp+b
				oriread=dna[position[i]:position[i]+length]
				hashvarpos1[varpos1].append(snp)					

			endcoord=position[i]+length
			readid=str(i)+":"+str(position[i])+":"+str(endcoord)
	stringmutnuc=""
	for key in hashvarpos1:
		stringmutnuc=""
		for base in bases:
			stringmutnuc+=base+":"+str(hashvarpos1[key].count(base))+" "
		hashaltallelefreq[key].append(stringmutnuc)	
	return readstring

def deletion(varpos,readstring,varsize):
	allvarpos,delreadpos=[],[]
	count=0
	bases= ['A', 'C', 'G', 'T']
	iter1=len(position)

	for i in range(iter1):
		varpos1=varpos-1
		length=lengthlist[i]
		varpos1=varpos-1
		if position[i]>=(varpos-length) and position[i]<=varpos-1:		#if position[i]>=(varpos-length+1) and position[i]<=varpos:
			endcoord=position[i]+length;
			readid=str(i)+":"+str(position[i])+":"+str(endcoord)
			#print varpos,readid
			hashofarray[varpos].append(readid)
			count+=1
		if position[i]>=(varpos-length+4) and position[i]+4<=varpos:
			allvarpos.append(position[i])
	
	numvarread=int(0.01*varrate*count)
	if len(allvarpos)<numvarread:
		numvarread=len(allvarpos)
		#print "adjustment:","alt:",numvarread, "Total:", count,"len(varpos):",len(allvarpos)
	
	for i in range(numvarread):
		tarpos=random.choice(allvarpos)
		allvarpos.remove(tarpos)
		delreadpos.append(tarpos)
	i=0
	#print "iter1 is", iter1, len(hashstpos1)
	while i < iter1:		
		if position[i] in delreadpos:
			npos=hashstpos1[i]
			length=lengthlist[i]
			rem=varpos-position[i]
			delreadpos.remove(position[i])
			stpos=npos+(varpos-position[i])
			if rem<varsize:
				inc=rem
			else:
				inc=varsize

			a=readstring[npos:stpos]
			b=readstring[stpos+inc:npos+length]
			endcoord=position[i]+length
			readid=str(i)+":"+str(position[i])+":"+str(endcoord)
			
			if len(a)<1 or len(b)<1:
				i-=1
				continue

			try:
				delidcount[readid]
			except:
				delidcount[readid]=0

			delidcount[readid]+=inc
			if delidcount[readid]>1:
				add=delidcount[readid]-1
				nextnuc=dna[position[i]+length:position[i]+length+inc]
				#print readid,nextnuc,"1"
			elif delidcount[readid]==1:
				add = 1
				nextnuc=dna[position[i]+length:position[i]+length+inc]
				#print readid,nextnuc,"2"				
			else:
				add=delidcount[readid]-1
				nextnuc=dna[position[i]+length+add:position[i]+length+add+inc]
				#print readid,nextnuc,"3"
			newread=a+b+nextnuc
			readstring=readstring[:npos]+newread+readstring[npos+length:]
			oriread=dna[position[i]:position[i]+length]
			#print readid,"\n", oriread,"\n",newread
			hashofarrayalt[varpos].append(readid)
		i+=1
	return readstring


def insertion(varpos,readstring,varsize,varnewallele):
	allvarpos,insreadpos,nbases=[],[],[]
	count=0
	bases= ['A', 'C', 'G', 'T']
	inscount=0
	hashinsnuc={}

	iter1=len(position)
	for i in range(iter1):
		varpos1=varpos-1
		length=lengthlist[i]
		if position[i]>=(varpos-length) and position[i]<=varpos-1:  #if   #if position[i]-1>=(varpos1-length) and position[i]<=varpos1:
			endcoord=position[i]+length
			readid=str(i)+":"+str(position[i])+":"+str(endcoord)
			hashofarray[varpos].append(readid)
			#print varpos1,readid
			count+=1
		if position[i]>=(varpos-length+1) and position[i]+1<=varpos:
			allvarpos.append(position[i])

	numvarread=int(0.01*varrate*count)
	if len(allvarpos)<numvarread:
		numvarread=len(allvarpos)
		#print "adjustment:","alt:",numvarread, "Total:", count,"len(varpos):",len(allvarpos)
	
	for i in range(numvarread):
		tarpos=random.choice(allvarpos)
		allvarpos.remove(tarpos)
		insreadpos.append(tarpos)

	i=0

	while i < iter1:	
		if position[i] in insreadpos:
			insreadpos.remove(position[i])
			npos=hashstpos1[i]
			rem=varpos-position[i]
			stpos=npos+rem
			length=lengthlist[i]
			oriread=dna[position[i]:position[i]+length]
			
			if rem<varsize:
				inc=rem
			else:
				inc=varsize
			a=readstring[npos:stpos]
			b=readstring[stpos:(npos+length)-inc]  ## -varsize
			nextnuc=dna[varpos:varpos+1]
			prvnuc=a[-1:]
			if len(b)>1 and len(a)>1:
				nextnuc=b[:1]				
			else:
				i-=1
				continue

			nbases = [x for x in bases if x != prvnuc]
			n1bases= [x for x in nbases if x != nextnuc]

			insertnuc=""
			if len(varnewallele) < 1:
				for k in range(inc):
					insertnuc+=random.choice(n1bases)
			else:
				insertnuc=varnewallele[:inc]
				if inc > len(varnewallele):
					print "Error: Insert size larger than insert fragment"
					sys.exit()

			newread=a+insertnuc+b
			readstring=readstring[:npos]+newread+readstring[npos+length:]
			endcoord=position[i]+length
			readid=str(i)+":"+str(position[i])+":"+str(endcoord)

			try:
				delidcount[readid]
			except:
				delidcount[readid]=0

			delidcount[readid]-=inc

			try:
				hashinsnuc[insertnuc]+=1
			except:
				hashinsnuc[insertnuc]=1	

			hashofarrayalt[varpos].append(readid)
			inscount+=1
			#print readid,npos,stpos,len(a),len(b),"\n", oriread,"\n",newread
		i+=1

	refins=count-inscount
	stringinsnuc=""
	for elem in hashinsnuc:
		stringinsnuc+=elem+":"+str(hashinsnuc[elem])+" "
	hashaltallelefreq[varpos].append(stringinsnuc)
	return readstring

def get_arguments():
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("input", help="Reference genome file")
	parser.add_argument("var", help="Tab delimited file with variant position")
	parser.add_argument("output", help="Output filename")
	parser.add_argument("-len", help="length of reads",nargs='?', const=50, default=50, type=int)
	parser.add_argument("-cov", help="coverage", nargs='?', const=100, default=100, type=int)
	parser.add_argument("-ser", help="SNP error rate", nargs='?', const=0.5, default=0.5, type=float)
	parser.add_argument("-der", help="Deletion error rate", nargs='?', const=0.5, default=0.5, type=float)
	parser.add_argument("-ier", help="Insertion error rate", nargs='?', const=0.5, default=0.5, type=float)
	parser.add_argument("-thr", help="Number of threads", nargs='?', const=4, default=4, type=int)
	parser.add_argument("-hsr", help="comma separated list of hotspot regions", nargs='*', default="WholeGenome")
	args = parser.parse_args()
	return (args.input, args.var, args.output, args.len, args.cov, args.ser, args.der, args.ier, args.thr, args.hsr)

param=get_arguments()
infile=param[0]
bed=param[1]
outfile=param[2]
length=param[3]
coverage=param[4]
errorsnp=param[5]
errordel=param[6]
errorins=param[7]
threads=param[8]
hotspots=param[9]

#print infile, bed, outfile, length, coverage, errorsnp, errordel, errorins, threads, hotspots


hashposlength,hashgenpos,hashrounditer,hashstpos,hashstpos1={},{},{},{},{}
lengthlist=[]

totlength=0
#print "totlenght is", totlength
delidcount,tothash, althash,hashvar,hashout={},{},{},{},{}
hashofarray = defaultdict(list)
hashofarrayalt = defaultdict(list)
hashaltallelefreq = defaultdict(list)
hashvarpos1=defaultdict(list)
varposition,hmp,varpos1,stlist,endlist,unalgn,varposex=[],[],[],[],[],[],[]
trinucs=['AAA','CCC','GGG','TTT','RRR','YYY','SSS','WWW','KKK','MMM','BBB','DDD','HHH','VVV','NNN']
global position
position=[]
readstring=""
dna=""

with open (infile) as file1:
	for line in file1:
		line=line.rstrip()
		if '>' in line:
			continue
		else:
			dna+=line

with open (bed) as bedfile:
	for line in bedfile:
		line=line.rstrip()
		varpos=int(line.split('\t')[1])
		varpos1.append(varpos)
		hashvar[varpos]=line

if hotspots is "WholeGenome" or not hotspots:
	stlist.append(0)
	endlist.append(len(dna)-length)
else:
	for hotspot in hotspots:
		#print hotspot
		try:
			stval=int(hotspot.split(":")[0].strip())		
			endval=int(hotspot.split(":")[1].strip())
			#print stval,endval
			stlist.append(stval)
			endlist.append(endval)
		except:
			print "Incorrect format for hotspot region"
			sys.exit()
totiter=0	
for stposval, endposval in zip(stlist,endlist):
	i=stposval
	while i<endposval-3:
		sub=dna[i:i+3]
		j=1
		if sub in trinucs:
			j=i
			switch=0
			while switch==0:
				aa=dna[j:j+1]
				if aa in sub:
					#print j, sub, aa
					hmp.append(j)
					j+=1
				else:
					switch=1
					j=1
		i+=j

	#print "hmp is", hmp
	dna1=dna[stposval:endposval]
	iter1=int((coverage*len(dna1))/length)
	
	numerror=int(errorsnp*0.01*coverage*len(dna1))
	numerrdel=int(errordel*0.01*coverage*len(dna1))
	numerrins=int(errorins*0.01*coverage*len(dna1))	

	readstringret=genReads(iter1,length,dna,readstring, stposval, endposval,totlength) #generate error and variation free reads
	totiter+=iter1
	readstring=readstringret[0]
	totlength=readstringret[1]
	varpos1.sort(reverse=True)
	for varpos in varpos1:
		line=hashvar[varpos]
		vartype=(line.split('\t')[0])
		varrate=float(line.split('\t')[2])
		varsize=int(line.split('\t')[3])
		try:
			varnewallele=line.split('\t')[4]	
		except:
			varnewallele=""
		if varpos>=stposval and varpos<=endposval:
			#print varpos,stposval,endposval 
			if (vartype) == 'SNP':
				varpos=varpos-1
				readstring=mutate(varpos,readstring,varsize,varnewallele)
				
				varposition.append(varpos+1)
				continue
			elif vartype == 'Del':
				#print varpos
				if varpos in hmp:
					print "Warning: Position",varpos,"is in homopolymeric region"
				readstring=deletion(varpos,readstring,varsize)
				varposition.append(varpos)
				continue
			elif vartype == 'Ins':
				readstring=insertion(varpos,readstring,varsize,varnewallele)
				varposition.append(varpos)
			else:
				print "Unrecognized variant type. Variant type must be SNP, Ins or Del"

readstring=errmutate(readstring, delidcount) #introduce substitution errors
readstring=errdeletion(readstring, delidcount)
readstring=errinsertion(readstring, delidcount)
iter1=totiter
f1=open(outfile,"w")
npos=0
for i in range (iter1):
	try:
		length=lengthlist[i]
		qualstring='5' * length
		end=position[i]+length
		mutread=readstring[npos:npos+length]
		string1="@"+str(i)+":"+str(position[i])+":"+str(end)+"\n"+mutread+"\n"+"+"+"\n"+qualstring+"\n"
		f1.write(string1)
		posst=position[i]+1
		npos+=length
	except:
		pass
i=0
'''
##To remove unaligned reads
os.system("bwa mem -L 100 -t %s %s %s -v1 > %s" %(threads, infile, outfile, outfile+".sam"))

with open (outfile+".sam") as samfile:
	for line in samfile:
		sam=line.split("\t")
		if sam[2] is "*" and sam[5] is "*":
			id1=sam[0]
			unalgn.append(id1)
			#continue
'''	
for key in hashofarray:
	totcount=0
	for item in hashofarray[key]:
		if item in unalgn:
			continue 
		else:
			totcount+=1
	tothash[key]=totcount

for key in hashofarrayalt:
	altcount=0
	for item in hashofarrayalt[key]:
		if item in unalgn:
			continue 
		else:
			altcount+=1
	althash[key]=altcount
	
for key in hashofarray:
	varposex.append(key)
	var=key
	try:
		tot=tothash[var]
		alt=althash[var]
		ref=tot-alt
		if tot<1:
			continue
		pc=(alt/float(tot))*100
		if (ref+alt)<100:
			#print "Warning: At pos",var,"total reads is less than 100"
			pass
		hashout[var]=str(var)+"\t"+str(ref)+"\t"+str(alt)+"\t"+str(pc)
	except:
		continue

varposex.sort()
print "Pos\tNo. of Ref\tNo. of Variant\t%Varaiant\tAlleles"
for var in varposex:
	try:
		print hashout[var]+"\t"+hashaltallelefreq[var][0]
	except:
		print hashout[var]

f1.close()
bedfile.close()
file1.close()
