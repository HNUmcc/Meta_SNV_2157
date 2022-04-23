#!/usr/bin/python
# -*- coding:utf-8 -*-

import os,sys
#import linecache
listfile=open(sys.argv[1],'r')
line = listfile.readline()
#num=len(line.strip().split('\t'))-2
num=len(line.strip().split('\t'))-1
#listfile.seek(0) #zhi zhen chong xin hui dao wen jian kai tou 
#num=len(linecache.getline('listfile',1).strip().split('\t'))-2
leval=float(num*0.5)
sum=0
outfile=open(sys.argv[2],'w')
outfile.write(line)
for line in listfile.readlines():
    sum=0
    abc=line.strip().split('\t')[1:]
#   abc=line.strip().split('\t')[2:]
    for i in range(0,num):
        if float(abc[i])>0:
           sum+=1
        else:
           pass
    if sum>=leval:
        outfile.write(line)
    else:
        pass
listfile.close()
outfile.close()
