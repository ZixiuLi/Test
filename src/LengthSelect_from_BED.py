#! /bin/usr/python
import sys

BED_file=sys.argv[1]
cutoff=sys.argv[2]

fi=open(BED_file,"r")
for line in fi:
        line=line.strip()
        element=line.split("\t")
        transcript=element[3]
        blocksize_infor=element[10]

        blocksize_list=blocksize_infor.split(",")

        transcript_length=0
        for i in range(0,len(blocksize_list)-1):
                transcript_length = transcript_length + int(blocksize_list[i])

        if transcript_length >= int(cutoff):
                print line
