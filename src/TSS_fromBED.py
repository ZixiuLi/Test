#! /usr/bin/python
import sys
import re

BED_file=sys.argv[1]

fi=open(BED_file,"r")
for line in fi:
        line=line.strip()
        element=line.split("\t")
        Chr=element[0]
        Startpos=element[1]
        Endpos=element[2]
        Transcript=element[3]
        Strand=element[5]

        if re.match("^chr",Chr): 
                if Strand == "+":
                        start = int(Startpos)    # Startpos: 0 based
                        end = int(Startpos)+1 
                elif Strand == "-":
                        start = int(Endpos)-1    # Endpos: 1 based
                        end = int(Endpos)

                print Chr+"\t"+str(start)+"\t"+str(end)+"\t"+Transcript+"\t"+".\t"+Strand





