#! /usr/bin/python
import sys
import re
import json

with open('config.json') as input_config:
        config = json.load(input_config)

chrNameLength_File = config['file_config']['chr_length']

BED_file = sys.argv[1]
Extend_length = int(sys.argv[2])

Chr_Length_Dict={}
fj=open(chrNameLength_File,"r")
for line in fj:
        line=line.strip()
        element=line.split("\t")
        Chr=element[0]
        Length=int(element[1])
        Chr_Length_Dict[Chr]=Length
fj.close()

fi=open(BED_file,"r")
for line in fi:
        line=line.strip()
        element=line.split("\t")
        Chr=element[0]
        Startpos=element[1]
        Endpos=element[2]
        Transcript=element[3]
        Strand=element[5]

        Chr_length=Chr_Length_Dict[Chr]
        if re.match("^chr",Chr):
                if Strand == "+":
                        if int(Startpos)-Extend_length < 0:
                                start=0
                        else:
                                start=int(Startpos)-Extend_length
                        if int(Startpos)+Extend_length > Chr_length:
                                end=Chr_length
                        else:
                                end=int(Startpos)+Extend_length
                elif Strand == "-":
                        if int(Endpos)+Extend_length > Chr_length:
                                end=Chr_length
                        else:
                                end=int(Endpos)+Extend_length
                        if int(Endpos)-Extend_length < 0:
                                start=0
                        else:
                                start=int(Endpos)-Extend_length
                        

                print Chr+"\t"+str(start)+"\t"+str(end)+"\t"+Transcript+"\t"+".\t"+Strand





