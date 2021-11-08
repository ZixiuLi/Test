#! /usr/bin/python
import sys

input_file=sys.argv[1]

fi=open(input_file,"r")
T=fi.readline().strip().split("(")[0].replace(">","")
fi.readline()
fi.readline()
Promoter=fi.readline().strip().replace(" promoter(s)  were predicted","")
fi.close()

print T+"\t"+Promoter
