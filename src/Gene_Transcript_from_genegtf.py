#! /usr/bin/python
import re
import sys

GTF_file=sys.argv[1]
Transcript_Gene_Dict={}

fi=open(GTF_file,"r")
for line in fi:
        line=line.strip()
        if re.match("^#",line):
                next
        else:
                element=line.split("\t")
                try:
                        GeneTranscript_Infor=element[8]
                        gene=""
                        transcript=""
                        for item in GeneTranscript_Infor.split('; '):
                                if re.search("^gene_id",item):
                                        gene=item.split("\"")[1]
                                elif re.search("^transcript_id",item):
                                        transcript=item.split("\"")[1]
                        Transcript_Gene_Dict[transcript]=gene
                except Exception:
                        pass
fi.close

for T in Transcript_Gene_Dict:
        G=Transcript_Gene_Dict[T]
        print G+"\t"+T
