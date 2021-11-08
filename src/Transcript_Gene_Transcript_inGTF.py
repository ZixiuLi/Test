#! /usr/bin/python
import sys
import re
from collections import defaultdict

coding_transcript_input=sys.argv[1]
gtf_input=sys.argv[2]

Coding_T_Dict={}
fj=open(coding_transcript_input,"r")
for line in fj:
        line=line.strip()
        T_ID=line.split("\t")[0]
        Coding_T_Dict[T_ID]=1
fj.close()

Coding_G_Dict={}
G_T_Dict=defaultdict(dict)
fi=open(gtf_input,"r")
for line in fi:
        line=line.strip()
        if re.match("^#",line):
                pass
        else:
                try:
                        element=line.split("\t")
                        Type=element[2]
                        Information=element[8]
                
                        if Type == "transcript":
                                Information_List=Information.split(";")
                                Gene=Information_List[0].replace("gene_id ","").replace("\"","").replace(" ","")
                                Transcript=Information_List[1].replace("transcript_id ","").replace("\"","").replace(" ","")
                        
                                G_T_Dict[Gene][Transcript]=1
                        
                                if Transcript in Coding_T_Dict:
                                        Coding_G_Dict[Gene]=1
                                else:
                                        pass
                except Exception:
                        pass
fi.close()

for gene in Coding_G_Dict:
        for transcript in G_T_Dict[gene]:
                print transcript
