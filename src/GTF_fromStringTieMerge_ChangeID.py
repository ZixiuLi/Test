#! /usr/bin/python
import sys
import re
from collections import OrderedDict

GeneID_Dict=OrderedDict()

gtf_file = sys.argv[1]
gtf_output = sys.argv[2]

fo = open(gtf_output, "w")

Index = int(0)
fi = open(gtf_file, "r")
for line in fi:
        line = line.strip()
        if re.match("^#", line):
                pass
        else:
                element = line.split("\t")
                try:    
                        Type = element[2]
                        if Type == "transcript":
                                GeneTranscript_Infor = element[8]
                                for item in GeneTranscript_Infor.split('; '):
                                        if re.match("^gene_id", item):
                                                gene_id=item.split("\"")[1]
                                
                                if re.match( "^MSTRG", gene_id):        
                                        if gene_id in GeneID_Dict:
                                                pass
                                        else:
                                                Index = Index + 1
                                                GeneID_Dict[gene_id] = "MSTRG." + str(Index)
                except:
                        pass
fi.close()

fi = open(gtf_file, "r")
for line in fi:
        line = line.strip()
        if re.match("^#", line):
                fo.write(line + "\n")
        else:
                element = line.split("\t")
                try:    
                        GeneTranscript_Infor = element[8].split('; ')
                        GeneID_ori = GeneTranscript_Infor[0].split("\"")[1]
                        TranscriptID_ori = GeneTranscript_Infor[1].split("\"")[1]
                        
                        if re.match("^MSTRG", GeneID_ori):
                                GeneID_change = GeneID_Dict[GeneID_ori]
                                GeneTranscript_Infor[0] = "gene_id \"" + GeneID_change + "\""
                                if re.match("^MSTRG", TranscriptID_ori):
                                        TranscriptID_change = GeneID_Dict[GeneID_ori] + "." + str(TranscriptID_ori.split(".")[2])
                                        GeneTranscript_Infor[1] = "transcript_id \"" + TranscriptID_change + "\""
                                GeneTranscript_Infor_output = "; ".join(GeneTranscript_Infor)
                                element.pop(8)
                                element_output = "\t".join(element)
                                fo.write(element_output + "\t" + GeneTranscript_Infor_output + "\n")
                        else:
                                fo.write(line + "\n")
                except:
                        pass
fi.close()

fo.close()
