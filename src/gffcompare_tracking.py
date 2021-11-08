#! /usr/bin/python
import sys
import json

Transcript_Strand = {}
Transcript_Type = {}

with open('config.json') as input_config:
        config = json.load(input_config)

Transcript_Type_File = config['file_config']['transcript_type']

f1 = open(Transcript_Type_File, "r")
for line in f1:
        line = line.strip()
        element = line.split("\t")
        Transcript = element[3]
        Type = element[4]
        Strand = element[5]

        Transcript_Strand[Transcript] = Strand
        Transcript_Type[Transcript] = Type
f1.close()

gffcompare_tracking_file = sys.argv[1]
lncRNA_BED_file = sys.argv[2]
divergent_file = sys.argv[3]

LncRNA_Strand_Dict = {}
f3 = open(lncRNA_BED_file, "r")
for line in f3:
        line = line.strip()
        element = line.split("\t")
        LncRNA = element[3]
        Strand = element[5]
        LncRNA_Strand_Dict[LncRNA] = Strand
f3.close()

LncRNA_Divergent = {}
f2 = open(divergent_file, "r")
for line in f2:
        line = line.strip()
        element = line.split("\t")
        LncRNA = element[0]
        if LncRNA in Transcript_Type:
                Type = Transcript_Type[LncRNA]
                print LncRNA+"\tdivergent\t"+Type
        else:
                print LncRNA+"\tdivergent\tNovel"
f2.close()

fi = open(gffcompare_tracking_file, "r")
for line in fi:
        line = line.strip()
        element = line.split("\t")
        GenomeLocation = element[3]
        LncRNA = element[4].split(":")[1].split("|")[1]
      
        if element[2] != "-":
                KnownTranscript = element[2].split("|")[1]
                KnownTranscript_Strand = Transcript_Strand[KnownTranscript]
        else:
                KnownTranscript_Strand = "NA"

        LncRNA_Strand = LncRNA_Strand_Dict[LncRNA]

        if (GenomeLocation == "e") or (GenomeLocation == "s") or (GenomeLocation == "p") or (GenomeLocation == "r"):
                print LncRNA+"\tNO_USE\tNovel"
        elif GenomeLocation == "=":
                Type = Transcript_Type[LncRNA]
                print LncRNA+"\t=\t"+Type 
        else:
                if GenomeLocation == "u":
                        GenomeLocation = "intergenic"
                elif GenomeLocation == "i":
                        if KnownTranscript_Strand == LncRNA_Strand:
                                GenomeLocation = "sense"
                        else: 
                                GenomeLocation = "antisense"
                elif GenomeLocation == "x":
                        GenomeLocation = "antisense"
                else:
                        GenomeLocation = "sense"
                print LncRNA+"\t"+GenomeLocation+"\tNovel"
fi.close()

