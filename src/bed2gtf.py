#! /usr/bin/python
import sys
BED12_Gene_file=sys.argv[1]


fi=open(BED12_Gene_file,"r")
for line in fi:
        line=line.strip()
        element=line.split("\t")
        Chr=element[0]
        Start=element[1]
        End=element[2]
        Transcript_id=element[3]
        Strand=element[5]
        Exon_count=element[9]
        Exon_size=element[10]
        Exon_relative_start=element[11]
        Gene_id=element[12]

        Exon_size_list=Exon_size.split(",")
        Exon_relative_start_list=Exon_relative_start.split(",")


        Start_gtf=str(int(Start)+1)
        print Chr+"\tHAVANA\ttranscript\t"+Start_gtf+"\t"+End+"\t.\t"+Strand+"\t.\tgene_id \""+Gene_id+"\"; transcript_id \""+Transcript_id+"\";"

        for i in range(0,int(Exon_count)):
                exon_start=int(Start_gtf)+int(Exon_relative_start_list[i])
                exon_end=int(exon_start)+int(Exon_size_list[i])-1
                print Chr+"\tHAVANA\texon\t"+str(exon_start)+"\t"+str(exon_end)+"\t.\t"+Strand+"\t.\tgene_id \""+Gene_id+"\"; transcript_id \""+Transcript_id+"\"; exon_number \""+str(i+1)+"\";"

fi.close()

