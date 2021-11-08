#! /usr/bin/python
import sys
import re

gtf_file=sys.argv[1]
TPM_file=sys.argv[2]
FPKM_file=sys.argv[3]
ReadCount_file=sys.argv[4]

Transcript_TPM_Dict={}
f_tpm=open(TPM_file,"r")
f_tpm.readline()
for tpm_line in f_tpm:
        tpm_line=tpm_line.strip()
        element=tpm_line.split("\t")
        transcript=element[0]
        tpm=element[4]
        Transcript_TPM_Dict[transcript]=tpm
f_tpm.close()

Gene_FPKM_Dict={}
f_fpkm=open(FPKM_file,"r")
for fpkm_line in f_fpkm:
        fpkm_line=fpkm_line.strip()
        element=fpkm_line.split("\t")
        gene=element[0]
        fpkm=element[1]
        Gene_FPKM_Dict[gene]=fpkm
f_fpkm.close()

Gene_ReadCount_Dict={}
f_rc=open(ReadCount_file,"r")
for rc_line in f_rc:
        rc_line=rc_line.strip()
        element=rc_line.split("\t")
        gene=element[0]
        readcount=element[1]
        Gene_ReadCount_Dict[gene]=readcount
f_rc.close()

print "transcript_id\tgene_id\tTPM\tFPKM/RPKM\tReadCounts"
f_gtf=open(gtf_file,"r")
for gtf_line in f_gtf:
        gtf_line=gtf_line.strip()

        if not re.match("^#",gtf_line):
                element=gtf_line.split("\t")
                try:
                        Type=element[2]
                        Infor=element[8]

                        if Type == "transcript":
                                Infor_list=Infor.split(";")
                                Gene=Infor_list[0].replace("gene_id ","").replace("\"","").replace(" ","")
                                Transcript=Infor_list[1].replace("transcript_id ","").replace("\"","").replace(" ","")
                        
                                try:
                                        TPM_value=Transcript_TPM_Dict[Transcript]
                                except Exception:
                                        TPM_value="NA"

                                try:
                                        FPKM_value=Gene_FPKM_Dict[Gene]
                                except Exception:
                                        FPKM_value="NA"
                        
                                try:
                                        ReadCount_value=Gene_ReadCount_Dict[Gene]
                                except Exception:
                                        ReadCount_value="NA"

                                print Transcript+"\t"+Gene+"\t"+TPM_value+"\t"+FPKM_value+"\t"+ReadCount_value
                except Exception:
                        pass
