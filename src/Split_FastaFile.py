#! /usr/bin/python
import sys

fasta_file=sys.argv[1]
output_dir_path=sys.argv[2]

fa_in = open(fasta_file,"r")
fa_Info = []
fa_Seq = []
fa_Num = -1
for line in fa_in.readlines():
        line = line.rstrip()
        if line[0] == ">":
                fa_Info.append(line)
                fa_Num = fa_Num + 1
                fa_Seq.append("")
        else:
                fa_Seq[fa_Num] = fa_Seq[fa_Num] + line

for i in range(fa_Num + 1):
        fa_out = open(output_dir_path+"/"+str(i)+".fa", "w")
        fa_out.write(fa_Info[i] + "\n")
        while len(fa_Seq[i]) > 60:
                fa_out.write(fa_Seq[i][:60] + "\n")
                fa_Seq[i] = fa_Seq[i][60:]
        else:
                fa_out.write(fa_Seq[i] + "\n")
