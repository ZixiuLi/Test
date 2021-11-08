import os
import re
import subprocess
import argparse
import logging
import datetime
import json

def HISAT2_Strandness_option(Strandness, SeqType):
        if Strandness == "first":
                if SeqType == "PE":
                        HISAT2_Strandness = "RF"
                elif SeqType == "SE":
                        HISAT2_Strandness = "R"
        elif Strandness == "second":
                if SeqType == "PE":
                        HISAT2_Strandness = "FR"
                elif SeqType == "SE":
                        HISAT2_Strandness = "F"
        return HISAT2_Strandness

def StringTie_Strandness_option(Strandness):
        if Strandness == "first":
                StringTie_Strandness = "rf"
        elif Strandness == "second":
                StringTie_Strandness = "fr"
        return StringTie_Strandness

def Kallisto_Strandness_option(Strandness):
        if Strandness == "first":
                Kallisto_Strandness = "rf-stranded"
        elif Strandness == "second":
                Kallisto_Strandness = "fr-stranded"
        return Kallisto_Strandness

def HTSeq_Strandness_option(Strandness):
        if Strandness == "first":
                HTSeq_Strandness = "reverse"
        elif Strandness == "second":
                HTSeq_Strandness = "yes"
        return HTSeq_Strandness

def rRNA_Proportion(BamFile, OutputDic):
        if not os.path.isdir(OutputDic):
                os.mkdir(OutputDic)

        with open('config.json') as input_config:
            config = json.load(input_config)
        
        rRNA_BED = config['file_config']['rRNA']

        cmd = '-i {BamFile} -r {rRNA_BED} -o {OutputDic}/rRNA >{OutputDic}/rRNA.txt'
        cmd_R = cmd.format(BamFile = BamFile, rRNA_BED = rRNA_BED, OutputDic = OutputDic)
        cmd_Run = subprocess.Popen(['split_bam.py'] + cmd_R.split())
        cmd_Run.communicate()

def HISAT2_Alignment(InputFile, OutputDic):
        with open('config.json') as input_config:
            config = json.load(input_config)
        
        HISAT2Index = config['file_config']['hisat2index']
        Known_splice_site = config['file_config']['known_splice_site']
       
        sambamba_view_path = "./sambamba_view.sh"
        sambamba_sort_path = "./sambamba_sort.sh"
        sed_1_path = "./sed_1.sh"
        sed_2_path = "./sed_2.sh"
        sed_3_path = "./sed_3.sh"

        if not os.path.isdir(OutputDic):
                os.mkdir(OutputDic)

        fi = open(InputFile,"r")
        for line in fi:
                line = line.strip()
                element = line.split("\t")
                SampleName = element[0]
                Strandness = element[1]
                SeqType = element[2]
                
                OutputDic_SampleName = OutputDic + "/" + SampleName
                if not os.path.isdir(OutputDic_SampleName):
                        os.mkdir(OutputDic_SampleName)

                HISAT2_Strandness = HISAT2_Strandness_option(Strandness,SeqType)
                HISAT2_novel_splice_site_path = OutputDic_SampleName + "/HISAT2_novel_splice_site.txt"
                
                if SeqType == "SE":
                        R1_file=element[3]
                        run_hisat2_1 = subprocess.Popen(['./centos_0.1.sif', 'hisat2', '-p', '10', '--dta', '-x', HISAT2Index,
                                                        '-U', R1_file,
                                                        '--add-chrname', '--rna-strandness', HISAT2_Strandness, '--fr', '--known-splicesite-infile',
                                                        Known_splice_site, '--novel-splicesite-outfile',HISAT2_novel_splice_site_path,
                                                        '--novel-splicesite-infile', HISAT2_novel_splice_site_path, '--seed', '168', '--phred33',
                                                        '--min-intronlen', '20', '--max-intronlen', '500000'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                elif SeqType == "PE":
                        R1_file=element[3]
                        R2_file=element[4]
                        run_hisat2_1 = subprocess.Popen(['./centos_0.1.sif', 'hisat2', '-p', '10', '--dta', '-x', HISAT2Index,
                                                        '-1', R1_file,'-2', R2_file,
                                                        '--add-chrname', '--rna-strandness', HISAT2_Strandness, '--fr', '--known-splicesite-infile',
                                                        Known_splice_site, '--novel-splicesite-outfile',HISAT2_novel_splice_site_path,
                                                        '--novel-splicesite-infile', HISAT2_novel_splice_site_path, '--seed', '168', '--phred33',
                                                        '--min-intronlen', '20', '--max-intronlen', '500000'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                run_hisat2_2 = subprocess.Popen([sambamba_view_path, '-h -S -f bam -t 10 /dev/stdin'], stdin=run_hisat2_1.stdout, stdout=subprocess.PIPE, )
                run_hisat2_3_str = '--sort-by-name -t 10 --tmpdir ' + OutputDic_SampleName + ' -o ' + OutputDic_SampleName + '/' + SampleName + '.sortedbyname.bam /dev/stdin'
                run_hisat2_3 = subprocess.Popen([sambamba_sort_path, run_hisat2_3_str], stdin=run_hisat2_2.stdout, )
                run_hisat2_3.communicate()
                
                stderr_value = run_hisat2_1.communicate()[1]
                HISAT2_log = OutputDic_SampleName + '/' + SampleName + '.HISAT2'
                file(HISAT2_log, 'w').write(stderr_value.encode('utf-8'))

                TimeNow = str(datetime.datetime.now())
                logging.info("Step 1: " + SampleName + " HISAT2 done " + TimeNow)

                run_sambamba_1_str = '-t 10 -h ' + OutputDic_SampleName + '/' + SampleName + '.sortedbyname.bam'
                run_sambamba_1 = subprocess.Popen([sambamba_view_path, run_sambamba_1_str], stdout=subprocess.PIPE, )
                run_sambamba_2 = subprocess.Popen([sed_1_path], stdin=run_sambamba_1.stdout, stdout=subprocess.PIPE, )
                run_sambamba_3 = subprocess.Popen([sed_2_path], stdin=run_sambamba_2.stdout, stdout=subprocess.PIPE, )
                run_sambamba_4 = subprocess.Popen([sed_3_path], stdin=run_sambamba_3.stdout, stdout=subprocess.PIPE, )
                run_sambamba_5_str = '-t 10 -f bam -S -o ' + OutputDic_SampleName + '/' + SampleName + '.sortedbyname.renamed.bam /dev/stdin'
                run_sambamba_5 = subprocess.Popen([sambamba_view_path, run_sambamba_5_str], stdin=run_sambamba_4.stdout)
                run_sambamba_5.communicate()
                os.system("rm -f " + OutputDic_SampleName + '/' + SampleName + '.sortedbyname.bam')

                run_sambamba_6_str = '-t 10 --tmpdir=' + OutputDic_SampleName + ' -o ' + OutputDic_SampleName + '/' + SampleName + '.sortedbycoord.bam ' + OutputDic_SampleName + '/' + SampleName + '.sortedbyname.renamed.bam'
                run_sambamba_6 = subprocess.Popen([sambamba_sort_path, run_sambamba_6_str])
                run_sambamba_6.communicate()
                os.system("rm -f " + OutputDic_SampleName + '/' + SampleName + '.sortedbyname.renamed.bam')
                
                TimeNow = str(datetime.datetime.now())
                logging.info("Step 2: " + SampleName + " sambamba done " + TimeNow)
        fi.close()

def StringTie_Assembly(InputFile, OutputDic):
        with open('config.json') as input_config:
            config = json.load(input_config)
        GTF_file = config['file_config']['StringTie_gtf']

        fi = open(InputFile,"r")
        for line in fi:
                line = line.strip()
                element = line.split("\t")
                SampleName = element[0]
                Strandness = element[1]
                SeqType = element[2]
                
                OutputDic_SampleName = OutputDic + "/" + SampleName
                if not os.path.isdir(OutputDic_SampleName):
                        os.mkdir(OutputDic_SampleName)

                StringTie_Strandness = StringTie_Strandness_option(Strandness)

                StringTie_GTF = OutputDic_SampleName + '/' + SampleName + '.stringtie.gtf'
                BamFile = OutputDic_SampleName + '/' + SampleName + '.sortedbycoord.bam'

                cmd = '--{StringTie_Strandness} -p 20 -G {GTF_file} -o {StringTie_GTF} -l {SampleName} -f 0 -m 200 -a 10 -j 1 -M 1 -g 50 {BamFile}'
                cmd_R = cmd.format(StringTie_Strandness = StringTie_Strandness, GTF_file = GTF_file, StringTie_GTF = StringTie_GTF, SampleName = SampleName, BamFile = BamFile)
                cmd_Run = subprocess.Popen(['./centos_0.1.sif', 'stringtie'] + cmd_R.split())
                cmd_Run.communicate()

                TimeNow = str(datetime.datetime.now())
                logging.info("Step 3: " + SampleName + " StringTie done " + TimeNow)
        fi.close()

def Strawberry_Assembly(InputFile, OutputDic):
        with open('config.json') as input_config:
            config = json.load(input_config)
        GTF_file = config['file_config']['Strawberry_gtf']

        fi = open(InputFile,"r")
        for line in fi:
                line = line.strip()
                element = line.split("\t")
                SampleName = element[0]
                Strandness = element[1]
                SeqType = element[2]
                
                OutputDic_SampleName = OutputDic + "/" + SampleName
                if not os.path.isdir(OutputDic_SampleName):
                        os.mkdir(OutputDic_SampleName)
                
                StringTie_Strandness = StringTie_Strandness_option(Strandness)

                Strawberry_GTF = OutputDic_SampleName + '/' + SampleName + '.strawberry.gtf'
                BamFile = OutputDic_SampleName + '/' + SampleName + '.sortedbycoord.bam'

                if os.path.exists(Strawberry_GTF):
                        os.remove(Strawberry_GTF)
                
                cmd = '--{StringTie_Strandness} -g {GTF_file} -o {Strawberry_GTF} -p 20 -m 0 -t 200 -s 10 -d 50 --no-quant --min-depth-4-transcript 0.1 {BamFile}'
                cmd_R = cmd.format(StringTie_Strandness = StringTie_Strandness, GTF_file = GTF_file, Strawberry_GTF = Strawberry_GTF, BamFile = BamFile)
                cmd_Run = subprocess.Popen(['./centos_0.1.sif', 'strawberry'] + cmd_R.split())
                cmd_Run.communicate()

                TimeNow = str(datetime.datetime.now())
                logging.info("Step 4: " + SampleName + " Strawberry done " + TimeNow)
        fi.close()

def Mkdir_MergeSampleNames(InputFile, OutputDic):
        SampleName_List=[]
        
        fi = open(InputFile,"r")
        for line in fi:
                line = line.strip()
                element = line.split("\t")
                SampleName = element[0]
                SampleName_List.append(SampleName)
        fi.close()

        SampleName_List_Value = "_".join(SampleName_List)
        OutputDic_SampleName_List_Value = OutputDic + '/' + SampleName_List_Value
                
        if not os.path.isdir(OutputDic_SampleName_List_Value):
                os.mkdir(OutputDic_SampleName_List_Value)
        
        TimeNow = str(datetime.datetime.now())
        logging.info("Step 5: Mkdir_MergeSampleNames function done " + TimeNow)

def GTF_fromStringTieMerge_ChangeID(GTF_ori, GTF_new):
                cmd_1 = '{GTF_ori}'
                cmd_1_R = cmd_1.format(GTF_ori = GTF_ori)
                
                fo = open(GTF_ori + '.sort', "w")
                cmd_1_Run = subprocess.Popen(['./gff3sort-master/gff3sort.pl'] + cmd_1_R.split(), stdout = fo)
                cmd_1_Run.communicate()
                fo.close()

                cmd_2 = 'src/GTF_fromStringTieMerge_ChangeID.py {GTF_ori}.sort {GTF_new}'
                cmd_2_R = cmd_2.format(GTF_ori = GTF_ori, GTF_new = GTF_new)
                cmd_2_Run = subprocess.Popen(['./centos_0.1.sif', 'python'] + cmd_2_R.split())
                cmd_2_Run.communicate()


def StringTie_Merge(InputFile, OutputDic):
        with open('config.json') as input_config:
            config = json.load(input_config)
        GTF_file = config['file_config']['StringTie_gtf']

        Mkdir_MergeSampleNames(InputFile, OutputDic)

        SampleName_List=[]
        GTF_List=[]
        fi = open(InputFile,"r")
        for line in fi:
                line = line.strip()
                element = line.split("\t")
                SampleName = element[0]
                SampleName_List.append(SampleName)

                OutputDic_SampleName = OutputDic + "/" + SampleName
                StringTie_GTF = OutputDic_SampleName + '/' + SampleName + '.stringtie.gtf' 
                Strawberry_GTF = OutputDic_SampleName + '/' + SampleName + '.strawberry.gtf'
                GTF_List.append(StringTie_GTF)
                GTF_List.append(Strawberry_GTF)
        fi.close()

        SampleName_List_Value = "_".join(SampleName_List)
        OutputDic_SampleName_List_Value = OutputDic + '/' + SampleName_List_Value
        GTF_List_Value = "\n".join(GTF_List)
        
        with open(OutputDic_SampleName_List_Value + '/GTF_MergeList.txt', 'w') as rsh:
                rsh.write(GTF_List_Value)
        
        Stringtie_merged_gtf = OutputDic_SampleName_List_Value + '/stringtie_merged.gtf'
        MergeList_File = OutputDic_SampleName_List_Value + '/GTF_MergeList.txt'

        cmd = '--merge -p 20 -G {GTF_file} -F 0 -g 0 -f 0 -i -m 0 -T 0 -c 0 -o {Stringtie_merged_gtf} {MergeList_File}'
        cmd_R = cmd.format(GTF_file = GTF_file, Stringtie_merged_gtf = Stringtie_merged_gtf, MergeList_File = MergeList_File)
        cmd_Run = subprocess.Popen(['./centos_0.1.sif', 'stringtie'] + cmd_R.split())
        cmd_Run.communicate()

        Stringtie_merged_gtf_changeID = OutputDic_SampleName_List_Value + '/stringtie_merged.changeID.gtf'
        GTF_fromStringTieMerge_ChangeID(Stringtie_merged_gtf, Stringtie_merged_gtf_changeID)
        
        TimeNow = str(datetime.datetime.now())
        logging.info("Step 6: StringTie_Merge function done " + TimeNow)
        
def GTF_RelatedFile_Generation(gtf_file):
        with open('config.json') as input_config:
            config = json.load(input_config)
        hg38_fasta = config['file_config']['genome_sequence']

        cmd_1 = 'src/Gene_Transcript_from_genegtf.py {gtf_file}'
        cmd_1_R = cmd_1.format(gtf_file = gtf_file)
        f1 = open(gtf_file + '.Gene_Transcript', "w")
        cmd_1_Run = subprocess.Popen(['./centos_0.1.sif', 'python'] + cmd_1_R.split(), stdout = f1)
        cmd_1_Run.communicate()
        f1.close()
        
        cmd_2 = 'src/gtf2bed.pl {gtf_file}'
        cmd_2_R = cmd_2.format(gtf_file = gtf_file)
        f2 = open(gtf_file + '.bed', "w")
        cmd_2_Run = subprocess.Popen(['./centos_0.1.sif', 'perl'] + cmd_2_R.split(), stdout = f2)
        cmd_2_Run.communicate()
        f2.close()

        cmd_3 = '{gtf_file} -g {hg38_fasta} -w {gtf_file}.fa'
        cmd_3_R = cmd_3.format(gtf_file = gtf_file, hg38_fasta = hg38_fasta)
        cmd_3_Run = subprocess.Popen(['./centos_0.1.sif', 'gffread'] + cmd_3_R.split())
        cmd_3_Run.communicate()
       
        cmd_4 = 'index -i {gtf_file}.kallisto_idx {gtf_file}.fa'
        cmd_4_R = cmd_4.format(gtf_file = gtf_file)
        cmd_4_Run = subprocess.Popen(['./centos_0.1.sif', 'kallisto'] + cmd_4_R.split())
        cmd_4_Run.communicate()
        
        cmd_5 = 'src/ExonLength_forGene_fromgtf.R {gtf_file} {gtf_file}.ExonLength'
        cmd_5_R = cmd_5.format(gtf_file = gtf_file)
        cmd_5_Run = subprocess.Popen(['./centos_0.2.sif', 'Rscript'] + cmd_5_R.split())
        cmd_5_Run.communicate()
       
        TimeNow = str(datetime.datetime.now())
        logging.info("Step 7: GTF_RelatedFile_Generation function done " + TimeNow)

def PureLoci_RelatedTo_FilterRNA(BED_file, Pseudogenes_included = "T"):
        with open('config.json') as input_config:
            config = json.load(input_config)
        config = config['file_config']

        if Pseudogenes_included == "T":
                Filtered_RNA = config['Six_types']
        elif Pseudogenes_included == "F":
                Filtered_RNA = config['Five_types']

        cmd_1 = 'intersect -a {BED_file} -b {Filtered_RNA} -s -v -wa'
        cmd_1_R = cmd_1.format(BED_file = BED_file, Filtered_RNA = Filtered_RNA)
        f1 = open(BED_file + '.FilterRNA.NotOverlap', "w")
        cmd_1_Run = subprocess.Popen(['./centos_0.1.sif', 'bedtools'] + cmd_1_R.split(), stdout = f1)
        cmd_1_Run.communicate()
        f1.close()

        cmd_2 = 'intersect -a {BED_file} -b {Filtered_RNA} -s -wa'
        cmd_2_R = cmd_2.format(BED_file = BED_file, Filtered_RNA = Filtered_RNA)
        f2 = open(BED_file + '.FilterRNA.Overlap', "w")
        cmd_2_Run = subprocess.Popen(['./centos_0.1.sif', 'bedtools'] + cmd_2_R.split(), stdout = f2)
        cmd_2_Run.communicate()
        f2.close()

        cmd_3 = 'intersect -a {BED_file}.FilterRNA.NotOverlap -b {BED_file}.FilterRNA.Overlap -s -v -wa'
        cmd_3_R = cmd_3.format(BED_file = BED_file)
        f3 = open(BED_file + '.FilterRNA.NotOverlap_Double', "w")
        cmd_3_Run = subprocess.Popen(['./centos_0.1.sif', 'bedtools'] + cmd_3_R.split(), stdout = f3)
        cmd_3_Run.communicate()
        f3.close()
    
        TimeNow = str(datetime.datetime.now())
        logging.info("Step 8: PureLoci_RelatedTo_FilterRNA function done " + TimeNow)
        return BED_file + '.FilterRNA.NotOverlap_Double'

def BED_2_fasta(BED_file):
        with open('config.json') as input_config:
            config = json.load(input_config)
        hg38_fasta = config['file_config']['genome_sequence']
        
        cmd_1 = 'getfasta -fi {hg38_fasta} -bed {BED_file} -name -split -s -fo {BED_file}.fa.tmp' # if "Segmentation fault" is reported, submit the job using bsub 
        cmd_1_R = cmd_1.format(BED_file = BED_file, hg38_fasta = hg38_fasta)
        cmd_1_Run = subprocess.Popen(['./centos_0.1.sif', 'bedtools'] + cmd_1_R.split())
        cmd_1_Run.communicate()

        cmd_2 = '-w 60 {BED_file}.fa.tmp'
        cmd_2_R = cmd_2.format(BED_file = BED_file)
        f2 = open(BED_file + '.fa', "w")
        cmd_2_Run = subprocess.Popen(['fold'] + cmd_2_R.split(), stdout = f2)
        cmd_2_Run.communicate()
        f2.close()

        TimeNow = str(datetime.datetime.now())
        logging.info("Step 9: BED_2_fasta function done " + TimeNow)
        return BED_file + '.fa'
        
def CPAT(FASTA_file, OutputDic):
        with open('config.json') as input_config:
            config = json.load(input_config)
        CPAT_hex = config['file_config']['CPAT_hex']
        CPAT_logitModel = config['file_config']['CPAT_logitmodel']
        
        cmd_1 = '--gene={FASTA_file} -o {OutputDic}/cpat.out --hex={CPAT_hex} --logitModel={CPAT_logitModel}'
        cmd_1_R = cmd_1.format(FASTA_file = FASTA_file, OutputDic = OutputDic, CPAT_hex = CPAT_hex, CPAT_logitModel = CPAT_logitModel)
        cmd_1_Run = subprocess.Popen(['./centos_0.1.sif', 'cpat.py'] + cmd_1_R.split())
        cmd_1_Run.communicate()
        
        cmd_2 = '{OutputDic}/cpat.out.r'
        cmd_2_R = cmd_2.format(OutputDic = OutputDic)
        cmd_2_Run = subprocess.Popen(['./centos_0.2.sif', 'Rscript'] + cmd_2_R.split())
        cmd_2_Run.communicate()
       
        cmd_3 = "awk 'NR!=1 && $6>=0.364 {print $1}' FS=\"\t\" OFS=\"\t\" " + OutputDic + "/cpat.out | cut -d\"(\" -f 1 | sort | uniq >" + OutputDic + "/cpat.out.CodingTranscript.id"
        os.system(cmd_3)
        
        TimeNow = str(datetime.datetime.now())
        logging.info("Step 10: CPAT function done " + TimeNow)
        return  OutputDic + "/cpat.out.CodingTranscript.id"

def LGC(FASTA_file, OutputDic):
        cmd_1 = '{FASTA_file} {OutputDic}/LGC.out'
        cmd_1_R = cmd_1.format(FASTA_file = FASTA_file, OutputDic = OutputDic)
        cmd_1_Run = subprocess.Popen(['./centos_0.3.sif', 'lgc-1.0.0.py'] + cmd_1_R.split())
        cmd_1_Run.communicate()
        
        cmd_2 = "awk '$5==\"Coding\" {print $1}' FS=\"\t\" OFS=\"\t\" " + OutputDic + "/LGC.out | cut -d\"(\" -f 1 | sort | uniq >" + OutputDic + "/LGC.out.CodingTranscript.id"
        os.system(cmd_2)

        TimeNow = str(datetime.datetime.now())
        logging.info("Step 11: LGC function done " + TimeNow)
        return OutputDic + "/LGC.out.CodingTranscript.id"

def PLEK(FASTA_file, OutputDic):
        cmd_1 = 'PLEK.py -fasta {FASTA_file} -out {OutputDic}/PLEK.out -thread 10'
        cmd_1_R = cmd_1.format(FASTA_file = FASTA_file, OutputDic = OutputDic)
        cmd_1_Run = subprocess.Popen(['./centos_0.1.sif'] + cmd_1_R.split())
        cmd_1_Run.communicate()

        cmd_2 = "awk '$1==\"Coding\" {print $3}' FS=\"\t\" OFS=\"\t\" " + OutputDic + "/PLEK.out | cut -d\">\" -f 2 | cut -d\"(\" -f 1 | sort | uniq >" + OutputDic + "/PLEK.out.CodingTranscript.id"
        os.system(cmd_2)
        
        TimeNow = str(datetime.datetime.now())
        logging.info("Step 12: PLEK function done " + TimeNow)
        return OutputDic + "/PLEK.out.CodingTranscript.id"

def CPPred(FASTA_file, OutputDic):
        with open('config.json') as input_config:
            config = json.load(input_config)

        Human_Hexamer = config['file_config']['CPPred_human_hexamer']
        Human_Range = config['file_config']['CPPred_human_range']
        Human_Model = config['file_config']['CPPred_human_model']
        
        CPPred_path = config['software_config']['CPPred_path']
        
        Current_path = os.getcwd()
        os.chdir(CPPred_path)
        
        cmd_1 = 'CPPred.py -i {FASTA_file} -hex {Human_Hexamer} -r {Human_Range} -mol {Human_Model} -spe Human -o {OutputDic}/CPPred.out'
        cmd_1_R = cmd_1.format(FASTA_file = FASTA_file, Human_Hexamer = Human_Hexamer, Human_Range = Human_Range, Human_Model = Human_Model, OutputDic = OutputDic)
        cmd_1_Run = subprocess.Popen([Current_path+'/centos_0.3.sif', 'python'] + cmd_1_R.split())
        cmd_1_Run.communicate()
       
        os.chdir(Current_path)

        cmd_2 = "awk '$(NF-1)==\"coding\" {print $1}' FS=\"\t\" OFS=\"\t\" " + OutputDic + "/CPPred.out | cut -d\"(\" -f 1 | sort | uniq >" + OutputDic + "/CPPred.out.CodingTranscript.id"
        os.system(cmd_2)

        TimeNow = str(datetime.datetime.now())
        logging.info("Step 13: CPPred function done " + TimeNow)
        return OutputDic + "/CPPred.out.CodingTranscript.id"

def Get_CodingTranscript_BED(gtf_file, BED_file, CodingTranscript_list):
        cmd_1 = 'src/Transcript_Gene_Transcript_inGTF.py {CodingTranscript_list} {gtf_file}'
        cmd_1_R = cmd_1.format(CodingTranscript_list = CodingTranscript_list, gtf_file = gtf_file)
        f1 = open(CodingTranscript_list + '.out', "w")
        cmd_1_Run = subprocess.Popen(['./centos_0.1.sif', 'python'] + cmd_1_R.split(), stdout = f1)
        cmd_1_Run.communicate()
        f1.close()
        
        cmd_2 = 'src/my_join.pl -a {BED_file} -b {CodingTranscript_list}.out -F 4 -f 1'
        cmd_2_R = cmd_2.format(BED_file = BED_file, CodingTranscript_list = CodingTranscript_list)
        f2 = open(CodingTranscript_list + '.out.bed.tmp', "w")
        cmd_2_Run = subprocess.Popen(['./centos_0.1.sif', 'perl'] + cmd_2_R.split(), stdout = f2)
        cmd_2_Run.communicate()
        f2.close()

        cmd_3 = "awk '$13!=\"\"' FS=\"\t\" OFS=\"\t\" " + CodingTranscript_list + ".out.bed.tmp | cut -f 1-12 >" + CodingTranscript_list + ".out.bed"
        os.system(cmd_3)

        os.system("rm -f " + CodingTranscript_list + '.out.bed.tmp')
        return CodingTranscript_list + ".out.bed"

def PureLoci_RelatedTo_CodingRNA(gtf_file, BED_file, FASTA_file, OutputDic):
        CPAT_CodingTranscript_list = CPAT(FASTA_file, OutputDic)
        CPAT_CodingTranscript_BED = Get_CodingTranscript_BED(gtf_file, BED_file, CPAT_CodingTranscript_list)

        LGC_CodingTranscript_list = LGC(FASTA_file, OutputDic)
        LGC_CodingTranscript_BED = Get_CodingTranscript_BED(gtf_file, BED_file, LGC_CodingTranscript_list)

        PLEK_CodingTranscript_list = PLEK(FASTA_file, OutputDic)
        PLEK_CodingTranscript_BED = Get_CodingTranscript_BED(gtf_file, BED_file, PLEK_CodingTranscript_list)
        
        CPPred_CodingTranscript_list = CPPred(FASTA_file, OutputDic)
        CPPred_CodingTranscript_BED = Get_CodingTranscript_BED(gtf_file, BED_file, CPPred_CodingTranscript_list)
        
        cmd_0 = '{CPAT_CodingTranscript_BED} {LGC_CodingTranscript_BED} {PLEK_CodingTranscript_BED} {CPPred_CodingTranscript_BED}'
        cmd_0_R = cmd_0.format(CPAT_CodingTranscript_BED = CPAT_CodingTranscript_BED, LGC_CodingTranscript_BED = LGC_CodingTranscript_BED, PLEK_CodingTranscript_BED = PLEK_CodingTranscript_BED, CPPred_CodingTranscript_BED = CPPred_CodingTranscript_BED)
        output_file = OutputDic + '/codingpotential_transcript.bed'
        f0 = open(output_file, "w")
        cmd_0_Run = subprocess.Popen(['cat'] + cmd_0_R.split(), stdout = f0)
        cmd_0_Run.communicate()
        f0.close()
        
        cmd_1 = 'intersect -a {BED_file} -b {CodingRNA} -s -v -wa'
        cmd_1_R = cmd_1.format(BED_file = BED_file, CodingRNA = output_file)
        f1 = open(BED_file + '.CodingRNA.NotOverlap', "w")
        cmd_1_Run = subprocess.Popen(['./centos_0.1.sif', 'bedtools'] + cmd_1_R.split(), stdout = f1)
        cmd_1_Run.communicate()
        f1.close()

        cmd_2 = 'intersect -a {BED_file} -b {CodingRNA} -s -wa'
        cmd_2_R = cmd_2.format(BED_file = BED_file, CodingRNA = output_file)
        f2 = open(BED_file + '.CodingRNA.Overlap', "w")
        cmd_2_Run = subprocess.Popen(['./centos_0.1.sif', 'bedtools'] + cmd_2_R.split(), stdout = f2)
        cmd_2_Run.communicate()
        f2.close()

        cmd_3 = 'intersect -a {BED_file}.CodingRNA.NotOverlap -b {BED_file}.CodingRNA.Overlap -s -v -wa'
        cmd_3_R = cmd_3.format(BED_file = BED_file)
        f3 = open(BED_file + '.CodingRNA.NotOverlap_Double', "w")
        cmd_3_Run = subprocess.Popen(['./centos_0.1.sif', 'bedtools'] + cmd_3_R.split(), stdout = f3)
        cmd_3_Run.communicate()
        f3.close()
      
        TimeNow = str(datetime.datetime.now())
        logging.info("Step 14: PureLoci_RelatedTo_CodingRNA function done " + TimeNow)

        return BED_file + '.CodingRNA.NotOverlap_Double'

def SelectLongRNA(BED_file):
        cmd = 'src/LengthSelect_from_BED.py {BED_file} 200'
        cmd_R = cmd.format(BED_file = BED_file)
        f0 = open(BED_file + '.Long', "w")
        cmd_Run = subprocess.Popen(['./centos_0.1.sif', 'python'] + cmd_R.split(), stdout = f0)
        cmd_Run.communicate()
        f0.close()
        
        TimeNow = str(datetime.datetime.now())
        logging.info("Step 15: SelectLongRNA function done " + TimeNow)
        return BED_file + '.Long'

def Kallisto(InputFile, OutputDic):
        SampleName_List = []
        
        fi = open(InputFile,"r")
        for line in fi:
                line = line.strip()
                element = line.split("\t")
                SampleName = element[0]
                SampleName_List.append(SampleName)
        fi.close()

        SampleName_List_Value = "_".join(SampleName_List)
        OutputDic_SampleName_List_Value = OutputDic + '/' + SampleName_List_Value
                
        fi = open(InputFile,"r")
        for line in fi:
                line = line.strip()
                element = line.split("\t")
                SampleName = element[0]
                Strandness = element[1]
                SeqType = element[2]
                
                OutputDic_SampleName = OutputDic + "/" + SampleName
                Kallisto_Strandness = Kallisto_Strandness_option(Strandness)

                if SeqType == "SE":
                        R1_file=element[3]
                        cmd_1 = 'quant --{Kallisto_Strandness} -t 20 --bias -i {OutputDic_SampleName_List_Value}/stringtie_merged.changeID.gtf.kallisto_idx -o {OutputDic_SampleName} --single -l 200 -s 20 -b 100 -g {OutputDic_SampleName_List_Value}/stringtie_merged.changeID.gtf {R1_file}'
                        cmd_1_R = cmd_1.format(Kallisto_Strandness = Kallisto_Strandness, OutputDic_SampleName_List_Value = OutputDic_SampleName_List_Value, OutputDic_SampleName = OutputDic_SampleName, R1_file = R1_file)
                        cmd_1_Run = subprocess.Popen(['./centos_0.1.sif', 'kallisto'] + cmd_1_R.split())
                        cmd_1_Run.communicate()

                elif SeqType == "PE":
                        R1_file=element[3]
                        R2_file=element[4]
                        cmd_1 = 'quant --{Kallisto_Strandness} -t 20 --bias -i {OutputDic_SampleName_List_Value}/stringtie_merged.changeID.gtf.kallisto_idx -o {OutputDic_SampleName} -b 100 -g {OutputDic_SampleName_List_Value}/stringtie_merged.changeID.gtf {R1_file} {R2_file}'
                        cmd_1_R = cmd_1.format(Kallisto_Strandness = Kallisto_Strandness, OutputDic_SampleName_List_Value = OutputDic_SampleName_List_Value, OutputDic_SampleName = OutputDic_SampleName, R1_file = R1_file, R2_file = R2_file)
                        cmd_1_Run = subprocess.Popen(['./centos_0.1.sif', 'kallisto'] + cmd_1_R.split())
                        cmd_1_Run.communicate()
                
                TimeNow = str(datetime.datetime.now())
                logging.info("Step 16: " + SampleName + " Kallisto done " + TimeNow)
        fi.close() 

def ReadCount(HISAT2_AlignSummary_file, condition = "Total"):
        number = 1        
        
        fi = open(HISAT2_AlignSummary_file, "r")
        for content in fi:
                if re.search("^Warning", content):
                        pass
                        number = number + 1
        fi.close()

        fi = open(HISAT2_AlignSummary_file, "r")
        for i in range(1,number):
                fi.readline()

        Total_ReadCount = fi.readline().strip().split(" ")[0]
        fi.readline()
        fi.readline()
        UniqueMapped_ReadCount = fi.readline().strip().split("(")[0].replace(" ","")
        fi.close()

        if condition == "Total":
                return Total_ReadCount
        elif condition == "Unique":
                return UniqueMapped_ReadCount
        else:
                pass

def FPKM(InputFile, OutputDic):
        SampleName_List = []
        
        fi = open(InputFile,"r")
        for line in fi:
                line = line.strip()
                element = line.split("\t")
                SampleName = element[0]
                SampleName_List.append(SampleName)
        fi.close()

        SampleName_List_Value = "_".join(SampleName_List)
        OutputDic_SampleName_List_Value = OutputDic + '/' + SampleName_List_Value
                
        fi = open(InputFile,"r")
        for line in fi:
                line = line.strip()
                element = line.split("\t")
                SampleName = element[0]
                Strandness = element[1]
                
                OutputDic_SampleName = OutputDic + "/" + SampleName
                HTSeq_Strandness = HTSeq_Strandness_option(Strandness)

                cmd_1 = '-f bam -r name --stranded={HTSeq_Strandness} -m union {OutputDic_SampleName}/{SampleName}.sortedbycoord.bam {OutputDic_SampleName_List_Value}/stringtie_merged.changeID.gtf'
                cmd_1_R = cmd_1.format(HTSeq_Strandness = HTSeq_Strandness, OutputDic_SampleName = OutputDic_SampleName, SampleName = SampleName, OutputDic_SampleName_List_Value = OutputDic_SampleName_List_Value)
                f1 = open(OutputDic_SampleName + '/HTSEQ_count.txt', "w")
                cmd_1_Run = subprocess.Popen(['./centos_0.1.sif', 'htseq-count'] + cmd_1_R.split(), stdout = f1)
                cmd_1_Run.communicate()
                f1.close()
                
                count = ReadCount(OutputDic_SampleName + "/" + SampleName + ".HISAT2", condition = "Total") 
                
                cmd_2 = 'src/rpkm.pl {OutputDic_SampleName}/HTSEQ_count.txt {OutputDic_SampleName_List_Value}/stringtie_merged.changeID.gtf.ExonLength {count}'
                cmd_2_R = cmd_2.format(OutputDic_SampleName = OutputDic_SampleName, OutputDic_SampleName_List_Value = OutputDic_SampleName_List_Value, count = count)
                f2 = open(OutputDic_SampleName + '/Gene_FPKM.txt', "w")
                cmd_2_Run = subprocess.Popen(['./centos_0.1.sif', 'perl'] + cmd_2_R.split(), stdout = f2)
                cmd_2_Run.communicate()
                f2.close()
                
                TimeNow = str(datetime.datetime.now())
                logging.info("Step 17: " + SampleName + " FPKM done " + TimeNow)
        fi.close()

def Cat_TPM_FPKM_ReadCount(InputFile, OutputDic):
        SampleName_List = []
        
        fi = open(InputFile,"r")
        for line in fi:
                line = line.strip()
                element = line.split("\t")
                SampleName = element[0]
                SampleName_List.append(SampleName)
        fi.close()

        SampleName_List_Value = "_".join(SampleName_List)
        OutputDic_SampleName_List_Value = OutputDic + '/' + SampleName_List_Value
               
        gtf_file = OutputDic_SampleName_List_Value + '/stringtie_merged.changeID.gtf'       

        fi = open(InputFile,"r")
        for line in fi:
                line = line.strip()
                element = line.split("\t")
                SampleName = element[0]
                
                OutputDic_SampleName = OutputDic + "/" + SampleName
                
                Kallisto_output = OutputDic_SampleName + '/abundance.tsv'
                FPKM_file = OutputDic_SampleName + '/Gene_FPKM.txt'
                HTSeq_output = OutputDic_SampleName + '/HTSEQ_count.txt'
                
                cmd_1 = 'src/Transcript_Gene_TPM_FPKM_ReadCounts.py {gtf_file} {Kallisto_output} {FPKM_file} {HTSeq_output}'
                cmd_1_R = cmd_1.format(gtf_file = gtf_file, Kallisto_output = Kallisto_output, FPKM_file = FPKM_file, HTSeq_output = HTSeq_output)
                f1 = open(OutputDic_SampleName + '/Transcript_TPM_FPKM_ReadCounts.txt', "w")
                cmd_1_Run = subprocess.Popen(['./centos_0.1.sif', 'python'] + cmd_1_R.split(), stdout = f1)
                cmd_1_Run.communicate()
                f1.close()
        fi.close()

        TimeNow = str(datetime.datetime.now())
        logging.info("Step 18: Cat_TPM_FPKM_ReadCount function done " + TimeNow)

def SelectExpress(InputFile, OutputDic):
        SampleName_List = []
        
        fi = open(InputFile,"r")
        for line in fi:
                line = line.strip()
                element = line.split("\t")
                SampleName = element[0]
                SampleName_List.append(SampleName)
        fi.close()

        SampleName_List_Value = "_".join(SampleName_List)
        OutputDic_SampleName_List_Value = OutputDic + '/' + SampleName_List_Value
       
        Expressed_TranscriptID_file = OutputDic_SampleName_List_Value + "/expressed_transcripts.id"
        os.system("rm -f " + Expressed_TranscriptID_file)
        
        BED_file = OutputDic_SampleName_List_Value + "/stringtie_merged.changeID.gtf.bed.FilterRNA.NotOverlap_Double.CodingRNA.NotOverlap_Double.Long"

        fi = open(InputFile,"r")
        for line in fi:
                line = line.strip()
                element = line.split("\t")
                SampleName = element[0]
                
                OutputDic_SampleName = OutputDic + "/" + SampleName
    
                cutoff = os.popen("./centos_0.2.sif Rscript src/FPKM_Cutoff.R {OutputDic_SampleName}/Transcript_TPM_FPKM_ReadCounts.txt {BED_file}".format(OutputDic_SampleName = OutputDic_SampleName, BED_file = BED_file)).read()
                cutoff = cutoff.split('\n')[-2].split()[-1]

                cmd_1 = ("awk 'NR!=1 && ($4>%s && $5>10) {print $1}' FS=\"\t\" OFS=\"\t\" " % cutoff) + OutputDic_SampleName + "/Transcript_TPM_FPKM_ReadCounts.txt >" + OutputDic_SampleName + "/expressed_transcripts.id"
                os.system(cmd_1)           
                
                cmd_2 = '{OutputDic_SampleName}/expressed_transcripts.id'
                cmd_2_R = cmd_2.format(OutputDic_SampleName = OutputDic_SampleName)
                f2 = open(Expressed_TranscriptID_file, "a")
                cmd_2_Run = subprocess.Popen(['cat'] + cmd_2_R.split(), stdout = f2)
                cmd_2_Run.communicate()
                f2.close()
        fi.close()
        
        cmd_3 = 'src/my_join.pl -a {BED_file} -b {Expressed_TranscriptID_file} -F 4 -f 1'
        cmd_3_R = cmd_3.format(BED_file = BED_file, Expressed_TranscriptID_file = Expressed_TranscriptID_file)
        f3 = open(BED_file + '.Expressed_tmp1', "w")
        cmd_3_Run = subprocess.Popen(['./centos_0.1.sif', 'perl'] + cmd_3_R.split(), stdout = f3)
        cmd_3_Run.communicate()
        f3.close()

        cmd_4 = "awk '$13!=\"\"' FS=\"\t\" OFS=\"\t\" " + BED_file + ".Expressed_tmp1 | cut -f 1-12 >" + BED_file + ".Expressed_tmp2"
        os.system(cmd_4)
        
        cmd_5 = 'src/my_join.pl -a {BED_file}.Expressed_tmp2 -b {OutputDic_SampleName_List_Value}/stringtie_merged.changeID.gtf.Gene_Transcript -F 4 -f 2'
        cmd_5_R = cmd_5.format(BED_file = BED_file, OutputDic_SampleName_List_Value = OutputDic_SampleName_List_Value)
        f5 = open(BED_file + '.Expressed_tmp3', "w")
        cmd_5_Run = subprocess.Popen(['./centos_0.1.sif', 'perl'] + cmd_5_R.split(), stdout = f5)
        cmd_5_Run.communicate()
        f5.close()

        cmd_6 = "awk '$1!~/chrKI/ && $1!~/^GL/ && $14!=\"\"' FS=\"\t\" OFS=\"\t\" " + BED_file + ".Expressed_tmp3 | cut -f 1-13 | sort | uniq >" + BED_file + ".Expressed"
        os.system(cmd_6)

        os.system("rm -f " + BED_file + ".Expressed_tmp1")
        os.system("rm -f " + BED_file + ".Expressed_tmp2")
        os.system("rm -f " + BED_file + ".Expressed_tmp3")

        TimeNow = str(datetime.datetime.now())
        logging.info("Step 19: SelectExpress funciton done " + TimeNow)
        return BED_file + ".Expressed"

def Genomic_Structure(BED_file):
        with open('config.json') as input_config:
            config = json.load(input_config)
       
        Protein_coding_TSS = config['file_config']['protein_coding_tss']
        GTF = config['file_config']['protein_coding_gtf']

        cmd_1 = 'src/TSS_fromBED.py {BED_file}'
        cmd_1_R = cmd_1.format(BED_file = BED_file)
        f1 = open(BED_file + ".TSS", "w")
        cmd_1_Run = subprocess.Popen(['./centos_0.1.sif', 'python'] + cmd_1_R.split(), stdout = f1)
        cmd_1_Run.communicate()
        f1.close()
       
        cmd_2 = '-a {Protein_coding_TSS} -b {BED_file}.TSS -w 2000 -Sm'
        cmd_2_R = cmd_2.format(BED_file = BED_file, Protein_coding_TSS = Protein_coding_TSS)
        f2 = open(BED_file + ".divergent.tmp", "w")
        cmd_2_Run = subprocess.Popen(['./centos_0.1.sif', 'windowBed'] + cmd_2_R.split(), stdout = f2)
        cmd_2_Run.communicate()
        f2.close()

        cmd_3 = "cut -f 10 " + BED_file + ".divergent.tmp |  sort | uniq >" + BED_file + ".divergent"
        os.system(cmd_3)

        cmd_4 = 'src/my_join.pl -a {BED_file} -b {BED_file}.divergent -F 4 -f 1'
        cmd_4_R = cmd_4.format(BED_file = BED_file)
        f4 = open(BED_file + '.non-divergent.tmp1', "w")
        cmd_4_Run = subprocess.Popen(['./centos_0.1.sif', 'perl'] + cmd_4_R.split(), stdout = f4)
        cmd_4_Run.communicate()
        f4.close()
        
        cmd_5 = "awk '$14==\"\"' FS=\"\t\" OFS=\"\t\" " + BED_file + ".non-divergent.tmp1 | cut -f 1-13 >" + BED_file + ".non-divergent.tmp2"
        os.system(cmd_5)

        os.system("./centos_0.1.sif python src/bed2gtf.py " + BED_file + ".non-divergent.tmp2 >" + BED_file + ".non-divergent.gtf") 
        os.system("echo " + BED_file + ".non-divergent.gtf >" + BED_file + ".non-divergent_gtf")

        cmd_6 = '-r {GTF} -i {BED_file}.non-divergent_gtf -o {BED_file}.non-divergent_gffcompare'
        cmd_6_R = cmd_6.format(GTF = GTF, BED_file = BED_file)
        cmd_6_Run = subprocess.Popen(['./centos_0.1.sif', 'gffcompare'] + cmd_6_R.split())
        cmd_6_Run.communicate()

        cmd_7 = 'src/gffcompare_tracking.py {BED_file}.non-divergent_gffcompare.tracking {BED_file} {BED_file}.divergent'
        cmd_7_R = cmd_7.format(BED_file = BED_file)
        f7 = open(BED_file + '.Genomic_Structure', "w")
        cmd_7_Run = subprocess.Popen(['./centos_0.1.sif', 'python'] + cmd_7_R.split(), stdout = f7)
        cmd_7_Run.communicate()
        f7.close()
        
        os.system("rm -f " + BED_file + '.divergent')
        os.system("rm -f " + BED_file + '.divergent.tmp')
        os.system("rm -f " + BED_file + '.non-divergent*')        
       
        return BED_file + '.Genomic_Structure'

def Single_Multi_Exon(BED_file):
        cmd_1 = "awk '$10==1 {print $4,\"SE\"}' FS=\"\t\" OFS=\"\t\" " + BED_file + " >" + BED_file + ".SE"
        cmd_2 = "awk '$10>1 {print $4,\"PE\"}' FS=\"\t\" OFS=\"\t\" " + BED_file + " >" + BED_file + ".PE"
        cmd_3 = "cat " + BED_file + ".SE " + BED_file + ".PE >" + BED_file + ".Single_Multi_Exon"
        os.system(cmd_1)
        os.system(cmd_2)
        os.system(cmd_3)
        
        os.system("rm -f " + BED_file + '.SE')
        os.system("rm -f " + BED_file + '.PE')
        
        return BED_file + '.Single_Multi_Exon'

def PromoterSignature(BED_file, TSSG_output, TSS_Extend_Length = 1000):
        with open('config.json') as input_config:
            config = json.load(input_config)
        hg38_fasta = config['file_config']['genome_sequence']
        
        if not os.path.isdir(TSSG_output):
                os.mkdir(TSSG_output)
        else:
                os.system("rm -rf " + TSSG_output)
                os.mkdir(TSSG_output)

        cmd_1 = 'src/Promoter_fromBED.py {BED_file} {TSS_Extend_Length}'
        cmd_1_R = cmd_1.format(BED_file = BED_file, TSS_Extend_Length = TSS_Extend_Length)
        f1 = open(BED_file + ".promoter", "w")
        cmd_1_Run = subprocess.Popen(['./centos_0.1.sif', 'python'] + cmd_1_R.split(), stdout = f1)
        cmd_1_Run.communicate()
        f1.close()

        cmd_2 = 'getfasta -fi {hg38_fasta} -bed {BED_file}.promoter -name -split -s -fo {BED_file}.promoter.fa.tmp'
        cmd_2_R = cmd_2.format(hg38_fasta = hg38_fasta, BED_file = BED_file)
        cmd_2_Run = subprocess.Popen(['./centos_0.1.sif', 'bedtools'] + cmd_2_R.split())
        cmd_2_Run.communicate()
       
        cmd_3 = 'fold -w 60 ' + BED_file +".promoter.fa.tmp >" + BED_file + ".promoter.fa"
        os.system(cmd_3)

        cmd_4 = 'src/Split_FastaFile.py {BED_file}.promoter.fa {TSSG_output}'
        cmd_4_R = cmd_4.format(BED_file = BED_file, TSSG_output = TSSG_output)
        cmd_4_Run = subprocess.Popen(['./centos_0.1.sif', 'python'] + cmd_4_R.split())
        cmd_4_Run.communicate()

        fa_list = os.listdir(TSSG_output)
        current_directory = os.getcwd()
        os.chdir(current_directory + "/TSSG") 
        RunTSSG_path = "./run_tssg.sh"
        for File in fa_list:
                cmd_5 = '{current_directory}/{TSSG_output}/{File} {current_directory}/{TSSG_output}/{File}.out'
                cmd_5_R = cmd_5.format(current_directory = current_directory, TSSG_output = TSSG_output, File = File)
                cmd_5_Run = subprocess.Popen([RunTSSG_path] + cmd_5_R.split())
                cmd_5_Run.communicate()
        
        os.chdir(current_directory)
        os.system("rm -f " + BED_file + ".PromoterSignature")
        for File in fa_list:
                if not re.search("out$", File):
                        cmd_6 = 'src/TSSG.py {current_directory}/{TSSG_output}/{File}.out'
                        cmd_6_R = cmd_6.format(current_directory = current_directory, TSSG_output = TSSG_output, File = File)
                        f6 = open(BED_file + ".PromoterSignature", "a")
                        cmd_6_Run = subprocess.Popen(['./centos_0.1.sif', 'python'] + cmd_6_R.split(), stdout = f6)
                        cmd_6_Run.communicate()
                        f6.close()
        
        os.system("rm -f " + BED_file + '.promoter*')
        os.system("rm -rf " + TSSG_output)
        
        return BED_file + ".PromoterSignature"

def Transcript_Length(BED_file):
        fo =open(BED_file + ".Transcript_Length", "w")
        fi = open(BED_file, "r")
        for line in fi:
                line=line.strip()
                element=line.split("\t")
                T=element[3]
                ExonSize=element[10]
                ExonSize_List=ExonSize.split(",")
                ExonSize_Value=[]
                for i in range(0,len(ExonSize_List)-1):
                        ExonSize_Value.append(int(ExonSize_List[i]))
                
                TranscriptLength=sum(ExonSize_Value)
                fo.write(T+"\t"+str(TranscriptLength)+"\n")
        fi.close()
        fo.close()

        return BED_file + ".Transcript_Length"

def Transcript_FPKM(InputFile, OutputDic):
        SampleName_List = []
        
        fi = open(InputFile,"r")
        for line in fi:
                line = line.strip()
                element = line.split("\t")
                SampleName = element[0]
                SampleName_List.append(SampleName)
        fi.close()

        SampleName_List_Value = "_".join(SampleName_List)
        OutputDic_SampleName_List_Value = OutputDic + '/' + SampleName_List_Value
       
        BED_file = OutputDic_SampleName_List_Value + "/lncRNA.bed"

        Transcript_FPKM_Dict = {}
        fi = open(InputFile,"r")
        for line in fi:
                line = line.strip()
                element = line.split("\t")
                SampleName = element[0]
                
                OutputDic_SampleName = OutputDic + "/" + SampleName
                Transcript_FPKM_file =  OutputDic_SampleName + "/Transcript_TPM_FPKM_ReadCounts.txt"
                
                fj = open(Transcript_FPKM_file, "r")
                fj.readline()
                for line_j in fj:
                        line_j = line_j.strip()
                        element_j = line_j.split("\t")
                        Transcript = element_j[0]
                        FPKM_value = float(element_j[3])
                        try:
                                Transcript_FPKM_Dict[Transcript].append(FPKM_value) 
                        except Exception:
                                Transcript_FPKM_Dict[Transcript] = []
                                Transcript_FPKM_Dict[Transcript].append(FPKM_value)
                fj.close()
        fi.close()        
        
        fo = open(BED_file + ".Transcript_FPKM", "w")
        fk = open(BED_file, "r")
        for line in fk:
                line = line.strip()
                element = line.split("\t")
                Transcript = element[3]
                Transcript_FPKM = Transcript_FPKM_Dict[Transcript]
                Transcript_FPKM_Mean = sum(Transcript_FPKM) / len(Transcript_FPKM)
                fo.write(Transcript+"\t"+str(Transcript_FPKM_Mean)+"\n")
        fk.close()
        fo.close()
        
        return BED_file + ".Transcript_FPKM"

def Features(InputFile, OutputDic):
        SampleName_List = []
        
        fi = open(InputFile,"r")
        for line in fi:
                line = line.strip()
                element = line.split("\t")
                SampleName = element[0]
                SampleName_List.append(SampleName)
        fi.close()

        SampleName_List_Value = "_".join(SampleName_List)
        OutputDic_SampleName_List_Value = OutputDic + '/' + SampleName_List_Value
       
        BED_file = OutputDic_SampleName_List_Value + "/lncRNA.bed"

        Genomic_Structure_output = Genomic_Structure(BED_file)
        Single_Multi_Exon_output = Single_Multi_Exon(BED_file)
        PromoterSignature_output = PromoterSignature(BED_file, "TSSG/input/" + SampleName_List_Value)
        Transcript_Length_output = Transcript_Length(BED_file)
        Transcript_FPKM_output = Transcript_FPKM(InputFile, OutputDic)
        
        f1 = open(Single_Multi_Exon_output, "r")
        Single_Multi_Exon_Dict = {}
        for line in f1:
                line = line.strip()
                element = line.split("\t")
                Transcript = element[0]
                Type = element[1]
                
                Type_new = "NA"
                if Type == "SE":
                        Type_new = "0"
                elif Type  == "PE":
                        Type_new = "1"
                Single_Multi_Exon_Dict[Transcript] = Type_new
        f1.close()

        f2 = open(Genomic_Structure_output, "r")
        Genomic_Struture_Dict = {}
        for line in f2:
                line = line.strip()
                element = line.split("\t")
                Transcript = element[0]
                Type = element[1]
                
                if Type != "NO_USE":
                        Type_new = "NA\tNA\tNA\tNA\tNA"
                        if Type == "divergent":
                                Type_new = "1\t0\t0\t0\t0"
                        elif Type == "intergenic":
                                Type_new = "0\t1\t0\t0\t0"
                        elif Type == "antisense":
                                Type_new = "0\t0\t1\t0\t0"
                        elif Type == "sense":
                                Type_new = "0\t0\t0\t1\t0"
                        elif Type == "=":
                                Type_new = "0\t0\t0\t0\t1"
                        
                        Genomic_Struture_Dict[Transcript] = Type_new
        f2.close()

        f3 = open(Transcript_FPKM_output, "r")
        FPKM_Dict = {}
        for line in f3:
                line = line.strip()
                element = line.split("\t")
                Transcript = element[0]
                FPKM = element[1]

                FPKM_Dict[Transcript] = FPKM
        f3.close()

        f4 = open(Transcript_Length_output, "r")
        Length_Dict = {}
        for line in f4:
                line = line.strip()
                element = line.split("\t")
                Transcript = element[0]
                Length = element[1]

                Length_Dict[Transcript] = Length
        f4.close()

        f5 = open(PromoterSignature_output, "r")
        PromoterSignature_Dict = {}
        for line in f5:
                line = line.strip()
                element = line.split("\t")
                Transcript = element[0]
                PromoterSig = element[1]
                
                if int(PromoterSig) >= 1:
                        PromoterSignature_Dict[Transcript] = "1"
                else:
                        PromoterSignature_Dict[Transcript] = "0"
        f5.close()

        fo = open(OutputDic_SampleName_List_Value + "/lncRNA.features", "w") 
        fo.write("ID\tSE_ME\tdivergent\tintergenic\tantisense\tsense_overlap\tknown\tFPKM\tTranscriptLength\tPromoter\n")
        for T in Genomic_Struture_Dict:
                SE_ME = Single_Multi_Exon_Dict[T]
                Genomic_Struture = Genomic_Struture_Dict[T]
                FPKM_value = FPKM_Dict[T]
                Size = Length_Dict[T]
                Promoter = PromoterSignature_Dict[T]
                output = T+"\t"+SE_ME+"\t"+str(Genomic_Struture)+"\t"+str(FPKM_value)+"\t"+str(Size)+"\t"+str(Promoter)
                fo.write(output+"\n")

def main():
        parser = argparse.ArgumentParser()
        parser.add_argument('--input_file')
        parser.add_argument('--output_dir')
        args = parser.parse_args()

        logging.basicConfig(filename = args.input_file + '.log', level = logging.INFO)
        logging.info("-----------------------------------------")
        logging.info(str(datetime.datetime.now()))
        logging.info("Input file: " + args.input_file)
        logging.info("Ouput directory: " + args.output_dir)
        #HISAT2_Alignment(args.input_file, args.output_dir)
        #StringTie_Assembly(args.input_file, args.output_dir)
        #Strawberry_Assembly(args.input_file, args.output_dir)
        #StringTie_Merge(args.input_file, args.output_dir)  
        
        SampleName_List = []
        fi = open(args.input_file,"r")
        for line in fi:
                line = line.strip()
                element = line.split("\t")
                SampleName = element[0]
                SampleName_List.append(SampleName)

        SampleName_List_Value = '_'.join(SampleName_List)
        GTF_input = args.output_dir + '/' + SampleName_List_Value + '/stringtie_merged.changeID.gtf'

        #GTF_RelatedFile_Generation(GTF_input)
        
        PureLoci_RelatedTo_FilterRNA_BED = PureLoci_RelatedTo_FilterRNA(GTF_input + '.bed')
        PureLoci_RelatedTo_FilterRNA_fasta = BED_2_fasta(PureLoci_RelatedTo_FilterRNA_BED) # bsub
        PureLoci_RelatedTo_CodingRNA_BED = PureLoci_RelatedTo_CodingRNA(GTF_input, PureLoci_RelatedTo_FilterRNA_BED, PureLoci_RelatedTo_FilterRNA_fasta, args.output_dir + '/' + SampleName_List_Value)
        SelectLongRNA(PureLoci_RelatedTo_CodingRNA_BED)

#        Kallisto(args.input_file, args.output_dir)
        FPKM(args.input_file, args.output_dir)
        Cat_TPM_FPKM_ReadCount(args.input_file, args.output_dir)
        
        SelectExpress_BED = SelectExpress(args.input_file, args.output_dir)

        os.system("cp " + SelectExpress_BED + " " + args.output_dir + '/' + SampleName_List_Value + "/lncRNA.bed")
        TimeNow = str(datetime.datetime.now()) 
        logging.info("Step 20: lncRNA.bed generation " + TimeNow)

        os.system("./centos_0.1.sif python src/bed2gtf.py " + args.output_dir + '/' + SampleName_List_Value + "/lncRNA.bed >" + args.output_dir + '/' + SampleName_List_Value + "/lncRNA.gtf")
        TimeNow = str(datetime.datetime.now())
        logging.info("Step 21: lncRNA.gtf generation " + TimeNow)
        logging.info("All 21 steps done (putative lncRNA identification pipeline finished !!!)")

        Features(args.input_file, args.output_dir)
        TimeNow = str(datetime.datetime.now())
        logging.info("Step 22: lncRNA.features generation (fetures detection done)" + TimeNow)
        logging.info("-----------------------------------------")
main()
