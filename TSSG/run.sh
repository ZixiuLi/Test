for sample in `ls myseq* | cut -d"." -f 1`
do
        echo ${sample}
        ./tssg ${sample}.fa output/${sample}.txt
done
