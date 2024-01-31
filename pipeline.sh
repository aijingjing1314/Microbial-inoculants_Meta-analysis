
    # Set work directoryand database
    wd=/b/Metaseqanalysis_Revised
    db=/c/EasyMicrobiome
    PATH=$PATH:${db}/win
    cd ${wd}


    # creat temp
    mkdir -p temp

    # metadata
    csvtk -t stat result/metadata_raw.txt
    cat -A result/metadata_raw.txt | head -n3
    sed 's/\r//' result/metadata_raw.txt > result/metadata.txt
    cat -A result/metadata.txt | head -n3

    # sequencing data
    # zless check sequencing data
    ls -sh seq/
    zless seq/Z1ACK1_1.fastq.gz | head -n4 
    zless seq/Z1ACK1_1.fastq | head | cut -c 1-60
    seqkit stat seq/Z1ACK1_1.fastq.gz
    seqkit stat seq/*.fastq.gz > result/seqkit.txt
    head result/seqkit.txt


    # rush reads merge and rename
    time tail -n+2 result/metadata.txt | cut -f 1 | \
     rush -j 3 "vsearch --fastq_mergepairs seq/{}_1.fastq.gz --reverse seq/{}_2.fastq.gz \
      --fastqout temp/{}.merged.fq --relabel {}."
    # check name
    head temp/`tail -n+2 result/metadata.txt | cut -f 1 | tail -n1`.merged.fq | grep ^@
    

    # integrate renamed reads
    cat temp/*.merged.fq > temp/all.fq
    #check
    ls -lsh temp/all.fq
    head -n 6 temp/all.fq|cut -c1-60


    # Cut primers and quality filter
    time vsearch --fastx_filter temp/all.fq \
      --fastq_stripleft 0 --fastq_stripright 0 \
      --fastq_maxee_rate 0.01 \
      --fastaout temp/filtered.fa
    # check
    head temp/filtered.fa


    # Dereplicate
    vsearch --derep_fulllength temp/filtered.fa \
      --minuniquesize 10 --sizeout --relabel Uni_ \
      --output temp/uniques.fa 
    # check
    ls -lsh temp/uniques.fa
    head -n 2 temp/uniques.fa


    # ASV Denoise: predict biological sequences and filter chimeras
    usearch -unoise3 temp/uniques.fa -minsize 250 \
      -zotus temp/zotus.fa
    sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/otus.fa
    head -n 2 temp/otus.fa

   
    # Reference-based chimera detect
    mkdir -p result/raw
    cp -f temp/otus.fa result/raw/otus.fa


    # Feature table create and filter
    time vsearch --usearch_global temp/filtered.fa \
      --db result/raw/otus.fa \
      --id 0.97 --threads 10 \
    	--otutabout result/raw/otutab.txt 
    sed -i 's/\r//' result/raw/otutab.txt
    head -n6 result/raw/otutab.txt | cut -f 1-6 |cat -A
    csvtk -t stat result/raw/otutab.txt


    # Species annotation.
    vsearch --sintax result/raw/otus.fa \
      --db ${db}/usearch/rdp_16s_v18.fa \
      --sintax_cutoff 0.1 \
      --tabbedout result/raw/otus.sintax 
    head result/raw/otus.sintax | cat -A
    sed -i 's/\r//' result/raw/otus.sintax

    # Remove plasmid and non-Bacteria
    wc -l result/raw/otutab.txt
    Rscript ${db}/script/otutab_filter_nonBac.R \
      --input result/raw/otutab.txt \
      --taxonomy result/raw/otus.sintax \
      --output result/otutab.txt\
      --stat result/raw/otutab_nonBac.stat \
      --discard result/raw/otus.sintax.discard
    wc -l result/otutab.txt
    cut -f 1 result/otutab.txt | tail -n+2 > result/otutab.id
    usearch -fastx_getseqs result/raw/otus.fa \
        -labels result/otutab.id -fastaout result/otus.fa
    awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}'\
        result/raw/otus.sintax result/otutab.id \
        > result/otus.sintax
        
    usearch -otutab_stats result/otutab.txt \
      -output result/otutab.stat
    cat result/otutab.stat


    # Normlize by subsample
    mkdir -p result/alpha
    Rscript ${db}/script/otutab_rare.R --input result/otutab.txt \
      --depth 10000 --seed 1 \
      --normalize result/otutab_rare.txt \
      --output result/alpha/vegan.txt
    usearch -otutab_stats result/otutab_rare.txt \
      -output result/otutab_rare.stat
    cat result/otutab_rare.stat



    # Species annotation
    cut -f 1,4 result/otus.sintax \
      |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g' \
      > result/taxonomy2.txt
    head -n3 result/taxonomy2.txt

    awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
      split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
      result/taxonomy2.txt > temp/otus.tax
    sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | \
      sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
      > result/taxonomy.txt
    head -n3 result/taxonomy.txt
    
    mkdir -p result/tax
    for i in p c o f g;do
      usearch -sintax_summary result/otus.sintax \
      -otutabin result/otutab_rare.txt -rank ${i} \
      -output result/tax/sum_${i}.txt
    done
    sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' result/tax/sum_*.txt
    # 列出所有文件
    wc -l result/tax/sum_*.txt
    head -n3 result/tax/sum_g.txt


