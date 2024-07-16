### Assembly pipeline

##### 1. Hifiasm:

```shell
#!/bin/bash
export PATH="/public3/home/sc53149/soft/hifiasm/hifiasm-master/hifiasm:$PATH"

hifiasm -o VF36.asm -t 64 VF36.HiFi.fastq -l 0
```

- ##### Ragtag:

  ```shell
  #!/bin/bash
  awk '/^S/{print ">"$2"\n"$3}' VF36.asm.bp.p_ctg.gfa | fold > VF36.asm.fa
  ragtag.py correct VF36.fa VF36.asm.fa
  ragtag.py scaffold VF36.fa ragtag_output/ragtag.correct.fasta
  ```

##### 2. HiCanu:

```shell
#!/bin/bash
export PATH="/public3/home/sc53149/D2A/sc40090/canu-2.1.1/bin:$PATH"

canu -d hicanu -p VF36_hicanu genomeSize=800m useGrid=false -maxMemory=256 -maxThreads=64 correctedErrorRate=0.01 -pacbio-hifi VF36.HiFi.fastq
```

##### 3. Flye:

```shell
#!/bin/bash
fastqfile=$1
outpwd=$2
ac=$3
genome_size=$4
minovlp=$5

flye --iterations 0 --pacbio-corr $fastqfile --threads 32 -o $outpwd --asm-coverage $ac --genome-size $genome_size -m $minovlp
```

##### 4. GPM manual editing



### SV loci calling by multiple software

##### 1. Mapping to the reference genome

- pbmm2:

  ```shell
  module load anaconda3/2019-03-new
  if [ "$1" == "-h"  ]; then
          echo "Usage: ./pbmm2.sh reference in.fastq ID out.bam"
          exit 0
  fi
  
  if [ ! -n "$1" ]; then
          echo "./pbmm2.sh reference in.fastq ID out.bam"
          exit 0
  fi
  
  ref=$1
  in=$2
  ID=$3
  out=$4
  pbmm2 align $ref $in $out  --unmapped --preset CCS -N 2 --sort --rg '@RG\tID:$ID\tSM:$ID\tLB:SL\tPL:PB' 
  ```

- NGMLR:

  ```shell
  #!/bin/bash
  export PATH="/public3/home/sc53149/soft/ngmlr/ngmlr-master/bin/ngmlr-0.2.8:$PATH"
  source /public3/soft/modules/module.sh
  module load samtools/1.9-wzm-public3
  
  ngmlr -r SL.fa -q VF36.HiFi.fastq.gz -t 64 -x pacbio --rg-sm VF36 --rg-lb SL --rg-pl PacBio -o ccs_ngmlr.all.bam
  samtools sort ccs_ngmlr.all.bam > ccs_ngmlr.sort.bam
  samtools index ccs_ngmlr.sort.bam
  ```

##### 2. SVs calling 

- PBSV

  ```shell
  #!/bin/bash
  fastqfile=$1
  repeatfile=$2
  out=$3
  ID=$4
  para=$5
  reference=$6
  
  pbsv discover $fastqfile --tandem-repeats $repeatfile $out/$ID.svsig.gz
  pbsv call --ccs -x $para $ref $ID.svsig.gz $out/$ID.var.vcf
  ```

- SVIM

  ```shell
  svim alignment $ID $bamfile
  samtools sort -@ 24 -o $out/$ID.sorted.bam $bamfile
  svim alignment --sequence_alleles  --read_names $out/result  $ID.sorted.bam $ref > $out/$ID.vcf
  ```

- cuteSV

  ```shell
  cuteSV VF36.bam SL.fa cuteSV.vcf /public3/home/sc53149/D2A/sc40090/work/VF36/SV/tmp --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --threads 64 --sample ccs_pbmm2 --min_support 10 --genotype
  ```

- Sniffles

  ```shell
  #step1:
  in=$1
  ID=$2
  out=$3
  node=$4
  #para=$5
  mkdir -p $out
  #sorted
  echo -e "#!/bin/bash \n" > sorted.$ID.sh
  echo "srun -n 1 -c 1 samtools sort -@ 24 -o $out/$ID.sorted.pbmm2.bam $in" >> sorted.$ID.sh
  sbatch -N 1 -p $node sorted.$ID.sh
  
  #step2:
  ref=$1
  ID=$2
  out=$3
  node=$4
  #para=$5
  mkdir -p $out
  
  echo -e "#!/bin/bash \n" > MD.$ID.sh
  echo "srun -n 1 -c 24 samtools calmd -u $out/$ID.sorted.pbmm2.bam $ref > $out/$ID.sorted.pbmm2.MD.bam" >> MD.$ID.sh
  sbatch -N 1 -p $node MD.$ID.sh
  
  #step3:
  ID=$1
  out=$2
  node=$3
  para=$4
  mkdir -p $out
  export PATH="/public1/home/sc40090/soft/T/sc40090/Sniffles-master/bin/sniffles-core-1.0.11:$PATH"
  echo -e "#!/bin/bash \n" > sniffles.$ID.sh
  echo "srun -n 1 -c 1 sniffles  --skip_parameter_estimation -s $para --ccs_reads -m $out/$ID.sorted.pbmm2.MD.bam -v $out/$ID.vcf" >> sniffles.$ID.sh
  sbatch -N 1 -p $node sniffles.$ID.sh
  ```

- SVision

  ```shell
  #!/bin/bash
  
  # conda activate svision-env
  export LD_LIBRARY_PATH=/public3/home/sc53149/D2A/sc40090/anaconda3/envs/svision-env/lib/python3.6/site-packages/pysam/:$LD_LIBRARY_PATH
  
  srun -n 1 -c 64 SVision -b VF36.all.filter.bam -m ./CNN/svision-cnn-model.ckpt.meta -n VF36 -s 5 -t 64 --min_sv_size 50 -g SL.fa -o ./svision_out
  ```

- Assemblytics

  - genome alignment

    ```shell
    srun -n 1 -c 64 nucmer --maxmatch -l 100 -c 500 /public3/home/sc53149/D2A/sc40090/work/VF36/reference/SL4/SL.fa /public3/home/sc53149/D2A/sc40090/work/VF36/assembly/arrow/VF36.fa -p VF36_SL -t 64
    ```

  - Call SV

    ```shell
    Assemblytics VF36_SL.delta.gz VF36_SL 10000 20 100000
    ```

- AnchorWave

  ```shell
  #step1:
  srun -n 1 -c 64 anchorwave gff2seq -i SL_gene_models.gff -r SL.fa -o cds.fa
  srun -n 1 -c 64 minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 SL.fa cds.fa > ref.sam
  srun -n 1 -c 64 minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 VF36.rechr.fa cds.fa > VF36.sam
  #step2:
  srun -n 1 -c 64 anchorwave genoAli -i SL_gene_models.gff -as cds.fa -r SL.fa -a VF36.sam -ar ref.sam -s VF36.rechr.fa -v VF36.vcf -n VF36.anchors -o VF36.maf -f VF36.f.maf > VF36.log
  ```

- SVMU

  - alignment by LASTZ

  ```shell
  #!/bin/bash
  if [ "$1" == "-h"  ]; then
          echo "Usage: ./call_SNP_InDel_HiFi_main.sh reference in.fastq ID out"
          exit 0
  fi
  
  if [ ! -n "$1" ]; then
          echo "./call_SNP_InDel_HiFi_main.sh reference in.fastq ID out"
          exit 0
  fi
  
  node=$1
  export PATH="/public3/home/sc53149/soft/sv_ccs/svmu-master:$PATH"
  #srun -n 1 -c 64 lastz /public3/home/sc53149/D2A/sc40090/work/VF36/reference/SL4/SL.fa[multiple] /public3/home/sc53149/D2A/sc40090/work/VF36/assembly/arrow/VF36.fa[multiple] --chain --format=general:name1,strand1,start1,end1,name2,strand2,start2,end2 --exact=20000 > VF36_SL_lastz.txt
  
  for i in {0..12}
  do
    echo -e "#!/bin/bash \n" > lastz_chr${i}.sh
    echo "srun -n 1 -c 64 lastz SL_chr${i}.fa[multiple] VF36_Chr${i}.fa[multiple] --chain --format=general:name1,strand1,start1,end1,name2,strand2,start2,end2  > VF36_SL_lastz.Chr${i}.txt & " >> lastz_chr${i}.sh
    echo "wait" >> lastz_chr${i}.sh
    sbatch -N 1 -p $node lastz_chr${i}.sh
  done
  ```

  - Call SVs

    ```shell
    #!/bin/bash
    export PATH="/public3/home/sc53149/soft/sv_ccs/svmu-master:$PATH"
    srun -n 1 -c 64 svmu VF36_SL.delta SL.fa VF36.fa I lastz.txt VF36_SL
    ```

- SyRI

  - alignment by minimap2:

    ```shell
    #!/bin/bash
    source /public3/soft/modules/module.sh
    module load samtools/1.9-wzm-public3
    
    minimap2 -ax asm5 --eqx SL.rechr.fa VF36.rechr.fa > VF36_SL.rechr.sam
    samtools view -b VF36_SL.rechr.sam > VF36_SL.rechr.bam
    
    ```

  - Call SVs

    ```shell
    #!/bin/bash
    source activate py35
    
    export PATH="/public3/home/sc53149/soft/sv_ccs/syri-master/syri/bin:$PATH"
    export LD_LIBRARY_PATH=/public3/home/sc53149/anaconda3/lib:$LD_LIBRARY_PATH
    #
    syri -c VF36_SL.rechr.bam -r SL.rechr.fa -q VF36.rechr0.fa -k -F B --nosnp
    ```



### SV visualization

##### 1. Samplot

```shell
#!/bin/bash
if [ "$1" == "-h"  ]; then
        echo "Usage: ./samplot.sh chr start end type node"
        exit 0
fi

if [ ! -n "$1" ]; then
        echo "./samplot.sh chr start end type node"
        exit 0
fi
c=$1#chromesome
s=$2#SV locus start position
e=$3#SV locus end position
t=$4#SV locus type
node=$5

export PATH="/public3/home/sc53149/D2A/sc40090/soft/samplot/bin:$PATH"

echo -e "#!/bin/bash \n" > ${c}_${s}_${e}_${t}.sh
echo "srun -n 1 -c 64 samplot plot \
    -n ngmlr pbmm2 \
    -b ngmlr.all.bam \
      pbmm2.all.bam \
    -o ./result/${c}_${s}_${e}_${t}.png \
    -c ${c} \
    -s ${s} \
    -e ${e} \
    -t ${t}" >> ${c}_${s}_${e}_${t}.sh
sbatch -N 1 -p $node ${c}_${s}_${e}_${t}.sh
```

##### 2. IGV

```shell
#!/bin/bash

# pwd
reference=$1
# pwd
snapshot=$2
# pwd
regionfile=$3

n=-1
for i in "$@";
do
  ((n++))
  align[$n]=$i
done

echo "new"
echo "genome $reference"
echo "snapshotDirectory $snapshot"
for i in "${align[@]}";
do
  echo "load $i"
done

echo "sort position"
echo "expand"

cat ${regionfile} | while read line ;
do

  str=$line
  OLD_IFS="$IFS"
  #IFS="\t"
  arr=($str)
  IFS="$OLD_IFS"


  chr=${arr[0]}
  region1=${arr[1]}
  region2=${arr[2]}
  type=${arr[3]}
  echo "goto ${chr}:${region1}-${region2}"
  echo "snapshot ${chr}:${region1}-${region2}-${type}.png"
done

####example for running:
#bash igv_batch.sh SL.fa /manual_check_igv/INS/ngmlr_cutesv/ ngmlr.cutesv.ins.pos.bed pbmm2.all.bam ngmlr.all.bam > ngmlr.cutesv.ins.sh
```



### SV loci merging

###### obtained SV regions in the reference genome

##### 1. code:

```shell
python svloci_merge.py -i sv.info.sort.bed -o sv.merge.txt -t 500
```

> sv.info.sort.bed: chr\tstart\tend\ttype\tcaller; filtered after visualization.

##### 2. manual inspection



### **Localization of SV regions in the query genome**

##### 1. Extension

###### > SV region in the reference genome extended by 500 bp upstream and downstream

```shell
python SV_custom.py -i sv.merge.txt -e SVregion.ex.txt -m expand -n 500
```

##### 2. Re-mapping

```shell
#!/bin/bash

#SBATCH -N 1
#SBATCH -n 64
#SBATCH --cpus-per-task=6
#SBATCH -p hebhcnormal01
#SBATCH --mem-per-cpu=2G

source ~/.bashrc
conda activate truvari
module load apps/samtools/1.9/gcc-7.3.1
module load apps/bedtools/2.29.1/gcc-7.3.1

dirname="/public/home/work/"
mkdir -p ${dirname}/run_code
mkdir -p ${dirname}/anchor_qrypos
outdir=${dirname}/anchor_qrypos

for chr in {1..12}
do
  inputname="chr_${chr}.txt"
  inputfile=${dirname}/expand_bychr/$inputname
  cat $inputfile | while read refchr refs refe type aligner
  do
    echo -e "#!/bin/bash" > \
    ${dirname}/run_code/chr${refchr}_${refs}_${refe}.${aligner}.sh
    echo "mkdir -p ${outdir}/chr${refchr}_${refs}_${refe}
    cd ${outdir}/chr${refchr}_${refs}_${refe}" >> \
    ${dirname}/run_code/chr${refchr}_${refs}_${refe}.${aligner}.sh

    if [ "$aligner" == "pbmm2"  ]; then
      echo "
      samtools view -hb /public/home/sunmiao/work/cui/complex_check/comsvcheck1031/INV/VF36.all.filter.bam ${refchr}:${refs}-${refe} > chr${refchr}_${refs}-${refe}.bam

      samtools fastq chr${refchr}_${refs}-${refe}.bam --reference /public/share/acirymd555/data/SL4.0/SL.fa > chr${refchr}_${refs}-${refe}.fq

      pbmm2 align /public/share/acirymd555/data/VF36/VF36.asm.rechr.fa chr${refchr}_${refs}-${refe}.fq chr${refchr}_${refs}-${refe}.re${aligner}.bam --preset CCS --sort --rg '@RG\tID:VF36\tSM:VF36\tLB:SL\tPL:PB'
      " >> ${dirname}/run_code/chr${refchr}_${refs}_${refe}.${aligner}.sh
    else
      echo "
      samtools view -hb /public/home/sunmiao/work/cui/complex_check/comsvcheck1031/INV/ccs_ngmlr.sort.filter.bam ${refchr}:${refs}-${refe} > chr${refchr}_${refs}-${refe}.bam

      samtools fastq chr${refchr}_${refs}-${refe}.bam --reference /public/share/acirymd555/data/SL4.0/SL.fa > chr${refchr}_${refs}-${refe}.fq

      ngmlr -r /public/share/acirymd555/data/VF36/VF36.asm.rechr.fa -q chr${refchr}_${refs}-${refe}.fq -x pacbio --rg-sm VF36 --rg-lb SL --rg-pl PacBio -o chr${refchr}_${refs}-${refe}.unsort.bam

      samtools sort chr${refchr}_${refs}-${refe}.unsort.bam > chr${refchr}_${refs}-${refe}.re${aligner}.bam
      rm chr${refchr}_${refs}-${refe}.unsort.bam
      " >> ${dirname}/run_code/chr${refchr}_${refs}_${refe}.${aligner}.sh
    fi
    echo "
    bedtools bamtobed -i chr${refchr}_${refs}-${refe}.re${aligner}.bam > chr${refchr}_${refs}-${refe}.bed

    python ${dirname}/SV_custom.py -b chr${refchr}_${refs}-${refe}.bed -rc ${refchr} -rs ${refs} -re ${refe} -q ${dirname}/ref2qry.pos.txt -m bedtopos
    " >> ${dirname}/run_code/chr${refchr}_${refs}_${refe}.${aligner}.sh
    sbatch ${dirname}/run_code/chr${refchr}_${refs}_${refe}.${aligner}.sh &
    wait
  done
  wait
done
```



### **Pre-filtering of SV regions containing INVs**

1. employed `mummer` mode to align SV regions in the reference genome with the entire corresponding chromosomes of the query genome.

```shell
#!/bin/bash

#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH -p hebhcnormal01
source ~/.bashrc
conda activate mummer
module load apps/samtools/1.9/gcc-7.3.1

out=$1
chr=$2
refchr=$3
refs=$4
refe=$5
reffa=$6
###alignment by chr
samtools faidx ${reffa} ${refchr}:${refs}-${refe} > ${out}/refseq/chr${refchr}_${refs}-${refe}.fa
##
mummer -mum -b -c ${out}/qryseq/VF36_${chr}.fa ${out}/refseq/chr${refchr}_${refs}-${refe}.fa > ${out}/mummer_result/chr${refchr}/chr${refchr}_${refs}_${refe}/${chr}.mums
mummer -mum -b -l 200 -c ${out}/qryseq/VF36_${chr}.fa ${out}/refseq/chr${refchr}_${refs}-${refe}.fa > ${out}/mummer_result/chr${refchr}/chr${refchr}_${refs}_${refe}/${chr}.filter200.mums
###
python mummerplot_format.py -c ${out}/mummer_result/chr${refchr}/chr${refchr}_${refs}_${refe}/${chr}.mums -p ${out}/mummer_result/chr${refchr}/chr${refchr}_${refs}_${refe}/${chr}.format.mums
python mummerplot_format.py -c ${out}/mummer_result/chr${refchr}/chr${refchr}_${refs}_${refe}/${chr}.filter200.mums -p ${out}/mummer_result/chr${refchr}/chr${refchr}_${refs}_${refe}/${chr}.filter200.format.mums
###visualization
cd ${out}/mummer_result/chr${refchr}/chr${refchr}_${refs}_${refe}/
mummerplot --postscript --prefix=${chr}.${aligner} ${out}/mummer_result/chr${refchr}/chr${refchr}_${refs}_${refe}/${chr}.format.mums
cd ${out}/mummer_result/chr${refchr}/chr${refchr}_${refs}_${refe}/
mummerplot --postscript --prefix=${chr}.filter200 ${out}/mummer_result/chr${refchr}/chr${refchr}_${refs}_${refe}/${chr}.filter200.format.mums
###
gnuplot ${out}/mummer_result/chr${refchr}/chr${refchr}_${refs}_${refe}/${chr}.gp
gnuplot ${out}/mummer_result/chr${refchr}/chr${refchr}_${refs}_${refe}/${chr}.filter200.gp

```

2. manual inspection

### **SV region alignment**

```shell
#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH -p hebhcnormal01

source ~/.bashrc
conda activate mummer
module load apps/samtools/1.9/gcc-7.3.1
module load apps/bedtools/2.29.1/gcc-7.3.1

dirname="/public/home/work/"
mkdir -p ${dirname}/mummer
outdir=${dirname}/mummer

cat ${dirname}/mummer/SVpos_info.txt  | while read refchr refs refe qrychr qrys qrye ex1 ex2 ex3 ex4
do
  mkdir -p ${outdir}/chr${refchr}_${refs}_${refe}
  echo -e "#!/bin/bash
  module load apps/samtools/1.9/gcc-7.3.1
  module load apps/R/3.6.3-gcc731
  " > ${outdir}/chr${refchr}_${refs}_${refe}/code_${refchr}.sh
  echo "
  cd ${outdir}/chr${refchr}_${refs}_${refe}" >> \
  ${outdir}/chr${refchr}_${refs}_${refe}/code_${refchr}.sh
  ###
  echo "
  block1=${ex1}
  block2=${ex2}
  block3=${ex3}
  block4=${ex4}
  num=0
  ref1=\`expr ${refs} - \${block1}\`
  ref2=\`expr ${refe} + \${block2}\`
  qry1=\`expr ${qrys} - \${block3}\`
  qry2=\`expr ${qrye} + \${block4}\`
  samtools faidx /public/share/acirymd555/data/SL4.0/SL.fa ${refchr}:\${ref1}-\${ref2} > chr${refchr}_${refs}-${refe}.ref.fa
  samtools faidx /public/share/acirymd555/data/VF36/VF36.asm.rechr.fa ${qrychr}:\${qry1}-\${qry2} > ${qrychr}_${qrys}-${qrye}.qry.fa

  if grep -q \"[Nn]\" chr${refchr}_${refs}-${refe}.ref.fa; then
    echo \"${rechr}:${refs}-${refe}\" >> refpos_N.txt
  else
    dnadiff -p dup\${num}kR ${qrychr}_${qrys}-${qrye}.qry.fa chr${refchr}_${refs}-${refe}.ref.fa
  fi
  " >> ${outdir}/chr${refchr}_${refs}_${refe}/code_${refchr}.sh
  sbatch ${outdir}/chr${refchr}_${refs}_${refe}/code_${refchr}.sh &
  wait
done
wait
```

This process is iterated until manual visualization identifies the UAS location accurately.

### Similarity calculated among copies

```shell
#!/bin/bash

#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH -p hebhcnormal01

source ~/.bashrc
conda activate mummer
module load apps/samtools/1.9/gcc-7.3.1

dirname="/public/home/work/"
cat $dirname/copies.info.txt | while read juge region1 region2 region3
#juge: more copies of the segment in this genome
#region1: the position of single copy
#region2: the position of copy1
#region3: the position of copy2
do
  if [ "$juge" == "qry"  ]; then #two or more copies of the segment in the qry genome
    mkdir -p $dirname/chr${region1}
    cd $dirname/chr${region1}
    samtools faidx SL.fa ${region1} > chr${region1}.fa
    samtools faidx VF36.asm.rechr.fa $region2 > $region2.fa
    samtools faidx VF36.asm.rechr.fa $region3 > $region3.fa
    for i in {$region2,$region3}
    do
      nucmer --mum -p chr${region1}_$i chr${region1}.fa ${i}.fa
      delta-filter -g chr${region1}_${i}.delta > chr${region1}_${i}.filter.delta
      show-snps -H -Clr -T chr${region1}_${i}.filter.delta > chr${region1}_${i}.smallvariant.txt
    done
  else
    mkdir -p $dirname/${region1}
    cd $dirname/${region1}
    samtools faidx VF36.asm.rechr.fa ${region1} > ${region1}.fa
    samtools faidx SL.fa $region2 > chr$region2.fa
    samtools faidx SL.fa $region3 > chr$region3.fa
    for i in {$region2,$region3}
    do
      nucmer --mum -p ${region1}_${i} ${region1}.fa chr${i}.fa
      delta-filter -g ${region1}_${i}.delta > ${region1}_${i}.filter.delta
      show-snps -H -Clr -T ${region1}_${i}.filter.delta > ${region1}_${i}.smallvariant.txt
    done
  fi
done
```

