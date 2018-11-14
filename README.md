# snATAC

### overall summary
check cluster results in different ranks, and then choose the robust one
work on TSCC
git clone https://github.com/YoungLeeBBS/snATAC

download example files (see step 0) to folder example/

```bash
## -- 1 -- ##
# run snATAC in different ranks (# of cell clusters)
for i in `seq 5 20`;
do echo $i;
echo -e "
#!/bin/bash

#PBS -q home
#PBS -N rank${i}
#PBS -l nodes=1:ppn=2,pmem=8gb,walltime=1:00:00
#PBS -V
#PBS -M xxx@gmail.com
#PBS -m abe
#PBS -A ren-group
#PBS -j oe

module load R

path0=/projects/ps-renlab/yangli/projects/CEMBA
n=1
p=`bc <<< ${i}*2`
prefix=example.r${i}.n\${n}

mkdir example/\$prefix/
python snATAC.nmf.lite.py -i example/tmp.repl1_CEMBA171212_4B.sparse.npz -x example/tmp.repl1_CEMBA171212_4B.xgi -y example/tmp.repl1_CEMBA171212_4B.ygi -o example/\$prefix/\$prefix -r ${i} -n \${n} -p 0.05 -c 1000 > example/\$prefix/\$prefix.log
Rscript snATAC.plotH.R -i example/\$prefix/\$prefix.H.mx -o example/\$prefix/\$prefix.H
Rscript snATAC.plotW.R -i example/\$prefix/\$prefix.W.mx -o example/\$prefix/\$prefix.W
python snATAC.nmf.stat.py -m example/\$prefix/\$prefix.npz -x example/\$prefix/\$prefix.xgi -y example/\$prefix/\$prefix.ygi --basis example/\$prefix/\$prefix.W.mx --coef example/\$prefix/\$prefix.H.mx -c 0.2 -o example/\$prefix/\$prefix
python snATAC.nmf.plot.py --normH example/\$prefix/\$prefix.normH --statH example/\$prefix/\$prefix.statH -p \${p} -o example/\$prefix/\$prefix
" | sed '1d' > qsub/snATAC.nmf.r${i}.qsub
qsub qsub/snATAC.nmf.r${i}.qsub -o log/snATAC.nmf.r${i}.log -e log/snATAC.nmf.r${i}.err
done;

## -- 2 -- ##
# after all the job finished, plot measurment (sparseness)
for i in cellSparse entropy; do echo $i >> example/sta.lite.header; done
cut -f 1 example/example.r5.n1/example.r5.n1.box.sta | sed "1 s/contributes/rank/g" > example/sta.box.lite.header

for i in `seq 5 20`;
do echo $i;
n=1
prefix=example.r${i}.n\${n}

# calculate cell sparseness and entropy using the statH file
Rscript snATAC.statBox.R -i example/$prefix/$prefix.statH -o example/$prefix/$prefix >> example/$prefix/$prefix.sta.txt
cat example/$prefix/$prefix.sta.txt | sed -e "s/^/example\t${i}\t/g" | paste - example/sta.lite.header > example/$prefix/$prefix.sta.tmp
sed '1d' example/$prefix/$prefix.box.sta | cut -f 2 | sed -e "1i ${i}" > example/$prefix/$prefix.box.contributes
sed '1d' example/$prefix/$prefix.box.sta | cut -f 3 | sed -e "1i ${i}" > example/$prefix/$prefix.box.sparseness
sed '1d' example/$prefix/$prefix.box.sta | cut -f 4 | sed -e "1i ${i}" > example/$prefix/$prefix.box.entropy
sed '1d' example/$prefix/$prefix.statH | sed -e "s/^/${s}\t${i}\t/g" > example/$prefix/$prefix.statH.tmp
done;
cat example/*/*.sta.tmp | sed -e "1i samples\tranks\tval\tstat" | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$4,$3}' > example/example.sta.txt
rm -rf $path0/01.nmf/${s}/*/*.sta.tmp
paste example/sta.box.lite.header example/*/*.box.contributes > example/example.contributes.sta.txt
paste example/sta.box.lite.header example/*/*.box.sparseness > example/example.sparseness.sta.txt
paste example/sta.box.lite.header example/*/*.box.entropy > example/example.entropy.sta.txt
rm -rf example/*/*.box.contributes example/*/*.box.sparseness example/*/*.box.entropy
cat example/*/*.statH.tmp | sed -e "1i samples\tranks\txgi\tindex\tclass0\tclass1\tcontributes\tsparseness\tentropy" > example/example.statH.sta.txt
rm -rf example/*/*.statH.tmp
Rscript snATAC.plotBox.R -i example/example.statH.sta.txt -o example/example.statH.sta
done;

## -- 3 -- ##
# plot Sankey diagram
for i in `seq 5 20`;
do echo $i;
n=1
prefix=example.r${i}.n\${n}
awk -v a=${i} 'BEGIN{FS=OFS="\t"}{print "r"a".c"$4*1}' example/$prefix/$prefix.statH > example/$prefix/$prefix.statH.tmp
done;

for j in `seq 5 19`;
do echo $j;
j1=`bc <<< ${j}+1`
echo ${j1}
prefix1=example.r${j}.n\${n}
prefix2=example.r${j1}.n\${n}
paste example/$prefix1/$prefix1.statH.tmp example/$prefix2/$prefix2.statH.tmp | sort | uniq -c | awk 'BEGIN{OFS="\t"}{print $2,$3,$1}' >> example/example.statH.sankey.summary
done;
rm -rf example/*/*.statH.tmp
sed -i -e "1i from\tto\tcount" example/example.statH.sankey.summary
done;

n=0
for i in `seq 3 20`;
do echo $i;
c=`bc <<< ${i}-1`
for j in `seq 0 ${c}`;
do echo $j;
echo -e "$n\tr${i}.c${j}\t${i}" >> example/example.statH.sankey.idx
n=`bc <<< ${n}+1`
done;
done;
sed -i -e "1i idx\tname\tcolor" example/example.statH.sankey.idx

python snATAC.nmf.sankey.py --data example/example.statH.sankey.summary  --idx example/example.statH.sankey.idx -o example/example.statH

## -- 4 -- ##
# extract read for each cell cluster and gerenate aggregate *.bam file
# for example choose rank = 15 after check the measurement box-plot and sankey plot 
python snATAC.nmf.bam.py --bam example/tmp.repl1_CEMBA171212_4B.sorted.bam --statH example/example.r15.n1/example.r15.n1.statH -o example/example.r15.n1/example.r15.n1
```

#### 0. download datasets

access to snATAC Cool ADMIN: http://epigenomics.sdsc.edu:8086/

download dataset: CEMBA171212_4B
	http://epigenomics.sdsc.edu:8086/download_page?dataset_name=CEMBA171212_4B&output_folder_name=output_17Q4_new

you will need "Sorted bam file", "Full sparse matrix", "Nuclei ID coordinate", "Peak ID coordinate" 

#### 1. run NMF on tf-idf normlization sparse matrix, plot normalized matrix H & W

snATAC.nmf.lite.py
```bash
python snATAC.nmf.py -h
usage: snATAC.nmf.py [-h] [-i INPUTF] [-x XGI] [-y YGI] [-r RANK] [-n N_RUN]
                      [-p PROB] [-c CT] [-o OUTPREFIX]

Run NMF using sklearn.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTF, --inputF INPUTF
                        input matrix in npz format
  -x XGI, --xgi XGI     input xgi index
  -y YGI, --ygi YGI     input ygi index
  -r RANK, --rank RANK  an integer for rank
  -n N_RUN, --n_run N_RUN
                        an integer for # of runs
  -p PROB, --prob PROB  an float of prob used for bins quantile filtering
  -c CT, --ct CT        an integer for # of read counts used for filtering
  -o OUTPREFIX, --outPrefix OUTPREFIX
                        output prefix
```

snATAC.plotH.R
```bash
Rscript snATAC.plotH.R -h
usage: snATAC.plotH.R [-h] -i INPUT -o OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input matrix
  -o OUTPUT, --output OUTPUT
                        output file prefix
```

example:
```bash
# run NMF
python snATAC.nmf.lite.py -i tmp.repl1_CEMBA171212_4B.sparse.npz \
			-x tmp.repl1_CEMBA171212_4B.xgi \
			-y tmp.repl1_CEMBA171212_4B.ygi \
			-o example \
			-r 15 -n 1 -p 0.05 -c 1000 > example.log

# draw heatmap
Rscript snATAC.plotH.R -i example.H.mx -o example.H
Rscript snATAC.plotW.R -i example.W.mx -o example.W
```

#### 2. calculate measurements and output matrix, coorfinates for each cell cluster

snATAC.nmf.stat.py
```bash
python snATAC.nmf.stat.py -h
usage: snATAC.nmf.stat.py [-h] [-m MATRIX] [-x XGI] [-y YGI] [--basis BASIS]
                           [--coef COEF] [-c CONTRIBUTE] [-o OUTPREFIX]

Run statistic for NMF using sklearn.

optional arguments:
  -h, --help            show this help message and exit
  -m MATRIX, --matrix MATRIX
                        matrix in npz format
  -x XGI, --xgi XGI     input xgi index
  -y YGI, --ygi YGI     input ygi index
  --basis BASIS         basis matrix W
  --coef COEF           coefficient matrix H
  -c CONTRIBUTE, --contribute CONTRIBUTE
                        an float cutoff for contribute of bins / cells
  -o OUTPREFIX, --outPrefix OUTPREFIX
                        output prefix
```

calculate cell sparseness and entropy using the statH file
```bash
Rscript snATAC.statBox.R -h
usage: snATAC.statBox.R [-h] -i INPUT -o OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input statH
  -o OUTPUT, --output OUTPUT
                        output file prefix
```

example:
```
python snATAC.nmf.stat.py -m example.npz \
			-x example.xgi \
			-y example.ygi \
			--basis example.W.mx -\
			--coef example.H.mx \
			-c 0.2 -o example

# calculate cell sparseness and entropy using the statH file
Rscript snATAC.statBox.R -i example.statH \
			-o example >> example.sta.txt
```

#### 3. calculate silhouette and plot tSNE using coefficient matrix H

snATAC.nmf.plot.py
```bash
python snATAC.nmf.plot.py -h
usage: snATAC.nmf.plot.py [-h] [--normH NORMH] [--statH STATH]
                           [-p PERPLEXITY] [-o OUTPREFIX]

Run NMF using sklearn.

optional arguments:
  -h, --help            show this help message and exit
  --normH NORMH         coefficient matrix H
  --statH STATH         input statH matrix
  -p PERPLEXITY, --perplexity PERPLEXITY
                        an int for perplexity
  -o OUTPREFIX, --outPrefix OUTPREFIX
                        output prefix
```

example:
```bash
python snATAC.nmf.plot.py --normH example.normH \
			--statH example.statH \
			-p 30 -o example
```


#### 4. extract reads barcode from bam file and generate aggregate bam files for each cell cluster

snATAC.nmf.bam.py
```bash
python snATAC.nmf.bam.py -h
usage: snATAC.nmf.bam.py [-h] [--bam BAM] [--statH STATH] [-o OUTPREFIX]

filter bam based on QNAMES

optional arguments:
  -h, --help            show this help message and exit
  --bam BAM             bam file
  --statH STATH         input statH matrix
  -o OUTPREFIX, --outPrefix OUTPREFIX
                        output prefix
```

example:
```bash
python snATAC.nmf.bam.py --bam tmp.repl1_CEMBA171212_4B.sorted.bam \
			--statH example.statH -o example
```


#### 5. determine which rank to choose

when you got cluster results for each rank (the # of clusters), use sparseness (in *.statH file) to determine which rank should choose.
you can calculate the median sparseness from column "sparseness" in *.statH file for each rank.

#### 6. draw Sankey Diagram

when you have cluster results for each rank, you can plot Sankey Diagram to see where cell goes for one rank to another.
Also, you can tell at which rank is the robust cluster result.

you need to create two files by using the information in column "class0" in *.statH file. 
I create them in directory example/
1. CEMBA171212_4B.statH.sankey.summary
2. CEMBA171212_4B.statH.sankey.idx

```bash
python sklearn.nmf.sankey.py -h
usage: sklearn.nmf.sankey.py [-h] [--data DATA] [--idx IDX] [-o OUTPREFIX]

plot sankey diagram

optional arguments:
  -h, --help            show this help message and exit
  --data DATA           summary of links for sankey diagram
  --idx IDX             index for each elemnets
  -o OUTPREFIX, --outPrefix OUTPREFIX
                        output prefix
```

example:

```bash
python sklearn.nmf.sankey.py --data example/CEMBA171212_4B.statH.sankey.summary
			--idx example/CEMBA171212_4B.statH.sankey.idx
			-o example/CEMBA171212_4B.statH
```




