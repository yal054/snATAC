# snATAC

#### 0. download datasets

access to snATAC Cool ADMIN: http://epigenomics.sdsc.edu:8086/

download dataset: CEMBA171212_4B
	http://epigenomics.sdsc.edu:8086/download_page?dataset_name=CEMBA171212_4B&output_folder_name=output_17Q4_new

you will need "Sorted bam file", "Full sparse matrix", "Nuclei ID coordinate", "Peak ID coordinate" 

#### 1. run NMF on tf-idf normlization sparse matrix, plot normalized matrix H & W

sklearn.nmf.py
```bash
python sklearn.nmf.py -h
usage: sklearn.nmf.py [-h] [-i INPUTF] [-x XGI] [-y YGI] [-r RANK] [-n N_RUN]
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

sklearn.plotH.R
```bash
Rscript sklearn.plotH.R -h
usage: sklearn.plotH.R [-h] -i INPUT -o OUTPUT

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
python sklearn.nmf.py -i tmp.repl1_CEMBA171212_4B.sparse.npz \
			-x tmp.repl1_CEMBA171212_4B.xgi \
			-y tmp.repl1_CEMBA171212_4B.ygi 
			-o example \
			-r 6 -n 30 -p 0.05 -c 1000 > example.log

# draw heatmap
Rscript sklearn.plotH.R -i example.H.mx -o example.H
Rscript sklearn.plotW.R -i example.W.mx -o example.W
```

#### 2. calculate measurements and output matrix, coorfinates for each cell cluster

sklearn.nmf.stat.py
```bash
python sklearn.nmf.stat.py -h
usage: sklearn.nmf.stat.py [-h] [-m MATRIX] [-x XGI] [-y YGI] [--basis BASIS]
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

example:
```
python sklearn.nmf.stat.py -m example.npz \
			-x example.xgi \
			-y example.ygi \
			--basis example.W.mx -\
			-coef example.H.mx \
			-c 0.2 -o example
```

#### 3. calculate silhouette and plot tSNE using coefficient matrix H

sklearn.nmf.plot.py
```bash
python sklearn.nmf.plot.py -h
usage: sklearn.nmf.plot.py [-h] [--normH NORMH] [--statH STATH]
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
python sklearn.nmf.plot.py --normH example.normH \
			--statH example.statH \
			-p 0.2 -o example
```


#### 4. extract reads barcode from bam file and generate aggregate bam files for each cell cluster

sklearn.nmf.bam.py
```bash
python sklearn.nmf.bam.py -h
usage: sklearn.nmf.bam.py [-h] [--bam BAM] [--statH STATH] [-o OUTPREFIX]

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
python sklearn.nmf.bam.py --bam tmp.repl1_CEMBA171212_4B.sorted.bam \
			--statH example.statH -o example
```






