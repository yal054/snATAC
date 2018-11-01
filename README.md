# snATAC

#### 0. download datasets
access to snATAC Cool ADMIN: http://epigenomics.sdsc.edu:8086/
download dataset: CEMBA171212_4B
	http://epigenomics.sdsc.edu:8086/download_page?dataset_name=CEMBA171212_4B&output_folder_name=output_17Q4_new
	you will need "Full sparse matrix", "Nuclei ID coordinate", "Peak ID coordinate" 

#### 1. run NMF on tf-idf normlization sparse matrix, plot normalized matrix H & W

sklearn.nmf.py
`
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
`

sklearn.plotH.R
`
Rscript sklearn.plotH.R -h
usage: sklearn.plotH.R [-h] -i INPUT -o OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input matrix
  -o OUTPUT, --output OUTPUT
                        output file prefix
`

example:
`
python sklearn.nmf.py -i tmp.repl1_CEMBA171212_4B.sparse.npz -x tmp.repl1_CEMBA171212_4B.xgi -y tmp.repl1_CEMBA171212_4B.ygi -o example -r 6 -n 30 -p 0.05 -c 1000 > example.log

`

#### 2. sklearn.nmf.stat.py
#### 3. sklearn.nmf.plot.py
#### 4. sklearn.nmf.bam.py
