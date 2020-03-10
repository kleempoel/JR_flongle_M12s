# Testing carnivore diet analysis on Oxford Nanopore Flongle

Execute on server

1. Merge all sequences in one file

`mkdir Data`

`cat fastq_pass/*.fastq > Data/fastq_pass.fastq`

File size is 269Mo

2. Calculate stats with seqkit, filter out long reads, convert to fasta

`seqkit fx2tab Data/fastq_pass.fastq -l -q -n -i -H -j 32 > JR_flongle_m12s_stats.txt`

`seqkit stats Data/fastq_pass.fastq`

file                  | format | type | num_seqs |     sum_len | min_len | avg_len | max_len

Data/fastq_pass.fastq |FASTQ  | DNA  |  295,234 | 109,867,923   |    84  |  372.1  |  4,367


`seqkit seq --min-len 100 --max-len 400 Data/fastq_pass.fastq > Data/fastq_pass_100-400bp.fastq`

`seqkit fq2fa Data/fastq_pass_100-400bp.fastq -o Data/fasta_pass.fasta`

`seqkit stats Data/fasta_pass.fasta`

**file**|**format**|**type**|**num\_seqs**|**sum\_len**|**min\_len**|**avg\_len**|**max\_len**
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
Data/fasta\_pass.fasta|FASTA|DNA|230,394|67,470,630|101|292.8|400


3.1. Cluster sequences with usearch, then sort by size and remove rare clusters

`../usearch/usearch -cluster_fast Data/fasta_pass.fasta -id 0.80 -threads 32 -sizeout -centroids Data/centroids.fasta -uc Data/clusters.uc -consout Data/jr_m12s_jm_consensus.fasta`

Seqs|230394 (230.4k)
Clusters|115250 (115.2k)
Max size|5989
Avg size|2.0
Min size|1
Singletons|102113 (102.1k), 44.3% of seqs, 88.6% of clusters
Max mem|730Mb
Time|24:47
Throughput|154.9 seqs/sec.

`../usearch/usearch -sortbysize Data/jr_m12s_jm_consensus.fasta -fastaout Data/jr_m12s_jm_consensus_min5.fasta -minsize 5`

`seqkit stats Data/jr_m12s_jm_consensus_min5.fasta`

**file**|**format**|**type**|**num\_seqs**|**sum\_len**|**min\_len**|**avg\_len**|**max\_len**
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
Data/jr\_m12s\_jm\_consensus\_min5.fasta|FASTA|DNA|4,296|1,232,430|178|286.9|336


https://drive5.com/usearch/manual/cmd_cluster_fast.html

https://drive5.com/usearch/manual/identity.html

3.2. Blastn against vertebrate mitochondrial database

`mkdir blast`

`blastn -db ../JR_SoilF1/db_ncbi_vrtmito/ncbi_vrtmito.fasta -query Data/jr_m12s_jm_consensus_min5.fasta -out blast/jr_flongle_m12s_clust_mito.txt -max_target_seqs 1 -perc_identity 80 -outfmt "6 qseqid sseqid sacc staxids evalue pident nident slen qstart qend length mismatch"`

3.3. Analysing Blastn results
[In this notebook](jr_flongle_m12s_blastn.ipynb). But not pretty, will work in R markdown in the future

3.4. Cutting adaptors
Our PCR products were designed for illumina sequencing, they contain adaptors at each end that we need to remove before processing in ecotag

**1st PCR**|**F R**|**Ref**|**Mean length target sequence**|**Barcode length w/o primers**|**Adaptor**|**N**|**Primer**|**Full Primer**
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
MiMammal-U (forward)|F|Ushio 2017|219|169|ACACTCTTTCCCTACACGACGCTCTTCCGATCT|NNNNNN|GGGTTGGTAAATTTCGTGCCAGC|ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNGGGTTGGTAAATTTCGTGCCAGC
MiMammal-U (reverse)|R|Ushio 2017| | |GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT|NNNNNN|CATAGTGGGGTATCTAATCCCAGTTTG|GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNNCATAGTGGGGTATCTAATCCCAGTTTG

Their reverse needs to be searched for as well, as I believe nanopores have no preferential direction.
Rev Forward: GCTGGCACGAAATTTACCAACCCNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
Rev Reverse: CAAACTGGGATTAGATACCCCACTATGNNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC

I started from the raw fastq, converted to fasta
seqkit fq2fa Data/fastq_pass.fastq -o Data/raw_pass.fasta

Forward first
cutadapt -o Data/fasta_pass_cutFOR.fasta Data/raw_pass.fasta -e 0.15 -g ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNGGGTTGGTAAATTTCGTGCCAGC...CAAACTGGGATTAGATACCCCACTATGNNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -j 4 --discard-untrimmed -O 20 > Data/fasta_pass_cutFOR.report.txt

Reverse first
cutadapt -o Data/fasta_pass_cutREV.fasta Data/raw_pass.fasta -e 0.15 -g GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNNCATAGTGGGGTATCTAATCCCAGTTTG...GCTGGCACGAAATTTACCAACCCNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -j 4 --discard-untrimmed -O 20 > Data/fasta_pass_cutREV.report.txt

seqkit seq --min-len 150 --max-len 300 -j 16 Data/fasta_pass_cutREV.fasta > Data/fasta_pass_cutREV_min150.fasta
seqkit seq --min-len 150 --max-len 300 -j 16 Data/fasta_pass_cutFOR.fasta > Data/fasta_pass_cutFOR_min150.fasta
seqkit seq -r -p -v -j 16 Data/fasta_pass_cutFOR_min150.fasta > Data/fasta_pass_cutFOR_min150_rc.fasta
cat Data/fasta_pass_cutFOR_min150_rc.fasta Data/fasta_pass_cutREV_min150.fasta > Data/fasta_pass_cut.fasta
seqkit stat -j 16 Data/*.fasta

We lose more than half of the reads through this step (fasta_pass_cut.fasta) compared to fasta_pass

**file**|**format**|**type**|**num\_seqs**|**sum\_len**|**min\_len**|**avg\_len**|**max\_len**
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
Data/fasta\_pass\_cutFOR.fasta|FASTA|DNA|104,470|18,944,808|3|181.3|2,352
Data/fasta\_pass\_cutFOR\_min150.fasta|FASTA|DNA|97,902|16,125,229|150|164.7|300
Data/fasta\_pass\_cutREV.fasta|FASTA|DNA|27,955|5,047,612|14|180.6|1,774
Data/fasta\_pass\_cutREV\_min150.fasta|FASTA|DNA|26,684|4,454,156|150|166.9|285
Data/fasta\_pass.fasta|FASTA|DNA|230,394|67,470,630|101|292.8|400
Data/raw\_pass.fasta|FASTA|DNA|295,234|109,867,923|84|372.1|4,367
Data/fasta\_pass\_cut.fasta|FASTA|DNA|124,586|20,579,385|150|165.2|300

3.5 Clustering on cut sequences
`../usearch/usearch -cluster_fast Data/fasta_pass_cut.fasta -id 0.80 -threads 32 -sizeout -centroids Data/centroids_cut.fasta -uc Data/clusters_cut.uc -consout Data/jr_m12s_cut_consensus.fasta`

Clustering works much better on these cut sequences

Seqs|124583 (124.6k)
Clusters|23342 (23.3k)
Max size|2805
Avg size|5.3
Min size|1
Singletons|16636 (16.6k), 13.4% of seqs, 71.3% of clusters
Max mem|425Mb
Time|02:15
Throughput|922.8 seqs/sec.

`../usearch/usearch -sortbysize Data/jr_m12s_cut_consensus.fasta -fastaout Data/jr_m12s_cut_consensus_min2.fasta -minsize 2`

`seqkit stats -j 16 Data/jr_m12s_cut_consensus_min2.fasta`

file|format|type|num\_seqs|sum\_len|min\_len|avg\_len|max\_len
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
Data/jr\_m12s\_cut\_consensus\_min5.fasta|FASTA|DNA|6,706|1,118,271|148|166.8|308


3.6. Ecotag
For comparison, I used ecotag to identify the taxa corresponding to OTUs

ecotag -d ../../db/EMBL_191121 -R ../../db/MiMammal-U/M12s_v1_clean_uniqID_uniq_length.fasta Data/jr_m12s_cut_consensus_min2.fasta > ecotag/jr_flongle_m12s_clust_ecotag.fasta
obitab -d -o ecotag/jr_flongle_m12s_clust_ecotag.fasta > ecotag/jr_flongle_m12s_clust_ecotag.txt
obiclean -r 1 -d 1 -H ecotag/jr_flongle_m12s_clust_ecotag.fasta > ecotag/jr_flongle_m12s_clust_ecotag_CleanH.fasta
obitab -d -o ecotag/jr_flongle_m12s_clust_ecotag_CleanH.fasta > ecotag/jr_flongle_m12s_clust_ecotag_CleanH.txt

Results from Ecotag on the clusters
> fam=data  %>% group_by(family_name) %>% tally()
**family\_name**|**n**
:-----:|:-----:
<chr>|<int>
1 Canidae|327
2 Cervidae|2
3 Columbidae|4
4 Cricetidae|117
5 Felidae|408
6 Leporidae|1
7 Muridae|1
8 NA|2946

> genus=data  %>% group_by(genus_name) %>% tally()
**genus\_name**|**n**
:-----:|:-----:
1 Canis|1
2 Neotoma|2
3 Patagioenas|2
4 Rattus|1
5 Urocyon|7
6 NA|3793

**scientific\_name**|**n**
:-----:|:-----:
<chr>|<int>
1 Arvicolinae|19
2 Canidae|319
3 Caniformia|153
4 Canis|1
5 Carnivora|1286
6 Cervidae|1
7 Columbidae|2
8 Cricetidae|91
9 Felidae|394
10 Feliformia|782

4.1. Blasting all sequences against mitochondrial database
Given that the clustering is not optimal, we can also check what the results look like if we blass all 230,394 sequences

`blastn -db ../JR_SoilF1/db_ncbi_vrtmito/ncbi_vrtmito.fasta -query Data/fasta_pass.fasta -out blast/jr_flongle_m12s_allseq_mito.txt -max_target_seqs 1 -perc_identity 80 -outfmt "6 qseqid sseqid sacc staxids evalue pident nident slen qstart qend length mismatch"`

`blastn -db ../JR_SoilF1/db_ncbi_vrtmito/ncbi_vrtmito.fasta -query Data/fasta_pass_cut.fasta -out blast/jr_flongle_m12s_allseqcut_mito.txt -max_target_seqs 1 -perc_identity 80 -num_threads 32 -outfmt "6 qseqid sseqid sacc staxids evalue pident nident slen qstart qend length mismatch"`



4.2. Blasting all sequences against curated M12s database