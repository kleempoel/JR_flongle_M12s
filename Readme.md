# Testing carnivore diet analysis on Oxford Nanopore Flongle

Execute on server

1. Merge all sequences in one file

`mkdir Data`

`cat fastq_pass/*.fastq > Data/fastq_pass.fastq`

File size is 269Mo

2. Calculate stats with seqkit, filter out long reads, convert to fasta

`seqkit fx2tab Data/fastq_pass.fastq -l -q -n -i -H -j 32 > JR_flongle_m12s_stats.txt`

`seqkit stats Data/fastq_pass.fastq`

file                   format  type  num_seqs      sum_len  min_len  avg_len  max_len
Data/fastq_pass.fastq  FASTQ   DNA    295,234  109,867,923       84    372.1    4,367


`seqkit seq --min-len 100 --max-len 400 Data/fastq_pass.fastq > Data/fastq_pass_100-400bp.fastq`

`seqkit fq2fa Data/fastq_pass_100-400bp.fastq -o Data/fasta_pass.fasta`

`seqkit stats Data/fasta_pass.fasta`

file                   format  type  num_seqs     sum_len  min_len  avg_len  max_len
Data/fasta_pass.fasta  FASTA   DNA    230,394  67,470,630      101    292.8      400


3. Cluster sequences with usearch, then sort by size and remove rare clusters

`../usearch/usearch -cluster_fast Data/fasta_pass.fasta -id 0.80 -threads 32 -sizeout -centroids Data/centroids.fasta -uc Data/clusters.uc -consout Data/jr_m12s_jm_consensus.fasta`

      Seqs  230394 (230.4k)
  Clusters  227232 (227.2k)
  Max size  54
  Avg size  1.0
  Min size  1
Singletons  225799 (225.8k), 98.0% of seqs, 99.4% of clusters
   Max mem  944Mb
      Time  16:24
Throughput  234.1 seqs/sec.


`../usearch/usearch -sortbysize Data/jr_m12s_jm_consensus.fasta -fastaout Data/jr_m12s_jm_consensus_min5.fasta -minsize 5`

`seqkit stats Data/jr_m12s_jm_consensus_min5.fasta`

file                                   format  type  num_seqs  sum_len  min_len  avg_len  max_len
Data/jr_m12s_jm_consensus_min10.fasta  FASTA   DNA         52   14,508      248      279      306


https://drive5.com/usearch/manual/cmd_cluster_fast.html

https://drive5.com/usearch/manual/identity.html

4. Blastn against vertebrate mitochondrial database

`mkdir blast`

`blastn -db ../JR_SoilF1/db_ncbi_vrtmito/ncbi_vrtmito.fasta -query Data/jr_m12s_jm_consensus_min10.fasta -out blast/jr_flongle_m12s_blast_mito.txt -max_target_seqs 1 -perc_identity 80 -outfmt "6 qseqid sseqid sacc staxids evalue pident nident slen qstart qend length mismatch"`

5. Analysing Blastn results
[In this notebook](jr_flongle_m12s_blastn.ipynb)
