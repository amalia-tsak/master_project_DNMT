For RNA-seq data the nf-core pipeline was used. 
For WGBS data the order of the scripts is the following: 
1)download
2)concatenate
3)chop
4)deduplicate_reads
5)trim galore
6)map which runs mapping.R script (map and extract meth)
7)run_CHH which runs run_CHH_all.R (if you want to extract methylation also for cytosines not in CpG context)

Quality is checked after different steps using fastqc and multiqc
