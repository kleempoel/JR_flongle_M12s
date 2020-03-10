library(tidyverse)
data=read_table2('blast/jr_flongle_m12s_allseqcut_mito.txt',col_names = FALSE)
colnames(data)=c('seqid', 'sseqid', 'sacc', 'staxids', 'evalue', 'pident', 'nident', 'slen', 'qstart', 'qend', 'length', 'mismatch')

#hist match length
hist(data$length,breaks = 100,col = 'black')

data=data %>% filter(length>150)


accessions=data %>% group_by(sacc) %>% summarise(count=n()) %>% arrange(-count)
library(taxize)
accessions$taxid=NA;accessions$common_name=NA;accessions$scientific_name=NA
accessions$Genus=NA;accessions$Family=NA;accessions$Order=NA;
accessions$Class=NA;accessions$Kingdom=NA;accessions$ref_name=NA;
for (iseq in 1:nrow(accessions)){
  temptax=genbank2uid(id = accessions$sacc[iseq])
  accessions$taxid[iseq]=as.numeric(unlist(temptax[[1]])[1])
  temptax2=unlist(id2name(accessions$taxid[iseq], db = "ncbi"))
  accessions$scientific_name[iseq]=temptax2[2]
  t2=tax_name(query = c(accessions$scientific_name[iseq]), 
              get = c("Genus","Family","Order","Class",'Kingdom'), db = "ncbi")
  t=unlist(sci2comm(scinames=accessions$scientific_name[iseq], db = "ncbi")) 
  if(is_empty(t)==F){
    accessions$common_name[iseq]=unlist(sci2comm(scinames=accessions$scientific_name[iseq], db = "ncbi"))
  }
  accessions$Genus[iseq]=t2$Genus
  accessions$Family[iseq]=t2$Family
  accessions$Order[iseq]=t2$Order
  accessions$Class[iseq]=t2$Class
  accessions$Kingdom[iseq]=t2$Kingdom
  accessions$ref_name[iseq]=attr(temptax[[1]],"name")
  rm(temptax,temptax2,t2)
}
print(which(is.na(accessions$Class)==T))
write_csv(accessions, 'blast/jr_flongle_m12s_allseqcut_mito_accessions.csv', na = "NA")

all_sp=accessions %>% group_by(common_name) %>% summarise(sum(count)) %>% arrange(-`sum(count)`)
all_fam=accessions %>% group_by(Family) %>% summarise(sum(count)) %>% arrange(-`sum(count)`)
all_gen=accessions %>% group_by(Genus) %>% summarise(sum(count)) %>% arrange(-`sum(count)`)


#get current species
current_sp='plateau mouse'
current_sacc=(accessions %>% filter(common_name==current_sp) %>% select(sacc))$sacc
current_blast=data %>% filter(sacc==current_sacc)

### Get fasta sequences
seq_file="Data/fasta_pass_cut.fasta"
out_file="test.fasta"
seq_list=current_blast$seqid
# seq_list=c('d7d1d81f-d61a-47db-8a7b-46486da2b53d','421511a2-37e3-4959-b3cd-ecffe715b783')
write('',file=out_file,append=FALSE)
for (iseq in 1:length(seq_list)){
  seq_name=seq_list[iseq]
  con = file(seq_file, "r")
  found_seq=0
  while(found_seq==0){
    line = readLines(con, n = 1)
    find_text=grep(pattern = seq_name,x = line)
    if (is_empty(find_text)==FALSE){
      print(paste('found',seq_name,iseq,'/',length(seq_list)),quote = FALSE)
      write(line,file=out_file,append=TRUE)
      found_seq=1
    }
    if(found_seq==1){
      find_nextseq=FALSE
      while(find_nextseq==FALSE){
        line = readLines(con, n = 1)
        test_line=grep(pattern = '>',x = line)
        if(is_empty(test_line)==TRUE){
          cat(line,file="output.txt",sep="\n",append=TRUE)
          write(line,file=out_file,append=TRUE)
        }else{
          find_nextseq=TRUE
        }
      }
    }
  }
  close(con)
}

