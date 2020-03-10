library(tidyverse)
data=read_table2('ecotag/jr_flongle_m12s_clust_ecotag.txt')

print(nrow(data))
#filter identity
data= data %>% filter(`best_identity:M12s_v1_clean_uniqID_uniq_length`>0.85)

#families
fam=data  %>% group_by(family_name) %>% tally()
print(fam)

genus=data  %>% group_by(genus_name) %>% tally()
print(genus)

sci_name=data  %>% group_by(scientific_name) %>% tally()
print(sci_name)
