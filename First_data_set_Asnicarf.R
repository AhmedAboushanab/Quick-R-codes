library(curatedMetagenomicData)
library(phyloseq)
library(ExperimentHub)
library(stringr)
Dataset <- curatedMetagenomicData("AsnicarF_2017.metaphlan_bugs_list.stool", dryrun = FALSE)
dataset2 <- mergeData(Dataset)
experimentData(dataset2)
Asnicarf <- pData(dataset2)
write.csv(Asnicarf, file = "Asnicarf_first_dataset.csv", row.names = F)
Asnicarf_results <- exprs(dataset2)
write.csv(Asnicarf_results, file = "Asnicarf_results_first_dataset.csv", row.names = T)
subdata <- read.csv("Asnicarf_relative abundance_first_dataset.csv")
healthy <- read.csv("subset_data_sampleID.csv")

#modifying column names in subdata
subdata2 <- list()
for(e in 1:ncol(subdata)){
  h = subdata[1 , e]
  if(str_detect(h, "AsnicarF_2017.metaphlan_bugs_list.stool:") == TRUE){
    new = gsub("AsnicarF_2017.metaphlan_bugs_list.stool:", "", h) #removes a pattern from file name 
    subdata2 = append(subdata2, new)
  }else{
    break
  }
}

print(subdata2)

#healthy data ISA
healthydata <- list()
for(f in 1:nrow(healthy)){
  for(e in 1:ncol(subdata)){
  g = healthy[f, 1]
  h = subdata[1 , e]
  if(str_detect(h, str(g)) == TRUE){
    print("true")
    healthydata <- append(healthydata, h)
  }else{
    print("false")
  }
 }
}

write.csv(unlist(Names_list), file = "study_Names.csv", row.names = F)

healthy[1,1]
dim(healthy)
dim(subdata)
