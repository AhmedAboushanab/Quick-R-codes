library(curatedMetagenomicData)
library(phyloseq)
library(ExperimentHub)
library(stringr)

#getting the dataset
Dataset <- curatedMetagenomicData("AsnicarF_2017.metaphlan_bugs_list.stool", dryrun = FALSE)
dataset2 <- mergeData(Dataset)
experimentData(dataset2)
Asnicarf <- pData(dataset2)
write.csv(Asnicarf, file = "Asnicarf_first_dataset.csv", row.names = F)
Asnicarf_results <- exprs(dataset2)
write.csv(Asnicarf_results, file = "Asnicarf_results_first_dataset.csv", row.names = T)

#filtering healthy samples data
subdata <- read.csv("DhakanDB_2019_relative abundance_dataset.csv", header = FALSE)
healthy <- read.csv("Backhed.csv", header = FALSE)

#modifying column names in subdata
subdata2 <- list()
for(e in 1:ncol(subdata)){
  h = subdata[1 , e]
  if(str_detect(h, "AsnicarF_2017.metaphlan_bugs_list.stool:") == TRUE){
    new = gsub("AsnicarF_2017.metaphlan_bugs_list.stool:", "", h) #removes a pattern from file name 
    subdata2 = append(subdata2, new)
  }else{
    print("FALSE")
  }
}

print(subdata2)
head(subdata)

#healthy data ISA
subdata3 <- read.csv("DhakanDB_2019_relative abundance_dataset.csv", header = FALSE)
write.csv(subdata3[1], file = "DhakanDB.csv", row.names = FALSE)

healthydata <- read.csv("DhakanDB.csv", header = FALSE)
for(f in 1:nrow(healthy)){
  for(e in 1:ncol(subdata)){
    g = healthy[f, 1]
    h = subdata[ , e]
    if(str_detect(h, g) == TRUE && !(h[1] %in% healthydata[1, ])){
      print("true")
      healthydata <- cbind(healthydata, h)
    }else{
      print("false")
    }
  }
}

write.csv(healthydata, file = "DhakanDB_healthy.csv", row.names = F)
# two conditions used to prevent repetition 
