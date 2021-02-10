#filtering datasets name
Names <- read.csv("subdata_study_names.csv")
Names_list <- list()
for(f in 1:nrow(Names)){
  g = Names[f, ]
  if(g %in% Names_list){
    print("there")
  }else{
    Names_list <- append(Names_list, g)
  }
}
write.csv(unlist(Names_list), file = "study_Names.csv", row.names = F)

