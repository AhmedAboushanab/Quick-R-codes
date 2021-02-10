library(curatedMetagenomicData)
View(combined_metadata)
data <- data.frame()
for(f in 1:nrow(combined_metadata)) {
  g = combined_metadata[f, ]
  if(combined_metadata$disease == "healthy") {
    print("healthy")
    data = rbind(data, g)
  } else if(combined_metadata$disease != "healthy") {
    print("diseased")
  } else {
    break
  }
}

#write(unlist(data), file = "data.csv")
write.csv(data, file = 'data2.csv')
#?write_csv
#?rbind
