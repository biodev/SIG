library(gdata)

# clean NA function
clean_na = function(data_set){
  
  for (i in 1:dim(data_set)[2]){
    #print(i)    
    if(sum(na.omit(data_set[,i] == "")) > 0){
      data_set[which(data_set[,i] == ""),i] <- NA
    }
    
  }
  return(data_set)
}

