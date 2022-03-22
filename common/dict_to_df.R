dict_to_df = function(dict) {
  
# adapted from   https://gist.github.com/erinshellman/7405574  
  
  require(plyr)
  df = data.frame()
  df_temp = list()
  store = list()
  
  for (i in 1:length(dict)) {
    
    # Split up the dictionary entry
    split = unlist(strsplit(dict[i], ','))
    split = gsub('\\{', '', split)
    split = gsub('\\}', '', split)
    values = unlist(strsplit(split, ':'))
    
    # Parse out what will be the df headers
    headers = values[seq(1, length(values), 2)]       
    headers = gsub('\"', '', headers) # Remove quotes
    headers = gsub(' ', '', headers)  # and whitespace
    headers = gsub('\'', '', headers)
    
    # Parse out what will be the df values
    row_values = values[seq(0, length(values), 2)]
    row_values = gsub('\"', '', row_values) # Remove quotes
    row_values = gsub(' ', '', row_values)  # and whitespace
    
    # Construct a dataframe with 1 row
    out = data.frame(t(row_values))
    colnames(out) = headers
    
    store[i] = list(out)
    
    if (i %% 1000 == 0) { print(round(i / length(dict), 2)) }
    
  }
  
  # rbind all the dataframes together into one dataframe 
  list_length = length(store)
  
  # If the dictionary is sufficiently large rbind will be slow
  # as all hell, so break the rbinding into multiple steps
  if (list_length >= 3000) {
    
    no_splits = round(list_length / 500)
    chunks = split(store, 1:no_splits)
    
    for (j in 1:no_splits) {
      
      df_temp[j] = list(rbind.fill(chunks[[j]]))
      
    }
    df = rbind.fill(df_temp)
    return(df)
  }
  
  else {
    
    df = rbind.fill(store)
    return(df)
  }
  
}

