# R Script
# Author: Adam W Hansen
# Date Created: Jan 26, 2018
# Date Last Modified: Apr 3, 2018

library(getPass)
library(RJDBC)

#USERNAME = readline(prompt = "Username: ")
USERNAME = getPass(msg = "Username: ")
USERNAME
PASSWORD = getPass(msg = "Password: ")
FILTER_PARAM = "polyphen_score"
MIN = 0
FINISH = 1.01
BIN_WIDTH = 0.05
SUBDIR_NAME = paste(c("./", as.character(Sys.Date()), "_", FILTER_PARAM, "_Default_", MIN, "-", FINISH, "_binwidth_", BIN_WIDTH, "_subqueries"), collapse = '')

get_binned_query_series = function(min=0, finish=10, binWidth=0.5){
  # Generate text for fixed query parameters
  preString = paste(
    c("select 
    count(distinct(sname)) as n_samples, gene_symbol, omim_gene_2016, omim_gene_2018
    from freeze2_3.variant_genotype_filtered_qualpass
    where impact = 'MODERATE'
    and (charge_af <= 0.0001 or charge_af is null)
    and (gnomad_af <= 0.0001 or gnomad_af is null)
    and (cast(gnomad_ext_af as float) <= 0.0001 or cast(gnomad_ext_af as float) is null)
    and (pubmed = '' or pubmed is null)
    and (project not in ('SASCHX', 'WGLTEST', 'WGLM', 'WGLBP', 'WGLCARDIO', 'WGLC', 'WGLHES'))
    and (existing_variation = '' or existing_variation is null)
    and (",
      FILTER_PARAM,
      " >= "),
    collapse = ''
  )
  midString = paste(
    c(" and ",
    FILTER_PARAM,
    " < "),
    collapse = ''
  )
  postString = 
    ")
    group by gene_symbol, omim_gene_2016, omim_gene_2018
    having count(distinct(sname)) >= 5
    order by omim_gene_2016 nulls first, omim_gene_2018 nulls last, count(distinct(sname)) desc"
  
  # Create empty vector to store queries
  queries = c()
  
  while (min+binWidth <= finish) {
    # Format query
    queryString = paste(c(preString, min, midString, min+binWidth, postString), collapse = '')
    # Append query to vector
    queries = append(queries, queryString)
    min = min + binWidth
  }
  
  # Return queries
  return(queries)
}

get_file_names = function(min=0, finish=10, binWidth=0.5){ 
  # Create empty vector to store file names
  fileNames = c()
  
  while (min+binWidth <= finish) {
    # Create curFileName
    curFileName = paste(c(min, '-', min+binWidth, '_', FILTER_PARAM, '_Default.tsv'), collapse = '')
    # Append curFileName to vector
    fileNames = append(fileNames, curFileName)
    min = min + binWidth
  }
  # Return file names
  return(fileNames)
}

get_greaterThanOrEqualTo_vals = function(min=0, finish=10, binWidth=0.5){
  # Create empty vector to store values
  vals = c()
  
  while (min+binWidth <= finish) {
    # Append min to vector
    vals = append(vals, min)
    min = min + binWidth
  }
  # Return vals
  return(vals)
}

get_lessThan_vals = function(min=0, finish=10, binWidth=0.5){
  # Create empty vector to store values
  vals = c()
  
  while (min+binWidth <= finish) {
    # Append max to vector
    vals = append(vals, min + binWidth)
    min = min + binWidth
  }
  # Return vals
  return(vals)
}

# Build query strings
queries = get_binned_query_series(MIN, FINISH, BIN_WIDTH)

# Build file names
fileNames = get_file_names(MIN, FINISH, BIN_WIDTH)

# Get min vals
minVals = get_greaterThanOrEqualTo_vals(MIN, FINISH, BIN_WIDTH)

# Get max vals
maxVals = get_lessThan_vals(MIN, FINISH, BIN_WIDTH)

# Make database connection
drv <- JDBC(driverClass = "com.cloudera.impala.jdbc41.Driver", classPath = list.files("/opt/JDBC/jars/ImpalaJDBC41/",pattern="jar$",full.names=T),identifier.quote="`")
impalaConnectionUrl <- paste("jdbc:impala://XXXX;AuthMech=4;SSLTrustStore=/opt/cloudera/security/jks/truststore.jks;SSLTrustStorePwd=XXXX;ssl=1;UID=", USERNAME, "@XXXX;PWD=", PASSWORD, sep="")
conn <- dbConnect(drv, impalaConnectionUrl)

# Create subdir to store individual queries
dir.create(SUBDIR_NAME)

# Initialize master data frame
allResults = NULL

# Execute queries
count = 1
for (query in queries){
  # Execute current query, get results
  curResult = dbGetQuery(conn, query)
  # Rename default col name "EXPR_0" as "n_samples
  names(curResult)[names(curResult)=="EXPR_0"] = "n_samples"
  
  if (length(curResult$gene_symbol > 0)){
    
    # Add two extra columns to current query (min and max query val)
    curResult$">=" = minVals[count]
    curResult$"<" = maxVals[count]
    
    # Append results to master data frame
    allResults = rbind(allResults, curResult)
  }
  
  # Write current query to unique file
  curFileName = fileNames[count]
  write.table(curResult, file=paste(c(SUBDIR_NAME, "/", curFileName), collapse = ''), sep="\t", quote=FALSE, row.names=FALSE)
  
  count = count + 1
}

# Write master data frame to file
write.table(allResults, file=paste(c(as.character(Sys.Date()), "_", FILTER_PARAM, "_Default_", MIN, "-", FINISH, "_binwidth_", BIN_WIDTH, "_all", ".tsv"), collapse = ''), sep="\t", quote=FALSE, row.names=FALSE)

