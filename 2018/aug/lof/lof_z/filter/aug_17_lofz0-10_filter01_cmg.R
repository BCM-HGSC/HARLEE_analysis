# R Script
# Author: Adam W Hansen
# Date Created: Jan 26, 2018
# Date Last Modified: Aug 17, 2018

library(getPass)
library(RJDBC)

#USERNAME = readline(prompt = "Username: ")
USERNAME = getPass(msg = "Username: ")
USERNAME
PASSWORD = getPass(msg = "Password: ")
FILTER_PARAM = "exac_lof_z"
MIN = 0
FINISH = 10.01
INCREMENT = 0.01
SUBDIR_NAME = paste(c("./", as.character(Sys.Date()), "_", FILTER_PARAM, "_Default_", MIN, "-", FINISH, "_increment_", INCREMENT, "_subqueries"), collapse = '')

get_binned_query_series = function(min=0, finish=10, increment=0.5){
  # Generate text for fixed query parameters
  preStringCases = paste(
    c("select 
    count(distinct(sname)) as n_samples, gene_symbol, omim_gene_2016, omim_gene_2018
    from 
    (select freeze2_3.variant_genotype_filtered_qualpass.* from freeze2_3.variant_genotype_filtered_qualpass
    left join freeze2_3.filtered_qualpass_maf
        on variant_genotype_filtered_qualpass.vid = filtered_qualpass_maf.vid
    where impact = 'HIGH'
    and f2_3_af < 0.01
    and (af < 0.001 or af is null)
    and (charge_af <= 0.0001 or charge_af is null)
    and (gnomad_af <= 0.0001 or gnomad_ext_af <= 0.0001 or (gnomad_af is null and gnomad_ext_af is null))
    and (pubmed = '' or pubmed is null)
    and (project = 'CMG')
    and (existing_variation = '' or existing_variation is null)
    and variant_genotype_filtered_qualpass.vid not like 'GL%'
    and (domains not like 'Low_complexity%' or domains is null)
    and exac_mu_syn is not null
    and (vr >=4 and vr >= 0.25*dp)
    and not (sname rlike '[Mm]other' or sname rlike '[Ff]ather' or sname rlike '[\\-_0-9][Mm]$' or sname rlike '[\\-_0-9][Ff]$' or sname rlike '[\\-_0-9][Dd]$' or sname rlike '[\\-_]0?[2-3]' or sname rlike '^[0-9]{4}[A-Z]{2}0[2-3]$' or sname like 'WPW%' or sname like 'HGSC%' or sname like 'NA%')
    and (",
      FILTER_PARAM,
      " >= "),
    collapse = ''
  )
  preStringControls = paste(
    c(")
    ) cases
    left anti join
    (select freeze2_3.variant_genotype_filtered_qualpass.vid from freeze2_3.variant_genotype_filtered_qualpass
    left join freeze2_3.filtered_qualpass_maf
        on variant_genotype_filtered_qualpass.vid = filtered_qualpass_maf.vid
    where impact = 'HIGH'
    and f2_3_af < 0.01
    and (af < 0.001 or af is null)
    and (charge_af <= 0.0001 or charge_af is null)
    and (gnomad_af <= 0.0001 or gnomad_ext_af <= 0.0001 or (gnomad_af is null and gnomad_ext_af is null))
    and (pubmed = '' or pubmed is null)
    and (project = 'CMG')
    and (existing_variation = '' or existing_variation is null)
    and variant_genotype_filtered_qualpass.vid not like 'GL%'
    and (domains not like 'Low_complexity%' or domains is null)
    and exac_mu_syn is not null
    and (vr >=4 and vr >= 0.25*dp)
    and (sname rlike '[Mm]other' or sname rlike '[Ff]ather' or sname rlike '[\\-_0-9][Mm]$' or sname rlike '[\\-_0-9][Ff]$' or sname rlike '[\\-_0-9][Dd]$' or sname rlike '[\\-_]0?[2-3]' or sname rlike '^[0-9]{4}[A-Z]{2}0[2-3]$' or sname like 'WPW%' or sname like 'HGSC%' or sname like 'NA%')
    and (",
      FILTER_PARAM,
      " >= "),
    collapse = ''
  )
  postString = ")
    ) controls
    on cases.vid = controls.vid
    group by gene_symbol, omim_gene_2016, omim_gene_2018
    having count(distinct(sname)) >= 5
    order by omim_gene_2016 nulls first, omim_gene_2018 nulls last, count(distinct(sname)) desc"
  
  # Create empty vector to store queries
  queries = c()
  
  while (round(min+increment, digits=3) <= finish) {
    # Format query
    queryString = paste(c(preStringCases, round(min, digits=3), preStringControls, round(min, digits=3), postString), collapse = '')
    # Append query to vector
    queries = append(queries, queryString)
    min = min + increment
  }
  
  # Return queries
  return(queries)
}

get_file_names = function(min=0, finish=10, increment=0.5){ 
  # Create empty vector to store file names
  fileNames = c()
  
  while (round(min+increment, digits=3) <= finish) {
    # Create curFileName
    curFileName = paste(c(round(min, digits=3), '_', FILTER_PARAM, '_Default_cmg.tsv'), collapse = '')
    # Append curFileName to vector
    fileNames = append(fileNames, curFileName)
    min = min + increment
  }
  # Return file names
  return(fileNames)
}

get_greaterThanOrEqualTo_vals = function(min=0, finish=10, increment=0.5){
  # Create empty vector to store values
  vals = c()
  
  while (round(min+increment, digits=3) <= finish) {
    # Append min to vector
    vals = append(vals, round(min, digits=3))
    min = min + increment
  }
  # Return vals
  return(vals)
}

get_lessThan_vals = function(min=0, finish=10, increment=0.5){
  # Create empty vector to store values
  vals = c()
  
  while (round(min+increment, digits=3) <= finish) {
    # Append max to vector
    vals = append(vals, round(min + increment, digits=3))
    min = min + increment
  }
  # Return vals
  return(vals)
}

# Build query strings
queries = get_binned_query_series(MIN, FINISH, INCREMENT)

# Build file names
fileNames = get_file_names(MIN, FINISH, INCREMENT)

# Get min vals
minVals = get_greaterThanOrEqualTo_vals(MIN, FINISH, INCREMENT)

# Get max vals
maxVals = get_lessThan_vals(MIN, FINISH, INCREMENT)

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
    
    # Add extra column to current query (min query val)
    curResult$">=" = minVals[count]
 
    # Append results to master data frame
    allResults = rbind(allResults, curResult)
  }
  
  # Write current query to unique file
  curFileName = fileNames[count]
  write.table(curResult, file=paste(c(SUBDIR_NAME, "/", curFileName), collapse = ''), sep="\t", quote=FALSE, row.names=FALSE)
  
  count = count + 1
}

# Write master data frame to file
write.table(allResults, file=paste(c(as.character(Sys.Date()), "_", FILTER_PARAM, "_Default_", MIN, "-", FINISH, "_increment_", INCREMENT, "_cmg", ".tsv"), collapse = ''), sep="\t", quote=FALSE, row.names=FALSE)

