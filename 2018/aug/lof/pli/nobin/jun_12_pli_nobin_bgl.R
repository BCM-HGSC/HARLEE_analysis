# R Script
# Author: Adam W Hansen
# Date Created: Jan 26, 2018
# Date Last Modified: Jun 12, 2018

library(getPass)
library(RJDBC)

#USERNAME = readline(prompt = "Username: ")
USERNAME = getPass(msg = "Username: ")
PASSWORD = getPass(msg = "Password: ")
FILTER_PARAM = "exac_pli"

get_query = function(){
  # Generate text for query
  query  = paste(c("select count(distinct(sname)) as n_samples, gene_symbol, omim_gene_2016, omim_gene_2018, ", FILTER_PARAM, " 
    from freeze2_3.variant_genotype_filtered_qualpass
    left join freeze2_3.filtered_qualpass_maf
                on variant_genotype_filtered_qualpass.vid = filtered_qualpass_maf.vid
    where impact = 'HIGH'
    and f2_3_af < 0.01
    and (af < 0.001 or af is null)
    and (charge_af <= 0.0001 or charge_af is null)
    and (gnomad_af <= 0.0001 or gnomad_ext_af <= 0.0001 or (gnomad_af is null and gnomad_ext_af is null))
    and (pubmed = '' or pubmed is null)
    and (project not in ('CMG', 'SASCHX', 'WGLTEST', 'WGLM', 'WGLBP', 'WGLCARDIO', 'WGLC', 'WGLHES'))
    and (existing_variation = '' or existing_variation is null)
    and variant_genotype_filtered_qualpass.vid not like 'GL%'
    and (domains not like 'Low_complexity%' or domains is null)
    and exac_mu_syn is not null
    and ", FILTER_PARAM ," is not null
    group by gene_symbol, omim_gene_2016, omim_gene_2018, ", FILTER_PARAM, " 
    having count(distinct(sname)) >= 5
    order by ", FILTER_PARAM,  " desc"), collapse="")
  return(query)
}

get_file_name = function(){ 
  # Return file name
  return(paste(c('nobin_', FILTER_PARAM, '_Default.tsv'), collapse = ''))
}

# Build query 
query = get_query()

# Build file name
fileName = get_file_name()

# Make database connection
drv <- JDBC(driverClass = "com.cloudera.impala.jdbc41.Driver", classPath = list.files("/opt/JDBC/jars/ImpalaJDBC41/",pattern="jar$",full.names=T),identifier.quote="`")
impalaConnectionUrl <- paste("jdbc:impala://XXXX;AuthMech=4;SSLTrustStore=/opt/cloudera/security/jks/truststore.jks;SSLTrustStorePwd=XXXX;ssl=1;UID=", USERNAME, "@XXXX;PWD=", PASSWORD, sep="")
conn <- dbConnect(drv, impalaConnectionUrl)

# Initialize master data frame
allResults = NULL

# Execute queries
result = dbGetQuery(conn, query)
# Rename default col name "EXPR_0" as "n_samples
names(result)[names(result)=="EXPR_0"] = "n_samples"
  
# Write result data frame to file
write.table(result, file=paste(c(as.character(Sys.Date()), "_", FILTER_PARAM, "_Default_nobin_bgl.tsv"), collapse = ''), sep="\t", quote=FALSE, row.names=FALSE)

