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

query =
    "select 
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
        and exac_pli >= 1.0
    ) cases
    left anti join 
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
        and (sname rlike '[Mm]other' or sname rlike '[Ff]ather' or sname rlike '[\\-_0-9][Mm]$' or sname rlike '[\\-_0-9][Ff]$' or sname rlike '[\\-_0-9][Dd]$' or sname rlike '[\\-_]0?[2-3]' or sname rlike '^[0-9]{4}[A-Z]{2}0[2-3]$' or sname like 'WPW%' or sname like 'HGSC%' or sname like 'NA%')
        and exac_pli >= 1.0
    ) controls
    on cases.vid=controls.vid
    group by gene_symbol, omim_gene_2016, omim_gene_2018
    having count(distinct(sname)) >= 5
    order by omim_gene_2016 nulls first, omim_gene_2018 nulls last, count(distinct(sname)) desc"
  
# Make database connection
drv <- JDBC(driverClass = "com.cloudera.impala.jdbc41.Driver", classPath = list.files("/opt/JDBC/jars/ImpalaJDBC41/",pattern="jar$",full.names=T),identifier.quote="`")
impalaConnectionUrl <- paste("jdbc:impala://XXXX;AuthMech=4;SSLTrustStore=/opt/cloudera/security/jks/truststore.jks;SSLTrustStorePwd=XXXX;ssl=1;UID=", USERNAME, "@XXXX;PWD=", PASSWORD, sep="")
conn <- dbConnect(drv, impalaConnectionUrl)

# Execute current query, get results
curResult = dbGetQuery(conn, query)
# Rename default col name "EXPR_0" as "n_samples
names(curResult)[names(curResult)=="EXPR_0"] = "n_samples"

# Add extra column to current query (min query val)
curResult$">=" = 1
 
# Write master data frame to file
write.table(curResult, file=paste(c(as.character(Sys.Date()), "_exac_pli_Default_1_fixed_cmg", ".tsv"), collapse = ''), sep="\t", quote=FALSE, row.names=FALSE)

