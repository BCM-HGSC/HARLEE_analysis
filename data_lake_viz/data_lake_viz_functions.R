library(tidyverse)
library(gridExtra)

# Publication-quality ggplot2 theme found online at `https://rpubs.com/Koundy/71792`

theme_Publication <- function(base_size=14, base_family="Helvetica") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(), 
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               #legend.position = "none",
               legend.direction = "horizontal",
               legend.key.size= unit(0.2, "cm"),
               legend.margin = margin(0),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
          ))
      
}

theme_Publication_legendright <- function(base_size=14, base_family="Helvetica") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(), 
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
          ))
      
}

theme_Publication_nolegend <- function(base_size=14, base_family="Helvetica") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(), 
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "none",
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
          ))
      
}

scale_fill_Publication <- function(...){
      library(scales)
      discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

scale_colour_Publication <- function(...){
      library(scales)
      discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

visualize <- function(data, omim_year, plot_title, score_name, x_limits=c(-0.05,1), x_breaks=seq(0,1,0.1), log=TRUE, only16_18=FALSE){
    colnames(data)[5] <- "greater_or_equal"
    colnames(data)[6] <- "less"
    data = format_for_viz(data, only16_18=only16_18)
    
    if (omim_year == 2016){
        column_to_plot = data$n_2016
    }
    if (omim_year == 2018){
        column_to_plot = data$n_2018
    }
    
    if (log==TRUE){
        ggplot(data, aes(greater_or_equal, column_to_plot, fill = has_omim_phenotype)) +
            geom_bar(stat="identity", position = "dodge") +
            geom_line(aes(color = has_omim_phenotype)) +
            scale_fill_brewer(palette = "Paired") +
            ylab("# Genes") +
            scale_y_log10() +
            xlab(paste(score_name, " (binned)")) +
            scale_x_continuous(breaks=x_breaks, limits=x_limits) +
            ggtitle(plot_title)
    }
    else {
        ggplot(data, aes(greater_or_equal, column_to_plot, fill = has_omim_phenotype)) +
            geom_bar(stat="identity", position = "dodge") +
            geom_line(aes(color = has_omim_phenotype)) +
            scale_fill_brewer(palette = "Paired") +
            ylab("# Genes") +
            xlab(paste(score_name, " (binned)")) +
            scale_x_continuous(breaks=x_breaks, limits=x_limits) +
            ggtitle(plot_title)
    }
}

calculate_if_discovery <- function(df, col_past_name, col_future_name, newcol_name){
    df[,newcol_name] = (is.na(df[,col_past_name]) & !is.na(df[,col_future_name]) )
    return(df)
}

format_for_viz <- function(data, binned=TRUE, simulate=FALSE, n_sims=3, simyear=2018, downsample_fraction=0.8, only16_18=FALSE){
    ## uncomment to Filter out any genes where OMIM was "corrected", or included a phenotype in 2016, but no phenotype for 2018. Rationale: don't know how to handle those, I suspect many are actually true disease genes.
    #data = filter(data, !(!is.na(omim_gene_2016) & is.na(omim_gene_2018)))
    if (binned==TRUE){
        if(only16_18==FALSE){
            grouped_data_2013 = data %>% count(greater_or_equal, !is.na(omim_gene_2013))
            colnames(grouped_data_2013)[2] = "has_omim_phenotype"
            colnames(grouped_data_2013)[3] = "n_2013"

            grouped_data_2014 = data %>% count(greater_or_equal, !is.na(omim_gene_2014))
            colnames(grouped_data_2014)[2] = "has_omim_phenotype"
            colnames(grouped_data_2014)[3] = "n_2014"
        }
        
        grouped_data_2016 = data %>% count(greater_or_equal, !is.na(omim_gene_2016))
        colnames(grouped_data_2016)[2] = "has_omim_phenotype"
        colnames(grouped_data_2016)[3] = "n_2016"

        grouped_data_2018 = data %>% count(greater_or_equal, !is.na(omim_gene_2018))
        colnames(grouped_data_2018)[2] = "has_omim_phenotype"
        colnames(grouped_data_2018)[3] = "n_2018"
        
        if(only16_18==FALSE){
            grouped_data_discoveries_2013to2014 = data %>% count(greater_or_equal, discovery_2013to2014==TRUE)
            colnames(grouped_data_discoveries_2013to2014)[2] = "discovery"
            colnames(grouped_data_discoveries_2013to2014)[3] = "n_discoveries_2013to2014"
            grouped_data_discoveries_2013to2014[is.na("n_discoveries_2013to2014"),"n_discoveries_2013to2014"] = 0
            grouped_data_discoveries_2013to2014 = filter(grouped_data_discoveries_2013to2014, discovery==TRUE)[,c("greater_or_equal", "n_discoveries_2013to2014")]

            grouped_data_discoveries_2013to2016 = data %>% count(greater_or_equal, discovery_2013to2016==TRUE)
            colnames(grouped_data_discoveries_2013to2016)[2] = "discovery"
            colnames(grouped_data_discoveries_2013to2016)[3] = "n_discoveries_2013to2016"
            grouped_data_discoveries_2013to2016[is.na("n_discoveries_2013to2016"),"n_discoveries_2013to2016"] = 0
            grouped_data_discoveries_2013to2016 = filter(grouped_data_discoveries_2013to2016, discovery==TRUE)[,c("greater_or_equal", "n_discoveries_2013to2016")]

            grouped_data_discoveries_2013to2018 = data %>% count(greater_or_equal, discovery_2013to2018==TRUE)
            colnames(grouped_data_discoveries_2013to2018)[2] = "discovery"
            colnames(grouped_data_discoveries_2013to2018)[3] = "n_discoveries_2013to2018"
            grouped_data_discoveries_2013to2018[is.na("n_discoveries_2013to2018"),"n_discoveries_2013to2018"] = 0
            grouped_data_discoveries_2013to2018 = filter(grouped_data_discoveries_2013to2018, discovery==TRUE)[,c("greater_or_equal", "n_discoveries_2013to2018")]

            grouped_data_discoveries_2014to2016 = data %>% count(greater_or_equal, discovery_2014to2016==TRUE)
            colnames(grouped_data_discoveries_2014to2016)[2] = "discovery"
            colnames(grouped_data_discoveries_2014to2016)[3] = "n_discoveries_2014to2016"
            grouped_data_discoveries_2014to2016[is.na("n_discoveries_2014to2016"),"n_discoveries_2014to2016"] = 0
            grouped_data_discoveries_2014to2016 = filter(grouped_data_discoveries_2014to2016, discovery==TRUE)[,c("greater_or_equal", "n_discoveries_2014to2016")]

            grouped_data_discoveries_2014to2018 = data %>% count(greater_or_equal, discovery_2014to2018==TRUE)
            colnames(grouped_data_discoveries_2014to2018)[2] = "discovery"
            colnames(grouped_data_discoveries_2014to2018)[3] = "n_discoveries_2014to2018"
            grouped_data_discoveries_2014to2018[is.na("n_discoveries_2014to2018"),"n_discoveries_2014to2018"] = 0
            grouped_data_discoveries_2014to2018 = filter(grouped_data_discoveries_2014to2018, discovery==TRUE)[,c("greater_or_equal", "n_discoveries_2014to2018")]
        }
        
        grouped_data_discoveries_2016to2018 = data %>% count(greater_or_equal, discovery_2016to2018==TRUE)
        colnames(grouped_data_discoveries_2016to2018)[2] = "discovery"
        colnames(grouped_data_discoveries_2016to2018)[3] = "n_discoveries_2016to2018"
        grouped_data_discoveries_2016to2018[is.na("n_discoveries_2016to2018"),"n_discoveries_2016to2018"] = 0
        grouped_data_discoveries_2016to2018 = filter(grouped_data_discoveries_2016to2018, discovery==TRUE)[,c("greater_or_equal", "n_discoveries_2016to2018")]
        
        if(only16_18==TRUE){
             grouped_data = grouped_data_2016 %>% left_join(grouped_data_2018, by=c("greater_or_equal", "has_omim_phenotype")) %>% left_join(grouped_data_discoveries_2016to2018, by="greater_or_equal")
        }
        else{
            grouped_data = grouped_data_2013 %>% left_join(grouped_data_2014, by=c("greater_or_equal", "has_omim_phenotype")) %>% left_join(grouped_data_2016, by=c("greater_or_equal", "has_omim_phenotype")) %>% left_join(grouped_data_2018, by=c("greater_or_equal", "has_omim_phenotype")) %>% left_join(grouped_data_discoveries_2013to2014, by="greater_or_equal") %>% left_join(grouped_data_discoveries_2013to2016, by="greater_or_equal") %>% left_join(grouped_data_discoveries_2013to2018, by="greater_or_equal") %>% left_join(grouped_data_discoveries_2014to2016, by="greater_or_equal") %>% left_join(grouped_data_discoveries_2014to2018, by="greater_or_equal") %>% left_join(grouped_data_discoveries_2016to2018, by="greater_or_equal")
            grouped_data$n_discoveries_2013to2014[is.na(grouped_data$n_discoveries_2013to2014)] <- 0
            grouped_data$n_discoveries_2013to2016[is.na(grouped_data$n_discoveries_2013to2016)] <- 0
            grouped_data$n_discoveries_2013to2018[is.na(grouped_data$n_discoveries_2013to2018)] <- 0
            grouped_data$n_discoveries_2014to2016[is.na(grouped_data$n_discoveries_2014to2016)] <- 0
            grouped_data$n_discoveries_2014to2018[is.na(grouped_data$n_discoveries_2014to2018)] <- 0
            grouped_data$n_discoveries_2016to2018[is.na(grouped_data$n_discoveries_2016to2018)] <- 0
            grouped_data_nonomim = grouped_data %>% filter(has_omim_phenotype==FALSE) %>% select(greater_or_equal, n_2013, n_2014, n_2016, n_2018, n_discoveries_2013to2014, n_discoveries_2013to2016, n_discoveries_2013to2018, n_discoveries_2014to2016,  n_discoveries_2014to2018, n_discoveries_2016to2018)

            grouped_data_nonomim$discovery_density_2013to2014 = grouped_data_nonomim$n_discoveries_2013to2014/grouped_data_nonomim$n_2014
            grouped_data_nonomim$discovery_density_2013to2016 = grouped_data_nonomim$n_discoveries_2013to2016/grouped_data_nonomim$n_2016
            grouped_data_nonomim$discovery_density_2013to2018 = grouped_data_nonomim$n_discoveries_2013to2018/grouped_data_nonomim$n_2018
            grouped_data_nonomim$discovery_density_2014to2016 = grouped_data_nonomim$n_discoveries_2014to2016/grouped_data_nonomim$n_2016
            grouped_data_nonomim$discovery_density_2014to2018 = grouped_data_nonomim$n_discoveries_2014to2018/grouped_data_nonomim$n_2018
            grouped_data_nonomim$discovery_density_2016to2018 = grouped_data_nonomim$n_discoveries_2016to2018/grouped_data_nonomim$n_2018
            grouped_data = left_join(grouped_data, (grouped_data_nonomim %>% select(greater_or_equal, discovery_density_2013to2014, discovery_density_2013to2016, discovery_density_2013to2018, discovery_density_2014to2016, discovery_density_2014to2018, discovery_density_2016to2018)), by="greater_or_equal")
            grouped_data = left_join(grouped_data, grouped_data %>% group_by(greater_or_equal) %>% summarize(total_genelist_size = sum(n_2016)), by="greater_or_equal")
        }
        
        
        if (simulate==TRUE){
            # Read OMIM table in
            if (simyear==2013){
                omim_raw = read.table("omim/omim_flat_02_07_2013_noheader.txt", header=FALSE, sep="|", quote = "", comment.char = "")
            }
            if (simyear==2014){
                omim_raw = read.table("omim/omim_flat_12_09_2014_noheader.txt", header=FALSE, sep="|", quote = "", comment.char = "")
            }
            if (simyear==2016){
                omim_raw = read.table("omim/omim_flat_08_18_2016_noheader.txt", header=FALSE, sep="\t", quote = "", comment.char = "")
            }
            if (simyear==2018){
                omim_raw = read.table("omim/omim_flat_01_08_2018_noheader.txt", header=FALSE, sep="\t", quote= "", comment.char = "")
            }
            colnames(omim_raw) = c("phenotype", "all_genes", "omim_id", "cyto_band", "gene_id")
            # For i in n_sims:
            set.seed(2)
            for (i in 1:n_sims){
                # create string paste('sim', i, sep=''), as col name prefix
                simstring = paste('sim', i, sep='')
                # Subsample 90% of OMIM rows
                omim_sub = sample_frac(omim_raw[,c("gene_id", "phenotype")], downsample_fraction )
                # Group OMIM by gene, concat OMIM annotations
                omim_sub_gr = omim_sub %>% group_by(gene_id) %>% summarise(test = toString(phenotype)) %>% ungroup()
                # Annotate gene hits with this OMIM info
                data$gene_symbol = as.character(data$gene_symbol)
                omim_sub_gr$gene_id = as.character(omim_sub_gr$gene_id)
                sub_data = left_join(data, omim_sub_gr, by=c("gene_symbol"="gene_id"))
                # Calculate "discovery" column
                sub_data$newcol = (is.na(sub_data$test) & !is.na(sub_data$omim_gene_2016))
                
                # Group as above, counting OMIM hits vs non-hits per gene bin; append to master data frame
                sub_grouped_data = sub_data %>% count(greater_or_equal, !is.na(test))
                colnames(sub_grouped_data)[2] = "has_omim_phenotype"
                colnames(sub_grouped_data)[3] = "n_sub"
                grouped_data = left_join(grouped_data, sub_grouped_data, by=c("greater_or_equal", "has_omim_phenotype"))
                colnames(grouped_data)[colnames(grouped_data) == "n_sub"] = paste(simstring, ".omim_gene", sep='')
                
                # Calculate n_discoveries as above; append to master data frame
                sub_grouped_data_discoveries = sub_data %>% count(greater_or_equal, newcol==TRUE)
                colnames(sub_grouped_data_discoveries)[2] = "discovery"
                colnames(sub_grouped_data_discoveries)[3] = "n_discoveries_sub"
                sub_grouped_data_discoveries[is.na("n_discoveries_sub"),"n_discoveries_sub"] = 0
                sub_grouped_data_discoveries = filter(sub_grouped_data_discoveries, discovery==TRUE)[,c("greater_or_equal", "n_discoveries_sub")]
                grouped_data = left_join(grouped_data, sub_grouped_data_discoveries, by="greater_or_equal")
                grouped_data$n_discoveries_sub[is.na(grouped_data$n_discoveries_sub)] = 0
                colnames(grouped_data)[colnames(grouped_data) =='n_discoveries_sub'] = paste(simstring, ".discovery" , sep='')
                
                # Calculate discovery density as above; append to master data frame
                sub_grouped_data_nonomim = sub_grouped_data %>% filter(has_omim_phenotype==FALSE) %>% left_join(sub_grouped_data_discoveries, by="greater_or_equal") %>% select(greater_or_equal, n_sub, n_discoveries_sub)
                sub_grouped_data_nonomim$discovery_density_sub = sub_grouped_data_nonomim$n_discoveries_sub/sub_grouped_data_nonomim$n_sub
                sub_grouped_data_discovery_density = select(sub_grouped_data_nonomim, c('greater_or_equal', 'discovery_density_sub'))
                grouped_data = left_join(grouped_data, sub_grouped_data_discovery_density, by="greater_or_equal")
                grouped_data$discovery_density_sub[is.na(grouped_data$discovery_density_sub)] = 0
                colnames(grouped_data)[colnames(grouped_data) =='discovery_density_sub'] = paste(simstring, ".discovery_density" , sep='')
            }
        }
    }
    if (binned==FALSE) {
        grouped_data_2013 = data %>% count(group_index, !is.na(omim_gene_2013))
        colnames(grouped_data_2013)[2] = "has_omim_phenotype"
        colnames(grouped_data_2013)[3] = "n_2013"
        
        grouped_data_2014 = data %>% count(group_index, !is.na(omim_gene_2014))
        colnames(grouped_data_2014)[2] = "has_omim_phenotype"
        colnames(grouped_data_2014)[3] = "n_2014"
        
        grouped_data_2016 = data %>% count(group_index, !is.na(omim_gene_2016))
        colnames(grouped_data_2016)[2] = "has_omim_phenotype"
        colnames(grouped_data_2016)[3] = "n_2016"

        grouped_data_2018 = data %>% count(group_index, !is.na(omim_gene_2018))
        colnames(grouped_data_2018)[2] = "has_omim_phenotype"
        colnames(grouped_data_2018)[3] = "n_2018"

        #grouped_data = dplyr::left_join(grouped_data_2016, grouped_data_2018, by=c("group_index", "has_omim_phenotype"))
        grouped_data = grouped_data_2013 %>% left_join(grouped_data_2014, by=c("group_index", "has_omim_phenotype")) %>% left_join(grouped_data_2016, by=c("group_index", "has_omim_phenotype")) %>% left_join(grouped_data_2018, by=c("group_index", "has_omim_phenotype"))
        grouped_data[is.na(grouped_data)] <- 0
    }
    return (grouped_data)
}

cutoff_calculation_prepare_df = function(nobin, reverse=FALSE, log_scale=FALSE){    
    if (reverse==TRUE){
        nobin <- nobin[nrow(nobin):1,]
        rownames(nobin) <- NULL
    }
    for (i in 1:nrow(nobin)){
        nobin$cutoff_gene_list_size_2018[i] <- sum(is.na(nobin$omim_gene_2018[1:i]))
        nobin$cutoff_discovery[i] <- sum(nobin$discovery_2013to2018[1:i])
        nobin$cutoff_discovery_density[i] <- nobin$cutoff_discovery[i] / nobin$cutoff_gene_list_size_2018[i]
    }
    return(nobin)
}

calculate_cutoff_nobin = function(nobin, min_gene_list_size=20){
    nobin_filtered <- filter(nobin, cutoff_gene_list_size_2018 >= min_gene_list_size)
    nobin_cutoff <- vector(mode="list", length=4)
    names(nobin_cutoff) <- c("density", "index", "value", "candidate_list_size_2018")
    nobin_cutoff$density = max(nobin_filtered$cutoff_discovery_density)
    nobin_cutoff$index = as.numeric(row.names(nobin[nobin$cutoff_discovery_density==nobin_cutoff$density & nobin$cutoff_gene_list_size_2018 >= min_gene_list_size,]))[1]
    nobin_cutoff$value = nobin[nobin_cutoff$index,5]
    nobin_cutoff$candidate_list_size_2018 = nobin[nobin_cutoff$index,"cutoff_gene_list_size_2018"]
    return(nobin_cutoff)
}

calculate_cutoff_filtered = function(grouped, min_gene_list_size=20, sim=FALSE, ddyearstart=2013, ddyearend=2018){
    grouped_non_omim_filtered = filter(grouped, n_2018 >= min_gene_list_size, has_omim_phenotype == FALSE)
    grouped_omim_filtered = semi_join(x=filter(grouped, has_omim_phenotype == TRUE), y=grouped_non_omim_filtered, by="greater_or_equal")
    grouped_cutoff <- vector(mode="list", length=4)
    if (sim==FALSE){
        names(grouped_cutoff) <- c("density", "value", "candidate_list_size", "index")
        if (ddyearend==2014){
            if (ddyearstart==2013){
                grouped_cutoff$density = max(grouped_omim_filtered$discovery_density_2013to2014)
                grouped_cutoff$value = grouped_omim_filtered[which(grouped_omim_filtered$discovery_density_2013to2014==grouped_cutoff$density)[1],1]
            }
            grouped_cutoff$candidate_list_size = grouped_non_omim_filtered[which(grouped_non_omim_filtered$greater_or_equal==grouped_cutoff$value),'n_2014']
        }
        if (ddyearend==2016){
            if (ddyearstart==2013){
                grouped_cutoff$density = max(grouped_omim_filtered$discovery_density_2013to2016)
                grouped_cutoff$value = grouped_omim_filtered[which(grouped_omim_filtered$discovery_density_2013to2016==grouped_cutoff$density)[1],1]
            }
            if (ddyearstart==2014){
                grouped_cutoff$density = max(grouped_omim_filtered$discovery_density_2014to2016)
                grouped_cutoff$value = grouped_omim_filtered[which(grouped_omim_filtered$discovery_density_2014to2016==grouped_cutoff$density)[1],1]
            }
            grouped_cutoff$candidate_list_size = grouped_non_omim_filtered[which(grouped_non_omim_filtered$greater_or_equal==grouped_cutoff$value),'n_2016']
        }
        if (ddyearend==2018){
            if (ddyearstart==2013){
                grouped_cutoff$density = max(grouped_omim_filtered$discovery_density_2013to2018)
                grouped_cutoff$value = grouped_omim_filtered[which(grouped_omim_filtered$discovery_density_2013to2018==grouped_cutoff$density)[1],1]
            }
            if (ddyearstart==2014){
                grouped_cutoff$density = max(grouped_omim_filtered$discovery_density_2014to2018)
                grouped_cutoff$value = grouped_omim_filtered[which(grouped_omim_filtered$discovery_density_2014to2018==grouped_cutoff$density)[1],1]
            }
            if (ddyearstart==2016){
                grouped_cutoff$density = max(grouped_omim_filtered$discovery_density_2016to2018)
                grouped_cutoff$value = grouped_omim_filtered[which(grouped_omim_filtered$discovery_density_2016to2018==grouped_cutoff$density)[1],1]
            }
            grouped_cutoff$candidate_list_size = grouped_non_omim_filtered[which(grouped_non_omim_filtered$greater_or_equal==grouped_cutoff$value),'n_2018']
        }
        grouped_cutoff$index = which(grouped$greater_or_equal == grouped_cutoff$value)[1]
    }
    if (sim==TRUE){
        names(grouped_cutoff) <- c("density", "value", "candidate_list_size", "index")
        grouped_cutoff$density = max(grouped_omim_filtered$discovery_density)
        grouped_cutoff$value = grouped_omim_filtered[which(grouped_omim_filtered$discovery_density==grouped_cutoff$density)[1],1]
        grouped_cutoff$candidate_list_size = grouped_non_omim_filtered[which(grouped_non_omim_filtered$greater_or_equal==grouped_cutoff$value),'n_2018']
        grouped_cutoff$index = which(grouped$greater_or_equal == grouped_cutoff$value)[1]
    }
    return(grouped_cutoff)
}

get_sim_cutoffs = function(grouped, N_SIMS=3){
    sim_cutoffs = vector(length = N_SIMS)
    for (i in 1:N_SIMS){
        simstring = paste0("sim", as.character(i), ".")
        n_string = paste0(simstring, 'omim_gene')
        discovery_string = paste0(simstring, 'discovery')
        discovery_density_string = paste0(simstring, 'discovery_density', sep='')
        sub_df = grouped[,c('greater_or_equal', 'has_omim_phenotype', n_string, 'n_2018', discovery_string, discovery_density_string, 'total_genelist_size')]
        colnames(sub_df) = c("greater_or_equal", "has_omim_phenotype", "omim_gene", "n_2018", "discovery", "discovery_density", "total_genelist_size")
        sub_cutoff = calculate_cutoff_filtered(sub_df, sim=TRUE)
        # Append sub cutoff to set of cutoffs
        sim_cutoffs[i] = sub_cutoff$value
    }
    return(sim_cutoffs)
}

get_real_cutoffs = function(grouped){
    real_cutoffs = vector(length = 6)
    names(real_cutoffs) <- c("from2013to2014", "from2013to2016", "from2013to2018", "from2014to2016", "from2014to2018", "from2016to2018")
    real_cutoffs$from2013to2014 = calculate_cutoff_filtered(grouped, ddyearstart=2013, ddyearend=2014)$value
    real_cutoffs$from2013to2016 = calculate_cutoff_filtered(grouped, ddyearstart=2013, ddyearend=2016)$value
    real_cutoffs$from2013to2018 = calculate_cutoff_filtered(grouped, ddyearstart=2013, ddyearend=2018)$value
    real_cutoffs$from2014to2016 = calculate_cutoff_filtered(grouped, ddyearstart=2014, ddyearend=2016)$value
    real_cutoffs$from2014to2018 = calculate_cutoff_filtered(grouped, ddyearstart=2014, ddyearend=2018)$value
    real_cutoffs$from2016to2018 = calculate_cutoff_filtered(grouped, ddyearstart=2016, ddyearend=2018)$value
    return(real_cutoffs)
}

import_nobin_data = function(fileString="", binSize=40, reverse=FALSE, prepare_for_cutoff_calculation=TRUE, calc_index=TRUE){
    # import
    nobin = read.table(fileString, header=TRUE, sep="\t", quote = "", comment.char = "")
    
    # Import 2014 OMIM, annotate, calculate discovery bool
    omim_raw_2014 = read.table("omim/omim_flat_12_09_2014_noheader.txt", header=FALSE, sep="|", quote = "", comment.char = "")
    colnames(omim_raw_2014) = c("omim_gene_2014_raw", "genes", "omim_id", "cyto_band", "gene_symbol")
    omim_raw_2014$gene_symbol = gsub(",", "",omim_raw_2014$gene_symbol)
    omim_2014 = select(omim_raw_2014, gene_symbol, omim_gene_2014_raw) %>% group_by(gene_symbol) %>% summarise(omim_gene_2014=paste(omim_gene_2014_raw, collapse=","))
    nobin = nobin %>% left_join(omim_2014)
    
    # Import 2013 OMIM, annotate, calculate discovery bool
    omim_raw_2013 = read.table("omim/omim_flat_02_07_2013_noheader.txt", header=FALSE, sep="|", quote = "", comment.char = "")
    colnames(omim_raw_2013) = c("omim_gene_2013_raw", "genes", "omim_id", "cyto_band", "gene_symbol")
    omim_raw_2013$gene_symbol = gsub(",", "",omim_raw_2013$gene_symbol)
    omim_2013 = select(omim_raw_2013, c("gene_symbol", "omim_gene_2013_raw")) %>% group_by(gene_symbol) %>% summarise(omim_gene_2013=paste(omim_gene_2013_raw, collapse=","))
    nobin = nobin %>% left_join(omim_2013)
    
    # calculate discovery bools
    nobin$discovery_2013to2014 = (is.na(nobin$omim_gene_2013) & !is.na(nobin$omim_gene_2014))
    nobin$discovery_2013to2016 = (is.na(nobin$omim_gene_2013) & !is.na(nobin$omim_gene_2016))
    nobin$discovery_2013to2018 = (is.na(nobin$omim_gene_2013) & !is.na(nobin$omim_gene_2018))
    nobin$discovery_2014to2016 = (is.na(nobin$omim_gene_2014) & !is.na(nobin$omim_gene_2016))
    nobin$discovery_2014to2018 = (is.na(nobin$omim_gene_2014) & !is.na(nobin$omim_gene_2018))
    nobin$discovery_2016to2018 = (is.na(nobin$omim_gene_2016) & !is.na(nobin$omim_gene_2018))
    
    if (calc_index==TRUE){
        # add group index
        nobin$group_index = ceiling(as.numeric(row.names(nobin))/binSize)
    }
    
    if (prepare_for_cutoff_calculation==TRUE){
        nobin = cutoff_calculation_prepare_df(nobin, reverse)
    }
    return(nobin)
}

group_nobin_data = function(nobin, binSize=40, plot=FALSE){
    # group
    nobinGrouped = format_for_viz(data=nobin, binned=FALSE)
    
    # remove "remainder" gene list with less than $binSize genes
    nobinGrouped = nobinGrouped[1:(nrow(nobinGrouped)-2),]
    
    # reverse order of group_index, such that lower indicies indicate a lower score
    rownames(nobinGrouped) <- NULL
    nobinGrouped=nobinGrouped[order(nrow(nobinGrouped):1),]
    rownames(nobinGrouped) <- NULL
    
    # Recalculate/assign group index
    nobinGrouped$group_index = ceiling(as.numeric(rownames(nobinGrouped))/2)
    
    if(plot==TRUE){
        ggplot(nobinGrouped, aes(group_index, n_2016, color = has_omim_phenotype)) + geom_point() + scale_fill_brewer(palette = "Paired") + ylab("# Genes") + xlab(paste(c("Group Index (bin size = ", as.character(binSize), " genes)"), sep=""))
    }
    
    # filter to only include TRUE count + ratios
    nobinGroupedRatios = filter(nobinGrouped, has_omim_phenotype==TRUE)
    nobinGroupedRatios$discoveries_2013to2014 = nobinGroupedRatios$n_2014 - nobinGroupedRatios$n_2013
    nobinGroupedRatios$discoveries_2013to2016 = nobinGroupedRatios$n_2016 - nobinGroupedRatios$n_2013
    nobinGroupedRatios$discoveries_2013to2018 = nobinGroupedRatios$n_2018 - nobinGroupedRatios$n_2013
    nobinGroupedRatios$discoveries_2014to2016 = nobinGroupedRatios$n_2016 - nobinGroupedRatios$n_2014
    nobinGroupedRatios$discoveries_2014to2018 = nobinGroupedRatios$n_2018 - nobinGroupedRatios$n_2014
    nobinGroupedRatios$discoveries_2016to2018 = nobinGroupedRatios$n_2018 - nobinGroupedRatios$n_2016
    nobinGroupedRatios$ratio_2013 = nobinGroupedRatios$n_2013/binSize
    nobinGroupedRatios$ratio_2014 = nobinGroupedRatios$n_2014/binSize
    nobinGroupedRatios$ratio_2016 = nobinGroupedRatios$n_2016/binSize
    nobinGroupedRatios$ratio_2018 = nobinGroupedRatios$n_2018/binSize
    
    return(nobinGroupedRatios)
}

fit_score_by_ratio_models = function(data){
    model.lo.2013 = loess(ratio_2013 ~ group_index, data, span=1)
    data$residual_2013_loess = model.lo.2013$residuals
    data$fit_2013_loess = model.lo.2013$fit

    model.lm.2013 = lm(ratio_2013 ~ group_index, data)
    data$residual_2013_lm = model.lm.2013$residuals
    data$fit_2013_lm = model.lm.2013$fit

    model.lo.2014 = loess(ratio_2014 ~ group_index, data, span=1)
    data$residual_2014_loess = model.lo.2014$residuals
    data$fit_2014_loess = model.lo.2014$fit

    model.lm.2014 = lm(ratio_2014 ~ group_index, data)
    data$residual_2014_lm = model.lm.2014$residuals
    data$fit_2014_lm = model.lm.2014$fit
    
    model.lo.2016 = loess(ratio_2016 ~ group_index, data, span=1)
    data$residual_2016_loess = model.lo.2016$residuals
    data$fit_2016_loess = model.lo.2016$fit

    model.lm.2016 = lm(ratio_2016 ~ group_index, data)
    data$residual_2016_lm = model.lm.2016$residuals
    data$fit_2016_lm = model.lm.2016$fit

    model.lo.2018 = loess(ratio_2018 ~ group_index, data, span=1)
    data$residual_2018_loess = model.lo.2018$residuals
    data$fit_2018_loess = model.lo.2018$fit

    model.lm.2018 = lm(ratio_2018 ~ group_index, data)
    data$residual_2018_lm = model.lm.2018$residuals
    data$fit_2018_lm = model.lm.2018$fit
    
    return(data)
}

# Multiple plot function (from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/)
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

plot_dd_year_combinations = function(grouped, scorename="greater_or_equal"){
    p0 = ggplot(filter(grouped, has_omim_phenotype==TRUE), aes(x=greater_or_equal)) + geom_line(aes(y=discovery_density_2013to2014), color="purple") + geom_line(aes(y=discovery_density_2014to2016), color="orange") + geom_line(aes(y=discovery_density_2016to2018), color="green") + xlab(scorename) + ylab("Discovery Density") + theme_Publication_nolegend()
    p1 = ggplot(filter(grouped, has_omim_phenotype==TRUE)) + geom_point(aes(x=discovery_density_2013to2014, y=discovery_density_2014to2016), alpha=0.1) + theme_Publication_nolegend() + xlab("DD 13-14") + ylab("DD 14-16")
    p2 = ggplot(filter(grouped, has_omim_phenotype==TRUE)) + geom_point(aes(x=discovery_density_2013to2014, y=discovery_density_2016to2018), alpha=0.1) + theme_Publication_nolegend() + xlab("DD 13-14") + ylab("DD 16-18")
    p3 = ggplot(filter(grouped, has_omim_phenotype==TRUE)) + geom_point(aes(x=discovery_density_2014to2016, y=discovery_density_2016to2018), alpha=0.1) + theme_Publication_nolegend() + xlab("DD 14-16") + ylab("DD 16-18")
    p4 = ggplot(filter(grouped, has_omim_phenotype==TRUE)) + geom_point(aes(x=discovery_density_2013to2014, y=discovery_density_2014to2018), alpha=0.1) + theme_Publication_nolegend() + xlab("DD 13-14") + ylab("DD 14-18")
    p5 = ggplot(filter(grouped, has_omim_phenotype==TRUE)) + geom_point(aes(x=discovery_density_2013to2016, y=discovery_density_2016to2018), alpha=0.1) + theme_Publication_nolegend() + xlab("DD 13-16") + ylab("DD 16-18")
    multiplot(p0, p1, p2, p3, p4, p5, cols=2)
}

get_cutoff_violin_plot = function(grouped, n_sims=3, xlabel, ylabel=""){
    grouped_real_cutoffs = get_real_cutoffs(grouped)
    grouped_sim_cutoffs = get_sim_cutoffs(grouped, n_sims)
    return(ggplot(data.frame(grouped_sim_cutoffs), aes(x=xlabel, y=grouped_sim_cutoffs)) + geom_violin() + xlab("") + ylab(ylabel) + geom_hline(yintercept = grouped_real_cutoffs$from2013to2014) + geom_hline(yintercept = grouped_real_cutoffs$from2013to2016) + geom_hline(yintercept = grouped_real_cutoffs$from2014to2016) + geom_hline(yintercept = grouped_real_cutoffs$from2014to2018) + geom_hline(yintercept = grouped_real_cutoffs$from2016to2018) + geom_hline(yintercept = grouped_real_cutoffs$from2013to2018, color="red") + theme_Publication())
}

plot_simulated_discovery = function(grouped, n_sims=3, score_plot_title){
    grouped = filter(grouped, has_omim_phenotype==FALSE)
    
    melted_df = data.frame(greater_or_equal=double(),
                 has_omim_phenotype=logical(),
                 omim_gene=character(),
                 n_2018_false=integer(),
                 discovery=integer(),
                 discovery_density=double(),
                 total_genelist_size=integer(),
                 sim_number=integer())
    
    for (i in 1:n_sims){
        cursim_discoveries = vector(length = n_sims)
        
        simstring = paste0("sim", as.character(i), ".")
        n_string = paste0(simstring, 'omim_gene')
        discovery_string = paste0(simstring, 'discovery')
        discovery_density_string = paste0(simstring, 'discovery_density', sep='')
        
        sub_df = grouped[,c('greater_or_equal', 'has_omim_phenotype', n_string, 'n_2018', discovery_string, discovery_density_string, 'total_genelist_size')]
        colnames(sub_df) = c("greater_or_equal", "has_omim_phenotype", "omim_gene", "n_2018_false", "discovery", "discovery_density", "total_genelist_size")
        sub_df$sim_number = i
        
        melted_df = rbind(melted_df, sub_df)

    }
    
    light_df = grouped %>% filter(n_discoveries_2013to2018!=0, !is.na(n_2018), !is.na(n_discoveries_2013to2018)) %>% select(greater_or_equal, n_discoveries_2013to2018, total_genelist_size, n_2018)
    colnames(light_df)[colnames(light_df)=="n_2018"] = 'n_2018_false'
    avg_df = melted_df %>% filter(discovery_density!=0) %>% group_by(greater_or_equal) %>% summarize(avg_sim_discovery_density = mean(discovery_density, na.rm=TRUE))
    avg_df = light_df %>% left_join(avg_df)
    melted_df = melted_df %>% left_join(light_df %>% select(greater_or_equal, n_discoveries_2013to2018)) %>% filter(discovery_density!=0)
    
    xval = (avg_df$avg_sim_discovery_density)
    yval = (avg_df$n_discoveries_2013to2018/(avg_df$n_2018_false)) 
    
    print(ggplot(melted_df, aes(x=discovery_density, y=n_discoveries_2013to2018/n_2018_false, color=greater_or_equal)) + xlab("Simulated Past Discovery Density (.8 x 2013)-2013") + ylab("Discovery Density 2013–2018") + geom_point(alpha=0.05) + theme_Publication_legendright() + labs(color = "Constraint\nCutoff"))
    
    equation_xpos = min(xval)+((range(xval)[2]-range(xval)[1])*.4)
    equation_ypos = min(yval)+((range(yval)[2]-range(yval)[1])*.925)

    print(ggplot(avg_df, aes(x=xval, y=yval, color=greater_or_equal)) + xlab(paste("Avg Simulated Past Discovery Density n=", as.character(N_SIMS), " (.8 x 2013)-2013", sep="")) + ylab("Discovery Density 2013–2018") + geom_point(alpha=0.1) + geom_smooth(method=lm, se=FALSE) + geom_text(x=equation_xpos, y=equation_ypos, label = lm_eqn(avg_df, x=xval, y=yval), parse = TRUE, color="black") + theme_Publication_legendright() + labs(color = "Constraint\nCutoff"))
    
    print(ggplot(melted_df, aes(x=greater_or_equal)) + geom_point(aes(y=melted_df$discovery_density), alpha=0.01) + xlab(score_plot_title) + ylab("Discovery Density 2013–2018 (red) or simulated (black)") + geom_line(data=grouped, aes(x=greater_or_equal, y=discovery_density_2013to2018), color="red") + theme_Publication())
}

# GET EQUATION AND R-SQUARED AS STRING
# SOURCE: http://goo.gl/K4yh

lm_eqn <- function(df, x, y){
    m <- lm(y ~ x, df)
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
         list(a = format(coef(m)[1], digits = 2), 
              b = format(coef(m)[2], digits = 2), 
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
}

import_process_and_simulate = function(inFileString, outFileString, n_sims){
    series_filter <- import_nobin_data(fileString=inFileString, prepare_for_cutoff_calculation = FALSE, calc_index = FALSE)
    colnames(series_filter)[5] = 'greater_or_equal'
    series_grouped <- format_for_viz(series_filter, simulate=TRUE, n_sims=n_sims)
    write.table(series_grouped, file=outFileString, sep="\t", row.names=FALSE)
}
