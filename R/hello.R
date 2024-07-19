#read in Env

#' Function to normalise column names and tidy up enviromental metadata
#' ready for use in PostGres DB
#' 
#' @param input_filename String containing filename of input dataset.  Expects 'csv'
#' @param output_filename String containing filename of data set to be written out as csv
#' 
#' @return None
#' 
#' @examples
#' get_clean_filter_enviromnetal_metadata('CS_envmeta.csv', 'CS_envmeta_filtered.csv')
#' 
#' @export
get_clean_filter_environmetal_metadata <- function(input_filename, output_filename){
     # only need 3 columns
    enviromental_metadata <- fread(input_filename,
                                select = c("avc_code", "avc", "pH"))

    #   # postgres requires unquoted names to be lower
    #   names(enviromental_metadata) <- base::tolower(names(enviromental_metadata))

    # perhaps faster?
    setnames(enviromental_metadata, 
            old = names(enviromental_metadata),
            new = base::tolower(names(enviromental_metadata))
            )

    # lose rows with no avc_code values
    na.omit(enviromental_metadata, cols="avc_code")

    # create a column named samples to hold row numbers and set as first column
    enviromental_metadata[, sample := row.names(enviromental_metadata)]
    setcolorder(enviromental_metadata, "sample", before=1)

    fwrite(enviromental_metadata, "output_filename")

}

#' NEW DESCRIPTION
#' 
#' @param input_filename String containing filename of input dataset.  Expects 'csv'
#' @param output_filename String containing filename of data set to be written out as csv
#' 
#' @return None
#' 
#' @examples
#' get_clean_filter_OTU_table('CS_envmeta.csv', 'CS_envmeta_filtered.csv')
#' 
#' @export
get_clean_filter_OTU_table <- function(input_filename, output_filename, sample_reads_limits=5000, occupancy_threshold=500){

    OTU_table <- fread(input_filename)
    # Remove samples with less than 5000 reads
    OTU_table <- OTU_table[rowSums(OTU_tab)>sample_reads_limits]
    
    column_names <- names(OTU_table)
    # set all values to logical boolean
    OTU_table[,(column_names) := lapply(.SD, as.logical), .SDcols = column_names]
    # convert booleans to 0 or 1
    OTU_table[,(column_names) := lapply(.SD, as.numeric), .SDcols = column_names]
    # remove columns where taxa do not meet occupancy threshold
    # https://github.com/Rdatatable/data.table/issues/795
    OTU_table[, OTU_table[, names(.SD), .SDcols = -(colSums(OTU_table>occupancy_threshold))]:=NULL]
    # normalise rows to sum to 1(what is the prc)
    OTU_table[, names(OTU_table) := .SD/rowSums(.SD)]
    fwrite(OTU_table, "output_filename")
}


# Workflow

## Step 1 Prepare tables for database
### Prepare  OTU tab and Env 

output_dir <- '/raid3/scratch/MEGshared/Fungi_Explorer_App/Fungi_Explorer_unpublished_app_experiments/Outputs'
input_dir <- '/raid3/scratch/MEGshared/Fungi_Explorer_App/Fungi_Explorer/Inputs/Data/'

OTU_input_file.csv <- 'basic_preprocessed/OTU_tab.csv'
OTU_output_file.csv <- '/occ_threshold_500/Supplementary/Tables_in_SQL/OTU_abund.csv'
enviromental_meta_input_file <- 'basic_preprocessed/Env.csv'
enviromental_meta_output_file <- '/occ_threshold_500//Supplementary/Tables_in_SQL/Env.csv'

get_clean_filter_enviromnetal_metadata(enviromental_meta_input_file, enviromental_meta_output_file)
get_clean_filter_OTU_table(OTU_input_file.csv, OTU_output_file.csv)





#' NEW DESCRIPTION
#' 
#' @param input_filename String containing filename of input dataset.  Expects 'csv'
#' @param output_filename String containing filename of data set to be written out as csv
#' 
#' @return None
#' 
#' @examples
#' get_clean_filter_OTU_table('CS_envmeta.csv', 'CS_envmeta_filtered.csv')
#' 
#' @export
get_clean_filter_OTU_table <- function(OTU_table_filename){

    OTU_table <- fread(OTU_table_filename)

    # Get the names of ASVs, and total the total abundance across the samples.
    ASV_names <- names(OTU_table)
    abundance <- colSums(OTU_table)
    sample_abundance <- data.table(hit=ASV_names, abundance=abundance)

    #get OTU abundance rank no.. where there are ties both values get the same rank i.e if two values that would be ranked 6 and 7 are the same they will both be ranked 6
    #contextualise the number by referencing the total amount of OTU/ASVs
    # Get rank as a fraction of ASVs, record as character
    abundance_stats[,abundance_rank:=paste(rank(abundance,ties.method="min"),ncol(OTU_table),sep='/'),]

    #get presence absence again for remaining taxa
    column_names <- names(OTU_table)
    # set all values to logical boolean
    OTU_table[,(column_names) := lapply(.SD, as.logical), .SDcols = column_names]
    # convert booleans to 0 or 1
    OTU_table[,(column_names) := lapply(.SD, as.numeric), .SDcols = column_names]












    fwrite(OTU_table, "output_filename")
}

### Prepare abundance_stats table

Get individual OTU stats to summarise abundance (rank) and occupancy (percentage and rank), 
these stats will form a table in the database (added to DB in step 2). In this code chunk 
saved to `r params$Output_dir`/Supplementary/Tables_in_SQL subfolder for reference.

**r**
```{r get individual OTU abundance and occupancy,eval=FALSE}
#lets get OTUs total abundance across all remaining samples
abundance_stats=data.frame(hit=colnames(OTU_tab_sub_occ_dec),abundance=colSums(OTU_tab_sub_occ_dec))
#get OTU abundance rank no.. where there are ties both values get the same rank i.e if two values that would be ranked 6 and 7 are the same they will both be ranked 6
#contextualise the number by referencing the total amount of OTU/ASVs
abundance_stats$abundance_rank=paste(rank(-abundance_stats$abundance,ties.method="min"),ncol(OTU_tab_sub_occ_dec),sep="/")
#get presence absence again for remaining taxa 
OTU_tab_sub_occ_dec_pa=(OTU_tab_sub_occ_dec !=0)*1
#get occupancy by summing cols using presence absense version of abundance table
abundance_stats$occupancy=colSums(OTU_tab_sub_occ_dec_pa)
#get occupancy as a percentage and add rank
abundance_stats$occupancy_proportion=paste(round(abundance_stats$occupancy/nrow(OTU_tab_sub_occ_dec_pa)*100,2),"% (Rank: ",rank(-round(abundance_stats$occupancy/nrow(OTU_tab_sub_occ_dec_pa)*100,2),ties.method="min"),"/",ncol(OTU_tab_sub_occ_dec_pa),")",sep="")
#remove unnecessary columns
abundance_stats=abundance_stats[,-c(2,4)]
write.csv(abundance_stats,paste0(Output_dir_with_occ,"/Supplementary/Tables_in_SQL/abundance_stats.csv"),row.names=FALSE)
####
#####

