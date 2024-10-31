# ---
# title: "Generate taxonomic explorer app inputs"
# author: "Briony Jones"
# date: "07/12/2022"
# output: 
#   html_document:
#     fig_caption: TRUE
#     toc: TRUE
#     toc_depth: 3
#     toc_float:
#       collapsed: FALSE
#       smooth_scroll: TRUE 
# params: 
#   OTU_tab_file: "/raid3/scratch/MEGshared/Fungi_Explorer_App/Fungi_Explorer/Inputs/Data/basic_preprocessed/OTU_tab.csv"
#   Tax_file: "/raid3/scratch/MEGshared/Fungi_Explorer_App/Fungi_Explorer/Inputs/Data/basic_preprocessed/Taxonomy.csv"
#   Fasta_file: "/raid3/scratch/MEGshared/Fungi_Explorer_App/Fungi_Explorer/Inputs/fasta/recreated_fasta/CS_Fungi.fasta"
#   Env_file: "/raid3/scratch/MEGshared/Fungi_Explorer_App/Fungi_Explorer/Inputs/Data/basic_preprocessed/Env.csv"
#   App_template_input_dir: "/raid3/scratch/MEGshared/Fungi_Explorer_App/Blank_app_framework/Generic_inputs"
#   Map_objs_input_dir : "/raid3/scratch/MEGshared/Fungi_Explorer_App/Fungi_Explorer/Inputs/UK_map_objs"
#   Output_dir: "/raid3/scratch/MEGshared/Fungi_Explorer_App/Fungi_Explorer_unpublished_app_experiments/Outputs"  
#   OTU_tab_occ_filter: 500
#   SQL_database_name: "molecular"
#   SQL_database_host: "connect-apps.ceh.ac.uk"
#   Schema_table_prefix: "fungal" 
#   Empty_database: FALSE
#   Use_occupancy_in_schema_table_names: TRUE
#   App_title: "Fungal Explorer"
#   Example_sequence: "CTACCTGATCCGAGGTCAACCTTGGTGCCGCCGGAGCGGGCTTGAGGGGGGTTTAGAGGCCGGATAGCCCGCAGGCTCCCGATGCGAGGCAGATGTTACTACGCAAAGGAAGGGCCCAACGGGTCCGCCACTGGTTTTCGGGGACTGCCTGGGCAGATCCCCAACGCCGGGCCACGGGGGCTCGAGGGTTGAAACGACGCTCGGACAGGCATGCCTCCCAGGATAC"
#   Info_text: "This is a useful description of my app"
#   External_link: "https://www.ceh.ac.uk/"
#   External_link_name: "More about UKCEH" 
  
# ---
# **r**
# ```{r setup, include=FALSE,eval=FALSE}
# knitr::opts_chunk$set(echo = TRUE)
# ```


## Summary
## Aim of script is to automate process of generating taxonomic explorer app linking marker gene sequences to environmental responses 
## (see ID-TaxER https://shiny-apps.ceh.ac.uk/ID-TaxER/ for example of similar existing app). Script directly adds to SQL database and generates blast database for app backend and modifies r shiny template file to produce front end. Preprocessed tables (otu/taxonomy etc) and map objects are also saved locally for reference. At the moment this script is suitable for UK data only as uses UK map outline objects to run, should be suitable for all CS molecular datasets (which was really my motivation for writing this). Due to complexity of script would recommend running code chunk by code chunk, rather than knitting it (atleast while it is still being tested). Apps for multiple taxonomic datasets can be generated using repeated runs of this script and added to the same database if they share the same environmental data.

### Script Input requirements:
# * OTU_tab_file should have OTU's/ASV's as column names and sample IDs as row
#   names (row names must be able to be matched to Env file row names- but do
#   not need to be in the same order)

# * Tax_file with two columns one with OTU/ASV names and the other with
#   taxonomic classification delimited by ";" (OTUs/ASVs must be able to be
#   matched to OTU_tab, but do not need to be in the same order)

#* Fasta_file with OTU representative sequences/ ASV sequences corresponding to
#  taxa in OTU_tab_file and Tax_file

#* Env_file should contain the following column names
#  "avc_code","avc","pH","eastings","northings, row names should be sample ID's
#  (these should be able to be matched to Env file but don't need to be in the
#  same order)

# * App_template_input_dir location of shiny app template 

# * OTU_tab_occ_filter: The minimum amount of samples an OTU should occur in to be included in app

#* Map_objs_input_dir: Directory with input map objects files
#  ''ukcoast_line.shp' and 'ukcoast1.shp'

# * Output_dir: Directory for all outputs

#* SQL_database_name: Name of database to input pre-processed tables and maps,
#  to be used as shiny app back end

# * SQL_database_host: Database host address for data inputs

#* Schema_table_prefix: This string will be added to database schema and table
#  names (e.g could name after taxonomic kingdom if multiple taxonomic kingdoms
#  sharing the same database)

#* Empty_database: TRUE/FALSE, if database is empty all database schemas will be
#  created; env table and uk mapping object will be added to the database
#  alongside OTU related data. If FALSE It will be assumed the database already
#  contains data for other taxonomic datasets and that the majority of schemas
#  will already exist.will have already been created and env table and mapping
#  object should already be exist in the database.

#* Use_occupancy_in_schema_table_names: TRUE/FALSE, if TRUE will add occupancy
#  into otu_attributes schema and into abund_tables table name, this is useful
#  if initially trialing different occupancy filtering in app (although will
#  likely want to go back and delete unnecessary schemas/ tables later).

#* App_title: Title to appear on the app header (not linked to future web
#  address)

#* Example_sequence: Sequence to be used as an example sequence for users when
#  exploring the app

# * Info_text: Text to feature on app describing purpose and function

# * External_link: External link you would like app to feature

# * External_link_name: Name to label external link

#*NOTE*-Additionally SQL_USER and SQL_PWD should be set in global environment to
#corresponding to the database username and password, have not included these as
#parameters to avoid hard coding DB credentials. These will be needed both to
#create the app and to subsequently run the app.


### Main Script Outputs
#All script outputs are stored either into `r params$Output_dir`/App or `r
#params$Output_dir`/Supplementary sub folders depending on if they are needed
#for the App to run or whether they are just saved for reference. Some outputs
#are also saved directly to the `r params$SQL_database_name` database.

#*Shiny front end file with ui and server code modified using rmarkdown
#parameters

#*Env file with just "avc_code","avc","pH" columns (these are the variables
#which will be used by the app itself) stored in `r params$SQL_database_name`
#database and saved to `r params$Output_dir`/Supplementary/Tables_in_SQL
#subfolder for reference.

#* Preprocessed OTU abundance table. Preprocessing includes removing samples
#  with a read number <5000 and OTUs with an occupancy of < `r
#  params$OTU_tab_occ_filter` as well as normalisation using decostands "total"
#  method (relative abundance). Table stored in `r params$SQL_database_name` and
#  saved to `r params$Output_dir`/Supplementary/Tables_in_SQL subfolder for
#  reference. 

#* Abundance_stats table with individual OTU abundance rank and occupancy
#  percentage/rank for all taxa meeting occupancy filter >= `r
#  params$OTU_tab_occ_filter`. Table stored in `r params$SQL_database_name` and
#  saved to `r params$Output_dir`/Supplementary/Tables_in_SQL subfolder for
#  reference.

#* Taxonomy table,split into
#  "kingdom","phylum","class","order","family","genus","species" fields. Table
#  filtered to only include taxa with an occupancy of >= `r
#  params$OTU_tab_occ_filter`. Table stored in `r params$SQL_database_name`
#  database and saved to `r params$Output_dir`/Supplementary/Tables_in_SQL
#  subfolder for reference.

#* Map objects per OTU/ASV (meeting occupancy filter >= `r
#  params$OTU_tab_occ_filter`) stored within `r params$SQL_database_name`
#  database and also saved to `r params$Output_dir`/Supplementary/Map_objects
#  subfolder for reference

#* Blast database for all sequences within `r params$Fasta_file` (meeting
#  occupancy filter >= `r params$OTU_tab_occ_filter`) saved to `r
#  params$Output_dir`/App/Blast_DB

## Step 0 initial setup

# Load R libraries

# **r**
# ```{R load libraries, eval=FALSE}
library(data.table)
library(vegan)
library(maptools)
library(rgdal)
library(ggplot2)
library(gstat)
library(parallel)
library(splitstackshape)
library(RPostgreSQL)
# ```


# Tweak output dir depending on occurrence cut off used
# Then pass this modified output dir to bash for command line steps

# **r**
# ```{R, eval=FALSE}
#lets add subdir to Output_dir specifying occupancy filter used
Output_dir_with_occ=paste0(params$Output_dir,"/occ_threshold_",params$OTU_tab_occ_filter)
#need to use this variable in bash chunk.. one way to do this is to export to bash like so
#https://github.com/yihui/knitr-examples/blob/master/061-bash-variable.md
Sys.setenv(Output_dir_with_occ = Output_dir_with_occ)
# ```

# Connect to specified PostgreSQL DB through R using global variables for authentication

# **r**


# ```{R, eval=FALSE}
drv=dbDriver("PostgreSQL")
conn<- dbConnect(drv, 
                 dbname = params$SQL_database_name,
                 host = params$SQL_database_host, 
                 port = '5432',
                 user=Sys.getenv('SQL_USER'),
                 password = Sys.getenv('SQL_PWD'))
# ```


## Step 1 Prepare tables for database

### Prepare  OTU tab and Env 

# Need to preprocess data, first will remove low read samples from OTU table.
# Then will remove OTUs with very low occupancy from OTU table (using `r
# params$OTU_tab_occ_filter` threshold) and normalise (decostand). Not going to
# filter samples from env based on OTU table read numbers in database -incase
# other taxonomic datasets added to database at a later date, to ensure all
# database OTU tables can be matched to env. Matching of samples between env and
# OTU table will be done within the app code itself.  
# Will filter out env samples with no AVC info for database as the purpose of
# the env info being there is to produce the habitat boxplots Will save modified
# OTU tab and env file to `r params$Output_dir`/Supplementary/Tables_in_SQL
# subfolder for reference.


## Step 1.1
#' Clean up enviromental metadata
#'
#' @description
#' Clean up enviromental metadata
#' 
#' @details
#' Takes in pre-prcoessed data and cleans it up:
#' - Removes unneeded fields
#' - rearranges into a format SQL is happier with
#' - Write new data to file
#'
#' @param input_file pre-processed intput file
#' @param output_file name to use for the output file
#' 
#' @return None
#'
#' @examples
#'  dir.create(paste0(Output_dir_with_occ,"/Supplementary/Tables_in_SQL"),
#'             showWarnings = FALSE,recursive=TRUE)
#' 
#' clean_environmental_metadata(params$Env_file # basic_preprocessed/Env.csv
#'                              paste0(Output_dir_with_occ,"/Supplementary/Tables_in_SQL/Env.csv"
#'
#' @note
#' Original code provided output variable: Env_for_SQL

clean_environmental_metadata <- function(input_file, output_file){
                                        #read in Env
    Env=data.frame(fread(input_file),row.names=1,check.names=FALSE)
                                        # dont need eastings and northings saved in database
                                        # as map objects made generated in advance
    Env_for_SQL=Env[,c("avc_code","avc","pH")]
                                        # filter out Env rows without avc code 
    Env_for_SQL=Env_for_SQL[-which(is.na(Env_for_SQL$avc_code)),]
                                        # Rearrange slightly so suitable for inserting into
                                        # SQL- e.g make sample a column(rather than rownames)
                                        # and change "pH" colname to "ph" as standard
                                        # postgres field names without double quotes have to be
                                        # lower case
                                        # https://deeplearning.lipingyang.org/2017/01/07/
                                        # postgresql-column-names-of-a-table-are-case-sensitive/
    Env_for_SQL=data.frame(sample=row.names(Env_for_SQL),Env_for_SQL[,1:2],ph=Env_for_SQL[,3])
  write.csv(Env_for_SQL,output_file)
}

## Step 1.2
#' Clean up the OTU table
#'
#' @description
#' Clean up the OTU table
#' 
#' @details
#' Filters out rows (OTU) by occupancy level
#' - Takes in Data from a raw OTU table
#' - Filters based on occupancy level
#'
#' @param input_file pre-processed OTU table
#' @param output_file name to use for the output file (filtered OTU table)
#' 
#' @return None
#'
#' @examples
#'
#' # Need to create dir first
#' dir.create(paste0(Output_dir_with_occ,"/Supplementary/Tables_in_SQL"),
#'            showWarnings = FALSE,recursive=TRUE)
#' 
#' clean_OTU_table(params$OTU_tab_file # basic_preprocessed/OTU_tab.csv,
#'                 Output_dir_with_occ,"/Supplementary/Tables_in_SQL/OTU_abund.csv")
#'
#' dir.create(paste0(Output_dir_with_occ,"/Supplementary/Tables_in_SQL"),
#'                    showWarnings = FALSE,recursive=TRUE)
#' 
#' @note
#' Original code provided output variable:  OTU_tab_sub_occ_dec

clean_OTU_table <- function(input_file, output_file, OTU_table_occupancy_filter=30){
    
                                        #read in OTU_tab
    OTU_tab=data.frame(fread(input_file),row.names=1,check.names=FALSE)
                                        #remove samples from OTU tab with reads less than 5000
    OTU_tab_sub<-OTU_tab[rowSums(OTU_tab)>5000,]
                                        # Convert OTU_tab_sub to presence and absence in order to
                                        # filter OTU_tab by OTU occupancy (i.e how many
                                        # samples an OTU is present in)
    OTU_tab_sub_pa=(OTU_tab_sub !=0)*1
                                        #remove taxa that do not meet occupancy threshold set in
                                        # OTU_table_occupancy_filter
    OTU_tab_sub_occ<-OTU_tab_sub[,which(colSums(OTU_tab_sub_pa)>=OTU_table_occupancy_filter)]
                                        #normalise OTU tab
    OTU_tab_sub_occ_dec=vegan::decostand(OTU_tab_sub_occ,method="total")
                                        #Write OTU_tab_sub to file
                                        #first make new subdir in our outdir(if doesnt already
                                        # exist) tO specify these are the tables that will be
                                        # stored in SQL 
    write.csv(OTU_tab_sub_occ_dec,output_file)
}

## Step 1.3
#' Prepare abundance_stats table
#'
#' @description
#' Prepare abundance_stats table
#' 
#' @details
#'
#' #Get individual OTU stats to summarise abundance (rank) and occupancy
#' (percentage and rank), these stats will form a table in the database (added to
#' DB in step 2). In this code chunk saved to `r
#' params$Output_dir`/Supplementary/Tables_in_SQL subfolder for reference.
#'
#' @param input_file filtered OTU table (created in step 1.2
#' @param output_file name to use for the output file, filtered abundace stats
#' 
#' @return None
#'
#' @examples
#'
#' This requires a variable "OTU_tab_sub_occ_dec" created in another functions
#' and written to a file in clean_OTU_table
#' get_abundance_stats(paste0(Output_dir_with_occ,"/Supplementary/Tables_in_SQL/OTU_abund.csv",
#'                     paste0(Output_dir_with_occ,"/Supplementary/Tables_in_SQL/abundance_stats.csv"))
#' 
#' @note
#' Original code provided output variable:  abundance_stats

get_abundance_stats <- function(input_file, output_file){

                                        # get the otu table
    OTU_tab_sub_occ_dec = data.frame(fread(input_file),row.names=1,check.names=FALSE)
    
                                        #lets get OTUs total abundance across all remaining samples
    abundance_stats=data.frame(hit=colnames(OTU_tab_sub_occ_dec),abundance=colSums(OTU_tab_sub_occ_dec))
                                        # get OTU abundance rank no.. where there are ties both values
                                        # get the same rank i.e if two values that would be ranked 6
                                        # and 7 are the same they will both be ranked 6
                                        # contextualise the number by referencing the total
                                        # amount of OTU/ASVs
    abundance_stats$abundance_rank=paste(rank(-abundance_stats$abundance,ties.method="min"),ncol(OTU_tab_sub_occ_dec),sep="/")
                                        #get presence absence again for remaining taxa 
    OTU_tab_sub_occ_dec_pa=(OTU_tab_sub_occ_dec !=0)*1
                                        #get occupancy by summing cols using presence absense version of abundance table
    abundance_stats$occupancy=colSums(OTU_tab_sub_occ_dec_pa)
                                        #get occupancy as a percentage and add rank
    abundance_stats$occupancy_proportion=paste(round(abundance_stats$occupancy/nrow(OTU_tab_sub_occ_dec_pa)*100,2),"% (Rank: ",rank(-round(abundance_stats$occupancy/nrow(OTU_tab_sub_occ_dec_pa)*100,2),ties.method="min"),"/",ncol(OTU_tab_sub_occ_dec_pa),")",sep="")
                                        #remove unnecessary columns
    abundance_stats=abundance_stats[,-c(2,4)]
    write.csv(abundance_stats,output_file,row.names=FALSE)
}

## Step 1.4
#' Prepare taxonomy table
#'
#' @description
#' Prepare abundance_stats table
#' 
#' @details
#' Read in taxonomy, split taxonomic fields by ";"and subset to only include OTUs
#' that meet occupancy threshold ,save to `r params$Output_dir`/Supplementary/Tables_in_SQL
#' subfolder for reference.
#'
#'
#' @param input_file filtered OTU table (created in step 1.2
#' @param output_file name to use for the output file, filtered abundace stats
#'
#' @param OTU_abund_filter_file filtered OTU table created in step 1.2 for input
#' @param input_file pre-processed taxonomy table to use as input
#' @param output_file name for output file, which is a filtered taxonomy table
#' 
#' @return None
#'
#' @examples
#'
#'
#' prepair_taxonomy_table(
#'     OTU_abund_filter_file = paste0(Output_dir_with_occ,"/Supplementary/Tables_in_SQL/OTU_abund.csv",
#'     input_file = params$Tax_file # basic_preprocessed/Taxonomy.csv,
#'     output_file = paste0(Output_dir_with_occ,"/Supplementary/Tables_in_SQL/Taxonomy.csv")
#' )
#'
#' @note
#' Original code provided output variable:
#' Taxonomy_filt (Taxonomy_Sort is now written as extra steps were needed.)
#' 
prepair_taxonomy_table <- function(input_file, output_file, OTU_abund_filter_file){
    
    Taxonomy <- read.csv(input_file)
                                        # data is two columns (OTU and taxa).
                                        # Split into multiple taxonomy columns (taxonomy_1, etc)
    Taxonomy<-as.data.frame(cSplit(indt=Taxonomy,splitCols=2,sep=";"))
                                        # rename taxonomy columns to correct ¿levels?
    colnames(Taxonomy)=c("hit","Kingdom","Phylum","Class","Order","Family","Genus","Species")
    row.names(Taxonomy)=Taxonomy$hit
    
                                        # rel Current code does not generate output as expected.
                                        # rel added this to ensure match
    Taxonomy[,'hit'] <- NULL  
                                        #filter to match OTU table
    OTU_abund_filter <- colnames(data.table::fread(OTU_abund_filter_file))
    Taxonomy_filt <- Taxonomy[OTU_abund_filter,]
    Taxonomy_Sort <- Taxonomy_filt[ order(row.names(Taxonomy_filt)),]
    
                                        # rel Current code does not generate output as expected.
                                        # rel added this to ensure match by cleaning up entries
    Taxonomy_Sort$Kingdom <- sub(".*__", "", Taxonomy_Sort$Kingdom)
    Taxonomy_Sort$Phylum <- sub(".*__", "", Taxonomy_Sort$Phylum)
    Taxonomy_Sort$Class <- sub(".*__", "", Taxonomy_Sort$Class)
    Taxonomy_Sort$Order <- sub(".*__", "", Taxonomy_Sort$Order)
    Taxonomy_Sort$Family <- sub(".*__", "", Taxonomy_Sort$Family)
    Taxonomy_Sort$Genus <- sub(".*__", "", Taxonomy_Sort$Genus)
    Taxonomy_Sort$Species <- sub(".*__", "", Taxonomy_Sort$Species)
    Taxonomy_Sort[is.na(Taxonomy_Sort)] <- ""
    write.csv(Taxonomy_Sort,output_file)
}


#####################################################################
#####################################################################

### ALL DB STUFF AND THEN APP NOW

#####################################################################
#####################################################################

## Steps 2.1, 2.2, and 2.3 have all been replaced with a single "Step 2"
## Step 2.1 used to setup the SQL db
## Step 2.2 modified the table name based on the abundance cut
## Step 2.3 populated the tables
## The single step 2 now creates the SQLite files and populates them in
## one function.


## Step 2
#' Put the filtered OTU, taxonomy, and abundance files into RSQLite files
#'
#' @description
#' Put the filtered OTU, taxonomy, and abundance files into RSQLite files
#' 
#' @details
#' Reads in the filtered OTU, taxonomy, and abundance files, and creates corresponding
#' tables in an RSQLite file
#'
#' @param input_file filtered OTU table
#' @param input_file filtered abundance table
#' @param input_file filtered taxonomy table
#' 
#' @return None
#'
#' @examples
#'
#'
#' prepair_taxonomy_table(
#'     OTU_abund_filter_file = paste0(Output_dir_with_occ,"/Supplementary/Tables_in_SQL/OTU_abund.csv",
#'     input_file = params$Tax_file # basic_preprocessed/Taxonomy.csv,
#'     output_file = paste0(Output_dir_with_occ,"/Supplementary/Tables_in_SQL/Taxonomy.csv")
#' )
#'
#' @note
#'
#' format_otu_for_Rsqlite(filtered_abundance_csv, filtered_taxonomy_csv, filtered_otu_csv)
#' 
#' 
format_otu_for_Rsqlite <- function(abundance_csv, taxonomy_csv, otu_csv){
    ## Lets try not splitting to begin with
    ## Current SQL command is:
    ## CREATE TABLE IF NOT EXISTS abund_tables.<prefix>_abund_2(hit character varying(30), "<column names each followed by '" numeric, "'>" numeric, CONSTRAINT '<prefix>_abund_2_pkey PRIMARY KEY("hit"))'
    ## Lets get rid of the prefix and double tables to:
    ## CREATE TABLE IF NOT EXISTS abund_tables(hit character varying(30), "<column names each followed by '" numeric, "'>" numeric, CONSTRAINT '<prefix>_abund_2_pkey PRIMARY KEY("hit"))'

    ## Lets use sprintf to get better formatting of the string.
    ## Fails with "CONSTRAINT", not sure on need, so removing.
    
    sql_command <- sprintf("create table if not exists abund_table (hit character varying(30), %s numeric, primary key ('hit'))",
                                     paste(colnames(abundance_csv), collapse=' numeric, ')   #colnames
                           )

    ## Create table, abundance
    abundance_db <- dbConnect(RSQLite::SQLite(), "abundance_db.sqlite")
    DBI::dbExecute(conn = abundance_db, statement = sql_command)
    ## Fill table
    DBI::dbWriteTable(abundance_db, "abund_table", abundance_csv, append = TRUE, row.names = FALSE)
    ## Disconnect, best practise?
    DBI::dbDisconnect(abundance_db)


    ## Create table, taxonomy

    ## "order" is an SQL command so need to escape it in quotes, easier to escape all
    ## fields in quotes
    
    sql_command_taxonomy <- sprintf("create table taxonomy_table (hit character varying (30), %s character varying (250), primary key ('hit'))",
                                    paste('"',colnames(abundance_csv),'"',
                                          collapse=' charater varying(250),',
                                          sep='')
                                    )

    ## Create table, taxonomy
    taxonomy_db <- dbConnect(RSQLite::SQLite(), "taxonomy_db.sqlite")
    DBI::dbExecute(conn = taxonomy_db, statement = sql_command)
    ## Fill table
    DBI::dbWriteTable(taxonomy_db, "taxonomy_table", taxonomy_csv, append = TRUE, row.names = FALSE)
    ## Disconnect, best practise?
    DBI::dbDisconnect(taxonomy_db)


    ## Create OTU table

    sql_command_otu <- sprintf("create table otu_table (hit character varying (30), %s character varying (30), primary key ('hit'))",
                                    paste('"',colnames(otu_csv),'"',
                                          collapse=' charater varying(250),',
                                          sep='')
                                    )

    ## Create table, otu
    otu_db <- dbConnect(RSQLite::SQLite(), "otu_db.sqlite")
    DBI::dbExecute(conn = otu_db, statement = sql_command)
    ## Fill table
    DBI::dbWriteTable(otu_db, "otu_table", otu_csv, append = TRUE, row.names = FALSE)
    ## Disconnect, best practise?
    DBI::dbDisconnect(otu_db)
    

    ## Create maps table


    sql_command_maps <- sprintf("create table otu_attributes_table (hit character varying (30), map_object bytea, primary key ('hit'))")
    

    ## Create table, maps
    maps_db <- dbConnect(RSQLite::SQLite(), "maps_db.sqlite")
    DBI::dbExecute(conn = maps_db, statement = sql_command)

    #Fill maps later

    ## Fill table
##    DBI::dbWriteTable(maps_db, "maps_table", maps_csv, append = TRUE, row.names = FALSE)
    ## Disconnect, best practise?
    DBI::dbDisconnect(maps_db)

                                        #Create empty maps table within otu_attributes_schema
    dbExecute(conn=conn,paste0('CREATE TABLE IF NOT EXISTS ', # SQL command
                               Schema_table_prefix_modified,  #table name start
                               '_otu_attributes.',
                               Schema_table_prefix_modified, # table name ends as '_maps'
                               '_maps(hit character varying(30),map_object bytea, CONSTRAINT '
                              ,Schema_table_prefix_modified, #primary key_name
                               '_maps_pkey PRIMARY KEY ("hit"));'
                               )
              )
                                        #will fill maps table in step 3 using save_otu_map function


    
}



## Step 3 Make map objects for DB

### 3.1 Map preparation

# **r**
# ```{R map prep, eval=FALSE}


#######
#######
#######
## rel notes ##
##  looks like this sets up a series of functions and variables to be used in a latter step to make a db.
## This last part maybe much harder to split up.   Need to focus on first parts and leave this until later.
## Perhaps it would be better to move to SQLite????  Perhaps some of the functions to setup, like grid bit
## below could be separated out??
############

## Inputs?
## var= OTU_tab_sub_occ_de
## file= /ukcoast1.shp
## file= /ukcoast_line.shp

map_prep <- function(){
                                        # first lets get env in the same order as otu table
    Env_sub = Env[row.names(OTU_tab_sub_occ_dec),]
                                        # Isnt this just a "test" so does nothing?
    identical(row.names(Env_sub),row.names(OTU_tab_sub_occ_dec))
                                        # read in some shapefiles (georeferenced maps of uk)
                                        # uk polygon
    uk.poly <- rgdal::readOGR(paste0(params$Map_objs_input_dir,"/ukcoast1.shp"))
                                        # outline of uk
    uk.line <- rgdal::readOGR(paste0(params$Map_objs_input_dir,"/ukcoast_line.shp"))

                                        # British National Grid System is square for UK.
                                        # Need to project maps long-lat onto this as
                                        # country side survey data is in easting/northings
                                        # grid system.
                                        # grid references of maps are in latitude longitude,
                                        # but CS data coords are in eastings/northings
                                        # some conversion trickery is required
                                        # before trickery we can define some conversion parameters
    ukgrid = "+init=epsg:27700"
    latlong = "+init=epsg:4326"
                                        # tell r what the current projections are of the maps
                                        # set the CRS (coordinate reference system) field
    uk.poly@proj4string = CRS(latlong)
                                        # set the CRS field
    uk.line@proj4string = CRS(latlong)
                                        # transform
                                        # transform coordinates to eastings northings
    uk.poly = spTransform(uk.poly, CRS(ukgrid))
                                        # transform coordinates to eastings northings
    uk.line = spTransform(uk.line, CRS(ukgrid))
    
                                        # check
    plot(uk.poly)
    plot(uk.line,add=T,col="red")
    
                                        # make a grid spanning uk, with which to interpolate
                                        # not sure what this does really
    g <- fortify(uk.poly)
    coordinates(g) = ~long+lat
    x.range <- as.integer(range(g@coords[,1]))
    y.range <- as.integer(range(g@coords[,2]))
                                        # specify spatial data as being gridded
    grd <- expand.grid(x=seq(from=x.range[1], to=x.range[2], by=5000),
                       y=seq(from=y.range[1], to=y.range[2], by=5000)
                       )

    coordinates(grd) <- ~ x+y

    gridded(grd) <- TRUE
                                        # check
    plot(grd)
    plot(uk.poly,add=T,col="red")
}


### 3.2 Generate maps per OTU


                                        # **save_otu_map** function  creates map object representing an individual OTUs/ ASVs geographical distribution this is saved to both the specified outputdir and database. Pngs of plotted map objects can also be saved to outdir for reference.

                                        # Arguments

                                        # * OTU_name- OTU/ASV name

                                        # * OTU_table- colnames should be OTUs, rows should be samples, row/sample order should be consistent with Env_table

                                        # * Env_table - table with environmental metadata including 'eastings' and 'northings',cols should be environmental variables, rows should be samples, row/sample order should be consistent with Env table

                                        # * Grid- Spatial pixels object for interpolation

                                        # * UK_poly- Spatial polygons dataframe
                                        # * UK_line- Spatial lines dataframe
                                        # * Eastings_col- Eastings column name within Env_table
                                        # * Northings_col- Northings column name within Env_table
                                        # * Conn- connection to postgres database
                                        # *Schema_table_prefix - prefix of otu_attributes schema and maps table
                                        # * Output_dir- Output directory for map objects
                                        # * Make_png- TRUE/FALSE if true png will be generated of map as well as map object

                                        # Example arguments for function testing
                                        # OTU_name="ASV_1"
                                        # OTU_table=OTU_tab_sub_occ_dec
                                        # Env_table=Env_sub
                                        # Grid=grd
                                        # UK_poly=uk.poly
                                        # UK_line=uk.line
                                        # Conn=conn
                                        # Schema_table_prefix=Schema_table_prefix_modified
                                        # Output_dir=paste0(Output_dir_with_occ,"/Supplementary/Map_objects")
                                        # Make_png=TRUE

                                        # **r**
                                        # ```{R maps function, eval=FALSE}

map_function <- function(){}
save_otu_map<-function(OTU_name,OTU_table,Env_table,Grid,UK_poly,UK_line,Conn,Schema_table_prefix,Output_dir,Make_png){
    otu_abund<-OTU_table[,OTU_name,drop=FALSE]
                                        #make dataframe with eastings and northings 
    dat<-cbind(Env_table[,c('eastings','northings')],otu_abund)
    colnames(dat)[3]="OTU"
    dat<-dat[complete.cases(dat), ]
    attach(dat)
                                        #specify the coordinates for the file with otu
    sp::coordinates(dat) = ~eastings+northings
                                        #interpolate
    spc.idw<-gstat::krige(OTU~1,dat,Grid)
    ukgrid = "+init=epsg:27700"
    spc.idw@proj4string<-sp::CRS(ukgrid)
                                        #now need to only show interpolation within UK coastline
    overlay.idw<-sp::over(spc.idw,UK_poly)
    newmap=spc.idw[complete.cases(overlay.idw),]
                                        #set heat map scale, want this to be dependent on OTU abundance i.e we want to show abundance contrast relative to that OTU
    minv=min(otu_abund)
    maxv=max(otu_abund)
    at=c(minv,signif(maxv/256,1),signif(maxv/128,1),signif(maxv/64,1),signif(maxv/32,1),signif(maxv/16,1),signif(maxv/8,),signif(maxv/4,2),signif(maxv/2,2),maxv)
                                        #save map object to outdir
    mapandinfo<-c(newmap,at)
    save(mapandinfo,file=paste(Output_dir,"/",OTU_name,".RData",sep=""))
                                        #Add object to database  
                                        #need to first get object into a form suitable for SQL
                                        #convert to stream of bytes
    ser_mapandinfo=serialize(mapandinfo,connection=NULL,ascii=TRUE)
                                        #convert to a form postgres will accept
    bytea_ser_mapandinfo<-RPostgreSQL::postgresqlEscapeBytea(ser_mapandinfo,con=Conn)
    dbSendQuery(Conn,paste0("INSERT INTO ", Schema_table_prefix,"_otu_attributes.",Schema_table_prefix,"_maps VALUES ('",OTU_name,"','",bytea_ser_mapandinfo,"')"))
                                        #save png vs if make_png=TRUE
    if (Make_png==TRUE){ 
        map_plot=sp::spplot(newmap["var1.pred"], at = at,sp.layout=list("sp.lines",UK_line,col="black",lwd=2))
        dir.create(paste0(Output_dir,"/png"), showWarnings = FALSE,recursive=TRUE)
        png(file=paste0(Output_dir,"/png/",OTU_name,"_plot.png"),width=200,height=400)
        print(map_plot)
        dev.off()
    }
}





# step 3.2b
                                        # ```
# Run save_otu_map function on all OTUs.

# Using rparallel to parallelise -need to reauthenticate database on all cpus working on.

# **r**
# ```{R maps parallelise, eval=FALSE}
maps_parallelise <- function(){
                                        #create outdir for map objects
    dir.create(paste0(Output_dir_with_occ,"/Supplementary/Map_objects"), showWarnings = FALSE,recursive=TRUE)
                                        #parallelise on 40 CPU
    cl<-makeCluster(40)
                                        #make r objs acceptable to all CPU
    clusterExport(cl,c("save_otu_map","OTU_tab_sub_occ_dec","Env_sub","grd","params","Output_dir_with_occ","Schema_table_prefix_modified","uk.poly","uk.line"))
    
    
                                        #need to connect to the database on all 40 cpu
    clusterEvalQ(cl, {
        library(RPostgreSQL)  
        drv=DBI::dbDriver("PostgreSQL")
        conn<- DBI::dbConnect(drv, 
                              dbname = params$SQL_database_name,
                              host = params$SQL_database_host, 
                              port = '5432',
                              user=Sys.getenv('SQL_USER'),
                              password = Sys.getenv('SQL_PWD'))
        
    })
    
    
    
    parSapply(cl,colnames(OTU_tab_sub_occ_dec), function(x) save_otu_map(OTU_name=x,OTU_table=OTU_tab_sub_occ_dec,Env_table=Env_sub,Grid=grd,UK_poly=uk.poly,UK_line=uk.line,Conn=conn,Schema_table_prefix=Schema_table_prefix_modified,Output_dir=paste0(Output_dir_with_occ,"/Supplementary/Map_objects"),Make_png=FALSE))
    stopCluster(cl)
    
                                        # ```
}

## Step 3 Make Blast database

###  3.3 Filter fasta file to contain OTU sequences that meet occupancy filter

#Need to filter sequences prior to making blast DB otherwise we will have hits
#in BLast DB with no supplementary info (i.e we dont want to keep sequences
#corresponding to OTU's/ ASV's that did not meet occupancy filter)
# Doing this using biopython.


make_blast_py <- function(){
# **python**
# ```{python filter fasta file, eval=FALSE}
from Bio import SeqIO
import os
#get OTUs we want to keep from OTU table
#python can access r variables within markdown and converts our r dataframe into a dictionary- the keys are the OTU names
OTU_tab_dict=r.OTU_tab_sub_occ_dec
#make output dir for filtered fasta
if not os.path.exists(r.Output_dir_with_occ+"/Supplementary/Filtered_sequences"):
  os.makedirs(r.Output_dir_with_occ+"/Supplementary/Filtered_sequences")
#only keep records in filtered otu tab (dictionary keys) 
#see http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec372 for reference
records=(record for record in SeqIO.parse(r.params["Fasta_file"],"fasta") if record.id in OTU_tab_dict.keys())
SeqIO.write(records, r.Output_dir_with_occ+"/Supplementary/Filtered_sequences/filtered_sequences.fasta", "fasta")
# ```
}

###  3.4 Make blast database

# Take filtered fasta and make blast database for back end of shiny app

make_blast_bash <- function(){
# **Bash:**
# ```{bash BlastDB, eval=FALSE}
#cmake makes dir and any parent dirs necessary
cmake -E make_directory $Output_dir_with_occ/App/Blast_DB
makeblastdb -in $Output_dir_with_occ/Supplementary/Filtered_sequences/filtered_sequences.fasta -dbtype nucl -out $Output_dir_with_occ/App/Blast_DB/Blast_DB
# ```
}




create_app_python <-function(){
#### 4 Create app front-end from template

# Using python to edit r code - find python easiest option for file handling

# **Python:**
# ```{python create app, eval=FALSE}
import re
#get blank template path using R App_template parameter  
blank_template=r.params["App_template_input_dir"]+"/Blank_taxonomic_explorer_app_with_functional_place_holders.R"
#open blank_template
with open(blank_template,'r') as blank_template_file:
#open output file  
  with open(r.Output_dir_with_occ+"/App/App.R",'w') as customised_file:
#loop over lines 
    for line in blank_template_file:
#empty array for any placeholder strings that need to be replaced within that line 
      occurence_strings=[]
#get all indices where placeholder  starts("<--")
      occurence_starts= [i.start() for i in re.finditer("<--",line)]
#if there are any beginnings of placeholders in the line  continue   
      if (len(occurence_starts)!=0):
#get all indices  of the end of place holder/s "-->"  
        occurence_ends= [i.start() for i in re.finditer("-->",line)]
#loop over number of start indices per line (as may be multiple placeholders per line)     
        for i in range(len(occurence_starts)):
#get whole place holder string  ,occurence_ends are the start indices of "-->" string so add 3
          occurence_string=line[occurence_starts[i]:occurence_ends[i]+3]
#if not in occurrence_strings array already append
          if (occurence_string not in occurence_strings):
            occurence_strings.append(occurence_string)
#define new line        
        new_line=line
#for every unique place holder string replace placeholder with relevant r parameter(exception being Schema_table_prefix_modified as it is not a parameter)   
        for occurence_string in occurence_strings:
          parameter=occurence_string.split("<--specified_")[1].split("-->")[0]
          if(parameter=="Schema_table_prefix"):
            new_line=new_line.replace(occurence_string,r.Schema_table_prefix_modified)
          else:
            new_line=new_line.replace(occurence_string,r.params[parameter])
#ok lets write the new modified line to our output file       
        customised_file.write(new_line)
#if no placeholder strings write original line to output file        
      else:  customised_file.write(line)      

# ```
# App should now be an executable 
}







































########

## Old functions for SQL, wnat to keep to hand for a bit for easier comparison.

########

### Step 2 Create Database tables, and populate with tables produced in Step 1



# Step 2.1

#### If this is a blank DB (e.g no other taxonomic datasets have previously been added) more setup will be needed
# Create all schemas and add env & plotting tools data

# ```{r setup database if no previous data added,eval=FALSE}
#Bit belt and braces with r and sql mechanisms to stop any accidental env/ plotting tools data duplication
#if Empty_database parameter set to TRUE

setup_R_db <- function(){
if(params$Empty_database==TRUE){
  #Create empty plotting_tools schema if doesnt already exist (this is just going to be used to store generic   plotting objects- currently just uk map outline)
  dbExecute(conn=conn,"CREATE SCHEMA IF NOT EXISTS plotting_tools;")
  #create map_tools table within plotting_tools schema if doesnt already exist 
  #make description a primary key -primary keys ensure that col contains unique values
  dbExecute(conn=conn,'CREATE TABLE IF NOT EXISTS plotting_tools.map_tools(description character varying(250),plot_object bytea,CONSTRAINT map_tools_pkey PRIMARY KEY (description));')
  #add map outline to table
  #first read in 
  uk.line<-rgdal::readOGR(paste0(params$Map_objs_input_dir,"/ukcoast_line.shp"))
  #convert to stream of bytes
  ser_uk.line=serialize(uk.line,connection=NULL,ascii=TRUE)
  #convert to a form postgres will accept
  bytea_ser_uk.line=RPostgreSQL::postgresqlEscapeBytea(ser_uk.line,con=conn)
  #populate plotting_tools.map_tools table
  #if 'map_outline' primary key already exists error will be produced 
  dbSendQuery(conn, paste0("INSERT Into plotting_tools.map_tools VALUES('map_outline','",bytea_ser_uk.line,"')"))
  #create empty env_attributes schema if doesnt already exist
  dbExecute(conn=conn,"CREATE SCHEMA IF NOT EXISTS env_attributes;")
  #create env table  within env_attributes schema
  #make sample primary key
  dbExecute(conn=conn,'CREATE TABLE IF NOT EXISTS env_attributes.env_attributes_all(sample character varying(250),avc_code numeric ,avc character varying(250),ph numeric,CONSTRAINT env_attributes_all_pkey PRIMARY KEY (sample));')
  #populate env table if any primary keys already exist errors will be produced
  append_cmd=sqlAppendTable(con=conn,table=Id(schema="env_attributes",table="env_attributes_all"), values =Env_for_SQL, row.names = FALSE )
  dbExecute(conn=conn,statement=append_cmd)
  #create abund_tables schema if doesnt already exist
  dbExecute(conn=conn,"CREATE SCHEMA IF NOT EXISTS abund_tables;")
  #create abund_table_descriptions table for recording abundance table preprocessing used
  #make abundance_table_name primary key
  dbExecute(conn=conn,"CREATE TABLE IF NOT EXISTS abund_tables.abund_table_descriptions(abund_table_name character varying(40),sample_read_filt character varying(40),otu_occupancy_filt character varying(40),normalisation character varying(40),further_notes character varying(80), CONSTRAINT abund_table_descriptions_pkey PRIMARY KEY (abund_table_name));")
  }
}

# ```


# Step 2.2

#### Create and populate tables specific to this taxonomic dataset 

# Create and populate otu_abund, taxonomy and abundance_stats tables.
# Create maps table, but populate in step 3.

# ```{r Add to database,eval=FALSE}

#modify schema_table_prefix if params$Use_occupancy_in_schema_table_names is TRUE
modify_schema_table_prefix <- function(){
if (params$Use_occupancy_in_schema_table_names==TRUE){
  Schema_table_prefix_modified=paste0(params$Schema_table_prefix,"_occ_",params$OTU_tab_occ_filter)
}else{Schema_table_prefix_modified=params$Schema_table_prefix}
}

#step 2.3

#make otu table structure going to transpose our dataframe to make samples
#columns and otus rows as sql tables can only have 1600 columns-apparently also
#want to include otu names (currently colnames/ rownames when transposed) as
#column called "hit" even with colnames being samples rows, are still too big so
#separating tables into two(column wise), can easily use a join to get full
#abundance table per ASV when needed in app


##Replaced with format_otu_for_Rsqlite
format_otu_for_sql <- function(){
    otu_tab_for_SQL_precursor=t(OTU_tab_sub_occ_dec)
                                        #split into two, also add row names to both (hit col)
                                        # i.e asv names so that the tables can be joined 
    half_cols = round(ncol(otu_tab_for_SQL_precursor)/2)
    otu_tab_for_SQL_1 = data.frame(hit=row.names(otu_tab_for_SQL_precursor),
                                   otu_tab_for_SQL_precursor[,1:half_cols],
                                   check.names=FALSE)
    otu_tab_for_SQL_2 = data.frame(hit=row.names(otu_tab_for_SQL_precursor),
                                   otu_tab_for_SQL_precursor[,(half_cols+1):ncol(otu_tab_for_SQL_precursor)],
                                   check.names=FALSE)
                                        #could have used sqlCreateTable command to create table
                                        # but constructed own command using colnames of
                                        # Otu_tab_for_SQL so I could specify data types 
    table_create_cmd = paste0('CREATE TABLE IF NOT EXISTS abund_tables.',
                              Schema_table_prefix_modified,
                              '_abund_1',
                              '(hit character varying(30), "',
                              paste(colnames(otu_tab_for_SQL_1)[2:ncol(otu_tab_for_SQL_1)],
                                    collapse='" numeric, "'
                                    ),
                              '" numeric, CONSTRAINT ',
                              Schema_table_prefix_modified,
                              '_abund_1_pkey PRIMARY KEY("hit"))'
                              )
    ## statement = <valid SQL>
    ## CREATE TABLE IF NOT EXISTS abund_tables.<prefix>_abund_2(hit character varying(30), "<column names each followed by '" numeric, "'>" numeric, CONSTRAINT '<prefix>_abund_2_pkey PRIMARY KEY("hit"))'
        
    dbExecute(conn=conn,statement=table_create_cmd)
    
    table_create_cmd = paste0('CREATE TABLE IF NOT EXISTS abund_tables.',
                              Schema_table_prefix_modified,
                              '_abund_2',
                              '(hit character varying(30), "',
                              paste(colnames(otu_tab_for_SQL_2)[2:ncol(otu_tab_for_SQL_2)],
                                    collapse='" numeric, "'
                                    ),
                              '" numeric, CONSTRAINT ',
                              Schema_table_prefix_modified,
                              '_abund_2_pkey PRIMARY KEY("hit"))'
                              )
    
    dbExecute(conn=conn,statement=table_create_cmd)
                                        #add data to abundance tables

    append_cmd = sqlAppendTable(con=conn,
                                table=Id(schema="abund_tables",
                                         table=paste0(Schema_table_prefix_modified,"_abund_1")),
                                values=otu_tab_for_SQL_1,
                                row.names=FALSE )
    
    dbExecute(conn=conn,statement=append_cmd)

    append_cmd = sqlAppendTable(con=conn,
                                table=Id(schema="abund_tables",
                                         table=paste0(Schema_table_prefix_modified,"_abund_2")),
                                values=otu_tab_for_SQL_2,
                                row.names=FALSE)
    
    dbExecute(conn=conn,statement=append_cmd)

##########
## Data descriptions table, ignore adding for now, lets get it wokring. ##
#########    
                                        # want to record how table was normalised/processed in db
                                        # will be useful if multiple otu abundance tables in the
                                        # same db- add to abund_table_descriptions table
    data_descriptions = data.frame(rbind(c(paste0(Schema_table_prefix_modified,"_abund_1"),
                                           ">5000 retained",
                                           paste0(">=",params$OTU_tab_occ_filter," retained"),
                                           "decostand total",
                                           paste0("first half of ",
                                                  Schema_table_prefix_modified,
                                                  "_abund table samples")
                                           ),
                                         c(paste0(Schema_table_prefix_modified,"_abund_2"),
                                           ">5000 retained",
                                           paste0(">=",params$OTU_tab_occ_filter," retained"),
                                           "decostand total",
                                           paste0("second half of ",
                                                  Schema_table_prefix_modified,
                                                  "_abund table samples")
                                           )
                                         )
                                   )
    
    colnames(data_descriptions) = c("abund_table_name",
                                    "sample_read_filt",
                                    "otu_occupancy_filt",
                                    "normalisation",
                                    "further_notes")
    
    append_cmd = sqlAppendTable(con=conn,
                                table=Id(schema="abund_tables",
                                         table="abund_table_descriptions"),
                                values=data_descriptions,
                                row.names=FALSE
                                )
    
    dbExecute(conn=conn,statement=append_cmd)
    
    ## ######
    ## end table descriptors, start taxonoomy table
    ## ######    
    
                                        # Create otu_attributes schema(specific to this
                                        # taxonomic dataset) if doesnt already exist 
    dbExecute(conn=conn,
              paste0("CREATE SCHEMA IF NOT EXISTS ",
                     Schema_table_prefix_modified,
                     "_otu_attributes;")
              )
                                        #create empty taxonomy table within otu_attributes
                                        # schema ,field name 'order' is in double quotes
                                        # as it is also a SQL command
    
    ## CREATE TABLE IF NOT EXISTS abund_tables(hit character varying(30), "<column names each followed by '" numeric, "'>" numeric, CONSTRAINT '<prefix>_abund_2_pkey PRIMARY KEY("hit"))'

    ##     paste0('CREATE TABLE IF NOT EXISTS otu_attributes_taxonomy(hit character varying(30),kingdom character varying(250),phylum character varying(250),class character varying(250),"order" character varying(250),family character varying(250),genus character varying(250),species character varying(250), CONSTRAINT taxonomy_pkey PRIMARY KEY ("hit"))'

    ## CREATE TABLE IF NOT EXISTS otu_attributes_taxonomy(hit character varying(30),kingdom character varying(250),phylum character varying(250),class character varying(250),"order" character varying(250),family character varying(250),genus character varying(250),species character varying(250), CONSTRAINT taxonomy_pkey PRIMARY KEY ("hit"))'


    dbExecute(conn=conn,
              paste0('CREATE TABLE IF NOT EXISTS ',
                     Schema_table_prefix_modified,
                     '_otu_attributes.',
                     Schema_table_prefix_modified,
                     '_taxonomy(hit character varying(30),kingdom character varying(250),phylum character varying(250),class character varying(250),"order" character varying(250),family character varying(250),genus character varying(250),species character varying(250), CONSTRAINT ',
                     Schema_table_prefix_modified,
                     '_taxonomy_pkey PRIMARY KEY ("hit"))')
              )

                                        #add data from data frame to empty taxonomy table
    append_cmd = sqlAppendTable(con=conn,
                                table=Id(schema=paste0(Schema_table_prefix_modified,
                                                                "_otu_attributes"
                                                       ),
                                         table=paste0(Schema_table_prefix_modified,
                                                      "_taxonomy"
                                                      )
                                                  ),
                                values=Taxonomy_filt,
                                row.names=FALSE
                                )

    dbExecute(conn=conn,statement=append_cmd)
                                        #create empty abundance stats table within
                                        # otu_attributes schema

##### Create new table, OTU

    ## CREATE TABLE IF NOT EXISTS abund_tables(hit character varying(30), "<column names each followed by '" numeric, "'>" numeric, CONSTRAINT '<prefix>_abund_2_pkey PRIMARY KEY("hit"))'

    
    dbExecute(conn=conn,paste0('CREATE TABLE IF NOT EXISTS ', #sql command
                               Schema_table_prefix_modified, ## table name start
                               '_otu_attributes.',
                               Schema_table_prefix_modified, #table name ends as _abundace_stats
                               '_abundance_stats(hit character varying(30),abundance_rank character varying(30),occupancy_proportion character varying(30), CONSTRAINT ',
                               Schema_table_prefix_modified,  # priary key name start
                               '_abundance_stats_pkey PRIMARY KEY ("hit"));'
                               )
              )


        sql_command_taxonomy <- sprintf("create table otu_table (hit character varying (30), %s character varying (30), primary key ('hit'))",
                                    paste('"',colnames(otu_csv),'"',
                                          collapse=' charater varying(250),',
                                          sep='')
                                    )


    
                                        #add data from data frame to empty abundance stats table  
    append_cmd = sqlAppendTable(con=conn,table=Id(schema=paste0(Schema_table_prefix_modified,
                                                                "_otu_attributes"),
                                                  table=paste0(Schema_table_prefix_modified,
                                                               "_abundance_stats")),
                                values=abundance_stats,
                                row.names=FALSE
                                )
    
    dbExecute(conn=conn,statement=append_cmd)

    ## CREATE TABLE IF NOT EXISTS abund_tables(hit character varying(30), "<column names each followed by '" numeric, "'>" numeric, CONSTRAINT '<prefix>_abund_2_pkey PRIMARY KEY("hit"))'


    
                                        #Create empty maps table within otu_attributes_schema
    dbExecute(conn=conn,paste0('CREATE TABLE IF NOT EXISTS ', # SQL command
                               Schema_table_prefix_modified,  #table name start
                               '_otu_attributes.',
                               Schema_table_prefix_modified, # table name ends as '_maps'
                               '_maps(hit character varying(30),map_object bytea, CONSTRAINT '
                              ,Schema_table_prefix_modified, #primary key_name
                               '_maps_pkey PRIMARY KEY ("hit"));'
                               )
              )
                                        #will fill maps table in step 3 using save_otu_map function
                                        # ```
}

