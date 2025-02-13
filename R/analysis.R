#' Microscope
#'
#' Generate taxonomic explorer app inputs
#'
#' 
#' @author "Briony Jones"
#' 
#' @import data.table (>= 1.15.4)
#' @import vegan (>= 2.6_8)
#' @import terra (>= 1.7_78) 
#' @import ggplot2 (>= 3.5.1)
#' @import gstat (>= 2.1_2)
#' @import parallel 
#' @import splitstackshape (>= 1.4.8)
#' @import RSQLite (>= 2.3)
#' @import DBI (>= 1.2.3
#' @import sp (>= 2.1)
#'
#' @description Aim of script is to automate process of generating taxonomic explorer app linking marker gene sequences to environmental responses (see ID-TaxER https://shiny-apps.ceh.ac.uk/ID-TaxER/ for example of similar existing app). Script directly adds to SQL database and generates blast database for app backend and modifies r shiny template file to produce front end. Preprocessed tables (otu/taxonomy etc) and map objects are also saved locally for reference. At the moment this script is suitable for UK data only as uses UK map outline objects to run, should be suitable for all CS molecular datasets (which was really my motivation for writing this). Due to complexity of script would recommend running code chunk by code chunk, rather than knitting it (atleast while it is still being tested). Apps for multiple taxonomic datasets can be generated using repeated runs of this script and added to the same database if they share the same environmental data.
#'



#' #@param OTU_tab_file should have OTU's/ASV's as column names and sample IDs as row
#   names (row names must be able to be matched to Env file row names- but do
#   not need to be in the same order)
#' #@param Tax_file with two columns one with OTU/ASV names and the other with
#   taxonomic classification delimited by ";" (OTUs/ASVs must be able to be
#   matched to OTU_tab, but do not need to be in the same order)
#' #@param Fasta_file with OTU representative sequences/ ASV sequences corresponding to
#  taxa in OTU_tab_file and Tax_file
#' #@param Env_file should contain the following column names
#  "avc_code","avc","pH","eastings","northings, row names should be sample ID's
#  (these should be able to be matched to Env file but don't need to be in the
#  same order)
#' #@param App_template_input_dir location of shiny app template 
#'
#' #@param OTU_tab_occ_filter: The minimum amount of samples an OTU should occur in to be included in app
#' #@param Map_objs_input_dir: Directory with input map objects files
#  ''ukcoast_line.shp' and 'ukcoast1.shp'
#' #@param Output_dir: Directory for all outputs
#' #@param SQL_database_name: Name of database to input pre-processed tables and maps,
#  to be used as shiny app back end
#' #@param SQL_database_host: Database host address for data inputs
#' #@param Schema_table_prefix: This string will be added to database schema and table
#  names (e.g could name after taxonomic kingdom if multiple taxonomic kingdoms
#  sharing the same database)
#' #@param Empty_database: TRUE/FALSE, if database is empty all database schemas will be
#  created; env table and uk mapping object will be added to the database
#  alongside OTU related data. If FALSE It will be assumed the database already
#  contains data for other taxonomic datasets and that the majority of schemas
#  will already exist.will have already been created and env table and mapping
#  object should already be exist in the database.
#' #@param Use_occupancy_in_schema_table_names: TRUE/FALSE, if TRUE will add occupancy
#  into otu_attributes schema and into abund_tables table name, this is useful
#  if initially trialing different occupancy filtering in app (although will
#  likely want to go back and delete unnecessary schemas/ tables later).
#' #@param App_title: Title to appear on the app header (not linked to future web
#  address)
#' #@param Example_sequence: Sequence to be used as an example sequence for users when
#  exploring the app
#' #@param Info_text: Text to feature on app describing purpose and function
#'


#' Step 0.1
#' Merge AVC with location from CountrySide survery.
#'
#' @description
#' Merge AVC with location from CountrySide survery.
#' 
#' @details
#' Filters out rows (OTU) by occupancy level
#' - Takes in Data from a raw OTU table
#' - Filters based on occupancy level
#'
#' @param AVC_data name of input file with AVC and PH data
#' @param CS_location_data name of input file to use for Countryside survey location data
#' @param CS_AVC_combined name of output file with combined data
#' 
#' @return None
#' @examples
#'
#' merge_AVC_location_data <- function('cs_avc_ph.csv', 'cs_location.csv', 'combined.csv')
#' 
#' @note
#' This used the 10km areas from country side survey
#' 
#' @export
merge_AVC_location_data <- function( AVC_data, CS_location_data, CS_AVC_combined){
    
    cs_avc = data.table::fread(AVC_data)
    cs_location = data.table::fread(CS_location_data,)

    ## "V1" means there were row names, we need to change this to a proprt column name
    data.table::setnames(cs_avc, old = "V1", new = "ID", skip_absent = TRUE)
    
    cs_avc_with_location <- merge(cs_avc, cs_location, by = 'ID')
    data.table::fwrite(cs_avc_with_location, CS_AVC_combined)
}



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
#' Clean up enviromental metadata. 
#' 
#' @details
#' Takes in pre-prcoessed data and cleans it up:
#' - Removes unneeded fields
#' - rearranges into a format SQL is happier with
#' - Write new data to file
#'
#' @param enviroment_data enviromental data file. should contain the following column names
#'  "avc_code","avc","pH","eastings","northings, row names should be sample ID's
#'  (these should be able to be matched to Env file but don't need to be in the
#'  same order)
#' @param filtered_env_data name to use for the outputed filtered env file
#' 
#' @return None
#' @examples
#' # Create directory for supplementary SQL tables
#' #dir.create(
#' #  'Supplementary/Tables_in_SQL,
#' #  showWarnings = FALSE, 
#' #  recursive = TRUE
#' #)
#' 
#' # Clean environmental metadata
#' #clean_environmental_metadata(
#' #  params$Env_file,
#' #  paste0(Output_dir_with_occ, "/Supplementary/Tables_in_SQL/Env.csv")
#' #)
#' 
#' @note
#' Original code provided output variable: Env_for_SQL
#' @export
clean_environmental_metadata <- function(enviroment_data,  filtered_env_data){

    # This function used to filter out location data, not sure why, seems needed.
    Env=data.table::fread(enviroment_data)

    Env_for_SQL=Env[-which(is.na(Env$avc_code)),]
                                        # Rearrange slightly so suitable for inserting into
                                        # SQL- e.g make sample a column(rather than rownames)
                                        # and change "pH" colname to "ph" as standard
                                        # postgres field names without double quotes have to be
                                        # lower case
                                        # https://deeplearning.lipingyang.org/2017/01/07/
                                        # postgresql-column-names-of-a-table-are-case-sensitive/
    data.table::fwrite(Env_for_SQL, filtered_env_data)
}

## Step 1.2
#' Clean up the OTU table
#'
#' @description
#' Clean up the OTU table.
#' 
#' @details
#' Filters out rows (OTU) by occupancy level
#' - Takes in Data from a raw OTU table
#' - Filters based on occupancy level
#'
#' @param OTU_file pre-processed OTU table
#' @param filtered_OTU_file name to use for the output file (filtered OTU table)
#' @param  OTU_table_occupancy_filter occupancy level to filter the table by
#' 
#' @return None
#' @examples
#' # Create directory for supplementary SQL tables
#' #Output_dir_with_occ = 'files'
#' #dir.create(
#' #  paste0(Output_dir_with_occ, "/Supplementary/Tables_in_SQL"),
#' #  showWarnings = FALSE, 
#' #  recursive = TRUE
#' #)
#' 
#' # Clean OTU table
#' #clean_OTU_table(
#' #  'OTU_tab_file',
#' #  paste0(Output_dir_with_occ, "/Supplementary/Tables_in_SQL/OTU_abund.csv")
#' #)
#'
#' @note
#' Original code provided output variable:  OTU_tab_sub_occ_dec
#' @export
clean_OTU_table <- function( OTU_file, filtered_OTU_file, OTU_table_occupancy_filter=30){
    
                                        #read in OTU_tab
    OTU_tab=data.frame(data.table::fread( OTU_file),row.names=1,check.names=FALSE)
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
    OTU_out <- data.table::data.table(rownames(OTU_tab_sub_occ_dec), OTU_tab_sub_occ_dec)
    data.table::setnames(OTU_out, 1, "ID")
    data.table::fwrite(OTU_out,filtered_OTU_file)
}

## Step 1.3
#' Prepare abundance_stats table
#'
#' @description
#' Prepare abundance_stats table.
#' 
#' @details
#' Get individual OTU stats to summarise abundance (rank) and occupancy
#' (percentage and rank), these stats will form a table in the database (added to
#' DB in step 2). In this code chunk saved to r
#' params$Output_dir /Supplementary/Tables_in_SQL subfolder for reference.
#'
#' @param filtered_OTU_file filtered OTU table (created in step 1.2)
#' @param abundance_stats_file name to use for the output file, filtered abundace stats
#' 
#' @return None
#'
#' @examples
#' # Requires OTU_tab_sub_occ_dec variable created in another function
#' # Get abundance statistics
#' #Output_dir_with_occ = 'files'
#' #get_abundance_stats(
#' #  paste0(Output_dir_with_occ, "/Supplementary/Tables_in_SQL/OTU_abund.csv"),
#' #  paste0(Output_dir_with_occ, "/Supplementary/Tables_in_SQL/abundance_stats.csv")
#' #)
#'
#' @note
#' Original code provided output variable:  abundance_stats
#' @export
get_abundance_stats <- function(filtered_OTU_file, abundance_stats_file){

                                        # get the otu table
    OTU_tab_sub_occ_dec = data.frame(data.table::fread(filtered_OTU_file),row.names=1,check.names=FALSE)
    
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
    data.table::fwrite(abundance_stats,abundance_stats_file,row.names=FALSE)
}

## Step 1.4
#' Prepare taxonomy table
#'
#' @description
#' Prepare abundance_stats table.
#' 
#' @details
#' Read in taxonomy, split taxonomic fields by ";"and subset to only include OTUs
#' that meet occupancy threshold ,save to r params$Output_dir/Supplementary/Tables_in_SQL
#' subfolder for reference.
#'
#' @param OTU_abund_filter_file filtered OTU table created in step 1.2 for input
#' @param taxonomy_file pre-processed taxonomy table to use as input
#' @param filtered_taxonomy_file name for output file, which is a filtered taxonomy table
#' 
#' @return None
#'
#' @examples
#' # Prepare taxonomy table
#' #Output_dir_with_occ = 'files'
#' #prepair_taxonomy_table(
#' #  OTU_abund_filter_file = paste0(Output_dir_with_occ, "/Supplementary/Tables_in_SQL/OTU_abund.csv")#,
#'  # input_file = 'Tax_file',
#'  # output_file = paste0(Output_dir_with_occ, "/Supplementary/Tables_in_SQL/Taxonomy.csv")
#' #)
#'
#' @note
#' Original code provided output variable:
#' Taxonomy_filt (Taxonomy_Sort is now written as extra steps were needed.)
#' 
#' @export
prepair_taxonomy_table <- function(taxonomy_file, filtered_taxonomy_file, OTU_abund_filter_file){
    
    Taxonomy <- read.csv(taxonomy_file)
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
    data.table::fwrite(Taxonomy_Sort,filtered_taxonomy_file)
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
#' Put the filtered OTU, taxonomy, and abundance files into RSQLite files.
#' 
#' @details
#' Reads in the filtered OTU, taxonomy, and abundance files, and creates corresponding
#' tables in an RSQLite file
#'
#' @param abundance_csv filtered OTU table
#' @param taxonomy_csv filtered abundance table
#' @param otu_csv filtered taxonomy table
#' 
#' @return None
#'
#' @examples
#' # Prepare taxonomy table
#' #Output_dir_with_occ = 'files'
#' #prepair_taxonomy_table(
#' #  OTU_abund_filter_file = paste0(Output_dir_with_occ, "/Supplementary/Tables_in_SQL/OTU_abund.csv"),
#'#   input_file = 'Tax_file.csv',
#' #  output_file = paste0(Output_dir_with_occ, "/Supplementary/Tables_in_SQL/Taxonomy.csv")
#' #)
#'
#' @note
#'
#' format_otu_for_Rsqlite(filtered_abundance_csv, filtered_taxonomy_csv, filtered_otu_csv)
#' 
#' 
#' @export
format_otu_for_Rsqlite <- function(filtered_abundance_csv, filtered_taxonomy_csv, filtered_otu_csv){

    print("Debug: Load files")
    
    abundance_csv <- read.csv(filtered_abundance_csv)
    taxonomy_csv <- read.csv(filtered_taxonomy_csv)
    otu_csv <- read.csv(filtered_otu_csv)
    
    ## Lets try not splitting to begin with
    ## Current SQL command is:
    ## CREATE TABLE IF NOT EXISTS abund_tables.<prefix>_abund_2(hit character varying(30), "<column names each followed by '" numeric, "'>" numeric, CONSTRAINT '<prefix>_abund_2_pkey PRIMARY KEY("hit"))'
    ## Lets get rid of the prefix and double tables to:
    ## CREATE TABLE IF NOT EXISTS abund_tables(hit character varying(30), "<column names each followed by '" numeric, "'>" numeric, CONSTRAINT '<prefix>_abund_2_pkey PRIMARY KEY("hit"))'

    ## Lets use sprintf to get better formatting of the string.
    ## Fails with "CONSTRAINT", not sure on need, so removing.


    sql_command <- sprintf("create table if not exists abund_table (hit character varying(30), %s numeric, primary key (hit))",
                                     paste(colnames(abundance_csv)[-1], collapse=' numeric, ')   #colnames
                           )

    ## Create table, abundance
    abundance_db <- DBI::dbConnect(RSQLite::SQLite(), "abundance_db.sqlite")
    DBI::dbExecute(conn = abundance_db, statement = sql_command)
    ## Fill table
    DBI::dbWriteTable(abundance_db, "abund_table", abundance_csv, append = TRUE, row.names = FALSE)
    ## Disconnect, best practise?
    DBI::dbDisconnect(abundance_db)


    ## Create table, taxonomy

    ## "order" is an SQL command so need to escape it in quotes, easier to escape all
    ## fields in quotes
    
    sql_command_taxonomy <- sprintf("create table taxonomy_table (hit character varying (30), %s character varying (250), primary key (hit))",
                                    paste('"',colnames(abundance_csv)[-1],'"',
                                          collapse=' charater varying(250),',
                                          sep='')
                                    )

    ## Create table, taxonomy
    taxonomy_db <- DBI::dbConnect(RSQLite::SQLite(), "taxonomy_db.sqlite")
    DBI::dbExecute(conn = taxonomy_db, statement = sql_command_taxonomy)
    ## Fill table
    DBI::dbWriteTable(taxonomy_db, "taxonomy_table", taxonomy_csv, append = TRUE, row.names = FALSE)
    ## Disconnect, best practise?
    DBI::dbDisconnect(taxonomy_db)



    ## OTU table has too many columns for SQL/SQLite so need to transpose
    transpose_otu_csv = t(otu_csv)
    ## convert back to a data.frame
    otu_csv = data.frame(hit=row.names(transpose_otu_csv), transpose_otu_csv, check.names=FALSE)

    ## Create OTU table
    sql_command_otu <- sprintf("create table otu_table (hit character varying (30), %s character varying (30), primary key (hit))",
                                    paste('"',colnames(otu_csv)[-1],'"',
                                          collapse=' charater varying(250),',
                                          sep='')
                                    )

    ## Create table, otu
    otu_db <- DBI::dbConnect(RSQLite::SQLite(), "otu_db.sqlite")
    DBI::dbExecute(conn = otu_db, statement = sql_command_otu)
    ## Fill table
    DBI::dbWriteTable(otu_db, "otu_table", otu_csv, append = TRUE, row.names = FALSE)
    ## Disconnect, best practise?
    DBI::dbDisconnect(otu_db)
    

    ## Create maps table
    ## Postgresql use bytea as the data type, but sqlite does not have this
    ## equvilant is "blob" so changing bytea to blob for this.
    ## https://jfaganuk.github.io/2015/01/12/storing-r-objects-in-sqlite-tables/
    sql_command_maps <- sprintf("create table otu_attributes_table (otu_name character varying (30), map_object blob, primary key (otu_name))")

    ## Create table, maps
    maps_db <- DBI::dbConnect(RSQLite::SQLite(), "maps_db.sqlite")
    DBI::dbExecute(conn = maps_db, statement = sql_command_maps)

    #Fill maps later

    ## Fill table
    ## DBI::dbWriteTable(maps_db, "maps_table", maps_csv, append = TRUE, row.names = FALSE)
    ## Disconnect, best practise?
    DBI::dbDisconnect(maps_db)

    ## will fill maps table in step 3 using save_otu_map function
   
}



## Step 3 Make map objects for DB

### 3.1 Map preparation

# **r**
# ```{R map prep, eval=FALSE}



## Step 3.1 Map preperation
#' Put the filtered OTU, taxonomy, and abundance files into RSQLite files
#'
#' @description
#' Step 3 Make map objects for DB.
#' 
#' @details
#' Step 3 Make map objects for DB
#'
#' @param ukcoast_poly shape file of uk polygon
#' @param ukcoast_line shape file of uk outline
#' @param uk_poly_converted shape file of uk polygon
#' @param uk_line_converted shape file of uk outline
#' @param uk_grid_converted shape file of uk outline
#'
#' @return None
#'
#' @examples
#' # Prepare map
#' #map_prep(
#' #  otu_table = "supplementary/Tables_in_SQL/OTU_abund.csv",
#' #  enviroment_data='country_side_survey/CS2007_Env_pH_AVC.csv'
#' #  ukcoast_poly = 'ukcoast1.shp',
#' #  ukcoast_line = 'ukcoast_line.shp'
#' #)
#' 
#' @note
#'
#' Running this as separate function, but perhaps this should be the start of another function?
#' 
#' abund_table is called by env name: OTU_tab_sub_occ_dec
#' looks like this sets up a series of functions and variables to be used in a latter step
#' to make a db. This last part maybe much harder to split up.   Need to focus on first
#' parts and leave this until later. Perhaps it would be better to move to SQLite????
#' Perhaps some of the functions to setup, like grid bit
#'
#' @export
map_prep <- function(ukcoast_poly,
                     ukcoast_line,
                     uk_poly_converted,
                     uk_line_converted,
                     uk_grid_converted){
    
    ## Explicitly define CRS
    ukgrid <- sf::st_crs(27700)  # British National Grid
    latlong <- sf::st_crs(4326) # WGS84 Lat/Long
    
    ## Read shapefiles, specifying CRS if needed
    uk_poly <- sf::st_read(ukcoast_poly)
    uk_poly <- sf::st_set_crs(uk_poly, latlong)
    uk_poly <- sf::st_transform(uk_poly, ukgrid)
    
    uk_line <- sf::st_read(ukcoast_line)
    uk_line <- sf::st_set_crs(uk_line, latlong)
    uk_line <- sf::st_transform(uk_line, ukgrid)
  
    ## Create interpolation grid based on UK polygon bounding box
    uk_poly_bbox <- sf::st_bbox(uk_poly)
    x_range <- c(uk_poly_bbox["xmin"], uk_poly_bbox["xmax"])
    y_range <- c(uk_poly_bbox["ymin"], uk_poly_bbox["ymax"])
  
    ## Create grid for interpolation
    grd <- expand.grid(
        x = seq(from = x_range[1], to = x_range[2], by = 5000),
        y = seq(from = y_range[1], to = y_range[2], by = 5000)
    )
  
    ## Convert grid to SpatVector
    grd_vect <- terra::vect(grd, geom = c("x", "y"), crs = ukgrid$input)
  
    ## Optional: Visualization
    terra::plot(terra::vect(uk_poly))
    terra::plot(terra::vect(uk_line), add = TRUE, col = "red")
    terra::plot(grd_vect, add = TRUE, col = "blue", pch = ".")
    
    ## Return list of spatial objects
#    return(list(
#        uk_poly = uk.poly, 
#        uk_line = uk.line, 
#        grid = grd_vect
#    ))

    sf::st_write(uk_poly, uk_poly_converted, driver = "ESRI Shapefile", overwrite = TRUE)
    sf::st_write(uk_line, uk_line_converted, driver = "ESRI Shapefile", overwrite = TRUE)
    terra::writeVector(grd_vect, uk_grid_converted)
}


### 3.2 Generate maps per OTU
#' Generate maps per OTU
#' 
#' @description
#' creates map object representing an individual.
#' OTUs/ ASVs geographical distribution.
#' 
#' @details
#' **save_otu_map** function  creates map object representing an individual
#' OTUs/ ASVs geographical distribution this is saved to both the specified
#' outputdir and database. Pngs of plotted map objects can also be saved to
#' outdir for reference.
#'
#' @param OTU_name Not sure what this is.
#' @param OTU_table filtered table
#' @param Env_table Env_sub #created in map_prep
#' @param Grid created in map_prep
#' @param UK_poly uk.poly #created in map_prep
#' @param UK_line uk.line #created in map_prep
#' @param Conn SQLdb, removed so correct
#' @param Schema_table_prefix removed for time being
#' @param Output_dir not used #abandoned in favour of dir in string
#' @param Make_png flag
#' #@param otu_table filtered OTU table (created in step 1.2)
#' #@param ukcoast_poly shape file of uk polygon
#' #@param ukcoast_line shape file of uk outline
#'
#' @return None
#'
#' @examples
#'
#' # Prepare map
#' #map_prep(
#' #  otu_table = "/Supplementary/Tables_in_SQL/OTU_abund.csv",
#' #  ukcoast_poly = 'ukcoast1.shp',
#' #  ukcoast_line = 'ukcoast_line.shp'
#' #)
#' 
#' @note
#' I think this might need merging with another function, can't remember
#' which or why
#'
#' @export
save_otu_map <- function(OTU_name,  # Not sure what this is.
                         OTU_table, # filtered table
                         Env_table, # Env_sub #created in map_prep
                         Grid,      # created in map_prep
                         UK_poly,   # uk.poly #created in map_prep
                         UK_line,   # uk.line #created in map_prep
                         Conn,      # SQLdb, removed so correct
                         Schema_table_prefix, # removed for time being
                         Output_dir, # again, abandoned in favour of dir in string
                         Make_png # flag
                         ){


                                        #read in OTU_tab
    OTU_table = data.frame(data.table::fread(otu_table),
                           row.names=1,
                           check.names=FALSE)

                                        #read in enviromental data
    env_data = data.frame(data.table::fread(enviroment_data),
                           row.names=1,
                           check.names=FALSE)
                                        # first lets get env in the same order as otu table

                                        # Env is Env.csv file read in as a dataframe which is csv
    
    Env_sub = env_data[row.names(OTU_table),]
                                        # Isnt this just a "test" so does nothing?
    identical(row.names(Env_sub),row.names(OTU_table))

    
    ## drop=FALSE protects against conversion to vector instead of dataframe
    ## OTU_name is column selection of some kind, need to find definition and
    ## reasoning
    otu_abund <- OTU_table[,OTU_name,drop=FALSE]
    
                                        #make dataframe with eastings and northings
                                        # binds to otu_abund as well?
    dat <- cbind(Env_table[,c('eastings','northings')],otu_abund)
    colnames(dat)[3] = "OTU"
    dat <- dat[complete.cases(dat), ]
    attach(dat)
                                        #specify the coordinates for the file with otu
    
    ## could we expans this to:
    ## sp::coordinates(dat) = ~dat$eastings+dat$northings  as attach puts columns in namespace?
    sp::coordinates(dat) = ~eastings+northings
                                        #interpolate
    
    spc.idw <- gstat::krige(OTU~1,dat,Grid)
    ukgrid = "+init=epsg:27700"
    spc.idw@proj4string <- sp::CRS(ukgrid)
                                        #now need to only show interpolation within UK coastline
    overlay.idw <- sp::over(spc.idw,UK_poly)
    newmap = spc.idw[complete.cases(overlay.idw),]
                                        #set heat map scale, want this to be dependent on OTU abundance i.e we want to show abundance contrast relative to that OTU
    minv = min(otu_abund)
    maxv = max(otu_abund)
    at = c(minv,
           signif(maxv/256,1),
           signif(maxv/128,1),
           signif(maxv/64,1),
           signif(maxv/32,1),
           signif(maxv/16,1),
           signif(maxv/8,),
           signif(maxv/4,2),
           signif(maxv/2,2),
           maxv
           )
                                        #save map object to outdir
    mapandinfo <- c(newmap,at)
    
    save(mapandinfo, file=paste(Output_dir,"/",OTU_name,".RData",sep=""))
    
                                        #Add object to database  
                                        #need to first get object into a form suitable for SQL
                                        #convert to stream of bytes
    
    ## Serialise converts all pointers and links to objects into one object
    ## so that it can be saved or transferred
    ser_mapandinfo = serialize(mapandinfo,connection=NULL,ascii=TRUE)
    
    ## will fill maps table in step 3 using save_otu_map function
    
    sql_command_maps <- sprintf("insert into maps_table values (?, ?)")
    DBI::dbExecute(maps_db, sql_command_maps, params = list(OTU_name, ser_mapandinfo))
    DBI::dbDisconnect(maps_db)
    
                   if (Make_png==TRUE){ 
                       map_plot = sp::spplot(newmap["var1.pred"],
                                             at = at,
                                             sp.layout=list("sp.lines",UK_line,col="black",lwd=2)
                                             )
                       
                       dir.create(paste0(Output_dir,"/png"),
                                  showWarnings = FALSE,
                                  recursive=TRUE)
                       
                       png(file=paste0(Output_dir,"/png/",OTU_name,"_plot.png"),
                           width=200, height=400)
                       print(map_plot)
                       dev.off()
                   }
}
    
    
    

### 3.2 Generate maps per OTU
#' Generate maps per OTU
#' 
#' @description
#' Using rparallel to parallelise -need to reauthenticate database on all cpus
#' working on.
#'
#' @details
#' Runs save_otu_map function on all OTUs.
#' Using rparallel to parallelise -need to reauthenticate database on all cpus
#' working on.
#'
#' @param Output_dir_with_occ blank
#' @param OTU_tab_sub_occ_dec blank
#' @param Env_sub blank
#' @param grd blank
#' @param uk.poly blank
#' @param uk.line blank
#' @param Schema_table_prefix_modified blank
#'
#' @return None
#'
#' @examples
#' # Prepare map
#' #map_prep(
#' #  otu_table = "/Supplementary/Tables_in_SQL/OTU_abund.csv",
#' #  ukcoast_poly = 'ukcoast1.shp',
#' #  ukcoast_line = 'ukcoast_line.shp'
#' #)
#' 
#' @note
#' I think this might need merging with another function, can't remember
#' which or why
#'
#' @export
maps_parallelise <- function(Output_dir_with_occ,
                             OTU_tab_sub_occ_dec,
                             Env_sub,
                             grd,
                             uk.poly,
                             uk.line,
                             Schema_table_prefix_modified
                             ){
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
                              password = Sys.getenv('SQL_PWD')
                              )
    }
    )
    
    
    
    parSapply(cl,colnames(OTU_tab_sub_occ_dec), function(x) save_otu_map(OTU_name=x,OTU_table=OTU_tab_sub_occ_dec,Env_table=Env_sub,Grid=grd,UK_poly=uk.poly,UK_line=uk.line,Conn=conn,Schema_table_prefix=Schema_table_prefix_modified,Output_dir=paste0(Output_dir_with_occ,"/Supplementary/Map_objects"),Make_png=FALSE))
    stopCluster(cl)
    
    ## ```
}



## Step 3 Make Blast database
###  3.3 Filter fasta file to contain OTU sequences that meet occupancy filter

## Need to filter sequences prior to making blast DB otherwise we will have hits
## in BLast DB with no supplementary info (i.e we dont want to keep sequences
## corresponding to OTU's/ ASV's that did not meet occupancy filter)
## Doing this using biopython.


# make_blast_py <- function(){
#    ## **python**
    ## ```{python filter fasta file, eval=FALSE}
#    from Bio import SeqIO
#    import os
    ## get OTUs we want to keep from OTU table
    ## python can access r variables within markdown and converts our r dataframe into a dictionary- the keys are the OTU names
#    OTU_tab_dict = r.OTU_tab_sub_occ_dec
    ## make output dir for filtered fasta
#    if not os.path.exists(r.Output_dir_with_occ+"/Supplementary/Filtered_sequences"):
#               os.makedirs(r.Output_dir_with_occ+"/Supplementary/Filtered_sequences")
    ## only keep records in filtered otu tab (dictionary keys) 
    ## see http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec372 for reference
#    records = (record for record in SeqIO.parse(r.params["Fasta_file"],"fasta") if record.id in OTU_tab_dict.keys())
#    SeqIO.write(records, r.Output_dir_with_occ+"/Supplementary/Filtered_sequences/filtered_sequences.fasta", "fasta")
    ## ```
#}

###  3.4 Make blast database
# Take filtered fasta and make blast database for back end of shiny app

#make_blast_bash <- function(){
    ## **Bash:**
    ## ```{bash BlastDB, eval=FALSE}
    ## cmake makes dir and any parent dirs necessary
#    cmake -E make_directory $Output_dir_with_occ/App/Blast_DB
#    makeblastdb -in $Output_dir_with_occ/Supplementary/Filtered_sequences/filtered_sequences.fasta -dbtype nucl -out $Output_dir_with_occ/App/Blast_DB/Blast_DB
    ## ```
#}




## create_app_python <-function(){
## #### 4 Create app front-end from template
    
##     ## Using python to edit r code - find python easiest option for file handling
    
##     ## **Python:**
##     ## ```{python create app, eval=FALSE}
##     import re
##     ## get blank template path using R App_template parameter  
##     blank_template=r.params["App_template_input_dir"]+"/Blank_taxonomic_explorer_app_with_functional_place_holders.R"
##     ## open blank_template
##     with open(blank_template,'r') as blank_template_file:
##                                          ## open output file  
##                                          with open(r.Output_dir_with_occ+"/App/App.R",'w') as customised_file:
## ### loop over lines 
##     for line in blank_template_file:
##                     ## empty array for any placeholder strings that need to be replaced within that line 
##                     occurence_strings=[]
##     ## get all indices where placeholder  starts("<--")
##     occurence_starts= [i.start() for i in re.finditer("<--",line)]
##     ## if there are any beginnings of placeholders in the line  continue   
##     if (len(occurence_starts)!=0):
##         ## get all indices  of the end of place holder/s "-->"  
##         occurence_ends= [i.start() for i in re.finditer("-->",line)]
##     ## loop over number of start indices per line (as may be multiple placeholders per line)     
##     for i in range(len(occurence_starts)):
##                  ## get whole place holder string  ,occurence_ends are the start indices of "-->" string so add 3
##                  occurence_string=line[occurence_starts[i]:occurence_ends[i]+3]
##     ## if not in occurrence_strings array already append
##     if (occurence_string not in occurence_strings):
##         occurence_strings.append(occurence_string)
##     ## define new line        
##     new_line=line
##     ## for every unique place holder string replace placeholder with relevant r parameter(exception being Schema_table_prefix_modified as it is not a parameter)   
##     for occurence_string in occurence_strings:
##                                 parameter=occurence_string.split("<--specified_")[1].split("-->")[0]
##     if(parameter=="Schema_table_prefix"):
##         new_line=new_line.replace(occurence_string,r.Schema_table_prefix_modified)
##     else:
##         new_line=new_line.replace(occurence_string,r.params[parameter])
##     ## ok lets write the new modified line to our output file       
##     customised_file.write(new_line)
##     ## if no placeholder strings write original line to output file        
##     else:  customised_file.write(line)      
    
##     ##  ```
##     ##  App should now be an executable 
## }
