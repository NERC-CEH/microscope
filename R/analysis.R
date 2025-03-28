#' Microscope: Generate taxonomic explorer app inputs
#'
#' @author "Briony Jones"
#'
#' @import data.table
#' @import vegan
#' @import terra
#' @import ggplot2
#' @import gstat
#' @import parallel
#' @import splitstackshape
#' @import RSQLite
#' @import DBI
#' @import sp
#' @import sf
#'
#' @description This script automates the process of generating a taxonomic explorer app 
#'   that links marker gene sequences to environmental responses. The application is similar to
#'   ID-TaxER (https://shiny-apps.ceh.ac.uk/ID-TaxER/). The script creates an SQL database 
#'   and a BLAST database for the app backend, and modifies an R Shiny template file to produce 
#'   the front end. Preprocessed tables (OTU/taxonomy) and map objects are also saved locally 
#'   for reference. Currently, this script is suitable for UK data only as it uses UK map outline
#'   objects. It should work for all CS (Countryside Survey) molecular datasets. Multiple 
#'   taxonomic datasets can be added to the same database if they share the same environmental data.
#'
#' @note Due to the complexity of this script, it's recommended to run it chunk by chunk
#'   rather than knitting the entire file at once, especially during testing.

#' Step 0.1
#' Merge AVC with location from CountrySide survey
#'
#' @description Combines AVC (Aggregate Vegetation Class) data with location information from
#'   the Countryside survey.
#' 
#' @param AVC_data Character string. Path to input file with AVC and pH data.
#' @param CS_location_data Character string. Path to input file with Countryside survey location data.
#' @param CS_AVC_combined Character string. Path where the combined output file will be saved.
#' 
#' @return None. The function writes the combined data to the specified output file.
#' 
#' @examples
#' merge_AVC_location_data("data/cs_avc_ph.csv", "data/cs_location.csv", "data/combined.csv")
#' 
#' @note This uses the 10km areas from Countryside survey.
#' 
#' @export
merge_AVC_location_data <- function(AVC_data, CS_location_data, CS_AVC_combined) {
    cs_avc <- data.table::fread(AVC_data)
    cs_location <- data.table::fread(CS_location_data)

    # "V1" means there were row names, we need to change this to a proper column name
    data.table::setnames(cs_avc, old = "V1", new = "ID", skip_absent = TRUE)
    
    cs_avc_with_location <- merge(cs_avc, cs_location, by = 'ID')
    data.table::fwrite(cs_avc_with_location, CS_AVC_combined)
}



## Step 1 Prepare tables for database

## Step 1.1
#' Clean up enviromental metadata
#'
#' @description Processes and formats environmental metadata for use in the taxonomic explorer app.
#' 
#' @details The function removes unnecessary fields and rearranges the data into a format 
#'   suitable for SQL. It removes samples with missing AVC codes.
#'
#' @param enviroment_data Character string. Path to environmental data file. Should contain 
#'   columns "avc_code", "avc", "pH", "eastings", "northings". Row names should be sample IDs
#'   that match the OTU table.
#' @param filtered_env_data Character string. Path where the filtered environmental data will be saved.
#' 
#' @return None. The function writes the filtered data to the specified output file.
#' 
#' @examples
#' # Create directory for supplementary SQL tables
#' dir.create("output/Supplementary/Tables_in_SQL", showWarnings = FALSE, recursive = TRUE)
#' 
#' # Clean environmental metadata
#' clean_environmental_metadata(
#'   "data/environmental_data.csv",
#'   "output/Supplementary/Tables_in_SQL/Env.csv"
#' )
#' 
#' @note
#' Original code provided output variable: Env_for_SQL
#' @export
clean_environmental_metadata <- function(enviroment_data, filtered_env_data) {
    # Read the environmental data
    Env <- data.table::fread(enviroment_data)

    # Remove rows with missing AVC codes
    Env_for_SQL <- Env[-which(is.na(Env$avc_code)), ]
    
    # Write filtered data to file
    data.table::fwrite(Env_for_SQL, filtered_env_data)
}

## Step 1.2
#' Clean up the OTU table
#'
#' @description Filters and normalizes an OTU (Operational Taxonomic Unit) table.
#' 
#' @details The function:
#'   1. Removes samples with low read counts (< 5000)
#'   2. Filters out OTUs that don't meet the specified occupancy threshold
#'   3. Normalizes the OTU table using total sum scaling
#'
#' @param OTU_file Character string. Path to the pre-processed OTU table file.
#' @param filtered_OTU_file Character string. Path where the filtered OTU table will be saved.
#' @param OTU_table_occupancy_filter Integer. Minimum number of samples an OTU must be present 
#'   in to be retained. Default is 30.
#' 
#' @return None. The function writes the filtered and normalized OTU table to the specified output file.
#'
#' @examples
#' # Clean OTU table with default occupancy filter (30)
#' clean_OTU_table(
#'   "data/raw_otu_table.csv",
#'   "output/Supplementary/Tables_in_SQL/OTU_abund.csv"
#' )
#'
#' # With custom occupancy filter
#' clean_OTU_table(
#'   "data/raw_otu_table.csv",
#'   "output/Supplementary/Tables_in_SQL/OTU_abund.csv",
#'   OTU_table_occupancy_filter = 20
#' )
#'
#' @export
clean_OTU_table <- function(OTU_file, filtered_OTU_file, OTU_table_occupancy_filter = 30, sample_read_number = 5000) {
    # Read in OTU table
    OTU_tab <- data.frame(data.table::fread(OTU_file), row.names = 1, check.names = FALSE)
    
    # Remove samples with reads less than 5000
    OTU_tab_sub <- OTU_tab[rowSums(OTU_tab) > sample_read_number, ]
    
    # Convert to presence/absence to filter by occupancy
    OTU_tab_sub_pa <- (OTU_tab_sub != 0) * 1
    
    # Remove taxa that do not meet occupancy threshold
    OTU_tab_sub_occ <- OTU_tab_sub[, which(colSums(OTU_tab_sub_pa) >= OTU_table_occupancy_filter)]
    
    # Normalize OTU table using total sum scaling
    OTU_tab_sub_occ_dec <- vegan::decostand(OTU_tab_sub_occ, method = "total")
    
    # Format and write filtered OTU table
    OTU_out <- data.table::data.table(rownames(OTU_tab_sub_occ_dec), OTU_tab_sub_occ_dec)
    data.table::setnames(OTU_out, 1, "ID")
    data.table::fwrite(OTU_out, filtered_OTU_file)
}

## Step 1.3
#' Generate abundance statistics for OTUs
#'
#' @description Calculates summary statistics for each OTU to characterize its abundance and occupancy.
#' 
#' @details The function computes:
#'   1. Total abundance of each OTU across all samples
#'   2. Abundance rank (where ties are given the same rank)
#'   3. Occupancy (number of samples an OTU is present in)
#'   4. Occupancy proportion as a percentage with rank
#'
#' @param filtered_OTU_file Character string. Path to the filtered OTU table file (created in step 1.2).
#' @param abundance_stats_file Character string. Path where the abundance statistics will be saved.
#' 
#' @return None. The function writes the abundance statistics to the specified output file.
#'
#' @examples
#' # Generate abundance statistics
#' get_abundance_stats(
#'   "output/Supplementary/Tables_in_SQL/OTU_abund.csv",
#'   "output/Supplementary/Tables_in_SQL/abundance_stats.csv"
#' )
#'
#' @note
#' Original code provided output variable:  abundance_stats
#' @export
get_abundance_stats <- function(filtered_OTU_file, abundance_stats_file) {
    # Read the filtered OTU table
    OTU_tab_sub_occ_dec <- data.frame(data.table::fread(filtered_OTU_file), row.names = 1, check.names = FALSE)
    
    # Get total abundance for each OTU
    abundance_stats <- data.frame(
        hit = colnames(OTU_tab_sub_occ_dec),
        abundance = colSums(OTU_tab_sub_occ_dec)
    )
    
    # Calculate abundance rank (ties get the same rank)
    abundance_stats$abundance_rank <- paste(
        rank(-abundance_stats$abundance, ties.method = "min"),
        ncol(OTU_tab_sub_occ_dec),
        sep = "/"
    )
    
    # Convert to presence/absence matrix
    OTU_tab_sub_occ_dec_pa <- (OTU_tab_sub_occ_dec != 0) * 1
    
    # Calculate occupancy
    abundance_stats$occupancy <- colSums(OTU_tab_sub_occ_dec_pa)
    
    # Calculate occupancy as percentage with rank
    abundance_stats$occupancy_proportion <- paste(
        round(abundance_stats$occupancy / nrow(OTU_tab_sub_occ_dec_pa) * 100, 2),
        "% (Rank: ",
        rank(-round(abundance_stats$occupancy / nrow(OTU_tab_sub_occ_dec_pa) * 100, 2), ties.method = "min"),
        "/",
        ncol(OTU_tab_sub_occ_dec_pa),
        ")",
        sep = ""
    )
    
    # Remove unnecessary columns
    abundance_stats <- abundance_stats[, -c(2, 4)]
    
    # Write to file
    data.table::fwrite(abundance_stats, abundance_stats_file, row.names = FALSE)
}

## Step 1.4
#' Prepare taxonomy table
#'
#' @description Processes and formats taxonomy data for the taxonomic explorer app.
#' 
#' @details The function:
#'   1. Reads in the taxonomy file
#'   2. Splits taxonomic classifications by ";" into separate columns
#'   3. Filters to include only OTUs that meet the occupancy threshold
#'   4. Formats taxonomy information by removing prefixes
#'
#' @param taxonomy_file Character string. Path to the pre-processed taxonomy file.
#' @param filtered_taxonomy_file Character string. Path where the filtered taxonomy data will be saved.
#' @param OTU_abund_filter_file Character string. Path to the filtered OTU table file to use for filtering.
#' 
#' @return None. The function writes the filtered taxonomy data to the specified output file.
#'
#' @examples
#' # Prepare taxonomy table
#' prepair_taxonomy_table(
#'   "data/raw_taxonomy.csv",
#'   "output/Supplementary/Tables_in_SQL/Taxonomy.csv",
#'   "output/Supplementary/Tables_in_SQL/OTU_abund.csv"
#' )
#'
#' @note
#' Original code provided output variable:
#' Taxonomy_filt (Taxonomy_Sort is now written as extra steps were needed.)
#' 
#' @export
prepair_taxonomy_table <- function(taxonomy_file, filtered_taxonomy_file, OTU_abund_filter_file) {
    # Read taxonomy file
    Taxonomy <- data.table::fread(taxonomy_file)
    
    # Split taxonomic classifications into separate columns
    Taxonomy <- splitstackshape::cSplit(indt = Taxonomy, splitCols = 2, sep = ";")
    
    # Rename taxonomy columns to the correct levels
    num_columns = length((Taxonomy))
    if (num_columns == 9){
        colnames(Taxonomy) <- c("hit", "Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species")
    } else if (num_columns == 8){
        colnames(Taxonomy) <- c("hit", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")                                                                                                        }

    row.names(Taxonomy) <- Taxonomy$hit
    
    # Filter to match OTU table
    OTU_abund_filter <- colnames(data.table::fread(OTU_abund_filter_file))
    Taxonomy_filt <- Taxonomy[hit %in% OTU_abund_filter]
    
    # Sort and format taxonomy
    Taxonomy_Sort <- Taxonomy_filt[order(row.names(Taxonomy_filt)), ]
    Taxonomy_Sort <- Taxonomy_Sort[, lapply(.SD, function(x) sub(".*__", "", x))]
    Taxonomy_Sort[is.na(Taxonomy_Sort)] <- ""
    
    # Write to file
    data.table::fwrite(Taxonomy_Sort, filtered_taxonomy_file)
}

#' Create SQLite database tables for OTU, taxonomy, and abundance data
#'
#' @description Generates SQLite database files and tables to store the filtered OTU, 
#'   taxonomy, and abundance data for the app backend.
#' 
#' @details The function:
#'   1. Creates separate SQLite database files for abundance, taxonomy, OTU, and map data
#'   2. Defines appropriate schema for each table
#'   3. Populates the tables with the filtered data
#'
#' @param filtered_abundance_csv Character string. Path to the filtered abundance statistics file.
#' @param filtered_taxonomy_csv Character string. Path to the filtered taxonomy file.
#' @param filtered_otu_csv Character string. Path to the filtered OTU table file.
#' @param ukcoast_line_shp Character string. Path to the uk coast shapefile.
#' 
#' @return None. The function creates SQLite database files with the appropriate tables.
#'
#' @examples
#' format_otu_for_Rsqlite(
#'   "output/Supplementary/Tables_in_SQL/abundance_stats.csv",
#'   "output/Supplementary/Tables_in_SQL/Taxonomy.csv",
#'   "output/Supplementary/Tables_in_SQL/OTU_abund.csv"
#'   "data/00_raw_data/uk_map_objs/ukcoast_line.shp"
#' )
#' 
#' @export
format_otu_for_Rsqlite <- function(filtered_abundance_csv, filtered_taxonomy_csv, filtered_otu_csv, ukcoast_line_shp, filtered_environmental_csv,
                                   molecular_db = "molecular_db.sqlite",
                                   maps_db = "maps_db.sqlite") {
    
    # Read the input files
    environmental_csv <- read.csv(filtered_environmental_csv)
    abundance_csv <- read.csv(filtered_abundance_csv)
    taxonomy_csv <- read.csv(filtered_taxonomy_csv)
    otu_csv <- read.csv(filtered_otu_csv)

    # Connect to environmental database and create table
    conn_molecular_db <- DBI::dbConnect(RSQLite::SQLite(), molecular_db)

    # SQL command to create table
    sql_command_env_table <- sprintf(
        "CREATE TABLE IF NOT EXISTS env_table (hit character varying(250), avc_code numeric, avc character varying(250), pH numeric, primary key (hit))"
    )

    DBI::dbExecute(conn = conn_molecular_db, statement = sql_command_env_table)

    #ensure first column is called "hit"
    colnames(environmental_csv)[1] <- "hit"

    # Fill table and disconnect
    DBI::dbWriteTable(conn_molecular_db, "env_table", environmental_csv[1:4], row.names = FALSE, overwrite = TRUE)

    # Create abundance table
    sql_command_abundance_table <- sprintf(
        "CREATE TABLE IF NOT EXISTS abund_table (hit character varying(30), %s numeric, primary key (hit))",
        paste(colnames(abundance_csv)[-1], collapse = ' numeric, ')
    )

    # Fill table and disconnect
    DBI::dbWriteTable(conn_molecular_db, "abund_table", abundance_csv, row.names = FALSE, overwrite = TRUE)

    
    ## Create table
    DBI::dbExecute(conn = conn_molecular_db, statement = sql_command_abundance_table)

    ## Create taxonomy table
    sql_command_taxonomy_table <- sprintf(
        "CREATE TABLE IF NOT EXISTS taxonomy_table (hit character varying (30), %s character varying (250), primary key (hit))",
        paste(
            '"', colnames(taxonomy_csv)[-1], '"',
            collapse = ' charater varying(250),',
            sep = ''
        )
    )
    
    # Connect to taxonomy database and create table
    DBI::dbExecute(conn = conn_molecular_db, statement = sql_command_taxonomy_table)
    
    # Fill table and disconnect
    DBI::dbWriteTable(conn_molecular_db, "taxonomy_table", taxonomy_csv, row.names = FALSE, overwrite = TRUE)
  
    ## Create OTU table
    sql_command_otu_table <- sprintf(
        "CREATE TABLE IF NOT EXISTS otu_table (hit character varying (30), %s character varying (30), primary key (hit))",
        paste(
            '"', colnames(otu_csv)[-1], '"',
            collapse = ' charater varying(250),',
            sep = ''
        )
    )

    ## OTU table preparation (transpose due to column limit)
    rownames(otu_csv) <- otu_csv[, 1]
    transpose_otu_csv <- t(otu_csv)
    transpose_otu_csv = as.data.frame(transpose_otu_csv, stringsAsFactors = FALSE)
    otu_csv <- data.frame(hit = row.names(transpose_otu_csv), transpose_otu_csv, check.names = FALSE)

    DBI::dbExecute(conn = conn_molecular_db, statement = sql_command_otu_table)

    
    ## Fill table
    DBI::dbWriteTable(conn_molecular_db, "otu_table", otu_csv, row.names = FALSE, overwrite = TRUE)


    ## Disconnect from database
    DBI::dbDisconnect(conn_molecular_db)


    #turns out that varchar is ignored by SQLlite, it only uses text.  Can be great to indicate to developers that we want short text though.
    # Create maps table
    sql_command_maps <- "CREATE TABLE IF NOT EXISTS maps_table (otu_name character varying (30), map_object blob, break_points text, primary key (otu_name))"
    
    # Connect to maps database and create table
    maps_db_conn <- DBI::dbConnect(RSQLite::SQLite(), maps_db)
    DBI::dbExecute(conn = maps_db_conn, statement = sql_command_maps)
    DBI::dbDisconnect(maps_db_conn)
   
    ## Maps table will be filled in step 3 using save_otu_map function


    conn <- DBI::dbConnect(RSQLite::SQLite(), maps_db)

    ## Create plotting_tools table if it doesn't already exist
    ## Note: SQLite doesn't have schemas like PostgreSQL, so we'll use table naming convention
    DBI::dbExecute(conn, "CREATE TABLE IF NOT EXISTS plotting_tools_map_tools (
           description TEXT PRIMARY KEY,
           plot_object BLOB
          );")
    
    ## Read in UK map outline using sf instead of rgdal
    uk.line <- sf::st_read(ukcoast_line_shp)

    ## Convert to stream of bytes
    ser_uk.line <- serialize(uk.line, connection=NULL, ascii = FALSE)

    ## Insert into the table using prepared statement via DBI
    ## This handles binary data properly without needing special escape functions
    DBI::dbExecute(
        conn, 
        "INSERT OR IGNORE INTO plotting_tools_map_tools (description, plot_object) VALUES (?, ?)",
        params = list('map_outline', list(ser_uk.line))
    )

    ## Close the connection
    DBI::dbDisconnect(conn)
    
}

## Step 3 Make map objects for DB

### 3.1 Map preparation

#' Prepare map data for visualization
#'
#' @description Prepares UK map data for use in OTU distribution visualization.
#' 
#' @details The function:
#'   1. Reads UK coastline shape files
#'   2. Sets the correct coordinate reference system
#'   3. Creates an interpolation grid based on the UK polygon bounding box
#'   4. Saves the processed spatial objects as shapefiles
#'
#' @param ukcoast_poly Character string. Path to UK polygon shapefile.
#' @param ukcoast_line Character string. Path to UK coastline shapefile.
#' @param uk_poly_converted Character string. Path where the processed UK polygon will be saved.
#' @param uk_line_converted Character string. Path where the processed UK coastline will be saved.
#' @param uk_grid_converted Character string. Path where the interpolation grid will be saved.
#'
#' @return None. The function saves the processed spatial objects to the specified output files.
#'
#' @examples
#' map_prep(
#'   "data/map/ukcoast1.shp",
#'   "data/map/ukcoast_line.shp",
#'   "output/map/uk_poly.shp",
#'   "output/map/uk_line.shp",
#'   "output/map/uk_grid.shp"
#' )
#' 
#' @export
map_prep <- function(ukcoast_poly,
                     ukcoast_line,
                     uk_poly_converted,
                     uk_line_converted,
                     uk_grid_converted) {
    
    # Explicitly define coordinate reference systems
    ukgrid <- sf::st_crs(27700)  # British National Grid
    latlong <- sf::st_crs(4326)  # WGS84 Lat/Long
    
    # Read UK polygon shapefile
    uk_poly <- sf::st_read(ukcoast_poly)
    uk_poly <- sf::st_set_crs(uk_poly, latlong)
    uk_poly <- sf::st_transform(uk_poly, ukgrid)
    
    # Read UK coastline shapefile
    uk_line <- sf::st_read(ukcoast_line)
    uk_line <- sf::st_set_crs(uk_line, latlong)
    uk_line <- sf::st_transform(uk_line, ukgrid)
  
    # Create interpolation grid based on UK polygon bounding box
    uk_poly_bbox <- sf::st_bbox(uk_poly)
    x_range <- c(uk_poly_bbox["xmin"], uk_poly_bbox["xmax"])
    y_range <- c(uk_poly_bbox["ymin"], uk_poly_bbox["ymax"])
  
    # Create grid for interpolation
    grd <- expand.grid(
        x = seq(from = x_range[1], to = x_range[2], by = 5000),
        y = seq(from = y_range[1], to = y_range[2], by = 5000)
    )
  
    # Convert grid to SpatVector
    grd_vect <- terra::vect(grd, geom = c("x", "y"), crs = ukgrid$input)
  
    # Optional: Visualization
    terra::plot(terra::vect(uk_poly))
    terra::plot(terra::vect(uk_line), add = TRUE, col = "red")
    terra::plot(grd_vect, add = TRUE, col = "blue", pch = ".")
    
    # Save processed spatial objects
    sf::st_write(uk_poly, uk_poly_converted, driver = "ESRI Shapefile", overwrite = TRUE)
    sf::st_write(uk_line, uk_line_converted, driver = "ESRI Shapefile", overwrite = TRUE)
    terra::writeVector(grd_vect, uk_grid_converted)
}


### 3.2 Generate maps per OTU
#' Generate map for an individual OTU
#' 
#' @description Creates a map object representing the geographical distribution of an individual OTU/ASV.
#' 
#' @details The function:
#'   1. Extracts abundance data for the specified OTU
#'   2. Performs kriging interpolation to predict values across the UK
#'   3. Clips the interpolation to the UK coastline
#'   4. Saves the map object to the maps database
#'   5. Optionally generates a PNG visualization
#'
#' @param OTU_name Character string. Name of the OTU to create a map for.
#' @param OTU_table_in Character string. Path to the filtered OTU table.
#' @param environment_data Character string. Path to the environmental data file.
#' @param Grid_file Character string. Path to the interpolation grid shapefile.
#' @param UK_poly_file Character string. Path to the UK polygon shapefile.
#' @param UK_line_file Character string. Path to the UK coastline shapefile.
#' @param Make_png Logical. Whether to generate a PNG visualization. Default is FALSE.
#'
#' @return None. The function saves the map object to the maps database and optionally 
#'   generates a PNG visualization.
#'
#' @examples
#' save_otu_map(
#'   OTU_name = "OTU1",
#'   OTU_table_in = "output/Supplementary/Tables_in_SQL/OTU_abund.csv",
#'   environment_data = "output/Supplementary/Tables_in_SQL/Env.csv",
#'   Grid_file = "output/map/uk_grid.shp",
#'   UK_poly_file = "output/map/uk_poly.shp",
#'   UK_line_file = "output/map/uk_line.shp",
#'   Make_png = TRUE
#' )
#'
#' @export
save_otu_map <- function(OTU_name,
                         OTU_table_in,
                         environment_data,
                         Grid_file,
                         UK_poly_file,
                         UK_line_file,
                         maps_db = "maps_db.sqlite",
                         Make_png = FALSE) {

  # Read input data
  env_table <- data.table::fread(environment_data)
  otu_table <- data.table::fread(OTU_table_in)
  env_table <- env_table[ID %in% otu_table$ID]

  # Convert to data frames
  OTU_table <- data.frame(otu_table, row.names = 1, check.names = FALSE)
  env_data <- data.frame(env_table, row.names = 1, check.names = FALSE)

  # Ensure env data is in the same order as OTU table
  Env_sub <- env_data[row.names(OTU_table), ]

  # Extract OTU abundance data
  otu_abund <- OTU_table[, OTU_name, drop = FALSE]
  dat <- cbind(Env_sub[, c('E_2_FIG_10KM', 'N_2_FIG_10KM')], otu_abund)
  colnames(dat)[3] <- "OTU"
  dat <- dat[complete.cases(dat), ]

  # Read spatial data
  Grid <- sf::st_read(Grid_file)
  UK_poly <- sf::st_read(UK_poly_file)
  UK_line <- sf::st_read(UK_line_file)

  # Convert data to sf object
  dat_sf <- sf::st_as_sf(dat, coords = c("E_2_FIG_10KM", "N_2_FIG_10KM"), crs = 27700)
  
  # Perform kriging interpolation
  spc.idw <- gstat::krige(OTU ~ 1, dat_sf, Grid)
  ukgrid <- 27700
  spc.idw_sf <- sf::st_as_sf(spc.idw)
  sf::st_crs(spc.idw_sf) <- ukgrid
  
  # Clip interpolation to UK coastline
  overlay.idw <- sf::st_intersection(spc.idw_sf, UK_poly)
  newmap <- overlay.idw

  # Define break intervals for visualization
  minv <- min(otu_abund)
  maxv <- max(otu_abund)
  at <- c(
    minv, 
    signif(maxv / 256, 1), 
    signif(maxv / 128, 1), 
    signif(maxv / 64, 1), 
    signif(maxv / 32, 1), 
    signif(maxv / 16, 1), 
    signif(maxv / 8, 1), 
    signif(maxv / 4, 2), 
    signif(maxv / 2, 2), 
    maxv
  )

  # Store breakpoints as JSON for db
  sf_json <- jsonlite::toJSON(sf::st_geometry(newmap))
  break_points_json <- jsonlite::toJSON(at, digits =NA) #NA prevents rounding

  # convert to list when serialising as SQLite sees serialised as many objects  
  serialized_sf <- list(serialize(newmap, NULL))
    
  # Connect to maps database
  maps_db_conn <- DBI::dbConnect(RSQLite::SQLite(), maps_db, synchronous = "normal")
  DBI::dbExecute(maps_db_conn, "PRAGMA busy_timeout = 5000;")

  # Insert map data into database
  sql_command_maps <- "insert or replace into maps_table (otu_name, map_object, break_points) values (?, ?, ?)"
  DBI::dbExecute(maps_db_conn, sql_command_maps, params = list(OTU_name, serialized_sf, break_points_json))
  DBI::dbDisconnect(maps_db_conn)

  # Generate PNG visualization if requested
  if (Make_png) {
  }
}

### 3.2 Generate maps per OTU
#' Generate maps for all OTUs in parallel
#' 
#' @description Processes all OTUs in the filtered table to create map objects for each one
#'   using parallel processing.
#'
#' @details The function:
#'   1. Sets up a parallel processing cluster
#'   2. Applies the save_otu_map function to each OTU
#'   3. Handles environment variables for proper functioning in parallel
#'
#' @param OTU_table_in Character string. Path to the filtered OTU table. 
#'   Default is 'data/02_processed_data/filtered_otu_table.csv'.
#' @param environment_data Character string. Path to the environmental data file.
#'   Default is 'data/01_pre-processed_data/cs_location_avc.csv'.
#' @param Grid_file Character string. Path to the interpolation grid shapefile.
#'   Default is 'data/02_processed_data/ukcoast_grid.shp'.
#' @param UK_poly_file Character string. Path to the UK polygon shapefile.
#'   Default is 'data/02_processed_data/ukcoast_poly.shp'.
#' @param UK_line_file Character string. Path to the UK coastline shapefile.
#'   Default is 'data/02_processed_data/ukcoast_line.shp'.
#' @param Make_png Logical. Whether to generate PNG visualizations. Default is FALSE.
#'
#' @return None. The function creates map objects for all OTUs and stores them in the maps database.
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
maps_parallelise <- function(
                             OTU_table_in = 'data/02_processed_data/filtered_otu_table.csv',
                             environment_data = 'data/01_pre-processed_data/cs_location_avc.csv',
                             Grid_file = 'data/02_processed_data/ukcoast_grid.shp',
                             UK_poly_file = 'data/02_processed_data/ukcoast_poly.shp',
                             UK_line_file = 'data/02_processed_data/ukcoast_line.shp',
                             Make_png = FALSE
                             ){
    run_save_otu_map <- function(OTU_name) {
        
        Sys.getenv("CONDA_PREFIX")
        Sys.setenv(CONDA_PREFIX="/data/conda/microscope/")
        proj_db_path <- file.path(Sys.getenv("CONDA_PREFIX"), "share", "proj")
        Sys.setenv(PROJ_LIB=proj_db_path)
        Sys.getenv("CONDA_PREFIX")
        
        microscope::save_otu_map(
                        OTU_name = OTU_name,
                        OTU_table_in = OTU_table_in,
                        environment_data = environment_data,
                        Grid_file = Grid_file,
                        UK_poly_file = UK_poly_file,
                        UK_line_file = UK_line_file,
                        Make_png = FALSE
                    )
    }
    
                                        # Create a cluster with 5 nodes
    cl <- parallel::makeCluster(40)
    
                                        # Export the required function to the workers
    parallel::clusterExport(cl, varlist = c("run_save_otu_map"))
    
    filtered_OTU_file='data/02_processed_data/filtered_otu_table.csv'
    filtered_OTU = data.table::fread(filtered_OTU_file)
    OTU_names = colnames(filtered_OTU)
    OTU_name = OTU_names[-1]
    
    results <- parallel::parSapply(cl, OTU_name, run_save_otu_map)
    parallel::stopCluster(cl)
    
}



## Step 3 Make Blast database
###  3.3 Filter fasta file to contain OTU sequences that meet occupancy filter

## Need to filter sequences prior to making blast DB otherwise we will have hits
## in BLast DB with no supplementary info (i.e we dont want to keep sequences
## corresponding to OTU's/ ASV's that did not meet occupancy filter)
## Doing this using biopython.

#' Make BLAST
#' @description blast
#' @export
make_blast_py <- function(filtered_OTU_file='data/02_processed_data/filtered_otu_table.csv',
                                 input_fasta_file='data/01_pre-processed_data/repseqs.fasta',
                                 filtered_fasta_file='data/02_processed_data/filtered_sequences.fasta'
                                ) {
    # Create a Python code block as a string
py_code <- '
from Bio import SeqIO
import pandas as pd
import os
# Get OTUs we want to keep from OTU table
otu_data = pd.read_csv(r.filtered_OTU_file, sep=None, engine="python")
otu_ids = otu_data.columns.values.tolist()

# Only keep records in filtered otu tab (dictionary keys)
records = (record for record in SeqIO.parse(r.input_fasta_file, "fasta") #r.input_fasta_file
           if record.id in otu_ids)
# Write filtered sequences to fasta file
SeqIO.write(records, r.filtered_fasta_file, "fasta")
    '
    assign("filtered_OTU_file", filtered_OTU_file, envir = .GlobalEnv)
    assign("input_fasta_file", input_fasta_file, envir = .GlobalEnv)
    assign("filtered_fasta_file", filtered_fasta_file, envir = .GlobalEnv)
    # Run the Python code
    reticulate::py_run_string(py_code)
}


###  3.4 Make blast database
# Take filtered fasta and make blast database for back end of shiny app

#' Make BLASTdb
#' @description blastdb
#' @export
make_blast_bash <- function(fasta_file, blast_db_out) {
  # Create the BLAST database directory
  system(paste("mkdir -p", shQuote(dirname(blast_db_out))), ignore.stdout = TRUE, ignore.stderr = TRUE)

  # Run makeblastdb
  Sys.setenv(PATH=paste0(Sys.getenv("PATH"), ":/data/conda/microscope/bin"))
  system(paste("makeblastdb -in", shQuote(fasta_file),
               "-dbtype nucl -out", shQuote(blast_db_out)))
}
