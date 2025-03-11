#' Microscope: Generate taxonomic explorer app inputs
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
#' @import DBI (>= 1.2.3)
#' @import sp (>= 2.1)
#' @import sf
#' @import tmap
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
clean_OTU_table <- function(OTU_file, filtered_OTU_file, OTU_table_occupancy_filter = 30) {
    # Read in OTU table
    OTU_tab <- data.frame(data.table::fread(OTU_file), row.names = 1, check.names = FALSE)
    
    # Remove samples with reads less than 5000
    OTU_tab_sub <- OTU_tab[rowSums(OTU_tab) > 5000, ]
    
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
    colnames(Taxonomy) <- c("hit", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
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
#' 
#' @return None. The function creates SQLite database files with the appropriate tables.
#'
#' @examples
#' format_otu_for_Rsqlite(
#'   "output/Supplementary/Tables_in_SQL/abundance_stats.csv",
#'   "output/Supplementary/Tables_in_SQL/Taxonomy.csv",
#'   "output/Supplementary/Tables_in_SQL/OTU_abund.csv"
#' )
#' 
#' @export
format_otu_for_Rsqlite <- function(filtered_abundance_csv, filtered_taxonomy_csv, filtered_otu_csv) {
    print("Debug: Load files")
    
    # Read the input files
    abundance_csv <- read.csv(filtered_abundance_csv)
    taxonomy_csv <- read.csv(filtered_taxonomy_csv)
    otu_csv <- read.csv(filtered_otu_csv)
    
    # Create abundance table
    sql_command <- sprintf(
        "create table if not exists abund_table (hit character varying(30), %s numeric, primary key (hit))",
        paste(colnames(abundance_csv)[-1], collapse = ' numeric, ')
    )
    
    # Connect to abundance database and create table
    abundance_db <- DBI::dbConnect(RSQLite::SQLite(), "abundance_db.sqlite")
    DBI::dbExecute(conn = abundance_db, statement = sql_command)
    
    # Fill table and disconnect
    DBI::dbWriteTable(abundance_db, "abund_table", abundance_csv, append = TRUE, row.names = FALSE)
    DBI::dbDisconnect(abundance_db)
    
    # Create taxonomy table
    sql_command_taxonomy <- sprintf(
        "create table taxonomy_table (hit character varying (30), %s character varying (250), primary key (hit))",
        paste(
            '"', colnames(taxonomy_csv)[-1], '"',
            collapse = ' charater varying(250),',
            sep = ''
        )
    )
    
    # Connect to taxonomy database and create table
    taxonomy_db <- DBI::dbConnect(RSQLite::SQLite(), "taxonomy_db.sqlite")
    DBI::dbExecute(conn = taxonomy_db, statement = sql_command_taxonomy)
    
    # Fill table and disconnect
    DBI::dbWriteTable(taxonomy_db, "taxonomy_table", taxonomy_csv, append = TRUE, row.names = FALSE)
    DBI::dbDisconnect(taxonomy_db)
    
    # OTU table preparation (transpose due to column limit)
    transpose_otu_csv <- t(otu_csv)
    otu_csv <- data.frame(hit = row.names(transpose_otu_csv), transpose_otu_csv, check.names = FALSE)
    
    # Create OTU table
    sql_command_otu <- sprintf(
        "create table otu_table (hit character varying (30), %s character varying (30), primary key (hit))",
        paste(
            '"', colnames(otu_csv)[-1], '"',
            collapse = ' charater varying(250),',
            sep = ''
        )
    )
    
    # Connect to OTU database and create table
    otu_db <- DBI::dbConnect(RSQLite::SQLite(), "otu_db.sqlite")
    DBI::dbExecute(conn = otu_db, statement = sql_command_otu)
    
    # Fill table and disconnect
    DBI::dbWriteTable(otu_db, "otu_table", otu_csv, append = TRUE, row.names = FALSE)
    DBI::dbDisconnect(otu_db)
    
    # Create maps table
    sql_command_maps <- "create table maps_table (otu_name character varying (30), map_object blob, primary key (otu_name))"
    
    # Connect to maps database and create table
    maps_db <- DBI::dbConnect(RSQLite::SQLite(), "maps_db.sqlite")
    DBI::dbExecute(conn = maps_db, statement = sql_command_maps)
    DBI::dbDisconnect(maps_db)
    
    # Maps table will be filled in step 3 using save_otu_map function
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
  dat <- cbind(env_data[, c('E_2_FIG_10KM', 'N_2_FIG_10KM')], otu_abund)
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

  # Combine map and interval data
  mapandinfo <- cbind(newmap, at)

  # Save map data locally if requested
  if (Make_png) {
    save(mapandinfo, file = paste0('data/02_processed_data/per_otu_map_data/', OTU_name, ".RData"))
  }
  
  # Serialize map data for database storage
  ser_mapandinfo <- serialize(mapandinfo, connection = NULL, ascii = FALSE)

  # Connect to maps database
  maps_db <- DBI::dbConnect(RSQLite::SQLite(), "maps_db.sqlite", synchronous = "normal")
  DBI::dbExecute(maps_db, "PRAGMA busy_timeout = 5000;")

  # Insert map data into database
  sql_command_maps <- "insert into maps_table (otu_name, map_object) values (?, ?)"
  DBI::dbExecute(maps_db, sql_command_maps, params = list(OTU_name, list(ser_mapandinfo)))
  DBI::dbDisconnect(maps_db)

  # Generate PNG visualization if requested
  if (Make_png) {
    tmap::tmap_mode("plot")
    options(tmap.raster.backend = "raster")
    
    # Create map plot using tmap
    map_plot <- tmap::tm_shape(newmap) + 
      tmap::tm_grid(breaks = at, col = "viridis", title = "Predicted Values") +
      tmap::tm_shape(UK_line) + tmap::tm_lines(col = "black", lwd = 2) +
      tmap::tm_layout(title = OTU_name, legend.outside = TRUE)

    # Create PNG directory if needed
    dir.create("output/png", showWarnings = FALSE, recursive = TRUE)

    # Save plot as PNG
    grDevices::png(file = paste0("output/png/", OTU_name, "_plot.png"), width = 200, height = 400)
    print(map_plot)
    grDevices::dev.off()
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
    print(OTU_name)
    
    results <- parallel::parSapply(cl, OTU_name, run_save_otu_map)
    parallel::stopCluster(cl)
    
}



## Step 3 Make Blast database
###  3.3 Filter fasta file to contain OTU sequences that meet occupancy filter

## Need to filter sequences prior to making blast DB otherwise we will have hits
## in BLast DB with no supplementary info (i.e we dont want to keep sequences
## corresponding to OTU's/ ASV's that did not meet occupancy filter)
## Doing this using biopython.


make_blast_py <- function(){
    ## **python**
    ```{python filter fasta file, eval=FALSE}
    from Bio import SeqIO
    import os
    ## get OTUs we want to keep from OTU table
    ## python can access r variables within markdown and converts our r dataframe into a dictionary- the keys are the OTU names
    OTU_tab_dict = r.OTU_tab_sub_occ_dec
    ## make output dir for filtered fasta
    if not os.path.exists(r.Output_dir_with_occ+"/Supplementary/Filtered_sequences"):
               os.makedirs(r.Output_dir_with_occ+"/Supplementary/Filtered_sequences")
    ## only keep records in filtered otu tab (dictionary keys) 
    ## see http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec372 for reference
    records = (record for record in SeqIO.parse(r.params["Fasta_file"],"fasta") if record.id in OTU_tab_dict.keys())
    SeqIO.write(records, r.Output_dir_with_occ+"/Supplementary/Filtered_sequences/filtered_sequences.fasta", "fasta")
    ```
}

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
