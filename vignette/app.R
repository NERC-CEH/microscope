# Load libraries with explicit namespace
library(shiny)
library(shinyjs)
library(shinydashboard)
library(ggplot2)
library(shinycssloaders)
library(pool)
library(DBI)
library(RSQLite)
library(sf)
library(DT)
library(bslib)
library(jsonlite)
library(splitstackshape)

devtools::source_url("https://github.com/NERC-CEH/UKCEH_shiny_theming/blob/main/theme_elements.R?raw=TRUE")

#theming with bslib
UKCEH_theme <- bs_theme(
  bg = "#fff",
  fg = "#292C2F",
  primary = "#0483A4",
  secondary = "#EAEFEC",
  success = "#37a635",
  info = "#34b8c7",
  warning = "#F49633",
  base_font = font_link(family = "Montserrat",href = "https://fonts.googleapis.com/css2?family=Montserrat:wght@400;600&display=swap")
)

#increase the font weight of the headings (h1, h2, h3, h4, h5, h6)
UKCEH_theme <- bs_add_variables(UKCEH_theme,
                                # low level theming
                                "headings-font-weight" = 600)

#titlePanel replacement
UKCEH_titlePanel <- function(title = "UKCEH Shiny app", windowTitle = title){
  
  div(
    img(src="https://www.ceh.ac.uk/sites/default/files/images/theme/ukceh_logo_long_720x170_rgb.png",style="height: 50px;vertical-align:middle;"),
    
    h2(  
      title,
      style ='vertical-align:middle; display:inline;padding-left:150px;'
    ),
    tagList(tags$head(tags$title(paste0(windowTitle," | UK Centre for Ecology & Hydrology")),
                      tags$link(rel="shortcut icon", href="https://brandroom.ceh.ac.uk/themes/custom/ceh/favicon.ico"))),
    style = "padding: 30px;"
  )
}


# Custom dashboard header function (to remove default buttons)

#============================================================
# UI DEFINITION
#============================================================

# Define UI for application that draws a histogram



ui <- fluidPage(
  theme = UKCEH_theme,
  
  UKCEH_titlePanel("Microscope"),
  
  tags$head(
    tags$style(HTML("
      .nav-tabs .nav-link.active {
        background-color: #2f7ece !important;
        color: #FFF !important;
      }
      .well {
        background-color: #c9e7f2;
        color: white;
      }
      .custom-box {
        border: 4px solid #c9e7f2; 
      }
    "))
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 2,
      style = "width: 250px; padding-left:0px; padding-right:0px;",
      
      fluidRow(
        style = 'padding-left:20px',
        actionButton("more_info_button", "More Information", 
                     icon = icon("info-circle"),  
                     style = "color: #fff; background-color: #2f7ece; border-color: #2f7ece", 
                     width = 200)
      ),
      fluidRow(
        style = 'padding-top:15px; padding-left:20px',
        actionButton("github_button", "  GitHub Page  ", 
                     onclick = "window.open('https://github.com/NERC-CEH/microscope', '_blank')", 
                     style = "color:#fff; background-color: #2f7ece; border-color: #2f7ece", 
                     width = 200)
      )
    ),
    
    mainPanel(
      width = 8,
      style = "padding-left:80px; padding-right:0px;",
      
      useShinyjs(),
      
      shinyjs::hidden(
        div(
          id = "more_info",
          tags$h4(tags$b("Summary")),
          shiny::tags$p("Re-implimentation of https://connect-apps.ceh.ac.uk/ID-TaxER/ using a modular workflow to build the data behind the app.
        ID-TaxER provides an interface to explore potential soil habitat preferences of bacterial taxa derived from 16S rRNA gene sequencing. Query sequences are blasted against a database of representative sequences of 97% OTUs obtained from a large soil survey conducted across Britain (the Countryside Survey). Each sequence in the database is linked to an additional trait matrix containing taxonomic assignments as well as environmentally derived information about that OTU (e.g pH or habitat preference). Results are displayed as an interactive table of hits with percentage match to a CS sequence, and associated taxonomy (greengenes). Upon selecting a hit, a plot of model fit to soil parameters is displayed indicating for example the pH optima of that taxon, as well as habitat preferences and spatial distribution (currently Britain only).")
          )
      ),
      
      br(),
      
      textInput(
        inputId = "mysequence", 
        label = "Please enter a sequence",
        value = "",
        width = 10000, 
        placeholder = ''
      ),
      
      br(),
      
      actionButton("blast", "Blast", style = "color: #fff; background-color:#2f7ece"),
      actionButton("clearInput", "Clear Input", style = "color: #fff; background-color:#34b8c7"),
      actionButton("exampleSequence", "Example Sequence", style = "color: #fff; background-color:#90c164"),
      
      selectInput(
        "db_choice", 
        "Choose SQLite Database:",
        choices = list(
          "16s" = "data/16s/",
          "18s" = "data/18s/"
        ),
        selected = "data/16s/"
      ),
      
      actionButton("connect_btn", "Connect to Database", style = "color: #fff; background-color:#2f7ece"),
      
      span(textOutput("Warning"), style = "color:red; font-size:17px"),
      
      shinyjs::hidden(
        div(
          id = "Results",
          hr(style = "border-color: #fff;"),
          
          fluidRow(
            column(
              width = 7,
              style = 'padding:0px;',
              hr(style = "border-color: #fff;"),
              HTML('<center><h4>Top Hits</h4></center>'),
              box(
                class = "custom-box",
                div(
                  DT::dataTableOutput("Main_output_table"),
                  style = "font-size: 75%"
                ),
                width = 600,
                height = 610
              )
            ),
            
            column(
              width = 5,
              style = 'padding-left:70px; padding-right:0px; padding-top:35px;',
              hr(style = "border-color: #fff;"),
              
              tabsetPanel(
                id = "plotTabset",
                tabPanel(
                  title = "GB Map",
                  box(
                    id = "Map_box",
                    withSpinner(plotOutput('OTU_map'), type = 7),
                    width = 350,
                    height = 425
                  )
                ),
                tabPanel(
                  title = "Habitats",
                  box(
                    id = "AVC_plotbox",
                    plotOutput('OTU_avc_boxplot'),
                    width = 350,
                    height = 425
                  )
                )
              ),
              
              fluidRow(
                HTML('<center><b>Additional Information (User Submitted)</b></center>'),
                box(
                  tags$style(HTML("
                    #Indicators {
                      height:64px;
                      overflow-y:scroll
                    }
                  ")),
                  width = 50,
                  height = 80,
                  div(
                    htmlOutput('Indicators'),
                    style = "font-size: 75%"
                  )
                )
              )
            )
          ),
          
          hr(style = "border-color: #fff;"),
          box(
            class = "custom-box",
            div(
              fluidRow(
                column(
                  width = 12,
                  HTML('<h4><center>Blast Output and Taxonomy</center></h4>'),
                  DT::dataTableOutput("BlastResults"),
                  DT::dataTableOutput("OTU_blast_match")
                )
              )
            ),
            style = "font-size: 75%",
            width = 950,
            height = 200
          )
        )
      ),
      
      hr(style = "border-color: #fff;"),
      div(
        style = "text-align: center;", 
        a("UKCEH Privacy Policy", href = "https://www.ceh.ac.uk/privacy-notice", class = "custom_link")
      )
    )
  )
)

#============================================================
# SERVER DEFINITION
#============================================================
server <- function(input, output, session) {
  db_conn_microscope <- reactiveVal(DBI::dbConnect(RSQLite::SQLite(), paste0("data/16s/molecular_db.sqlite")))
  db_conn_maps <- reactiveVal(DBI::dbConnect(RSQLite::SQLite(), paste0("data/16s/maps_db.sqlite"), flags = RSQLite::SQLITE_RO))
  blast_path <- reactiveVal('data/16s/')
  # AACATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTATAAACAACTTTATACTGTGAAACTGCGTACAGCTCATTAAATCAGTCATTATTTATTTGATGGTTTCTTACTTACATGGATACCCGTAGTAATTCTAGAGCTAATACATGCAAAAAATCCCGACTTTTGAAGGGATGTATTTATTAGATAAAAAACCAATGCGTCCTTCGGGGCGGTTTGTGGTGATTCATAATAACTGATCGAATCGCATGGCTTTGCCGGCGATAGTCCACCTAAGTTTTTGACCTATCAGCTAGACGGTAGGGTATTGTCCTACCGTGGCATTGACGGGTGACGGGGAATTAGGGTTTGATTCCGGAGAGGGAGCCTGAGAGACAGCTACCATGTCTACGGACAGCAACAGGCCCGCAAATTGTCCAATTCCAACATATCGGAGAGACAGTGAAAATTA
  
  observeEvent(input$connect_btn, {
    
    # Clear the sequence information
    shinyjs::reset("mysequence")
    shinyjs::hide("Results")       # Hide results panel
    
    
    if (!is.null(db_conn_microscope())) {
      dbDisconnect(db_conn_microscope())
      db_conn_microscope(NULL)  # Reset connection
    }
    if (!is.null(db_conn_maps())) {
      dbDisconnect(db_conn_maps())
      db_conn_maps(NULL)  # Reset connection
    }
    
    # Get the selected database path
    selected_path <- paste0(input$db_choice)
    
    # Update the path_to_db reactive value
    blast_path(selected_path)
    
    conn_microscope <- DBI::dbConnect(RSQLite::SQLite(), paste0(input$db_choice,"molecular_db.sqlite"))
    
    maps_db <- DBI::dbConnect(RSQLite::SQLite(), paste0(input$db_choice,"maps_db.sqlite"), flags = RSQLite::SQLITE_RO)
    
    db_conn_microscope(conn_microscope)
    db_conn_maps(maps_db)
  })
  
  
  # Functions to call env and uk_line dynamically depending on data source.  
  env <- function() {
    query <- "SELECT * FROM env_table;"
    return(DBI::dbGetQuery(db_conn_microscope(), query))
  }
  
  uk_line <- function() {
    query <- "SELECT plot_object FROM plotting_tools_map_tools where description= 'map_outline';"
    maps_data <- DBI::dbGetQuery(db_conn_maps(), query)
    return(unserialize(maps_data$plot_object[[1]]))
  }
  
  example_sequence <- function() {
    switch(blast_path(),
           "data/16s/" = 'ACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGTGATGCAAGTCTGGTGTGAAATCTCGGGGCTCAACTCCGAAATTGCACCGGATACTGCGTGACTCGAGGACTGTAGAGGAGATCGGAATTCACGGTGTAGCAGTGAAATGCGTAGATATCGTGAGGAAGACCAGTTGCGAAGGCGGATCTCTGGGCAGTTCCTGACACTGAGGCACGAAGGCCAGGGGAGCAAACGGG',
           "data/18s/" = 'AACATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTATAAACAACTTTATACTGTGAAACTGCGTACAGCTCATTAAATCAGTCATTATTTATTTGATGGTTTCTTACTTACATGGATACCCGTAGTAATTCTAGAGCTAATACATGCAAAAAATCCCGACTTTTGAAGGGATGTATTTATTAGATAAAAAACCAATGCGTCCTTCGGGGCGGTTTGTGGTGATTCATAATAACTGATCGAATCGCATGGCTTTGCCGGCGATAGTCCACCTAAGTTTTTGACCTATCAGCTAGACGGTAGGGTATTGTCCTACCGTGGCATTGACGGGTGACGGGGAATTAGGGTTTGATTCCGGAGAGGGAGCCTGAGAGACAGCTACCATGTCTACGGACAGCAACAGGCCCGCAAATTGTCCAATTCCAACATATCGGAGAGACAGTGAAAATTA',
           'ACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGTGATGCAAGTCTGGTGTGAAATCTCGGGGCTCAACTCCGAAATTGCACCGGATACTGCGTGACTCGAGGACTGTAGAGGAGATCGGAATTCACGGTGTAGCAGTGAAATGCGTAGATATCGTGAGGAAGACCAGTTGCGAAGGCGGATCTCTGGGCAGTTCCTGACACTGAGGCACGAAGGCCAGGGGAGCAAACGGG')
  }
  
  
  # Make Results visible when blast button is clicked
  shinyjs::onclick("blast", shinyjs::show(id = "Results", anim = TRUE))
  
  # Toggle more info section
  shinyjs::onclick("more_info_button", shinyjs::toggle(id = "more_info", anim = TRUE))
  
  # Blast function definition
  make.comparison <- function(query) {
    # Check if query isn't empty
    if (query != "") { 
      # Create blast command
      
      cmd <- paste('echo -e ">Query seq\n', query, '"', 
                   '|/data/conda/microscope/bin/blastn -db ',blast_path(),'Blast_DB/Blast_DB -num_alignments 20 -evalue 0.001 -outfmt 7',sep='')
      
      # Run system command and capture output
      blast_capture <- system(paste("/bin/bash -c", shQuote(cmd)), intern = TRUE)
      
      # Check if there are hits
      if (blast_capture[4] != "# 0 hits found") {
        # Set output switch to on
        output_switch <- "on"
        
        # Process blast output
        blast_capture <- blast_capture[6:(length(blast_capture) - 1)]
        blast_capture_df <- as.data.frame(splitstackshape::cSplit(
          as.data.frame(blast_capture), "blast_capture", sep = "\t"))
        
        # Remove first column with query ID
        blast_capture_df <- blast_capture_df[, -1]
        
        
        # SQL query to get taxonomy and abundance data
        SQL_command <- paste("
          SELECT taxonomy_table.*, abund_table.abundance_rank, abund_table.occupancy_proportion
          FROM taxonomy_table
          JOIN abund_table ON taxonomy_table.hit = abund_table.hit
          WHERE taxonomy_table.hit IN ('", paste(blast_capture_df$blast_capture_02, collapse = "', '"), "');", sep = "")
        
        relevant_tax_and_stats <- DBI::dbGetQuery(db_conn_microscope(), SQL_command)
        
        # Merge blast results with taxonomy data
        blast_relevant_tax_and_stats <- merge(blast_capture_df, relevant_tax_and_stats, by = 1)
        
        # Clean up and format the data
        blast_relevant_tax_and_stats <- blast_relevant_tax_and_stats[, c(2, 1, 12:20)]
        colnames(blast_relevant_tax_and_stats)[1:2] <- c("identity", "hit")
        
        # Order like original blast output
        row.names(blast_relevant_tax_and_stats) <- blast_relevant_tax_and_stats[, 2]
        blast_relevant_tax_and_stats <- blast_relevant_tax_and_stats[blast_capture_df$blast_capture_02, ]
        
        # Add column names to full blast output
        colnames(blast_capture_df) <- c("subject id", "% identity", "alignment length", 
                                        "mismatches", "gap opens", "q.start", "q.end", 
                                        "s.start", "s.end", "evalue", "bit score")
        
        # Get abundance data for OTUs
        SQL_command <- paste0("
          SELECT *
          FROM otu_table
          WHERE otu_table.hit IN('", paste(blast_capture_df$`subject id`, collapse = "', '"), "');")
        
        #OTU_abund <- DBI::dbGetQuery(db_conn_otu(), SQL_command)#
        OTU_abund <- DBI::dbGetQuery(db_conn_microscope(), SQL_command)
        
        # Remove duplicate columns and format data
        OTU_abund <- OTU_abund[, unique(colnames(OTU_abund))]
        row.names(OTU_abund) <- OTU_abund[, 1]
        OTU_abund <- OTU_abund[, -1]
        
        # Order same as blast output
        OTU_abund <- OTU_abund[blast_capture_df$`subject id`, ]
        
        # Transpose
        OTU_abund <- data.frame(t(OTU_abund))
        
        # Return all data and output switch
        return(list(
          tax_and_stats = blast_relevant_tax_and_stats,
          full_blast_output = blast_capture_df,
          OTU_abundance = OTU_abund,
          output_switch = "on"
        ))
      } else {
        # No hits found
        return(list("output_switch" = "off"))
      }
    } else {
      # Empty query
      return(list("output_switch" = "off"))
    }
  }
  
  # Run blast when button is clicked
  run_sequence <- shiny::eventReactive(input$blast, {
    make.comparison(input$mysequence)
  })
  
  # Main output table
  output$Main_output_table <- DT::renderDataTable({
    run_sequence_output <- run_sequence()
    if (run_sequence_output$output_switch == "on") {
      run_sequence_output$tax_and_stats
    }
  }, 
  selection = list(mode = 'single', selected = 1),
  options = list(scrollX = TRUE, pageLength = 7, dom = 'tp'),
  rownames = FALSE)
  
  # OTU Map plot
  output$OTU_map <- shiny::renderPlot({
    run_sequence_output <- run_sequence()
    if (run_sequence_output$output_switch == "on") {
      # Get selected row
      s <- input$Main_output_table_rows_selected
      
      # If row has been selected
      if (length(s)) {
        # Get OTU from run_sequence output
        OTU <- run_sequence_output$tax_and_stats[s, 2]
        par(mar = c(4, 4, 1, 4))
        
        # Get relevant map
        otu_data <- DBI::dbGetQuery(db_conn_maps(), 
                                    paste("SELECT map_object, break_points FROM maps_table WHERE otu_name='", toString(OTU), "';", sep = ""))
        
        # Convert from bytea to R object
        #r_object_map <- unserialize(unescape_bytea_map$map_object[[1]])
        serialized_points <- otu_data$map_object[[1]]
        db_break_points <- otu_data$break_points[1]
        r_break_points <- jsonlite::fromJSON(db_break_points)
        r_sf_map <- unserialize(serialized_points)
        
        # Temporary fix to solve issue in break points not serialising correctly
        minv <- r_break_points[1]#min(otu_abund)                                                                                                                                                                                           
        maxv <- r_break_points[10]#max(otu_abund)                                                                                                                                                                                           
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
        r_break_points <- at
        # temp fix ^^^ remove
        
        # Plot the map
        plot(r_sf_map["var1.pred"], main = "OTU Abundance Map", breaks = unique(r_break_points), pch=20, cex = 0.8)
        #        plot(r_object_map["var1.pred"], main = "OTU Abundance Map", breaks = unlist(at), pch=20)
      }
    }
  })
  
  # OTU Habitat boxplot
  output$OTU_avc_boxplot <- shiny::renderPlot({
    run_sequence_output <- run_sequence()
    if (run_sequence_output$output_switch == "on") {
      # Get selected row
      s <- input$Main_output_table_rows_selected
      
      # If row has been selected
      if (length(s)) {
        # Get OTU from run_sequence output
        abund <- run_sequence_output$OTU_abundance[, s, drop = FALSE]
        # Merge with environmental data
        abund_env <- merge(abund, env(), by.x = 0, by.y = 1)
        
        # Order AVC factor levels to represent habitat gradient
        ord <- stats::reorder(abund_env$avc, X = as.numeric(abund_env$avc_code), FUN = mean)
        
        # Convert abundance column to numeric (this fixes the error)
        abund_env[, 2] <- as.numeric(abund_env[, 2])
        
        # Set plot margins
        par(mar = c(11, 5.8, 2, 0.4) + 0.1)
        bxp_stats <- graphics::boxplot(abund_env[, 2] ~ ord, plot = FALSE)
        max_whisker <- max(bxp_stats$stats[5,])  # Row 5 contains upper whiskers
        
        # Create boxplot
        graphics::boxplot(abund_env[, 2] ~ ord, 
                          las = 2, 
                          ylab = "", 
                          xlab = "", 
                          ylim = c(0, 1.1 * max_whisker), 
                          cex = 0.5, 
                          outline = FALSE)
        
        # Add title
        title(ylab = paste("Relative Abundance (", colnames(abund_env)[2], ")", sep = ""), 
              line = 4.5, 
              cex.lab = 0.9)
      }
    } else {
      plot.new()
    }
  })
  
  # Blast match table
  output$OTU_blast_match <- DT::renderDataTable({
    run_sequence_output <- run_sequence()
    
    if (run_sequence_output$output_switch == "on") {
      # Get selected row
      s <- input$Main_output_table_rows_selected
      
      # If row has been selected
      if (length(s)) {
        run_sequence_output$full_blast_output[s, ]
      }
    }
  }, 
  selection = 'none',
  options = list(dom = 't'),
  rownames = FALSE,
  colnames = c("subject id", "% identity", "alignment length", "mismatches", 
               "gap opens", "q.start", "q.end", "s.start", "s.end", 
               "evalue", "bit score"))
  
  # Clear input button event
  shiny::observeEvent(input$clearInput, {
    shinyjs::reset("mysequence")
  })
  
  # Example sequence button event
  shiny::observeEvent(input$exampleSequence, {
    shiny::updateTextInput(session, "mysequence", 
                           value = example_sequence()#"ACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGTGATGCAAGTCTGGTGTGAAATCTCGGGGCTCAACTCCGAAATTGCACCGGATACTGCGTGACTCGAGGACTGTAGAGGAGATCGGAATTCACGGTGTAGCAGTGAAATGCGTAGATATCGTGAGGAAGACCAGTTGCGAAGGCGGATCTCTGGGCAGTTCCTGACACTGAGGCACGAAGGCCAGGGGAGCAAACGGG")
    )
  }
  )
  
  # Warning message
  output$Warning <- shiny::renderText({ 
    run_sequence_output <- run_sequence()
    if (run_sequence_output$output_switch == "off") {
      "No Hits Found!"
    }
  })
  

  # Close connections when app terminates
  shiny::onSessionEnded(function() {
    if (!is.null(isolate(db_conn_microscope()))) {
      DBI::dbDisconnect(isolate(db_conn_microscope()))
    }
    if (!is.null(isolate(db_conn_maps()))) {
      DBI::dbDisconnect(isolate(db_conn_maps()))
    }
  })
}

# Run the application 
shiny::shinyApp(ui = ui, server = server)


