

#============================================================================================================================================
#load libraries

library(DT)
library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinycssloaders)
library(bslib)


#rather than reading this in have just copied and pasted so I could edit the left padding for the title

#devtools::source_url("https://github.com/NERC-CEH/UKCEH_shiny_theming/blob/main/theme_elements.R?raw=TRUE")

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




# Define UI for application that draws a histogram
ui <- fluidPage(
  
  #theme
  theme = UKCEH_theme,
  
  # Application title
  #div(style = "text-align: center;",
  UKCEH_titlePanel("ID-TaxER - Identification of Taxa & Environment Responses"),
  #),
  #change tab selection color
  #change side panel colour
  #define box class with blue outline  
  tags$head(
    tags$style(HTML("
      /* tabsetPanel tab selection */
      .nav-tabs .nav-link.active {
       background-color: #2f7ece !important;
        color: #FFF !important;
       }
       /* side panel colour */
      .well {
        background-color: #c9e7f2;
        color: white;
      }
       /* coloured box outline */
       .custom-box {
      border: 4px solid #c9e7f2; 
     }
     
   
    "))
  ),
  
  sidebarLayout(
    sidebarPanel(width = 2, 
                 # Allow for custom CSS
                 style = "width: 250px; padding-left:opx;padding-right:0px;",
                 
                 #here I want the buttons to not be too close vertically so have used fluid rows
                 fluidRow(style='padding-left:20px',actionButton("more_info_button", "More Information",icon=icon("info-circle"),  style="color: #fff; background-color: #2f7ece; border-color: #2f7ece",width=200)),
                 fluidRow(style='padding-top:15px; padding-left:20px',actionButton("github_button", "  GitHub Page  ",onclick ="window.open('https://github.com/brijon/ID-TaxER-flat-files', '_blank')", style="color:#fff; background-color: #2f7ece; border-color: #2f7ece",width=200)),
            
                 
                 
                 # actionButton("more_info_button", "More Information",   style="color: #fff; background-color: #2f7ece;
                 # border-color: #2f7ece",width=200),
                 # hr(),
                 ##github button
                 #actionButton("github_button", "  GitHub Page  ",onclick ="window.open('https://github.com/brijon/ID-TaxER-flat-files', '_blank')",
                 #            style="color:#fff; background-color: #2f7ece; border-color: #2f7ece",width=200)
                 #end of sidebarPanel
                 #width=2),
    ),
    
    mainPanel(width = 8, 
              style="padding-left:80px; padding-right:0px;",
              #Hidden section about app
              useShinyjs(),
              #tags$h is heading, tags$p is paragraph etc , tags$b is bold etc                
              shinyjs::hidden(div(id="more_info",tags$h4(tags$b("Summary")),
                                  tags$p( 
                                    "ID-TaxER provides an interface to explore potential soil habitat preferences of bacterial taxa derived from 16S rRNA gene sequencing. Query sequences are blasted against a database of representative sequences of 97% OTUs obtained from a large soil survey conducted across Britain (the Countryside Survey). Each sequence in the database is linked to an additional trait matrix containing taxonomic assignments as well as environmentally derived information about that OTU (e.g pH or habitat preference). Results are displayed as an interactive table of hits with percentage match to  a CS sequence, and associated taxonomy (greengenes). Upon selecting a hit, a plot of model fit to soil parameters is displayed indicating for example the pH optima of that taxon, as well as habitat preferences and spatial distribution (currently Britain only).  
		"),
                                  tags$h4(tags$b("Limitations")),
                                  tags$p("The database encompasses the V3-V4 region of the 16SrRNA, amplified with 341f/806r primers. Queries which do not cover this region will obviously give incorrect results, and additionally taxa poorly amplified with these primers will be under represented. Importantly this tool is based on homology mapping to a short portion of the conserved 16S rRNA gene, and so all the usual limitations apply regarding accuracy of taxonomic (and habitat preference) assignment.  It is therefore for research purposes only."),
                                  tags$h4(tags$b("Ongoing work")),
                                  
                                  tags$p("A",tags$a(href="https://github.com/brijon/ID_TaxER-Custom-Database-for-DADA2","github",target="_blank")," page has been set up for batch sequence querying using Dada2. We are exploring options to allow users to upload their own ecological trait information to the trait matrix (e.g if a sequence with high homology to a CS sequence comes up in user experiments as an indicator of warming, drought, plant species X etc then it would be useful to capture this information in the trait matrix).
		We also have similar ITS and 18S datasets which could be developed in a similar portal if enough interest."),
                                  tags$h4(tags$b("Contact:")),
                                  tags$p("
		Briony Jones", tags$b("(brijon@ceh.ac.uk)")),
                                  
                                  tags$p("Rob Griffiths",tags$b(" (rig@ceh.ac.uk)
                ")) 
                                  #end of more info divider                    
              )
              #end of more info hidden javascript section
              ),
              br(),
              #enter sequence box      
              textInput(inputId = "mysequence",label="Please enter a sequence",
                        value="",
                        width = 10000, placeholder = ''),
              br(),
              #various buttons               
              actionButton("blast", "Blast",style="color: #fff; background-color:#2f7ece"),
              actionButton("resetSequence", "Clear Input",style="color: #fff; background-color:#34b8c7"),
              actionButton("exampleSequence", "Example Sequence",style="color: #fff; background-color:#90c164"),
              #area where any warning messages appear                
              span(textOutput("Warning"),style="color:red;font-size:17px  "),
              
              #hidden results section            
              shinyjs::hidden(div(id="Results", hr(style = "border-color: #fff;"),
                                  #===================================================================================================================================================================================================================   
                                  #hits table
                                  fluidRow(column(width=7,style='padding:0px;',hr(style = "border-color: #fff;"),HTML('<center><h4>Top Hits</h4></center>'),box(class = "custom-box",div(DT::dataTableOutput("blastout"),style = "font-size: 75%"),width=600,height=610)),
                                           #===================================================================================================================================================================================================================                  
                                           #define collumn ie portion of interface to the right for plots etc                
                                           column(width=5,style='padding-left:70px; padding-right:0px;padding-top:35px;',hr(style = "border-color: #fff;"),
                                                  tabsetPanel(id="plotTabset",
                                                              #===================================================================================================================================================================================================================   
                                                              #first plot HOF                
                                                              tabPanel(title="pH Model",box(id="plotbox",plotOutput('modelplot'), width=350,height=425
                                                                                            #end of box 
                                                              )
                                                              #end of plot tab panel
                                                              ),
                                                              #===================================================================================================================================================================================================================   
                                                              #LOESS plot              
                                                              tabPanel(title="pH LOESS",box(id="plotbox2",plotOutput('loessplot'), width=350,height=425
                                                                                            #end of box 
                                                              )
                                                              #end of plot tab panel
                                                              ),
                                                              #===================================================================================================================================================================================================================   
                                                              #Map                
                                                              tabPanel(title="GB Map",box(id="Map_box",withSpinner(plotOutput('map'),type=7),width=350,height=425
                                                                                          #end of box 
                                                              )
                                                              #		end of plot tab panel
                                                              ),                  
                                                              #===================================================================================================================================================================================================================   
                                                              #AVC plot
                                                              tabPanel(title="Habitats",
                                                                       box(id="AVC_plotbox",plotOutput('AVC_box_plot'),width=350,height=425
                                                                           # end of box
                                                                       )
                                                                       # end of plot tab panel
                                                              )
                                                              #end of plot tab panel set
                                                  )
                                                  #===================================================================================================================================================================================================================   
                                                  #indicators text box
                                                  ,fluidRow( HTML('<center><b>Additional Information (User Submitted)</b></center>'),
                                                             box(
                                                               #want to make scrollable box               
                                                               tags$style(HTML("
                                                                     #Indicators {
                                                                     height:64px;
                                                                     overflow-y:scroll
                                                                    }
                                                                   ")),
                                                               width=50,
                                                               height=80,  
                                                               div(htmlOutput('Indicators'),style = "font-size: 75%")
                                                               #end of box               
                                                             )
                                                             #end of fluid row             
                                                  )
                                                  #===================================================================================================================================================================================================================   
                                                  
                                                  
                                                  #end of column
                                           )  
                                           #===================================================================================================================================================================================================================                               
                                           #blast output and taxonomy
                                           ,hr(style = "border-color: #fff;"),box(class = "custom-box",div(fluidRow(column(width=12,HTML('<h4><center>Blast Output and Taxonomy</center></h4>'),DT::dataTableOutput("BlastResults"),DT::dataTableOutput("OTU.Taxon")
                                                                                                                           #end of box
                                           ))
                                           #end of column
                                           ,style = "font-size: 75%")
                                           #end of fluid row
                                           ,width=950,height=200 )
                                           #=====================================================================================================================================================================================================================               
                                           
                                  )
                                  
              )
              #end of hidden java section      
              ),#link to privacy privacy
              hr(style = "border-color: #fff;"),
              img(src="logos_narrow.svg"),
              hr(style = "border-color: #fff;"),
              div(style = "text-align: center;", a("UKCEH Privacy Policy", href = "https://www.ceh.ac.uk/privacy-notice",class="custom_link")),
              
              #end of main panel 
              
    )
    
    
    #end of sidebarlayout
  )
  #end of fluid page   
)

