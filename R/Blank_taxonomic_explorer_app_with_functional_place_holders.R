#load libraries
library(shiny)
library(shinyjs)
#This library enables us to make a dashboardPage and use box() function
#could change the type of shiny app e.g to fluidPage etc, but box functions will not work (so need to remember to take those out)
library(shinydashboard)
library(ggplot2)
library(shinycssloaders)
#for Csplit
library(splitstackshape)
#for database connection
library(RPostgreSQL)
library(pool)
library(sp)

#hack to get rid of  dashboard button
mydashboardHeader <- function(..., title = NULL, disable = FALSE,title.navbar=NULL, .list = NULL) {
  items <- c(list(...), .list)
  tags$header(class = "main-header",
  style = if (disable) "display: none;",
  span(class = "logo", title),
  tags$nav(class = "navbar navbar-static-top", role = "navigation",
  span(shiny::icon("bars"), style = "display:none;"), 
  title.navbar,
  div(class = "navbar-custom-menu",
  tags$ul(class = "nav navbar-nav",items)
  )))
}

#UI PAGE===================================================================================================================================================================================== 

# Define UI for application that draws a histogram
ui <- dashboardPage(
#Define header title with a bit of styling   
  mydashboardHeader(title=tags$p("",style="font-size: 20px; font-family:Calibri; background-color: #337ab7;")),
  #mydashboardHeader(title=tags$p("",style="font-size: 20px; font-family:Calibri; background-color: #AAFF00;")),
  dashboardSidebar(width=230,
###NEED A LOGO PLACE HOLDER HERE
#This is where logo image will be displayed
  img(height=231.141,width=230),      
#more info button
  actionButton("more_info_button", "More Information",   style="color:#097969; background-color: #fff;border-color: #2e6da4",width=200),
#button to take you to relevant external link
  actionButton("External_link_button", "<--specified_External_link_name-->",onclick ="window.open('<--specified_External_link-->', '_blank')",
  style="color:#097969; background-color: #fff; border-color: #2e6da4",width=200)
  ), 
  dashboardBody(
#Hack to add title in the center of the dashboard header (ln 41-56) see https://stackoverflow.com/questions/45176030/add-text-on-right-of-shinydashboard-header
#other code in this section can be used to customise colours background and rest of dashboard     https://stackoverflow.com/questions/31711307/how-to-change-color-in-shiny-dashboard
  tags$head(tags$style(HTML(
###COULD DO COLOUR CUSTOMISATION-MAYBE HAVE GROUPS OF COLOURS TO INSERT
  '.myClass {
  font-size: 30px;
  line-height: 50px;
  text-align: center;
  font-family: "Calibri";
  padding: 0 130px;
  overflow: hidden;
  color: white;
  }
 .content-wrapper{
  background-color: #F7FCFA;
 }
 
  /* logo */
  .skin-blue .main-header .logo {
  background-color: #097969;}
  /* logo when hovered*/
  .skin-blue .main-header .logo:hover {
  background-color: #097969;}
  /* navbar (rest of the header) */
  .skin-blue .main-header .navbar {
  background-color: #097969;}   
  /* main sidebar */
  .skin-blue .main-sidebar {background-color: #03483E;}
  
  '))),
  #add heading
  tags$script(HTML('
  $(document).ready(function() {
  $("header").find("nav").append(\'<span class="myClass"> <--specified_App_title--></span>\');
  })
  ')),
  #Colour of selection for data tables
  tags$style(HTML('table.dataTable tr.selected td, table.dataTable td.selected {background-color: #8BDACE !important;}')),
#MORE INFO
#  Hidden more information section using shinyjs
#  useShinyjs 'function must be called from a Shiny app's UI in order for all other shinyjs' https://www.rdocumentation.org/packages/shinyjs/versions/2.1.0/topics/useShinyjs
useShinyjs(),
#have used some shiny html tags here
#tags$h is heading, tags$p is paragraph etc , tags$b is bold etc more on this here https://shiny.rstudio.com/articles/tag-glossary.html            
  shinyjs::hidden(div(id="more_info",tags$h4(tags$b("Summary"),style="color:#000080;font-family:Calibri;"),tags$p("<--specified_Info_text-->",style="color:#000080;font-family:Calibri;"))),
  br(),
##INPUTS    
#Sequence Input    
  textInput(inputId = "mysequence",label="Please enter a sequence",
  value="",width = 10000, placeholder = ''),
#add empty line
  br(),
#buttons
#blast sequence
  actionButton("blast", "Blast",style="color: #fff; background-color:#DF4242"),
#clear sequence input
  actionButton("clearInput", "Clear Input",style="color: #fff; background-color:#097969"),
#enter example sequence
  actionButton("exampleSequence", "Example Sequence",style="color: #fff; background-color:#74A5E1"),
  br(),
#area where any warning messages appear                
  span(textOutput("Warning"),style="color:red;font-size:17px"),
  br(),
### OUTPUTS
# 'This function must be called from a Shiny app's UI in order for all other shinyjs' https://www.rdocumentation.org/packages/shinyjs/versions/2.1.0/topics/useShinyjs
# useShinyjs(),
#hidden results section
  shinyjs::hidden(div(id="Results", hr(),
#ok lets add some structure for the outputs (we dont just want plots and tables stacked), here im using fluid rows and columns (although many other options for laying out app)
#create fluid row
#first column for main output table
  fluidRow(column(width=8,hr(),HTML('<center><h4>Top Hits</h4></center>'),
#box basically creates a white box around the output
#only works in shiny dashboards not fluid page etc
#some dataTable aesthetic options (e.g number of rows to display/pagination of large tables, colnames displayed etc) are controlled from renderDataTable in server.R                 
  box(DT::dataTableOutput('Main_output_table'),width=500,height=425)
  #end of collumn
  ),
#second row for outputs
  column(width=4,
  hr(),
#Tabset for plots for neatness
#'Tabsets are useful for dividing output into multiple independently viewable sections.' https://www.rdocumentation.org/packages/shiny/versions/1.7.3/topics/tabsetPanel
  tabsetPanel(id="plotTabset",
#Tab plot 1  
#Include loading spinner for map
  tabPanel(title="GB Maps",withSpinner(plotOutput('OTU_map')),width=350,height=425),
#Tab plot 2
  tabPanel(title="Habitats",plotOutput('OTU_avc_boxplot'),width=350,height=425)
#end of tabset
  )
#end of column
  ),
#end of fluid row
  ), 
  fluidRow(column(width=12,HTML('<h4><center>Blast Output </center></h4>'),
  box(DT::dataTableOutput("OTU_blast_match"),width=1000
#end of box
  )
#end of column
    )
    #end of fluid row
    )
    #end of div function
    )
    #emd of shinyjs::hidden
    )
    #end of dashboard body
    )
    #end of dashboard
    )

               
    
# Define server logic required to draw a histogram
###POSTGRES###
server <- function(input, output,session) {
#postgres database connection
  con <- dbPool(
    drv = RPostgreSQL::PostgreSQL(max.con=140),
    dbname = '<--specified_SQL_database_name-->',
    host = '<--specified_SQL_database_host-->',
    port='5432',
    user=Sys.getenv('SQL_USER'),
    password = Sys.getenv('SQL_PWD')
  )  
#ShinyJS command to make output visible
    shinyjs::onclick("blast",shinyjs::show(id="Results",anim=TRUE))  
#ShinyJS command to get more info        	
    onclick("more_info_button",toggle(id="more_info",anim=TRUE))
#get all sample environment information going to use this for habitat box plots
    SQL_command=paste("select * from env_attributes.env_attributes_all;")
#env consists of sample avc_code (habitat code), avc (habitat description) and pH
    env <- dbGetQuery(con, SQL_command)    
#get map outline to use later
    #for purpose of loading serialised objects 
    dbGetQuery(con, "set standard_conforming_strings to 'on'")
    #get uk map outline to plot uk mapping objects  
    SQL_command=paste("select plot_object from plotting_tools.map_tools where description= 'map_outline';")
    uk.line <-unserialize(postgresqlUnescapeBytea( dbGetQuery(con, SQL_command))) 
    
# Blast function    
    make.comparison <- function(query ){
      #check query isnt empty
      if (query!=""){ 
        #blast_command for aligning sequences returns top 20 hits
        #using echo as a way to pass sequence not in a file to blastn #https://www.biostars.org/p/17265/      
        cmd<- paste('echo -e ">Query seq\n',query,'"', '|blastn -db Blast_DB/Blast_DB -num_alignments 20 -evalue 0.001 -outfmt  7')
        #run system command and capture output
        #shQuote Quotes a string to be passed to an operating system shell
        blast_capture<- system(paste("/bin/bash -c", shQuote(cmd)),intern=TRUE)
        #check there are hits
        #first few lines are not results if there are 0 hits this will be displayed on fourth line
        if(blast_capture[4]!="# 0 hits found"){
          #this variable will be used later to identify that output should be displayed (as hits have been returned from blast command)      
          output_switch <-"on"
          #remove all descriptive lines of output that are not results (first five lines and last line of output)
          blast_capture<-blast_capture[6:(length(blast_capture)-1)]
          #split output by tabs
          blast_capture_df=as.data.frame(cSplit(as.data.frame(blast_capture),"blast_capture",sep="\t"))
#get rid of first column with query ID  (ln 147)        
          blast_capture_df=blast_capture_df[,-1]
          #now lets start getting some wider information about these taxa stored in the SQL database         
          #using a join(similar to r merge) to get information in taxonomy and abundance_stats tables using WHERE statement to get relevant ASV/OTUs from blast output
          #tx and abs are SQL aliases for taxonomic and abundance_stats tables to reduce length of command
          SQL_command=paste("SELECT tx.*, abs.abundance_rank, abs.occupancy_proportion FROM <--specified_Schema_table_prefix-->_otu_attributes.<--specified_Schema_table_prefix-->_taxonomy tx JOIN <--specified_Schema_table_prefix-->_otu_attributes.<--specified_Schema_table_prefix-->_abundance_stats abs ON tx.hit = abs.hit WHERE tx.hit IN ('",paste(blast_capture_df$blast_capture_02,collapse="', '"),"');",sep="")
          relevant_tax_and_stats=dbGetQuery(con, SQL_command)
          blast_relevant_tax_and_stats=merge(blast_capture_df,relevant_tax_and_stats,by=1)
          #remove columns we dont want for this table
          blast_relevant_tax_and_stats=blast_relevant_tax_and_stats[,c(2,1,12:20)]
          #change first two colnames 
          colnames(blast_relevant_tax_and_stats)[1:2]=c("identity","hit")
          #now order like original blast output
          row.names(blast_relevant_tax_and_stats)=blast_relevant_tax_and_stats[,2]
          blast_relevant_tax_and_stats=blast_relevant_tax_and_stats[blast_capture_df$blast_capture_02,]
          #add colnames to full blast output as we want to return this too
          colnames(blast_capture_df)= c("subject id","% identity","alignment length","mismatches","gap opens","q.start","q.end","s.start","s.end","evalue","bit score")
          #get abundance data for all otus here , could do this at point when rows are selected (e.g in output$AVC_box_plot code) but as a slow command this would slow things down as navigating the app which would be frustrating
          #could consider getting map data within this function too
          #need to excute join command as data in two parts
          SQL_command=paste0("SELECT * FROM abund_tables.<--specified_Schema_table_prefix-->_abund_1 JOIN abund_tables.<--specified_Schema_table_prefix-->_abund_2 ON abund_tables.<--specified_Schema_table_prefix-->_abund_1.hit=abund_tables.<--specified_Schema_table_prefix-->_abund_2.hit WHERE abund_tables.<--specified_Schema_table_prefix-->_abund_1.hit IN('",paste(blast_capture_df$`subject id`,collapse="', '"),"');")
          OTU_abund=dbGetQuery(con,SQL_command)
          #two hit columns need to remove one
          OTU_abund=OTU_abund[,unique(colnames(OTU_abund))]
          #make remaining hit collumn rownames
          row.names(OTU_abund)=OTU_abund[,1]
          #remove first column
          OTU_abund=OTU_abund[,-1]
          #order in the same order as blast output for consistency and in order to access correct otu data in output$AVC_box_plot code
          OTU_abund= OTU_abund[blast_capture_df$`subject id`,] 
          #transpose 
          OTU_abund=data.frame(t(OTU_abund))
          #return blast_relevant_tax_and_stats, full blast output and output_switch , signalling output should be displayed
          return(list(tax_and_stats=blast_relevant_tax_and_stats,full_blast_output=blast_capture_df,OTU_abundance=OTU_abund,output_switch="on"))
        }
        else{
          #if query has no hits, switch variable is given value off to identify that output should not be displayed         
          return(list("output_switch"="off"))       
        }
      }
      else{
        #if query empty,switch variable is given value off to identify that output should not be displayed     
        return(list("output_switch"="off"))
        
      }
    }    
#make blast run when you click blast button using eventReactive Function
#"The function eventReactive() is used to compute a reactive value that only updates in response to a specific event."-https://campus.datacamp.com/courses/building-web-applications-with-shiny-in-r/reactive-programming-3?ex=11
    run_sequence<-eventReactive(input$blast,{make.comparison(input$mysequence)})
    

#renderDataTable selection mode='single' ensures only one row can be selected at a time
#selected=1 means first row will be selected as default
#in options scrollX = TRUE means horixontal scrolling is enabled 
#pageLength=7 refers to how many rows to show per page    
# dom refers to which elements of datatable to include , here dom='tp', refers to table and pagination control (other options include a search bar etc)
#see https://datatables.net/reference/option/dom
#row.names not included and colnames set here too (dont need to do this if happy with original dataframes rownames )
    output$Main_output_table=DT::renderDataTable({
      run_sequence_output=run_sequence()
#if output_switch on      
      if(run_sequence_output$output_switch=="on"){
        run_sequence_output$tax_and_stats
      }
      },selection = list(mode='single',selected=1),options=list(scrollX=TRUE,pageLength=7,dom='tp'),rownames=FALSE)

    
    
    output$OTU_map=renderPlot({
      run_sequence_output<-run_sequence()
#if output_switch on      
      if ( run_sequence_output$output_switch=="on"){
#get selected row       
        s=input$Main_output_table_rows_selected
#if row has been selected        
        if(length(s)){
#get OTU from  run_sequence output          
          OTU=run_sequence_output$tax_and_stats[s,2]
          par(mar = c(4, 4, 1, 4))
          dbGetQuery(con, "set standard_conforming_strings to 'on'")         
#get relevant map       
          unescape_bytea_map=dbGetQuery(con, paste("SELECT map_object FROM <--specified_Schema_table_prefix-->_otu_attributes.<--specified_Schema_table_prefix-->_maps WHERE hit='",toString(OTU),"';",sep=""))
#convert from bytea to r object          
          bytea_map<-postgresqlUnescapeBytea(unescape_bytea_map) 
          r_object_map<-unserialize( bytea_map)
#plot       		
          spplot(r_object_map[[1]]["var1.pred"],at=unlist(r_object_map[2:11]),xlab=toString(OTU),ylab.pos=c(5,10,100),sp.layout=list("sp.lines",uk.line,lwd=2,col="black"))
        }
      }   
     })
    output$OTU_avc_boxplot = renderPlot({    
      run_sequence_output<-run_sequence()
      if ( run_sequence_output$output_switch=="on"){
        #get selected row       
        s=input$Main_output_table_rows_selected
        #if row has been selected        
        if(length(s)){
          #get OTU from  run_sequence output     
          abund=run_sequence_output$OTU_abundance[,s,drop=FALSE]
        #if otu has occupancy of  greater than 30
         # if(length(which(abund!=0))>20){
            #merge with env
            abund_env=merge(abund,env,by.x=0,by.y=1)
             #order avc factor levels so x axis represents habitat gradient
            ord<-reorder(abund_env$avc,X=as.numeric(abund_env$avc_code),FUN=mean)
            par(mar=c(11,5.8,2,0.4)+0.1)
            #outline =TRUE show outliers
            boxplot(abund_env[,2]~ord,las=2,ylab="",xlab="",ylim=c(0,1.1*max(abund_env[,2])),cex=0.5,outline=TRUE)
            title(ylab=paste("Relative Abundance (",colnames(abund_env)[2],")",sep=""), line=4.5, cex.lab=0.9)
           # }
          #else{
          #  par(mar = c(0,0,0,0))
          #  empty_plot=plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
          #  empty_plot+ text(x = 0.8, y = 0.8, paste("Not enough data to generate box plot"),
          #  cex = 1, col = "red", family="sans", font=1, adj=1)
         # }
         }
      }else{plot.new()}
    })

#selection='none' means no rows can be selected     
# options dom='t' is telling renderDataTable just to show table and not include fancy extras
#No rownames
    output$OTU_blast_match=DT::renderDataTable({
      run_sequence_output<-run_sequence()
#if output_switch on      
      if (run_sequence_output$output_switch=="on"){
#get selected row       
       s=input$Main_output_table_rows_selected
#if row has been selected        
       if (length(s)){
#get relevant         
          run_sequence_output$full_blast_output[s,]
        }
      }
    },selection='none',options=list(dom='t'),rownames=FALSE,colnames= c("subject id","% identity","alignment length","mismatches","gap opens","q.start","q.end","s.start","s.end","evalue","bit score"))
  
#if clear input button is pressed, query is set back to orignal value i.e blank  https://www.rdocumentation.org/packages/shinyjs/versions/2.1.0/topics/reset 
    observeEvent(input$clearInput, {
      reset("mysequence")
    })
###SEQUENCE###    
#if example sequence is button is pressed, example sequence entered into query box
    observeEvent(input$exampleSequence, {
      updateTextInput(session,"mysequence",value="<--specified_Example_sequence-->")
    })     
#display warning if switch set to off  
    output$Warning <- renderText({ 
      run_sequence_output<-run_sequence()
      if(run_sequence_output$output_switch=="off"){
        paste("No Hits Found!")
      }
      
    })
    
        
}

# Run the application 
shinyApp(ui = ui, server = server)
