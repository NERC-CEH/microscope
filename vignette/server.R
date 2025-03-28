
#Pointing R to where libraries are 
#.libPaths("/home/brijon/R/x86_64-pc-linux-gnu-library/3.3/")
library(shiny)
#Data table
library(DT)
#to connect to SQL
library(RPostgreSQL)
library(pool)
#for pH HOF model loading
library(eHOF)
library(reshape)

library(vegan)
#library(maptools)
library(RColorBrewer)
library(ggplot2)
library(gstat)
library(sp)

#============================================================================================================================================================== 
#Define server

server <- function(input, output,session) {
	
#connect to postgres 
      
	con <- dbPool(
 	 drv = RPostgreSQL::PostgreSQL(max.con=140),
 	 dbname = "molecular",
 	 host = "connect-apps.ceh.ac.uk",
 	 user = getOption("trait_matrix_userid"),
	 password = getOption("trait_matrix_password"),
 #time out is clearly not working
 #no of milliseconds way to small anyway... not sure if function changed since I fist wrote this
 #not a major issue as connections closed each day on server anyway... but need to sort this if we do end up having a lot of users 
         idleTimeout=1800

	)
#=====================================================================================================================================
#shiny JS commands
#Get info about app , shiny js allows the information to appear and dissapear when clicking more information
           	
  shinyjs::onclick("more_info_button",shinyjs::toggle(id="more_info",anim=TRUE))

#get results
  shinyjs::onclick("blast",shinyjs::show(id="Results",anim=TRUE))
#====================================================================================================================================  
#Get AVC data 

  SQL_command=paste("select * from env_attributes.env_attributes;")
  env <- dbGetQuery(con, SQL_command)
#Get mapping outline
  dbGetQuery(con, "set standard_conforming_strings to 'on'")
  SQL_command=paste("select plot_object from plotting_tools.map_tools where description= 'map_outline';")
  uk.line <-unserialize(postgresqlUnescapeBytea( dbGetQuery(con, SQL_command)))
 	
  rownames(env)<-env[,1]
  env[env=="NA"]<-NA
  
#==============================================================================================================================================================
#Define function that runs blast ,arguments are database and query sequence

  make.comparison <- function(db, query ){
   
#if query isnt empty then pass to blast
   if (query!=""){    
#blast_command 
   cmd <- paste('echo -e ">Name\n',query,'"', '|blastn' ,'-db',
'blast_db/repseqs_ID_TaxER', '-num_descriptions',20,'-num_alignments',20,'-evalue',0.001, '-outfmt', 7)
   
#run system command and Capture Output have to specify bash because server running bourne shell... 
    blast_capture<- system(paste("/bin/bash -c", shQuote(cmd)),intern=TRUE)
#blast_capture<-blast_capture[which(duplicated(str_split_fixed(blast_capture, fixed("\t"),n=12)[,2],incomparables=FALSE)==FALSE)]
#if there are hits
    if(blast_capture[4]!="# 0 hits found"){
     output_switch <-"on"
#remove first and last lines of output (descriptive not hits)
      blast_capture<-blast_capture[-length(blast_capture)]
      blast_capture<-blast_capture[6:25]
#some dont have 20 hits! so nas appear....
       blast_capture<-blast_capture[!is.na(blast_capture)]
    
#==============================================================================================================================================================    
    #make empty arrays/dataframes
      all_blst_output<-c()
      OTUlist<-c()

#empty blast output matrix
      m <- matrix(0, ncol = 9)
 
#convert to dataframe
      ph.model<-data.frame(m)  
#remove row of 0s ..have to do this dont delete(again!)
      ph.model<-ph.model[-1,]
#tax dataframe
      m<-matrix(0, ncol=8)
      tax<-data.frame(m)
#remove first rows of 0's      
      tax<-tax[-1,]
 #get col names for abundance data
      myQuery <- "SELECT * FROM information_schema.columns WHERE table_schema = 'otu_abund' AND table_name = 'otu_abund_ur'"
      samp_colnames <- dbGetQuery(con, myQuery)[,4]
 #make empty matrix with sample names as colnames     
      m<-matrix(0,ncol=1007)
      abund<-data.frame(m)
      abund<-abund[-1,]
      colnames(abund)<-samp_colnames

#==============================================================================================================================================================   
#fill arrays and df
      for (blast_hit in blast_capture){
#get OTU
        OTU=strsplit(blast_hit,"\t")[[1]][2]
# look up OTU in postgres to get model info
        SQL_command=paste("SELECT bm.hit, bm.hof_model, bm.optimum1, bm.optimum2, bm.reclass1, bm.reclass2, abnd.abundance_rank, abnd.occupancy_proportion FROM bacterial_otu_attributes.bacterial_ph_ehof_models bm
INNER JOIN bacterial_otu_attributes.bacterial_abundance_stats abnd ON bm.hit = abnd.hit WHERE bm.hit='",toString(OTU),"';",sep="")
        ph.model.row <- dbGetQuery(con, SQL_command)
        ph.model.row$Blast_Percentage_Identity<-strsplit(blast_hit,"\t")[[1]][3]
#Attach to ph.model df
        ph.model<-rbind(ph.model,ph.model.row[,c(9,1:8)])
        ph.model<- subset(ph.model, !duplicated(ph.model$hit))
        row.names(ph.model)<-ph.model$hit
 #look up otu in postgres to get tax info	     
        SQL_command=paste("SELECT * FROM bacterial_otu_attributes.bacterial_taxonomy WHERE hit='",toString(OTU),"';",sep="")
        tax.row<-dbGetQuery(con,SQL_command)
#Attatch to tax df
        tax<-rbind(tax,tax.row)
        tax<- subset(tax, !duplicated(tax$hit))
#look up otu in postgres to get abundance info    
        SQL_command=paste("SELECT * FROM abund_tables.bacterial_abund WHERE hit='",toString(OTU),"';",sep="")
        abund.row<-dbGetQuery(con,SQL_command)
#Attatch to abund df      
        abund<-rbind(abund,abund.row)
        abund<- subset(abund, !duplicated(abund$hit))
     
        all_blst_output<-c(all_blst_output,strsplit(blast_hit,"\t"))
   
        OTUlist=c(OTUlist,OTU)
    
      }
#==============================================================================================================================================================   
#remove replicates
      OTUlist<-unique(OTUlist)
#make blast output into nice df for the GUI 
      all_blst_output_df<-t(as.data.frame(all_blst_output))[,2:12]
      row.names(all_blst_output_df)<-all_blst_output_df[,1]
#make colnames
      colnames(all_blst_output_df)<-c('subject id',' % identity', 'alignment length', 'mismatches', 'gap opens', 'q. start', 'q. end', 's. start', 's. end', 'evalue', 'bit score')
      all_blst_output_df<-as.data.frame(all_blst_output_df[OTUlist,])
#assign rownames to tax df	
      row.names(tax)<-tax$hit
#assign colnames to taxonomy df
      colnames(tax)=c("subject id","kingdom","phylum","class","order","family","genus","species")
      tax<-as.data.frame(tax)
      ph.model<-ph.model[!duplicated(ph.model[order(-as.numeric(ph.model$Blast_Percentage_Identity)),]),]
    
      colnames(ph.model)<-c("Blast Percentage Identity","CS OTU hit","pH HOF model","pH Optimum 1","pH Optimum 2","Model description","pH Class","Abundance rank","Occupancy")
      ph.model<-as.data.frame(ph.model)
      all_blst_output_df<-all_blst_output_df[row.names(ph.model),]
   
      abund<-t(abund)
      colnames(abund)<-abund[1,]
      abund<-abund[-1,]
      return(list(ph.model,all_blst_output_df,tax,abund,output_switch))
    }
    else{
         
         output_switch <-"off"
       
         return(output_switch)
        
    }}
    else{
    
       output_switch <-"off"
   
      
        return(output_switch)

   }}
 #===================================================================================================================================================================================== 
 #make blast run when you click blast button should probably merge this with the java script function that reveals output for neatness!!!
 
   run_sequence<-eventReactive(input$blast,{make.comparison("/home/brijon/repseqs_db", input$mysequence) })
  

 #=====================================================================================================================================================================================
 #Define table ,selection mode is one at a time (dont want to show several plots at once etc), automatically select top row otherwise other outputs
# will be empty when page loads which is just ugly :)
  output$blastout<-DT::renderDataTable({
   
   all_mod_output<-run_sequence()
   if(all_mod_output[[length(all_mod_output)]]=="on"){
     ph_mod_output<-all_mod_output[[1]]
   }
   

   }, selection = list(mode='single',selected=1),options=list(scrollX=TRUE,pageLength=7,dom='tp'),rownames=FALSE, colnames = c("Blast Percentage Identity","CS OTU hit","pH HOF model","pH Optimum 1","pH Optimum 2","Model description","pH Class","Abundance rank","Occupancy")
  ) 
  

#================================================================================================================================================================================================
#HOF plot code :)
  output$modelplot = renderPlot({
    all_mod_output<-run_sequence()
#selected rows 
     if(all_mod_output[[length(all_mod_output)]]=="on"){
     
      
       s=input$blastout_rows_selected
#get selected rows in df
       if (length(s)){
         OTU=all_mod_output[[2]][s,1]
         par(mar = c(4, 4, 1, 4))
        dbGetQuery(con, "set standard_conforming_strings to 'on'")
        ph.mod.obj=dbGetQuery(con, paste("SELECT mod_object FROM bacterial_otu_attributes.bacterial_ph_ehof_models WHERE hit='",toString(OTU),"';",sep="")) 
#unserialize 
         ph.mod.obj<-postgresqlUnescapeBytea(ph.mod.obj)
         ph.mod.obj<- unserialize(ph.mod.obj)
      
#Own custom plot with relative abundance as modifying plot hof object is a nightmare :)
         mod_choice<-Para(ph.mod.obj)$model
         fitted=ph.mod.obj$models[mod_choice][[1]][9] 
         mod_stats=as.data.frame(cbind(ph.mod.obj$y,ph.mod.obj$x,fitted$fitted))
         mod_stats=mod_stats[order(mod_stats$V2),]
         mod_col=c("black","red","#2BF33F","#3554EE","#895A3B")
         names(mod_col)=c("I","II","III","IV","V")
         plot(mod_stats$V2,mod_stats$V1,ylim=c(0,quantile(mod_stats$V1,0.999)),cex=0.7,pch=19,xlab="pH",ylab=paste("Number of reads ","(",toString(OTU),")",sep=""))
         lines(mod_stats$V2,mod_stats$V3,ylim=c(0,quantile(mod_stats$V1,0.999)),lwd=2,col=mod_col[mod_choice])
     
      
      
       }
     }
     
     
   else{
         plot.new()
    }
#end of render plot   
  })


#================================================================================================================================================================================================
#LOESS plot code 
  output$loessplot = renderPlot({
    all_mod_output<-run_sequence()
#selected rows 
   if(all_mod_output[[length(all_mod_output)]]=="on"){
       s=input$blastout_rows_selected
      if(length(s)){
         par(mar = c(4, 4, 1, 4)) 
         abund <-all_mod_output[[4]]
         OTU.spc<-as.data.frame(as.numeric(abund[,s]),row.names=row.names(abund))
         OTU.spc_pH<-merge(OTU.spc,env,by=0)
#drop unecessary collumns
         OTU.spc_pH<-OTU.spc_pH[,-c(3:5)]
#order df by ph
         OTU.spc_pH<-OTU.spc_pH[order(OTU.spc_pH[,3]),]
         loessMod50 <- loess(OTU.spc_pH[,2]~OTU.spc_pH[,3], span=0.5)
         smoothed50 <- predict(loessMod50,se=TRUE) 
         plot(OTU.spc_pH[,3],OTU.spc_pH[,2],ylim=c(min(OTU.spc_pH[,2]),quantile(OTU.spc_pH[,2],0.999)),cex=0.7,pch=19,xlab="pH",ylab=paste("Relative Abundance (",colnames(abund)[s],")",sep=""),col="#878273")
         lines(smoothed50$fit, x=OTU.spc_pH[,3],ylim=c(min(OTU.spc_pH[,2]),quantile(OTU.spc_pH[,2],0.999)), col="black",lwd=2)
      }
    }     
  })
      
      
      
      
#================================================================================================================================================================================================
#user submitted section (indicators)

    output$Indicators=renderText({
         all_mod_output<-run_sequence()
         if(all_mod_output[[length(all_mod_output)]]=="on"){
         s=input$blastout_rows_selected
#get selected rows in df
         if (length(s)){
                 OTU=all_mod_output[[2]][s,1]
                 SQL_command=paste("select * from bacterial_otu_attributes.bacterial_user_submitted_attributes WHERE hit='",toString(OTU),"';",sep="")
                 indicator.tab<- dbGetQuery(con, SQL_command)
                
#if there are indicators for this otu                
                if (nrow(indicator.tab)>=1){
#empty variable that will be printed in indicators text box                
                  indicators=""
#loop over difference datasets in reference collumn
                  for (study in unique(indicator.tab$reference)){
#get indicator subsets  for each dataset                  
                    study_subset=indicator.tab[which(indicator.tab$reference==study),]
#sort by identity and pvalue
                    study_subset=study_subset[order(-study_subset$identity,study_subset$pvalue),]
                    indicator=paste(study_subset[1,4],"<b>  p value: </b>",study_subset[1,5],"\n <b> Source: </b>",study_subset[1,6],sep="") 
                    #indicators=paste(indicators,indicator,sep="\n")
#paste with line breaks
                    indicators=paste(indicator,indicators,sep='<br/>')
                  }
                  indicators=paste("<center> <b>",OTU,"</center> </b> <center>",indicators,"</center>",sep="")
   		            HTML(indicators)
                 }else{
                 noindicators= paste("<center> <b>",OTU," </b> </center> \n <center> No user submitted information.</center>",sep="")
                 HTML(noindicators)}
         }
         }

  })





#===================================================================================================================================================================================================================   
#Map

  output$map=renderPlot({
    all_mod_output<-run_sequence()
  
    if (all_mod_output[[length(all_mod_output)]]=="on"){
       	s=input$blastout_rows_selected
       if(length(s)){
	      	OTU=all_mod_output[[2]][s,1]
        
          par(mar = c(4, 4, 1, 4))
          dbGetQuery(con, "set standard_conforming_strings to 'on'")
          map.obj=dbGetQuery(con, paste("SELECT map_object FROM bacterial_otu_attributes.bacterial_maps WHERE hit='",toString(OTU),"';",sep=""))
          raw_retrieved_map<-postgresqlUnescapeBytea(map.obj) 
       		object_retrieved_map<-unserialize(raw_retrieved_map)
          spplot(object_retrieved_map[[1]]["var1.pred"],at=unlist(object_retrieved_map[2:11]),xlab=toString(OTU),ylab.pos=c(5,10,100),sp.layout=list("sp.lines",uk.line,lwd=2,col="black"))
       	}
    }
})

#===================================================================================================================================================================================================================   

    output$AVC_box_plot = renderPlot({    
    all_mod_output<-run_sequence()
#selected rows 
     if(all_mod_output[[length(all_mod_output)]]=="on"){
       s=input$blastout_rows_selected
       if(length(s)){
        abund <-all_mod_output[[4]]
         OTU.spc<-as.data.frame(as.numeric(abund[,s]),row.names=row.names(abund))
         numb_otus<-specnumber(OTU.spc)
         if (sum(numb_otus)>150){

          ord<-reorder(env$avc,X=as.numeric(env$avc_code),FUN=mean)
          par(mar=c(11,5.8,2,0.4)+0.1)
          boxplot(OTU.spc[,1]~ord,las=2,ylab="",cex=0.5,outline=FALSE,col="white",xlab="")
          title(ylab=paste("Relative Abundance (",colnames(abund)[s],")",sep=""), line=4.5, cex.lab=0.9)
      
         }else{
      
          par(mar = c(0,0,0,0))
          empty_plot=plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
          empty_plot+ text(x = 0.8, y = 0.8, paste("Not enough data to generate box plot"),
          cex = 1, col = "red", family="sans", font=1, adj=1)

         }
       }
     }
    else{
      plot.new()
     

    }
   })

#============================================================================================================================================
#Taxonomy
  
    output$OTU.Taxon<-DT::renderDataTable({
      all_mod_output<-run_sequence()
      if(all_mod_output[[length(all_mod_output)]]=="on"){
        s=input$blastout_rows_selected
        if (length(s)) {
           as.data.frame(all_mod_output[[3]])[s,]
        
        }
      }
#no selection and just tables no extras like search(dom=t)
    },selection='none',options=list(dom='t'),rownames=FALSE)
    

#============================================================================================================================================  
#blast output
    

   output$BlastResults<-DT::renderDataTable({
      #dat_raw
     all_mod_output<-run_sequence()
       if(all_mod_output[[length(all_mod_output)]]=="on"){
         s=input$blastout_rows_selected
      
         if (length(s)){
     
          all_mod_output[[2]][s,]
     #end of second if statement    
         }
   #end of first if statement
     }
     #end of render datatable 
    },selection='none',options=list(dom='t'),rownames=FALSE)
    
    

    
    output$Warning <- renderText({ 
      all_mod_output<-run_sequence()
      if(all_mod_output[[length(all_mod_output)]]!="on"){
       
        paste("No Hits Found!")
 
        }
     
    }  )
   
    observeEvent(input$resetSequence, {
       # updateTextInput(session,"mysequence",value="")
        reset("mysequence")
      })


    observeEvent(input$exampleSequence, {
        updateTextInput(session,"mysequence",value="ACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGTGATGCAAGTCTGGTGTGAAATCTCGGGGCTCAACTCCGAAATTGCACCGGATACTGCGTGACTCGAGGACTGTAGAGGAGATCGGAATTCACGGTGTAGCAGTGAAATGCGTAGATATCGTGAGGAAGACCAGTTGCGAAGGCGGATCTCTGGGCAGTTCCTGACACTGAGGCACGAAGGCCAGGGGAGCAAACGGG")
       # reset("mysequence")
    })

   
}
#=================================================================================================================================================


     
