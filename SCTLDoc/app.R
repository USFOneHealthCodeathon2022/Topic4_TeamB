library(shiny)
library(leaflet)
library(phyloseq)
library(dada2)

## contains physeq_16S; a phyloseq object with count table, metadata, and taxonomy information bundled together
load('reference_data.RData')

subset_samples <- function(physeq_16S, criteria) {
  ## create new set of samples for plotting based on interactive selection criteria
  subset_samples(physeq_16S, sample_type == "TissueSlurry" | sample_type == "Mucus" | sample_type == "TissueSlurry_Skleton" |
                   sample_type=="Seawater" | sample_type=="Sediment")
}

plot_map <- function(ps) {
  ## generate plot for panel 1
  ## https://rstudio.github.io/leaflet/shiny.html
  filtDat <- sample_data(ps)
  filtDat$latitude <- as.numeric(filtDat$latitude)
  filtDat$longitude <- as.numeric(filtDat$longitude)
  filtDat <- filtDat[!(is.na(filtDat$latitude) | is.na(filtDat$longitude))]
  leaflet() %>%
    addCircleMarkers(
      data = filtDat,
      lat = ~latitude,
      lng = ~longitude
    )
  ##data = map_click(),
}

plot_taxonomy <- function(ps) {
  ## generate plot for panel 2
  sum_ps <- ps  %>%
    tax_glom(taxrank = "Family") %>%
    transform_sample_counts(function(x) {x/sum(x)}) %>% # Transform to rel. abundance
    psmelt()
  
  sum_ps$tissue_type <- factor(sum_ps$tissue_type, 
                               levels = c("AH","DU", "DL"))
  ggplot(subset(sum_ps, Abundance>0.5 ), 
         aes(x=sample_type, y=Abundance, fill=Order)) + 
    geom_bar(stat="identity", position="fill") +
    facet_grid(tissue_type~., scales = "free") +
    scale_fill_manual(values=c("#56B4E9","#CBD588","#5F7FC7", "orange","#DA5724","#CD9BCD",
                               "gray80", "#AD6F3B", "#673770","#D14285", "#652926","#8569D5", 
                               "#5E738F","#D1A33D", "#8A7C64","lightsalmon","aquamarine4",
                               "lightblue4", "lightpink", "ivory4","royalblue4", "darkorchid", 
                               "palevioletred1", "#56B4E9","#CBD588","yellow2","#5F7FC7", "orange","#DA5724",
                               "#CD9BCD", "gray80",
                               "#AD6F3B", "#673770","#D14285", "#652926","#8569D5", "#5E738F",
                               "#56B4E9","#CBD588","#5F7FC7", "orange","#DA5724","#CD9BCD", "blue", "red")) +
    guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.70, ncol=1)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme_classic()
}

plot_ordination <- function(ps) {
  ## generate plot for panel 3
  ps_clr <- microbiome::transform(ps, 'clr')
  psr_clr.ord <- ordinate(ps_clr, "RDA", "euclidean")
  plot_ordination(ps_clr, psr_clr.ord,
                  color="tissue_type",
                  axes = c(1,2)) +
    theme_classic() +
    theme(legend.text = element_text(size=6),
          legend.title = element_text(size=7)) +
    theme(axis.title.x = element_text(size = 7)) +
    theme(axis.title.y = element_text(size = 7)) +
    theme(axis.text.x = element_text(size = 6)) +
    theme(axis.text.y = element_text(size = 6)) +
    scale_color_manual(values=c("#007f00", "Navy", "#800020"),
                       labels=c(AH="Apparently Healthy",
                                DU="Diseased Unaffected", 
                                DL="Diseased Lesion")) +
    stat_ellipse()
}

plot_boxes <- function() {
  ## generate plot for panel 4
  
  
}

# increase max R-Shiny user-input file size from 5 MB to 3 GB
options(shiny.maxRequestSize = 3 * 1024 ^ 3)

# the ui object has all the information for the user-interface
ui <- fluidPage(
  fluidRow(
    fileInput("user_reads", "Upload merged fastqs",
              multiple = FALSE,
              accept = c("*.fastq","*.fastq.gz")),
    fileInput("user_disease", "Upload disease annotations",
              multiple = FALSE,
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv"))
  ),
  fluidRow(
    column(6,leafletOutput("panel1", height = 300)),
    column(6,plotOutput("panel2", height = 300))
  ),
  fluidRow(
    column(6,plotOutput("panel3", height = 300)),
    column(6,plotOutput("panel4", height = 300))
  )
)

server <- function(input, output) {
  p <- reactiveValues(panel1=NULL,panel2=NULL,panel3=NULL,panel4=NULL)
  observeEvent(input$user_disease, {
    #user_counts <- makeSequenceTable(dada(input$user_reads$datapath, pool='pseudo', selfConsist=TRUE, err=NULL))
    #user_phyloseq <- phyloseq(otu_table(user_counts, taxa_are_rows = T), sample_data(read.table(input$user_disease$datapath)))
    #merged_data <- merge_phyloseq(physeq_16S,user_phyloseq)
    merged_data <- physeq_16S
    subsetted_data <- merged_data
    #criteria <- list() ## create a set of subsetting criteria based on user input
    #subsetted_data <- subset_samples(merged_data, criteria)
    p$panel1 <- plot_map(merged_data)
    p$panel2 <- plot_taxonomy(subsetted_data)
    p$panel3 <- plot_ordination(subsetted_data)
    p$panel4 <- plot_boxes(subsetted_data)
  })
  output$panel1 <- renderLeaflet(if(is.null(p$panel1)) return() else {p$panel1})
  output$panel2 <- renderPlot(if(is.null(p$panel2)) return() else {p$panel2})
  output$panel3 <- renderPlot(if(is.null(p$panel3)) return() else {p$panel3})
  output$panel4 <- renderPlot(if(is.null(p$panel4)) return() else {p$panel4})
}

shinyApp(ui, server)
