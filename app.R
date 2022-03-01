library(shiny)

subset_samples <- function(physeq_16S, criteria) {
  ## create new set of samples for plotting based on interactive selection criteria
  subset_samples(physeq_16S, sample_type == "TissueSlurry" | sample_type == "Mucus" | sample_type == "TissueSlurry_Skleton" |
                   sample_type=="Seawater" | sample_type=="Sediment")
}

plot_map <- function() {
  ## generate plot for panel 1
  ## https://rstudio.github.io/leaflet/shiny.html
  leaflet() %>%
    addCircleMarkers(
      data = map_click(),
      lat = ~latitude,
      lng = ~longitude
    )
}

plot_taxonomy <- function(ps, taxrank) {
  ## generate plot for panel 2
  sum_ps <- ps  %>%
    tax_glom(taxrank = "Family") %>%
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
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
    guides(fill = guide_legend(keywidth = 0.5, , keyheight =.70, ncol=1)) +
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
    theme(legend.text =element_text(size=6),
          legend.title=element_text(size=7)) +
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


# increase max R-Shiny user-input file size from 5 to 30 MB
options(shiny.maxRequestSize = 30 * 1024 ^ 2)

# define fonts for plot
f1 <- list(family = "Arial, sans-serif",
           size = 24,
           color = "black")

f2 <- list(family = "Arial, sans-serif",
           size = 20,
           color = "black")

f3 <- list(family = "Arial, sans-serif",
           size = 16,
           color = "black")

# set plot margins
m <- list(
  l = 20,
  r = 20,
  b = 10,
  t = 100,
  pad = 4
)

# the ui object has all the information for the user-interface
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      plotOutput("panel2")
    ),
    mainPanel(
      plotOutput("panel3")
    )
))

server <- function(input, output, session) {
  #taxrank load this as pre-computed object merged with user-input data
  #physeq_16S load this as pre-computed object merged with user-input data
  criteria <- list() ## create a set of subsetting criteria based on user input
  subsetted_data <- subset_samples(physeq_16S, criteria)
  output$panel2 <- plot_taxonomy(subsetted_data, taxrank)
  output$panel3 <- plot_ordination(subsetted_data)
}

shinyApp(ui, server)
