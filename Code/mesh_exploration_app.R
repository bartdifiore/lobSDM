#----------------------------------
## Libraries and preliminaries
#----------------------------------
library(tidyverse)
library(shiny)
library(bslib)
library(sdmTMB)
library(ggplot2)

#----------------------------------
## Data
#----------------------------------
all_mod_data<- readRDS(here::here("Data/Derived/all_model_data.rds"))

# Focus on just lobster, and during GLORYs time series
year_min <- 1993
year_max <- 2023

red_mod_data <- all_mod_data |>
  filter(life_class == "juvenile") |>
  filter(between(year, year_min, year_max)) 

summary(red_mod_data)

# Still some weird NA biomass values to figure out, dropping those for now
mod_data<- red_mod_data |>
  drop_na(total_biomass)
summary(mod_data)

# What the heck is going on with that 2014 value?
t<- mod_data[which.max(mod_data$total_biomass),]
plot(mod_data$total_biomass)

# Doesn't seem like there is anyway that point is real.
mod_data<- mod_data[-which.max(mod_data$total_biomass),]
plot(mod_data$total_biomass)


mod_data <- mod_data |>
  sdmTMB::add_utm_columns(ll_names = c("longitude", "latitude"), units = "km") 

mod_data_unique<- mod_data |>
    distinct(longitude, latitude)

#----------------------------------
## App
#----------------------------------

ui <- page_sidebar(
  title = "sdmTMB mesh explorer",
  sidebar = sidebar(
    # Mesh method selection
    selectInput("mesh_type", "Mesh creation method",
      choices = c("cutoff", "cutoff_search", "kmeans"),
      selected = "cutoff"),
    
    # Mesh parameters that change based on method
    conditionalPanel(
      condition = "input.mesh_type == 'kmeans'",
      numericInput("n_knots", "Number of knots", 
        value = 100, min = 10, max = 500, step = 10)
    ),
    conditionalPanel(
      condition = "input.mesh_type != 'kmeans'",
      numericInput("cutoff", "Distance cutoff", 
        value = 600, min = 5, max = 600, step = 5)
    ),
    
    # Add a download button
    downloadButton("downloadMesh", "Download Mesh")
  ),

  # Layout stuff
  # layout_columns(
  #   fill = FALSE,
  #   value_box(
  #     title = "Number of Vertices",
  #     value = textOutput("n_vertices")
  #   ),
  #   value_box(
  #     title = "Number of Triangles",
  #     value = textOutput("n_triangles")
  #   )
  # ),
  card(
    card_header("Mesh Visualization"),
    plotOutput("mesh_plot", height = "600px")
  ),
    # Main panel with mesh plot and statistics
  card(
    card_header("Mesh Statistics"),
    tableOutput("mesh_stats")
  )
)

server <- function(input, output) {
  mesh <- reactive({
    # Force reactivity to key parameters
    input$mesh_type
    input$n_knots
    input$cutoff
    
    # Set parameters
    mesh_args <- list(
      data = mod_data,
      xy_cols = c("longitude", "latitude"),
      type = input$mesh_type
    )
    
    # Add additional params based on mesh type
    if (input$mesh_type == "kmeans") {
      mesh_args$n_knots = input$n_knots
    } else {
      mesh_args$cutoff <- input$cutoff
    }
    
    do.call(make_mesh, mesh_args)
  })

  output$mesh_plot <- renderPlot({
    req(mesh())

    plot(mesh()$mesh,
      main = "",
      asp = 1
    )
    # points(mod_data_unique$longitude, mod_data_unique$latitude, col = "#2b8cbe", pch = 16, size = 0.5, alpha = 0.5)
  })

  # Output mesh statistics
  output$mesh_stats <- renderTable({
    req(mesh())

    data.frame(
      Statistic = c("Number of vertices", "Number of triangles"),
      Value = c(
        length(mesh()$mesh$loc[, 1]),
        nrow(mesh()$mesh$graph$tv)
      )
    )
  })

  # Download handler
  output$downloadMesh <- downloadHandler(
    filename = function() {
      paste0("mesh_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
    },
    content = function(file) {
      saveRDS(mesh(), file)
    }
  )
}

shinyApp(ui, server)
