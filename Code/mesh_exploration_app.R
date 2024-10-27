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
  title = "Mesh Explorer for sdmTMB",
  sidebar = sidebar(
    numericInput("cutoff", "Mesh cutoff",
      value = 10000, min = 10, max = 200, step = 10
    ),
    downloadButton("downloadMesh", "Download Mesh")
    # numericInput("n_knots", "Number of knots",
    #   value = 50, min = 10, max = 200, step = 10
    # ),
    # numericInput("type", "Type (1 = cutoff, 2 = kmeans)",
    #   value = 1, min = 1, max = 2, step = 1
    # )
  ),
  layout_columns(
    fill = FALSE,
    value_box(
      title = "Number of Vertices",
      value = textOutput("n_vertices")
    ),
    value_box(
      title = "Number of Triangles",
      value = textOutput("n_triangles")
    )
  ),
  card(
    card_header("Mesh Visualization"),
    plotOutput("mesh_plot", height = "600px")
  )
)

server <- function(input, output) {
  mesh <- reactive({
    sdmTMB::make_mesh(mod_data, 
    xy_cols = c("longitude", "latitude"), 
    type = "cutoff", 
    cutoff = input$cutoff, 
    fmesher_func = fmesher::fm_mesh_2d_inla)
  })

  output$n_vertices <- renderText({
    length(mesh()$mesh$n)
  })

  output$n_triangles <- renderText({
    nrow(mesh()$mesh$graph$tv)
  })

  output$mesh_plot <- renderPlot({
    current_mesh <- mesh()
    plot(current_mesh$mesh,
      main = "",
      asp = 1
    )
    # points(mod_data_unique$longitude, mod_data_unique$latitude, col = "#2b8cbe", pch = 16, size = 0.5, alpha = 0.5)
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
