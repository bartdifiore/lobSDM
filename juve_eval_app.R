library(shiny)
library(bslib)
library(ggplot2)
library(DHARMa)
library(dplyr)

options(shiny.autoreload = TRUE)

# Safely read in the saved model fits with error handling
safely_read_rds <- function(file_path) {
  tryCatch(
    readRDS(file_path),
    error = function(e) NULL
  )
}

fit_base <- safely_read_rds("Juve_Base.rds")
fit_sp <- safely_read_rds("Juve_Sp.rds")
fit_spst <- safely_read_rds("Juve_SpST.rds")

model_list <- list("Baseline" = fit_base, "Spatial" = fit_sp, "Spatio-Temporal" = fit_spst)

# Check if models were loaded successfully
if (all(sapply(list(fit_base, fit_sp, fit_spst), is.null))) {
  stop("Could not load any model files. Please ensure the RDS files are in the correct directory.")
}

ui <- fluidPage(
  title = "Juvenile lobster SDM evaluation",
  sidebar = sidebar(
    selectInput("model_select", "Select Model:", 
                choices = names(model_list)),
    radioButtons("plot_type", "Diagnostic Plot:",
                choices = c("Residuals vs. Fitted" = "resid_fitted",
                          "QQ Plot" = "qq",
                          "Scale-Location" = "scale_loc",
                          "DHARMa Residuals" = "dharma")),
    hr(),
    h4("Model Comparison"),
    tableOutput("model_metrics")
  ),
  
  layout_columns(
    fill = FALSE,
    # value_box(
    #   title = "Deviance Explained",
    #   value = textOutput("dev_explained"),
    #   showcase = bsicons::bs_icon("bar-chart")
    # ),
    value_box(
      title = "AICc",
      value = textOutput("aicc"),
      showcase = bsicons::bs_icon("graph-up")
    )
  ),
  card(
    full_screen = TRUE,
    card_header("Diagnostic Plots"),
    plotOutput("diagnostic_plot", height = "500px")
  )
)

server <- function(input, output, session) {
  selected_model <- reactive({
    model_list[[input$model_select]]
  })
  
  output$diagnostic_plot <- renderPlot({
    req(input$plot_type)
    
    model <- selected_model()
    
    if(input$plot_type == "dharma") {
      dharma_residuals <- DHARMa::createDHARMa(
        simulatedResponse = t(model$family$simulate_response(model$tmb_obj$env$last.par.best)),
        observedResponse = model$data[[model$response]],
        integerResponse = model$family$family %in% c("poisson", "nbinom2")
      )
      plot(dharma_residuals)
    } else {
      resids <- residuals(model)
      fitted_vals <- fitted(model)
      
      if(input$plot_type == "resid_fitted") {
        ggplot(data.frame(fitted = fitted_vals, resids = resids), 
               aes(x = fitted, y = resids)) +
          geom_point() +
          geom_hline(yintercept = 0, linetype = "dashed") +
          theme_minimal() +
          labs(x = "Fitted values", y = "Residuals",
               title = "Residuals vs Fitted")
      } else if(input$plot_type == "qq") {
        ggplot(data.frame(resids = resids), aes(sample = resids)) +
          stat_qq() + stat_qq_line() +
          theme_minimal() +
          labs(title = "Normal Q-Q Plot")
      } else {
        ggplot(data.frame(fitted = fitted_vals, 
                         std_resids = sqrt(abs(resids))), 
               aes(x = fitted, y = std_resids)) +
          geom_point() +
          theme_minimal() +
          labs(x = "Fitted values", 
               y = expression(sqrt("|Standardized residuals|")),
               title = "Scale-Location Plot")
      }
    }
  })
  
  output$model_metrics <- renderTable({
    data.frame(
      Model = names(model_list),
      AICc = sapply(model_list, AIC, k = 2)
      # Dev_Explained = sapply(model_list, function(m) {
      #   sprintf("%.1f%%", 100 * m$deviance_explained)
      # })
    )
  })
  
  # output$dev_explained <- renderText({
  #   sprintf("%.1f%%", 100 * selected_model()$deviance_explained)
  # })
  
  output$aicc <- renderText({
    sprintf("%.1f", AIC(selected_model(), k = 2))
  })
}

# For standalone use
shinyApp(ui, server)

