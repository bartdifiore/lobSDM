# Scaling/unscaling function to facilitate model convergence
set.seed(13)
x <- rnorm(100)

scaled_x <- scale(x)
unscaled_x <- as.numeric((scaled_x * attr(scaled_x, "scaled:scale")) + attr(scaled_x, "scaled:center"))
all.equal(x, unscaled_x)

scaled<- function(x, center, scale){
  (x - center) / scale
}

unscale <- function(scaled_x, center, scale) {
  if (is.null(attr(scaled_x, "scaled:scale")) == F) {
    # (scaled_x * sd) + m
    (scaled_x * attr(scaled_x, "scaled:scale")) + attr(scaled_x, "scaled:center")
  }
  if (is.null(attr(scaled_x, "scaled:scale")) == T) {
    (scaled_x * scale) + center
  }
}

unscale_aja <- function(scaled_x, orig_mean, orig_sd) {
  (scaled_x * orig_sd) + orig_mean
}

# Function to plot smooth effects
plot_smooths<- function(mod_fit, n_preds = 100, y_lab = "Predicted Biomass", rescale_means = NULL, rescale_sds = NULL){
  # Get smooth terms from the model
  formula_terms <- as.character(mod_fit$smoothers$labels)  # Right-hand side of formula
  
  # Clean up to get just the variable names
  smooth_terms <- gsub("\\)", "", gsub("s\\(", "", formula_terms))
  
  # Create prediction data frame for all smooth terms
  all_preds <- lapply(smooth_terms, function(term) {
    # Generate term string
    term_string <- paste0(term, paste0(" [n=", n_preds, "]"))
    
    # Pass to ggpredict
    pred <- ggpredict(mod_fit, terms = term_string)
    
    # Unscale for plotting?
    pred$smooth_term <- term
    if(!is.null(rescale_means) & !is.null(rescale_sds)){
      pred$x_raw <- unscale_aja(pred$x, rescale_means[[gsub("_scaled", "", term)]], rescale_sds[[gsub("_scaled", "", term)]])
    }
    pred$smooth_term <- term
    return(pred)
  })
  
  # Combine all predictions
  pred_data <- data.frame(bind_rows(all_preds))
  
  # Create the plot
  p <- ggplot() +
    geom_ribbon(data = pred_data, aes(x = x_raw, ymin = conf.low, ymax = conf.high, fill = smooth_term), alpha = 0.1, color = NA) +
    geom_line(data = pred_data, aes(x = x_raw, y = predicted, color = smooth_term), linewidth = 1) +
    labs(
      x = "Predictor value",
      y = y_lab
    ) +
    theme(
      text = element_text(size = 14),
      legend.position = "bottom"
    ) +
    facet_wrap(~ gsub("_scaled", "", smooth_term), scales = "free_x") +
    theme_bw()
  
  return(p)
}