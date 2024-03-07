library(shiny)
library(ggplot2)
library(dplyr)
library(readr)
library(ggpubr)
library(pracma)
library(shinydashboard) # Make sure to load the shinydashboard package

# Define UI
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .js-irs-0 {
        border-color: #5bc0de; /* Change slider border color */
      }
      .shiny-input-container {
        margin-bottom: 20px;
      }
      .box-title {
        font-size: 16px; /* Adjust title size */
        text-align: center; /* Center-align the box titles */
      }
      .box-subtitle {
        font-size: 12px; /* Adjust subtitle size */
        color: #777; /* Adjust subtitle color */
      }
    "))
  ),
  
  titlePanel("Dynamic Visualization of Expected Cost Differences"),
  
  fluidRow(
    column(width = 3,  # Adjust the width as necessary
           wellPanel(
             h3("Model Choice Panel", style = "margin-top: 0;"),
           selectInput("attribute_ratio", "Size of Minority Group?)",
                       choices = c("5% of total pop" = 0.05,
                                   "15% of total pop" = 0.15,
                                   "25% of total pop" = 0.25), width = '100%'),
           selectInput("p_pos_space", "How Many Positive Cases?",
                       choices = c("5% positives" = 0.05,
                                   "20% positives" = 0.20,
                                   "35% positives" = 0.35,
                                   "50% positives" = 0.50), width = '100%'),
            selectInput("ev_pos_a", "Minority Group Disease Presence",
                                   choices = c("63%" = 0.63, "72%" = 0.72, "80%" = 0.80), 
                                   selected = 0.80,
                                   width = '100%'),
            selectInput("ev_pos_b", "Majority Group Disease Presence",
                                   choices = c("63%" = 0.63, "72%" = 0.72, "80%" = 0.80), 
                                   selected = 0.80,
                                   width = '100%'),
            sliderInput("tau_max", "Max Threshold Considered", min = 0, max = 1, value = 0.5, step = 0.1,  width = '100%')
           
    )),
    
    column(width = 9,
           wellPanel(
             h3("Model Choice Distribution", style = "margin-top: 0;"),
             plotOutput("plot_freq_tot", height = "450px")
           )
    )
  ),
  
  # Control panel on the left
  fluidRow(
    column(3,
           wellPanel(
             h3("Relative Costs of FP's and FN's", style = "margin-top: 0;"),
             p("In medical testing, like cancer screenings, false positives (FP) and false negatives (FN) have significant implications. A false positive leads to unnecessary stress and medical procedures, while a false negative might result in delayed treatment for serious conditions."),
             p("Adjust the sliders to set the relative costs. Typically, FN's are costlier than FP's, especially in critical screenings where early detection is vital. These costs represent not just financial aspects but also the impact on patient outcomes and healthcare efficiency."),
             sliderInput("cost_fn", "Cost of False Negative's (FNs)", min = 0, max = 8, value = 4, step = 2, width = '100%'),
             sliderInput("cost_fp", "Cost of False Positive's (FPs)", min = 0, max = 8, value = 4, step = 2, width = '100%')
           )
    ),
    # First row of visualizations on the right
    column(9,
           wellPanel(
             h3("Expected Costs Over Different Thresholds"),
             h4("Area Under Expected Cost for AUROC in Red, and AUPRC in Blue", style = "color: #666;"),
             fluidRow(
               box(title = "Minority Group (A)", status = "primary", solidHeader = TRUE, plotOutput("exp_cost_plot_a"), width = 6),
               box(title = "Majority Group (B)", status = "primary", solidHeader = TRUE, plotOutput("exp_cost_plot_b"), width = 6)
             )
           )
    )
  ),
  
  # Second row of visualizations
  fluidRow(
    column(3,  # Text explaining the concepts in detail
           wellPanel(
           h3("Understanding AUEC", style = "margin-top: 0;"),
           h4("AUEC stands for Area Under the Expected Cost Curve, summarizing the model's cost efficiency over all thresholds."),
           p("The Expected Cost Curve plots the combined cost of false positives and negatives at each threshold, highlighting the trade-off between sensitivity and specificity."),
           p("Differences in AUEC between AUROC and AUPRC models indicate which is more cost-effective against a backdrop of varied false positive/negative costs."),
           p("A red heatmap shows where AUROC-optimized models have higher costs, while blue suggests AUPRC models are less cost-efficient in those cases."),
           p("The choice of color intensity on the heatmap underscores the magnitude of the cost difference, guiding decision-makers in model selection.")
    )),    
    column(9,
           wellPanel(
             h3("Difference in Area Under Expected Cost"),
             h4("Red Means Higher Values for AUROC and Blue Higher for AUPRC", style = "color: #666;"),
             h4("Threshold is chosen after Max Threshold Considered", style = "color: #666;"),
             fluidRow(
               box(title = "Minority Group (A)", status = "primary", solidHeader = TRUE, plotOutput("auec_plot_a"), width = 6),
               box(title = "Minority Group (B)", status = "primary", solidHeader = TRUE, plotOutput("auec_plot_b"), width = 6)
             )
           )
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  reactive_thresh_max <- reactive({
    suffix <- switch(as.character(input$tau_max),
                     "0.1" = "10",
                     "0.2" = "20",
                     "0.3" = "30",
                     "0.4" = "40",
                     "0.5" = "50",
                     "0.6" = "60",
                     "0.7" = "70",
                     "0.8" = "80",
                     "0.9" = "90",
                     "1" = "100")
    list(a = paste0("auec_mean_diff_", suffix, "_a"), 
         b = paste0("auec_mean_diff_", suffix, "_b"))
  })
  
  cost_thresh_data <- reactive({
    cost_thresh_data <- read_csv("shiny_data/thresholds_exp_costs_shiny_v2.csv")
    
    cost_thresh_data_tot <- cost_thresh_data %>%
      filter(
             attribute_ratio == input$attribute_ratio,
             p_pos == input$p_pos_space,
             ev_pos_a == input$ev_pos_a,
             ev_pos_b == input$ev_pos_b,
             cost_fn == as.numeric(input$cost_fn),
             cost_fp == as.numeric(input$cost_fp)) %>% 
      group_by(model_choice) %>% 
      slice(1)
    
    cost_thresh_data_a <- cost_thresh_data %>%
      filter(subgroup == "A",
             attribute_ratio == input$attribute_ratio,
             p_pos == input$p_pos_space,
             threshold <= input$tau_max,
             ev_pos_a == input$ev_pos_a,
             ev_pos_b == input$ev_pos_b,
             cost_fn == as.numeric(input$cost_fn),
             cost_fp == as.numeric(input$cost_fp))
    
    cost_thresh_data_b <- cost_thresh_data %>%
      filter(subgroup == "B",
             attribute_ratio == input$attribute_ratio,
             p_pos == input$p_pos_space,
             threshold <= input$tau_max,
             ev_pos_a == input$ev_pos_a,
             ev_pos_b == input$ev_pos_b,
             cost_fn == as.numeric(input$cost_fn),
             cost_fp == as.numeric(input$cost_fp))    
    
    return(list(cost_thresh_data_tot = cost_thresh_data_tot, cost_thresh_data_a = cost_thresh_data_a, cost_thresh_data_b = cost_thresh_data_b))
  })
  
  auec_diff_data <- reactive({
    auec_diffs <- read_csv("shiny_data/auec_data_shiny_v3.csv")
    
    auec_diffs <- auec_diffs %>% filter(attribute_ratio == input$attribute_ratio,
                                        p_pos_space == input$p_pos_space,
                                        ev_pos_a == input$ev_pos_a,
                                        ev_pos_b == input$ev_pos_b) %>% 
      mutate(cost_fn = round(cost_fn),
             cost_fp = round(cost_fp)
      )
    
    return(auec_diffs)
  })
  
  freq_data <- reactive({
    freq_data <- read_csv("shiny_data/simulated_probs_shiny.csv")
    
    freq_data <- freq_data %>% filter(attribute_ratio == input$attribute_ratio,
                                      p_pos_space == input$p_pos_space,
                                      ev_pos_a == input$ev_pos_a,
                                      ev_pos_b == input$ev_pos_b)
    return(freq_data)
  })
  
  output$plot_freq_tot <- renderPlot({
    req(freq_data())  # Make sure the reactive expression is ready
    
    freq_data_tot <- freq_data()
    
    freq_data_tot$positive = factor(ifelse(freq_data_tot$probability >= input$tau_max, "Positive", "Negative"))
    
    p_freq <- ggplot(freq_data_tot, aes(x = probability, fill = positive, color = positive)) + 
      geom_histogram(bins = 100, alpha = 0.9) +
      labs(x = "Predicted Risk", y = "Frequency",
           color = "", fill = "") +
      scale_color_manual(values = c("Positive" = "black", "Negative" = "black")) +
      scale_fill_manual(values = c("Positive" = "#9b2226", "Negative" = "#669bbc")) +
      scale_x_continuous(labels=scales::percent) +
      facet_grid(model_choice ~ subgroup, labeller = labeller(subgroup = c(a = "Minority (A)", b = "Minority (B)"))) +
      theme_minimal() +
      theme(
        legend.text = element_text(size = 14, color = "black", face = "bold"),
        axis.title.x = element_text(size = 14, color = "black", face = "bold"),
        axis.title.y = element_text(size = 14, color = "black", face = "bold"),
        axis.text.x = element_text(size = 14, color = "black", face = "bold"),
        axis.text.y = element_text(size = 14, color = "black", face = "bold"),
        legend.position = "bottom",
        legend.key.width = unit(1, "cm"),
        strip.text.x = element_text(size = 16, color = "black", face = "bold"),
        strip.text.y = element_text(size = 16, color = "black", face = "bold")
      )
    
    
    plot(p_freq)
  })
  

  
  output$exp_cost_plot_a <- renderPlot({
    req(cost_thresh_data())  # Make sure the reactive expression is ready
    
    cost_thresh_data <- cost_thresh_data()
    
    cost_thresh_data_a <- cost_thresh_data$cost_thresh_data_a
    
    tt_max_auroc_a = cost_thresh_data_a %>% filter(model_choice == "AUROC",
                                                   threshold <= input$tau_max)
    
    tt_max_auprc_a = cost_thresh_data_a %>% filter(model_choice == "AUPRC",
                                                   threshold <= input$tau_max)
    tt_max_auroc_a <- tt_max_auroc_a[order(tt_max_auroc_a$threshold), ]
    tt_max_auprc_a <- tt_max_auprc_a[order(tt_max_auprc_a$threshold), ]
    
    auec_auroc_a <- pracma::trapz(tt_max_auroc_a$threshold, tt_max_auroc_a$expected_cost)
    auec_auprc_a <- pracma::trapz(tt_max_auprc_a$threshold, tt_max_auprc_a$expected_cost)
    
    # Create the plot
    p_cost_thresh_data_a <- ggplot(cost_thresh_data_a, aes(x = threshold, y = expected_cost, color = model_choice)) +
      geom_area(aes(fill = model_choice), position = position_dodge(), alpha = 0.7) +  # Adjust alpha for visibility
      geom_line() +
      scale_color_manual(values = c("AUROC" = "black", "AUPRC" = "black"),
                         labels = c(paste0("AUROC (AUEC: ", round(auec_auroc_a, 2), ")"),
                                    paste0("AUPRC (AUEC: ", round(auec_auprc_a, 2), ")"))) +  # Assign colors for each model_choice
      scale_fill_manual(values = c("AUROC" = "#9b2226", "AUPRC" = "#669bbc"), 
                        labels = c(paste0("AUROC (AUEC: ", round(auec_auroc_a, 2), ")"),
                                   paste0("AUPRC (AUEC: ", round(auec_auprc_a, 2), ")"))) +   
      labs(x = "Threshold",
           y = "Expected Cost",
           color = "",
           fill = "") +
      theme_minimal() +
      theme(legend.text = element_text(size = 14, color = "black", face = "bold"),
            axis.title.x = element_text(size = 14, color = "black", face = "bold"),
            axis.title.y = element_text(size = 14, color = "black", face = "bold"),
            axis.text.x = element_text(size = 14, color = "black", face = "bold"),
            axis.text.y = element_text(size = 14, color = "black", face = "bold"),
            legend.position = "bottom",
            legend.key.width = unit(1.5, "cm"))
    plot(p_cost_thresh_data_a)
  })
  
  output$exp_cost_plot_b <- renderPlot({
    req(cost_thresh_data())  # Make sure the reactive expression is ready
    
    cost_thresh_data <- cost_thresh_data()
    cost_thresh_data_b <- cost_thresh_data$cost_thresh_data_b
    
    tt_max_auroc_b = cost_thresh_data_b %>% filter(model_choice == "AUROC",
                                                   threshold <= input$tau_max)
    
    tt_max_auprc_b = cost_thresh_data_b %>% filter(model_choice == "AUPRC",
                                                   threshold <= input$tau_max)
    tt_max_auroc_b <- tt_max_auroc_b[order(tt_max_auroc_b$threshold), ]
    tt_max_auprc_b <- tt_max_auprc_b[order(tt_max_auprc_b$threshold), ]
    
    auec_auroc_b <- pracma::trapz(tt_max_auroc_b$threshold, tt_max_auroc_b$expected_cost)
    auec_auprc_b <- pracma::trapz(tt_max_auprc_b$threshold, tt_max_auprc_b$expected_cost)
    
    # Create the plot
    p_cost_thresh_data_b<- ggplot(cost_thresh_data_b, aes(x = threshold, y = expected_cost, color = model_choice)) +
      geom_area(aes(fill = model_choice), position = position_dodge(), alpha = 0.7) +  # Adjust alpha for visibility
      geom_line() +
      scale_color_manual(values = c("AUROC" = "black", "AUPRC" = "black"),
                         labels = c(paste0("AUROC (AUEC: ", round(auec_auroc_b, 2), ")"),
                                    paste0("AUPRC (AUEC: ", round(auec_auprc_b, 2), ")"))) +  # Assign colors for each model_choice
      scale_fill_manual(values = c("AUROC" = "#9b2226", "AUPRC" = "#669bbc"),
                        labels = c(paste0("AUROC (AUEC: ", round(auec_auroc_b, 2), ")"),
                                   paste0("AUPRC (AUEC: ", round(auec_auprc_b, 2), ")"))) +  # Assign fill colors
      labs(x = "Threshold",
           y = "Expected Cost",
           color = "",
           fill = "") +
      theme_minimal() +
      # annotate("text", x = Inf, y = Inf, label = paste("AUEC AUROC:", round(auec_auroc_b, 2)), vjust = 3, hjust = 1.1) +
      # annotate("text", x = Inf, y = Inf, label = paste("AUEC AUPRC:", round(auec_auprc_b, 2)), vjust = 1, hjust = 1.1) +
      theme(legend.text = element_text(size = 14, color = "black", face = "bold"),
            axis.title.x = element_text(size = 14, color = "black", face = "bold"),
            axis.title.y = element_text(size = 14, color = "black", face = "bold"),
            axis.text.x = element_text(size = 14, color = "black", face = "bold"),
            axis.text.y = element_text(size = 14, color = "black", face = "bold"),
            legend.position = "bottom",
            legend.key.width = unit(1.5, "cm"))
    
    
    plot(p_cost_thresh_data_b)
  }) 
  
  output$auec_plot_a <- renderPlot({
    req(auec_diff_data())  # Make sure the reactive expression is ready
    auec_diff_data_a <- auec_diff_data()
    column_name_a <- reactive_thresh_max()$a
    
    # Ensure the dynamically selected column exists in the dataframe
    stopifnot(column_name_a %in% names(auec_diff_data_a))
    
    # Use the dynamic column name for the fill aesthetic
    p_auec_diff_data_a <- ggplot(auec_diff_data_a, aes(x = as.factor(cost_fn), y = as.factor(cost_fp), fill = .data[[column_name_a]])) +
      geom_tile() +
      scale_fill_gradient2(low = "#669bbc", high = "#9b2226", 
                           n.breaks = 4) +
      labs(x = "Relative Cost of False Negative", y = "Relative Cost of False Positive",
           fill = "") +
      theme_classic() +
      theme(legend.text = element_text(size = 14, color = "black", face = "bold"),
            axis.title.x = element_text(size = 14, color = "black", face = "bold"),
            axis.title.y = element_text(size = 14, color = "black", face = "bold"),
            axis.text.x = element_text(size = 14, color = "black", face = "bold"),
            axis.text.y = element_text(size = 14, color = "black", face = "bold"),
            legend.position = "bottom",
            legend.key.width = unit(1.5, "cm"))
    
    plot(p_auec_diff_data_a)
  })
  
  output$auec_plot_b <- renderPlot({
    req(auec_diff_data())  # Make sure the reactive expression is ready
    auec_diff_data_b <- auec_diff_data()
    column_name_b <- reactive_thresh_max()$b
    
    # Ensure the dynamically selected column exists in the dataframe
    stopifnot(column_name_b %in% names(auec_diff_data_b))
    # Get the range of values for the fill aesthetic
    value_range <- range(auec_diff_data_b[[column_name_b]], na.rm = TRUE)
    
    # Create three breaks, including min and max
    breaks <- pretty(value_range, n = 3)
    
    p_auec_diff_data_b <- ggplot(auec_diff_data_b, aes(x = as.factor(cost_fn), y = as.factor(cost_fp), fill = .data[[column_name_b]])) +
      geom_tile() +
      scale_fill_gradient2(low = "#669bbc", high = "#9b2226",
                           n.breaks = 4) +
      labs(x = "Relative Cost of False Negative", y = "Relative Cost of False Positive",
           fill = "") +
      theme_classic() +
      theme(legend.text = element_text(size = 14, color = "black", face = "bold"),
            axis.title.x = element_text(size = 14, color = "black", face = "bold"),
            axis.title.y = element_text(size = 14, color = "black", face = "bold"),
            axis.text.x = element_text(size = 14, color = "black", face = "bold"),
            axis.text.y = element_text(size = 14, color = "black", face = "bold"),
            legend.position = "bottom",
            legend.key.width = unit(1.5, "cm"))
    
    plot(p_auec_diff_data_b)
  })
  
}

shinyApp(ui = ui, server = server)