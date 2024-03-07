## Visualizations Overview

Shiny app link:

[https://fnh30i-lasse-hansen.shinyapps.io/AUCIAYN_shiny_v2/](https://fnh30i-lasse-hansen.shinyapps.io/AUCIAYN_shiny_v2/)

The app integrates various visualizations to aid in the understanding and analysis of expected costs related to medical testing. Each visualization is designed to offer insights into different aspects of the cost-benefit analysis of diagnostic tests.

### Model Choice Distribution

This visualization showcases the distribution of predicted risk scores, segmented by the chosen model (e.g., AUROC vs. AUPRC) and subgroup (e.g., Minority Group A vs. Majority Group B). It provides a histogram that helps in understanding how the risk scores are spread across different probabilities, highlighting the frequency of scores that lead to a positive or negative classification.

#### Interpretation:

- **Histograms**: Observe the distribution of risk scores. A well-calibrated model should ideally show a clear distinction between the low-risk (closer to 0) and high-risk (closer to 1) predictions.
- **Color Coding**: Different colors represent the model choice, aiding in quick visual comparison between the outcomes of different modeling approaches.

### Expected Costs Over Different Thresholds

This section presents a comparative view of the expected costs for two major groups (A and B), usually representing different population segments or conditions. The costs are plotted over a range of thresholds, providing insights into how sensitive the cost outcomes are to changes in the decision threshold.

#### Interpretation:

- **Area Plots**: The shaded area under the curve illustrates the expected cost at each threshold, helping to identify the threshold that minimizes the cost.
- **AUROC vs. AUPRC**: Comparing these plots can highlight which model performs better in terms of cost-efficiency across different thresholds.

### Difference in Area Under Expected Cost

These visualizations delve into the cost differences between models optimized for AUROC and AUPRC. They utilize heatmaps to represent the cost efficiency of each model at various combinations of false-positive and false-negative costs.

#### Interpretation:

- **Heatmaps**: Each cell represents the cost difference for a specific combination of FP and FN costs. Colors indicate whether one model (AUROC or AUPRC) demonstrates a higher expected cost than the other.
- **Color Intensity**: The intensity of the color corresponds to the magnitude of the cost difference, with red indicating areas where the AUROC model is costlier and blue where the AUPRC model is less efficient.

Feel free to explore these interactive visualizations to gain comprehensive insights into the expected costs associated with different diagnostic thresholds and model choices. Adjust the sliders and inputs to see how varying parameters affect the visualized outcomes, providing a deeper understanding of the underlying cost dynamics in medical decision-making.
