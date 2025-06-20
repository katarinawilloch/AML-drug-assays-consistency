# A multi-center study on the consistency of drug sensitivity testing in patients with AML

This repository contains the code and analysis for the project:  
**"A multi-center study on the consistency of drug sensitivity testing in patients with AML"**  


---

## 🧰 Project Overview

This project performs statistical modeling and generates diagnostic plots to evaluate model assumptions and performance. It is written in R using packages like `lme4`, `ggplot2`, and `dplyr`.

### 🔍 Key Features
- Cleans and processes input data
- Perfomes some exploratory analysis
- Fits mixed-effects models using `lmer()`
- Produces diagnostic plots for residuals, random effects, and influential points

---
## 🔖 Data file Structure
   - Modelling.R contains code for creating the lmer() models and the fresh vs frozen comparisons. It also uses functions in model_funtions.R to create the diagnostic plots. 
   - Explore.R has code to create the heatmaps and ppca plots. 
   - The python files have code to clean the initial downloaded datasets and get the common names using the pubchem API. 
   - DSS_calculation.R has code to calculate DSS values from the raw downloaded dose responses. 
   - AUC_Andrea_way_calculation.R has code to calculare rAUC from the raw downloaded dose responses. 


---
## 🚀 How to Run

1. **Clone the repository**:
   ```bash
   git clone https://github.com/yourusername/your-repo.git
   cd your-repo

##

## 🖼️ Supplementary Diagnostic Plots

The following image shows an example of a diagnodtic plot for the Time until Usage LMM fitted on transformed DSS2 values. It helps assess normality, homoscedasticity, and influential observations. 
The rest of the supplementary diagnostic plots following the article can be found in [Supplementary Diagnostic Plots](https://github.com/katarinawilloch/AML-drug-assays-consistency/tree/main/Diagnostic%20plots)


![Residual Diagnostic Plot](https://github.com/katarinawilloch/AML-drug-assays-consistency/blob/main/Diagnostic%20plots/DSS2/time_until_sample_usage_final_plot.png)
