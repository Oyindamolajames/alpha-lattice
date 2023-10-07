# Readme for `lattice.R` R Function

## Overview

The `lattice_anv.R` script includes a function that allows you to conduct analysis of variance (ANOVA) on a collection of traits. Additionally, it utilizes the `lme4` package to compute Best Linear Unbiased Predictions (BLUPs) and access variance components.

This Readme serves as a guide for effectively utilizing this R function. Please ensure you have already loaded the essential libraries (`dplyr`, `readxl`, and `lme4`) before proceeding.

## Getting Started

To use the function, download the repository and arrange your code in a manner similar to the `script.r` file included in the repository. This will help ensure that you can easily incorporate and execute the `lattice_anv.R` script and the `analyse_lattice` function within your R project.

1. **Load Required Libraries**:

   Make sure you have the `dplyr` and `readxl` libraries installed and loaded in your R environment before using the `lattice_anv.R` script. You can load them using the following commands:

   ```R
   library(dplyr)
   library(readxl)
   ```

2. **Source the `lattice_anv.R` Script**:

   You need to source the `lattice_anv.R` script, which contains the `analyse_lattice` function. Make sure you provide the correct path to the script if it's not in your working directory:

   ```R
   source("path/to/lattice_anv.R")
   ```

## Data Preparation

Before using the `analyse_lattice` function, you need to load your data into R. You can load your data using the `read_excel` function from the `readxl` library or any other method that suits your data format.

```R
# Load your data here
data <- read_excel("path/to/your/data.xlsx")
```

### Checking Data

Always check the head of your data to ensure that it has been loaded correctly and that it contains the necessary columns for analysis.

```R
# Check data head
head(data)
```

## Using the `analyse_lattice` Function

The `analyse_lattice` function is used to perform the analysis of variance (ANOVA) on a set of traits. It takes the following parameters:

- `traits`: A subset of your data containing the columns for the traits you want to analyze.
- `replicates`: A vector or column from your data that represents the replicates.
- `blocks`: A vector or column from your data that represents the blocks.
- `genotypes`: A vector or column from your data that represents the genotypes.

Here's an example of how to use the `analyse_lattice` function:

```R
# Example of how to use the function
result <- analyse_lattice(data[4:20], data$Rep, data$Block, data$Genotype)
```

## Accessing the Output

The results of the analysis are stored in the `result` object. You can access the following components:

- `result$combined_results`: Combined results from the analysis.
- `result$variance`: Variance components.
- `result$mean`: Mean values.
- `result$blup`: Best linear unbiased predictions (BLUPs).

You can save these components to separate CSV files for further analysis or reporting. Here's an example of how to save them to a "result" folder:

```R
# Create a "result" folder if it doesn't exist
if (!dir.exists("result")){
  dir.create("result")
}

# Save the components to CSV files
write_csv(result$combined_results, "result/anvtable.csv")
write_csv(result$variance, "result/variance_components.csv")
write_csv(result$mean, "result/mean.csv")
write_csv(result$blup, "result/blups.csv")
```

## Reporting ANOVA Results

The ANOVA results can be printed to a text file for easier inspection. You can do this by using the `sink` function to redirect the output to a text file:

```R
# Redirect ANOVA results to a text file
sink("result/anova.txt")
result$anova_result
sink()
```

## Support and Issues

If you encounter any issues or have questions about using the `lattice_anv.R` script or the `analyse_lattice` function, please feel free to contact the author at oyindamolajames@gmail.com.

Please make sure to customize the script and usage instructions to fit the specific details of your data and analysis.
