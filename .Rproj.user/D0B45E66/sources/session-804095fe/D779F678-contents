library(dplyr)
library(readxl)
source("lattice_anv.R")

## Load your data here
data <-

# check data head
head(data)

# how to use the function
# the 4:20 is the coulumn number of the traits

result <- analyse_lattice(data[4:20], data$Rep, data$Block, data$Genotype)

# Access the output from result folder

if (!dir.exists("result")){
  dir.create("result")
}

write_csv(result$combined_results, "result/anvtable.csv")
write_csv(result$variance, "result/variance_components.csv")
write_csv(result$mean, "result/mean.csv")
write_csv(result$blup, "result/blups.csv")

sink("result/anova.txt")
result$anova_result
sink()

# If you have any issues, send a mail to oyindamolajames@gmail.com
