#' Function to analyse alpha lattice design
#'
#' @param traits
#' @param rep
#' @param blk
#' @param entry
#'
#' @return
#' @export
#'
#' @examples
analyse_lattice <- function(traits, rep, blk, entry) {
library(lme4)
library(tidyverse)
    # rep <- data_2021$rep
    # blk <- data_2021$inc_block
    # entry <- data_2021$entry
    # traits <<- data_2021[8:10]

  rep <- as.factor(rep)
  entry <- as.factor(entry)
  blk <- as.factor(blk)

  # Check to ensure data is of right data type
  if (!all(sapply(traits, is.double))) {
    print("traits should be of type numeric")
    return(NULL)
  }

  data_mod <- data.frame(rep = rep,
                         genotype = entry,
                         blk = blk)
  data_set <-  cbind(data_mod, traits)

  # Create an empty list to store the model results
  model_results <- list()
  model_anv <- list()

  # get trait columns
  trait_columns <- colnames(traits)

  # Loop through each trait column and run the model
  for (trait in trait_columns) {
    # Create the formula for the model
    formula <- as.formula(paste(trait,
                                "~ rep + blk%in%rep + genotype"))

    # Fit the model
    fix_mod <- aov(formula, data = data_set)

    # Perform ANOVA
    anv <- anova(fix_mod)


    # Extract mean squares column and level of significance
    mean_squares <- round(anv[, "Mean Sq"], 2)

    significance <- anv %>%
      mutate(significance = case_when(
        `Pr(>F)` < 0.01 ~ "**", # If p-value is less than 0.01, assign "**"
        `Pr(>F)` < 0.05 ~ "*", # If p-value is less than 0.05 but
        #greater than or equal to 0.01, assign "*"
        TRUE ~ ""  # For all other cases (p-values >= 0.05), assign an empty string ""
      ))
    significance <- significance$significance


    # Concatenate mean square and significance into a single column
    ms_sig <- paste(mean_squares, significance)

    # Create a new dataframe for the ANOVA results
    result_df <- data.frame(Trait = trait, Mean_Square = ms_sig)

    # Store the dataframe in the list
    model_results[[trait]] <- result_df

    model_anv[[trait]] <- anv
  }



  # Combine all the ANOVA results into a single dataframe
  combined_results <- do.call(cbind,
                              lapply(model_results,
                                     function(x) x[, 2]))
  # Add the df and sov column to the anova
  names <- rownames(model_anv[[trait]])
  df <- model_anv[[trait]][1]$Df
  names[1] <- 'rep'
  names[2] <- 'Gen'
  names[3] <- "blk (rep)"
  names[4] <- 'error'

  # convert combined_results to a dataframe
  combined_results <- as.data.frame(combined_results)

  combined_results <- combined_results |>
    mutate(sov = names, df = df) |>
    select(sov, df, everything())

  # Print the combined results
  knitr::kable(combined_results)


  # Print the results for each trait
  for (trait in trait_columns) {
    print(paste("Results for", trait))
    print(model_anv[[trait]])
    print("\n")
  }

  ## Means
  mean_1 <- data_set |>
    group_by(genotype) |>
    summarize(across(3:(ncol(data_set) - 1), mean))


# Variance components -----------------------------------------------------

  # Create data output to contain blup
  data <- data_set
  DataOutput <- data.frame(genotype = unique(data_set$genotype))

  # Get column numbers of columns with data
  colnum <- 4:ncol(data_set)

  # Get the column names of the dataframe (excluding non-numeric columns if needed)
  trait_columns <- colnames(data_set)[sapply(data_set, is.numeric)]

  # Create an empty list to store the model results
  model_results <- list()

  # Variance Component
  # Create Dataframe for variance components
  DataVarComp <- data.frame()
  DataVarCompOutput <- data.frame()
  HeritabilityData <- data.frame()

  #this empty dataframe is for dropped variance components
  drops <- c("var1","var2","sdcor")

  for (trait in trait_columns) {
    # Create the formula for the model
    formula <- as.formula(paste(trait,
                                "~ (1|blk%in%rep) + (1|genotype) + (1|rep)"))

    # Fit the model
    mix_mod <- lmer(formula, data = data_set)
    summary(mix_mod)
    blup = coef(mix_mod)$genotype
    colnames(blup) <- trait
    blup$genotype <- rownames(blup)
    # Bind the column to the output
    DataOutput <- merge(DataOutput,blup, by="genotype", all=T)

    varComp<-as.data.frame(VarCorr(mix_mod,comp="vcov"))

    varComp <- varComp[ , !(names(varComp) %in% drops)]
    varComp$Trait <- trait
    DataVarComp <- rbind(DataVarComp,varComp)


    model_results[[trait]] <- mix_mod
  }

  DataVarCompOutput <- reshape(DataVarComp,
                               idvar = "Trait",
                               timevar = "grp",
                               direction = "wide")
  #reshape our variance components dataframe so that we can run the heritability script
  HeritabilityData <- ( (DataVarCompOutput[,2]) / ( (DataVarCompOutput[,2]) + (DataVarCompOutput[,5]) )) *100

  #create function to identify maximum value in a column
  colMax <- function(data) sapply(data, max, na.rm = TRUE)
  #create function to identify minimum value in a column
  colMin <- function(data) sapply(data, min, na.rm = TRUE)
  #create function to identify mean value in a column
  colMean <- function(data) sapply(data, mean, na.rm = TRUE)
  #create function to identify median value in a column
  colMedian <- function(data) sapply(data, median, na.rm = TRUE)
  #create function to identify standard dev value in a column
  colStdev <- function(data) sapply(data, sd, na.rm = TRUE)


  df <- data
  #summary statistics
  DataColMax <- colMax(df[,colnum])
  DataColMin <- colMin(df[,colnum])
  DataColMean <- colMean(df[,colnum])
  DataColMedian <- colMedian(df[,colnum])
  DataColStdev <- colStdev(df[,colnum])
  GCV <- (sqrt(DataVarCompOutput[,2])) / ((DataColMean)) * 100
  PCV <- (sqrt(DataVarCompOutput[,2]) + sqrt(DataVarCompOutput[,5])) / ((DataColMean)) * 100
  CV <- DataColStdev / DataColMean * 100

  DataVarCompOutput <- cbind(DataVarCompOutput, HeritabilityData,
                             GCV, PCV, CV, DataColMin, DataColMax, DataColMean,
                             DataColMedian, DataColStdev)


  result = list()

  result$variance <-  DataVarCompOutput
  result$blup <- DataOutput

  result$combined_results = combined_results
  result$anova_result = model_anv
  result$mean = mean_1

  return(result)
}

analyse_lattice_comb <- function(traits, rep, blk, env, entry) {
  library(lme4)
  library(tidyverse)
  # rep <- data_combined$rep
  # blk <- data_combined$inc_block
  # entry <- data_combined$accession
  # env <- data_combined$year

  #traits <<- data_combined[8:10]

  rep <- as.factor(rep)
  entry <- as.factor(entry)
  blk <- as.factor(blk)
  env <- as.factor(env)

  # Check to ensure data is of right data type
  if (!all(sapply(traits, is.double))) {
    print("traits should be of type numeric")
    return(NULL)
  }

  data_mod <- data.frame(rep = rep,
                         genotype = entry,
                         blk = blk,
                         env = env)
  data_set <-  cbind(data_mod, traits)

  # Create an empty list to store the model results
  model_results <- list()
  model_anv <- list()

  # get trait columns
  trait_columns <- colnames(traits)

  # Loop through each trait column and run the model
  for (trait in trait_columns) {
    # Create the formula for the model
    formula <- as.formula(paste(trait,
                                "~ genotype*env + rep%in%env +  blk%in%rep%in%env"
                                ))

    # Fit the model
    fix_mod <- aov(formula, data = data_set)

    # Perform ANOVA
    anv <- anova(fix_mod)


    # Extract mean squares column and level of significance
    mean_squares <- round(anv[, "Mean Sq"], 2)

    significance <- anv %>%
      mutate(significance = case_when(
        `Pr(>F)` < 0.01 ~ "**", # If p-value is less than 0.01, assign "**"
        `Pr(>F)` < 0.05 ~ "*", # If p-value is less than 0.05 but
        #greater than or equal to 0.01, assign "*"
        TRUE ~ ""  # For all other cases (p-values >= 0.05), assign an empty string ""
      ))
    significance <- significance$significance


    # Concatenate mean square and significance into a single column
    ms_sig <- paste(mean_squares, significance)

    # Create a new dataframe for the ANOVA results
    result_df <- data.frame(Trait = trait, Mean_Square = ms_sig)

    # Store the dataframe in the list
    model_results[[trait]] <- result_df

    model_anv[[trait]] <- anv
  }



  # Combine all the ANOVA results into a single dataframe
  combined_results <- do.call(cbind,
                              lapply(model_results,
                                     function(x) x[, 2]))
  # Add the df and sov column to the anova
  names <- rownames(model_anv[[trait]])
  df <- model_anv[[trait]][1]$Df
  # names[1] <- 'rep'
  # names[2] <- 'Gen'
  # names[3] <- "blk (rep)"
  # names[4] <- 'error'

  # convert combined_results to a dataframe
  combined_results <- as.data.frame(combined_results)

  combined_results <- combined_results |>
    mutate(sov = names, df = df) |>
    select(sov, df, everything())

  # Print the combined results
  knitr::kable(combined_results)


  # Print the results for each trait
  for (trait in trait_columns) {
    print(paste("Results for", trait))
    print(model_anv[[trait]])
    print("\n")
  }

  ## Means
  mean_1 <- data_set |>
    group_by(genotype) |>
    summarize(across(4:(ncol(data_set) - 1), mean))


  # Variance components -----------------------------------------------------

  # Create data output to contain blup
  data <- data_set
  DataOutput <- data.frame(genotype = unique(data_set$genotype))

  # Get column numbers of columns with data
  colnum <- 5:ncol(data_set)

  # Get the column names of the dataframe (excluding non-numeric columns if needed)
  trait_columns <- colnames(data_set)[sapply(data_set, is.numeric)]

  # Create an empty list to store the model results
  model_results <- list()

  # Variance Component
  # Create Dataframe for variance components
  DataVarComp <- data.frame()
  DataVarCompOutput <- data.frame()
  HeritabilityData <- data.frame()

  #this empty dataframe is for dropped variance components
  drops <- c("var1","var2","sdcor")

  for (trait in trait_columns) {
    # Create the formula for the model
    formula <- as.formula(paste(trait,
                                "~ (1|env) + (1|blk%in%rep) + (1|genotype) + (1|genotype:env)"
                                ))

    # Fit the model
    mix_mod <- lmer(formula, data = data_set)
    summary(mix_mod)
    blup = coef(mix_mod)$genotype
    colnames(blup) <- trait
    blup$genotype <- rownames(blup)
    # Bind the column to the output
    DataOutput <- merge(DataOutput,blup, by="genotype", all=T)

    varComp<-as.data.frame(VarCorr(mix_mod,comp="vcov"))

    varComp <- varComp[ , !(names(varComp) %in% drops)]
    varComp$Trait <- trait
    DataVarComp <- rbind(DataVarComp,varComp)


    model_results[[trait]] <- mix_mod
  }

  DataVarCompOutput <- reshape(DataVarComp,
                               idvar = "Trait",
                               timevar = "grp",
                               direction = "wide")
  #reshape our variance components dataframe so that we can run the heritability script
  # Broad sense heritability script
  # Taken from Gioia et al. 2017
  nloc = 2
  HeritabilityData <- ((DataVarCompOutput[,3])) / (((DataVarCompOutput[,3])) + (((DataVarCompOutput[,6])) / (nloc))) *100
  #create function to identify maximum value in a column
  colMax <- function(data) sapply(data, max, na.rm = TRUE)
  #create function to identify minimum value in a column
  colMin <- function(data) sapply(data, min, na.rm = TRUE)
  #create function to identify mean value in a column
  colMean <- function(data) sapply(data, mean, na.rm = TRUE)
  #create function to identify median value in a column
  colMedian <- function(data) sapply(data, median, na.rm = TRUE)
  #create function to identify standard dev value in a column
  colStdev <- function(data) sapply(data, sd, na.rm = TRUE)


  df <- data
  #summary statistics
  DataColMax <- colMax(df[,colnum])
  DataColMin <- colMin(df[,colnum])
  DataColMean <- colMean(df[,colnum])
  DataColMedian <- colMedian(df[,colnum])
  DataColStdev <- colStdev(df[,colnum])
  GCV <- (sqrt(DataVarCompOutput[,3])) / ((DataColMean)) * 100
  PCV <- (sqrt(DataVarCompOutput[,3]) + sqrt(DataVarCompOutput[,6])) / ((DataColMean)) * 100
  CV <- DataColStdev / DataColMean * 100

  DataVarCompOutput <- cbind(DataVarCompOutput, HeritabilityData,
                             GCV, PCV, CV, DataColMin, DataColMax, DataColMean,
                             DataColMedian, DataColStdev)


  result = list()

  result$variance <-  DataVarCompOutput
  result$blup <- DataOutput

  result$combined_results = combined_results
  result$anova_result = model_anv
  result$mean = mean_1

  return(result)
}

