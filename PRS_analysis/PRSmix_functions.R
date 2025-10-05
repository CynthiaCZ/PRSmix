# function to call PRSmix
call_PRSmix <- function(df, outcome_name, outcome_binary=FALSE, covar_list, prs_list, coh, cat_covar_list) {
  in_dir <- "path/to/data/PRSmix_in_files"
  out_dir <- "path/to/data/PRSmix_out_files"
  
  # get phenotype
  pheno <- df[, c("sid", outcome_name)]
  colnames(pheno)[1] <- "IID"
  fwrite(pheno, file = paste0(in_dir, "/pheno_", coh, "_", outcome_name, ".txt"), sep = "\t")
  
  # get covariates
  covariate <- df[, c("sid", covar_list)]
  colnames(covariate)[1] <- "IID"
  covariate <- na.omit(covariate)
  fwrite(covariate, file = paste0(in_dir, "/covariate_", coh, ".txt"), sep = "\t")
  
  # get PRS
  pgs <- df[, c("sid", prs_list)]
  names(pgs)[-1] <- paste0(names(pgs)[-1], "_SUM")
  names(pgs)[1] <- "IID"
  fwrite(pgs, file = paste0(in_dir, "/pgs_", coh, ".txt"), sep = "\t")
  
  # write list of PRSs
  writeLines(prs_list, paste0(in_dir, "/trait_specific_score.txt"))
  
  # PRSmix
  combine_PRS(
    pheno_file = paste0(in_dir, "/pheno_", coh, "_", outcome_name, ".txt"),
    covariate_file = paste0(in_dir, "/covariate_", coh, ".txt"),
    score_files_list = paste0(in_dir, "/pgs_", coh, ".txt"),
    trait_specific_score_file = paste0(in_dir, "/trait_specific_score.txt"),
    pheno_name = outcome_name,
    isbinary = outcome_binary,
    out = paste0(out_dir, "/", coh, "_", outcome_name),
    IID_pheno = "IID",
    covar_list = covar_list,
    cat_covar_list = cat_covar_list
  )
}


# function to scale coefficients with results from PRSmix
scale_PRSmix <- function(df, outcome_name, coh) {
  out_dir <- "path/to/data/PRSmix_out_files"
  weight_file <- paste0(out_dir, "/", coh, "_", outcome_name, "_power.0.95_pthres.0.05_weight_PRSmixPlus.txt")
  df[[paste0("mix_prs_", outcome_name)]] <- 0
  
  # df to store what prss are selected and the weights
  weights_table <- data.frame(prs = character(), weight = character(), stringsAsFactors = FALSE)
  
  # Check if the weight file exists
  if (file.exists(weight_file)) {
    # File exists, read the weights and process as usual
    weights <- read.table(weight_file, header = TRUE)
    
    for (i in 1:nrow(weights)) {
      prs_name <- weights$topprs[i]
      weight <- weights$ww[i]
      df[[paste0("mix_prs_", outcome_name)]] <- df[[paste0("mix_prs_", outcome_name)]] + df[[prs_name]] * weight
      
      weights_table <- rbind(weights_table, data.frame(prs = prs_name, weight = weight))
    }
    
    # if there is no weight file
  } else {
    best_acc <- read.table(paste0(out_dir, "/", coh, "_", outcome_name, "_best_acc.txt"), header = TRUE, sep = "\t", na.strings = "")
    best_prs <- best_acc$pgs
    df[[paste0("mix_prs_", outcome_name)]] <- df[[best_prs]]
    # weight is 1 when only 1 prs is selected
    weights_table <- data.frame(prs = best_prs, weight = 1)
  }
  
  df[[paste0("mix_prs_", outcome_name)]] <- scale(df[[paste0("mix_prs_", outcome_name)]])
  names(weights_table)[2] <- outcome_name
  
  return(list(df = df, weights_table = weights_table))
}

# function to build linear models
build_model <- function(df, outcome_name, prs=NULL, pcs=NULL, covars=NULL, model_type = "lm", adjustment=NULL) {
  
  # NULL values will be ignored
  formula_components <- c(prs, pcs, covars, adjustment)
  formula_str <- paste0(outcome_name, " ~ ", paste0(formula_components, collapse = " + "))
  
  # build model from formula
  if (model_type == "lm") {
    model <- lm(as.formula(formula_str), data = df)
  } else if (model_type == "nb") {
    model <- glm.nb(as.formula(formula_str), data = df)
  } else if (model_type == "binomial") {
    model <- glm(as.formula(formula_str), data = df, family = binomial())
  } else if (model_type == "poisson") {
    model <- glm(as.formula(formula_str), data = df, family = poisson())
  } else {
    stop("Unsupported model type. Use 'lm', 'nb', 'binomial', or 'poisson'")
  }
  
  return(model)
}

# function to call linear models
call_model <- function(df, outcome_name, pcs, covars, model_type = "lm", adjustment = NULL, coh) {

  out_dir <- "path/to/data/PRSmix_out_files"
  
  weight_file <- paste0(out_dir, "/", coh, "_", outcome_name, "_power.0.95_pthres.0.05_weight_PRSmixPlus.txt")
  
  results <- data.frame(prs = character(), beta = numeric(), ci_lower = numeric(), ci_upper = numeric(), 
                        beta_ci = character(), se = numeric(), stringsAsFactors = FALSE)
  
  # Check if the weight file exists
  if (file.exists(weight_file)) {
    weights <- read.table(weight_file, header = TRUE)
    # get prss that were selected, and add mix_prs
    prss <- weights$topprs
    if (length(prss) > 1) {
      prss <- c(prss, paste0("mix_prs_", outcome_name))
    }

    for (prs in prss) {
      model <- build_model(df, outcome_name, prs, pcs, covars, model_type, adjustment)
      beta <- coef(model)[prs]
      ci <- confint(model, prs, level = 0.95)
      beta_ci <- sprintf("%.3f (%.3f, %.3f)", beta, ci[1], ci[2])
      se <- summary(model)$coefficients[prs, "Std. Error"]
      pval <- summary(model)$coefficients[prs, 4]
      
      prs <- ifelse(grepl("^mix_prs_", prs), "mix_prs", prs) #rename
      results <- rbind(results, data.frame(prs = prs, beta = beta, ci_lower = ci[1], ci_upper = ci[2], 
                                           beta_ci = beta_ci, se = se, pval = pval))
    }
    
    # if there is no weight file, use the best one
  } else {
    best_acc <- read.table(paste0(out_dir, "/", coh, "_", outcome_name, "_best_acc.txt"), header = TRUE, sep = "\t", na.strings = "")
    prs <- best_acc$pgs
    model <- build_model(df, outcome_name, prs, pcs, covars, model_type, adjustment)
    beta <- coef(model)[prs]
    ci <- confint(model, prs, level = 0.95)
    beta_ci <- sprintf("%.3f (%.3f, %.3f)", beta, ci[1], ci[2])
    se <- summary(model)$coefficients[prs, "Std. Error"]
    pval <- summary(model)$coefficients[prs, 4]
    
    results <- data.frame(prs = prs, beta = beta, ci_lower = ci[1], ci_upper = ci[2], 
                          beta_ci = beta_ci, se = se, pval = pval)
    
  }
  colnames(results)[colnames(results) == "beta_ci"] <- outcome_name
  return(results)
}

# function to compare ROC with DeLong p-value, test for model miscalibration using Hosmer-Lemeshow tests, and plot
compare_roc <- function(df, outcome_binary, prs, pcs, covars, adjustment) {
  # clinical covariates
  model1 <- build_model(df, outcome_binary, prs=NULL, pcs=NULL, covars, model_type='binomial', adjustment)
  pred_prob1 <- predict(model1, type = "response")
  roc_covar <- roc(df[[outcome_binary]], pred_prob1)
  
  # PRS-mix + PCs
  model2 <- build_model(df, outcome_binary, prs, pcs, covars=NULL, model_type='binomial', adjustment)
  pred_prob2 <- predict(model2, type = "response")
  roc_PRSmix <- roc(df[[outcome_binary]], pred_prob2)
  
  # full model of clinical covariates + PRSmix + PCs
  model3 <- build_model(df, outcome_binary, prs, pcs, covars, model_type='binomial', adjustment)
  pred_prob3 <- predict(model3, type = "response")
  roc_full <- roc(df[[outcome_binary]], pred_prob3)
  
  # compare with DeLong p-value
  print(roc.test(roc_covar, roc_PRSmix))
  print(roc.test(roc_covar, roc_full))
  print(roc.test(roc_PRSmix, roc_full))
  
  # test for model miscalibration using Hosmer-Lemeshow tests
  hoslem.test1 <- hoslem.test(df[[outcome_binary]], pred_prob1)
  hoslem.test2 <- hoslem.test(df[[outcome_binary]], pred_prob2)
  hoslem.test3 <- hoslem.test(df[[outcome_binary]], pred_prob3)
  
  print(hoslem.test1)
  print(hoslem.test2)
  print(hoslem.test3)
  
  # plot curves
  plot(roc_covar, col = "blue", main = paste(outcome_binary, "predicted by", prs), cex.main = 0.8)
  plot(roc_PRSmix, col = "red", add = TRUE)
  plot(roc_full, col = "green", add = TRUE)
  
  legend("bottomright", legend = c("Clinical Covariates", "PRS-mix + PCs", "Full Model"), 
         col = c("blue", "red", "green"), lwd = 2)
  
  # return auc value of PRSmix model
  return(as.numeric(auc(roc_PRSmix)))
}

compare_roc_5_models <- function(df, outcome_binary, mix_prs, prs_to_compare, pcs, covars, adjustment=NULL, plot_title) {
  # clinical covariates
  model1 <- build_model(df, outcome_binary, prs=NULL, pcs=NULL, covars, model_type='binomial', adjustment)
  pred_prob1 <- predict(model1, type = "response")
  roc_covar <- roc(df[[outcome_binary]], pred_prob1)
  
  # PRSmix + PCs
  model2 <- build_model(df, outcome_binary, prs=mix_prs, pcs, covars=NULL, model_type='binomial', adjustment)
  pred_prob2 <- predict(model2, type = "response")
  roc_PRSmix <- roc(df[[outcome_binary]], pred_prob2)
  
  # clinical covariates + PRSmix + PCs
  model3 <- build_model(df, outcome_binary, prs=mix_prs, pcs, covars, model_type='binomial', adjustment)
  pred_prob3 <- predict(model3, type = "response")
  roc_covar_PRSmix <- roc(df[[outcome_binary]], pred_prob3)
  
  # PRS to compare + PCs
  model4 <- build_model(df, outcome_binary, prs=prs_to_compare, pcs, covars=NULL, model_type='binomial', adjustment)
  pred_prob4 <- predict(model4, type = "response")
  roc_prs_to_compare <- roc(df[[outcome_binary]], pred_prob4)
  
  # clinical covariates + PRS to compare + PCs
  model5 <- build_model(df, outcome_binary, prs=prs_to_compare, pcs, covars, model_type='binomial', adjustment)
  pred_prob5 <- predict(model5, type = "response")
  roc_covar_prs_to_compare <- roc(df[[outcome_binary]], pred_prob5)
  
  auc_table <- data.frame(
    Model = c(
      "Clinical Covariates",
      "PRS-mix + PCs",
      "Clinical + PRS-mix + PCs",
      "FEV1/FVC PRS + PCs",
      "Clinical + FEV1/FVC PRS + PCs"
    ),
    AUC = c(
      signif(auc(roc_covar), 3),
      signif(auc(roc_PRSmix), 3),
      signif(auc(roc_covar_PRSmix), 3),
      signif(auc(roc_prs_to_compare), 3),
      signif(auc(roc_covar_prs_to_compare), 3)
    )
  )
  
  # test for model miscalibration using Hosmer-Lemeshow tests
  cat("clinical covariates")
  print(hoslem.test(df[[outcome_binary]], pred_prob1))
  cat("PRSmix + PCs")
  print(hoslem.test(df[[outcome_binary]], pred_prob2))
  cat("clinical covariates + PRSmix + PCs")
  print(hoslem.test(df[[outcome_binary]], pred_prob3))
  cat("PRS to compare + PCs")
  print(hoslem.test(df[[outcome_binary]], pred_prob4))
  cat("clinical covariates + PRS to compare + PCs")
  print(hoslem.test(df[[outcome_binary]], pred_prob5))
  
  
  # compare with DeLong p-value
  print(roc.test(roc_PRSmix, roc_prs_to_compare))
  print(roc.test(roc_covar_PRSmix, roc_covar_prs_to_compare))
  print(roc.test(roc_covar_PRSmix, roc_covar))
  print(roc.test(roc_covar_prs_to_compare, roc_covar))
  
  pval1 <- roc.test(roc_PRSmix, roc_prs_to_compare)$p.value
  pval2 <- roc.test(roc_covar_PRSmix, roc_covar_prs_to_compare)$p.value
  pval3 <- roc.test(roc_covar_PRSmix, roc_covar)$p.value
  pval4 <- roc.test(roc_covar_prs_to_compare, roc_covar)$p.value
  
  pval_table <- data.frame(
    comparison = c("PRS-mix + PCs vs FEV1/FVC PRS + PCs", 
                   "Clinical + PRS-mix + PCs vs Clinical + FEV1/FVC PRS + PCs",
                   "Clinical + PRS-mix + PCs vs Clinical",
                   "Clinical + FEV1/FVC PRS + PCs vs Clinical"),
    p_value = c(signif(pval1, 3), signif(pval2, 3), signif(pval3, 3), signif(pval4, 3))
  )
  
  # plot curves
  plot(roc_covar, col = "#006600", cex.main = 1.2, font.main = 1, cex.lab = 1.2)
  plot(roc_PRSmix, col = "#FF6600", add = TRUE)
  plot(roc_covar_PRSmix, col = "#CC0000", add = TRUE)
  plot(roc_prs_to_compare, col = "#0033FF", add = TRUE)
  plot(roc_covar_prs_to_compare, col = "#00CCFF", add = TRUE)
  
  labels <- list(
    bquote("Clinical Covariates (AUC = " * .(signif(auc(roc_covar), 3)) * ")"),
    bquote(PRS[plain(multi)] * " + PCs (AUC = " * .(signif(auc(roc_PRSmix), 3)) * ")"),
    bquote("Clinical + " * PRS[plain(multi)] * " + PCs (AUC = " * .(signif(auc(roc_covar_PRSmix), 3)) * ")"),
    bquote(PRS[plain(ratio)] * " + PCs (AUC = " * .(signif(auc(roc_prs_to_compare), 3)) * ")"),
    bquote("Clinical + " * PRS[plain(ratio)] * " + PCs (AUC = " * .(signif(auc(roc_covar_prs_to_compare), 3)) * ")")
  )
  
  legend("bottomright", legend = labels, col = c("#006600", "#FF6600", "#CC0000", "#0033FF", "#00CCFF"), lwd = 2, cex = 0.95, y.intersp = 0.85, inset = 0)
  text(1.25, 1, labels = plot_title, adj = c(0, 1), cex = 1.2)
  
  return(list(auc_table = auc_table, pval_table = pval_table))
}
# function to plot ROC only
plot_roc <- function(df, outcome_binary, prs, pcs, covars, adjustment, plot_title) {
  # clinical covariates
  model1 <- build_model(df, outcome_binary, prs=NULL, pcs=NULL, covars, model_type='binomial', adjustment)
  pred_prob1 <- predict(model1, type = "response")
  roc_covar <- roc(df[[outcome_binary]], pred_prob1)
  
  # PRS-mix + PCs
  model2 <- build_model(df, outcome_binary, prs, pcs, covars=NULL, model_type='binomial', adjustment)
  pred_prob2 <- predict(model2, type = "response")
  roc_PRSmix <- roc(df[[outcome_binary]], pred_prob2)
  
  # full model of clinical covariates + PRSmix + PCs
  model3 <- build_model(df, outcome_binary, prs, pcs, covars, model_type='binomial', adjustment)
  pred_prob3 <- predict(model3, type = "response")
  roc_full <- roc(df[[outcome_binary]], pred_prob3)
  
  # plot curves
  plot(roc_covar, col = "blue", main = plot_title, cex.main = 1.2, font.main = 1)
  plot(roc_PRSmix, col = "red", add = TRUE)
  plot(roc_full, col = "green", add = TRUE)
  
  legend("bottomright", legend = c("Clinical Covariates", "PRS-mix + PCs", "Full Model"), 
         col = c("blue", "red", "green"), lwd = 2)
}
# 
# ################## test in other cohorts functions #####################
# 
# function to call scale other cohorts prs with copdgene PRSmix weights (no training)
scale_PRSmix_other_coh <- function(df, outcome_name_copd, coh) {
  out_dir <- "path/to/data/PRSmix_out_files"
  
  # scale coefficients with results from PRSmix
  weight_file <- paste0(out_dir, "/", coh, "_", outcome_name_copd, "_power.0.95_pthres.0.05_weight_PRSmixPlus.txt")
  mix_prs_name <- paste0("mix_prs_", outcome_name_copd)
  df[[mix_prs_name]] <- 0

  # Check if the weight file exists
  if (file.exists(weight_file)) {
    # File exists, read the weights and process as usual
    weights <- read.table(weight_file, header = TRUE)

    for (i in 1:nrow(weights)) {
      prs_name <- weights$topprs[i]
      weight <- weights$ww[i]
      df[[mix_prs_name]] <- df[[mix_prs_name]] + df[[prs_name]] * weight
    }

    # if there is no weight file
  } else {
    best_acc <- read.table(paste0(out_dir, "/", coh, "_", outcome_name_copd, "_best_acc.txt"), header = TRUE, sep = "\t", na.strings = "")
    best_prs <- best_acc$pgs
    df[[mix_prs_name]] <- df[[best_prs]]
  }

  df[[mix_prs_name]] <- scale(df[[mix_prs_name]])

  return(df)
}


# function to call linear models on scale other cohorts prs with copdgene PRSmix weights
call_model_other_coh <- function(df, outcome_name_copd, outcome_name_other_coh, pcs, covars, model_type = "lm", adjustment = NULL, coh, additional_prs = NULL) {
  out_dir <- "path/to/data/PRSmix_out_files"
  
  weight_file <- paste0(out_dir, "/", coh, "_", outcome_name_copd, "_power.0.95_pthres.0.05_weight_PRSmixPlus.txt")
  mix_prs_name <- paste0("mix_prs_", outcome_name_copd)

  results <- data.frame(prs = character(), beta = numeric(), ci_lower = numeric(), ci_upper = numeric(),
                        beta_ci = character(), se = numeric(), stringsAsFactors = FALSE)

  # Check if the weight file exists
  if (file.exists(weight_file)) {
    weights <- read.table(weight_file, header = TRUE)
    # get prss that were selected, and add mix_prs
    prss <- weights$topprs
    prss <- c(prss, mix_prs_name, additional_prs)

    for (prs in prss) {
      model <- build_model(df, outcome_name_other_coh, prs, pcs, covars, model_type, adjustment)
      beta <- coef(model)[prs]
      ci <- confint(model, prs, level = 0.95)
      beta_ci <- sprintf("%.3f (%.3f, %.3f)", beta, ci[1], ci[2])
      se <- summary(model)$coefficients[prs, "Std. Error"]
      pval <- summary(model)$coefficients[prs, 4]

      # prs <- ifelse(grepl("^mix_prs_", prs), "mix_prs", prs) #rename
      results <- rbind(results, data.frame(prs = prs, beta = beta, ci_lower = ci[1], ci_upper = ci[2],
                                           beta_ci = beta_ci, se = se, pval = pval), row.names = NULL)
    }

    # if there is no weight file, use the best one
  } else {
    best_acc <- read.table(paste0(out_dir, "/", coh, "_", outcome_name_copd, "_best_acc.txt"), header = TRUE, sep = "\t", na.strings = "")
    prs <- best_acc$pgs
    model <- build_model(df, outcome_name_other_coh, prs, pcs, covars, model_type, adjustment)
    beta <- coef(model)[prs]
    ci <- confint(model, prs, level = 0.95)
    beta_ci <- sprintf("%.3f (%.3f, %.3f)", beta, ci[1], ci[2])
    se <- summary(model)$coefficients[prs, "Std. Error"]
    pval <- summary(model)$coefficients[prs, 4]

    results <- data.frame(prs = prs, beta = beta, ci_lower = ci[1], ci_upper = ci[2],
                          beta_ci = beta_ci, se = se, pval = pval, row.names = NULL)
  }
  # colnames(results)[colnames(results) == "beta_ci"] <- paste0(outcome_name_other_coh, "_weighted_by_", outcome_name_copd)
  # names(results)[2] <- paste0(outcome_name_other_coh, "_weighted_by_", outcome_name_copd)
  return(results)
}