#' Call train from caret on a trainTest set
#'
#' This function retruns a caret model of named method using trainTest set data.
#'
#' @param set Set of train and test data returned from trainTest.
#' @param method Avaliable methods in caret::train().
#' @return A fit.caret model
#' @export
caretTrain <- function(set, method = 'adaboost'){
  fit.caret <- caret::train(
    make.names(Label) ~ .,
    data       = set$Train,
    method     = method,
    trControl  = caret::trainControl(method = "cv", number = 5, classProbs = TRUE),
    tuneLength = 1,
    importance = 'impurity'
  )
  return(fit.caret)
}

#' Perform feature importance analysis using Lime
#'
#' This function returns a feature importance analysis of caret model using LIME.
#'
#' @param fit.caret Fit model resulting from caretTrain.
#' @param set Set of train and test data returned from trainTest.
#' @param permutation Number of permutations in LIME
#' @return A LIME explanation for caret
#' @export
caretLime <- function(fit.caret, set, permutation = 100){
  explainer_caret <- lime::lime(set$Test, fit.caret, n_bins = 5)
  explanation_caret <- lime::explain(
    x               = set$Test,
    explainer       = explainer_caret,
    n_permutations  = permutation,
    dist_fun        = "gower",
    kernel_width    = .75,
    n_features      = 5,
    feature_select  = "auto",
    labels          = as.factor('1')
  )
  return(explanation_caret)
}

#' Get feature importance with LIME
#'
#' This function returns a feature importance analysis of caret model using LIME.
#'
#' @param method Avaliable methods in caret::train().
#' @param set Set of train and test data returned from trainTest.
#' @param permutation Number of permutations in LIME
#' @return A LIME explanation for caret
#' @export
lbpLIME <- function(set, method = 'adaboost', permutation = 100){
  model <- caretTrain(set, method)
  out <- caretLime(model, set, permutation)
  return(out)
}

#' Plot LIME feature importance
#'
#' This function plots the feature importance from LIME
#'
#' @param lbpLIME Output of pbpLIME
#' @return A feature importance plot
#' @export
plotLIME <- function(lbpLIME){
  out <- lime::plot_features(lbpLIME[1:min(length(lbpLIME), 50)])
  return(out)
}

#' Predict results using adaboostTrain model and featureMatrix
#'
#' This fucntion gives boosting prediction results using adaboostTrain model
#' and featureMatrix
#'
#' @param model Output of adaboostTrain
#' @param featureMatrix Output of featureMatrix. Do not enter for negative in
#' featureMatrix function. Test samples have label 1 as place holder.
#' @return A boosting prediction
#' @export
predict <- function(model, featureMatrix = NA){
  adaboostModel <- model$model
  tryCatch({
    adabag::boosting()
  }, error = function(e) e)
  options(warn=-1)
  if(is.na(featureMatrix)){
    out <- stats::predict(adaboostModel, attr(model, 'Test'))
  }else{
    out <- stats::predict(adaboostModel, featureMatrix)
  }
  options(warn=0)
  return(out)
}
