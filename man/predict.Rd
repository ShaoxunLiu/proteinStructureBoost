% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PSB_functions_2.R
\name{predict}
\alias{predict}
\title{Predict results using adaboostTrain model and featureMatrix}
\usage{
predict(model, featureMatrix = NA)
}
\arguments{
\item{model}{Output of adaboostTrain}

\item{featureMatrix}{Output of featureMatrix. Do not enter for negative in
featureMatrix function. Test samples have label 1 as place holder.}
}
\value{
A boosting prediction
}
\description{
This fucntion gives boosting prediction results using adaboostTrain model
and featureMatrix
}
