% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PSB_functions.R
\name{adaboostModel}
\alias{adaboostModel}
\title{Train a adaboost model using feature Matrix}
\usage{
adaboostModel(
  set = NA,
  Train = NA,
  Test = NA,
  no.Aug = NA,
  positive.label = 1,
  negative.label = 0,
  iteration = 10,
  mfinal = 100,
  mfinal.sequence = NA,
  threshold = 0.5,
  verbose = FALSE
)
}
\arguments{
\item{set}{Set of train and test data.frames, output of trainTest.}

\item{Train}{The train data.frame. Only enter when set is not provided.}

\item{Test}{The test data.frame. Only enter when set is not provided.}

\item{no.Aug}{Numebr of augmentations plus the original. Only enter when set
is not provided.}

\item{positive.label}{The positive label of the data set (used in
featureMatrix function).}

\item{negative.label}{The negative label of the data set (used in
featureMatrix function).}

\item{iteration}{Number of iterations in training round.}

\item{mfinal}{The mfinal value in boosting function. See more in
?adabog::boosting}

\item{mfinal.sequence}{The list of mfinal values used. Length must be the
same as iteration. Only enter when mfinal is not provided.}

\item{threshold}{The proportion of positive votes to be seen as consensus.}

\item{verbose}{If true, prints the error of each iteration.}
}
\value{
A vector of voted results.
}
\description{
This function takes the data set generated by featureMatrix and processed by
trainTest and trains for a optimized adaboost model.
}
