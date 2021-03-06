% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fcr_package.R
\docType{package}
\name{fcr-package}
\alias{fcr-package}
\alias{dynfcr}
\title{Dynamic prediction in functional concurrent regression with sparse functional covariates}
\description{
This package contains the functions for fitting dynamic functional concurrent regression
with sparse data.
}
\section{Notation}{


Let $y_{ij}$ denote some outcome measured at
$t_{ij}$ on the functional domain (e.g. time) for subject $i$ at observation
$j$.
We focus on fitting models of the form

\deqn{
 y_{ij} = f_0(t_{ij}) + f_1(t_{ij})X_{ij} + \cdots + b_i(t_{ij}) + \epsilon_{ij}
}
}

\section{Estimation}{


Estimation is performed using an iterative procedure described in Leroux et. al (2017).
Initially, a model is fit without $b_i(t_{ij})$. Using the residuals from this initial fit,
the covariance function is estimated. The model is then re-fit using this covariance function.
This procedure can be iterated as many times as desired.
}

\references{
Leroux A, Xiao L, Crainiceanu C, Checkley W (2017).
Dynamic prediction in functional concurrent regression with an application to child growth.
}
