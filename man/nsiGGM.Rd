\name{nsiGGM}
\alias{nsiGGM}
\title{Node-structured Integrative Gaussian Graphical Model}
\usage{nsiGGM(Y, lambda1, lambda2, lambda3, params.in_gr,
              params.out_gr, params.dummy, params.verbose,
              penalty = "group", rho = 1, penalize.diagonal = FALSE,
              maxiter = 50, tol = 1e-5, warm = NULL,
              truncate = 1e-5, is.Naive.Group = TRUE)
}
\arguments{Please note an example below.}
\value{A list contains information on the final model}
\description{
This function develops an estiamation and variable selection algorithm for
node-structured integrative gaussian graphical model (nsiGGM).
}
\details{
We propose a novel statistical framework called the
``node-structured joint Gaussian graphical model",
for fitting joint Gaussian graphical model simultaneously with informative pathways
consistently expressed across multiple studies.
With an application to simulated and breast cancer genomic data,
the proposed model is found to be superior in efficiently capturing transcriptional
modules pre-defined by pathway database.
}
\author{
SungHwan Kim, Jae-Hwan Jhong, JungJun Lee, Ja-Yong Koo, ByungYong Lee, SungWon Han
}
\examples{
library(nsiGGM)
library(doMC)
library(mnormt)
library(igraph)
registerDoMC(10)

# True precision matrix
path1 = system.file("extdata", "TR_PR_mat.Rdata", package = "nsiGGM")
TR_PR_mat = get(load(path1))
# Three data sets generated from the Barabasi Albert algorithm
Y <- list()
n.S <- 3
for(i in 1:n.S){
  Y[[i]] <- rmnorm(100, mean = rep(0,dim(TR_PR_mat)[1]), solve(TR_PR_mat))
}

# Node-group index list that contatins a single pathway geneset in the 449th slot in the list.
path2 = system.file("extdata", "out_gr.Rdata", package = "nsiGGM")
out_gr <- get(load(path2))

# Study-group index that captures all involved three studies.
in_gr = list(c(1,2,3))

params.in_gr = in_gr
params.out_gr = out_gr
params.verbose = 1;

lambda1 <- 0.05
lambda2 <- 0.05
lambda3 <- 0.05

rst_nsJGL <- nsiGGM(Y, lambda1, lambda2, lambda3,  params.in_gr = params.in_gr,
                    params.out_gr = params.out_gr, params.dummy = 1,
                    params.verbose = params.verbose, rho = 1,
                    penalize.diagonal = FALSE, maxiter = 50,
                    tol = 1e-5, warm = NULL, truncate = 1e-5)

tmp_theta <- rst_nsJGL$theta
}

