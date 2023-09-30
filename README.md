gbggm: Generalized Bayesian Gaussian Graphical Models
=============

This package includes some implementations of Bayesian Gaussian Graphical Models based on the generalized approach proposed by Franco et al. (under review). Currently, there are 6 models implemented:
* One-parameters regularization (i.e., models that control the regularization with a single parameter): "normal", "laplace", "logistic", and "cauchy";
* Two-parameters regularization (i.e., models that control the regularization with two parameter: a "regularization" parameter per se and an "heavy-tailedness" parameter, which may be useful when there are "outlier" correlations): "t", and "lomax".

The user can also decide if they will estimate a sparse or non-sparse network. The estimation of the BGGMs is done with the `bggm` function.

# Installation #

Using the 'remotes' package:

    install.packages("remotes")
    remotes::install_github("vthorrf/gbggm")
