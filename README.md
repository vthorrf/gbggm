gbggm: Generalized Bayesian Gaussian Graphical Models
=============

This package includes some implementations of Bayesian Gaussian Graphical Models based on the generalized approach proposed by Franco et al. (under review). Currently, there are 8 models implemented:
* Priors with one parameter for regularization (i.e., models that control the regularization with a single parameter): `"normal"`, `"laplace"`, `"logistic"`, `"cauchy"`, and `"hypersec"`. Respectively, these values set a normal, laplace, logistic, Cauchy, or hyperbolic secant as prior distributions;
* Priors with two parameters for regularization (i.e., models that control the regularization with two parameters: a "regularization" parameter per se and an "heavy-tailedness" parameter, which may be useful when there are "outlier" correlations): `"t"`, `"lomax"`, and `"NEG"`. Respectively, these values set a t, double lomax, or normal-exponential-gamma as prior distributions.

The user can also decide if they will estimate a sparse or non-sparse network. The estimation of the BGGMs is done with the `bggm` function.

# Installation #

Using the 'remotes' package:

    install.packages("remotes")
    remotes::install_github("vthorrf/gbggm")
