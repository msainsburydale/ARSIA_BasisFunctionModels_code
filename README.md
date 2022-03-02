# ARSIA review of basis-function models

As part of our review of basis-function models in spatial statistics (see [here](https://www.annualreviews.org/doi/abs/10.1146/annurev-statistics-040120-020733) for the published manuscript, and see [here](https://arxiv.org/abs/2202.03660) for an arXiv version of the manuscript), we compared the spatial predictions obtained using the `R` packages `FRK` (Zammit-Mangion and Cressie, 2021; Sainsbury-Dale et al., 2021), `INLA` (Rue et al., 2009), `LatticeKrig` (Nychka, 2015), `mgcv` (Wood, 2017), and `gstat` (Pebesma, 2004), as well as a `C++` implementation of the multi-resolutional approximation (MRA; Huang et al., 2019).  The data consists of sea-surface temperatures (SST) in the Brazil-Malvinas confluence zone, shown below. This is Figure 4 of the manuscript.

![Figure 4: Sea-surface temperature data](/img/global_and_training_data.png?raw=true)

The predictions and standard errors obtained using each package are shown below. This is Figure 5 of the manuscript.

![Figure 5: Predictions and standard errors](/img/DEM_results.png?raw=true)


## Dependencies

The package versions used by this repo are:
- `FRK` version 2.0.1,  
- `LatticeKrig` version 8.4,
- `mgcv` version 1.8.31,
- `gstat` version 2.0.5,
-  `MRA` from [this](https://github.com/hhuang90/MRA_For_NCAR_Technical_Note) repo. 

Note that it is particularly important that `FRK` is up-to-date, as the code will not work for any version less than 2.0.0.

## Instructions

To reproduce the results of the manuscript, please download this repository (see [here](https://superuser.com/a/1309684) for steps to download a repository). Open SST_analysis.R; this is the controlling script for the analysis. The code populates the img/ and results/ directories, the contents of which are either used in the paper or by subsequent code.

Within SST_analysis.R, first load the required packages, and then enter the path to the directory containing SST_analysis.R in the DIRECTORY variable. The first stages of this script consist of data pre-processing and visualisation, yielding Figure 4 in the manuscript.

Next, the script loads the model fitting and prediction functions, which are kept in modelling_functions/, and predicts the SST at the testing locations. The data frame containing the testing data, predictions, and prediction standard errors is saved as results/df_test.csv. One may load this data frame directly to skip model fitting and prediction (commands for loading df_test.csv are included in the SST_analysis.R script). Using df_test, out-of-sample diagnostics are produced, as given in Table 1 of the manuscript.

Finally, using the previously fitted model objects, the script generates predictions over the spatial domain, D, yielding Figure 5 of the manuscript. These results are saved in results/grid_over_D.csv. Again, one may load this data frame directly to skip model fitting and prediction, and commands for doing so are provided in the script.


#### References

* Huang H, Blake LR, Hammerling DM. 2019. Pushing the limit: A hybrid parallel implementation of the multi-resolution approximation for massive data. arXiv preprint arXiv:1905.00141
* Nychka D, Bandyopadhyay S, Hammerling D, Lindgren F, Sain S. 2015. A multiresolution Gaussian process model for the analysis of large spatial datasets. Journal of Computational and Graphical Statistics 24:579–599
* Pebesma EJ. 2004. Multivariable geostatistics in S: the gstat package. Computers \& Geosciences 30:683–691
* Rue H, Martino S, Chopin N. 2009. Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations. Journal of the Royal Statistical Society 71:319–392
* Sainsbury-Dale M, Zammit-Mangion A, Cressie N. 2021. Modelling, fitting, and prediction with non-Gaussian spatial and spatio-temporal data using FRK. arXiv:2110.02507
* Wood SN. 2017. Generalized Additive Models: An Introduction with R. Boca Raton, FL: Chapman \& Hall/CRC, 2nd ed.
* Zammit-Mangion A, Cressie N. 2021. FRK: An R package for spatial and spatio-temporal prediction
with large datasets. Journal of Statistical Software 98(4):1–48
