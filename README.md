# ARSIA review of basis-function models

As part of our review of basis-function models in spatial statistics, we compared the spatial predictions obtained using the `R` packages `FRK` (Zammit-Mangion, 2021), `INLA` (Rue et al., 2009), `LatticeKrig` (Nychka, 2015), `mgcv` (Wood, 2017), and `gstat` (Pebesma, 2004), as well as a `C++` implementation of the multi-resolutional approximation (Huang et al., 2019).  The data consists of sea-surface temperatures in the Brazil-Malvinas confluence zone, shown below. This is Figure 4 of the manuscript. 

![Figure 4: Sea-surface temperature data](/img/global_and_training_data.png?raw=true)

The predictions and standard errors obtained using each package are shown below. This is Figure 5 of the manuscript.

![Figure 5: Predictions and standard errors](/img/DEM_results.png?raw=true)


## Instructions

To reproduce the results please download this repository. Open SST_analysis.R; this is the controlling script for the entire analysis. The code populates the img/ and results/ directories, the contents of which are either used in the paper or by subsequent code. 

Within SST_analysis.R, first load the required packages, and then enter the path to this repository in the DIRECTORY variable. The first stages of this script consist of data pre-processesing and visualisation, yielding Figure 4 in the manuscript. 

Next, it loads the model fitting and prediction functions, which are kept in ./modelling_functions/, and produces predictions the testing locations. The data frame containing the testing data, predictions, and prediction standard errors, is saved in results/df_test.csv. One may load this data frame directly if one wishes to skip model fitting and prediction (commands for loading df_test.csv are included in the SST_analysis.R script). Using df_test, out-of-sample diagonstics are produced, as given in Table 1 of the manuscript.

Finally, using the previously fitted model objects, the script generates predictions over the spatial domain, D, yielding Figure 5 of the manuscript. These results are saved in ./results/grid_over_D.csv. Again, one may load this data frame directly if one wishes to skip model fitting and prediction, and commands for doing so are provided in the script. 


#### References

* Huang H, Blake LR, Hammerling DM. 2019. Pushing the limit: A hybrid parallel implementation of the multi-resolution approximation for massive data. arXiv preprint arXiv:1905.00141
* Nychka D, Bandyopadhyay S, Hammerling D, Lindgren F, Sain S. 2015. A multiresolution Gaussian process model for the analysis of large spatial datasets. Journal of Computational and Graphical Statistics 24:579–599
* Pebesma EJ. 2004. Multivariable geostatistics in S: the gstat package. Computers \& Geosciences 30:683–691
* Rue H, Martino S, Chopin N. 2009. Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations. Journal of the Royal Statistical Society 71:319–392
* Wood SN. 2017. Generalized Additive Models: An Introduction with R. Boca Raton, FL: Chapman \& Hall/CRC, 2nd ed.
* Zammit-Mangion A, Cressie N. 2021. FRK: An R package for spatial and spatio-temporal prediction
with large datasets. Journal of Statistical Software In press