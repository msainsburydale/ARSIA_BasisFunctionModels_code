# ARSIA review of basis-function models

As part of our review of basis-function models in spatial statistics, we compared the spatial predictions obtained using the `R` packages `FRK`, `INLA`, `LatticeKrig`, `mgcv`, and `gstat`, as well as a `C++` implementation of the multi-resolutional approximation (Huang et al., 2019).  The data consists of sea-surface temperatures in the Brazil-Malvinas confluence zone, shown in the following.

![Sea-surface temperature data](/img/global_and_training_data.png?raw=true)

The predictions and standard errors are shown in the following.

![Predictions and standard errors](/img/DEM_results.png?raw=true)


## Instructions

To reproduce the results please download this repository. Open SST_analysis.R; this is the controlling script for the entire analysis. The code populates the img/ and results/ directories, the contents of which are either used in the paper or by subsequent code. 

WithinSST_analysis.R, first load the required packages, and then enter the path to this repository in the DIRECTORY variable. The first stages of this script consist of data pre-processesing and visualisation, yielding Figure 4 in the manuscript. 

Next, it loads the model fitting and prediction functions, which are kept in ./modelling_functions/, and produces predictions the testing locations. The data frame containing the testing data, predictions, and prediction standard errors, is saved in results/df_test.csv. One may load this data frame directly if one wishes to skip model fitting and prediction (commands for loading df_test.csv are included in the SST_analysis.R script). Using df_test, out-of-sample diagonstics are produced, as given in Table 1 of the manuscript.

Finally, using the previously fitted model object, the script generates predictions over the spatial domain, D, yielding Figure 5 of the manuscript. These results are saved in ./results/grid_over_D.csv. Again, one may load this data frame directly if one wishes to skip model fitting and prediction, and commands for doing so are provided in the script. 


#### References

* Huang H, Blake LR, Hammerling DM. 2019. Pushing the limit: A hybrid parallel implementation of the multi-resolution approximation for massive data. arXiv preprint arXiv:1905.00141