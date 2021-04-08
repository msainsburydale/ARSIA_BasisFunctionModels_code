# ARSIA review of basis-function models

As part of our review of basis-function models in spatial statistics, we compared the spatial predictions obtained using the `R` packages `FRK`, `INLA`, `LatticeKrig`, `mgcv`, and `gstat`, as well as a `C++` implementation of the multi-resolutional approximation (Huang et al., 2019).  The data consists of sea-surface temperatures in the Brazil-Malvinas confluence zone, shown in the following.

![Sea-surface temperature data](/img/global_and_training_data2.png?raw=true)

The predictions and standard errors are shown in the following.

![Predictions and standard errors](/img/DEM_results.png?raw=true)



* Huang H, Blake LR, Hammerling DM. 2019. Pushing the limit: A hybrid parallel implementation of the multi-resolution approximation for massive data. arXiv preprint arXiv:1905.00141