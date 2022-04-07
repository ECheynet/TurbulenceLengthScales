# A note on the calculation of the length scales of wind turbulence

Matlab functions to compute the integral and turbulence length scales in the atmospheric boundary layer.

[![View Estimation of the length scales of wind turbulence on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/108944-estimation-of-the-length-scales-of-wind-turbulence)


## Summary

The Matlab functions used in this small toolbox are used to estimate some of the turbulence length scales in the atmospheric boundary layer [1]. The velocity fluctuations used should be stationary. The assumption of homogeneity is required for the study of the crosswind turbulence length scales. The study of these turbulence length scales can be used for wind engineering applications or fundamental turbulence analysis.  The toolbox includes two methods to compute the integral length scales: one method relying on the autocovariance function and another using the von Kármán spectrum.

The random error (or statistical uncertainty) is quite large for the integral length scale of the along-wind component. So its application for wind engineering application should be considered with caution [1]


## Content

The submission file contains:

- The function Lyz.m, which computes the crosswind turbulence length scale, either in the y-direction or the z-direction
- The function Ly.m, which computes the integral length scale using the auto-covariance function
- The function fitVK.m computes the integral length scale using a least-square fit of the von Kármán spectrum to the estimated power spectral densities of the velocity fluctuations.
- An example file using Lyz.m Crosswind_turbulence_length_scale.mlx
- An example file using Ly.m and fitVK.mExample_integral_length_scales.mlx
- A data file data.mat used in the example files.


## References

[1] Cheynet, E. (2016). Wind-induced vibrations of a suspension bridge: A case study in full-scale. PhD thesis. University of Stavanger, Norway
