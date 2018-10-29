Analytical kernel vine copula Toolbox for Matlab
=============================

Toolbox for canonical vine copula trees with kernel copula.


Description
-----------



Demonstration
-------------

The script `demo_kernel.m` demonstrates how to apply the Toolbox. 

Main functions
---------

demo_kernel.m

Fit_vCopula.m      - Fit and evaluate kernel vine-Copula to the data

KernelCop.m        - Fit and evaluate single copula 

prep_vine.m        - order the nodes of the c or nc -vine according to the tau-algorithm. 

predict_response.m - Predict mle and em estimates for the response.

Plot_vCopula.m     - Plot various structures of the nodes of the vine-copula

kernelpdf.m        - 1d marginal pdf's

kernelcdf.m        - 1d marginal cdf's

#rand_vCopula.m     - vine-copula random generation

DenseNaiv.m        -Naive kernel density

Kernel_LL.m        -Analytical Local likelihood density


---
mixedvinefit      - Mixed copula vine estimates

mixedvinepdf      - Mixed copula vine probability density function

mixedvinernd      - Mixed copula vine random numbers

mixedvineentropy  - Mixed copula vine entropy estimate

mixedvineinfo     - Mixed copula vine mutual information estimate

mixedgaussfit     - Mixed copula vine estimates with Gaussian copula

marginfit         - Univariate margin estimates

marginpdf         - Univariate margin probability density function

margincdf         - Univariate margin cumulative distribution function

margininv         - Inverse of univariate margin CDF

copulafit         - Copula parameter estimates

copulapdf         - Copula probability density function

copulacdf         - Copula cumulative distribution function

copulaccdf        - Copula conditional cumulative distribution function

copulaccdfinv     - Inverse of copula conditional CDF

copularnd         - Copula random numbers



License
-------
