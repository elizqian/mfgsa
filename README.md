# Multifidelity global sensitivity analysis

This code implements the multifidelity approach to estimating variance and sensitivity indices as described in:

Qian, E., Peherstorfer, B., O'Malley, D., Vesselinov, V., and Willcox, K. 
[Multifidelity Monte Carlo Estimation of Variance and Sensitivity Indices](https://www.dropbox.com/s/y77c42t9po52384/QPOVW_mfgsa_juq2018.pdf?dl=0).
SIAM/ASA Journal on Uncertainty Quantification, 6(2):683-706, 2018.

Two examples are provided:
* The `ishigami` folder contains everything needed to run the numerical experiments for the analytical Ishigami function of the above paper. 
* The `CDR` folder contains a `samples.mat` file with pre-computed yA, yB, and yC samples and bootstraps from these samples to compute the sensitivity index estimates. 

## Related papers 
* Peherstorfer, B., Willcox, K., and Gunzburger, M. [Optimal model management for multifidelity Monte Carlo estimation](https://pehersto.engr.wisc.edu/preprints/multi-fidelity-monte-carlo-peherstorfer-willcox-gunzburger.pdf).
SIAM Journal on Scientific Computing, 38(5):A3163-A3194, 2016. ([MFMC GitHub](https://github.com/pehersto/mfmc))
