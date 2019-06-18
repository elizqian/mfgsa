# Multifidelity global sensitivity analysis

This code implements the multifidelity approach to estimating variance and sensitivity indices as described in:

*  Qian, E., Peherstorfer, B., O'Malley, D., Vesselinov, V., and Willcox, K. 
[Multifidelity Monte Carlo Estimation of Variance and Sensitivity Indices](https://www.dropbox.com/s/y77c42t9po52384/QPOVW_mfgsa_juq2018.pdf?dl=0).
SIAM/ASA Journal on Uncertainty Quantification, 6(2):683-706, 2018.

The file `main.m` computes multifidelity variance and sensitivity index estimates for the Ishigami function example considered in the paper.

## Related papers 
* Peherstorfer, B., Willcox, K., and Gunzburger, M. [Optimal model management for multifidelity Monte Carlo estimation](https://pehersto.engr.wisc.edu/preprints/multi-fidelity-monte-carlo-peherstorfer-willcox-gunzburger.pdf).
SIAM Journal on Scientific Computing, 38(5):A3163-A3194, 2016. ([MFMC GitHub](https://github.com/pehersto/mfmc))