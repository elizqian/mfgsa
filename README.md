# Multifidelity global sensitivity analysis

## Repository overview

This code implements the multifidelity approach to estimating variance and global sensitivity indices as described in:

Qian, E., Peherstorfer, B., O'Malley, D., Vesselinov, V., and Willcox, K. 
[Multifidelity Monte Carlo Estimation of Variance and Sensitivity Indices](https://www.dropbox.com/s/y77c42t9po52384/QPOVW_mfgsa_juq2018.pdf?dl=0).
SIAM/ASA Journal on Uncertainty Quantification, 6(2):683-706, 2018.

Two examples are provided:
* The `ishigami` folder contains everything needed to run the numerical experiments for the analytical Ishigami function of the above paper. 
* The `CDR` folder contains a `samples.mat` file with pre-computed y<sub>A</sub>, y<sub>B</sub>, and y<sub>C</sub> samples and bootstraps from these samples to compute the sensitivity index estimates. 

## Global sensitivity analysis
When a model has uncertain inputs, the model output is also uncertain. Variance-based global sensitivity analysis quantifies the relative influence of each of these uncertain inputs on the output by dividing the total variance into percentages of the variance due to each of the inputs and due to interactions between inputs. For example, the Ishigami function has 3 random inputs, Z<sub>1</sub>, Z<sub>2</sub>, and Z<sub>3</sub>, each distributed uniformly on the interval [-&#120587;,&#120587;]:


## Related papers 
* Peherstorfer, B., Willcox, K., and Gunzburger, M. [Optimal model management for multifidelity Monte Carlo estimation](https://pehersto.engr.wisc.edu/preprints/multi-fidelity-monte-carlo-peherstorfer-willcox-gunzburger.pdf).
SIAM Journal on Scientific Computing, 38(5):A3163-A3194, 2016. ([MFMC GitHub](https://github.com/pehersto/mfmc))

## Contact
Please feel free to contact [Elizabeth Qian](http://www.elizabethqian.com/) with any questions about this repository.
