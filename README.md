# MIDestimator
"Maximally-Informative-Dimensions" (MID) estimator for the Linear-Nonlinear-Poisson (LNP) model

**Description:** Estimates the parameters of an LNP model (a set of linear filters and a nonlinear function) from a stimulus and
 spike train using the *maximally informative dimension* (MID) estimator (introduced in Sharpee *et al* 2004). MID is a nonparametric maximum-likelihood estimator for the parameters of a linear-nonlinear-Poisson (LNP) neural response model.

**Relevant publications:**

* [Sharpee, Rust, & Bialek, *Neural Computation* 2004](http://www.mitpressjournals.org/doi/abs/10.1162/089976604322742010) (original publication).
* Williamson, Sahani, & Pillow, *PLoS Comp Bio* 2015. [[abs](http://pillowlab.princeton.edu/pubs/abs_Williamson15_PLoSCB.html),
      [pdf](http://pillowlab.princeton.edu/pubs/Williamson_etal_plosCB2015.pdf), 
      [link](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004141)]  
      - Shows equivalence of MID to maximum-likelihood estimator for LNP

Download
==========

* **Command line**: clone the repository from github:
```git clone git@github.com:pillowlab/MIDestimator.git```
* **Browser**:  download zipped archive:  [MIDestimator-master.zip](https://github.com/pillowlab/MIDestimator/archive/master.zip)


Usage
=====

* Launch matlab and cd into the directory containing the code
 (e.g. `cd code/MIDestimator/`).

* Examine the tutorial scripts (`tutorial1_MID_tfilter.m`, ...) for detailed instructions and applications to simulated datasets.


Notes
=====

Primary differences between our implementation and that of Sharpee et
al 2004:

* This implementation parametrizes the nonlinearity with smooth RBF
  (or CBF) basis functions (followed by a point nonlinearity) instead
  of square histogram bins, which makes the log-likelihood
  differentiable and therefore easier to ascend.

* This implementation uses standard quasi-Newton methods to find a
  local optimum, as opposed to simulated annealing.
