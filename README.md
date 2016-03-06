# MID_LNPestimator
Maximally Informative Dimensions estimator: nonparametric estimator
for LNP model 

**Description:** Estimates a set of linear filters from a stimulus and
 spike train using the *maximally informative dimension* (MID)
 estimator (introduced in Sharpee *et al* 2004), which is a
 nonparametric maximum-likelihood estimator for the parameters of a
 linear-nonlinear-Poisson (LNP) neural response model.

**Relevant publication:**
Williamson, Sahani & Pillow, *PLoS Comp Bio*
2015. [[abs](http://pillowlab.princeton.edu/pubs/abs_Williamson15_PLoSCB.html),
      [pdf](http://pillowlab.princeton.edu/pubs/Williamson_etal_plosCB2015.pdf)]

Download
==========

* **Command line**: clone the repository from github (e.g., ```git
  clone git@github.com:pillowlab/MID_nonparLNP.git``` )
* **Browser**:  download zipped archive:  [MID_nonparLNP-master.zip](https://github.com/pillowlab/MID_nonparLNP/archive/master.zip)


Usage
=====

* Launch matlab and cd into the directory containing the code
 (e.g. `cd code/MIDforLNP/`).

* Examine the script `test.m` for a line-by-line tutorial on how to
use the code contained in this package, which goes through several
simulated examples.


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
