# LNPfitting
Linear-Nonlinear-Poisson (LNP) model fitting via maximum likelihood, aka Maximally-Informative-Dimensions (MID).

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
```git clone git@github.com:pillowlab/LNPfitting.git```
* **Browser**:  download zipped archive:  [LNPfitting-master.zip](https://github.com/pillowlab/LNPfitting/archive/master.zip)


Usage
=====

* Launch matlab and cd into the directory containing the code
 (e.g. `cd code/LNPfitting/`).

* Examine the demo scripts for annotated example analyses of simulated
datasets: 
	*  `demo1_1temporalfilter.m` - estimate a single temporal filter
    using STA, iSTAC, GLM with exponential nonlinearity, and using
    ML / MID estimation of filter and a nonlinearity parametrized
    by radial basis functions. Shows validation using R^2 of true
	filter and single-spike information (log-likelihood).
	*  `demo2_2temporalfilters.m` - similar, but for two-filter LNP
	models; shows how to plot inferred 2D nonlinearity.
	* `demo3_3temporalfilters` - compares iSTAC, and ML/MID estimates with cylindrical basis
	functions (CBFs) and with radial basis functions (RBFs) for
	parametrizing the nonlinearity.
	* `demo4_1spacetimefilter.m` - similar to demo1 but for 2D (space x
	time) binary white noise stimulus.
	* `demo5_5spacetimefilters.m` - similar to demo3, but recovers 5
	space-time filters. Compares iSTAC, and ML/MID with CBFs and RBFs.



Notes
=====

Primary differences between our implementation and that of Sharpee et
al 2004:

* Parametrizes the nonlinearity with smooth RBF
  (or CBF) basis functions (followed by a rectifying point
  nonlinearity to keep spike rates positive) instead
  of square histogram bins. This makes the log-likelihood
  differentiable and therefore easier to ascend.

* Uses standard quasi-Newton optimization methods (fminunc or minFunc)
  instead of simulated annealing. (Thus no guarantee that it finds
  global optimum, but runtime is substantially faster and performance
  is high in simulated and real datasets considered so far).
