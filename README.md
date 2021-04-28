# Robust Learning of Fixed-Structure Bayesian Networks
A MATLAB implementation of [Robust Learning of Fixed-Structure Bayesian Networks](https://arxiv.org/abs/1606.07384).

Explanation of Files
===

Our experiments (`filter` directory)
---
* `filter_bn_sparse.m`: Synthetic experiments (Section 4.1).  Evaluate the performance of different algorithms (MLE, RANSAC, Filtering) when the graph structure of the ground-truth Bayes net is a random tree or a random graph.
* `filter_bn_alarm.m`: Semi-synthetic experiments (Section 4.2).  Evaluate the performance of different algorithms when the ground-truth Bayes net is the ALARM network (see below for more details).
* `empirical_bn.m`: Compute the empirical conditional probabilities of a known-structure Bayes net.
* `dtv_bn.m`: Estimate the total variation distance between two Bayes nets that have the same structure (Appendix C.2).

The ALARM Bayes net (`alarm` directory)
---
* `alarm.bif`: The ALARM network from [The ALARM Monitoring System: A Case Study with Two Probabilistic Inference Techniques for Belief Networks](https://link.springer.com/chapter/10.1007/978-3-642-93437-7_28) in BIF format, downloaded from https://www.bnlearn.com/bnrepository.
* `parseBIF.cpp`: Parse `alarm.bif` and convert the multi-valued ALARM network into binary-valued.  The binary Bayes net is used in `filterBN_ALARM.m`.

Reference
===
This repository is an implementation of the paper [Robust Learning of Fixed-Structure Bayesian Networks](https://arxiv.org/abs/1703.00893) which appeared in NeurIPS 2018, authored by [Yu Cheng](https://homepages.math.uic.edu/~yucheng/), [Ilias Diakonikolas](http://www.iliasdiakonikolas.org/), [Daniel M. Kane](https://cseweb.ucsd.edu/~dakane/), and [Alistair Stewart](http://www.alistair-stewart.com/).

```
@inproceedings{ChengDKS18,
  author    = {Yu Cheng and
               Ilias Diakonikolas and
               Daniel Kane and
               Alistair Stewart},
  title     = {Robust Learning of Fixed-Structure {B}ayesian Networks},
  booktitle = {Proceedings of the 32nd Conference on Neural Information Processing Systems (NeurIPS)},
  pages     = {10304--10316},
  year      = {2018}
}
```
