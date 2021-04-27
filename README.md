# Robust Learning of Fixed-Structure Bayesian Networks
A MATLAB implementation of [Robust Learning of Fixed-Structure Bayesian Networks](https://arxiv.org/abs/1606.07384).

Explanation of Files
===

Our algorithms
---


The ALARM Bayes net (`alarm` directory)
---
* `alarm.bif`: The ALARM Bayes net from [The ALARM Monitoring System: A Case Study with Two Probabilistic Inference Techniques for Belief Networks](https://link.springer.com/chapter/10.1007/978-3-642-93437-7_28) in BIF format, downloaded from https://www.bnlearn.com/bnrepository/.
* `parseBIF.cpp`: Parse `alarm.bif` and convert the multi-valued ALARM Bayes net into a binary one.  The binary Bayes net is used in `filterBN_ALARM.m`.

Reference
===
This repository is an implementation of the paper [Robust Learning of Fixed-Structure Bayesian Networks](https://arxiv.org/abs/1703.00893) appeared in NeurIPS 2018, authored by [Yu Cheng](https://homepages.math.uic.edu/~yucheng/), [Ilias Diakonikolas](http://www.iliasdiakonikolas.org/), [Daniel M. Kane](https://cseweb.ucsd.edu/~dakane/), and [Alistair Stewart](http://www.alistair-stewart.com/).

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
