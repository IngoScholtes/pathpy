<img src="https://github.com/IngoScholtes/pathpy/blob/master/pathpy_logo.png" width="300" alt="pathpy logo" />

# Introduction

`pathpy` is an OpenSource python package for the modeling and analysis of pathways and temporal networks
using **higher-order** and **multi-order** graphical models.

The package is specifically tailored to analyze sequential data which capture multiple observations of short, independent paths
observed in an underlying graph topology. Examples for such data include user click streams in information networks,
biological pathways, or traces of information propagating in social media. Unifying the analysis of pathways and temporal networks,
`pathpy` also supports the extraction of time-respecting paths from time-stamped network data. It extends (and will eventually supersede)
the package [`pyTempnets`](https://github.com/IngoScholtes/pyTempNets).

`pathpy` facilitates the analysis of temporal correlations in such sequential data. It uses a principled model selection
technique to infer higher-order graphical representations that capture both topological and temporal
characteristics of time-resolved relational data. It specifically allows to answer the question whether a (first-order) network
abstraction of such data is justified, or whether higher-order network abstractions are needed.

The theoretical foundation of this package, **higher-order network models**, has been developed in the following research works:

1. I Scholtes: [When is a network a network? Multi-Order Graphical Model Selection in Pathways and Temporal Networks](https://arxiv.org/abs/1702.05499), arXiv:1702.05499
2. I Scholtes, N Wider, A Garas: [Higher-Order Aggregate Networks in the Analysis of Temporal Networks: Path structures and centralities](http://dx.doi.org/10.1140/epjb/e2016-60663-0), The European Physical Journal B, 89:61, March 2016
3. I Scholtes, N Wider, R Pfitzner, A Garas, CJ Tessone, F Schweitzer: [Causality-driven slow-down and speed-up of diffusion in non-Markovian temporal networks](http://www.nature.com/ncomms/2014/140924/ncomms6024/full/ncomms6024.html), Nature Communications, 5, September 2014
4. R Pfitzner, I Scholtes, A Garas, CJ Tessone, F Schweitzer: [Betweenness preference: Quantifying correlations in the topological dynamics of temporal networks](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.110.198701), Phys Rev Lett, 110(19), 198701, May 2013

<img src="https://github.com/IngoScholtes/pathpy/blob/master/multiorder.png" width="300" alt="Illustration of Multi-Order Model" />

# Download and installation

The module is written in pure python. It has no platform-specific dependencies and should thus work on all platforms. It builds on `numpy` and `scipy`. The latest version of `pathpy` can be installed by typing:

`> pip install git+git://github.com/IngoScholtes/pathpy.git`

# Tutorial

A [comprehensive educational tutorial](https://ingoscholtes.github.io/pathpy/tutorial.html) which shows how you can use `pathpy` to analyze your data sets is [available online](https://ingoscholtes.github.io/pathpy/tutorial.html).
Moreover, a tutorial which illustrates the abstraction of **higher-order networks** in the modeling of dynamical processes in temporal networks is [available here](https://www.sg.ethz.ch/team/people/ischoltes/research-insights/temporal-networks-demo/). The
latter tutorial is based on the predecessor library `pyTempNets` but most of its features have already been included in `pathpy`.

# Documentation

The code is fully documented via docstrings which are accessible through python's built-in help system. Just type `help(SYMBOL_NAME)` to see the documentation of a class or method. A [reference manual is available here](https://ingoscholtes.github.io/pathpy/hierarchy.html).

# Acknowledgements

The research behind this data analysis framework was funded by the Swiss State Secretariat for Education, Research and Innovation [(Grant C14.0036)](https://www.sg.ethz.ch/projects/seri-information-spaces/). The development of this package was generously supported by the [MTEC Foundation](http://www.mtec.ethz.ch/research/support/MTECFoundation.html) in the context of the project [The Influence of Interaction Patterns on Success in Socio-Technical Systems: From Theory to Practice](https://www.sg.ethz.ch/projects/mtec-interaction-patterns/).

# Contributors

[Ingo Scholtes](http://www.ingoscholtes.net) (project lead, development)
Roman Cattaneo (development)
Nicolas Wider (testing)

# Copyright

(c) Copyright ETH Zürich, Chair of Systems Design, 2015-2017
