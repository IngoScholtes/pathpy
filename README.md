# Introduction

`pathpy` is an OpenSource python package for the modeling and analysis of pathways and temporal networks
using higher-order and multi-order graphical models. 

The package is specifically tailored to analyze sequential data multiple observations of short, independent paths 
observed in a graph topology. Examples for such data include user click streams in information networks,
biological pathways, or traces of information propagation in social media. It also supports the extraction 
of pathways from data on temporal networks, such as time-stamped social interactions or contact patterns.

`pathpy` facilitates the analysis of temporal correlations in such data, using principled model selection 
techniques to infer higher-order graphical representatiosn that capture both topological and temporal 
characteristics of time-resolved relational data.

This methods implemented in this package have been introduced in the following research works: 

1. I Scholtes: [When is a network a network? Multi-Order Graphical Model Selection in Pathways and Temporal Networks](http://arxiv.org/), arXiv preprint
2. I Scholtes, N Wider, A Garas: [Higher-Order Aggregate Networks in the Analysis of Temporal Networks: Path structures and centralities](http://dx.doi.org/10.1140/epjb/e2016-60663-0), The European Physical Journal B, 89:61, March 2016
3. I Scholtes, N Wider, R Pfitzner, A Garas, CJ Tessone, F Schweitzer: [Causality-driven slow-down and speed-up of diffusion in non-Markovian temporal networks](http://www.nature.com/ncomms/2014/140924/ncomms6024/full/ncomms6024.html), Nature Communications, 5, September 2014
4. R Pfitzner, I Scholtes, A Garas, CJ Tessone, F Schweitzer: [Betweenness preference: Quantifying correlations in the topological dynamics of temporal networks](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.110.198701), Phys Rev Lett, 110(19), 198701, May 2013

The module is written in pure python, has no platform-specific dependencies and should thus work on all platforms. it builds on numpy and scipy. 

The latest development version can be installed by typing:

`> pip install git+git://github.com/IngoScholtes/pathpy.git`

## Tutorial

`pathpy` is the successor of `pyTempNets`. A detailed educational tutorial illustrating its theoretical foundation is [available online](https://www.sg.ethz.ch/team/people/ischoltes/research-insights/temporal-networks-demo/).

## Acknowledgements

The research behind this framework was funded by the Swiss State Secretariat for Education, Research and Innovation [(Grant C14.0036)}(https://www.sg.ethz.ch/projects/seri-information-spaces/). The development of this package was generously supported by the [MTEC Foundation](http://www.mtec.ethz.ch/research/support/MTECFoundation.html) in the context of the project [The Influence of Interaction Patterns on Success in Socio-Technical Systems: From Theory to Practice](https://www.sg.ethz.ch/projects/mtec-interaction-patterns/).

## Contributors

Ingo Scholtes (project lead, development)
Roman Cattaneo (development)

## Copyright

(c) Copyright ETH Zürich, Chair of Systems Design, 2015-2017