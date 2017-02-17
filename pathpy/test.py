# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 11:59:22 2015
@author: Ingo Scholtes

(c) Copyright ETH Zurich, Chair of Systems Design, 2015-2017
"""

import pathpy as pp
import numpy as _np

#########################
# TEST PATHWAY MODELING #
#########################

# Example without second-order correlations
paths = pp.Paths()

paths.addPath('a,c')
paths.addPath('b,c')
paths.addPath('c,d')
paths.addPath('c,e')

for k in range(3):
    paths.addPath('a,c,d')
    paths.addPath('b,c,e')
    paths.addPath('b,c,d')
    paths.addPath('a,c,e')

m = pp.MultiOrderModel(paths, maxOrder=2)
assert m.estimateOrder(paths) == 1, "Error, wrongly detected higher-order correlations"

# Example with second-order correlations
paths = pp.Paths()

paths.addPath('a,c')
paths.addPath('b,c')
paths.addPath('c,d')
paths.addPath('c,e')

for k in range(3):
    paths.addPath('a,c,d')
    paths.addPath('b,c,e')   

m = pp.MultiOrderModel(paths, maxOrder=2)
assert m.estimateOrder(paths) == 2, "Error, did not detect second-order correlations"

x = list(map(str, _np.random.choice(range(10), 100000)))
ms = pp.MarkovSequence(x)
assert ms.estimateOrder(maxOrder=2, method='BIC') == 1, "Error, wrongly detected higher-order correlations"
assert ms.estimateOrder(maxOrder=2, method='AIC') == 1, "Error, wrongly detected higher-order correlations"

#########################
# TEST TEMPORAL NETWORK #
#########################

t = pp.TemporalNetwork()
# Path of length two
t.addEdge("c", "e", 1)
t.addEdge("e", "f", 2)

# Path of length two
t.addEdge("a", "e", 3)
t.addEdge("e", "g", 4)

# Path of length two
t.addEdge("c", "e", 5)
t.addEdge("e", "f", 6)

# Path of length two
t.addEdge("a", "e", 7)
t.addEdge("e", "g", 8)

# Path of length two
t.addEdge("c", "e", 9)
t.addEdge("e", "f", 10)

# The next two edges continue the previous path to ( c-> e-> f-> e -> b )
t.addEdge("f", "e", 11)
t.addEdge("e", "b", 12)

# This is an isolated edge (i.e. path of length one)
t.addEdge("e", "b", 13)

# And multiple edges in a single time step, which should result 
# in two paths of length two
t.addEdge("g", "e", 14)
t.addEdge("c", "e", 14)
t.addEdge("e", "f", 15)

# Path of length two
t.addEdge("b", "e", 16)
t.addEdge("e", "g", 17)

# Path of length two
t.addEdge("c", "e", 18)
t.addEdge("e", "f", 19)

# Path of length two
t.addEdge("c", "e", 20)
t.addEdge("e", "f", 21)

# Extract (time-respecting) paths
paths = pp.Paths.fromTemporalNetwork(t, delta=1)

print("Test network has", paths.ObservationCount(), "time-respecting paths")
assert paths.ObservationCount() == 13, "Extracted wrong number of time-respecting paths"

# TODO: Compute betweenness preference of nodes
bw = 1.2954618442383219

assert bw == 1.2954618442383219