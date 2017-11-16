# -*- coding: utf-8 -*-
"""
    pathpy is an OpenSource python package for the analysis of sequential data on pathways and temporal networks using higher- and multi order graphical models

    Copyright (C) 2016-2017 Ingo Scholtes, ETH Zürich

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact the developer:
    
    E-mail: ischoltes@ethz.ch
    Web:    http://www.ingoscholtes.net
"""

import collections as _co
import bisect as _bs
import itertools as _iter

import numpy as _np 

import scipy.sparse as _sparse
import scipy.sparse.linalg as _sla
import scipy.linalg as _la
import scipy as _sp

from pathpy.Log import Log
from pathpy.Log import Severity


class EmptySCCError(Exception):
    """
    This exception is thrown whenever a non-empty strongly 
    connected component is needed, but we encounter an empty one
    """
    pass


class HigherOrderNetwork:
    """
    Instances of this class capture a k-th-order representation of path statistics. Path statistics 
    can originate from pathway data, temporal networks, or from processes observed on top of a network topology.
    """


    def __init__(self, paths, k=1, separator='-', nullModel=False, method='FirstOrderTransitions', lanczosVecs=15, maxiter=1000):
        """
        Generates a k-th-order representation based on the given path statistics.

        @param paths: An instance of class Paths, which contains the path statistics to be 
            used in the generation of the k-th order representation 

        @param k: The order of the network representation to generate. For the default case of 
            k=1, the resulting representation corresponds to the usual (first-order) aggregate network, 
            i.e. links connect nodes and link weights are given by the frequency of each interaction. For 
            k>1, a k-th order node corresponds to a sequence of k nodes. The weight of a k-th order link 
            captures the frequency of a path of length k.

        @param separator: The separator character to be used in higher-order node names.

        @param nullModel: For the default value False, link weights are generated based on the statistics of 
            paths of length k in the underlying path statistics instance. If True, link weights are generated 
            from the first-order model (k=1) based on the assumption of independent links (i.e. corresponding) 
            to a first-order Markov model.

        @param method: specifies how the null model link weights in the k-th order model are calculated. 
            For the default method='FirstOrderTransitions', the weight w('v_1-v_2-...v_k', 'v_2-...-v_k-v_k+1') of 
            a k-order edge is set to the transition probability T['v_k', 'v_k+1'] in the first order network.
            For method='KOrderPi' the entry pi['v1-...-v_k'] in the stationary distribution of the 
            k-order network is used instead.
        """

        assert nullModel == False or (nullModel and k>1)

        assert method == 'FirstOrderTransitions' or method == 'KOrderPi', 'Error: unknown method to build null model'

        assert len(paths.paths.keys())>0 and max(paths.paths.keys())>=k, 'Error: constructing a model of order k requires paths of at least length k'
        
        ## The order of this HigherOrderNetwork
        self.order = k

        ## The paths object used to generate this instance
        self.paths = paths

        ## The nodes in this HigherOrderNetwork 
        self.nodes = []

        ## The separator character used to label higher-order nodes. 
        ## For separator '-', a second-order node will be 'a-b'.
        self.separator = separator

        ## A dictionary containing edges as well as edge weights
        self.edges = _co.defaultdict( lambda: _np.array([0,0]) )

        ## A dictionary containing the sets of successors of all nodes
        self.successors = _co.defaultdict( lambda: set() )

        ## A dictionary containing the sets of predecessors of all nodes
        self.predecessors = _co.defaultdict( lambda: set() )

        if k>1: 
            # For k>1 we need the first-order network to generate the null model
            # and calculate the degrees of freedom

            # For a multi-order model, the first-order network is generated multiple times!
            # TODO: Make this more efficient
            g1 = HigherOrderNetwork(paths, k=1)
            A = g1.getAdjacencyMatrix(includeSubPaths = True, weighted=False, transposed=True)

        if not nullModel:
            # Calculate the frequency of all paths of
            # length k, generate k-order nodes and set
            # edge weights accordingly
            for p in paths.paths[k]:
                # For a 0-order model, we generate a dummy start node
                if k==0: 
                    v = 'start'
                    w = p[0]
                else:
                    # Generate names of k-order nodes v and w
                    v = p[0]
                    w = p[1]
                    for l in range(1, k):
                        v = v + separator + p[l]
                        w = w + separator + p[l+1]
                self.nodes.append(v)
                self.nodes.append(w)

                # as edge weights of the k-th order model, we sum the
                # occurrence of paths of length k as subpath and longest path
                self.edges[(v,w)] += paths.paths[k][p]
                self.successors[v].add(w)
                self.predecessors[w].add(v)
            self.nodes = list(set(self.nodes)) # removes duplicates, does however randomize the order
            # Note: For all sequences of length k which (i) have never been observed, but (ii) 
            #       do actually represent paths of length k in the first-order network, we 
            #       may want to include some 'escape' mechanism along the lines of (Cleary and Witten 1994)                        

        else:
            # generate the *expected* frequencies of all possible 
            # paths based on independently occurring (first-order) links
            
            # generate all possible paths of length k 
            # based on edges in the first-order network
            possiblePaths = list(g1.edges.keys())

            E_new = list()
            for e1 in possiblePaths:
                for e2 in g1.edges:
                    if e1[-1] == e2[0]:
                        p = e1 + (e2[1],)
                        E_new.append(p)
            possiblePaths = E_new

            # validate that the number of unique generated paths corresponds to the sum of entries in A**k            
            assert (A**k).sum() == len(possiblePaths), 'Expected ' + str((A**k).sum()) + ' paths but got ' + str(len(possiblePaths))
            
            if method == 'KOrderPi':
                # compute stationary distribution of a random walker in the k-th order network
                g_k = HigherOrderNetwork(paths, k=k, separator = separator, nullModel = False)
                pi_k = HigherOrderNetwork.getLeadingEigenvector(g_k.getTransitionMatrix(includeSubPaths=True), normalized=True, lanczosVecs=lanczosVecs, maxiter=maxiter)
            else:
                # A = g1.getAdjacencyMatrix(includeSubPaths=True, weighted=True, transposed=False)
                T = g1.getTransitionMatrix(includeSubPaths=True)

            # assign link weights in k-order null model
            for p in possiblePaths:
                v = p[0]   
                w = p[1]
                # add k-order nodes and edges
                for l in range(1, k):
                    v = v + separator + p[l]
                    w = w + separator + p[l+1]
                self.nodes.append(v)
                self.nodes.append(w)

                # NOTE: under the null model's assumption of independent events, we
                # have P(B|A) = P(A ^ B)/P(A) = P(A)*P(B)/P(A) = P(B)
                # In other words: we are encoding a k-1-order Markov process in a k-order
                # Markov model and for the transition probabilities T_AB in the k-order model 
                # we simply have to set the k-1-order probabilities, i.e. T_AB = P(B)
                
                # Solution A: Use entries of stationary distribution, 
                # which give stationary visitation frequencies of k-order node w
                if method == 'KOrderPi':
                    self.edges[(v,w)] = _np.array( [ 0, pi_k[ g_k.nodes.index(w) ] ] )

                # Solution B: Use relative edge weight in first-order network
                # Note that A is *not* transposed
                # self.edges[(v,w)] = A[(g1.nodes.index(p[-2]),g1.nodes.index(p[-1]))] / A.sum()

                # Solution C: Use transition probability in first-order network
                # Note that T is transposed (!)
                elif method == 'FirstOrderTransitions':
                    p_vw = T[ ( g1.nodes.index(p[-1]), g1.nodes.index(p[-2]) ) ]
                    self.edges[(v,w)] =  _np.array( [0, p_vw] )  

                # Solution D: calculate k-path weights based on entries of squared k-1-order adjacency matrix

                # Note: Solution B and C are equivalent
                self.successors[v].add(w)
            self.nodes = list(set(self.nodes))

        # Compute degrees of freedom of models
        if k==0:
            # for a zero-order model, we just fit node probabilities (excluding the special 'start' node) 
            # Since probabilities must sum to one, the effective degree of freedom is one less than the number of nodes
            # This holds for both the paths and the ngrams model
            self.dof_paths = self.vcount() - 2
            self.dof_ngrams = self.vcount() - 2
        else:
            # for a first-order model, self is the first-order network
            if k==1:
                g1 = self
                A = g1.getAdjacencyMatrix(includeSubPaths = True, weighted=False, transposed=True)
           
            # Degrees of freedom in a higher-order ngram model
            s = g1.vcount()

            ## The degrees of freedom of the higher-order model, under the ngram assumption
            self.dof_ngrams = (s**k)*(s-1)
                 
            # For k>0, the degrees of freedom of a path-based model depend on 
            # the number of possible paths of length k in the first-order network.
            # Since probabilities in each row must sum to one, the degrees 
            # of freedom must be reduced by one for each k-order node 
            # that has at least one possible transition.

            # (A**k).sum() counts the number of different paths of exactly length k
            # based on the first-order network, which corresponds to the number of 
            # possible transitions in the transition matrix of a k-th order model. 
            paths_k = (A**k).sum()

            # For the degrees of freedom, we must additionally consider that 
            # rows in the transition matrix must sum to one, i.e. we have to 
            # subtract one degree of freedom for every non-zero row in the (null-model)
            # transition matrix. In other words, we subtract one for every path of length k-1 
            # that can possibly be followed by at least one edge to a path of length k

            # This can be calculated by counting the number of non-zero elements in the 
            # vector containing the row sums of A**k
            non_zero = _np.count_nonzero((A**k).sum(axis=0))

            ## The degrees of freedom of the higher-order model, under the paths assumption
            self.dof_paths = paths_k - non_zero


    def vcount(self):
        """ Returns the number of nodes """
        return len(self.nodes)


    def ecount(self):
        """ Returns the number of links """
        return len(self.edges)


    def totalEdgeWeight(self):
        """ Returns the sum of all edge weights """
        if len(self.edges)>0:
            return sum(self.edges.values())
        else:
            return _np.array([0,0])


    def HigherOrderNodeToPath(self, node):
        """
        Helper function that transforms a node in a
        higher-order network of order k into a corresponding 
        path of length k-1. For a higher-order node 'a-b-c-d' 
        this function will return ('a','b','c','d')

        @param node: The higher-order node to be transformed to a path.
        """

        return tuple(node.split(self.separator))


    def pathToHigherOrderNodes(self, path, k=None):
        """
        Helper function that transforms a path into a sequence of k-order nodes 
        using the separator character of the HigherOrderNetwork instance 

        Consider an example path (a,b,c,d) with a separator string '-'
        For k=1, the output will be the list of strings ['a', 'b', 'c', 'd']
        For k=2, the output will be the list of strings ['a-b', 'b-c', 'c-d']
        For k=3, the output will be the list of strings ['a-b-c', 'b-c-d']
        etc. 

        @param path: the path tuple to turn into a sequence of higher-order nodes 

        @param k: the order of the representation to use (default: order of the HigherOrderNetwork instance)
        """
        if k == None:
            k = self.order
        assert len(path)>k, 'Error: Path must be longer than k'

        if k == 0 and len(path)==1:
            return ['start', path[0]]

        nodes = []
    
        for s in range(0, len(path)-k+1):
            if s>len(path)-1 or s<0:
                print(path)
                print(s)
                print(k)
            v = path[s]
            for l in range(1,k):
                v = v + self.separator + path[s+l]
            nodes.append(v)

        return nodes


    def getNodeNameMap(self):
        """
        Returns a dictionary that can be used to map 
        node nodes to matrix/vector indices
        """

        name_map = {}
        for idx,v in enumerate(self.nodes):
            name_map[v] = idx
        return name_map


    def getDoF(self, assumption="paths"):
        """
        Calculates the degrees of freedom (i.e. number of parameters) of 
        this k-order model. Depending on the modeling assumptions, this either
        corresponds to the number of paths of length k in the first-order network 
        or to the number of all possible k-grams. The degrees of freedom of a model 
        can be used to assess the model complexity when calculating, e.g., the 
        Bayesian Information Criterion (BIC).

        @param assumption: if set to 'paths', for the degree of freedon calculation in the BIC, 
            only paths in the first-order network topology will be considered. This is 
            needed whenever we are interested in a modeling of paths in a given network topology.
            If set to 'ngrams' all possible n-grams will be considered, independent of whether they 
            are valid paths in the first-order network or not. The 'ngrams' and the 'paths' assumption 
            coincide if the first-order network is fully connected.
        """
        assert assumption == 'paths' or assumption == 'ngrams', 'Error: Invalid assumption'
        
        if assumption == 'paths':            
            return self.dof_paths
        else:            
            return self.dof_ngrams   


    def getDistanceMatrix(self):
        """
        Calculates shortest path distances between all pairs of 
        higher-order nodes using the Floyd-Warshall algorithm.
        """

        Log.add('Calculating distance matrix in higher-order network (k = ' + str(self.order) + ') ...', Severity.INFO)

        dist = _co.defaultdict( lambda: _co.defaultdict( lambda: _np.inf ) )

        for v in self.nodes:
            dist[v][v] = 0

        for e in self.edges:
            dist[e[0]][e[1]] = 1

        for v in self.nodes:
            for w in self.nodes:
                for k in self.nodes:
                    if dist[v][w] > dist[v][k] + dist[k][w]:
                        dist[v][w] = dist[v][k] + dist[k][w]

        Log.add('finished.', Severity.INFO)

        return dist


    def getShortestPaths(self):
        """
        Calculates all shortest paths between all pairs of 
        higher-order nodes using the Floyd-Warshall algorithm.
        """
        
        Log.add('Calculating shortest paths in higher-order network (k = ' + str(self.order) + ') ...', Severity.INFO)

        dist = _co.defaultdict( lambda: _co.defaultdict( lambda: _np.inf ) )
        shortest_paths = _co.defaultdict( lambda: _co.defaultdict( lambda: set() ) )

        for e in self.edges:
            dist[e[0]][e[1]] = 1
            shortest_paths[e[0]][e[1]].add(e)
       
        for v in self.nodes:
            for w in self.nodes:
                if v != w:
                    for k in self.nodes:
                        if dist[v][w] > dist[v][k] + dist[k][w]:
                            dist[v][w] = dist[v][k] + dist[k][w]
                            shortest_paths[v][w] = set()
                            for p in list(shortest_paths[v][k]):
                                for q in list(shortest_paths[k][w]):
                                    shortest_paths[v][w].add(p+q[1:])
                        elif dist[v][w] == dist[v][k] + dist[k][w]:
                            for p in list(shortest_paths[v][k]):
                                for q in list(shortest_paths[k][w]):
                                    shortest_paths[v][w].add(p+q[1:])
        
        for v in self.nodes:
            dist[v][v] = 0
            shortest_paths[v][v].add((v,))

        Log.add('finished.', Severity.INFO)

        return shortest_paths


    def getDistanceMatrixFirstOrder(self):
        """
        Projects a distance matrix from a higher-order to 
        first-order nodes, while path lengths are calculated 
        based on the higher-order topology
        """

        dist = self.getDistanceMatrix()
        dist_first = _co.defaultdict( lambda: _co.defaultdict( lambda: _np.inf ) )

        # calculate distances between first-order nodes based on distance in higher-order topology 
        for vk in dist:
            for wk in dist[vk]:
                v1 = self.HigherOrderNodeToPath(vk)[0]
                w1 = self.HigherOrderNodeToPath(wk)[-1]
                if dist[vk][wk] + self.order-1 < dist_first[v1][w1]:
                    dist_first[v1][w1] = dist[vk][wk] + self.order - 1

        return dist_first


    def ClosenessCentrality(self):
        """ 
        Calculates the closeness centralities of all nodes.
        If the order of the higher-order network is larger than one 
        centralities calculated based on the higher-order 
        topology will automatically be projected back to first-order 
        nodes.
        """

        dist_first = self.getDistanceMatrixFirstOrder()
        node_centralities = _co.defaultdict( lambda: 0 )   
        
        Log.add('Calculating closeness centralities (k = ' + str(self.order) + ') ...', Severity.INFO)             

        # calculate closeness values
        for v1 in dist_first:
            for w1 in dist_first[v1]:
                if v1 != w1 and dist_first[v1][w1] < _np.inf:
                    node_centralities[v1] += 1.0 / dist_first[v1][w1]
        
        # assign centrality zero to nodes not occurring on higher-order shortest paths
        nodes = self.paths.getNodes()
        for v in nodes:
            node_centralities[v] += 0

        Log.add('finished.', Severity.INFO)

        return node_centralities    


    def EvCent(self, projection='scaled', includeSubPaths=True):
        """
        Calculates the eigenvector centralities of higher-order nodes. If 
        the order of the HigherOrderNetwork is larger than one, the centralities
        will be projected to the first-order nodes. 

        @param projection: Indicates how the projection from k-th-order nodes (v1, v2, ... , v{k-1})
            shall be performed. For the method 'all', the eigenvector centrality of the higher-order node 
            will be added to *all* first-order nodes on the path corresponding to the higher-order node. For 
            the method 'last', the centrality of the higher-order node will only be assigned to *last* 
            first-order node v{k-1}. For the method 'scaled' (default), the eigenvector centrality of higher-order 
            nodes will be assigned proportionally to first-order nodes, i.e. each of the three nodes in the 
            third-order node (a,b,c) will receive one third of the eigenvector centrality of (a,b,c).
        @param includeSubPaths: whether or not to include subpath statistics in the calculation (default True)
        """
        A = self.getAdjacencyMatrix(includeSubPaths=includeSubPaths, weighted=False, transposed=True)

        # calculate leading eigenvector of A
        w, v = _sla.eigs(A, k=1, which="LM", ncv=13 )
        
        v = v.reshape(v.size,)
       
        higher_order_evcent = dict(zip(self.nodes, map(_np.abs, v)))

        # project evcent of higher-order nodes to first-order network
        first_order_evcent = _co.defaultdict( lambda: 0.0 )

        # sum evcent values based on higher-order nodes 
        # and normalize the result
        for v in self.nodes:
            # turns node a-b-c in path tuple (a,b,c)
            p = self.HigherOrderNodeToPath(v)
            if projection == 'all':
                # assign evcent of higher-order node to all first-order nodes
                for x in p:                    
                    first_order_evcent[x] += higher_order_evcent[v]
            elif projection == 'scaled':
                for x in p:                    
                    first_order_evcent[x] += higher_order_evcent[v] / float( len(p) )
            elif projection == 'last':
                # assign evcent of higher-order node to last first-order node
                first_order_evcent[p[-1]] += higher_order_evcent[v]
            elif projection == 'first':
                # assign evcent of higher-order node to last first-order node
                first_order_evcent[p[0]] += higher_order_evcent[v]

        # for scaled, values sum to one anyway
        if projection != 'scaled':
            for v in first_order_evcent:
                first_order_evcent[v] /= sum(first_order_evcent.values())

        Log.add('finished.', Severity.INFO)

        return first_order_evcent

        return v



    def PageRank(self, alpha=0.85, maxIterations=100, convergenceThres=1.0e-6, projection='scaled', includeSubPaths=True):
        """
        Calculates the PageRank of higher-order nodes based on a 
        power iteration. If the order of the higher-order network is larger than one,
        the PageRank calculated based on the higher-order
        topology will automatically be projected back to first-order 
        nodes.

        @param projection: Indicates how the projection from k-th-order nodes (v1, v2, ... , v{k-1})
            shall be performed. For the method 'all', the pagerank value of the higher-order node 
            will be added to *all* first-order nodes on the path corresponding to the higher-order node. For 
            the method 'last', the PR value of the higher-order node will only be assigned to *last* 
            first-order node v{k-1}. For the method 'scaled' (default), the PageRank of higher-order 
            nodes will be assigned proportionally to first-order nodes, i.e. each of the three nodes in the 
            third-order node (a,b,c) will receive one third of the PageRank of (a,b,c).
        @param includeSubpaths: whether or not to use subpath statistics in the PageRank calculation
        """

        assert projection == 'all' or projection == 'last' or projection == 'first' or projection == 'scaled', 'Invalid projection method'

        Log.add('Calculating PageRank in ' + str(self.order) + '-th order network...', Severity.INFO)

        higher_order_PR = _co.defaultdict( lambda: 0 )

        n = float(len(self.nodes))

        assert n>0, "Number of nodes is zero"

        # entries A[s,t] give directed link s -> t
        A = self.getAdjacencyMatrix(includeSubPaths=includeSubPaths, weighted=False, transposed=False)

        # sum of outgoing node degrees
        row_sums = _sp.array(A.sum(axis=1)).flatten()
        
        row_sums[row_sums != 0] = 1.0 / row_sums[row_sums != 0]
        d = _sp.where(row_sums == 0)[0]

        Q = _sparse.spdiags(row_sums.T, 0, *A.shape, format='csr')
        A = Q * A

        p = _sp.array([1.0 / n] * int(n))
        pr = p
       
        # Power iteration
        for i in range(maxIterations):
            last = pr
            pr = alpha * (pr * A + sum(pr[d]) * p) + (1 - alpha) * p

            if _sp.absolute(pr - last).sum() < n * convergenceThres:
                higher_order_PR = dict(zip(self.nodes, map(float, pr)))
                break

        if self.order == 1:
            return higher_order_PR

        # project PageRank of higher-order nodes to first-order network
        first_order_PR = _co.defaultdict( lambda: 0.0 )

        # sum PageRank values based on higher-order nodes 
        # and normalize the result
        for v in self.nodes:
            # turns node a-b-c in path tuple (a,b,c)
            p = self.HigherOrderNodeToPath(v)
            if projection == 'all':
                # assign PR of higher-order node to all first-order nodes
                for x in p:                    
                    first_order_PR[x] += higher_order_PR[v]
            elif projection == 'scaled':
                for x in p:
                    # each node on e.g. a 4-th-order path a-b-c-d receives one fourth of the 
                    # PageRank value, to ensure that the resulting first-order PageRank sums 
                    # to one
                    first_order_PR[x] += higher_order_PR[v] / float( len(p) )
            elif projection == 'last':
                # assign PR of higher-order node to last first-order node
                first_order_PR[p[-1]] += higher_order_PR[v]
            elif projection == 'first':
                # assign PR of higher-order node to last first-order node
                first_order_PR[p[0]] += higher_order_PR[v]

        # for projection method 'scaled', the values sum to one anyway
        if projection != 'scaled':
            for v in first_order_PR:
                first_order_PR[v] /= sum(first_order_PR.values())

        # assign centrality zero to nodes not occurring in higher-order PR
        nodes = self.paths.getNodes()
        for v in nodes:
            first_order_PR[v] += 0

        Log.add('finished.', Severity.INFO)

        return first_order_PR


    def HigherOrderPathToFirstOrder(self, path):
        """
        Maps a path in the higher-order network 
        to a path in the first-order network. As an 
        example, the second-order path ('a-b', 'b-c', 'c-d')
        of length two is mapped to the first-order path ('a','b','c','d')
        of length four. In general, a path of length l in a network of 
        order k is mapped to a path of length l+k-1 in the first-order network. 

        @param path: The higher-order path that shall be mapped to the first-order network
        """
        p1 = self.HigherOrderNodeToPath(path[0])
        for x in path[1:]:
            p1 += (self.HigherOrderNodeToPath(x)[-1],)
        return p1


    def BetweennessCentrality(self, normalized=False):
        """ 
        Calculates the betweenness centralities of all nodes.
        If the order of the higher-order network is larger than one 
        centralities calculated based on the higher-order 
        topology will automatically be projected back to first-order 
        nodes.

        @param normalized: If set to True, betweenness centralities of 
            nodes will be scaled by the maximum value (default False)
        """

        shortest_paths = self.getShortestPaths()
        node_centralities = _co.defaultdict( lambda: 0 )

        shortest_paths_firstorder =  _co.defaultdict( lambda:  _co.defaultdict( lambda: set() ) )

        Log.add('Calculating betweenness centralities (k = ' + str(self.order) + ') ...', Severity.INFO)

        for sk in shortest_paths:
            for dk in shortest_paths:
                s1 = self.HigherOrderNodeToPath(sk)[0]
                d1 = self.HigherOrderNodeToPath(dk)[-1]
                
                # we consider a path in a k-th order network 
                # connecting first-order node s1 to d1
                for pk in shortest_paths[sk][dk]:
                     # convert k-th order path to first-order path and add
                    shortest_paths_firstorder[s1][d1].add(self.HigherOrderPathToFirstOrder(pk))


        for s1 in shortest_paths_firstorder:
            for d1 in shortest_paths_firstorder[s1]:                               
                for p1 in shortest_paths_firstorder[s1][d1]:
                    # increase betweenness centrality of all intermediary nodes
                    # on path from s1 to d1
                    for v in p1[1:-1]:
                        if s1 != v != d1:
                            #print('node ' + x + ': ' + str(1.0 / len(shortest_paths[vk][wk])))
                            node_centralities[v] += 1.0 / (len(shortest_paths_firstorder[s1][d1]) + self.order-1)
                            #else:
                            #    node_centralities[v] += 1.0
        if normalized:
            m = max(node_centralities.values())
            for v in node_centralities:
                node_centralities[v] /= m

        # assign centrality zero to nodes not occurring on higher-order shortest paths
        nodes = self.paths.getNodes()
        for v in nodes:
            node_centralities[v] += 0

        Log.add('finished.', Severity.INFO)

        return node_centralities


    def reduceToGCC(self):
        """ 
        Reduces the higher-order network to its 
        largest (giant) strongly connected component 
        (using Tarjan's algorithm)
        """

        # nonlocal variables (!)        
        index = 0
        S = []
        indices = _co.defaultdict( lambda: None )
        lowlink = _co.defaultdict( lambda: None )        
        onstack = _co.defaultdict( lambda: False )

        # Tarjan's algorithm
        def strong_connect(v):
            nonlocal index
            nonlocal S
            nonlocal indices
            nonlocal lowlink            
            nonlocal onstack
                        
            indices[v] = index
            lowlink[v] = index
            index += 1
            S.append(v)
            onstack[v] = True

            for w in self.successors[v]:
                if indices[w] == None:
                    strong_connect(w)
                    lowlink[v] = min(lowlink[v], lowlink[w])
                elif onstack[w]:
                    lowlink[v] = min(lowlink[v], indices[w])
            
            # Generate SCC of node v
            component = set()
            if lowlink[v] == indices[v]:
                while True:
                    w = S.pop()
                    onstack[w] = False
                    component.add(w)
                    if v==w:
                        break
            return component

        # Get largest strongly connected component 
        components = _co.defaultdict( lambda: set() )
        max_size = 0
        max_head = None
        for v in self.nodes:
            if indices[v] == None:
                components[v] = strong_connect(v)
                if len(components[v]) > max_size:
                    max_head = v
                    max_size = len(components[v])
        
        scc = components[max_head]

        # Reduce higher-order network to SCC
        for v in list(self.nodes):
            if v not in scc:
                self.nodes.remove(v)
                del self.successors[v]

        for (v,w) in list(self.edges):
            if v not in scc or w not in scc:
                del self.edges[(v,w)]
                    

    def summary(self):
        """ 
        Returns a string containing basic summary statistics 
        of this higher-order graphical model instance
        """

        summary = 'Graphical model of order k = ' + str(self.order)
        summary += '\n'
        summary += 'Nodes:\t\t\t\t' +  str(self.vcount()) + '\n'
        summary += 'Links:\t\t\t\t' + str(self.ecount()) + '\n'
        summary += 'Total weight (sub/longest):\t' + str(self.totalEdgeWeight()[0]) + '/' + str(self.totalEdgeWeight()[1]) + '\n'
        return summary


    def __str__(self):
        """
        Returns the default string representation of 
        this graphical model instance
        """
        return self.summary()


    def degrees(self, includeSubPaths=True, weighted=True, mode="OUT"):
        """
        Returns the (weighted) degrees of nodes in the higher-order network

        @param weighted: If true, calculates the sum of weights for each node. If false, the 
            number of links is calculated

        @param mode: either "IN", "OUT", or "TOTAL" 
        """
        degrees = [0]*self.vcount()
        for v in self.nodes:
            for u,w in self.edges:
                if (mode == "OUT" and u == v) or (mode == "IN" and w == v) or (mode == "TOTAL" and (u==v or w==v)):
                    if weighted:
                        if includeSubPaths:
                            degrees[self.nodes.index(v)] += self.edges[(u,w)].sum()
                        else: 
                            degrees[self.nodes.index(v)] += self.edges[(u,w)][1]
                    else:
                        if includeSubPaths:
                            degrees[self.nodes.index(v)] += 1
                        else: 
                            if self.edges[(u,w)][1]>0:
                                degrees[self.nodes.index(v)] += 1
        return degrees


    def getAdjacencyMatrix(self, includeSubPaths = True, weighted = True, transposed=False):
        """
        Returns a sparse adjacency matrix of the higher-order network. By default, the entry 
            corresponding to a directed link source -> target is stored in row s and column t
            and can be accessed via A[s,t].
    
        @param includeSubPaths: if set to True, the returned adjacency matrix will 
            account for the occurrence of links of order k (i.e. paths of length k-1)
            as subpaths

        @param weighted: if set to False, the function returns a binary adjacency matrix.
          If set to True, adjacency matrix entries will contain the weight of an edge.
      
        @param transposed: whether to transpose the matrix or not.
        """        
    
        row = []
        col = []
        data = []
           
        if transposed:
            for s,t in self.edges:
                row.append(self.nodes.index(t))
                col.append(self.nodes.index(s))
        else:
            for s,t in self.edges:
                row.append(self.nodes.index(s))
                col.append(self.nodes.index(t))

        # create array with non-zero entries
        if not weighted:
            data = _np.ones(len(self.edges.keys()))
        else:
            if includeSubPaths:
                data = _np.array([float(x.sum()) for x in self.edges.values()])
            else:
                data = _np.array([float(x[1]) for x in self.edges.values()])

        return _sparse.coo_matrix( (data, (row, col)), shape=(self.vcount(), self.vcount()) ).tocsr()



    def getTransitionMatrix(self, includeSubPaths=True):
        """
        Returns a (transposed) random walk transition matrix 
        corresponding to the higher-order network.

        @param includeSubpaths: whether or not to include subpath statistics in the 
            transition probability calculation (default True)
        """
        row = []
        col = []
        data = []
        D = self.degrees(includeSubPaths=includeSubPaths, weighted=True, mode='OUT')
        for (s,t) in self.edges:
            # either s->t has been observed as a longest path, or we are interested in subpaths as well

            # the following makes sure that we do not accidentially consider zero-weight edges (automatically added by default_dic)
            if (self.edges[(s,t)][1] > 0) or (includeSubPaths and self.edges[(s,t)][0]>0):
                row.append(self.nodes.index(t))
                col.append(self.nodes.index(s))
                if includeSubPaths:
                    count = self.edges[(s,t)].sum()
                else:
                    count = self.edges[(s,t)][1]
                #print((s,t))
                #print(D[self.nodes.index(s)])
                #print(count)
                assert D[self.nodes.index(s)]>0, 'Encountered zero out-degree for node ' + str(s) + ' while weight of link (' + str(s) +  ', ' + str(t) + ') is non-zero.'
                prob = count / D[self.nodes.index(s)]
                if prob < 0 or prob > 1:
                    tn.Log.add('Encountered transition probability outside [0,1] range.', Severity.ERROR)
                    raise ValueError()
                data.append( prob )
    
        data = _np.array(data)
        data = data.reshape(data.size,)

        return _sparse.coo_matrix( (data, (row, col)), shape=(self.vcount(), self.vcount()) ).tocsr()


    @staticmethod
    def getLeadingEigenvector(A, normalized=True, lanczosVecs = 15, maxiter = 1000):
        """Compute normalized leading eigenvector of a given matrix A.

        @param A: sparse matrix for which leading eigenvector will be computed
        @param normalized: wheter or not to normalize. Default is C{True}
        @param lanczosVecs: number of Lanczos vectors to be used in the approximate
            calculation of eigenvectors and eigenvalues. This maps to the ncv parameter 
            of scipy's underlying function eigs. 
        @param maxiter: scaling factor for the number of iterations to be used in the 
            approximate calculation of eigenvectors and eigenvalues. The number of iterations 
            passed to scipy's underlying eigs function will be n*maxiter where n is the 
            number of rows/columns of the Laplacian matrix.
        """

        if _sparse.issparse(A) == False:
            raise TypeError("A must be a sparse matrix")

        # NOTE: ncv sets additional auxiliary eigenvectors that are computed
        # NOTE: in order to be more confident to find the one with the largest
        # NOTE: magnitude, see https://github.com/scipy/scipy/issues/4987
        w, pi = _sla.eigs( A, k=1, which="LM", ncv=lanczosVecs, maxiter=maxiter)
        pi = pi.reshape(pi.size,)
        if normalized:
            pi /= sum(pi)
        return pi


    def getLaplacianMatrix(self, includeSubPaths=True):
        """
        Returns the transposed Laplacian matrix corresponding to the higher-order network.

        @param includeSubpaths: Whether or not subpath statistics shall be included in the 
            calculation of matrix weights
        """   
    
        T = self.getTransitionMatrix(includeSubPaths)
        I  = _sparse.identity( self.vcount() )

        return I-T


    def getEigenValueGap(self, includeSubPaths=True, lanczosVecs = 15, maxiter = 20):
        """
        Returns the eigenvalue gap of the transition matrix.

        @param includeSubPaths: whether or not to include subpath statistics in the 
            calculation of transition probabilities.
        """
    
        #NOTE to myself: most of the time goes for construction of the 2nd order
        #NOTE            null graph, then for the 2nd order null transition matrix   
    
        Log.add('Calculating eigenvalue gap ... ', Severity.INFO)

        # Build transition matrices
        T = self.getTransitionMatrix(includeSubPaths)
    
        # Compute the two largest eigenvalues
        # NOTE: ncv sets additional auxiliary eigenvectors that are computed
        # NOTE: in order to be more confident to actually find the one with the largest
        # NOTE: magnitude, see https://github.com/scipy/scipy/issues/4987
        w2 = _sla.eigs(T, which="LM", k=2, ncv=lanczosVecs, return_eigenvectors=False, maxiter = maxiter)
        evals2_sorted = _np.sort(-_np.absolute(w2))
        
        Log.add('finished.', Severity.INFO)
    
        return _np.abs(evals2_sorted[1])    


    def getFiedlerVectorSparse(self, normalized = True, lanczosVecs = 15, maxiter = 20):
        """Returns the (sparse) Fiedler vector of the higher-order network. The Fiedler 
        vector can be used for a spectral bisectioning of the network.
     
        Note that sparse linear algebra for eigenvalue problems with small eigenvalues 
        is problematic in terms of numerical stability. Consider using the dense version
        of this method in this case. Note also that the sparse Fiedler vector might be scaled by 
        a factor (-1) compared to the dense version.
          
        @param normalized: whether (default) or not to normalize the fiedler vector.
          Normalization is done such that the sum of squares equals one in order to
          get reasonable values as entries might be positive and negative.
        @param lanczosVecs: number of Lanczos vectors to be used in the approximate
            calculation of eigenvectors and eigenvalues. This maps to the ncv parameter 
            of scipy's underlying function eigs. 
        @param maxiter: scaling factor for the number of iterations to be used in the 
            approximate calculation of eigenvectors and eigenvalues. The number of iterations 
            passed to scipy's underlying eigs function will be n*maxiter where n is the 
            number of rows/columns of the Laplacian matrix.
        """
    
        # NOTE: The transposed matrix is needed to get the "left" eigenvectors
        L = self.getLaplacianMatrix()

        # NOTE: ncv sets additional auxiliary eigenvectors that are computed
        # NOTE: in order to be more confident to find the one with the largest
        # NOTE: magnitude, see https://github.com/scipy/scipy/issues/4987
        maxiter = maxiter*L.get_shape()[0]
        w = _sla.eigs( L, k=2, which="SM", ncv=lanczosVecs, return_eigenvectors=False, maxiter=maxiter )
    
        # compute a sparse LU decomposition and solve for the eigenvector 
        # corresponding to the second largest eigenvalue
        n = L.get_shape()[0]
        b = _np.ones(n)
        evalue = _np.sort(_np.abs(w))[1]
        A = (L[1:n,:].tocsc()[:,1:n] - _sparse.identity(n-1).multiply(evalue))
        b[1:n] = A[0,:].toarray()
    
        lu = _sla.splu(A)
        b[1:n] = lu.solve(b[1:n])

        if normalized:
            b /= np.sqrt(np.inner(b, b))
        return b


    def getFiedlerVectorDense(self):
        """
         Returns the (dense)Fiedler vector of the higher-order network. The Fiedler 
         vector can be used for a spectral bisectioning of the network.             
        """
    
        # NOTE: The Laplacian is transposed for the sparse case to get the left
        # NOTE: eigenvalue.
        L = self.getLaplacianMatrix()
        # convert to dense matrix and transpose again to have the untransposed
        # laplacian again.
        w, v = _la.eig(L.todense().transpose(), right=False, left=True)

        return v[:,_np.argsort(_np.absolute(w))][:,1]



    def getAlgebraicConnectivity(self, lanczosVecs = 15, maxiter = 20):
        """
        Returns the algebraic connectivity of the higher-order network.    
        
        @param lanczosVecs: number of Lanczos vectors to be used in the approximate
            calculation of eigenvectors and eigenvalues. This maps to the ncv parameter 
            of scipy's underlying function eigs. 
        @param maxiter: scaling factor for the number of iterations to be used in the 
            approximate calculation of eigenvectors and eigenvalues. The number of iterations 
            passed to scipy's underlying eigs function will be n*maxiter where n is the
            number of rows/columns of the Laplacian matrix.         
        """
    
        Log.add('Calculating algebraic connectivity ... ', Severity.INFO)

        L = self.getLaplacianMatrix()
        # NOTE: ncv sets additional auxiliary eigenvectors that are computed
        # NOTE: in order to be more confident to find the one with the largest
        # NOTE: magnitude, see https://github.com/scipy/scipy/issues/4987
        w = _sla.eigs( L, which="SM", k=2, ncv=lanczosVecs, return_eigenvectors=False, maxiter = maxiter )
        evals_sorted = _np.sort(_np.absolute(w))

        Log.add('finished.', Severity.INFO)

        return _np.abs(evals_sorted[1])
