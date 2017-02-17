# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 11:49:39 2015
@author: Ingo Scholtes, Roman Cattaneo

(c) Copyright ETH Zürich, Chair of Systems Design, 2015-2017
"""
import numpy as _np
import sys as _sys
import collections as _co
import bisect as _bs
import datetime as _dt

from pathpy.Log import Log
import pathpy.Paths


class TemporalNetwork:
    """ This class represents a sequence of time-stamped edges.
       Instances of this class can be used to generate path statistics 
       based on the time-respecting paths resulting from a given maximum
       time difference between consecutive time-stamped edges.
    """

    def __init__(self, tedges = None):
        """
        Constructor generating a temporal network instance
        
        @param tedges: an optional list of (possibly unordered time-stamped) links 
            from which to construct a temporal network instance        
        """

        self.tedges = []
        nodes_seen = _co.defaultdict( lambda:False )
        self.nodes = []

        # Generate index structures which help to efficiently extract time-respecting paths

        # A dictionary storing all time-stamped links, indexed by time-stamps
        self.time = _co.defaultdict( lambda: list() )

        # A dictionary storing all time-stamped links, indexed by time and target node
        self.targets = _co.defaultdict( lambda: dict() )

        # A dictionary storing all time-stamped links, indexed by time and source node 
        self.sources = _co.defaultdict( lambda: dict() )

        # A dictionary storing time stamps at which links (v,*;t) originate from node v
        self.activities = _co.defaultdict( lambda: list() )

        # A dictionary storing sets of time stamps at which links (v,*;t) originate from node v
        # Note that the insertion into a set is much faster than repeatedly checking whether 
        # an element already exists in a list!
        self.activities_sets = _co.defaultdict( lambda: set() )

        # An ordered list of time-stamps
        self.ordered_times = []

        self.tedges = []

        if tedges is not None:
            Log.add('Building index data structures ...')

            for e in tedges:
                self.activities_sets[e[0]].add(e[2])
                self.time[e[2]].append(e)
                self.targets[e[2]].setdefault(e[1], []).append(e)
                self.sources[e[2]].setdefault(e[0], []).append(e)
                if not nodes_seen[e[0]]:
                    nodes_seen[e[0]] = True
                if not nodes_seen[e[1]]:
                    nodes_seen[e[1]] = True
            self.tedges = tedges
            self.nodes = list(nodes_seen.keys())

            Log.add('Sorting time stamps ...')

            self.ordered_times = sorted(self.time.keys())
            for v in self.nodes:
                self.activities[v] = sorted(self.activities_sets[v])
            Log.add('finished.')



    def readFile(filename='', sep=',', timestampformat="%s", maxlines=_sys.maxsize):
        """ Reads time-stamped links from a file and returns a new instance 
            of the class TemporalNetwork
        """
        assert (filename != ''), 'Empty filename given'
        
        # Read header
        with open(filename, 'r') as f:
            tedges = []
            twopaths = []
        
            header = f.readline()
            header = header.split(sep)

            # If header columns are included, arbitrary column orders are supported
            time_ix = -1
            source_ix = -1
            mid_ix = -1
            weight_ix = -1
            target_ix = -1
            for i in range(len(header)):
                header[i] = header[i].strip()
                if header[i] == 'node1' or header[i] == 'source':
                    source_ix = i
                elif header[i] == 'node2' or header[i] == 'target':
                    target_ix = i
                elif header[i] == 'time' or header[i] == 'timestamp':
                    time_ix = i

            assert (source_ix >= 0 and target_ix >= 0), "Detected invalid header columns: %s" % header

            if time_ix<0:
                Log.add('No time stamps found in data, assuming consecutive links', Severity.WARNING)
        
            Log.add('Reading time-stamped links ...')

            line = f.readline()
            n = 1 
            while line and n <= maxlines:
                fields = line.rstrip().split(sep)
                try:
                    if time_ix >=0:
                        timestamp = fields[time_ix]            
                        if timestamp.isdigit():
                            t = int(timestamp)
                        else:
                            x = _dt.datetime.strptime(timestamp, "%Y-%m-%d %H:%M")
                            t = int(time.mktime(x.timetuple()))
                    else:
                        t = n                
                    if t>=0:
                        tedge = (fields[source_ix], fields[target_ix], t)
                        tedges.append(tedge)
                    else:
                        Log.add('Ignoring negative timestamp in line ' + str(n+1) + ': "' + line.strip() + '"', Severity.WARNING)
                except (IndexError, ValueError):
                    Log.add('Ignoring malformed data in line ' + str(n+1) + ': "' +  line.strip() + '"', Severity.WARNING)
                line = f.readline()
                n += 1
        # end of with open()

        return TemporalNetwork(tedges = tedges)



    def filterEdges(self, edge_filter):
        """Filter time-stamped edges according to a given filter expression. 

        @param edge_filter: an arbitrary filter function of the form filter_func(v, w, time) that 
            returns True for time-stamped edges that shall pass the filter, and False for all edges that 
            shall be filtered out.
        """

        Log.add('Starting filtering ...', Severity.INFO)
        new_t_edges = []

        for (v,w,t) in self.tedges:
            if edge_filter(v,w,t):
                new_t_edges.append((v,w,t))

        Log.add('finished. Filtered out ' + str(self.ecount() - len(new_t_edges)) + ' time-stamped edges.', Severity.INFO)

        return TemporalNetwork(tedges=new_t_edges)


    def addEdge(self, source, target, ts):
        """Adds a directed time-stamped edge (source,target;time) to the temporal network. To add an undirected 
            time-stamped link (u,v;t) at time t, please call addEdge(u,v;t) and addEdge(v,u;t).
        
        @param source: naem of the source node of a directed, time-stamped link
        @param target: name of the target node of a directed, time-stamped link
        @param ts: (integer) time-stamp of the time-stamped link
        """
        e = (source, target, ts)
        self.tedges.append(e)
        if source not in self.nodes:
            self.nodes.append(source)
        if target not in self.nodes:
            self.nodes.append(target)

        # Add edge to index structures
        self.time[ts].append(e)
        self.targets[ts].setdefault(target, []).append(e)
        self.sources[ts].setdefault(source, []).append(e)

        if ts not in self.activities[source]:
            self.activities[source].append(ts)
            self.activities[source].sort()

        # Reorder time stamps
        self.ordered_times = sorted(self.time.keys())       


    def vcount(self):
        """
        Returns the total number of different vertices active across the whole evolution of the 
        temporal network. This number corresponds to the number of nodes in the (first-order) 
        time-aggregated network.
        """

        return len(self.nodes)

        
    def ecount(self):
        """
        Returns the number of time-stamped edges (u,v;t)
        """

        return len(self.tedges)


    def getObservationLength(self):
        """
        Returns the length of the observation time.
        """

        return max(self.ordered_times)-min(self.ordered_times)
    

    def getInterEventTimes(self):
        """
        Returns a numpy array containing all time differences between any 
        two consecutive time-stamped links (involving any node)
        """

        timediffs = []
        for i in range(1, len(self.ordered_times)):
            timediffs += [self.ordered_times[i] - self.ordered_times[i-1]]
        return _np.array(timediffs)


    def getInterPathTimes(self):
        """
        Returns a dictionary which, for each node v, contains all time differences 
        between any time-stamped link (*,v;t) and the next link (v,*;t') (t'>t)
        in the temporal network
        """

        interPathTimes = defaultdict( lambda: list() )
        for e in self.tedges:
            # Get target v of current edge e=(u,v,t)
            v = e[1]
            t = e[2]

            # Get time stamp of link (v,*,t_next) with smallest t_next such that t_next > t
            i = _bs.bisect_right(self.activities[v], t)
            if i != len(self.activities[v]):
                interPathTimes[v].append(self.activities[v][i]-t)
        return interPathTimes


    def summary(self):
        """
        Returns a string containing basic summary statistics of this temporal network
        """

        summary = ''

        summary += 'Nodes:\t\t\t' +  str(self.vcount()) + '\n'
        summary += 'Time-stamped links:\t' + str(self.ecount()) + '\n'
        if self.vcount()>0:
            summary += 'Links/Nodes:\t\t' + str(self.ecount()/self.vcount()) + '\n'
        else:
            summary += 'Links/Nodes:\t\tN/A\n'
        if len(self.ordered_times)>1:
            summary += 'Observation period:\t[' + str(min(self.ordered_times)) + ', ' + str(max(self.ordered_times)) + ']\n'
            summary += 'Observation length:\t' + str(max(self.ordered_times) - min(self.ordered_times)) + '\n'
            summary += 'Time stamps:\t\t' + str(len(self.ordered_times)) + '\n'

            d = self.getInterEventTimes()    
            summary += 'Avg. inter-event dt:\t' + str(_np.mean(d)) + '\n'
            summary += 'Min/Max inter-event dt:\t' + str(min(d)) + '/' + str(max(d)) + '\n'
        
        return summary       


    def __str__(self):
        """
        Returns the default string representation of 
        this temporal network instance
        """
        return self.summary()
   

    def ShuffleEdges(self, l=0, with_replacement=False):        
        """
        Generates a shuffled version of the temporal network in which edge statistics (i.e.
        the frequencies of time-stamped edges) are preserved, while all order correlations are 
        destroyed. The shuffling procedure randomly reshuffles the time-stamps of links.
        
        @param l: the length of the sequence to be generated (i.e. the number of time-stamped links.
            For the default value l=0, the length of the generated shuffled temporal network will be 
            equal to that of the original temporal network. 
        """

        tedges = []        
        
        timestamps = [e[2] for e in self.tedges]
        edges = list(self.tedges)
        
        if l==0:
            l = len(self.tedges)
        for i in range(l):
            
            if with_replacement:
            # Pick random link
                edge = edges[_np.random.randint(0, len(edges))]
                # Pick random time stamp
                time = timestamps[_np.random.randint(0, len(timestamps))]
            else:
                # Pick random link
                edge = edges.pop(_np.random.randint(0, len(edges)))
            # Pick random time stamp
                time = timestamps.pop(_np.random.randint(0, len(timestamps)))            
            
            # Generate new time-stamped link
            tedges.append( (edge[0], edge[1], time) )

        # Generate temporal network
        t = TemporalNetwork(tedges=tedges)

        # Fix node order to correspond to original network
        t.nodes = self.nodes
            
        return t    


    def GetTemporalBetweenness(t, delta=1, normalized=False):
        """
        Calculates the temporal betweenness centralities of all nodes based on the shortest 
        time-respecting paths with a maximum waiting time of delta. This function returns a 
        numpy array of temporal betweenness centrality values of nodes.
    
        @param t: the temporal network for which temporal closeness centralities will be computed    
        @param delta: the maximum time difference used in the time-respecting path definition (default 1).
            Note that this parameter is independent from the delta used internally for the extraction 
            of two-paths by the class TemporalNetwork
        @param normalized: whether or not to normalize centralities by dividing each value by the total 
            number of shortest time-respecting paths.
        """

        bw = np.array([0]*len(t.nodes))
        S = 0

        minD, minPaths = Distances.GetMinTemporalDistance(t, delta=1, collect_paths=True)

        for v in t.nodes:
            for w in t.nodes:
                for p in minPaths[v][w]:
                    for i in range(1,len(p)-1):
                        bw[t.nodes.index(p[i][0])] += 1
                        S+=1
        return bw


    def GetTemporalBetweennessInstantaneous(t, start_t=0, delta=1, normalized=False):
        """
        Calculates the temporal betweennness values of 
        all nodes fir a given start time start_t in an empirical temporal network t.
        This function returns a numpy array of (temporal) betweenness centrality values. 
        The ordering of these values corresponds to the ordering of nodes in the vertex 
        sequence of the igraph first order time-aggregated network. A mapping between node names
        and array indices can be found in Utilities.firstOrderNameMap().
    
        @param t: the temporal network for which temporal betweenness centralities will be computed
        @param start_t: the start time for which to consider time-respecting paths (default 0). This is 
            important, since any unambigious definition of a shortest time-respecting path between
            two nodes must include the time range to be considered (c.f. Holme and Saramäki, Phys. Rep., 2012)
        @param delta: the maximum waiting time used in the time-respecting path definition (default 1)
            Note that this parameter is independent from the delta used internally for the extraction of two-paths
            by the class TemporalNetwork
        @param normalized: whether or not to normalize the temporal betweenness centrality values by
            dividing by the number of all shortest time-respecting paths in the temporal network.
        """

        bw = np.array([0]*len(t.nodes))

        # First calculate all shortest time-respecting paths starting at time start_t
        D, paths = Distances.GetTemporalDistanceMatrix(t, start_t, delta, collect_paths=True)

        # Compute betweenness scores of all nodes based on shortest time-respecting paths
        k=0
        for u in t.nodes:
            for v in t.nodes:
                if u != v:
                    for p in paths[u][v]:
                        for i in range(1, len(p)-1):
                            bw[t.nodes.index(p[i][0])] += 1
                            k+=1

        # Normalize by dividing by the total number of shortest time-respecting paths
        if normalized:
            bw = bw/k
        return bw


    def GetTemporalCloseness(t, delta=1):
        """
        Calculates the temporal closeness centralities of all nodes based on the minimal 
        shortest time-respecting paths with a maximum time difference of delta. This function 
        returns a numpy array of average (temporal) closeness centrality values of nodes.
    
        @param t: the temporal network for which temporal closeness centralities will be computed   
        @param delta: the maximum waiting time used in the time-respecting path definition (default 1)            
        """

        cl = np.array([0.]*len(t.nodes))

        minD, minPaths = Distances.GetMinTemporalDistance(t, delta, collect_paths=False)

        for u in t.nodes:
            for v in t.nodes:
                if u!= v:
                    cl[t.nodes.index(v)] += 1./minD[t.nodes.index(u), t.nodes.index(v)]
        return cl


    def GetTemporalClosenessInstantaneous(t, start_t=0, delta=1):
        """
        Calculates the temporal closeness values of all nodes for a given start time start_t.
        This function returns a numpy array of (temporal) closeness centrality values.        
    
        @param t: the temporal network for which temporal closeness centralities will be computed
        @param start_t: the start time for which to consider time-respecting paths (default 0). This is 
            important, since any unambigious definition of a shortest time-respecting path between
            two nodes must include the time range to be considered (c.f. Holme and Saramäki, Phys. Rep., 2012)
        @param delta: the maximum time difference time used in the time-respecting path definition (default 1)
            Note that this parameter is independent from the delta used internally for the extraction of two-paths
            by the class TemporalNetwork.
        """
    
        closeness = np.array([0.]*len(t.nodes))

        # Calculate all shortest time-respecting paths
        D, paths = Distances.GetTemporalDistanceMatrix(t, start_t, delta, collect_paths=False)

        # Get a mapping between node names and matrix indices
        name_map = Utilities.firstOrderNameMap( t )

        # Calculate closeness for each node u, by summing the reciprocal of its 
        # distances to all other nodes. Note that this definition of closeness centrality 
        # is required for directed networks that are not strongly connected. 
        for u in t.nodes:
            for v in t.nodes:
                if u!=v:
                    closeness[name_map[u]] += 1./D[name_map[v], name_map[u]]

        return closeness