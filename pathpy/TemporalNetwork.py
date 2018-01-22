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

import numpy as _np
import sys as _sys
import collections as _co
import bisect as _bs
import datetime as _dt
import time as _t

from pathpy.Log import Log
from pathpy.Log import Severity
import pathpy.Paths


class TemporalNetwork:
    """ This class represents a sequence of time-stamped edges.
       Instances of this class can be used to generate path statistics 
       based on the time-respecting paths resulting from a given maximum
       time difference between consecutive time-stamped edges.
    """

    def __init__(self, tedges = None):
        """
        Constructor that generates a temporal network instance. 
        
        @param tedges: an optional list of (possibly unordered time-stamped) links 
            from which to construct a temporal network instance. For the default value None        
            an empty temporal network will be created.
        """

        ## A list of time-stamped edges of this temporal network
        self.tedges = []        

        ## A list of nodes of this temporal network
        self.nodes = []

        ## A dictionary storing all time-stamped links, indexed by time-stamps
        self.time = _co.defaultdict( lambda: list() )

        ## A dictionary storing all time-stamped links, indexed by time and target node
        self.targets = _co.defaultdict( lambda: dict() )

        ## A dictionary storing all time-stamped links, indexed by time and source node 
        self.sources = _co.defaultdict( lambda: dict() )

        ## A dictionary storing time stamps at which links (v,*;t) originate from node v
        self.activities = _co.defaultdict( lambda: list() )

        ## A dictionary storing sets of time stamps at which links (v,*;t) originate from node v
        ## Note that the insertion into a set is much faster than repeatedly checking whether 
        ## an element already exists in a list!
        self.activities_sets = _co.defaultdict( lambda: set() )

        ## An ordered list of time-stamps
        self.ordered_times = []        

        nodes_seen = _co.defaultdict( lambda:False )

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


    @staticmethod
    def readFile(filename, sep=',', timestampformat="%Y-%m-%d %H:%M", maxlines=_sys.maxsize):
        """ Reads time-stamped links from a file and returns a new instance 
            of the class TemporalNetwork. The file is assumed to have a header 

                source target time 

            where columns can be in arbitrary order and separated by arbitrary characters. 
            Each time-stamped link must occur in a separate line and links are assumed to be
            directed.
             
            The time column can be omitted and in this case all links are assumed to occur 
            in consecutive time stamps (that have a distance of one). Time stamps can be simple 
            integers, or strings to be converted to UNIX time stamps via a custom timestamp format. 
            For this, the python function datetime.strptime will be used. 

            @param sep: the character that separates columns 
            @param filename: path of the file to read from
            @param timestampformat: used to convert string timestamps to UNIX timestamps. This parameter is 
                ignored, if the timestamps are digit types (like a simple int).
            @param maxlines: limit reading of file to certain number of lines, default sys.maxsize

        """
        assert (filename != ''), 'Empty filename given'
        
        # Read header
        with open(filename, 'r') as f:
            tedges = []
        
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
                        # if the timestamp is a number, we use this 
                        if timestamp.isdigit():
                            t = int(timestamp)
                        else:   # if it is a string, we use the timestamp format to convert it to a UNIX timestamp                                
                            x = _dt.datetime.strptime(timestamp, timestampformat)
                            t = int(_t.mktime(x.timetuple()))
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
        
        @param source: name of the source node of a directed, time-stamped link
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
        Returns the number of vertices in the temporal network. 
        This number corresponds to the number of nodes in the (first-order) 
        time-aggregated network.
        """

        return len(self.nodes)

        
    def ecount(self):
        """
        Returns the number of time-stamped edges (u,v;t) in the temporal network.
        This number corresponds to the sum of link weights in the (first-order)
        time-aggregated network.
        """

        return len(self.tedges)


    def getObservationLength(self):
        """
        Returns the length of the observation time in time units.
        """

        return max(self.ordered_times)-min(self.ordered_times)
    

    def getInterEventTimes(self):
        """
        Returns an array containing all time differences between any 
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

        interPathTimes = _co.defaultdict( lambda: list() )
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
        this temporal network instance.
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
        @param with_replacement: Whether or not the sampling should be with replacement (default False)
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


    def exportUnfoldedNetwork(self, filename):
        """
        Generates a tex file that can be compiled to a time-unfolded 
        representation of the temporal network.

        @param filename: the name of the tex file to be generated.
        """    
        
        import os as _os

        output = []
            
        output.append('\\documentclass{article}\n')
        output.append('\\usepackage{tikz}\n')
        output.append('\\usepackage{verbatim}\n')
        output.append('\\usepackage[active,tightpage]{preview}\n')
        output.append('\\PreviewEnvironment{tikzpicture}\n')
        output.append('\\setlength\PreviewBorder{5pt}%\n')
        output.append('\\usetikzlibrary{arrows}\n')
        output.append('\\usetikzlibrary{positioning}\n')
        output.append('\\begin{document}\n')
        output.append('\\begin{center}\n')
        output.append('\\newcounter{a}\n')
        output.append("\\begin{tikzpicture}[->,>=stealth',auto,scale=0.5, every node/.style={scale=0.9}]\n")
        output.append("\\tikzstyle{node} = [fill=lightgray,text=black,circle]\n")
        output.append("\\tikzstyle{v} = [fill=black,text=white,circle]\n")
        output.append("\\tikzstyle{dst} = [fill=lightgray,text=black,circle]\n")
        output.append("\\tikzstyle{lbl} = [fill=white,text=black,circle]\n")

        last = ''
            
        for n in _np.sort(self.nodes):
            if last == '':
                output.append("\\node[lbl]                     (" + n + "-0)   {$" + n + "$};\n")
            else:
                output.append("\\node[lbl,right=0.5cm of "+last+"-0] (" + n + "-0)   {$" + n + "$};\n")
            last = n
            
        output.append("\\setcounter{a}{0}\n")
        output.append("\\foreach \\number in {"+ str(min(self.ordered_times))+ ",...," + str(max(self.ordered_times)+1) + "}{\n")
        output.append("\\setcounter{a}{\\number}\n")
        output.append("\\addtocounter{a}{-1}\n")
        output.append("\\pgfmathparse{\\thea}\n")
        
        for n in  _np.sort(self.nodes):
            output.append("\\node[v,below=0.3cm of " + n + "-\\pgfmathresult]     (" + n + "-\\number) {};\n")
        output.append("\\node[lbl,left=0.5cm of " + _np.sort(self.nodes)[0] + "-\\number]    (col-\\pgfmathresult) {$t=$\\number};\n")
        output.append("}\n")
        output.append("\\path[->,thick]\n")
        i = 1
        
        for ts in self.ordered_times:
            for edge in self.time[ts]:
                output.append("(" + edge[0] + "-" + str(ts) + ") edge (" + edge[1] + "-" + str(ts + 1) + ")\n")
                i += 1                                
        output.append(";\n")
        output.append(
    """\end{tikzpicture}
    \end{center}
    \end{document}""")
        
        # create directory if necessary to avoid IO errors
        directory = _os.path.dirname( filename )
        if directory != '':
            if not _os.path.exists( directory ):
                _os.makedirs( directory )

        with open(filename, "w") as tex_file:
            tex_file.write(''.join(output))