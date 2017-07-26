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
import collections as _co
import bisect as _bs
import itertools as _iter

import scipy.sparse as _sparse
import scipy.misc as _misc
import scipy.sparse.linalg as _sla
import scipy.linalg as _la

from scipy.stats import chi2

from pathpy.Log import Log
from pathpy.Log import Severity

_np.seterr(all='warn')


class MarkovSequence:
    """ Instances of this class can be used to fit
        standard higher-order Markov models for
        sequences generated from concatenated paths """

    def __init__(self, sequence):
        """
        Generates a Markov model for a sequence, given
        as a single list of strings
        """

        ## The sequence to be modeled
        self.sequence = sequence

        ## The transition probabilities of higher-order Markov chains
        self.P = {}

        ## the set of states of higher-order Markov chains
        self.states = {}
        self.states[1] = set(sequence)


    def fitMarkovModel(self, k=1):
        """ Generates a k-th order Markov model
            for the underlying sequence
        """

        # TODO: Add support for k=0

        assert len(self.sequence)>0, "Error: Empty sequence"

        # MLE fit of transition probabilities
        self.P[k] = _co.defaultdict( lambda:  _co.defaultdict( lambda: 0.0 )  )

        Log.add('Fitting Markov model with order k = ' + str(k))

        # Generate initial memory prefix
        mem = (())
        for s in self.sequence[:k]:
            mem += (s,)

        # count state transitions
        for s in self.sequence[k:]:
            self.P[k][mem][s] += 1.0

            # shift memory by one element
            mem = mem[1:] + (s,)

        # normalize transitions
        for m in self.P[k]:
            S = float(sum(self.P[k][m].values()))
            for s in self.P[k][m]:
                self.P[k][m][s] /= S
        Log.add('finished.')


    def getLikelihood(self, k=1, log=True):
        """
        Returns the likelihood of the sequence
        assuming a k-th order Markov model
        """

        if k not in self.P:
            self.fitMarkovModel(k)

        L = 0

         # Generate initial prefix
        mem = (())
        for s in self.sequence[:k]:
            mem += (s,)

        for s in self.sequence[k:]:
            L += _np.log(self.P[k][mem][s])

            # shift memory by one element
            mem = mem[1:] + (s,)

        if log:
            return L
        else:
            return _np.exp(L)


    def getBIC(self, k=1, m=1):
        """ Returns the Bayesian Information Criterion
            assuming a k-th order Markov model """

        if k not in self.P:
            self.fitMarkovModel(k)

        if m not in self.P:
            self.fitMarkovModel(m)

        L_k = self.getLikelihood(k, log=True)
        L_m = self.getLikelihood(m, log=True)

        s = len(self.states[1])
        n = len(self.sequence)-k

        # the transition matrix of a first-order model with s states has s**2 entries, subject to the
        # constraint that entries in each row must sum up to one (thus effectively reducing
        # the degrees of freedom by a factor of s, i.e. we have s**2-s**1. Generalizing this to order k,
        # we arrive at s**k * (s-1) = s**(k+1) - s**k derees of freedom
        bic = _np.log(n) * (s**k - s**m) * (s-1) - 2.0 * (L_k-L_m)

        return bic


    def getAIC(self, k=1, m=1):
        """ Returns the Aikake Information Criterion
            assuming a k-th order Markov model """

        if k not in self.P:
            self.fitMarkovModel(k)

        if m not in self.P:
            self.fitMarkovModel(m)

        L_k = self.getLikelihood(k, log=True)
        L_m = self.getLikelihood(m, log=True)

        s = len(self.states[1])
        n = len(self.sequence)

        aic = 2 * (s**k - s**m) * (s-1) - 2.0 * (L_k - L_m)

        return aic


    def estimateOrder(self, maxOrder, method='BIC'):
        """ Estimates the optimal order of a Markov model
            based on Likelihood, BIC or AIC """

        assert method == 'BIC' or method == 'AIC' or method == 'Likelihood', "Error: Expecting method 'AIC', 'BIC' or 'Likelihood'"

        values = []
        orders = []

        # We need k < m for the BIC and AIC calculation, which
        # is why we only test up to maxOrder - 1
        for k in range(1, maxOrder):
            if k not in self.P:
                self.fitMarkovModel(k)

            orders.append(k)

            if method == 'AIC':
                values.append(self.getAIC(k, maxOrder))
            elif method == 'BIC':
                values.append(self.getBIC(k, maxOrder))
            elif method == 'Likelihood':
                values.append(self.getLikelihood(k, log=True))

        if method == 'Likelihood':
            values.append(self.getLikelihood(maxOrder, log=True))
            orders.append(maxOrder)

            # return order at which likelihood is maximized
            return orders[_np.argmax(values)]
        else:
            # return order at which BIC/AIC are minimized
            return orders[_np.argmin(values)]
