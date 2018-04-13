# -*- coding: utf-8 -*-
#
# Copyright (c) 2017, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Calculate properties of nuclei based on data parsed by cclib."""

import logging

import numpy

from cclib.method.calculationmethod import Method


class Nuclear(Method):
    """A container for methods pertaining to atomic nuclei."""

    def __init__(self, data, progress=None, loglevel=logging.INFO, logname="Log"):

        self.required_attrs = ('natom','atomcoords','atomnos')

        super(Nuclear, self).__init__(data, progress, loglevel, logname)

    def __str__(self):
        """Return a string representation of the object."""
        return "Nuclear"

    def __repr__(self):
        """Return a representation of the object."""
        return "Nuclear"

    def repulsion_energy(self):
        """Return the nuclear repulsion energy."""
        CstBohr2Ang = 0.52917721092 # !!!
        nre = 0.0
        for i in range(data.natom):
            ri = data.atomcoords[-1][i]
            zi = data.atomnos[i]
            for j in range(i + 1, data.natom):
                rj = data.atomcoords[-1][j]
                zj = data.atomnos[j]
                d = np.linalg.norm(ri - rj)
                nre += zi * zj / d
        return nre * CstBohr2Ang  # !!!


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
