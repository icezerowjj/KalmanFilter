# -*- coding: utf-8 -*-
"""
Created on Sat Nov 26 08:54:41 2016

@author: zerow
"""

from pandas.io.data import DataReader
from numpy import *

if __name__ == "__main__":
    secs = ['EWA', 'EWC']
    data = DataReader(secs, 'yahoo', '2010-1-1', '2014-8-1')['Adj Close']
    obs_mat = vstack([data.EWA, ones(data.EWA.shape)]).T[:, newaxis]

