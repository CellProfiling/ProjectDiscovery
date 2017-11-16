#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 13:27:15 2017

@author: devinsullivan
"""

#This is a basic function for saving the high resolution figure for the
#scatter plot of number of tasks vs hamming score (Fig 3)

import plotly.plotly as py

fig = py.get_figure("https://plot.ly/~dpsulliv/24/")
fig.layout.yaxis.range = [0,1.1]
fig.layout.xaxis.ticks = [1,2,3,4,5,6]
fig.layout.xaxis.range = [1,6.1]
fig.data.append(Scatter(x=(1,6),y=(0.5631,0.5631),mode='lines',line=Line(color='rgba(0,0,0)',width=1.0)))

#data = py.get_figure("https://plot.ly/~dpsulliv/24/").get_data()

py.image.save_as(fig, filename='tasks_vs_hamming.png',scale=3)