# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 15:41:27 2022

@author: danie
"""

import bgheatmaps as bgh


values = dict(  # scalar values for each region
   CTX=2608,
    CNU=280,
    IB=1636,
    MB=0,
    VL=0,
    cc =408,
    HIP= 1000
    
)


f = bgh.heatmap(
    values,
    position=(8000),  # when using a named orientation you can pass a single value!
    orientation="frontal",  # 'frontal' or 'sagittal', or 'horizontal' or a tuple (x,y,z)
    title="",
    vmin=0,
    vmax=10000,
    atlas_name = "allen_mouse_25um",
    cmap='jet',
    format="2D",
)

f.show()



