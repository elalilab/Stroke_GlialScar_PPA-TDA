# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 13:24:51 2022

@author: danie
"""

# brainrender script for plotting cells from ClearMap 

# Author: 	Luke Hammond
# Cellular Imaging | Zuckerman Institute, Columbia University
# Date:	1st October, 2021

import numpy as np

import brainrender
from brainrender import Scene
from brainrender.actors import Points, PointsDensity
from brainrender.actors import Volume

from rich import print
from myterial import orange
from pathlib import Path
import pandas as pd

#if making videos
from brainrender import VideoMaker, Animation

#Modify Settings
brainrender.settings.WHOLE_SCREEN = False  # make the rendering window be smaller
brainrender.settings.SHOW_AXES = False
brainrender.settings.SHADER_STYLE = "plastic"  #[cartoon, metallic, plastic, shiny, glossy]
brainrender.settings.DEFAULT_ATLAS = "allen_mouse_25um" 
brainrender.settings.DEFAULT_CAMERA = "three_quarters"
    #[sagittal, sagittal2, frontal, top, top_side, three_quarters]
brainrender.settings.SCREENSHOT_SCALE = 1
brainrender.settings.ROOT_ALPHA = 0.2  # transparency of the overall brain model's actor'


print(f"[{orange}]Running example: {Path(__file__).name}")


# Provide path to ClearMap cells csv file that you wish to plot
NeuN = "D:/Research/Project_GlialTopology/3.DataAnalysis/Exp2-Gfap,NeuN,Iba1/Results/Coordinates/ECM_Exp2_M05_30D_NeuN_Coordinates.csv"
Gfap = "D:/Research/Project_GlialTopology/3.DataAnalysis/Exp2-Gfap,NeuN,Iba1/Results/Coordinates/ECM_Exp2_M05_30D_Gfap_Coordinates.csv"
Iba1 = "D:/Research/Project_GlialTopology/3.DataAnalysis/Exp2-Gfap,NeuN,Iba1/Results/Coordinates/ECM_Exp2_M05_30D_Iba1_Coordinates.csv"



#Define brain regions - leave array empty if plotting all cells
#CP = 672 SNr = GPe = GPi = PF = 930 PPN = STN = STR =477


def read_in_cells(filename):
    #this function reads in ClearMap csv files and can be used to filter to region at the same time
    #expects csv created using modified cellmap protocol, with region IDs and acronyms included
    #use pandas to read in specific columns from BrainJ csv output file
    cells = pd.read_csv(filename, usecols=[3,5,6,7,8], names = ['Ac','ID','Z','X','Y'], #, usecols=col_list)
                          skiprows = [0]) #skip header row

    print(cells.head())    
    Y = cells['X'].tolist()
    X = cells['Z'].tolist()
    Z = cells['Y'].tolist()
   # ID = cells['ID'].tolist()

    #If orientation was flipped during processing in ClearMap (e.g. slicing 
    # orientation=(1,-2,3) vs orientation=(1,2,3)) then reverse the order of that axis such as:
    #Z = Z[::-1]
    
    #data is in voxels, so multiply by atlas resolution (10micronXYZ)
    AtlasRes = 1
    X = [element * AtlasRes for element in X]
    Y = [element * AtlasRes for element in Y]
    Z = [element * AtlasRes for element in Z]
    pts = [[x, y, z] for x, y, z in zip(X, Y, Z)]
    return np.vstack(pts)

coordinates1 = read_in_cells(NeuN)
coordinates2 = read_in_cells(Gfap)
coordinates3 = read_in_cells(Iba1)

#Create the scene

scene = Scene(title="Cells")

#add in the relevant brain regions - use acronyms as below
scene.add_brain_region("CP", color = "lightblue", alpha=0.5)
scene.add_brain_region("HIP", color= "#07EEF7", alpha=0.5)
scene.add_brain_region("TH", color= "lightgreen", alpha=0.5)



#scene.add(Points(coordinates1, radius = 30, colors="#000000", alpha=0.5))
#scene.add(Points(coordinates2, radius = 50, colors="black", alpha=0.3))
#scene.add(Points(coordinates3, radius = 70, colors= "#CB1E21", alpha=1))

coordinates1B = coordinates1 
coordinates1B [:, 2] = -coordinates1B[:, 2] 

coordinates2B = coordinates2 
coordinates2B [:, 2] = -coordinates2B[:, 2] 

coordinates3B = coordinates3 
coordinates3B [:, 2] = -coordinates3B[:, 2] 


#coordinates3B [:, 2] = coordinates3B[:, 2] - 12500 

scene.add(PointsDensity(coordinates1B, radius = 1000, name="NeuN+"))
#scene.add(PointsDensity(coordinates2B, radius = 1000, name="Gfap+"))
#scene.add(PointsDensity(coordinates3B, radius = 1000, name="Iba1+"))

# render
scene.content
scene.render(zoom=0.9)

#scene.export("D:/Research/Project_GlialTopology/3.DataAnalysis/Exp2-Gfap,NeuN,Iba1/Plots/test.html")

#vm = VideoMaker(scene, "./Movies", "Movie1")
# Create video
# azimuth = X #elevation = Y #roll= X
# e.g. 1 = 1 degree per frame
#vm.make_video(duration=3, azimuth=5, fps=15)

#anim = Animation(scene, "./examples", "vid3")
#anim.add_keyframe(0, camera="top", zoom=1.3)
#anim.add_keyframe(1, camera="sagittal", zoom=3)
#anim.add_keyframe(2, camera="frontal", zoom=0.8)
#anim.add_keyframe(3, camera="frontal", )

#anim.make_video( duration=3, fps=15)

