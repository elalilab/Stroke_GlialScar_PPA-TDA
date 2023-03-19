#
# Author: Dhananjay Bhaskar <dhananjay.bhaskar@yale.edu>
# Usage: srun --pty -t 16:00:00 --gpus a100:1 --cpus-per-gpu 16 --mem=128G --partition gpu zsh
#

import os, glob, copy, random

import numpy as np
import pandas as pd

import ripserplusplus as rpp_py

from gtda.pipeline import Pipeline
from gtda.homology import VietorisRipsPersistence
from gtda.diagrams import BettiCurve

from sklearn.ensemble import RandomForestClassifier

import matplotlib.pyplot as plt

### helper methods

def ripser2gtda(dgm, max_dim):
    
    diags = []
    for dim in range(max_dim+1):
        num_pts = len(dgm[dim])
        pers_diag = np.zeros((num_pts, 3))
        for idx in range(num_pts):
            pers_diag[idx,0] = dgm[dim][idx][0]
            pers_diag[idx,1] = dgm[dim][idx][1]
            pers_diag[idx,2] = dim
        diags.append(copy.deepcopy(pers_diag))
        
    return(np.vstack(diags))    

##### Ref: https://colab.research.google.com/drive/1addhGqN3ZE1mIn4L6jQnnkVs7_y__qSE?usp=sharing
##### Credit: Florent Poux
def grid_subsampling(points, voxel_size):

  nb_vox = np.ceil((np.max(points, axis=0) - np.min(points, axis=0))/voxel_size)
  non_empty_voxel_keys, inverse, nb_pts_per_voxel = np.unique(((points - np.min(points, axis=0)) // voxel_size).astype(int), axis=0, return_inverse=True, return_counts=True)
  idx_pts_vox_sorted = np.argsort(inverse)
  voxel_grid = {}
  grid_barycenter,grid_candidate_center = [],[]
  last_seen = 0

  for idx,vox in enumerate(non_empty_voxel_keys):
    voxel_grid[tuple(vox)] = points[idx_pts_vox_sorted[last_seen:last_seen+nb_pts_per_voxel[idx]]]
    grid_barycenter.append(np.mean(voxel_grid[tuple(vox)],axis=0))
    grid_candidate_center.append(voxel_grid[tuple(vox)][np.linalg.norm(voxel_grid[tuple(vox)]-np.mean(voxel_grid[tuple(vox)],axis=0),axis=1).argmin()])
    last_seen += nb_pts_per_voxel[idx]

  return grid_candidate_center
    
### params

cell_type = "Iba1"
max_dim = 1
grid_voxel_size = 120  #150

csv_files = glob.glob(os.path.join("data", cell_type, "*.csv"))
betti_curve = BettiCurve(n_bins=100, n_jobs=-1)

TDA_data = []

for csv_fname in csv_files:
    
    fname_parts = os.path.splitext(os.path.basename(csv_fname))[0].split("_")
    ctype = fname_parts[-1]
    dpi = fname_parts[-2]
    mid = fname_parts[-3]
    
    print(csv_fname + ":\t MouseID: " + mid[1:] + "\t Days Post Injury: " + dpi[:-1])
        
    #### load and sample point cloud

    df = pd.read_csv(csv_fname)
    coords_Iba = df[['X', 'Y', 'Z']].to_numpy()
    print("Total number of microglia (Iba1) : " + str(coords_Iba.shape[0]))
    pcd_Iba = np.array(grid_subsampling(coords_Iba, grid_voxel_size))
    print("Sampled " + str(pcd_Iba.shape[0]) + " microglia out of " + str(coords_Iba.shape[0]))

    dgm = rpp_py.run("--dim " + str(max_dim) + " --format point-cloud", pcd_Iba)
    gtda_dgm_Iba = ripser2gtda(dgm, max_dim)
    bc_Iba = betti_curve.fit_transform([gtda_dgm_Iba])

    df = pd.read_csv(os.path.join("data", "Gfap", f"ECM_Exp2_{mid}_{dpi}_Gfap.csv"))
    coords_Gfap = df[['X', 'Y', 'Z']].to_numpy()
    print("Total number of astrocytes (Gfap) : " + str(coords_Gfap.shape[0]))
    pcd_Gfap = np.array(grid_subsampling(coords_Gfap, grid_voxel_size))
    print("Sampled " + str(pcd_Gfap.shape[0]) + " astrocyctes out of " + str(coords_Gfap.shape[0]))

    dgm = rpp_py.run("--dim " + str(max_dim) + " --format point-cloud", pcd_Gfap)
    gtda_dgm_Gfap = ripser2gtda(dgm, max_dim)
    bc_Gfap = betti_curve.fit_transform([gtda_dgm_Gfap])

    df = pd.read_csv(os.path.join("data", "NeuN", f"ECM_Exp2_{mid}_{dpi}_NeuN.csv"))
    coords_NeuN = df[['X', 'Y', 'Z']].to_numpy()
    print("Total number of neurons (NeuN) : " + str(coords_NeuN.shape[0]))
    pcd_NeuN = np.array(grid_subsampling(coords_NeuN, grid_voxel_size))
    print("Sampled " + str(pcd_NeuN.shape[0]) + " neurons out of " + str(coords_NeuN.shape[0]))

    dgm = rpp_py.run("--dim " + str(max_dim) + " --format point-cloud", pcd_NeuN)
    gtda_dgm_NeuN = ripser2gtda(dgm, max_dim)
    bc_NeuN = betti_curve.fit_transform([gtda_dgm_NeuN])

    coords_Iba_Gfap = np.vstack([coords_Iba, coords_Gfap])
    pcd_Iba_Gfap = np.array(grid_subsampling(coords_Iba_Gfap, grid_voxel_size))
    print("Sampled " + str(pcd_Iba_Gfap.shape[0]) + " microglia & astrocyctes out of " + str(coords_Iba_Gfap.shape[0]))
    
    dgm = rpp_py.run("--dim " + str(max_dim) + " --format point-cloud", pcd_Iba_Gfap)
    gtda_dgm_Iba_Gfap = ripser2gtda(dgm, max_dim)
    bc_Iba_Gfap = betti_curve.fit_transform([gtda_dgm_Iba_Gfap])

    coords_Iba_NeuN = np.vstack([coords_Iba, coords_NeuN])
    pcd_Iba_NeuN = np.array(grid_subsampling(coords_Iba_NeuN, grid_voxel_size))
    print("Sampled " + str(pcd_Iba_NeuN.shape[0]) + " microglia & neurons out of " + str(coords_Iba_NeuN.shape[0]))

    dgm = rpp_py.run("--dim " + str(max_dim) + " --format point-cloud", pcd_Iba_NeuN)
    gtda_dgm_Iba_NeuN = ripser2gtda(dgm, max_dim)
    bc_Iba_NeuN = betti_curve.fit_transform([gtda_dgm_Iba_NeuN])

    coords_Gfap_NeuN = np.vstack([coords_Gfap, coords_NeuN])
    pcd_Gfap_NeuN = np.array(grid_subsampling(coords_Gfap_NeuN, grid_voxel_size))
    print("Sampled " + str(pcd_Gfap_NeuN.shape[0]) + " astrocytes & neurons out of " + str(coords_Gfap_NeuN.shape[0]))

    dgm = rpp_py.run("--dim " + str(max_dim) + " --format point-cloud", pcd_Gfap_NeuN)
    gtda_dgm_Gfap_NeuN = ripser2gtda(dgm, max_dim)
    bc_Gfap_NeuN = betti_curve.fit_transform([gtda_dgm_Gfap_NeuN])

    coords_Iba_Gfap_NeuN = np.vstack([coords_Iba, coords_Gfap, coords_NeuN])
    pcd_Iba_Gfap_NeuN = np.array(grid_subsampling(coords_Iba_Gfap_NeuN, grid_voxel_size))
    print("Sampled " + str(pcd_Iba_Gfap_NeuN.shape[0]) + " cells out of " + str(coords_Iba_Gfap_NeuN.shape[0]))

    dgm = rpp_py.run("--dim " + str(max_dim) + " --format point-cloud", pcd_Iba_Gfap_NeuN)
    gtda_dgm_Iba_Gfap_NeuN = ripser2gtda(dgm, max_dim)
    bc_Iba_Gfap_NeuN = betti_curve.fit_transform([gtda_dgm_Iba_Gfap_NeuN])

    save_data = {
		"Iba_ptcloud" : pcd_Iba,
		"Gfap_ptcloud" : pcd_Gfap,
		"NeuN_ptcloud" : pcd_NeuN,
        "Iba_Gfap_ptcloud" : pcd_Iba_Gfap,
        "Iba_NeuN_ptcloud" : pcd_Iba_NeuN,
        "Gfap_NeuN_ptcloud" : pcd_Gfap_NeuN,
        "Iba_Gfap_NeuN_ptcloud" : pcd_Iba_Gfap_NeuN,

		"Iba_dgm" : gtda_dgm_Iba,
		"Gfap_dgm" : gtda_dgm_Gfap,
		"NeuN_dgm" : gtda_dgm_NeuN,
        "Iba_Gfap_dgm" : gtda_dgm_Iba_Gfap,
        "Iba_NeuN_dgm" : gtda_dgm_Iba_NeuN,
        "Gfap_NeuN_dgm" : gtda_dgm_Gfap_NeuN,
        "Iba_Gfap_NeuN_dgm" : gtda_dgm_Iba_Gfap_NeuN,

		"Iba_bc" : bc_Iba,
		"Gfap_bc" : bc_Gfap,
		"NeuN_bc" : bc_NeuN,
        "Iba_Gfap_bc" : bc_Iba_Gfap,
        "Iba_NeuN_bc" : bc_Iba_NeuN,
        "Gfap_NeuN_bc" : bc_Gfap_NeuN,
        "Iba_Gfap_NeuN_bc" : bc_Iba_Gfap_NeuN
    }

    if not os.path.exists("results_pairs"):
        os.makedirs("results_pairs")
    ofilename = mid + "_" + dpi + "_TDA.npy"
    np.save(os.path.join("results_pairs", ofilename), save_data)    
