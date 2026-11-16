# Identification and Orientation Extraction of Shallow Fractures on High and Steep Rock Slopes Using Point Cloud Data

## 1. Overview

This repository provides a point cloud-based method for identifying shallow fractures on high and steep rock slopes and extracting their orientations. The method combines digital image processing techniques with 3D point cloud analysis.

The workflow consists of:

1. Point cloud preprocessing: voxel downsampling, vegetation removal, KD-tree construction, and average point spacing estimation;
2. Adaptive side length: determination of the characteristic scale L_i for each point by searching for the peak of the local geometric response;
3. Semantic segmentation criterion: construction of the fracture feature value Nf, which replaces the image grey value;
4. 3D edge detection: dual-threshold suppression and adaptive hysteresis tracking;
5. DBSCAN grouping and segment-level filtering: removal of isolated clusters and misidentified clusters;
6. Orientation extraction: segment-wise PCA and Mean Shift clustering;
7. Accuracy evaluation: comparison with field-measured orientations.

---

## 2. Quick Start (How to Run)

To quickly test this algorithm, follow these steps:

1. Clone or download this repository to your local machine.
2. Open the `Main.m` file in MATLAB.
3. In the code, locate the point cloud reading section (around line 8) and input the absolute path of the test file:
   ```matlab
   % ---------- 1. Load point cloud ----------
   ptCloud = pcread('your_file_path\Test.pcd'); ;

## 3. Requirements
MATLAB (Developed and tested on R2019b)

Computer Vision Toolbox (for point cloud processing and KNN search)

Statistics and Machine Learning Toolbox (for GMM, K-means, and DBSCAN)

## 4. File Structure
main.m: Main script that runs the complete workflow.

computeG.m: Computes local fracture geometry response Gi and adaptive side length Li.

computeS.m: Computes directional continuity response Si.

gmmDoubleThreshold.m: Dual thresholds based on Gaussian mixture model (GMM).

doubleOtsuThreshold.m: Dual thresholds based on three-class Otsu.

kmeansDoubleThreshold.m: Dual thresholds based on K-means (K = 3).

maxEntropyThreshold.m: Dual thresholds based on maximum entropy.

hysteresisTrackingSpatial.m: Adaptive hysteresis tracking.

dbscanAndOrientations.m: DBSCAN grouping and segment-wise orientation extraction.

Extractionfractureinformation.m: Segment-wise PCA and Mean Shift clustering.

interpolateToFull.m: Mapping from downsampled to full point cloud.

## 5. Input Data
Format: .pcd or .ply point cloud file.

Requirement: Data fusion and vegetation removal should be completed beforehand.

## 6. Outputs
Li: Adaptive side length of each point.

Nf: Fracture feature value of each point.

fracture_mask: Fracture mask after hysteresis tracking.

cluster_results: Orientation and segment data for each fracture cluster.

## 7. License
This project is licensed under the MIT License - see the LICENSE file for details.
