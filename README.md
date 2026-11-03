# Fractures
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

## 2. File Structure

main.m
Main script that runs the complete workflow

computeG.m
Computes the local fracture geometry response G_i and adaptive side length L_i

computeS.m
Computes the directional continuity response S_i; each cell of the 3x3 grid is treated as a segment

gmmDoubleThreshold.m
Dual thresholds based on a three-component Gaussian mixture model

doubleOtsuThreshold.m
Dual thresholds based on three-class Otsu

kmeansDoubleThreshold.m
Dual thresholds based on K-means (K = 3)

maxEntropyThreshold.m
Dual thresholds based on maximum entropy

hysteresisTrackingSpatial.m
Adaptive hysteresis tracking

dbscanAndOrientations.m
DBSCAN grouping, segment-wise orientation extraction, and segment-level D-value filtering

Extractionfractureinformation.m
Segment-wise PCA and Mean Shift clustering for orientation extraction

interpolateToFull.m
Mapping from the downsampled point cloud to the full point cloud

---

## 3. Input Data

Format: .ply point cloud file;
Requirement: data fusion and vegetation removal should already be completed;
Recommended point density: 1000-2000 points per square metre;
Coordinate unit: metre.

---

## 4. Parameters

voxelSize
Voxel grid size
Determined by the average point spacing

Pn
Front direction of the point cloud model
Determined by the slope aspect

L_max
Maximum fracture aperture
Determined by field investigation

alpha
Steepness of the Sigmoid function
Set to 8

bandwidth_ori
Angular bandwidth of Mean Shift
Set to 15 degrees

minPts
Minimum number of points in a cell
Set to 10

---

## 5. Workflow

### 5.1 Point cloud preprocessing

ptCloud = pcread('your_point_cloud.ply');
ptCloudDown = pcdownsample(ptCloud, 'gridAverage', voxelSize);

### 5.2 Computation of G_i and L_i

[G, L_i] = computeG(ptCloudDown, Pn, q, alpha, L_max);

L_i is determined by the peak of the local geometric response. G_i is obtained by combining the curvature response and the centroid distance response, each normalized by a robust z-score, and multiplied by the directional constraint.

### 5.3 Computation of S_i

mask = true(size(ptCloudDown.Location, 1), 1);
[S, normals] = computeS(ptCloudDown, L_i, mask, thetaThr);

For each candidate point, a 3x3 grid is constructed on the local tangent plane. Each cell is treated as a segment. PCA is performed on the points within each cell to obtain its extension direction. The directional consistency response between a surrounding cell and the central cell is defined as:

w_j = 1 - theta_j / 90

S_i is the average of w_j over all valid cells:

S_i = (1 / n_valid) * sum(w_j)

Cells with insufficient points are excluded from n_valid.

### 5.4 Fusion of Nf

Nf = G .* S;

### 5.5 Dual-threshold segmentation

[B, T, S_T_idx, S_B_idx, S_cand_idx] = gmmDoubleThreshold(Nf);

Points with Nf > T are classified as FCPs, points with B < Nf <= T as SFCPs, and points with Nf <= B are suppressed.

### 5.6 Hysteresis tracking

[fracture_mask, ~] = hysteresisTrackingSpatial(ptCloudDown, Nf, L_i, S_T_idx, S_cand_idx, B);

### 5.7 DBSCAN grouping and orientation extraction

[cluster_results, idx_crack_clean, labels_clean] = dbscanAndOrientations(ptCloudDown, idx_crack, L_i, bandwidth_ori);

Each DBSCAN cluster corresponds to one continuous fracture. Within each cluster, the point cloud is divided into segments by KNN neighbourhoods. PCA is performed on each segment to obtain its normal vector. Mean Shift clustering is then applied to all segment normals, and the centre of the largest cluster is taken as the orientation of the fracture.

### 5.8 Orientation comparison

The main direction of each cluster in cluster_results is converted into dip direction and dip angle, and compared with field measurements.

---

## 6. Outputs

L_i
Adaptive side length of each point

Nf
Fracture feature value of each point

fracture_mask
Fracture mask after hysteresis tracking

idx_crack_clean
Fracture point indices after DBSCAN and segment-level filtering

cluster_results
Orientation and segment data for each fracture cluster

azimuths_dips
Dip direction and dip angle of each segment

---

## 7. Notes

1. mask setting: if no pre-screening is applied, mask should be all ones; if pre-screening is retained, mask should correspond to the candidate point indices.
2. L_max: it must be determined by field investigation. It serves both as the upper bound of the adaptive side length and as the basis for computing numNeighbors.
3. Pn: it is determined by the slope aspect and affects both the directional constraint and the orientation calculation.
4. DBSCAN parameters: eps is set to three times the median nearest-neighbour distance, and minPts is set to 10.
5. Segment-level filtering: a cluster is discarded as misidentified if it contains fewer than two fracture-like segments. Each segment is a KNN neighbourhood, not a single point.

---

## 8. Citation

If you use this method, please cite:

Liu Y, Wang W, Xie S, He B, Xiao H, Wu L. A point cloud technique for identifying fractures in unstable rock masses. (Under review)

---

## 9. Contact

Corresponding author: Lizhou Wu, Ph.D., Professor
Email: lzuw@cqjtu.edu.cn
Affiliation: Chongqing Jiaotong University
