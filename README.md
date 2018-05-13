# BundleAdjustmentWith_G-N
Bundle Adjustment is a least-squares objective function for reprojection errors. Researchers demonstrate that Bundle Adjustment can be implemented by Ceres and g2o and can be used for PnP pose estimation. In this question, I write a simple Gaussian Newton method , and use Bundle Adjustment to optimize the pose.  

The whole demo is tested in Ubuntu Platorforms.

# Mathematical derivation
<div align=center>  
  
![](https://github.com/TianQi-777/BundleAdjustmentWith_G-N/blob/master/Images/G-N.png)
</div>

# Data description
**p2d.txt**:Pixel coordinates in reprojection.  
**Data storageform**:u v

**p3d.txt**:Map point in reprojection. 
**Data storage form**:x y z  

fx=520.9 fy=521.0 cx=325.1 cy=249.7

# Additional Prerequisites for this demo
**Sophus**  
We use [Sophus](https://github.com/strasdat/Sophus) for Lie groups commonly used for 2d and 3d geometric problems. 
Dowload and install instructions can be found at: https://github.com/strasdat/Sophus.

## Build and Run
```
cd XX/XX(include GN-BA.cpp ,p2d.txt ,p3d.txt and CMakeLists.txt)  
mkdir build  
cd build  
cmake ..  
make -j2  
./GN-BA
```
