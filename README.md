# Flower2d
Geodetic tool package for 2-dimensional explorations of fault geometry and slip rates written in Python programming language. It can be utilized flexibly for plotting topography, InSAR and GPS data profiles and modeling tectonic deformations. It provides various tectonic geometries, such as ramps, flower structures, pop-up, or bookshelf.


To download the package
=============
```git clone https://github.com/simondaout/Flower2d.git```

To update the software: 
```git pull```

* All scripts are written in python and do not need installation. Just add them to your PATH. 

Requirements
=============
This project needs the following external components:
 * Python (2.7)
 * NumPy
 * SciPy 
 * PyMC
 * getopt
 * matplotlib 


Organisation
=============
This project contains the following folders:
 
* src/profiles: Python tool package to plot profiles and estimate residual ramps across InSAR and GPS data. 

* src/flower: Python inversion tool based on PyMC to invert for fault geometries and slip rates.

* examples/haiyuam: profile and inversion example  with a simple ramp-decollement structure or a flower structure geometry. All input parameters are difined in `tianzhuwestopti.py` and `tianzhuproin.py`
To run the profile: `plotPro.py tianzhuproin.py`
To run the inversion: `optimize.py tianzhuwestopti.py`

* example/socal: profile and inversion examples. Input file for profile is `gabrielproin.py`. Input file for inversion is `gabrielflower.py`. Example of inversion with a more complex fault geometries with several imbricated structures. 
To run the profile: `plotPro.py gabrielproin.py`
To run the inversion: `optimize.py gabrielflower.py`

:memo: Note: all data need to be projected in local UTM coordinates. Read readme.txt files in each directories for insar data, gps data or gmt files projection examples.

 
 References
============

  * [Along-strike variations of the partitioning of convergence across the Haiyuan fault system detected by InSAR, S. Daout; R. Jolivet; C. Lasserre; M.-P. Doin; S. Barbot; P. Tapponnier; G. Peltzer; A. Socquet; J. Sun; GJI 2016 205 (1): 536-547 doi: 10.1093/gji/ggw028](https://academic.oup.com/gji/article-abstract/205/1/536/2594860/Along-strike-variations-of-the-partitioning-of)

  * [Constraining the Kinematics of Metropolitan Los Angeles Faults with a Slip-Partitioning Model., S. Daout; S. Barbot, G. Peltzer, M.-P Doin, Z. Liu, R. Jolivet; GRL 2016](http://onlinelibrary.wiley.com/doi/10.1002/2016GL071061/full)

