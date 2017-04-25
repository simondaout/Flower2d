# protools
Geodetic tool package for 2-dimensional explorations of fault geometry and slip rates written in Python programming language. It can be utilized flexibly for plotting topography, InSAR and GPS data and modeling deformation sources. It provides various tectonic sources, such as decollement, ramp, pop-up, bookshelf, or flower structures.

 Requirements
=============
This project needs the following external components:
 * Python (2.7)
 * NumPy
 * SciPy 
 * Pymc
 * getopt
 * matplotlib 

 Source tree
============

  * src/flower: python tools based on pymc to invert for fault geometries and slip rates. The inversion code aims to explore several typical fault structures defined as follow in modelopti.py : 

- maindecol (decollement as main/first structure):

    - - - - - surface  - - - - - - - - - - - - - - - -
   |  
 H |       1   <--
   |  +------------------------
     (x,y)      --> 

-  mainflower (flower as main/first structure):

  - - - - - surface- - - - - - - - - - - - - - - -
H |
  +         + 
   \ºalpha / ºbeta
    \ 2   / 3
     \   /
      \ /           1   <--
       +------------------------
	  (x,y)             -->        
                                       
with:
1: the décollement
2: the ramp
3: the back-thrust


- ramp (secondary structure):

  - - - - - surface- - - - - - - - - - - - - - - -
H |
  +   
   \ºalpha  
    \ 1    
     \   
      \  
	   +	

-  flower structure (secodary structure):

 -+ - - - - +- - - - surface- - - - - - - - - - - -
   \ºalpha / ºbeta
    \ 1   / 2
     \   /
      \ /           
       +

Each segment is defined by a name, a strike-slip and dip-slip values, a dip angle, and a depth. The main segment (usually the decollement) needs to be position by its East, North coordinates, while the secondary segments are positioned by their horizontal distances to the main segment (parameter D). 

The code computes the surface displacements due to these edge dislocations or half-infinite strike-slip dislocations from [Segall 1996] equations in modelopti.py. Note that the surface displacements due to an half-infite strike-slip source is independant of its dip angle. The code imposes the conservation of the strike-slip and dip-slip motion along the various segments of the fault system.

The code supports GPS and InSAR data (class defined in networkopti.py).

The code can compute the full covariance matrix of the data vector following [Sudhaus and Jonsson, 2009] method (defined in networkopti.py).

The code saves all plaussible posterior models in the output directory and provide several plots of the results (defined in plot2d.py). 

  * src/profiles: python tool package to plot profiles across InSAR and GPS data. The code is made of these main classes:

- class network defined in network2d.py: reads InSAR or/and GPS data

- class fault2d defined in model2d.py : defines 2-dimensional fault in Eeast, West, strike coordinate 

- class prof defined in model2d.py: define position, azimuth and size of the profile to be plot

- class topo, seismi, moho defined in model2d.py: read data in x,y,z format for plots

- class gmtfiles defined in readgmt.py: read files in gmt format for plots

  * examples/haiyuam: inversion example  with a simple ramp-decollement structure or a flower structure geometry
  * example/socal: inversion example with a more complex fault geometries with several imbricated structures
  * example/atf: profile examples

Note: all data need to be projected in local UTM coordinates. Read readme.txt files in each directories for insar data, gps data or gmt files projection examples.
 
 References
============

  * Along-strike variations of the partitioning of convergence across the Haiyuan fault system detected by InSAR, S. Daout; R. Jolivet; C. Lasserre; M.-P. Doin; S. Barbot; P. Tapponnier; G. Peltzer; A. Socquet; J. Sun; GJI 2016 205 (1): 536-547 doi: 10.1093/gji/ggw028

  * Constraining the Kinematics of Metropolitan Los Angeles Faults with a Slip-Partitioning Model., S. Daout; S. Barbot, G. Peltzer, M.-P Doin, Z. Liu, R. Jolivet; GRL 2016.

