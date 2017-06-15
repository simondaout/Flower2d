# protools
Geodetic tool package for 2-dimensional explorations of fault geometry and slip rates written in Python programming language. It can be utilized flexibly for plotting topography, InSAR and GPS data and modeling deformation sources. It provides various tectonic sources, such as decollement, ramp, flower structure, pop-up, or bookshelf.

 Requirements
=============
This project needs the following external components:
 * Python (2.7)
 * NumPy
 * SciPy 
 * PyMC
 * getopt
 * matplotlib 

 Source tree
============

  * src/flower: python inversion tool based on PyMC to invert for fault geometries and slip rates. The inversion code aims to explore several typical fault structures defined as follow in `modelopti.py` : 

```
    - maindecol() (decollement as main structure):


    - - - - - surface  - - - - - - - - - - - - - - - -
   |  
 H |       1   <--
   |  +------------------------
     (x,y)      --> 


    - mainflower() (flower as main structure):


   - - - - - surface- - - - - - - - - - - - - - - -
 H |
   +         + 
    \ºgamma / ºbeta
     \ 2   / 3
      \   /
       \ /           1   <--
        +------------------------
  (x,y)  \ ºalpha       -->        
          \                             
      
	 with:
      1: the décollement
      2: the ramp
      3: the back-thrust

    - ramp() (secondary structure):

  - - - - - surface- - - - - - - - - - - - - - - -
H |
  +   
   \ºgamma  
    \ 1    
     \   
      \  
	   +	
    
    -  flower() (secodary structure):


 -+ - - - - +- - - - surface- - - - - - - - - - - -
   \ºgamma / ºbeta
    \ 1   / 2
     \   /
      \ /           
       +  

    - ...

```

Each segment is defined by a name, a strike-slip and dip-slip values, a dip angle, and a depth. The main segment (usually the decollement) needs to be position by its East and North coordinates, while the secondary segments are positioned by their horizontal distances to the main segment (parameter D). 

The code computes the surface displacements due to these edge dislocations or half-infinite strike-slip dislocations from [Segall 1996] equations in `modelopti.py`. Note that the surface displacements due to an half-infite strike-slip source is independant of its dip angle. The code imposes the conservation of the strike-slip and dip-slip motion along the various segments of the fault system such as:

```
	- the deep-seated strike-slip rate (i.e far-field strike-slip motion imposed by the slip 
  on the decollement) is superior or equal to the sum of the strike-slip rates on all segments:
					

	- the dip slip rate on each structures is imposed by the geometry of the fault system:


    V
     2        sin(beta - gamma + alpha)
   --- =       ----------------------
    V             sin(beta - gamma)
     1


    V
     3               sin(alpha)
   --- =       ----------------------
    V             sin(beta - gamma)
     1

```   

Each prior distributions are defined in the input file by their sigma value and their distributions (refered to `modelopti.py`). Several prior distributions are available (default :'Unif' for Uniform distribution, 'Gaus': Gaussian distribution, 'Logn': Log Normal distribution).

The code saves all plaussible posterior models in the output directory and provide several plots of the results (defined in `plot2d.py`). 


The code supports GPS and InSAR data (class defined in `networkopti.py`).

The code can compute the full covariance matrix of the data vector following [Sudhaus and Jonsson, 2009] method (defined in `networkopti.py`).


  * src/profiles: python tool package to plot profiles across InSAR and GPS data. The code is made of these main classes:

    - class `network()` defined in `network2d.py`: reads InSAR or/and GPS data
    - class `prof()` defined in `model2d.py`: define position, azimuth and size of the profile to be plot. It also defined the type offigure you desire:
		- default: scatter plot 
		- type = 'std': plot mean and standard deviation insar
		- type = 'distscale': scatter plot with color scale fonction of the profile-parallel distance		
		- type = 'stdscat': plot mean and standard deviation	

    - class `topo()`, `seismi()`, `moho()` defined in `model2d.py`: read data in x,y,z format for plots
    - class `gmtfiles()` defined in `readgmt.py`: read files in gmt format for plots

  * examples/haiyuam: inversion example  with a simple ramp-decollement structure or a flower structure geometry. All input parameters are difined in `tianzhuwestopti.py`.
  
  * example/socal: profile and inversion examples. Input file for profile is `gabrielproin.py`. Input file for inversion is `gabrielflower.py`. Example of inversion with a more complex fault geometries with several imbricated structures. 

:memo: Note: all data need to be projected in local UTM coordinates. Read readme.txt files in each directories for insar data, gps data or gmt files projection examples.
 
 References
============

  * [Along-strike variations of the partitioning of convergence across the Haiyuan fault system detected by InSAR, S. Daout; R. Jolivet; C. Lasserre; M.-P. Doin; S. Barbot; P. Tapponnier; G. Peltzer; A. Socquet; J. Sun; GJI 2016 205 (1): 536-547 doi: 10.1093/gji/ggw028](https://academic.oup.com/gji/article-abstract/205/1/536/2594860/Along-strike-variations-of-the-partitioning-of)

  * [Constraining the Kinematics of Metropolitan Los Angeles Faults with a Slip-Partitioning Model., S. Daout; S. Barbot, G. Peltzer, M.-P Doin, Z. Liu, R. Jolivet; GRL 2016](http://onlinelibrary.wiley.com/doi/10.1002/2016GL071061/full)

