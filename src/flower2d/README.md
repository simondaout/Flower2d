# optimize.py 

Python script based on PyMC to invert for fault geometries and slip rates. The inversion code aims to explore several typical fault structures.
 
Profile:
=============
Everything is defined relatively to the center of the profile and in the 2-dimension in the along-profile distance (class defined in `modelopti.py`).

```
Parameters profile: 
    name: name profile
    x,y: reference point: everything is refered to this point 
    l,w: length, width profile
    strike: strike profile (default: defined by inversion class)
    proj=[east, notth, up]: average LOS projection into east, north, up used for plots [default: None]
    type:  std - plot mean and standard deviation InSAR;
    distscale - scatter plot with color scale function of the profile-parallel distance;
    stdscat - plot scatter + standar deviation. 
```

Data:
=============
The code supports GPS and InSAR data (class defined in `networkopti.py`).

```
Parameters notwork: 
    network: name input text file
    reduction: reduction name for plots
    wdir: relative path input file
    dim: 1=InSAR, 2 or 3=GPS
    weight: weight for inversion (default: 1)
    scale: scaling value input data (default: 1)
    errorfile: optional error file  (default: None)
    los: if True then read look angle in the 4th column of InSAR text file (default: False)
    head: if True then read heading angle in the 5th column of InSAR text file (default: False)
    av_los: average los angle value (eg.: 23 as for desc Envisat) (default: None)
    av_heading: average heading angle value (eg.: -76 as for desc Envisat) (default: None)
    cov=[sill, sigm, lamb]: covariance parameters, default: None,None,None
    mask: optional mask (default: None)
    base: uncertainties for reference frame
    if InSAR data, give uncertainties for cst and linear term of the ramp, default: [10, 0.1]
    if GPS data, give uncertainties for each components, default: [50,50,50]
    color: plot option, default: 'black' 
    samp: subsample option, default:1 
    perc: cleaning outliers option (default: 100)
    lmin,lmax: min max options for plots
    width: size scatter point for plots (default: 3.)
```

Inversion:
=============

```
Defined your inversion parameters: 
    structure: define all structures (see bellow)
    strike: azimuth inversion profile (ie perpandicular to the main fault)
    profile: profile class 
    fullcov      :   if True, estimates a full covariance matrix for InSAR data
                     if False, returns a diagonal matrix with sigmad on the diagonal
    maskcov: give a mask for covariance estimation on sub-sampled data in perpendicular 
    distance to the center of profile km eg. [0,20,40,50] (default: None)
    rampcov: remove ramp before covariance estimation. eg. lin, quad, cub (default: lin)
    name: give a name (Optional)
    depthmax: maximum depth for plots (Optional)
```

Tectonic structures are defined as follow: 

```
    - mainfault (Defined main segment):


    - - - - - surface  - - - - - - - - - - - - - - - -
   |  
 w |       1   <--
   |  +------------------------
     (D)      --> 

Free parameters: ss (strike-slip), short (shortening), dip (dip-angle),
w (locking-depth), D (horizontal distance to the center of profile), L (length) Default: D=0 (fault located in the center of the profile), L=660 (half-infinite segement)

Uncertainties: sigmass, sigmashort, sigmaw, sigmaD, sigmaL must be given in the same units that parameters. Default: sigmaD=0, sigmaL=0 (paramters not explored).

Distributions: distss,distshort,distH,distL,distD,distdip defined prior distributions. Default: Unif (Uniform distributions)


    - mainflower() (flower structure as the main structure):


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

Feee parameters: ss, short, w, sigmaw, dip, D, L (as for mainfault) 
H2 (vertical distance of segement 2 to the main fault 1), ss2 (strike-slip segment 2), D2 (horizontal distance segment 2 to main fault 1)
H3 (vertical distance of segement 3 to the main fault 1), ss3 (strike-slip segment 2), D3 (horizontal distance segment 2 to main fault 1)

    - ramp() (secondary fault branching to the previous segment):

  - - - - - surface- - - - - - - - - - - - - - - -
H |
  +   
   \ºgamma  
    \ 1    
     \   
      \  
	   +	

Free parameters: ss, D, H
    
    -  flower() (secodary structure branching to the previous segment):


 -+ - - - - +- - - - surface- - - - - - - - - - - -
   \ºgamma / ºbeta
    \ 1   / 2
     \   /
      \ /           
       +  

    - ...

Feee parameters: 
H1 (vertical distance of segement 1 to the previous fault), ss1 (strike-slip segment 1), D1 (horizontal distance segment 1 to main previous fault)
H2 (vertical distance of segement 2 to the  fault 1), ss2 (strike-slip segment 2), D2 (horizontal distance segment 2 to  fault 1)
H3 (vertical distance of segement 3 to the  fault 1), ss3 (strike-slip segment 2), D3 (horizontal distance segment 2 to  fault 1)

```

The code computes the surface displacements due to these edge dislocations or half-infinite strike-slip dislocations from [Segall 1996] equations in `modelopti.py`. Note that the surface displacements due to an half-infite strike-slip source is independant of its dip angle. 

The code imposes the conservation of the strike-slip and dip-slip motion along the various segments of the fault system such as:

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

The code can compute the full covariance matrix of the data vector following [Sudhaus and Jonsson, 2009] method (defined in `networkopti.py`).

Plots:
=============

The code saves all plaussible posterior models in the output directories and provide several plots of the results (defined in `plot2d.py`). 

