# Flower2d

Geodetic tool package for 2-dimensional explorations of fault geometry and slip rates, written in Python 3. It can be used to plot topography, InSAR and GPS data profiles and to model tectonic deformations. It provides various tectonic geometries such as ramps, flower structures, pop-up, and bookshelf faults.


Installation
=============

```bash
git clone https://github.com/simondaout/Flower2d.git
cd Flower2d
```

To update:
```bash
git pull
```

All scripts are standalone Python 3 modules — no installation required. Add `src/flower2d/` and `src/profiles/` to your `PYTHONPATH` or `PATH`.


Requirements
=============

| Package | Version |
|---------|---------|
| Python  | ≥ 3.9   |
| NumPy   | ≥ 1.22  |
| SciPy   | ≥ 1.9   |
| PyMC    | ≥ 4.0   |
| pytensor | ≥ 2.0  |
| arviz   | ≥ 0.14  |
| matplotlib | ≥ 3.5 |
| pyproj  | ≥ 3.0   |

Install dependencies with conda:
```bash
conda install -c conda-forge pymc pytensor arviz numpy scipy matplotlib pyproj
```


Branches
=============

| Branch   | Description |
|----------|-------------|
| `main`   | Python 3 — current development |
| `python2`| Legacy Python 2.7 code (read-only archive) |


Organisation
=============

```
flower2d/
├── src/
│   ├── flower2d/    ← Bayesian MCMC inversion of fault geometries
│   └── profiles/    ← Profile plotting and ramp estimation tools
└── examples/
    ├── haiyuam/     ← Ramp-decollement / flower structure example
    └── socal/       ← Southern California multi-structure example
```

### src/flower2d — MCMC inversion

Inverts InSAR and/or GPS profiles for fault geometry and slip rates using PyMC (Metropolis sampler).

Key files:
- `optimize.py` — main entry point; defines the PyMC model and runs the sampler
- `modelopti.py` — fault geometry classes (`mainfault`, `ramp`, `flower`, `popup`, `bookshelf`, …)
- `networkopti.py` — data loading classes (`network`: InSAR `dim=1`, GPS `dim=2/3`)
- `plot2d.py` — posterior plot routines
- `readgmt.py` — GMT segment file loader

### src/profiles — Profile tools

Plots data profiles and estimates residual ramps. Can be used independently of the inversion.

Key files:
- `plotPro.py` — main entry point
- `network2d.py` — data loading (same interface as `networkopti.network`)
- `model2d.py` — fault position overlays
- `readgmt.py` — GMT segment file loader
- `atan_fit.py`, `tanh_fit.py` — ramp-fitting utilities


Data format
=============

**InSAR** (`dim=1`): space-separated columns

| Columns | Description |
|---------|-------------|
| `x y los` | average LOS angles given via `av_los` / `av_heading` |
| `x y los look` | per-pixel incidence angle (col 4); `los=True` |
| `x y los look head` | per-pixel incidence + heading (cols 4–5); `los=True, head=True` |

- `x`, `y`: coordinates in km (local UTM) or degrees (if `utm_proj` is set)
- `los`: LOS displacement (mm/yr or mm)
- Nodata: pixels equal to `0.0`, `>9990`, or the value given in `nodata=` are masked automatically

**GPS** (`dim=2` or `dim=3`): station list file

```
# name  lon  lat
STAX    120.5  35.2
...
```

Individual station files in `wdir/reduction/<name>`:
- dim=2: `date  east  north  sigma_e  sigma_n`
- dim=3: `date  east  north  up  sigma_e  sigma_n  sigma_up`

Data can be provided in WGS84 (lon/lat) by setting `utm_proj=<EPSG>` (e.g. `utm_proj=32647`). The code reprojects automatically using `pyproj`.


network class — key parameters
=============

```
network(network, reduction, wdir, dim,
        weight=1., scale=1., errorfile=None, perc=100,
        los=False, head=False, av_heading=None, av_los=None, nodata=None,
        mask=None, cov=None, base=None,
        color='black', samp=1, lmin=None, lmax=None, width=0.5,
        utm_proj=None, ref=None)

network    : input text file name
reduction  : label used for output plots and GPS sub-directory
wdir       : path to input files
dim        : 1=InSAR, 2=GPS (E+N), 3=GPS (E+N+Up)
weight     : data weight (default: 1)
scale      : unit scaling factor (default: 1)
errorfile  : optional per-pixel uncertainty file (x y sigma)
perc       : percentile threshold for outlier removal (default: 100 = no removal)
los        : if True, read incidence angle from col 4 (default: False)
head       : if True, read heading angle from col 5 (default: False)
av_los     : average incidence angle, e.g. 23 for descending Envisat
av_heading : average heading angle, e.g. -76 for descending Envisat
nodata     : nodata sentinel value to mask (default: None)
mask       : [ymin, ymax] — exclude profile-perpendicular distance range
cov        : [sill, sigm, lamb] — exponential covariance parameters
base       : ramp uncertainties. InSAR: [sigma_cst, sigma_lin] (default: [10, 0.1])
             GPS: [sigma_E, sigma_N] or [sigma_E, sigma_N, sigma_Up] (default: [50, 50, 50])
samp       : spatial subsampling factor (default: 1)
utm_proj   : EPSG code for UTM reprojection (e.g. 32647)
ref        : [lon, lat] reference origin for UTM projection
```


Examples
=============

**Haiyuan fault** (`examples/haiyuam/`): ramp-decollement or flower structure geometry.

```bash
cd examples/haiyuam
plotPro.py tianzhuproin.py       # plot profiles
optimize.py tianzhuwestopti.py   # run MCMC inversion
```

**Southern California** (`examples/socal/`): imbricated multi-structure geometry.

```bash
cd examples/socal
plotPro.py gabrielproin.py
optimize.py gabrielflower.py
```


References
=============

- [Daout et al. (2016), *Along-strike variations of the partitioning of convergence across the Haiyuan fault system detected by InSAR*, GJI 205(1): 536–547](https://academic.oup.com/gji/article-abstract/205/1/536/2594860)

- [Daout et al. (2016), *Constraining the Kinematics of Metropolitan Los Angeles Faults with a Slip-Partitioning Model*, GRL](http://onlinelibrary.wiley.com/doi/10.1002/2016GL071061/full)
