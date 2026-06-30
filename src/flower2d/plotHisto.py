#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plotHisto.py — Plot MCMC posterior histograms from flower2d optimize.py output.

Figure 1: direct sampled posteriors
  (ZT short, depth, D, L, dip | Ramp1 H, D | Ramp2 H, D)
Figure 2: derived quantities
  (Ramp1 dip, L, top-depth | Ramp2 dip, L, top-depth)

Usage
-----
    python3 plotHisto.py <stat_dir> [output_dir] [--filter-mode] [--winit FLOAT]

    stat_dir     : directory containing the .txt trace files
    output_dir   : where to save the PDFs (defaults to stat_dir)
    --filter-mode: keep only samples within the FWHM of the primary KDE peak
    --winit FLOAT: initial ZT depth (km) used during inversion (default: 10.0)
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import arviz as az
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks


# ======================================================================================================================
# KDE mode finder (same as plot_histo.py)
# ======================================================================================================================

def primary_mode_window(samples, n_grid=2000):
    """Return (mode, lo, hi) — position and FWHM window of the primary KDE peak."""
    kde  = gaussian_kde(samples, bw_method='silverman')
    x    = np.linspace(samples.min(), samples.max(), n_grid)
    y    = kde(x)

    peaks, _ = find_peaks(y)
    if len(peaks) == 0:
        return x[np.argmax(y)], x[0], x[-1]

    primary  = peaks[np.argmax(y[peaks])]
    mode     = x[primary]
    half_max = y[primary] / 2

    left_idx  = np.where(y[:primary] < half_max)[0]
    right_idx = np.where(y[primary:] < half_max)[0]
    lo = x[left_idx[-1]]           if len(left_idx)  else x[0]
    hi = x[primary + right_idx[0]] if len(right_idx) else x[-1]
    return mode, lo, hi


def filter_to_mode(samples):
    _, lo, hi = primary_mode_window(samples)
    sel = samples[(samples >= lo) & (samples <= hi)]
    return sel if len(sel) >= 10 else samples


# ======================================================================================================================
# Histogram helper
# ======================================================================================================================

def histo(ax, post, label, xlabel='', color='steelblue', filter_mode=False):
    """Plot a posterior histogram with mean/mode and 95% HDI."""
    post = np.asarray(post, dtype=float)
    post = post[np.isfinite(post)]
    if len(post) == 0:
        ax.set_title(f"{label}\n(no data)")
        return

    mode, _, _ = primary_mode_window(post)

    # Full distribution (light)
    ax.hist(post, bins=40, density=True, histtype='stepfilled',
            color=color, alpha=0.25)

    if filter_mode:
        sel    = filter_to_mode(post)
        suffix = " (mode)"
    else:
        sel    = post
        suffix = ""

    # Mode-selected distribution (solid)
    ax.hist(sel, bins=40, density=True, histtype='stepfilled',
            color=color, alpha=0.55)

    mean   = sel.mean()
    lo, hi = az.hdi(sel, hdi_prob=0.95)

    ax.axvline(mean, color='crimson', lw=1.5)
    ax.axvline(lo,   color='crimson', lw=1.0, linestyle='--', alpha=0.7)
    ax.axvline(hi,   color='crimson', lw=1.0, linestyle='--', alpha=0.7)
    ax.axvline(mode, color='black',   lw=1.0, linestyle=':',  alpha=0.8)
    ax.axvspan(lo, hi, alpha=0.10, color='crimson')

    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_title(
        f'{label}{suffix}\n{mean:.3g}  [{lo:.3g} – {hi:.3g}]',
        fontsize=8
    )
    ax.tick_params(labelsize=8)
    ax.set_yticks([])


# ======================================================================================================================
# Argument parsing
# ======================================================================================================================

args        = sys.argv[1:]
filter_mode = '--filter-mode' in args
args        = [a for a in args if a != '--filter-mode']

winit = 10.0
if '--winit' in args:
    idx   = args.index('--winit')
    winit = float(args[idx + 1])
    args  = args[:idx] + args[idx + 2:]

if len(args) < 1:
    print(__doc__)
    sys.exit(1)

wdir   = args[0]
outdir = args[1] if len(args) > 1 else wdir
os.makedirs(outdir, exist_ok=True)


# ======================================================================================================================
# Load traces
# ======================================================================================================================

def load(fname):
    path = os.path.join(wdir, fname)
    if not os.path.exists(path):
        alt = path.replace('_D.txt', '_L.txt')
        if os.path.exists(alt):
            print(f'[warn] using legacy file {os.path.basename(alt)} for {fname}')
            path = alt
        else:
            print(f'[error] file not found: {path}')
            sys.exit(1)
    return np.loadtxt(path)


ZT_short = load('ZT_short.txt')
ZT_H     = load('ZT_H.txt')
ZT_D     = load('ZT_D.txt')
ZT_L     = load('ZT_L.txt')
ZT_dip   = load('ZT_dip.txt')
R1_H     = load('Ramp1_H.txt')
R1_D     = load('Ramp1_D.txt')
R2_H     = load('Ramp2_H.txt')
R2_D     = load('Ramp2_D.txt')

ZT_depth = winit - ZT_H


# ======================================================================================================================
# Figure 1: direct posteriors
# ======================================================================================================================

fig1, axes = plt.subplots(3, 3, figsize=(13, 9), constrained_layout=True)
axes = axes.flatten()

histo(axes[0], ZT_short,  label='ZT shortening',      xlabel='mm/yr', color='steelblue',     filter_mode=filter_mode)
histo(axes[1], ZT_depth,  label='ZT depth',            xlabel='km',    color='steelblue',     filter_mode=filter_mode)
histo(axes[2], ZT_D,      label='ZT D',                xlabel='km',    color='steelblue',     filter_mode=filter_mode)
histo(axes[3], ZT_L,      label='ZT length',           xlabel='km',    color='steelblue',     filter_mode=filter_mode)
histo(axes[4], ZT_dip,    label='ZT dip',              xlabel='°',     color='steelblue',     filter_mode=filter_mode)
histo(axes[5], R1_H,      label='Ramp1 H (vertical)',  xlabel='km',    color='coral',         filter_mode=filter_mode)
histo(axes[6], R1_D,      label='Ramp1 D (horiz.)',    xlabel='km',    color='coral',         filter_mode=filter_mode)
histo(axes[7], R2_H,      label='Ramp2 H (vertical)',  xlabel='km',    color='mediumseagreen', filter_mode=filter_mode)
histo(axes[8], R2_D,      label='Ramp2 D (horiz.)',    xlabel='km',    color='mediumseagreen', filter_mode=filter_mode)

suffix = ' — mode-filtered' if filter_mode else ''
fig1.suptitle(f'flower2d — direct posteriors{suffix}', fontsize=11)
out1 = os.path.join(outdir, 'histo_direct.pdf')
fig1.savefig(out1, dpi=150)
print(f'Saved: {out1}')


# ======================================================================================================================
# Figure 2: derived quantities
# ======================================================================================================================

R1_dip   = np.rad2deg(np.arctan2(R1_H, -R1_D))
R1_L     = np.sqrt(R1_D**2 + R1_H**2)
R1_depth = ZT_depth - R1_H

R2_dip   = np.rad2deg(np.arctan2(R2_H, -R2_D))
R2_L     = np.sqrt(R2_D**2 + R2_H**2)
R2_depth = ZT_depth - R2_H

fig2, axes2 = plt.subplots(2, 3, figsize=(13, 6), constrained_layout=True)
axes2 = axes2.flatten()

histo(axes2[0], R1_dip,   label='Ramp1 dip',          xlabel='°',  color='coral',          filter_mode=filter_mode)
histo(axes2[1], R1_L,     label='Ramp1 L (down-dip)', xlabel='km', color='coral',          filter_mode=filter_mode)
histo(axes2[2], R1_depth, label='Ramp1 top depth',    xlabel='km', color='coral',          filter_mode=filter_mode)
histo(axes2[3], R2_dip,   label='Ramp2 dip',          xlabel='°',  color='mediumseagreen', filter_mode=filter_mode)
histo(axes2[4], R2_L,     label='Ramp2 L (down-dip)', xlabel='km', color='mediumseagreen', filter_mode=filter_mode)
histo(axes2[5], R2_depth, label='Ramp2 top depth',    xlabel='km', color='mediumseagreen', filter_mode=filter_mode)

for lbl, trace in [('Ramp1', R1_dip), ('Ramp2', R2_dip)]:
    print(f'{lbl} dip   — Mean: {trace.mean():.2f}°   HDI: {az.hdi(trace, hdi_prob=0.95)}')
for lbl, trace in [('Ramp1', R1_L), ('Ramp2', R2_L)]:
    print(f'{lbl} L     — Mean: {trace.mean():.2f} km  HDI: {az.hdi(trace, hdi_prob=0.95)}')
for lbl, trace in [('Ramp1', R1_depth), ('Ramp2', R2_depth)]:
    print(f'{lbl} depth — Mean: {trace.mean():.2f} km  HDI: {az.hdi(trace, hdi_prob=0.95)}')

fig2.suptitle(f'flower2d — derived ramp geometry{suffix}', fontsize=11)
out2 = os.path.join(outdir, 'histo_derived.pdf')
fig2.savefig(out2, dpi=150)
print(f'Saved: {out2}')

plt.show()
