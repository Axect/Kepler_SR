from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scienceplots

# Import netCDF file
ncfile = './kepler.nc'
data = Dataset(ncfile)
var = data.variables

# Prepare Data to Plot
t = var['t'][:]
x = var['x'][:]
y = var['y'][:]  
vx = var['vx'][:]
vy = var['vy'][:]
E = var['E'][:]
L = var['L'][:]
g = var['g'][:]

g_unique = np.unique(g)

df = pd.DataFrame({
    't': t,
    'x': x,
    'y': y,
    'vx': vx,
    'vy': vy,
    'E': E,
    'L': L,
    'g': g
})

# Plot params
pparam = dict(
    xlabel = r'$x$',
    ylabel = r'$y$',
    xscale = 'linear',
    yscale = 'linear',
)

# Plot
with plt.style.context(["science", "nature"]):
    fig, ax = plt.subplots()
    ax.autoscale(tight=True)
    ax.set(**pparam)
    for grp in g_unique:
        ax.plot(df.x[df.g == grp], df.y[df.g == grp], label=grp)
    ax.legend()
    fig.savefig('pos.png', dpi=300, bbox_inches='tight')

# Plot
with plt.style.context(["science", "nature"]):
    fig, ax = plt.subplots()
    ax.autoscale(tight=True)
    ax.set(**pparam)
    for grp in g_unique:
        ax.plot(df.vx[df.g == grp], df.vy[df.g == grp], label=grp)
    ax.legend()
    fig.savefig('vel.png', dpi=300, bbox_inches='tight')

# Plot
with plt.style.context(["science", "nature"]):
    fig, ax = plt.subplots()
    ax.autoscale(tight=True)
    ax.set(**pparam)
    for grp in g_unique:
        ax.plot(df.t[df.g==grp], df.x[df.g == grp], label=f"x_{grp}")
        ax.plot(df.t[df.g==grp], df.y[df.g == grp], label=f"y_{grp}")
    ax.legend()
    fig.savefig('pos_t.png', dpi=300, bbox_inches='tight')

# Plot
with plt.style.context(["science", "nature"]):
    fig, ax = plt.subplots()
    ax.autoscale(tight=True)
    ax.set(**pparam)
    for grp in g_unique:
        ax.plot(df.t[df.g==grp], df.vx[df.g == grp], label=f"vx_{grp}")
        ax.plot(df.t[df.g==grp], df.vy[df.g == grp], label=f"vy_{grp}")
    ax.legend()
    fig.savefig('vel_t.png', dpi=300, bbox_inches='tight')

# Plot
with plt.style.context(["science", "nature"]):
    fig, ax = plt.subplots()
    ax.autoscale(tight=True)
    ax.set(**pparam)
    for grp in g_unique:
        ax.plot(df.t[df.g==grp], df.E[df.g == grp], label=f"{grp}")
    ax.legend()
    fig.savefig('energy_t.png', dpi=300, bbox_inches='tight')

# Plot
with plt.style.context(["science", "nature"]):
    fig, ax = plt.subplots()
    ax.autoscale(tight=True)
    ax.set(**pparam)
    for grp in g_unique:
        ax.plot(df.t[df.g==grp], df.L[df.g == grp], label=f"{grp}")
    ax.legend()
    fig.savefig('angular_momentum_t.png', dpi=300, bbox_inches='tight')
