#!/usr/bin/env python
"""Read MAPS array location file and plot the UV distribution
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('array_file', type=str, help='array location file')
parser.add_argument('lat', type=float,
                    help='latitude of the center of the array [degree]')
parser.add_argument('ha', type=float,
                    help='hour angle of the pointing center [hour]')
parser.add_argument('dec', type=float,
                    help='declination of the pointing center [degree]')
parser.add_argument('wavelenght', type=float,
                    help='wavelenght of the observation [meter]')
parser.add_argument('-o', '--output', type=str, metavar='filename.type',
                    help='save the plot to a file')
args = parser.parse_args()

# convert arguments to radian
h0 = args.ha * 15. * np.pi / 180.0
d0 = args.dec * np.pi / 180.0
l0 = args.lat * np.pi / 180.0

E, N = np.genfromtxt(args.array_file, usecols=(0, 1), unpack=True)
x = - N * np.sin(l0)
y = E
z = N * np.cos(l0)

xx1, xx2 = np.meshgrid(x, x)
yy1, yy2 = np.meshgrid(y, y)
zz1, zz2 = np.meshgrid(z, z)

bx = np.ravel(xx2 - xx1)
by = np.ravel(yy2 - yy1)
bz = np.ravel(zz2 - zz1)

u = (np.sin(h0) * bx + np.cos(h0) * by) / args.wavelenght
v = (-np.sin(d0) * np.cos(h0) * bx + np.sin(d0) * np.sin(h0) * by
     + np.cos(d0) * bz) / args.wavelenght

lim = (np.min((u, v)) * 1.10, np.max((u, v)) * 1.10)
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(1, 1, 1, adjustable='box', aspect=1.0)
ax.scatter(u, v, s=1, c='k', marker='.')
ax.set_xlim(lim)
ax.set_ylim(lim)
ax.set_xlabel('$U\;[\lambda={:.2f}\,m]$'.format(args.wavelenght))
ax.set_ylabel('$V\;[\lambda={:.2f}\,m]$'.format(args.wavelenght))
plt.title('$UVwave\;(ha={0:.2f}h,\,\delta={1:.2f}d)$'.format(args.ha, args.dec))
if args.output is not None:
    plt.savefig(args.output)
plt.show()
