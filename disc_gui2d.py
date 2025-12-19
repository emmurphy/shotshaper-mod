# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 15:53:30 2022

@author: 2913452
"""

from shotshaper.projectile import DiscGolfDisc
import matplotlib.pyplot as pl
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider, Button, TextBox
import numpy as np
from shotshaper.transforms import T_21
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import argparse
import pint
import math

#create Parser
parser = argparse.ArgumentParser(description="Launch Shotshaper with specified disc data file and units.")

# Add an argument for unit selection with default value 'imperial'
parser.add_argument('--units', type=str, default='metric', choices=['metric', 'imperial'],
                    help='Unit system to use (metric or imperial). Default is metric. All calculations are done in metric and converted if desired.')

# Modify the argument for disc data file to expect a single disc name
parser.add_argument('disc', help='Name of the disc data file.')

#Mass of Disc
parser.add_argument('--mass', type=float, default=170.0, help='Mass of disc in grams, default is 170g')

# Parse the arguments
args = parser.parse_args()

name = args.disc
units = args.units

def disc_vertices(attitude):
    a,nt = attitude.shape
    nvert = 40
    r = np.linspace(0, 2*np.pi, nvert)
    xd = np.cos(r)
    yd = np.sin(r)
    zd = np.zeros(nvert)
    discoutline = np.vstack((xd, yd, zd)).T
    p=0
    while p < nvert:
        discoutline[p] = np.matmul(T_21(attitude[:,0]), discoutline[p]) # Convert outline from disc coords to ground coords
        p += 1
    d = [list(zip(discoutline[:,0],discoutline[:,1],discoutline[:,2]))]
    
    return d

#name = 'cd5'

#convert to kg
mass = args.mass / 1000.0;

ureg = pint.UnitRegistry()

if units == 'imperial':
    velo = 1.0 * ureg.miles / ureg.hour
    distance = 1.0 * ureg.feet

else:
    velo = 1.0 * ureg.meter / ureg.second
    distance = 1.0 * ureg.meter
    
rotation = 1.0 * ureg.rpm
angle = 1.0 * ureg.degree

d = DiscGolfDisc(name, mass=mass)
speed = 24
omega = d.empirical_spin(speed)
spinRate_rpm = d.empirical_spin(speed) * 60.0 / ( 2.0 * math.pi)
z0 = 1.3
pos = np.array((0,0,z0))
pitch = 10.0
nose = -3.0
roll = 0.0  
yaw = 0
adjust_axes = False

s = d.shoot(speed=speed, omega=omega, pitch=pitch, position=pos, nose_angle=nose, roll_angle=roll,yaw=yaw)

x,y,z = s.position * ureg.meter
    
if units == 'imperial':
    x = x.to( ureg.feet )
    y = y.to( ureg.feet )
    z = z.to( ureg.feet )
        
x = x.magnitude
y = y.magnitude
z = z.magnitude

# Creating figure
fig = pl.figure(1,figsize=(13, 6), dpi=80)
ax1 = pl.subplot(2,3,2)
ax2 = pl.subplot(2,3,5)
ax3 = pl.subplot(1,3,3)
fac = 1.6
ylim = fac*max(abs(min(y)),abs(max(y)))
ax1.axis((min(x),fac*max(x),-ylim,ylim))
ax2.axis((min(x),fac*max(x),min(z),fac*max(z)))
ax3.axis((-ylim,ylim,min(z),fac*max(z)))
        
ax3.invert_xaxis()

l1, = ax1.plot(x,y,lw=2)
ax1.set_xlabel('Distance ({})'.format( str(distance.units) ) )
ax1.set_ylabel('Drift ({})'.format( str(distance.units) ) )
l2, = ax2.plot(x,z,lw=2)
ax2.set_xlabel('Distance ({})'.format( str(distance.units) ) )
ax2.set_ylabel('Height ({})'.format( str(distance.units) ) )

l3, = ax3.plot(y,z,lw=2)
ax3.set_xlabel('Drift ({})'.format( str(distance.units) ) )
ax3.set_ylabel('Height ({})'.format( str(distance.units) ) )

xax = 0.07
ax4  = pl.axes([xax, 0.80, 0.25, 0.03], facecolor='lightgrey')
ax5  = pl.axes([xax, 0.75, 0.25, 0.03], facecolor='lightgrey')
ax6  = pl.axes([xax, 0.70, 0.25, 0.03], facecolor='lightgrey')
#ax7  = pl.axes([xax, 0.65, 0.25, 0.03], facecolor='lightgrey')
ax8  = pl.axes([xax, 0.65, 0.25, 0.03], facecolor='lightgrey')
ax9  = pl.axes([xax, 0.60, 0.25, 0.03], facecolor='lightgrey')
ax11 = pl.axes([xax, 0.55, 0.25, 0.03], facecolor='lightgrey')

if units == 'imperial':
    s1 = Slider(ax=ax4, label='Speed (mph)', valmin=1,  valmax=150, valinit=speed)
else:
    s1 = Slider(ax=ax4, label='Speed (m/s)', valmin=1,  valmax=70, valinit=speed)
    
s2 = Slider(ax=ax5, label='Roll (deg)',   valmin=-90, valmax=90, valinit=roll)
s3 = Slider(ax=ax6, label='Pitch (deg)', valmin=-10,   valmax=50, valinit=pitch)
s5 = Slider(ax=ax8, label='Nose (deg)',   valmin=-25, valmax=25, valinit=nose)
s7 = Slider(ax=ax9, label='Mass (g)',   valmin=140, valmax=200, valinit=mass)
s6 = Slider(ax=ax11, label='Spin (RPM)',   valmin=1, valmax=3000, valinit=spinRate_rpm)

def update(x):
    
    if units == 'imperial':
        speed = s1.val * ureg.miles / ureg.hours
        speed = speed.to( ureg.meter / ureg.second )
        speed = speed.magnitude
    else:
        speed = s1.val
    roll = s2.val
    pitch = s3.val
    nose = s5.val
    
    #convert to kg
    mass = s7.val / 1000
    d = DiscGolfDisc(name,mass=mass)
    
    #convert RPM to rad/sec
    omega = ( s6.val / 60.0 ) * 2 * math.pi
    s = d.shoot(speed=speed, omega=omega, pitch=pitch, position=pos, nose_angle=nose, roll_angle=roll)
    x,y,z = s.position * ureg.meter
    
    if units == 'imperial':
        x = x.to( ureg.feet )
        y = y.to( ureg.feet )
        z = z.to( ureg.feet )
        
    x = x.magnitude
    y = y.magnitude
    z = z.magnitude
        
  
    
    l1.set_xdata(x)
    l1.set_ydata(y)
    l2.set_xdata(x)
    l2.set_ydata(z)
    l3.set_xdata(y)
    l3.set_ydata(z)
    
    if adjust_axes:
        ax1.axis((min(x),max(x),min(y),max(y)))
        ax2.axis((min(x),max(x),min(z),max(z)))
        ax3.axis((min(y),max(y),min(z),max(z)))
        
    fig.canvas.draw_idle()

s1.on_changed(update)
s2.on_changed(update)
s3.on_changed(update)
s5.on_changed(update)
s6.on_changed(update)
s7.on_changed(update)

pl.show()
    
