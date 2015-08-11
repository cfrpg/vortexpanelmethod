#!/usrp/bin/python3

""" plotting.py

Created: 8/11/2015
Author: Michel Robijns

This file is part of vortexpanelmethod which is released under the MIT license.
See the file LICENSE or go to http://opensource.org/licenses/MIT for full
license details.

TODO: Add description
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.spines as spn
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from vpm import vortex_panel_method


def main():
    pass


if __name__ == '__main__':
    main()


def plot_cp(x, cp, x_airfoil, y_airfoil, use_TeX=False):
    """ TODO: Add description
    """
    
    blue = [0.118, 0.506, 0.984, 1]
    #cyan = [0.129, 0.996, 0.996, 1]
    #green = [0.118, 0.988, 0.133, 1]
    #orange = [0.992, 0.596, 0.118, 1]
    #red = [0.988, 0.047, 0.082, 1]
    
    plt.figure(num='Pressure Distribution', figsize=(6, 6), dpi=100)
    
    if use_TeX:
        plt.rc('font', **{'family': 'serif', 'serif': ['Palatino']})
        plt.rc('text', usetex=True)
    
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 0.2])
    gs.update(hspace=0.1) # Set the spacing between axes
    
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    
    ax1.plot(x, cp, color=blue)
    ax1.set_ylabel('$C_p$ $[-]$', fontsize=16)
    ax1.set_xlabel('$x/c$ $[-]$', fontsize=16)
    ax1.set_xlim(0, 1)
    ax1.set_ylim(1, -6)
    ax1.set_xticks(np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]))
    ax1.set_yticks(np.linspace(1, -6, 8))
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_position('zero')
    ax1.spines['bottom'].set_position('zero')
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')
    ax1.xaxis.labelpad = 20
    
    ax2.plot(x_airfoil, y_airfoil, color='k')
    ax2.set_xlim(-0.01, 1)
    ax2.set_ylim(-0.11, 0.11)
    ax2.axis('off')
    
    #plt.savefig('name.pdf', format = 'pdf', dpi = 1000, bbox_inches='tight')
    plt.show()


def plot_polar(coordinates, alpha_begin=-10, alpha_end=20, samples=31, use_TeX=False):
    """ TODO: Add description
    """
    
    blue = [0.118, 0.506, 0.984, 1]
    #cyan = [0.129, 0.996, 0.996, 1]
    #green = [0.118, 0.988, 0.133, 1]
    #orange = [0.992, 0.596, 0.118, 1]
    red = [0.988, 0.047, 0.082, 1]
    
    alpha_range = np.linspace(alpha_begin, alpha_end, samples)
    
    cl = np.zeros(np.size(alpha_range))
    cm = np.zeros(np.size(alpha_range))
    
    for i in range(np.size(alpha_range)):
        x, cp, cl[i], cm[i] = vortex_panel_method(alpha_range[i], coordinates)
    
    plt.figure(num='Lift Coefficient and Moment Coefficient versus Angle of Attack', figsize=(6, 6), dpi=100)
    ax = plt.gca()
    
    if use_TeX:
        plt.rc('font', **{'family': 'serif', 'serif': ['Palatino']})
        plt.rc('text', usetex=True)
    
    # Use no labels for minor ticks
    yLocator = MultipleLocator(0.1)
    majorFormatter = FormatStrFormatter('%d')
    xLocator = MultipleLocator(1)
    ax.yaxis.set_minor_locator(yLocator)
    ax.xaxis.set_minor_locator(xLocator)
    
    # Plot the data
    plt.plot(alpha_range, cl, color=blue, label=r'$C_l$')
    plt.plot(alpha_range, cm, color=red, label=r'$C_m$')
    plt.legend(loc='center right', frameon=False)
    plt.ylabel('$C_l / C_m$ $[-]$', fontsize=16)
    plt.xlabel(r'$\alpha$ $[deg]$', fontsize=16)
    plt.xlim(-8, 20)
    plt.ylim(-1.0, 2.5)
    plt.xticks(np.linspace(-8, 20, 8))
    plt.yticks(np.array([-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5]))
    
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # Only show ticks on the left and bottom spines
    ax.spines['left'].set_position('zero')
    ax.spines['bottom'].set_position('zero')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.xaxis.labelpad = 20
    
    # Hide the tick label at the origin
    yticks = ax.get_yticklabels()
    yticks[2].set_visible(False)
    xticks = ax.xaxis.get_major_ticks()
    xticks[2].set_visible(False)
    
    #plt.savefig('name.pdf', format = 'pdf', dpi = 1000, bbox_inches='tight')
    plt.show()
