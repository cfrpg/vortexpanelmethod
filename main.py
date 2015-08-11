#!/usrp/bin/python3

""" main.py

Created: 7/28/2015
Author: Michel Robijns

This file is part of vortexpanelmethod which is released under the MIT license.
See the file LICENSE or go to http://opensource.org/licenses/MIT for full
license details.

TODO: Add description
"""

from vpm import vortex_panel_method
import numpy as np
from naca import NACA4
from plotting import plot_cp
from plotting import plot_polar


def main():
    # This function runs the vortex panel method on a NACA 2412
    
    # Generate airfoil coordinates
    airfoil_coordinates = NACA4('2412', 50)
        
    x_airfoil = airfoil_coordinates[:, 0]
    y_airfoil = airfoil_coordinates[:, 1]
    
    # Set the angle of attack
    alpha = 10

    # Run the vortex panel method
    x, cp, cl, cm = vortex_panel_method(alpha, airfoil_coordinates)
    
    # Print the values of cl and cm in the terminal 
    print('C_l = ', cl)
    print('C_m = ', cm)
    
    # Plot the pressure distribution
    plot_cp(x, cp, x_airfoil, y_airfoil)
    
    # Plot the lift coefficient and moment coefficient versus angle of attack
    plot_polar(airfoil_coordinates)


if __name__ == "__main__":
    main()
