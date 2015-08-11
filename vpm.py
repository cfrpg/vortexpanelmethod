#!/usr/bin/python3

""" vpm.py

Created: 12/28/2014
Author: Michel Robijns

This file is part of vortexpanelmethod which is released under the MIT license.
See the file LICENSE or go to http://opensource.org/licenses/MIT for full
license details.

TODO: Add description
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt
import vpmhelperfunctions as pm


def main():
    pass


if __name__ == '__main__':
    main()


def vortex_panel_method(alpha, coordinates):
    """ TODO: Add description
    """
    
    # Convert the angle of attack from degrees to radians
    alpha = np.radians(alpha)
    
    # Compute the velocity unit vector
    Q = np.array([cos(alpha), sin(alpha)]);
    
    # N represents the number of airfoil coordinates
    N = np.max(np.shape(coordinates))
    
    # M represents the number of panels
    M = N - 1
    
    # Change the panel sequence to clockwise starting at the trailing edge.
    # This is the opposite of the convention for airfoil data files.
    coordinates = np.flipud(coordinates)
    
    # Generate two matrices containing the coordinates of the respective edges
    # of each individual panel. Each matrix has two columns containing x and y
    # coordinates, respectively. The n-th row of the matrix corresponds to the
    # n-th panel.
    panelBegins, panelEnds = pm.generatePanels(coordinates)
    
    # Compute the angle of each panel with respect to the horizontal.
    panelAngles = pm.computePanelAngles(panelBegins, panelEnds)
    
    # Compute each panel's normal vector
    normalVectors = pm.computeNormalVectors(panelAngles)
    
    # Compute each panel's tangent vector
    tangentVectors = pm.computeTangentVectors(panelAngles)
    
    # Compute the collocation (midpoint) of each panel
    panelCollocations = pm.computePanelCollocations(panelBegins, panelEnds)
    
    # Preallocate the vectors and matrices in which the results will be stored
    panelLengths = np.zeros((M, 1), dtype=float)
    A = np.zeros((N, N), dtype=float)
    B = np.zeros((M, N), dtype=float)
    RHS = np.zeros((N, 1), dtype=float)
    
    # Create the influence matrices
    for i in np.arange(M):
        for j in np.arange(M):
            # We are now looking at each individual panel. Convert the
            # global reference frame to a local panel reference frame.
            # The beginning of a panel is always the origin.
            
            x2 = panelEnds[j, 0] - panelBegins[j, 0]
            z2 = panelEnds[j, 1] - panelBegins[j, 1]
            end = np.array([x2, z2])
            x = panelCollocations[i, 0] - panelBegins[j, 0]
            z = panelCollocations[i, 1] - panelBegins[j, 1]
            collocation = np.array([x, z])
                        
            # Convert the local reference frame to horizontal because
            # ONLY THEN equations 11.44 and 11.45 are valid
            end = pm.rotate2DVector(end, panelAngles[j])
            collocation = pm.rotate2DVector(collocation, panelAngles[j])
            
            x2 = end[0]
            z2 = end[1]
            x = collocation[0]
            z = collocation[1]
                        
            # Store panel lengths
            if i == 0:
                panelLengths[j] = x2
            
            # Compute velocities
            if (i == j):
                V1 = np.array([-0.5 * (x - x2) / x2, -1 / (2 * pi)])
                V2 = np.array([0.5 * x / x2, 1 / (2 * pi)])
            else:
                V1 = pm.velocityDueToLinearVortexPanelBeginEdge(x, x2, z)
                V2 = pm.velocityDueToLinearVortexPanelEndEdge(x, x2, z)
            
            # Convert the local reference frame back to the global
            # reference frame
            V1 = pm.rotate2DVector(V1, -panelAngles[j])
            V2 = pm.rotate2DVector(V2, -panelAngles[j])
            
            if j == 0:
                A[i, 0] = np.dot(V1, normalVectors[i, :])
                tempA = np.dot(V2, normalVectors[i, :])
                B[i, 0] = np.dot(V1, tangentVectors[i, :])
                tempB = np.dot(V2, tangentVectors[i, :])
            elif j == M - 1:
                A[i, M - 1] = np.dot(V1, normalVectors[i, :]) + tempA;
                A[i, N - 1] = np.dot(V2, normalVectors[i, :])
                B[i, M - 1] = np.dot(V1, tangentVectors[i, :]) + tempB;
                B[i, N - 1] = np.dot(V2, tangentVectors[i, :])
            else:
                A[i, j] = np.dot(V1, normalVectors[i, :]) + tempA
                tempA = np.dot(V2, normalVectors[i, :])
                B[i, j] = np.dot(V1, tangentVectors[i, :]) + tempB
                tempB = np.dot(V2, tangentVectors[i, :])
        
        RHS[i] = np.dot(-Q, normalVectors[i, :]);
    
    # Add the Kutta condition
    A[N - 1, 0] = 1
    A[N - 1, N - 1] = 1
    RHS[N - 1] = 0
    
    # Solve the system of equations
    vortexStrength = np.linalg.solve(A, RHS)

    cl = 0
    cp = np.zeros((M,1), dtype=float)
    mx = 0
    mz = 0
    
    unit_l = pm.rotate2DVector(Q, -pi / 2)
    
    unit_z = np.array([0, 1])
    unit_x = np.array([1, 0])
        
    for i in np.arange(M):
        # Compute the sum of rows in B
        sum = 0
        for j in np.arange(M):
            sum += B[i, j] * vortexStrength[j]

        tangentialVelocity = sum + np.dot(Q, tangentVectors[i, :])
        
        cp[i] = 1 - tangentialVelocity**2
        cl += np.dot(cp[i] * panelLengths[i] * -1 * normalVectors[i, :], unit_l)
        
        p = cp[i] * panelLengths[i] * -1 * normalVectors[i, :]
        
        mx += np.dot(p, -unit_z) * (panelCollocations[i, 0] - 0.25)
        
        if panelCollocations[i, 1] < 0:
            mz += np.dot(p, -unit_x) * panelCollocations[i, 1]
        else:
            mz += np.dot(p, unit_x) * panelCollocations[i, 1]
    
    cl *= cos(alpha)
    cm = mx + mz
            
    return panelCollocations[:, 0], cp, cl, cm
