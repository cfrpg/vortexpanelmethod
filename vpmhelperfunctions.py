#!/usr/bin/python3

""" vpmhelperfunctions.py

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


def main():
    pass


if __name__ == '__main__':
    main()


def generatePanels(coordinates):
    """ TODO: Add a description 
    """
    
    M = np.max(np.shape(coordinates)) - 1
    
    panelBegins = np.zeros((M, 2), dtype=float)
    panelEnds = np.zeros((M, 2), dtype=float)
    
    for i in np.arange(M):
        panelBegins[i, 0] = coordinates[i, 0];
        panelEnds[i, 0] = coordinates[i + 1, 0];
        panelBegins[i, 1] = coordinates[i, 1];
        panelEnds[i, 1] = coordinates[i + 1, 1];
    
    return panelBegins, panelEnds


def computePanelAngles(panelBegins, panelEnds):
    """ TODO: Add a description 
    """
    
    M = np.max(np.shape(panelBegins))
    
    panelAngles = np.zeros((M, 1), dtype=float)
    
    for i in np.arange(M):
        dz = panelEnds[i, 1] - panelBegins[i, 1];
        dx = panelEnds[i, 0] - panelBegins[i, 0];
        panelAngles[i] = atan2(dz, dx);
        
    return panelAngles


def computeNormalVectors(panelAngles):
    """ TODO: Add a description 
    """
    
    M = np.max(np.shape(panelAngles))
    
    normalVectors = np.zeros((M, 2), dtype=float)
    
    for i in np.arange(M):
        nx = -sin(panelAngles[i]);
        ny = cos(panelAngles[i]);
        normalVectors[i, :] = [nx, ny]
        
    return normalVectors


def computeTangentVectors(panelAngles):
    """ TODO: Add a description 
    """
    
    M = np.max(np.shape(panelAngles))
    
    tangentVectors = np.zeros((M, 2), dtype=float)
    
    for i in np.arange(M):
        tx = cos(panelAngles[i]);
        ty = sin(panelAngles[i]);
        tangentVectors[i, :] = [tx, ty]
        
    return tangentVectors


def computePanelCollocations(panelBegins, panelEnds):
    """ TODO: Add a description 
    """
    
    M = np.max(np.shape(panelBegins))
    
    panelCollocations = np.zeros((M, 2), dtype=float)
    
    for i in np.arange(M):
        cx = (panelEnds[i, 0] - panelBegins[i, 0]) / 2 + panelBegins[i, 0];
        cy = (panelEnds[i, 1] - panelBegins[i, 1]) / 2 + panelBegins[i, 1];
        panelCollocations[i, :] = [cx, cy]
        
    return panelCollocations


def rotate2DVector(vector, angle):
    """ TODO: Add a description 
    """
    
    rotationMatrix = np.array([[cos(angle), -sin(angle)],
                               [sin(angle), cos(angle)]])
    
    rotatedVector = np.dot(vector, rotationMatrix);
    
    return rotatedVector


def normalVelocityDueToVortex(gamma, x, x1, x2, z, z1, z2):
    """ TODO: Add a description 
    """
    
    w = -gamma / (4 * pi) * log(((x - x1)**2 + (z - z1)**2) / ((x - x2)**2 +
        (z - z2)**2))
    
    return w


def tangentialVelocityDueToVortex(gamma, x, x1, x2, z, z1, z2):
    """ TODO: Add a description 
    """
    
    u = gamma / (2 * pi) * (atan2((z - z2), (x - x2))
        - atan2((z - z1), (x - x1)))
    
    return u


def velocityDueToLinearVortexPanelBeginEdge(x, x2, z):
    """ TODO: Add a description 
    """
        
    r1 = sqrt(x**2 + z**2)
    r2 = sqrt((x - x2)**2 + z**2)
    
    theta1 = atan2(z, x)
    theta2 = atan2(z, (x - x2))
    
    u1 = -1 * (z * log(r2 / r1) + x * (theta2 - theta1) - x2 * (theta2 - theta1)) / (2 * pi * x2)
    w1 = -1 * ((x2 - z * (theta2 - theta1)) - x * log(r1 / r2) + x2 * log(r1 / r2)) / (2 * pi * x2)
    
    return np.array([u1, w1])


def velocityDueToLinearVortexPanelEndEdge(x, x2, z):
    """ TODO: Add a description 
    """
    
    r1 = sqrt(x**2 + z**2)
    r2 = sqrt((x - x2)**2 + z**2)
    
    theta1 = atan2(z, x)
    theta2 = atan2(z, (x - x2))
    
    u2 = (z * log(r2 / r1) + x * (theta2 - theta1)) / (2 * pi * x2)
    w2 = ((x2 - z * (theta2 - theta1)) - x * log(r1 / r2)) / (2 * pi * x2)
    
    return np.array([u2, w2])
