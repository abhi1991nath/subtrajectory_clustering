import numpy as np
import sys

from base import *
from distanceUtils import *


class trajGrid:
    """ A class representing a simple grid to index trajectory points.
    
    Atrributes:
        data { (int,int): {int: [int]} }: dictionary mapping a grid cell
                to a dictionary containing traj points lying in the cell.
        startLat (float): smallest y-value of indexed points.
        startLon (float): smallest x-value of indexed points.
        delta (float): size of grid cell.
        numX (int): number of cells in the x-direction.
        numY (int): number of cells in the y-direction.
    """
    
    def __init__(self, data, startLat, startLon, delta, numX, numY):
        """ Simple initialization of class. """
        self.data = data
        self.startLat = startLat
        self.startLon = startLon
        self.delta = delta
        self.xCells = numX
        self.yCells = numY


def gridData(trajs, simpTrajs, delta):
    """ Create an index of trajectory points.
    
    Args:
        trajs ({int : traj}): dictionary of traj objects.
        simpTajs ([simpleTraj]): list of simplified trajectories whose points
                are to be indexed.
        delta (float): size of grid cell.
        
    Returns:
        trajGrid object containing the points feom the simplified trajectories.
    """
    
    # Compute the x and y extent of the points.
    
    simpTraj = simpTrajs[0]
    trID = simpTraj.trajID
    p = trajs[trID].pts[simpTraj.indices[0]]
    xMin, yMin = p.lon, p.lat
    xMax, yMax = xMin, yMin

    for simpTraj in simpTrajs:
        trID = simpTraj.trajID
        for index in simpTraj.indices:
            p = trajs[trID].pts[index]
            xMin, xMax = min(xMin, p.lon), max(xMax, p.lon)
            yMin, yMax = min(yMin, p.lat), max(yMax, p.lat)


    # Initialize grid which is a dictionary mapping from the cell (xIndex, yIndex)
    # to the contents of the cell.
    # grid is of the form {(xIndex, yIndex) : {trajID : [indices of points in simpTraj]}}.
    
    xCells = int(np.ceil((xMax-xMin)/delta))
    yCells = int(np.ceil((yMax-yMin)/delta))
    grid = {}
    for simpTraj in simpTrajs:
        trID = simpTraj.trajID
        for index in simpTraj.indices:
            p = trajs[trID].pts[index]
            xCell = int(np.floor((p.lon - xMin)/delta))
            yCell = int(np.floor((p.lat - yMin)/delta))
            # Make sure cell exists.
            if grid.has_key((xCell, yCell)) is False:
                grid[(xCell, yCell)] = {}
            # Add trajID to cell if not there.
            if grid[(xCell, yCell)].has_key(trID) is False:
                grid[(xCell, yCell)][trID] = []
            grid[(xCell, yCell)][trID].append(index)
    
    return trajGrid(grid, yMin, xMin, delta, xCells, yCells)
