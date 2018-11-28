import cPickle as pickle
import numpy as np
import sys

from base import *
from distanceUtils import *
from frechet import *
from grid import *


def inSquare(center, point, r):
    """ Check if a point lies inside an axis-aligned square.
    
    Args:
        center (pt): center of the square.
        point (pt): point to be checked.
        r (float): half side lengtg of square.
        
    Returns:
        True if point lies inside the square, False otherwise. 
    """
    
    return abs(center.lat - point.lat) <= r and abs(center.lon - point.lon) <= r


def findNeighboringPts(grid, point):
    """ Find points in grid that are 'close' to the input point.
    
    Compute all points in the (at most 9) neighboring cells of the grid cell
    containing the input point (inluding the input point's cell).
    
    Args:
        grid (trajGrid): grid containing trajectory points.
        point (pt): input point.
        
    Returns:
        Dictionay of the form {int : [int]}, mapping trajID to the indices of 
        its points that lie in the neighboring cells of input point's cell.
    """
    
    # Find grid cell containing input point.
    xCell = int(np.floor((point.lon - grid.startLon)/grid.delta))
    yCell = int(np.floor((point.lat - grid.startLat)/grid.delta))

    possibleTrajs = {}
    for i in range(-1, 2):
        for j in range(-1, 2):
            if 0 <= xCell + i < grid.xCells and 0 <= yCell + j < grid.yCells: # Boundary check
                # If there are no pts in the cell it won't be in data dictionary.
                if grid.data.has_key((xCell+i, yCell+j)):
                    for trajID in grid.data[(xCell+i, yCell+j)]:
                        # If this is the first time you're adding points from this trajectory.
                        if possibleTrajs.has_key(trajID) is False:
                            possibleTrajs[trajID] = []
                        possibleTrajs[trajID].extend(grid.data[(xCell+i, yCell+j)][trajID])
    return possibleTrajs


def computeDistances(trajs, simpTraj, pth, r, possStarts, possEnds, upth, pathDic):
    """ Compute all subtrajectories, with both endpoints coming from a given set of
        points, that have small Frechet distance to a pathlet.
        
        Args:
            trajs ({int : traj}): dictionary mapping IDs to trajectories.
            simpTraj (simpleTraj): simplified trajectory whose subtrajectories are
                    to be checked for proximity to input pathlet.
            pth (simpleTraj): input simplified pathlet.
            r (float): guess on the Frechet distance.
            possStarts ([int]): indices of possible starting positions for subtrajectories
                    of simpTraj.
            possEnds ([int]): indices of possible ending positions for subtrajectories
                    of simpTraj.
            upth (pathlet): original unsimplified pathlet.
            pathDic ({(pathlet, subTraj) : float}): dictionary mapping pathlet-subtrajectory
                    pair to the approximate Frechet distance b/w them.
                    
        Returns:
            Nothing is returned. However, pathDic is updated to store the newly computed
            distances.
            
    """
    
    if len(possStarts)==0 or len(possEnds)==0:
        return
        
    possStarts.sort()
    possEnds.sort()

    # Extract pathlet points.
    pathletPoints = []
    for i in xrange(len(pth.indices)):
        currentIndex = pth.indices[i]
        pathletPoints.append(trajs[pth.trajID].pts[currentIndex])

    subTrajMatches = []
    # Keeps track of already assigned subtrajectories.
    startSkips = []
    
    # Indices of points in simpTraj.
    indices = simpTraj.indices

    # Compute possible assignable subtrajectories to the pathlet.
    curBound = possStarts[0]
    
    # Loop through all possible starting and ending points
    for i in xrange(len(possStarts)):
        for j in xrange(len(possEnds)-1, -1 , -1):
            if possEnds[j] < possStarts[i]:
                break
            if possEnds[j] <= curBound:
                break
                
            # If (upth, straj) is already present in pathDic, do not add it.
            if pathDic.has_key((upth, subTraj(simpTraj.trajID, (possStarts[i], possEnds[j])))) is True:
                curBound = possEnds[j]
                break
      
            # Pull out the actual points of the "simplified" subtrajectory.
            simpTrajPoints = []
            for k in xrange(len(indices)):
                if possStarts[i] <= indices[k] <= possEnds[j]:
                    simpTrajPoints.append(trajs[simpTraj.trajID].pts[indices[k]])
                
            # Check if subtrajectory is within distance 2r to pathlet, if yes then add it to pathDic.
            if semiContFrechetDec(simpTrajPoints, pathletPoints, 2*r) is True:
                pathDic[(upth, subTraj(simpTraj.trajID, (possStarts[i], possEnds[j])))] = r
                curBound = possEnds[j]
                break


def process(trajs, rmin, rmax):
    """ Compute approximate Frechet distances b/w subtrajectories and pathlets.
    
    These distances are used as inputs to the greedy algorithm. We avoid computing
    all pairwise distances by overlaying a grid over the points, and discarding
    pathlet-subtrajectiry pairs whose endpoints do not lie in neighboring cells.
    
    Args:
        trajs ({int: traj}): dictionary of trajectories.
        rmin, rmax (float): lower and upper bounds on the Frechet distances.
        
    Returns:
        A dictionary of the form {(pathlet, subTraj) : float} storing the distance
        b/w pathlet and subtrajectory.
    """

    r = rmin
    pathDic = {}
    pathlets = []

    # Create "unsimplified" canonical pathlets.
    
    for trID, tr in trajs.iteritems():
        for b in canonise(0, len(tr.pts) - 1):
            pathlets.append(pathlet(trID, (b[0], b[1])))

    # Compute pairs with distance O(r) (we lose a constant factor
    # due to simplification and exponential search on the possible
    # Frechet distances, but that's OK in the grand scheme of things.
    
    while r <= rmax:
        fSimpTrajs, bSimpTrajs = {}, {}
        
        # (r/2)-simplify trajectories, store in list.
        for trID, tr in trajs.iteritems():
            fTraj = fSimplify(tr, r/2, 0, len(tr.pts)-1)
            bTraj = bSimplify(tr, r/2, 0, len(tr.pts)-1)
            fSimpTrajs[trID] = fTraj
            bSimpTrajs[trID] = bTraj
            
        # Create grid of resolution r for fSimpTrajs and bSimpTrajs.
        fGrid = gridData(trajs, list(fSimpTrajs.values()), r)
        bGrid = gridData(trajs, list(bSimpTrajs.values()), r)
        
        # Find candidate subtrajectories for each pathlet. These are
        # subtrajectories whose endpoints are close to the pathlet's endpoints.
        
        for pth in pathlets:
            # Simplify pathlet.
            forwardPathlet = fSimplify(trajs[pth.trajID], r/2, pth.bounds[0], pth.bounds[1])
            backwardPathlet = bSimplify(trajs[pth.trajID], r/2, pth.bounds[0], pth.bounds[1])

            # Forward simplified subtrajectories with endpoints close to those of pth.
            startFNeighbors = findNeighboringPts(fGrid, trajs[pth.trajID].pts[pth.bounds[0]])
            endFNeighbors = findNeighboringPts(fGrid, trajs[pth.trajID].pts[pth.bounds[1]])
            possibleFSubTrajs = list(set(startFNeighbors.keys()) & set(endFNeighbors.keys()))
            
            # Backward simplified subtrajectories with endpoints close to those of pth.
            startBNeighbors = findNeighboringPts(bGrid, trajs[pth.trajID].pts[pth.bounds[0]])
            endBNeighbors = findNeighboringPts(bGrid, trajs[pth.trajID].pts[pth.bounds[1]])
            possibleBSubTrajs = list(set(startBNeighbors.keys()) & set(endBNeighbors.keys()))

            # Record distance between pathlet and itself to be 0.
            pairedTraj = subTraj(pth.trajID, pth.bounds)
            if pathDic.has_key((pth, pairedTraj)) is False:
                pathDic[(pth, pairedTraj)] = 0
                
            # Compute distance between forwardPathlet and possibleFSubTrajs.
            for trajID in possibleFSubTrajs:
                simpTraj = fSimpTrajs[trajID]
                computeDistances(trajs, simpTraj, forwardPathlet, r, startFNeighbors[trajID], endFNeighbors[trajID], pth, pathDic)
                
            # Compute distance between backwardPathlet and possibleBSubTrajs.
            for trajID in possibleBSubTrajs:
                simpTraj = bSimpTrajs[trajID]
                computeDistances(trajs, simpTraj, backwardPathlet, r, startBNeighbors[trajID], endBNeighbors[trajID], pth, pathDic)
                
        r*=2
        
    return pathDic

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "Wrong command. See README for more details."
        
    else:    
        trajs = pickle.load(open(sys.argv[1] + ".CleanedUp.p", "rb"))
        rmin, rmax = float(sys.argv[2]), float(sys.argv[3])
        distPairs = process(trajs, rmin, rmax)
#        for key, value in distPairs.iteritems():
#            print key[0], key[1], value
        pickle.dump(distPairs, open(sys.argv[1] +".distPairs1.sc.p", "wb"))
#        print len(trajs.keys())
#        print len(distPairs.keys())
#        # Plot the pathlets
#        pths = [x[0] for x in distPairs.keys()]
#        plotPathlets(pths, trajs)
#        pickle.dump(distPairs, open(sys.argv[1] +".distPairs1.p", "wb"))
