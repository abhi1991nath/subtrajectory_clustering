import sys

from base import *
from frechet import sqDist


class simpleTraj:
    """ Class representing a 'simplified' (sub)trajectory.
    
    Instead of storing the actual points, stores the indices of points left
    over from the parent (sub)trajectory after simplification.
    
    Attributes:
        trajID: ID of parent (sub)trajectory.
        indices: indices of points left over from the parent trajectory. Indices
                 are based on the list of points sorted by timestamp.
    """
    
    def __init__(self):
        """ Basic initialization. """
        self.trajID = -1
        self.indices = []


def fSimplify(uTraj, tau, start, end):
    """ Simplify a (sub)trajectory in the forward direction.
    
    Traverse the points in ascending order of timestamps, only retaining a point if
    it is more than a certain distance away from the last retained point. The first
    and last points are always retained.
    
    Args:
        uTraj (traj): trajectory to be simplified (points are sorted by timestamp).
        tau (float): distance for simplification.
        start (int): starting index of trajectory point to be simplified.
        end (int): ending index of trajectory point to be simplified.
        
    Returns:
        A forward-simplified version of the subtrajectory lying between start and end.
    """
    
    sTraj = simpleTraj()
    sTraj.trajID = uTraj.pts[0].trajID
    sTraj.indices.append(start)
    curPoint = start
    for i in xrange(start+1, end+1):
        if i == end:
            sTraj.indices.append(i)
            continue
        if sqDist(uTraj.pts[curPoint], uTraj.pts[i]) >= tau**2 :
            sTraj.indices.append(i)
            curPoint = i

    return sTraj


def bSimplify(uTraj, tau, start, end):
    """ Simplify a (sub)trajectory in the reverse direction.
    
    Traverse the points in descending order of timestamps, only retaining a point if
    it is more than a certain distance away from the last retained point. The first
    and last points are always retained.
    
    Args:
        uTraj (traj): trajectory to be simplified (points are sorted by timestamp).
        tau (float): distance for simplification.
        start (int): starting index of trajectory point to be simplified.
        end (int): ending index of trajectory point to be simplified.
        
    Returns:
        A backward-simplified version of the subtrajectory lying between start and end.
    """
    
    sTraj = simpleTraj()
    sTraj.trajID = uTraj.pts[0].trajID
    sTraj.indices.append(end)
    curPoint = end
    for i in xrange(end-1, start-1, -1):
        if i == start:
            sTraj.indices.append(i)
            continue
        if sqDist(uTraj.pts[curPoint], uTraj.pts[i]) >= tau**2:
            sTraj.indices.append(i)
            curPoint = i
            
    sTraj.indices.reverse()
    return sTraj
    

def canonise(i,j):
    """ Return canonical intervals of the input interval.
    
    Args:
        i (int): left endpoint of input interval (inclusive).
        j (int): right endpoint of input interval (inclusive).
        
    Returns:
        List of canonical intervals. Each interval is of the form (start,end)
        where start and end are integers.
    """
    
    if j-i <= 8:  # Hard-coded lower bound on length of canonical pathlets.
        return []
    retval = []
    retval.append((i, j))
    midpoint = (i + j)/2
    retval.extend(canonise(i, midpoint))
    retval.extend(canonise(midpoint + 1, j))
    return retval


# test canonise function
if __name__ == "__main__":
    print canonise(0, 16)
