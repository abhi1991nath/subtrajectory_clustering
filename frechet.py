import numpy as np
import pickle
import sys

from base import *
from sets import Set


def sqDist(pt1, pt2):
    """ Compute squared Euclidean distance b/w points.
    
    Args:
        pt1, pt2 (pt): Input points.
        
    Returns:
        Squared Euclidean distance b/w pt1 and pt2.
    """
    
    return (pt1.lat - pt2.lat)**2 + (pt1.lon - pt2.lon)**2


def distPtSegment(p, seg):
    """ Compute the distance b/w a point and a line segment.
    
    The distance is defined as the shortest over all distances b/w the point
    and all the points of the segment.
    
    Args:
        p (pt): Input point.
        seg (pt,pt): Endpoints of input segment.
        
    Returns:
        The distance b/w p and seg.
    """
    
    (q1, q2) = seg
    x,y,x1,y1,x2,y2 = p.lon, p.lat, q1.lon, q1.lat, q2.lon, q2.lat
    
    if x1==x2 and y1==y2: # Degenerate segment.
        return float(np.sqrt(sqDist(p,q1)))
        
    if x1==x2: # Vertical segment.
        if y1 <= y <= y2 or y2 <= y <= y1:
            return abs(x-x1)
        else:
            return min(float(np.sqrt(sqDist(p,q1))), float(np.sqrt(sqDist(p,q2))))
            
    elif y1==y2: # Horizontal segment.
        if x1 <= x <= x2 or x2 <= x <= x1:
            return abs(y-y1)
        else:
            return min(float(np.sqrt(sqDist(p,q1))), float(np.sqrt(sqDist(p,q2))))
            
    else:
        x,y,x2,y2 = x-x1,y-y1,x2-x1,y2-y1 # Translate so that (x1,y1) is at the origin.
        m = y2/x2
        c = y + x/m
        x3, y3 = c/(m+1/m), m*c/(m+1/m) # Projection of (x,y) on line passing through origin and (x2,y2).
        if x2*x3 + y2*y3 >= 0: # (x3,y3) between origin and (x2,y2).
            return float(np.sqrt((x3-x)**2 + (y3-y)**2))
        else:
            return min(float(np.sqrt(sqDist(p,q1))), float(np.sqrt(sqDist(p,q2))))
            

def frechetDec(trajA, trajB, delta):
    """ Decide if the discrete Frechet distance b/w trajectories is at most delta.
    
    Uses the classic dynamic programming algorithm to find the optimal correspondence.
    
    Args:
        trajA, trajB ([pt]): input trajectories, represented as list of pt objects.
        delta (float): guess on the discrete Frechet distance.
        
    Returns:
        True if the discrete Frechet distance b/w trajA and trajB is at most delta,
        False otherwise.
    """
    
    ptQueue = [(0, 0)]
    visited = Set([])
    while len(ptQueue) > 0:
        current = ptQueue.pop(0)
        if current in visited:
            continue
        if current == (len(trajA) - 1, len(trajB) - 1):
            return True
        visited.add(current)
        i = current[0]
        j = current[1]
        #bounds check, add points which are within delta (using squared distance)
        if i + 1 < len(trajA):
            if sqDist(trajA[i + 1], trajB[j]) <= delta**2:
                ptQueue.append((i+1, j))
        if j + 1 < len(trajB):
            if sqDist(trajA[i], trajB[j + 1]) <= delta**2:
                ptQueue.append((i, j+1))
        if i + 1 < len(trajA) and j + 1 < len(trajB):
            if sqDist(trajA[i + 1], trajB[j + 1]) <= delta**2:
                ptQueue.append((i+1, j+1))
    return False


def semiContFrechetDec(trajA, trajB, delta):
    """ Decide if the semi-contiuous Frechet distance b/w trajectories is at most delta.
    
    The semi-continuous Frechet distance is defined similarly to its discrete cousin,
    except that the correspondences are defined b/w points on one trajectory and segments
    on the other. The cost of a correspondence pair is calculated as the distance b/w the
    point and segment. The decision procedure is based on dynamic programming.
    
    Args:
        trajA, trajB ([pt]): input trajectories, represented as list of pt objects.
        delta (float): guess on the discrete Frechet distance.
        
    Returns:
        True if the semi-continuous Frechet distance b/w trajA and trajB is at most delta,
        False otherwise.
    """

    trajA.append(trajA[-1])
    trajB.append(trajB[-1])
    
    ptQueue = [(0,0)]
    visited = set()
    while len(ptQueue) > 0:
        current = ptQueue.pop(0)
        if current in visited:
            continue
        if current == (len(trajA)-1, len(trajB)-1):
            return True
        visited.add(current)
        i, j = current[0], current[1]
        if i+1 < len(trajA):
            if j != 0:
                seg1, seg2 = (trajA[i], trajA[i+1]), (trajB[j-1], trajB[j])
                if max(distPtSegment(trajA[i+1],seg2), distPtSegment(trajB[j], seg1)) <= delta:
                    ptQueue.append((i+1,j))
            else:
                if sqDist(trajA[i+1],trajB[j]) <= delta**2:
                    ptQueue.append((i+1,j))
        if j+1 < len(trajB):
            if i != 0:
                seg1, seg2 = (trajA[i-1], trajA[i]), (trajB[j], trajB[j+1])
                if max(distPtSegment(trajB[j+1], seg1), distPtSegment(trajA[i], seg2)) <= delta:
                    ptQueue.append((i,j+1))
            else:
                if sqDist(trajA[i], trajB[j+1]) <= delta**2:
                    ptQueue.append((i,j+1))
        if i+1 < len(trajA) and j+1 < len(trajB):
            seg1, seg2 = (trajA[i], trajA[i+1]), (trajB[j], trajB[j+1])
            if max(distPtSegment(trajA[i+1], seg2), distPtSegment(trajB[j+1], seg1)) <= delta:
                ptQueue.append((i+1,j+1))
    return False
   

## test decision functions
if __name__ == "__main__":
    trajs = pickle.load(open("manhattan_grid_100_100_5.txt.CleanedUp.p", "rb"))
#    for trajID in trajs:
#        tr = trajs[trajID]
#        pts = densify(tr.pts, 0.0001)
#        for i in xrange(len(pts) - 1):
#            p1, p2 = pts[i], pts[i+1]
#            if float(np.sqrt(sqDist(p1,p2))) > 0.0001:
#                print "Error %f" % float(np.sqrt(sqDist(p1,p2)))
        
    #trajs = readTrajsFromTxtFile(sys.argv[1])
    k = list(trajs.keys())
    traj1, traj2 = trajs[k[0]], trajs[k[1]]
    print semiContFrechetDec(traj1.pts, traj2.pts, Decimal(sys.argv[1]))
