import sys
import cPickle as pickle

from base import *
from heapdict import heapdict


def preprocessGreedy(trajs, distPairs):
    """ Compute pre-requisite data structures for the greedy algorithm.

        Args:
            trajs ({int : traj}): dict mapping ID to traj objects.
            distPairs ({(pathlet, int) : [(subTraj, float)]): dict mapping a
                    pathlet-trajID pair to a list of subtraj-float pairs, where the
                    subtraj belong to the resp. traj, and the float values are computed
                    Frechet distances.
            
        Returns:
            A 4-tuple (strajCov, ptStraj, strajPth, trajCov), where
            strajCov ({subtraj : int}) : dict storing coverage (#points) in all subtrajs.
                    in distPairs.
            ptStraj ({pt : {subtraj}}) : dict storing for each point, the set of subtrajs.
                    in distPairs containing it.
            strajPth ({subtraj: [pathlet]) : dict storing for each pathlet in distPairs,
                    the list of pathlets associated with it.
            trajCov ({int : int}) : dict storing the #points in each trajectory.
    """
    
    strajCov, ptStraj, strajPth = {}, {}, {}
    for key, value in distPairs.iteritems():
        pth, trID = key[0], key[1]
        for i in xrange(len(value)):
            straj, dist = value[i][0], value[i][1]
            strajCov[straj] = straj.bounds[1] - straj.bounds[0] + 1
            if strajPth.has_key(straj):
                strajPth[straj].append(pth)
            else:
                strajPth[straj] = [pth]
            for j in xrange(straj.bounds[0], straj.bounds[1]+1):
                p = trajs[trID].pts[j]
                if ptStraj.has_key(p):
                    ptStraj[p].add(straj)
                else:
                    ptStraj[p] = set([straj])
                    
    trajCov = {}  
    for trID, tra in trajs.iteritems():
        trajCov[trID] = len(tra.pts)
        for i in xrange(len(tra.pts)):
            p = tra.pts[i]
            if ptStraj.has_key(p) is False: # this means p can't be assigned to any pathlet
                ptStraj[p] = set()
    
    return (strajCov, ptStraj, strajPth, trajCov)


def processPoint(p, ptStraj, strajCov, strajPth, trajs, trajCov, distPairs, numUnprocessedPts, queue):
    """ Process a point picked in an interation of the greedy algorithm.
    
        Args:
            p (pt): Point to be processed.
            ptStraj ({pt : {subtraj}}) : dict storing for each point, the set of subtrajs
                    in distPairs containing it.
            strajCov ({subtraj : int}) : dict storing coverage (#points) in all subtrajs
                    in distPairs.
            strajPth ({subtraj: [pathlet]) : dict storing for each subtraj in distPairs,
                    the list of pathlets associated with it.
            trajs ({int : traj}) : dict mapping IDs to trajectories.
            trajCov ({int : int}) : dict storing the #points in each trajectory.
            distPairs ({(pathlet, int) : [(subTraj, float)]): dict mapping a
                    pathlet-trajID pair to a list of subtraj-float pairs, where the
                    subtraj belong to the resp. traj, and the float values are computed
                    Frechet distances.
            numUnprocessedPts (int) : no. of points left to be processed.
            queue : priority queue
            
        Returns:
            Set of the form {(pathlet, int}) containing pathlet-trajID pairs that the point
            can be assigned to.
    """
    
    retVal = set()
    for straj in ptStraj[p]:
        trID = straj.trajID
        strajCov[straj] = strajCov[straj]-1
        for i in xrange(len(strajPth[straj])):
            pth = strajPth[straj][i]
            retVal.add((pth, trID))
    ptStraj[p] = None # this is also marking that p is processed.
    numUnprocessedPts[0] = numUnprocessedPts[0] - 1

    # Check if we also have singleton sets.
    if queue is not None:
        trajCov[p.trajID] -= 1
        # Change priority of trajID to zero if its coverage is 0.
        if trajCov[p.trajID] == 0:
            queue[p.trajID] = 0
    return retVal


def processSubtraj(straj, strajCov, trajs, trajCov, ptStraj, strajPth, distPairs, numUnprocessedPts, queue):
    """ Process the points of a subtrajectory in an iteration of the greedy algorithm.
    
        Args:
            straj (subtraj): subtrajectory whose points are to be processed.
            strajCov ({subtraj : int}) : dict storing coverage (#points) in all subtrajs
                    in distPairs.
            trajs ({int : traj}) : dict mapping IDs to trajectories.
            trajCov ({int : int}) : dict storing the #points in each trajectory.
            ptStraj ({pt : {subtraj}}) : dict storing for each point, the set of subtrajs
                    in distPairs containing it.           
            strajPth ({subtraj: [pathlet]) : dict storing for each subtraj in distPairs,
                    the list of pathlets associated with it.
            distPairs ({(pathlet, int) : [(subTraj, float)]): dict mapping a
                    pathlet-trajID pair to a list of subtraj-float pairs, where the
                    subtraj belong to the resp. traj, and the float values are computed
                    Frechet distances.
            numUnprocessedPts (int) : no. of points left to be processed.
            queue : priority queue
            
        Returns:
            Set of the form {(pathlet, int}) containing pathlet-trajID pairs that the point
            can be assigned to.
    """
    trID = straj.trajID
    tr = trajs[trID]
    retVal = set()
    for i in xrange(straj.bounds[0], straj.bounds[1]+1):
        p = tr.pts[i]
        # Check if point p is unprocessed.
        if ptStraj[p] is not None:
            retVal = retVal.union(processPoint(p, ptStraj, strajCov, strajPth, trajs, trajCov, distPairs, numUnprocessedPts, queue))
    if strajCov[straj] != 0:
        print "Error!! Coverage should have been 0, instead of %d" %strajCov[straj]
        
    return retVal


def processTraj(trID, ptStraj, strajCov, strajPth, trajs, trajCov, distPairs, numUnprocessedPts, queue, unassignedPts):
    """ Process the unprocessed points of a trajectory.
    
        This function is called when an interation of the greedy algorithm decides to leave the unprocessed
        points of the trajectory unassigned.
        
        Args:
            trID (int): ID of the traj whose points are to be processed.
            ptStraj ({pt : {subtraj}}) : dict storing for each point, the set of subtrajs
                    in distPairs containing it.
            strajCov ({subtraj : int}) : dict storing coverage (#points) in all subtrajs
                    in distPairs.
            strajPth ({subtraj: [pathlet]) : dict storing for each subtraj in distPairs,
                    the list of pathlets associated with it.
            trajs ({int : traj}) : dict mapping IDs to trajectories.
            trajCov ({int : int}) : dict storing the #points in each trajectory.
            distPairs ({(pathlet, int) : [(subTraj, float)]): Dictionary mapping a
                    pathlet-trajID pair to a list of subtraj-float pairs, where the
                    subtraj belong to the resp. traj, and the float values are computed
                    Frechet distances.
            numUnprocessedPts (int) : no. of points left to be processed.
            queue : priority queue
                
        Returns:
            Set of the form {(pathlet, int}) containing pathlet-trajID pairs that the point
            can be assigned to.
    """
    points = trajs[trID].pts
    retVal = set()
    for i in xrange(len(points)):
        p = points[i]
        # Check if p is unprocessed.
        if ptStraj[p] is not None:
            # Add p to the list of unassigned points.
            unassignedPts.append(p)
            retVal = retVal.union(processPoint(p, ptStraj, strajCov, strajPth, trajs, trajCov, distPairs, numUnprocessedPts, queue))
            
    if trajCov[trID] != 0:
        print "Error!! Coverage should have been zero."
        
    return retVal
        

def computeCovCostRatio(trajStrajDist, c1, c3, strajCov):
    """ Compute coverage-cost ratio for a pathlet.
    
        Args:
            trajStrajDist {int : (subtraj, float)} : dict containing a subtrajectory for
                    a trajectory ID alongwith the distance from the pathlet.
            c1, c3 (float): parameters of the greedy algorithm.
            strajCov ({subtraj : int}) : dict storing coverage (#points) of subtrajs.
            
        Returns:
            The coverage-cost ratio of the pathlet.
    """
    curCov, curCost = 0, c1
    for k, v in trajStrajDist.iteritems():
        straj, dist = v[0], v[1]
        if straj is None:
            continue
        curCov += strajCov[straj]
        curCost += 1.0*c3*dist
    if curCov == 0:
        return 0
    else:
        return (1.0*curCov)/(1.0*curCost)


def optStrajAdvancedHelper(strajDists, strajCov, r, c3):
    """ Help in computing the subtraj of a trajectory with the maximum coverage-cost ratio.
    
        For a guess 'r' on the coverage-cost ratio, computes the subtraj s that maximizes the quantity
        (coverage - c3*r*distance), and if it is > 0, returns the (s,distance), else returns (None, None).
        
        Args:
            strajDists ([(subtraj, float)]): list of subtraj-distance pairs from a single trajectory.
            strajCov ({subtraj : int}): dict storing coverage (#points) of subtrajs.
            r (float): guess on coverage-cost ratio.
            c3 : parameter of greedy algorithm.
            
        Returns:
            (subtraj, float) or (None,None).
    """
    temp, stra, dista = 0, None, None
    for i in xrange(len(strajDists)):
        straj, dist = strajDists[i][0], strajDists[i][1]
        cov = strajCov[straj]
        if temp < cov - c3*r*dist :
            temp = cov - c3*r*dist
            stra, dista = straj, dist
    
    return (stra, dista)


def computeOptStrajsAdvanced(pth, distPairs, pthOptStrajs, strajCov, c1, c3, m, affectedTrajs):
    """ Compute subtrajectories for a pathlet with optimal coverage-cost ratio.
    
        Args:
            pth (pathlet) : pathlet of concern.
            distPairs ({(pathlet, int) : [(subTraj, float)]): dict mapping a
                    pathlet-trajID pair to a list of subtraj-float pairs, where the
                    subtraj belong to the resp. traj, and the float values are computed
                    Frechet distances.
            pthOptStrajs ({pathlet : {int : (subtraj, float)}}): dict storing for a pathlet,
                    the (at most one) optimal subtraj. from  each traj. along with the distance.
                    This is updated by the function.
            strajCov ({subtraj : int}) : dict storing coverage (#points) of subtrajs.
            c1, c3 (float) : parameters of the greedy algorithm.
            m (int) : total no. of points.
            affectedTrajs [int] : list of trajIDs with newly covered points.
            
        Returns:
            Nothin, but updates pthOptStrajs.
    """
    optStrajs = pthOptStrajs[pth]
    ret, temp = {}, {}
    if c3 == 0:
        # Just try to maximize coverage when c3 = 0.
        for trID in affectedTrajs:
            strajDists = distPairs[(pth, trID)]
            # Helper's result does not depend on r when c3 = 0.
            (straj,dist) = optStrajAdvancedHelper(strajDists, strajCov, 1, c3)
            ret[trID] = (straj, dist)
    else:
        # r is a guess on the coverage-cost ratio.
        summation = 0 # Quantity to check when to stop iterating through values of r.
        rmin = 1.0/(c1 + c3)
        r, rmax = rmin, m*rmin
        while r <= rmax:
            for trID in affectedTrajs:
                strajDists = distPairs[(pth, trID)]
                (straj,dist) = optStrajAdvancedHelper(strajDists, strajCov, r, c3)
                temp[trID] = (straj, dist)
                if straj is not None:
                    summation += (strajCov[straj] - c3*r*dist) # summation += cov - c3*r*d.
                    
            if summation < c1*r:
                break
                
            else:
                ret = temp # Replace ret with current subtrajs.
                temp = {}
                summation = 0
            r *= 2
    
    # Update pthOptStrajs.
    pthOptStrajs[pth] = ret


def runGreedy(trajs, distPairs, strajCov, ptStraj, strajPth, trajCov, c1, c2, c3):
    """ Run the greedy algorithm for pathlet cover.
    
    At each step, the algorithm either chooses to leave a point unassigned, or picks
    a pathlet and a set of subtrajectories assigned to the pathlet (at most one from
    each subtrajectory) depending on whichever has the highest coverage-cost ratio.
    The points that are covered by the sets picked up in each greedy step are said to
    be "processed".
    
    Args:
        trajs ({int : traj}): dict mapping ID to traj objects.
        distPairs ({(pathlet, int) : [(subTraj, float)]): dict mapping a
                    pathlet-trajID pair to a list of subtraj-float pairs, where the
                    subtraj belong to the resp. traj, and the float values are computed
                    Frechet distances.
        strajCov ({subtraj : int}) : dict storing coverage (#points) in all subtrajs
                    in distPairs.
        ptStraj ({pt : {subtraj}}) : dict storing for each point, the set of subtrajs
                    in distPairs containing it.
        strajPth ({subtraj: [pathlet]) : dict storing for each subtraj in distPairs,
                the list of pathlets associated with it.
        trajCov ({int : int}) : dict storing the #points in each trajectory.
        c1,c2,c3 (float): parameters of the greedy algorithm.
        
    Returns:
        Pathlet assignments and unassigned points as determined by the greedy algorithm,
        alongwith other relevant info about the pathlets picked.
    """
    
    # Initialize coverage-cost ratios for each pathlet.
    # pthOptCovCost is a dict mapping a pathlet to its optimal coverage-cost ratio.
    # pthOptStrajs is a dict of the form {pathlet : {int : (subtraj, dist)}} mapping
    # a pathlet to the optimal assignment of subtrajectories to it, alongwith the Frechet
    # distances. Both these dicts evolve as the greedy algorithm progresses.
    
    pthOptCovCost, pthOptStrajs = {}, {}
    
    # Build skeleton of pthOptStrajs.
    for key, value in distPairs.iteritems():
        pth, trID, strajDists = key[0], key[1], value
        if pthOptStrajs.has_key(pth) is False:
            pthOptStrajs[pth] = {}
        pthOptStrajs[pth][trID] = (None, None)

    # Compute pthOptStrajs.
    for key, value in distPairs.iteritems():
        pth = key[0]
        affectedTrajs = pthOptStrajs[pth].keys()
        computeOptStrajsAdvanced(pth, distPairs, pthOptStrajs, strajCov, c1, c3, len(ptStraj), affectedTrajs)
    
    # Compute pthOptCovCost from pthOptStrajs.
    for key, value in pthOptStrajs.iteritems():
        pth = key
        pthOptCovCost[pth] = computeCovCostRatio(value, c1, c3, strajCov)
        if pthOptCovCost[pth] == 0:
            print "Error"
      
    # Initialize a max priority queue of pathlets ordered by coverage cost ratios.
    queue1 = heapdict()
    for pth, ccratio in pthOptCovCost.iteritems():
        queue1[pth] = -1.0*ccratio # Need to negate, since heapdict is a min-heap.
    
    # Initialize a priority queue of trajs, with coverage to cost ratio of a singleton set, i.e., |T|/c2.
    queue2 = heapdict()
    for trID, cov in trajCov.iteritems():
        queue2[trID] = -(1.0*cov)/(1.0*c2)
    
    # Initial pathlet assignments, it is of the form {pathlet : [subtraj]} and stores the 
    # list of subtraj assigned to the pathlets 
    pthAssignments = {}
    pthStats = [] # this is of the form [(pathlet,ccRatio,iter,fracThickness)], where iter was
                  # the iteration of the greedy algorithm when pathlet was picked,
                  # and fracThickness is the proportion of points of a traj assigned to pth,
                  # summed over all trajs. Note that the same pathlet can be picked multiple
                  # times in different iterations.

    # Unassigned points
    unassignedPts = []

    count = 0
    numUnprocessedPts = [len(ptStraj)]
    while numUnprocessedPts[0] > 0:
        print "num of points is %d" %numUnprocessedPts[0]
        x1, x2 = queue1.peekitem(), queue2.peekitem()

        # Set of the form {(pth, trID)} whose coverage changes after new points are
        # processed in the current iteration.
        affectedPths = set()
        
        # If the most beneficial pathlet is more beneficial than leaving a point unassigned.
        if x1[1] <= x2[1]:
            x = queue1.popitem()
            pth = x[0]
            fracThickness = 0
            for trID, pair in pthOptStrajs[pth].iteritems():
                straj = pair[0]
                if straj is None : # possible when all strajs of trID in pth's kitty have zero coverage
                    continue
                if pthAssignments.has_key(pth):
                    pthAssignments[pth].append(straj)
                else:
                    pthAssignments[pth] = [straj]
                    
                fracThickness += 1.0*strajCov[straj]/len(trajs[straj.trajID].pts)
                
                # Process the straj, and also return affected pathlets, i.e., pathlets whose optimal
                # coverage-cost ratio changes since some points that could have been assigned to it
                # have now been processed.
                affectedPths = affectedPths.union(processSubtraj(straj, strajCov, trajs, trajCov, ptStraj, strajPth, distPairs, numUnprocessedPts, queue2))
                
            pthStats.append((pth,pthOptCovCost[pth],count,fracThickness))

        # If more beneficial to leave a point unassigned than picking a pathlet.
        else:            
            x = queue2.popitem()
            trID = x[0]
            # Process the unassigned point, and also return pathlets whose optimal coverage-cost ratio
            # changes.
            affectedPths = affectedPths.union(processTraj(trID, ptStraj, strajCov, strajPth, trajs, trajCov, distPairs, numUnprocessedPts, queue2, unassignedPts))
        
        # Update coverage-cost ratio of affected pathlets.
        affectedPathlets = {path for (path, traID) in affectedPths}
        affectedTrajs = {traID for (path, traID) in affectedPths}

        for path in affectedPathlets:
            computeOptStrajsAdvanced(path, distPairs, pthOptStrajs, strajCov, c1, c3, len(ptStraj), affectedTrajs)
            pthOptCovCost[path] = computeCovCostRatio(pthOptStrajs[path], c1, c3, strajCov)
            queue1[path] = -1.0*pthOptCovCost[path]
            
        count += 1
        
    return (pthAssignments, pthStats, unassignedPts)

# test above function
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print "Wrong command. See README for details."
        
    else:
        print "Loading trajectories ..."
        trajs = pickle.load(open(sys.argv[1] + ".CleanedUp.p", "rb"))
        c1, c2, c3 = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
        #for trajID in trajs:
            #print "TrajID %d : %d points " %(trajID, len(trajs[trajID].pts))
        ret = pickle.load(open(sys.argv[1] + ".pre.sc.p", "rb"))
        print "Loading prerequisite data structures ..."
        distPairs = pickle.load(open(sys.argv[1] + ".distPairs2.sc.p", "rb"))
        (strajCov, ptStraj, strajPth, trajCov) = ret
        #print strajPth
        #validate(strajCov, ptStraj,  trajs)
        #testProcessSubtraj(strajCov, trajs, None, ptStraj, strajPth, distPairs)
        print "Running greedy algorithm ..."
        retVal = runGreedy(trajs, distPairs, strajCov, ptStraj, strajPth, trajCov, c1, c2, c3)
        print "Finished."
        pickle.dump(retVal, open(sys.argv[1] +".out.sc.p", "wb"))

