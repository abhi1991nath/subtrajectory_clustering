# Python script to run all modules in one go, without generating intermediate files.
import cPickle as pickle

from base import *
from distance import *
from greedy import *
from ioUtils import *

if __name__ == "__main__":
    #if len(sys.argv) != 7:
    #    print "Wrong command. Correct command is python allInOne.py fileName rmin rmax c1 c2 c3."
        
        print "Loading trajectories ..."
        trajs = readTrajsFromTxtFile("data/small.txt")
        
        rmin, rmax = 0.5, 2
        
        print "Computing Frechet distances ..."
        distPairs1 = process(trajs, rmin, rmax)
        
        #distPairs1 is of form {(pth, straj):dist}, change it to distPairs2 of the form {(pth, trajID):[(straj, dist)]}
        distPairs2 = {}
        for k,v in distPairs1.iteritems():
            pth, trID, dist, straj = k[0], k[1].trajID, v, k[1]
            if distPairs2.has_key((pth, trID)):
                distPairs2[(pth, trID)].append((straj,dist))
            else:
                distPairs2[(pth, trID)] = [(straj,dist)]
                
        print "Computing prerequisite data structures ..."
        (strajCov, ptStraj, strajPth, trajCov) = preprocessGreedy(trajs, distPairs2)
        
        c1,c2,c3 = 1,1,1
        
        print "Running greedy algorithm ..."
        retVal = runGreedy(trajs, distPairs2, strajCov, ptStraj, strajPth, trajCov, c1, c2, c3)
