class pt(object):
    """ A class representing a point in a trajectory.
    
    Attributes:
        lat (float): latitude.
        lon (float): longitude.
        trajID (int): ID of the point's trajectory.
        t (float): timestamp associated with the point.
    """
    
    def __init__(self):
        """ Default constructor with dummy initialization values. """
        self.lat = 0.0
        self.lon = 0.0
        self.trajID = -1
        self.t = -1.0
        
    def __hash__(self):
        """ Computes a hash so that pt objects can be used as dict keys. """
        return hash((self.lat, self.lon, self.trajID, self.t))
        
    def __eq__(self, other):
        """ Define == operator for pt objects. """
        return (self.lat==other.lat and self.lon==other.lon and self.trajID==other.trajID and self.t==other.t)
        
    def __str__(self):
        """ Return string to be output while printing a pt object. """
        return "Point TrajID %d ; lat-long (%f,%f); time %f" % (self.trajID, self.lat, self.lon, self.t)
    


class traj(object):
    """ A class representing a trajectory. 
    
    Attributes:
        pts: list of points in the trajectory.
    """
    
    def __init__(self):
        """ Initialize trajectory with empty list of points. """
        self.pts = []

    def addPt(self, lat, lon, trajID, t):
        """ Add a pt to the trajectory.
        
        Args:
            lat (float): latitude of point.
            lon (float): longitude of point.
            trajID (int): trajID of the point (all points of a traj will have same ID).
        """
        p = pt()
        p.lat = lat
        p.lon = lon
        p.trajID = trajID
        p.t = t
        self.pts.append(p)

    def sortPts(self):
        """ Sort points of trajectory in ascending order of timestamp. """
        self.pts = sorted(self.pts, key = lambda x: x.t)
        

class pathlet(object):
    """ A class representing a pathlet (essentially a subtrajectory).
    
        Similar to traj, however avoids storing the points explicitly by referring to the
        trajectory to which the points belong, along with the starting and ending indices
        of the points in the list of points of the trajectory, sorted by timestamp.
        
        Attributes:
            trajID (int): ID of the trajectory to which points of the pathlet belong.
            bounds (int,int): start and end indices of points in the timestamp-sorted list of
                              trajectory points.
    """
    
    def __init__(self, trajID, bounds):
        """ Initialize pathlet with points from a trajectory.
        
            Args:
                trajID (int): ID of trajectory.
                bounds (int,int): start and end indices.
        """
        self.trajID = trajID
        self.bounds = bounds
        
    def __str__(self):
        """ Return string to be output while printing a pathlet object. """
        return "Pathlet TrajID %d ; bounds (%d, %d)" % (self.trajID, self.bounds[0], self.bounds[1])
        
    def __eq__(self, other):
        """ Define == operator for pathlet objects. """
        return (self.trajID==other.trajID and self.bounds[0]==other.bounds[0] and self.bounds[1]==other.bounds[1])
        
    def __hash__(self):
        """ Define a hash function so that pathlets can be used as keys in a dict. """
        return hash((self.trajID, self.bounds[0], self.bounds[1]))
        


class subTraj(object):
    """ A class representing a subtrajectory.
    
    Exactly identical to a pathlet class. Defined as a class of its own for conceptual reasons.
    """
    
    def __init__(self, trajID, bounds):
        self.trajID = trajID
        self.bounds = bounds
        
    def __str__(self):
        return "Subtraj TrajID %d ; bounds (%d, %d)" % (self.trajID, self.bounds[0], self.bounds[1])
        
    def __eq__(self, other):
        return (self.trajID==other.trajID and self.bounds[0]==other.bounds[0] and self.bounds[1]==other.bounds[1])
        
    def __hash__(self):
        return hash((self.trajID, self.bounds[0], self.bounds[1]))

