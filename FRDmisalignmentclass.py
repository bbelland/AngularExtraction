class FRDmisalignmentextract(FRDlist,thetalist,philist):
    
    
    def __init__(self,FRDlist,thetalist,philist):
        """Generates an Angular Misalignment Extraction Object"""
        
        #Generate a function to accept inputs other than self
        self.FRDlist = FRDlist
        self.FRDerrlist = FRDerrlist
        self.thetalist = thetalist
        self.philist = philist
        
        #Run checks to verify the format of FRDlist, thetalist, philist
        
        self.FRDcheck()
        self.thetacheck()
        self.phicheck()
        
        #Then test for angular misalignment.
        
    def FRDcheck(self):
        
        if not (isinstance(self.FRDlist,list) or isinstance(self.FRDlist,np.ndarray)):
            #Error
            raise Exception('FRDlist should be a np.ndarray or a list. FRDlist\'s current type is {}'.format(type(self.FRDlist))
            
        elif isinstance(self.FRDlist,list):
            self.FRDlist = np.array(self.FRDlist)
    
    def thetacheck(self):
            
        if not (isinstance(self.thetalist,list) or isinstance(self.thetalist,np.ndarray) or isinstance(self.thetalist,np.ndarray)):
            #Error
            raise Exception('thetalist should be a np.ndarray, a list, or a float. thetalist\'s current type is {}'.format(type(self.thetalist))
            
        elif isinstance(self.thetalist,np.float):
            #Need to set thetalist as a constant list of length FRDlist if only one value input
            self.thetalist = self.thetalist*np.array(np.ones(np.shape(self.FRDlist)))
            
        elif isinstance(self.thetalist,list):
            self.thetalist = np.array(self.thetalist)
            
    def phicheck(self):
           
        if not (isinstance(self.philist,list) or isinstance(self.philist,np.ndarray) or isinstance(self.philist,np.ndarray)):
            #Error
            raise Exception('philist should be a np.ndarray, a list, or a float. philist\'s current type is {}'.format(type(self.philist))
            
        elif isinstance(self.philist,np.float):
            #Need to set philist as a constant list of length philist if only one value input
            self.philist = self.philist*np.array(np.ones(np.shape(self.philist)))
            
        elif isinstance(self.philist,list):
            self.philist = np.array(self.philist)   
                            

    #What I get: FRD as a function of theta and phi position
    #What I need to extract: Angular misalignment
    #In *theory* this could be doable. In a sweep of theta (say), one goes from a maximum angular misalignment to a minimum angular misalignment. That's the first obvious test. You know that you're near alignment if a scan over all theta positions yields the same angle (FRD). 
                            
    #If I assume that I'm given (or find) the theta, phi position with maximum angular misalignment, then I know that the angular misalignments there are generally adding constructively. Then I can use the shifts by pi in theta and phi to deduce more information. (Of course, directly extracting angular misalignment from the images would be ideal, but this is a good first step).