import numpy as np

from scipy.interpolate import interp1D,interp2D
from scipy.optimize import root

def sigFRDtosigTOT(Asquared):
    return 1+ (1+Asquared)/(2+1.5*Asquared)*Asquared

def angularmisalignment_position(baseFRD,sigTOT):
              
    squaredsigTOTtobaseFRD = sigTOT**2./baseFRD**2  
              
    def functionmin(Asquared,squaredfrac):
        return sigFRDtosigTOT(Asquared)-squaredfrac
              
    Asquared = root(functionmin, x0=[12.*(squaredsigTOTtobaseFRD-1.)/7],args=(squaredsigTOTtobaseFRD))#The smallest angular misalignment in the extreme, highest true sigFRD case. This should return 0 or else something weird is going on.                               
        
    misalignment = np.sqrt(Asquared)*baseFRD   
              
    return misalignment

class FRDmisalignmentextract(FRDarray,thetaarray,phiarray):
    
    
    def __init__(self,FRDarray,thetaarray,phiarray):
        """Generates an Angular Misalignment Extraction Object"""
        
        #Generate a function to accept inputs other than self
        self.FRDarray = FRDarray
        # self.FRDerrarray = FRDerrarray
        self.thetaarray = thetaarray
        self.phiarray = phiarray
        
        #Run checks to verify the format of FRDarray, thetaarray, phiarray
        
        self.FRDcheck()
        self.thetacheck()
        self.phicheck()
        
        if not self.phiscan and not self.thetascan:
            raise Exception('At least one of theta or phi should cover a range of values')
        
        #Then test for angular misalignment.
        
    def FRDcheck(self):
        
        if not (isinstance(self.FRDarray,list) or isinstance(self.FRDarray,np.ndarray)):
            #Error
            raise Exception('FRDarray should be a np.ndarray or a list. FRDarray\'s current type is {}'.format(type(self.FRDarray))
            
        elif isinstance(self.FRDarray,list):
            self.FRDarray = np.array(self.FRDarray)
    
    def thetacheck(self):
                            
        #Still need to check that length is valid relative to FRD array.
            
        if not (isinstance(self.thetaarray,list) or isinstance(self.thetaarray,np.ndarray) or isinstance(self.thetaarray,np.ndarray)):
            #Error
            raise Exception('thetaarray should be a np.ndarray, a list, or a float. thetaarray\'s current type is {}'.format(type(self.thetaarray))
            
        elif isinstance(self.thetaarray,np.float):
            #Need to set thetaarray as a constant list of length FRDarray if only one value input
            self.thetaarray = self.thetaarray*np.array(np.ones(np.shape(self.FRDarray)))
            self.thetascan = False #If we know theta is held at one value, we can make interpolation easier

            
        elif isinstance(self.thetaarray,list):
            self.thetaarray = np.array(self.thetaarray)
            if np.max(np.mod(self.thetaarray,2*np.pi)) - np.min(np.mod(self.thetaarray,2*np.pi)) > 0:
                self.thetascan = True #There is a range of angles in theta
            else:
                self.thetascan = False #Only one effective angle in theta

            #Check length. Could force an easier way to input or correct the harder way to input
                            
        else: #If thetaarray is a np.ndarray
            if np.max(np.mod(self.thetaarray,2*np.pi)) - np.min(np.mod(self.thetaarray,2*np.pi)) > 0:
                self.thetascan = True #There is a range of angles in theta
            else:
                self.thetascan = False #Only one effective angle in theta
            
    def phicheck(self):
                            
        #Still need to check that length is valid relative to FRD array.
           
        if not (isinstance(self.phiarray,list) or isinstance(self.phiarray,np.ndarray) or isinstance(self.phiarray,np.ndarray)):
            #Error
            raise Exception('phiarray should be a np.ndarray, a list, or a float. phiarray\'s current type is {}'.format(type(self.phiarray))
            
        elif isinstance(self.phiarray,np.float):
            #Need to set phiarray as a constant list of length phiarray if only one value input
            self.phiarray = self.phiarray*np.array(np.ones(np.shape(self.phiarray)))
            self.phiscan = False #If we know theta is held at one value, we can make interpolation easier
            
        elif isinstance(self.phiarray,list):
            self.phiarray = np.array(self.phiarray)
            if np.max(np.mod(self.phiarray,2*np.pi)) - np.min(np.mod(self.phiarray,2*np.pi)) > 0:
                self.phiscan = True #There is a range of angles in theta
            else:
                self.phiscan = False #Only one effective angle in theta
        
        else: #If phiarray is a np.ndarray
            if np.max(np.mod(self.phiarray,2*np.pi)) - np.min(np.mod(self.phiarray,2*np.pi)) > 0:
                self.phiscan = True #There is a range of angles in theta
            else:
                self.phiscan = False #Only one effective angle in theta
                            

    #What I get: FRD as a function of theta and phi position
    #What I need to extract: Angular misalignment
    #In *theory* this could be doable. In a sweep of theta (say), one goes from a maximum angular misalignment to a minimum angular misalignment. That's the first obvious test. You know that you're near alignment if a scan over all theta positions yields the same angle (FRD). 
                            
    #If I assume that I'm given (or find) the theta, phi position with maximum angular misalignment, then I know that the angular misalignments there are generally adding constructively. Then I can use the shifts by pi in theta and phi to deduce more information. (Of course, directly extracting angular misalignment from the images would be ideal, but this is a good first step).
                            
    #It's not strictly impossible for there to always be some angular misalignment at any Cobra position (most obviously where focal plane-Cobra body angle is nonzero but with zero misalignment between Cobra body and Cobra arm, Cobra arm and fiber). So the code has to find the maximum angle, find the minimum angle, find the difference between the angles as much as possible to create a model, and then output the model.
                            
    def maxanglecheck(self):
        phiarray = self.phiarray
        thetaarray = self.thetaarray
        FRDarray = self.FRDarray
                            
        firstguessindex = np.argwhere(FRDarray == np.max(FRDarray)) #Not great because it gets confused if there's any noise in FRD extraction (there will be!) Needs to more smartly find the maximum.
                            
        phi_maxFRD = phiarray[firstguessindex]
        theta_maxFRD = thetaarray = thetaarray[firstguessindex]
        
        if self.phiscan and self.thetascan: #A.K.A there's a 2D 
            FRDfunction = interp2D(thetaarray,phiarray,FRDarray) #Could use RectBivariateSpline instead. Might be worth investigating
                            
        elif self.phiscan and not self.thetascan:
            FRDfunction = interp1D(phiarray,FRDarray)
                            
        else: #if self.thetascan and not self.phiscan:
            FRDfunction = interp1D(thetaarray,FRDarray)
        
        #phimax, thetamax corresponds to angle of the maximum misalignment. Presumably, this corresponds to where all of the angular misalignments constructively add together. At the very least, some interpolation would be better.
                            
        #Next we can learn about one angle by taking thetamax and taking the closest value to pi away
                            
        theta_anglefrommaxFRD = np.mod(thetaarray-theta_maxFRD, 2*np.pi) #This is the most important value to ascertain, since it should always be available and 
        theta_piawayfrommaxFRD_condition = np.abs(theta_anglefrommaxFRD-np.pi) #These are angular distances from pi away
        theta_piawayfrommaxFRD_anglediff = np.min(theta_piawayfrommaxFRD_condition)
                            
        print('Angular difference between the maximum FRD angle and pi away in data points is 
              
        #Again, obviously interpolation would be better.
                            
    def AngularMisalignmentRange(self,minfiberFRD):
        minSigTot = np.min(self.FRDarray)
        maxSigTot = np.max(self.FRDarray)
        
    def firstguess_AngMis(self):
              
        lowestFRD = np.min(self.FRDarray)      
              
        minFRD = 0.06 #Could be an input based on the fiber characteristic(s) in the system.
        maximum_minimumFRD = lowestFRD #Obviously there cannot be any lower FRD than the minimum sigTOT of the system, neglecting errors
              
        minmisalignment_max = angularmisalignment_position(minFRD,lowestFRD)
        minmisalignment_min = angularmisalignment_position(lowestFRD,lowestFRD) #Should be zero.
                    
        highestFRD = np.max(self.FRDarray)      
            
        maxmisalignment_max = angularmisalignment_position(minFRD,highestFRD)
        maxmisalignment_min = angularmisalignment_position(lowestFRD,highestFRD)
              
        return [[minmisalignment_min, maxmisalignment_min],[minmisalignment_max, maxmisalignment_max]]
                      
    def guessbestbaseFRD(self):
        
        if self.thetascan == 1:
            thetaangle = self.guessbestthetaangle() #Try to probe the 2pi radians of the theta rotation to feed into the angular misalignment model
        if self.phiscan == 1:
            phiangle = self.guessbestphiangle() #Try to probe the pi radians of the phi rotation to feed into my angular misalignment model
              
        Cobraangle = self.guessbestCobraangle() #Try to determine the constant angular offset
              
        bestmaxangle = Cobraangle+phiangle+Cobraangle #Not really, since the phi angle may not align with the Cobra+theta angle, but.
        
        functionFRD(sigFRDsq,sigTOT,a):
              return sigTOT**2 - sigFRDsq*(1+(1+a**2/sigFRDsq)/(2+1.5*a**2/sigFRDsq))
        
        bestbaseFRDsquared = root(functionFRD, x0=[sigTOT],args=(np.max(self.FRDarray),bestmaxangle))
              
        return np.sqrt(bestbaseFRDsquared)
              
    def guessbestthetaangle(self):
        pass
              
    def guessbestphiangle(self):
        pass
              
    def guessbestCobraangle(self):
        pass
