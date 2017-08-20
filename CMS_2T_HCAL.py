from heppy.papas.detectors.detector import Detector, DetectorElement
import heppy.papas.detectors.material as material
from heppy.papas.detectors.geometry import VolumeCylinder
import math
import heppy.statistics.rrandom as random

class ECAL(DetectorElement):

    def __init__(self):
        volume = VolumeCylinder('ecal', 1.55, 2.1, 1.30, 2. )
        mat = material.Material('ECAL', 8.9e-3, 0.275)
        self.eta_crack = 1.479
        self.emin = {'barrel':0.3, 'endcap':1.}
        self.eres = {'barrel':[4.22163e-02, 1.55903e-01, 7.14166e-03], 'endcap':[-2.08048e-01, 3.25097e-01, 7.34244e-03]}
        self.eresp = {'barrel':[1.00071, -9.04973, -2.48554], 'endcap':[9.95665e-01, -3.31774, -2.11123]}
        super(ECAL, self).__init__('ecal', volume,  mat)

    def energy_resolution(self, energy, eta=0.):
        part = 'barrel'
        if abs(eta)>1.479 and abs(eta)<3.0:
            part = 'endcap'
        stoch = self.eres[part][0] / math.sqrt(energy)
        noise = self.eres[part][1] / energy
        constant = self.eres[part][2]
        return math.sqrt( stoch**2 + noise**2 + constant**2) 

    def energy_response(self, energy, eta=0):
        part = 'barrel'
        if abs(eta)>self.eta_crack:
            part = 'endcap'
        return self.eresp[part][0]/(1+math.exp((energy-self.eresp[part][1])/self.eresp[part][2])) #using fermi-dirac function : [0]/(1 + exp( (energy-[1]) /[2] ))

    def cluster_size(self, ptc):
        pdgid = abs(ptc.pdgid())
        if pdgid==22 or pdgid==11:
            return 0.04
        else:
            return 0.07

    def acceptance(self, cluster):
        energy = cluster.energy
        eta = abs(cluster.position.Eta())
        if eta < self.eta_crack:
            return energy>self.emin['barrel']
        elif eta < 2.93:
            return energy>self.emin['endcap'] and cluster.pt>0.2
        else:
            return False

    def space_resolution(self, ptc):
        pass
    
class HCAL(DetectorElement):

    def __init__(self):
        volume = VolumeCylinder('hcal', 4.3, 5.3, 1.9, 2.85 )
        # not sure about X0 and lambda_i, but these don't matter anyway
        mat = material.Material('HCAL', 0.018, 0.17)
        # resolution from CLIC CDR Fig. 6.11, 1st hypothesis
        self.eres = [0.60, 0., 0.025]        
        super(HCAL, self).__init__('ecal', volume, mat)

    def energy_resolution(self, energy, eta=0.):
        stoch = self.eres[0] / math.sqrt(energy)
        noise = self.eres[1] / energy
        constant = self.eres[2]
        return math.sqrt( stoch**2 + noise**2 + constant**2)

    def energy_response(self, energy, eta=0):
        return 1.0
    
    def cluster_size(self, ptc):
        '''returns cluster size in the HCAL
        
        25 cm for CLIC, c.f. CLIC CDR Fig. 6.12
        Value taken to get a 80% chance of separating two showers.
        '''
        return 0.25

    def acceptance(self, cluster):
        energy = cluster.energy
        eta = abs(cluster.position.Eta())
        if eta < 2.76:  #TODO: check this value
            return energy>1.
        else:
            return False
    
    def space_resolution(self, ptc):
        pass


    
class Tracker(DetectorElement):
    #TODO acceptance and resolution 
    #depend on the particle type
    
    def __init__(self):
        volume = VolumeCylinder('tracker', 1.29, 1.99)
        mat = material.void
        super(Tracker, self).__init__('tracker', volume,  mat)

    def acceptance(self, track):
        # return False
        pt = track.p3() .Pt()
        eta = abs(track.p3() .Eta())
        if eta < 1.35 and pt>0.5:
            return random.uniform(0,1)<0.95
        elif eta < 2.5 and pt>0.5:
            return random.uniform(0,1)<0.9 
        else:
            return False

    def resolution(self, track):
        # TODO: depends on the field
        pt = track.p3() .Pt()
        return 1.1e-2

    

class Field(DetectorElement):

    def __init__(self, magnitude):
        self.magnitude = magnitude
        volume = VolumeCylinder('field', 2.9, 3.6)
        mat = material.void
        super(Field, self).__init__('tracker', volume,  mat)

class BeamPipe(DetectorElement):
    '''Beam pipe is not used in the simulation at the moment, so no need to define it.'''

    def __init__(self):
        #Material Seamless AISI 316 LN, External diameter 53 mm, Wall thickness 1.5 mm (hors cms) X0 1.72 cm
        #in CMS, radius 25 mm (?), tchikness 8mm, X0 35.28 cm : berylluim
        factor = 1.0
        volume = VolumeCylinder('beampipe', 2.5e-2*factor+0.8e-3, 1.98, 2.5e-2*factor, 1.9785 )
        mat = material.Material('BeamPipe', 35.28e-2, 0)
        super(BeamPipe, self).__init__('beampipe', volume, mat)

        
class CMS(Detector):
        
    def electron_acceptance(self, track):
        return track.p3() .Mag() > 5 and abs(track.p3() .Eta()) < 2.5

    def electron_resolution(self, ptc):
        return 0.1 / math.sqrt(ptc.e())
            
    def muon_acceptance(self, track):
        return track.p3() .Pt() > 5 and abs(track.p3() .Eta()) < 2.5
            
    def muon_resolution(self, ptc):
        return 0.02 
    
    def __init__(self):
        super(CMS, self).__init__()
        self.elements['tracker'] = Tracker()
        self.elements['ecal'] = ECAL()
        self.elements['hcal'] = HCAL()
        self.elements['field'] = Field(2.)
        self.elements['beampipe'] = BeamPipe()

cms = CMS()
