from detector import Detector, DetectorElement
import material as material
from geometry import VolumeCylinder
import math
import heppy.statistics.rrandom as random

class ECAL(DetectorElement):

    def __init__(self):
        depth = 0.25
        inner_radius = 2.15
        inner_z = 2.6
        nX0 = 23  #CLIC CDR, page 70, value for CLIC_ILD
        nLambdaI = 1  # ibid
        outer_radius = inner_radius + depth
        outer_z = inner_z + depth
        X0 = depth / nX0
        lambdaI = depth / nLambdaI
        volume = VolumeCylinder('ecal', outer_radius, outer_z, inner_radius, inner_z)
        mat = material.Material('ECAL', X0, lambdaI)
        # todo: recompute
        self.eta_junction = volume.inner.eta_junction()
        # as for ILD (thresholds chosen by Mogens)
        self.emin = {'barrel':0.5, 'endcap':0.5}
        # CLIC CDR p.123. adding a noise term of 1%
        self.eres = {'barrel':[0.167, 0.010, 0.011]}
        super(ECAL, self).__init__('ecal', volume,  mat)

    def energy_resolution(self, energy, eta=0.):
        part = 'barrel'
        stoch = self.eres[part][0] / math.sqrt(energy)
        noise = self.eres[part][1] / energy
        constant = self.eres[part][2]
        return math.sqrt( stoch**2 + noise**2 + constant**2) 

    def energy_response(self, energy, eta=0):
        return 1
    
    def cluster_size(self, ptc):
        '''just guessing numbers (from Mogens, as in ILD).'''
        pdgid = abs(ptc.pdgid())
        if pdgid==22 or pdgid==11:
            return 0.015
        else:
            return 0.045

    def acceptance(self, cluster):
        energy = cluster.energy
        eta = abs(cluster.position.Eta())
        if eta < self.eta_junction:
            return energy>self.emin['barrel']
        elif eta < 2.76:  #TODO check this value
            return energy>self.emin['endcap']
        else:
            return False

    def space_resolution(self, ptc):
        pass
    
class HCAL(DetectorElement):

    def __init__(self):
        volume = VolumeCylinder('hcal', 2.9, 3.6, 1.9, 2.6 )
        mat = material.Material('HCAL', None, 0.17)
        self.eta_crack = 1.3
        self.eres = {'barrel':[0.8062, 2.753, 0.1501], 'endcap':[6.803e-06, 6.676, 0.1716]}
        self.eresp = {'barrel':[1.036, 4.452, -2.458], 'endcap':[1.071, 9.471, -2.823]}
        super(HCAL, self).__init__('ecal', volume, mat)

    def energy_resolution(self, energy, eta=0.):
        part = 'barrel'
        if abs(eta)>self.eta_crack:
            part = 'endcap'
        # stoch = self.eres[part][0] / math.sqrt(energy)
        # noise = self.eres[part][1] / energy
        # constant = self.eres[part][2]
        stoch = 1.1 / math.sqrt(energy)
        noise = 0
        constant = 0.09
        return math.sqrt( stoch**2 + noise**2 + constant**2)

    def energy_response(self, energy, eta=0):
        part = 'barrel'
        if abs(eta)>self.eta_crack:
            part = 'endcap'
        return self.eresp[part][0]/(1+math.exp((energy-self.eresp[part][1])/self.eresp[part][2])) #using fermi-dirac function : [0]/(1 + exp( (energy-[1]) /[2] ))

    def cluster_size(self, ptc):
        return 0.2

    def acceptance(self, cluster):
        energy = cluster.energy
        eta = abs(cluster.position.Eta())
        if eta < self.eta_crack :
            if energy>1.:
                return random.uniform(0,1)<(1/(1+math.exp((energy-1.93816)/(-1.75330))))
            else:
                return False
        elif eta < 3. : 
            if energy>1.1:
                if energy<10.:
                    return random.uniform(0,1)<(1.05634-1.66943e-01*energy+1.05997e-02*(energy**2))
                else:
                    return random.uniform(0,1)<(8.09522e-01/(1+math.exp((energy-9.90855)/-5.30366)))
            else:
                return False
        elif eta < 5.:
            return energy>7.
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
    
    
    def jet_energy_correction(self, jet):
        '''The factor roughly corresponds to the raw PF jet response in CMS,
        which is around 90%. The factor was checked in the reconstruction
        of Z->jj in papas.
        '''
        return 1.1
    
    def __init__(self):
        super(CMS, self).__init__()
        self.elements['tracker'] = Tracker()
        self.elements['ecal'] = ECAL()
        self.elements['hcal'] = HCAL()
        self.elements['field'] = Field(2.0)
        self.elements['beampipe'] = BeamPipe()

cms = CMS()
