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
        return 1.
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
        
    def electron_acceptance(self, ptc):
        """Delphes parametrization
        https://github.com/delphes/delphes/blob/master/cards/delphes_card_CMS.tcl
        96d6bcf 
        """
        rnd = random.uniform(0, 1)
        if ptc.pt() < 10.:
            return False
        else:
            eta = abs(ptc.eta())
            if eta < 1.5:
                return rnd < 0.95
            elif eta < 2.5:
                return rnd < 0.85
            else:
                return False

    def electron_resolution(self, ptc):
        # return 0.1 / math.sqrt(ptc.e())
        return 0.03
            
    def muon_acceptance(self, ptc):
        """Delphes parametrization
        https://github.com/delphes/delphes/blob/master/cards/delphes_card_CMS.tcl
        96d6bcf 
        """        
        rnd = random.uniform(0, 1)
        eta = abs(ptc.eta())        
        if ptc.pt() < 10.:
            return False
        elif eta < 2.4:
            return rnd < 0.95
        else:
            return False
            
    def muon_resolution(self, ptc):
        """Delphes parametrization
        
          # resolution formula for muons
          set ResolutionFormula {
                         (abs(eta) <= 0.5) * (pt > 0.1) * sqrt(0.01^2 + pt^2*1.0e-4^2) +
                         (abs(eta) > 0.5 && abs(eta) <= 1.5) * (pt > 0.1) * sqrt(0.015^2 + pt^2*1.5e-4^2) +
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1) * sqrt(0.025^2 + pt^2*3.5e-4^2)}
        """
        rnd = random.uniform(0, 1)
        eta = abs(ptc.eta())
        cstt = None
        vart = None
        if eta < 0.5:
            cstt, vart = 0.01, 1e-4
        elif eta < 1.5:
            cstt, vart = 0.015, 1.5e-4
        else:
            cstt, vart = 0.025, 3.5e-4
        res = math.sqrt(cstt**2 + vart**2)
        return res
    
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
        self.elements['field'] = Field(2.)
        self.elements['beampipe'] = BeamPipe()

cms = CMS()
