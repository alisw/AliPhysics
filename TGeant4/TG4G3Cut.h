// $Id$
// Category: global

#ifndef TG4_G3_CUT_H
#define TG4_G3_CUT_H

enum TG4G3Cut
{
  kCUTGAM, // gammas 
               // G4 particles: "gamma"         
               // G3 default value: 0.001 GeV
  kCUTELE, // electrons        
               // G4 particles: "e-"         
               // ?? positrons
               // G3 default value: 0.001 GeV
  kCUTNEU, // neutral hadrons  
               // G4 particles: of type "baryon", "meson", "nucleus"   
	       //               with zero charge     
               // G3 default value: 0.01 GeV
  kCUTHAD, // charged hadrons  
               // G4 particles: of type "baryon", "meson", "nucleus"        
	       //               with non-zero charge     
               // G3 default value: 0.01 GeV
  kCUTMUO, // muons            
               // G4 particles: "mu+", "mu-"         
               // G3 default value: 0.01 GeV
  kBCUTE,  // electron bremsstrahlung         
               // G4 particles: "gamma"         
               // G3 default value: CUTGAM
  kBCUTM,  // muon and hadron bremsstrahlung  
               // G4 particles: "gamma"         
               // G3 default value: CUTGAM
  kDCUTE,  // delta-rays by electrons 
               // G4 particles: "e-"         
               // G3 default value: 10**4
  kDCUTM,  // delta-rays by muons        
               // G4 particles: "e-"         
               // G3 default value: 10**4
  // the following cut is not yet implemented
  kPPCUTM, // direct pair production by muons 
               // G4 particles: ??         
               // G3 default value: 0.01 GeV
  kNoG3Cuts    
};

#endif //TG4_G3_CUT_H
