// $Id$
// Category: physics

#ifndef TG3_FLAG_H
#define TG3_FLAG_H

enum TG3Flag
{
  kPAIR, // pair production       
             // G3 default value: 1
             // G4 processes: G4GammaConversion
             //               (?? G4MuPairProduction/G4IMuPairProduction)
  kCOMP, // Compton scattering    
             // G3 default value: 1
             // G4 processes: G4ComptonScattering 
  kPHOT, // photo electric effect 
             // G3 default value: 1
             // G4 processes: G4PhotoElectricEffect
  kPFIS, // photofission          
             // G3 default value: 0
             // G4 process: ??
  kDRAY, // delta-ray              
             // G3 default value: 2
	     // CHECK:
             // G4 processes: G4eIonisation/G4IeIonization,
	     //               G4eIonisationPlus (??)
             //               G4MuIonisation/G4IMuIonization, 
	     //               G4hIonisation/G4IhIonisation
	     // !! G4 treats delta rays in different way
  kANNI, // annihilation          
             // G3 default value: 1
             // G4 processes: G4eplusAnnihilation/G4IeplusAnnihilation
	     // only for e+ 
  kBREM, // bremsstrahlung        
             // G3 default value: 1
             // G4 processes: G4eBremsstrahlung/G4IeBremsstrahlung,
	     //               G4eBremsstrahlungPlus (??),  
	     //               G4MuBremsstrahlung/G4IMuBremsstrahlung
	     // only for e-/e+; mu+/mu- 
  kHADR, // hadronic process      
             // G3 default value: 1
             // ??
  kMUNU, // muon nuclear interaction 
             // G3 default value: 0
	     // G4 processes: G4MuNuclearInteraction
  kDCAY, // decay                 
             // G3 default value: 1
	     // G4 process: G4Decay
  kLOSS, // energy loss           
             // G3 default value: 2
             // G4 processes: G4eIonisation/G4IeIonization,
	     //               G4eIonisationPlus (??)
             //               G4MuIonisation/G4IMuIonization, 
	     //               G4hIonisation/G4IhIonisation
  kMULS, // multiple scattering   
             // G3 default value: 1
	     // G4 process: G4MultipleScattering/G4IMultipleScattering
	     // all charged particles
/* to be added	 
  kCKOV  // Cerenkov photon generation
           // G3 default value: 0
	   // G4 process: G4Cerenkov
	   //             + light photon absorption processes (??which)
	   // all charged particles  
  kRAYL, // Rayleigh scattering
           // G3 default value: 0	     
	   // G4 process: ?? G4OpRayleigh (check)
  kLABS, // light photon absorption
         // it is turned on when Cerenkov process is turned on
	 // --> may be removed from the enum
           // G3 default value: 0	     
	   // G4 process: ?? G4PhotoAbsorption, G4OpAbsorption (check)
  kSYNC, // synchrotron radiation in magnetic field	   
           // G3 default value: 0	     
	   // G4 process: ?? G4SynchrotronRadiation (check)
*/
  kNoG3Flags
};

enum TG3FlagValue {
// in G3 the control flag values meaning can be different for
// different processes, but for most of them is:
//   0  process is not activated
//   1  process is activated WITH generation of secondaries
//   2  process is activated WITHOUT generation of secondaries
// if process does not generate seconadaries => 1 same as 2
//
// Exceptions:
//   MULS:  also 3
//   LOSS:  also 3, 4 
//   RAYL:  only 0,1
//   HADR:  may be > 2
//
  kUnset      = -1, 
  kInActivate = 0, 
  kActivate   = 1,
  kActivate2  = 2   
}; 

#endif //TG3_FLAG_H

