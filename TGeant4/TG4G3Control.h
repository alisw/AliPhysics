// $Id$
// Category: global
//
// Author: I. Hrivnacova
//
// Enum TG4G3Cut
// -------------
// Enumeration for G3 types of physics processes controls.
// The G4 physics processes and particles which the process
// control is applied to are indicated in the comments.
// The process control is mapped to the physics processes
// (using TG4ProcessControlMap) at the moment of process creating
// by the physics constructor; the physics contructor type
// is indicated by the "Physics:" comment.

#ifndef TG4_G3_CONTROL_H
#define TG4_G3_CONTROL_H

enum TG4G3Control
{
  kPAIR, // pair production       
             // G3 default value: 1
             // G4 processes: G4GammaConversion,
             //               G4MuPairProduction/G4IMuPairProduction
	     //               G4LowEnergyGammaConversion
	     // Particles: gamma, mu
	     // Physics:   EM  
	     
  kCOMP, // Compton scattering    
             // G3 default value: 1
             // G4 processes: G4ComptonScattering, 
	     //               G4LowEnergyCompton,
	     //               G4PolarizedComptonScattering
	     // Particles: gamma
	     // Physics:   EM  
	     
  kPHOT, // photo electric effect 
             // G3 default value: 1
             // G4 processes: G4PhotoElectricEffect
	     //               G4LowEnergyPhotoElectric
	     // Particles: gamma
	     // Physics:   EM  

  kPFIS, // photofission          
             // G3 default value: 0
             // G4 process: ??
	     //
	     // Particles: gamma
	     // Physics:   ??  
	     
  kDRAY, // delta-ray              
             // G3 default value: 2
	     // !! G4 treats delta rays in different way
             // G4 processes: G4eIonisation/G4IeIonization,
             //               G4MuIonisation/G4IMuIonization, 
	     //               G4hIonisation/G4IhIonisation
	     // Particles: charged 
	     // Physics:   EM  

  kANNI, // annihilation          
             // G3 default value: 1
             // G4 processes: G4eplusAnnihilation/G4IeplusAnnihilation
	     // Particles: e+ 
	     // Physics:   EM  
	     
  kBREM, // bremsstrahlung        
             // G3 default value: 1
             // G4 processes: G4eBremsstrahlung/G4IeBremsstrahlung,
	     //               G4MuBremsstrahlung/G4IMuBremsstrahlung,
	     //               G4LowEnergyBremstrahlung
	     //               
	     // Particles: e-/e+; mu+/mu- 
	     // Physics:   EM  
	     
  kHADR, // hadronic process      
             // G3 default value: 1
             // G4 processes: all defined by TG4PhysicsConstructorHadron 
	     //               
	     // Particles: hadrons 
	     // Physics:   Hadron
	     
  kMUNU, // muon nuclear interaction 
             // G3 default value: 0
	     // G4 processes: G4MuNuclearInteraction
	     //
	     // Particles: mu
	     // Physics:   Not set
	     
  kDCAY, // decay                 
             // G3 default value: 1
	     // G4 process: G4Decay
	     //
	     // Particles: all which decay is applicable for
	     // Physics:   General
	     
  kLOSS, // energy loss           
             // G3 default value: 2
             // G4 processes: G4eIonisation/G4IeIonization,
             //               G4MuIonisation/G4IMuIonization, 
	     //               G4hIonisation/G4IhIonisation
	     //
	     // Particles: charged 
	     // Physics:   EM  
	     
  kMULS, // multiple scattering   
             // G3 default value: 1
	     // G4 process: G4MultipleScattering/G4IMultipleScattering
	     //
	     // Particles: charged 
	     // Physics:   EM  

  kCKOV, // Cerenkov photon generation
             // G3 default value: 0
	     // G4 process: G4Cerenkov
	     //            
	     // Particles: charged  
	     // Physics:   Optical  
	   
  kRAYL, // Rayleigh scattering
             // G3 default value: 0	     
	     // G4 process: G4OpRayleigh
	     //            
	     // Particles: optical photon  
	     // Physics:   Optical  
	     
  kLABS, // light photon absorption
             // it is turned on when Cerenkov process is turned on
             // G3 default value: 0	     
	     // G4 process: G4OpAbsorption, G4OpBoundaryProcess
	     //
	     // Particles: optical photon  
	     // Physics:   Optical  

  kSYNC, // synchrotron radiation in magnetic field	   
             // G3 default value: 0	     
	     // G4 process: G4SynchrotronRadiation
	     //
	     // Particles: ??
	     // Physics:   Not set  

  kNoG3Controls
};

enum TG4G3ControlValue 
{
// in G3 the process control values meaning can be different for
// different processes, but for most of them is:
//   0  process is not activated
//   1  process is activated WITH generation of secondaries
//   2  process is activated WITHOUT generation of secondaries
// if process does not generate secondaries => 1 same as 2
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

#endif //TG4_G3_CONTROL_H

