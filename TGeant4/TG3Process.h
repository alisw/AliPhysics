// $Id$
// Category: physics

#ifndef TG3_PROCESS_H
#define TG3_PROCESS_H

//  NAMEC    List  of possible  mechanisms  for  step size  limitation  *
//           filled in GINIT :                                          *
//     DATA MEC/'NEXT','MULS','LOSS','FIEL','DCAY','PAIR','COMP','PHOT' *
//    +        ,'BREM','DRAY','ANNI','HADR','ECOH','EVAP','FISS','ABSO' *
//    +        ,'ANNH','CAPT','EINC','INHE','MUNU','TOFM','PFIS','SCUT' *
//    +        ,'RAYL','PARA','PRED','LOOP','NULL','STOP'/              *
//

enum TG3Process
{
         // process description           //LMEC value
                                          //in /GCTRAK/ common   
  kNEXT,
  kMULS, // multiple scattering   
  kLOSS, // energy loss    
  kFIEL,     
  kDCAY, // decay                 
  kPAIR, // pair production       
  kCOMP, // Compton scattering    
  kPHOT, // photo electric effect 
  kBREM, // bremsstrahlung        
  kDRAY, // delta-ray              
  kANNI, // positron annihilation          
  kHADR, // hadronic process      

  kECOH,
  kEVAP,
  kFISS,
  kABSO,
  kANNH,
  kCAPT,
  kEINC,
  kINHE,
  kMUNU, // muon nuclear interaction 
  kTOFM
  kPFIS, // photofission          
  kSCUT,
  kRAYL, // Rayleigh scattering;          
  kPARA,
  kPRED,
  kLOOP,
  kNULL,
  kSTOP,
  
  kLABS, // light photon absorption;      code: 101
  kLREF, // photon boundary effects;      code: 102
  kCKOV  // Cerenkov photon generation;   code: 105    
  kREFL, // photon reflection;            code: 106
  kREFR, // photon refraction;            code: 107

  kNoG3Process
};

#endif //TG3_PROCESS_H
