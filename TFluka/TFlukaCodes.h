#ifndef TFLUKACODES
#define TFLUKACODES 
typedef enum {
    kNoProcess         =   0,
    kKASKAD            =   1,  // KASKAD code first digit
    kKASKADelarecoil   =  10,  // elastic interaction recoil
    kKASKADinelarecoil =  11,  // inelastic interaction recoil   
    kKASKADstopping    =  12,  // stopping particle
    kKASKADpseudon     =  13,  // pseudo-neutron deposition
    kKASKADescape      =  14,  // escape
    kKASKADtimekill    =  15,  // time kill 
    kKASKADboundary    =  19,  // boundary crossing
    kKASKADnelint      = 100,  // elastic   interaction 
    kKASKADinelint     = 101,  // inelastic interaction 
    kKASKADdecay       = 102,  // particle decay  
    kKASKADdray        = 103,  // delta ray  generation 
    kKASKADpair        = 104,  // pair production
    kKASKADbrems       = 105,  // bremsstrahlung
    kEMFSCO            =   2,  // EMFSCO code first digit
    kEMSCOlocaledep    =  20,  // local energy deposition (i.e. photoelectric)
    kEMFSCOstopping1   =  21,  // below user-defined cut-off
    kEMFSCOstopping2   =  22,  // below user cut-off
    kEMFSCOescape      =  23,  // escape  
    kEMFSCOtimekill    =  24,  // time kill
    kEMFSCOboundary    =  29,  // boundary crossing
    kEMFSCObrems       = 208,  // bremsstrahlung
    kEMFSCOmoller      = 210,  // Moller
    kEMFSCObhabha      = 212,  // Bhabha
    kEMFSCOanniflight  = 214,  // in-flight annihilation
    kEMFSCOannirest    = 215,  // annihilation at rest
    kEMFSCOpair        = 217,  // pair production
    kEMFSCOcompton     = 219,  // Compton scattering
    kEMFSCOphotoel     = 221,  // photoelectric effect 
    kEMFSCOrayleigh    = 225,  // Rayleigh scattering  
    kKASNEU            =   3,  // KASNEU code first digit
    kKASNEUtargrecoil  =  30,  // target recoil
    kKASNEUstopping    =  31,  // neutron below threshold
    kKASNEUescape      =  32,  // escape 
    kKASNEUtimekill    =  33,  // time kill
    kKASNEUboundary    =  39,  // boundary crossing
    kKASNEUhadronic    = 300,  // neutron interaction
    kKASHEA            =   4,  // KASHEA code first digit
    kKASHEAescape      =  40,  // escape
    kKASHEAtimekill    =  41,  // time kill 
    kKASHEAboundary    =  49,  // boundary crossing
    kKASHEAdray        = 400,  // delta ray generation
    kKASOPH            =   5,  // KASOPH code first digit
    kKASOPHabsorption  =  50,  // optical photon absorption 
    kKASOPHescape      =  51,  // escape 
    kKASOPHtimekill    =  52,  // time kill
    kKASOPHrefraction  =  59   // boundary crossing (i.e. refraction)
}
FlukaProcessCode_t;	  

typedef enum {
    kNoCaller       =  0,
    kEEDRAW         =  2,      // Stepping called from eedraw 
    kENDRAW         =  3,      // Stepping called from endraw      
    kMGDRAW         =  4,      // Stepping called from mgdraw 
    kSODRAW         =  5,      // Stepping called from sodraw 
    kUSDRAW         =  6,      // Stepping called from usdraw 
    kBXEntering     = 11,      // Stepping called from bxdraw (entering track) 
    kBXExiting      = 12,      // Stepping called from bxdraw (exiting  track) 
    kMGResumedTrack = 40,      // Stepping called from mgdraw (resumed  track) 
    kUSTCKV         = 50       // Stepping called from ustckv 
}
FlukaCallerCode_t;	  

typedef enum {
    kFLUKAcodemin  = -  6,     // minimum particle code used by FLUKA
    kFLUKAcodemax  =  250,     // maximum particle code used by FLUKA
    kFLUKAoptical  = -  1,     // code for optical photon
    kFLUKAelectron =    3,     // electron
    kFLUKApositron =    4,     // positron
    kFLUKAphoton   =    7,     // photon
    kFLUKAmuplus   =   10,     // mu+
    kFLUKAmuminus  =   11      // mu-
}
FlukaParticleCode_t;	  

#endif //TFLUKACODE
