#ifndef TFLUKACODES
#define TFLUKACODES 
//
// Enumeration of the constants for the PDG particle IDs.
//

typedef enum {
    kNoProcess         =   0,
    kKASKAD            =   1,  // any KASKAD code / 100
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
    kEMFSCO            =   2,
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
    kKASNEU            =   3,
    kKASNEUtargrecoil  =  30,  // target recoil
    kKASNEUstopping    =  31,  // neutron below threshold
    kKASNEUescape      =  32,  // escape 
    kKASNEUtimekill    =  33,  // time kill
    kKASNEUboundary    =  39,  // boundary crossing
    kKASNEUhadronic    = 300,  // neutron interaction
    kKASHEA            =   4,
    kKASHEAescape      =  40,  // escape
    kKASHEAtimekill    =  41,  // time kill 
    kKASHEAboundary    =  49,  // boundary crossing
    kKASHEAdray        = 400,  // delta ray generation
    kKASOPH            =   5,
    kKASOPHabsorption  =  50,  // optical photon absorption 
    kKASOPHescape      =  51,  // escape 
    kKASOPHtimekill    =  52,  // time kill
    kKASOPHrefraction  =  59   // boundary crossing (i.e. refraction)
}
FlukaProcessCode_t;	  

typedef enum {
    kEEDRAW = 2, kENDRAW = 3, kMGDRAW = 4, kSODRAW = 5, kUSDRAW = 6,
    kBXEntering = 11, kBXExiting = 12,
    kMGResumedTrack = 40,
    kUSTCKV = 50
}
FlukaCallerCode_t;	  

typedef enum {
    kFLUKAoptical  = -1,
    kFLUKAelectron = 3,
    kFLUKApositron = 4,
    kFLUKAphoton   = 7,
    kFLUKAmuplus   = 10,
    kFLUKAmuminus  = 11
}
FlukaParticleCode_t;	  

#endif //TFLUKACODE
