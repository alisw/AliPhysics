// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Enum TG4G3ParticleWSP
// ---------------------
// The particles types which a special process (cuts, controls)
// is applicable for.

#ifndef TG4_G3_PARTICLE_WSP_H
#define TG4_G3_PARTICLE_WSP_H

enum TG4G3ParticleWSP
// particle with special process
{
  kGamma,           // kPAIR, kCOMP, kPHOT, kPHIS
  kElectron,        // kDRAY, kBREM, kMULS, kLOSS
  kEplus,           // kDRAY, kBREM, kMULS, kLOSS, kANNI
  kNeutralHadron,   // kHADR 
  kChargedHadron,   // kDRAY, kMULS, kLOSS, kHADR,
  kMuon,            // kDRAY, kBREM, kMULS, kLOSS, kMUNU
  kAny,             // kDCAY
  kNofParticlesWSP
};
   
#endif //TG4_G3_PARTICLE_WSP_H

