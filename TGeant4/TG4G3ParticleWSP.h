// $Id$
// Category: physics

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

