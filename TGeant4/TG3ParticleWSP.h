// $Id$
// Category: physics

#ifndef TG3_PARTICLE_WSP_H
#define TG3_PARTICLE_WSP_H

enum TG3ParticleWSP
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
   
#endif //TG3_PARTICLE_WSP_H

