///
/// \file  AliFemtoTrio.h
/// \author Jeremi Niedziela

#ifndef AliFemtoTrio_H
#define AliFemtoTrio_H

#include "AliFemtoParticle.h"
#include "AliFemtoTypes.h"

class AliFemtoTrio {
public:
  AliFemtoTrio();
  ~AliFemtoTrio();

  // track Gets:
  inline AliFemtoParticle* Track1(){return fTrack1;}
  inline AliFemtoParticle* Track2(){return fTrack2;}
  inline AliFemtoParticle* Track3(){return fTrack3;}
  
  // track Sets:
  inline void SetTrack1(AliFemtoParticle *track){fTrack1 = track;}
  inline void SetTrack2(AliFemtoParticle *track){fTrack2 = track;}
  inline void SetTrack3(AliFemtoParticle *track){fTrack3 = track;}

  double MInv(); // returns invariant mass of the trio
  
private:
  AliFemtoParticle* fTrack1;  // first particle in the trio
  AliFemtoParticle* fTrack2;  // second particle in the trio
  AliFemtoParticle* fTrack3;  // third particle in the trio
};

#endif
