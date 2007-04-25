////////////////////////////////////////////////////////////////////////////////
/// AliFemtoHiddenInfo - pure virtual base class for the hidden info         ///
/// Hidden info stores additional information, which is not in a standard    ///
/// track. 
////////////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoHiddenInfo_hh
#define AliFemtoHiddenInfo_hh

#include "Infrastructure/AliFemtoTypes.h"

class AliFemtoHiddenInfo{

public:
  AliFemtoHiddenInfo(){/* no-op */};
  virtual ~AliFemtoHiddenInfo(){/* no-op */};

// !!! MANDATORY !!!
// --- Copy the hidden info from AliFemtoTrack to AliFemtoParticle
  virtual AliFemtoHiddenInfo* getParticleHiddenInfo() const =0;
  virtual AliFemtoHiddenInfo* clone() const;

};
//_______________________________________
inline AliFemtoHiddenInfo* AliFemtoHiddenInfo::clone() const{
  // return exact copy of this hidden info
  return getParticleHiddenInfo();
}

#endif
