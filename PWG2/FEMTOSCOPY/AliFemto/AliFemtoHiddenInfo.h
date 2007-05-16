////////////////////////////////////////////////////////////////////////////////
/// AliFemtoHiddenInfo - pure virtual base class for the hidden info         ///
/// Hidden info stores additional information, which is not in a standard    ///
/// track. 
////////////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoHiddenInfo_hh
#define AliFemtoHiddenInfo_hh

#include "AliFemtoTypes.h"

class AliFemtoHiddenInfo{

 public:
  AliFemtoHiddenInfo(){/* no-op */};
  virtual ~AliFemtoHiddenInfo(){/* no-op */};
  
  // !!! MANDATORY !!!
  // --- Copy the hidden info from AliFemtoTrack to AliFemtoParticle
  virtual AliFemtoHiddenInfo* Clone() const;
  
 protected:
  virtual AliFemtoHiddenInfo* GetParticleHiddenInfo() const =0;

};
//_______________________________________
inline AliFemtoHiddenInfo* AliFemtoHiddenInfo::Clone() const{
  // return exact copy of this hidden info
  return GetParticleHiddenInfo();
}

#endif
