////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoModelGlobalHiddenInfo - the hidden info for model calculations     //
// Stores information needed for the weight generation - the true             //
// simulated momenta, freeze-out coordinates from model and particle PID      //
// and global creation point                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOMODELGLOBALHIDDENINFO_H
#define ALIFEMTOMODELGLOBALHIDDENINFO_H

#include <TH1D.h>
#include "AliFemtoTypes.h"
#include "AliFemtoThreeVector.h"
#include "AliFemtoLorentzVector.h"
#include "AliFemtoHiddenInfo.h"
#include "AliFemtoModelHiddenInfo.h"

class AliFemtoModelGlobalHiddenInfo : public AliFemtoModelHiddenInfo{

public:
  AliFemtoModelGlobalHiddenInfo();
  AliFemtoModelGlobalHiddenInfo(const AliFemtoModelGlobalHiddenInfo &aInfo);
  virtual ~AliFemtoModelGlobalHiddenInfo();

  AliFemtoModelGlobalHiddenInfo& operator=(const AliFemtoModelGlobalHiddenInfo& aInfo);

  AliFemtoThreeVector   *GetGlobalEmissionPoint() const;
  void                   SetGlobalEmissionPoint(const AliFemtoThreeVector& aPos);
  void                   SetGlobalEmissionPoint(Double_t aRx, Double_t aRy, Double_t aRz);

// !!! MANDATORY !!!
// --- Copy the hidden info from AliFemtoTrack to AliFemtoParticle
  virtual AliFemtoHiddenInfo* Clone() const;
  
 protected:
  virtual AliFemtoHiddenInfo* GetParticleHiddenInfo() const;

  AliFemtoThreeVector   *fGlobalEmissionPoint;
};
//_______________________________________
inline AliFemtoHiddenInfo* AliFemtoModelGlobalHiddenInfo::Clone() const{
  // return exact copy of this hidden info
  return GetParticleHiddenInfo();
}

#endif
