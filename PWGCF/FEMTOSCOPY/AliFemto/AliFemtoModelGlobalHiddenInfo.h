///
/// \file AliFemtoModelGlobalHiddenInfo.h
///

#ifndef ALIFEMTOMODELGLOBALHIDDENINFO_H
#define ALIFEMTOMODELGLOBALHIDDENINFO_H

#include "AliFemtoModelHiddenInfo.h"

/// \class AliFemtoModelGlobalHiddenInfo
/// \brief The hidden info for model calculations
///
/// Stores information needed for the weight generation - the true simulated
/// momenta, freeze-out coordinates from model and particle PID and global
/// creation point
///
class AliFemtoModelGlobalHiddenInfo : public AliFemtoModelHiddenInfo {
public:
  AliFemtoModelGlobalHiddenInfo();
  AliFemtoModelGlobalHiddenInfo(const AliFemtoModelGlobalHiddenInfo &aInfo);
  virtual ~AliFemtoModelGlobalHiddenInfo();

  AliFemtoModelGlobalHiddenInfo& operator=(const AliFemtoModelGlobalHiddenInfo& aInfo);

  AliFemtoThreeVector *GetGlobalEmissionPoint() const;
  void SetGlobalEmissionPoint(const AliFemtoThreeVector& aPos);
  void SetGlobalEmissionPoint(Double_t aRx, Double_t aRy, Double_t aRz);

  virtual AliFemtoHiddenInfo* Clone() const;  /// Crate a copy of all hidden information

protected:

  virtual AliFemtoHiddenInfo* GetParticleHiddenInfo() const;

  AliFemtoThreeVector *fGlobalEmissionPoint;
};
//_______________________________________
inline AliFemtoHiddenInfo* AliFemtoModelGlobalHiddenInfo::Clone() const {
  // return exact copy of this hidden info
  return GetParticleHiddenInfo();
}

#endif
