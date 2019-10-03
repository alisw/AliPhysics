#ifndef ALIANALYSISTASKONTHEFLYALIAODDIMUON_H
#define ALIANALYSISTASKONTHEFLYALIAODDIMUON_H

#include "AliAnalysisTaskSE.h"

///
/// Task to create (on-the-fly) from AODs the now deprecated AliAODDimuon object
///
/// It is only intended to ease the transition to the direct usage of muon tracks
/// instead of relying on this pre-computed object.
///
/// \author L. Aphecetche
///

class TClonesArray;

class AliAnalysisTaskOnTheFlyAliAODDimuon : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskOnTheFlyAliAODDimuon();
  virtual ~AliAnalysisTaskOnTheFlyAliAODDimuon();
  
  virtual void UserExec(Option_t* opt);

private:
  TClonesArray* fDimuons; //!<! (transient) storage for AliAODDimuon objects
  
  ClassDef(AliAnalysisTaskOnTheFlyAliAODDimuon,1) /// on the fly creation of (deprecated) dimuon object
};

#endif