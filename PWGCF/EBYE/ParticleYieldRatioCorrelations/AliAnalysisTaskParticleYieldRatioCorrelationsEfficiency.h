#ifndef AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency_H
#define AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency_H

class TList;
class TH1F;
class TH1D;
class TH2D;
class TH3D;
class TF1;
class AliPIDResponse;

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include "AliEventCuts.h"
class AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency();
  AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency(const char *name);
  virtual ~AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

private:
  AliAODEvent *fAOD;            //! input event
  AliPIDResponse *fPIDResponse; //! pid response object
  TList *fOutputList;           //! output list

  AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency(const AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency &);            // not implemented
  AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency &operator=(const AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency &); // not implemented

  ClassDef(AliAnalysisTaskParticleYieldRatioCorrelationsEfficiency, 1);
};

#endif
