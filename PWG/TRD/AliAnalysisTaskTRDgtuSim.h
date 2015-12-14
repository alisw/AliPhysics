#ifndef ALIANALYSISTASKTRDGTUSIM_H
#define ALIANALYSISTASKTRDGTUSIM_H

#include "AliAnalysisTaskSE.h"

class TH1;
class AliTRDgtuSim;

class AliAnalysisTaskTRDgtuSim : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskTRDgtuSim(const char *name = "");
  virtual ~AliAnalysisTaskTRDgtuSim();

  virtual Bool_t Notify();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(const Option_t*);

  void Check(Int_t label, Int_t labelRef);

  Int_t GetDeltaAlpha() const { return fDeltaAlpha; }
  Int_t GetDeltaY() const { return fDeltaY; }
  Int_t GetTrackletLabel() const { return fTrackletLabel; }
  Int_t GetLabel() const { return fLabel; }

  void SetDeltaAlpha(Int_t deltaAlpha) { fDeltaAlpha = deltaAlpha; }
  void SetDeltaY(Int_t deltaY) { fDeltaY = deltaY; }
  void SetTrackletLabel(Int_t trackletLabel) { fTrackletLabel = trackletLabel; }
  void SetLabel(Int_t label) { fLabel = label; }

 protected:
  // output objects
  TList  *fOutputList;

  TH1 *fHistStat;
  TH1 *fHistDeltaA;
  TH1 *fHistDeltaB;
  TH1 *fHistDeltaC;

  AliTRDgtuSim *fGtuSim;

  Int_t fTrackletLabel; // required label of the tracklets used for GTU
			// raw tracklets (-2) by default
  Int_t fLabel; // label used for re-simulated GTU tracks
  Int_t fDeltaY; // y window for GTU
  Int_t fDeltaAlpha; // alpha window for GTU
  Bool_t fLimitNoTracklets; // enable limitation on tracklet number
  Int_t fMaxNoTracklets; // maximum no of tracklets (if enabled)

 private:
  AliAnalysisTaskTRDgtuSim(const AliAnalysisTaskTRDgtuSim &rhs);
  AliAnalysisTaskTRDgtuSim& operator=(const AliAnalysisTaskTRDgtuSim &rhs);

  ClassDef(AliAnalysisTaskTRDgtuSim, 1);
};

#endif
