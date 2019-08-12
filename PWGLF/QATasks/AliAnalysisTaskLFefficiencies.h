#ifndef __AliAnalysisTaskLFefficiencies__
#define __AliAnalysisTaskLFefficiencies__

#include "AliAnalysisTaskSE.h"

#include <string>

#include <AliEventCuts.h>
#include <AliPID.h>

class TF1;
class TH1D;
class TH2F;
class TH3D;
class AliFlowTrackCuts;
class AliAODTrack;
class AliVVertex;
class AliPIDResponse;
class TList;

class AliAnalysisTaskLFefficiencies : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskLFefficiencies(TString taskname = "LFefficienciesTask");
  virtual ~AliAnalysisTaskLFefficiencies();

  static bool  HasTOF(AliVTrack *t);

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);

  AliEventCuts  fEventCut;

  bool fUseMCtruthParams;
  static const std::string fPosNeg[2];
  static const int fNcuts;
  static const std::string fCutNames[8];
private:
  AliAnalysisTaskLFefficiencies (const AliAnalysisTaskLFefficiencies &source);
  AliAnalysisTaskLFefficiencies &operator=(const AliAnalysisTaskLFefficiencies &source);

  TList* fOutputList;                                    //!<! Output list
  TH1D* fNumberOfRecoPrimaryTracks;                      //!<! Number of recostructed primary tracks per event
  TH3D* fGeneratedYPhiPt[AliPID::kSPECIESC][2];          //!<! Generated particles
  TH3D* fReconstructedYPhiPt[AliPID::kSPECIESC][2][8];   //!<! Reconstructed particles vs y, Phi and pT, {FB4,FB5,FB5+PID TPC, FB5 + TOF matching, FB5 + PID TOF, FB5 + TOF matching - TOF mismatch, FB5 + TOF matching - TOF mismatch + TOF PID}
  TH3D* fGeneratedEtaPhiPt[AliPID::kSPECIESC][8];        //!<! Generated particles in the eta
  TH3D* fReconstructedEtaPhiPt[AliPID::kSPECIESC][2][8]; //!<! Reconstructed particles vs eta, Phi and pT, {FB4,FB5,FB5+PID TPC, FB5 + TOF matching, FB5 + TOF matching - TOF mismatch, FB5 + TOF matching - TOF mismatch + TOF PID}

  TH2D* fNsigmaTOFvsPt[AliPID::kSPECIESC][2];  //!<! N sigma distribution for tracks passing FB5 + hasTOF without mismatch;

  /// \cond CLASSDEF
  ClassDef(AliAnalysisTaskLFefficiencies, 1);
  /// \endcond
};


#endif /* defined(__AliAnalysisTaskLFefficiencies__) */
