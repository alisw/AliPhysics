#ifndef ALIANALYSISTASKTRACKINGSYSPROPAGATION_H
#define ALIANALYSISTASKTRACKINGSYSPROPAGATION_H

////////////////////////////////////////////////////////////////////////////
// AliAnalysisTask for Tracking Systematics (matching efficiency          //
// + tracking efficiency) porpagation at the D-meson level                //
// (including daughter's kinematics)                                      //
////////////////////////////////////////////////////////////////////////////

/* $Id$ */

class TString;
class TH1F;
class TH2F;

#include "AliAnalysisTaskSE.h"
#include "AliRDHFCuts.h"

class AliAnalysisTaskTrackingSysPropagation : public AliAnalysisTaskSE {
 public:
  enum DecChannel {kDplustoKpipi,kD0toKpi,kDstartoKpipi,kDstoKKpi};

  AliAnalysisTaskTrackingSysPropagation();
  AliAnalysisTaskTrackingSysPropagation(DecChannel ch, AliRDHFCuts* cuts, TH1F *HistMESys, TH1F *HistTrEffSys);
  virtual ~AliAnalysisTaskTrackingSysPropagation();
    
  virtual void   UserCreateOutputObjects();
  virtual void   Init();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(Option_t *);

  void SetAODMismatchProtection(Int_t opt=1) {fAODProtection=opt;}
  void SetMaximumPt(Double_t maxpt) {fMaxPt = maxpt;}
    
  DecChannel GetDecayChannel()const {return fDecayChannel;}

    
 private:
  AliAnalysisTaskTrackingSysPropagation(const AliAnalysisTaskTrackingSysPropagation &source);
  AliAnalysisTaskTrackingSysPropagation& operator=(const AliAnalysisTaskTrackingSysPropagation &source);
    
    
  TString fPartName;            // string for particle name
  TList *fOutput;               //! tlist with output
    
  TH1F *fHistNEvents;           //! histo with number of events
  TH1F *fHistMESyst;            /// histo with match. eff. systematics vs pt (need to be passed as input)
  TH1F *fHistTrEffSyst;         /// histo with track. eff. systematics vs pt (need to be passed as input)
  TH2F *fhPtDauVsD;             //! histo with Pt daughters vs pt candidate
  TH2F *fhSystMatchEffD;        //! histo with systematic uncertainty on the candidate
    
  DecChannel fDecayChannel;     //identify the decay channel
    
  Int_t fPDGcode;
  Int_t fAODProtection;         /// flag to activate protection against AOD-dAOD mismatch.
    
  Double_t fMaxPt;              /// max pt in the outputs histos
  AliRDHFCuts* fAnalysisCuts;   /// cuts
    
  ClassDef(AliAnalysisTaskTrackingSysPropagation, 2);
};

#endif
