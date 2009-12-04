#ifndef AliAnalysisTaskITSTrackingCheck_cxx
#define AliAnalysisTaskITSTrackingCheck_cxx

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysiTaskITSTrackingCheck
// AliAnalysisTask to extract from ESD tracks the information on the
// ITS tracking efficiency and resolutions
//
// Author: A.Dainese, andrea.dainese@pd.infn.it
//*************************************************************************

class TNtuple;
class TParticle;
class TH1F;
class AliESDEvent;
class AliESDVertex;
class AliESDfriend;

#include "AliAnalysisTask.h"

class AliAnalysisTaskITSTrackingCheck : public AliAnalysisTask 
{
 public:

  AliAnalysisTaskITSTrackingCheck(const char *name = "AliAnalysisTaskITSTrackingCheck");
  virtual ~AliAnalysisTaskITSTrackingCheck(); 
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  Bool_t         GetReadMC() const { return fReadMC; }
  void           SetReadMC(Bool_t flag=kTRUE) { fReadMC=flag; }
  void           SetReadRPLabels(Bool_t flag=kTRUE) { fReadRPLabels=flag; }
  void           SetFillNtuples(Bool_t flag=kTRUE) { fFillNtuples=flag; }
  void           SetUseITSSAforNtuples(Bool_t flag=kTRUE) { fUseITSSAforNtuples=flag; }

  
 protected:
  Bool_t       fReadMC; // read Monte Carlo
  Bool_t       fReadRPLabels; // read MC labels from ITS.RecPoints
  Bool_t       fFillNtuples; // fill expert ntuples
  Bool_t       fUseITSSAforNtuples; // fill expert ntuples with ITSSA tracks
  AliESDEvent  *fESD;    // ESD object
  AliESDfriend *fESDfriend; // ESD friend object
  TList        *fOutput; //! list send on output slot 0
  TH1F         *fHistNtracks; //! output hist
  TH1F         *fHistNclsITSMI; //! output hist
  TH1F         *fHistNclsITSSA; //! output hist
  TH1F         *fHistNclsITSSAInAcc; //! output hist
  TH1F         *fHistClusterMapITSMI; //! output hist
  TH1F         *fHistClusterMapITSMIok; //! output hist
  TH1F         *fHistClusterMapITSMIbad; //! output hist
  TH1F         *fHistClusterMapITSMIskipped; //! output hist
  TH1F         *fHistClusterMapITSMIoutinz; //! output hist
  TH1F         *fHistClusterMapITSMInorefit; //! output hist
  TH1F         *fHistClusterMapITSMInocls; //! output hist
  TH1F         *fHistClusterMapITSMIokoutinzbad; //! output hist
  TH1F         *fHistClusterMapITSSA; //! output hist
  TH1F         *fHistClusterMapITSSAok; //! output hist
  TH1F         *fHistClusterMapITSSAbad; //! output hist
  TH1F         *fHistClusterMapITSSAskipped; //! output hist
  TH1F         *fHistClusterMapITSSAoutinz; //! output hist
  TH1F         *fHistClusterMapITSSAnorefit; //! output hist
  TH1F         *fHistClusterMapITSSAnocls; //! output hist
  TH1F         *fHistClusterMapITSSAokoutinzbad; //! output hist
  TH1F         *fHistClusterMapITSSAInAcc; //! output hist
  TH1F         *fHistClusterMapITSSAokInAcc; //! output hist
  TH1F         *fHistClusterMapITSSAbadInAcc; //! output hist
  TH1F         *fHistClusterMapITSSAskippedInAcc; //! output hist
  TH1F         *fHistClusterMapITSSAoutinzInAcc; //! output hist
  TH1F         *fHistClusterMapITSSAnorefitInAcc; //! output hist
  TH1F         *fHistClusterMapITSSAnoclsInAcc; //! output hist
  TH1F         *fHistClusterMapITSSAokoutinzbadInAcc; //! output hist
  TH1F         *fHistClusterMapModuleITSSAokInAcc; //! output hist
  TH1F         *fHistClusterMapModuleITSSAbadInAcc; //! output hist
  TH1F         *fHistClusterMapModuleITSSAnoclsInAcc; //! output hist
  TH1F         *fHistPhiTPCInAcc; //! output hist
  TH1F         *fHistPtTPC; //! output hist
  TH1F         *fHistPtTPCInAcc; //! output hist
  TH1F         *fHistPtITSMI2; //! output hist
  TH1F         *fHistPtITSMI3; //! output hist
  TH1F         *fHistPtITSMI4; //! output hist
  TH1F         *fHistPtITSMI5; //! output hist
  TH1F         *fHistPtITSMI6; //! output hist
  TH1F         *fHistPtITSMISPD; //! output hist
  TH1F         *fHistPtITSMI2InAcc; //! output hist
  TH1F         *fHistPtITSMI3InAcc; //! output hist
  TH1F         *fHistPtITSMI4InAcc; //! output hist
  TH1F         *fHistPtITSMI5InAcc; //! output hist
  TH1F         *fHistPtITSMI6InAcc; //! output hist
  TH1F         *fHistPtITSMISPDInAcc; //! output hist
  TH1F         *fHistPtITSMIokbadoutinz6; //! output hist
  TH1F         *fHistPtITSMIokbadoutinz6InAcc; //! output hist
  TH1F         *fHistPhiITSMIokbadoutinz6InAcc; //! output hist
  TNtuple      *fNtupleESDTracks; //! output ntuple
  TNtuple      *fNtupleITSAlignExtra; //! output ntuple
  TNtuple      *fNtupleITSAlignSPDTracklets; //! output ntuple
  Int_t         fCountsPerPtBin[10]; // track per pt bin

 private:    

  AliAnalysisTaskITSTrackingCheck(const AliAnalysisTaskITSTrackingCheck&); // not implemented
  AliAnalysisTaskITSTrackingCheck& operator=(const AliAnalysisTaskITSTrackingCheck&); // not implemented
  

  Int_t NumberOfITSClustersMC(Int_t label) const;
  Double_t ParticleImpParMC(TParticle *part,AliESDVertex *vert,Double_t bzT) const;
  Bool_t SelectPt(Double_t pt);

  ClassDef(AliAnalysisTaskITSTrackingCheck,1); // ITS tracks analysis
};

#endif
