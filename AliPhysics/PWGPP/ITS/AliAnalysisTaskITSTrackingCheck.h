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
class AliESDtrackCuts;
class AliTriggerConfiguration;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskITSTrackingCheck : public AliAnalysisTaskSE 
{
 public:
    AliAnalysisTaskITSTrackingCheck();
    AliAnalysisTaskITSTrackingCheck(const char *name);
  virtual ~AliAnalysisTaskITSTrackingCheck(); 
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  Bool_t         GetReadMC() const { return fReadMC; }
  void           SetReadMC(Bool_t flag=kTRUE) { fReadMC=flag; }
  void           SetUsePhysSel() { fUsePhysSel=kTRUE; }
  void           SetRequireSPDvtx(Bool_t flag=kTRUE) { fRequireSPDvtx=flag; }
  void           SetRejPileupSPD() { fRejPileupSPD=kTRUE; }
  void           SetReadRPLabels(Bool_t flag=kTRUE) { fReadRPLabels=flag; }
  void           SetFillNtuples(Bool_t flag=kTRUE) { fFillNtuples=flag; }
  void           SetUseITSSAforNtuples(Bool_t flag=kTRUE) { fUseITSSAforNtuples=flag; }
  void           SetESDtrackCutsTPC(AliESDtrackCuts *c) { fESDtrackCutsTPC=c; }
  void           SetESDtrackCutsITSTPC(AliESDtrackCuts *c) { fESDtrackCutsITSTPC=c; }
  void           SetMultiplicityRange(Int_t min,Int_t max) { fMinMult=min; fMaxMult=max; }
  void           SetCheckSDDIsIn(Int_t check=0) { fCheckSDDIsIn=check; }
  void           SetTriggerClass(TString tclass="") { fTriggerClass=tclass; }
  void           SetOCDBPath(TString path="") { fOCDBPath=path; }
  void           SetNITSLayers(Int_t n=6) { fNITSLayers=n; }
  void           SetUsePtBinsForBox() { fUsePtBinsForBox=kTRUE; }

 protected:
  Int_t        fNITSLayers; // number of layers (6 for current, 7 for upgrade)
  Bool_t       fReadMC; // read Monte Carlo
  Bool_t       fReadRPLabels; // read MC labels from ITS.RecPoints
  Bool_t       fFillNtuples; // fill expert ntuples
  Bool_t       fUseITSSAforNtuples; // fill expert ntuples with ITSSA tracks
  Bool_t       fUsePhysSel; // use AliPhysicsSelection
  Bool_t       fRequireSPDvtx; // check for SPD vtx to be reconstructed
  Bool_t       fRejPileupSPD; // reject pileup events based on SPD vertex
  Int_t        fCheckSDDIsIn; // check for ITSSDD in the trigger cluster: 0 no check; +1 only wSDD; -1 only woSDD
  AliESDEvent  *fESD;    // ESD object
  Int_t        fMinMult; // minimum multiplicity
  Int_t        fMaxMult; // maximum multiplicity
  TString      fTriggerClass; // trigger class to be inspected
  AliTriggerConfiguration *fTrigConfig; // trigger configuration (read from OCDB)
  Bool_t       fUsePtBinsForBox; // to use special pt binning
  TString      fOCDBPath; // to the OCDB
  TList        *fOutput; //! list send on output slot 0
  TH1F         *fHistNEvents; //! output hist
  TH1F         *fHistNEventsFrac; //! output hist
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
  TH1F         *fHistClusterMapModuleITSMIokInAcc; //! output hist
  TH1F         *fHistClusterMapModuleITSMIbadInAcc; //! output hist
  TH1F         *fHistClusterMapModuleITSMInoclsInAcc; //! output hist
  TH1F         *fHistNClustersMapModule; //! output hist
  TH1F         *fHistZatSPDouter0ok; //! output hist
  TH1F         *fHistZatSPDouter1ok; //! output hist
  TH1F         *fHistZatSPDouter2ok; //! output hist
  TH1F         *fHistZatSPDouter3ok; //! output hist
  TH1F         *fHistZatSPDouter0notok; //! output hist
  TH1F         *fHistZatSPDouter1notok; //! output hist
  TH1F         *fHistZatSPDouter2notok; //! output hist
  TH1F         *fHistZatSPDouter3notok; //! output hist
  TH1F         *fHistxlocSDDok; //! output hist
  TH1F         *fHistzlocSDDok; //! output hist
  TH2F         *fHistxlocVSmodSDDok; //! output hist
  TH1F         *fHistxlocSDDall; //! output hist
  TH1F         *fHistzlocSDDall; //! output hist
  TH1F         *fHistPhiTPCInAcc; //! output hist
  TH1F         *fHistEtaTPCInAcc; //! output hist
  TH1F         *fHistPtTPC; //! output hist
  TH1F         *fHistPtTPCInAcc; //! output hist
  TH1F         *fHistPtTPCInAccTOFbc0; //! output hist
  TH1F         *fHistPtTPCInAccwSDD; //! output hist
  TH1F         *fHistPtTPCInAccTOFbc0wSDD; //! output hist
  TH1F         *fHistPtTPCInAccwoSDD; //! output hist
  TH1F         *fHistPtTPCInAccTOFbc0woSDD; //! output hist
  TH1F         *fHistPtTPCInAccMCtwoSPD; //! output hist
  TH1F         *fHistPtTPCInAccMConeSPD; //! output hist
  TH2F         *fHistdEdxVSPtTPCInAcc; //! output hist
  TH2F         *fHistdEdxVSPtITSTPCsel; //! output hist
  TH2F         *fHistTPCclsVSPtTPCInAcc; //! output hist
  TH2F         *fHistTPCclsVSPtITSMISPDInAcc; //! output hist
  TH2F         *fHistPtVSphiTPCInAcc; //! output hist
  TH1F         *fHistPtTPCInAccNoTRDout; //! output hist
  TH1F         *fHistPtTPCInAccNoTOFout; //! output hist
  TH1F         *fHistPtTPCInAccWithPtTPCAtInnerWall; //! output hist
  TH1F         *fHistPtTPCInAccWithPtTPCAtVtx; //! output hist
  TH2F         *fHistDeltaPtTPC; //! output hist
  TH1F         *fHistPtTPCInAccP; //! output hist
  TH1F         *fHistPtTPCInAccS; //! output hist
  TH1F         *fHistPtTPCInAccPfromStrange; //! output hist
  TH1F         *fHistPtTPCInAccSfromStrange; //! output hist
  TH1F         *fHistPtTPCInAccSfromMat; //! output hist
  TH1F         *fHistPtITSMI2; //! output hist
  TH1F         *fHistPtITSMI3; //! output hist
  TH1F         *fHistPtITSMI4; //! output hist
  TH1F         *fHistPtITSMI5; //! output hist
  TH1F         *fHistPtITSMI6; //! output hist
  TH1F         *fHistPtITSMI7; //! output hist
  TH1F         *fHistPtITSMISPD; //! output hist
  TH1F         *fHistPtITSMIoneSPD; //! output hist
  TH1F         *fHistPtITSMItwoSPD; //! output hist
  TH1F         *fHistPtITSMI2InAcc; //! output hist
  TH1F         *fHistPtITSMI3InAcc; //! output hist
  TH1F         *fHistPtITSMI4InAcc; //! output hist
  TH1F         *fHistPtITSMI5InAcc; //! output hist
  TH1F         *fHistPtITSMI6InAcc; //! output hist
  TH1F         *fHistPtITSMI7InAcc; //! output hist
  TH1F         *fHistPtITSMISPDInAcc; //! output hist
  TH1F         *fHistPtITSMIoneSPDInAcc; //! output hist
  TH1F         *fHistPtITSMItwoSPDInAcc; //! output hist
  TH1F         *fHistPtITSMI2InAccTOFbc0; //! output hist
  TH1F         *fHistPtITSMI3InAccTOFbc0; //! output hist
  TH1F         *fHistPtITSMI4InAccTOFbc0; //! output hist
  TH1F         *fHistPtITSMI5InAccTOFbc0; //! output hist
  TH1F         *fHistPtITSMI6InAccTOFbc0; //! output hist
  TH1F         *fHistPtITSMISPDInAccTOFbc0; //! output hist
  TH1F         *fHistPtITSMIoneSPDInAccTOFbc0; //! output hist
  TH1F         *fHistPtITSMI2InAccwSDD; //! output hist
  TH1F         *fHistPtITSMI3InAccwSDD; //! output hist
  TH1F         *fHistPtITSMI4InAccwSDD; //! output hist
  TH1F         *fHistPtITSMI5InAccwSDD; //! output hist
  TH1F         *fHistPtITSMI6InAccwSDD; //! output hist
  TH1F         *fHistPtITSMISPDInAccwSDD; //! output hist
  TH1F         *fHistPtITSMIoneSPDInAccwSDD; //! output hist
  TH1F         *fHistPtITSMI2InAccTOFbc0wSDD; //! output hist
  TH1F         *fHistPtITSMI3InAccTOFbc0wSDD; //! output hist
  TH1F         *fHistPtITSMI4InAccTOFbc0wSDD; //! output hist
  TH1F         *fHistPtITSMI5InAccTOFbc0wSDD; //! output hist
  TH1F         *fHistPtITSMI6InAccTOFbc0wSDD; //! output hist
  TH1F         *fHistPtITSMISPDInAccTOFbc0wSDD; //! output hist
  TH1F         *fHistPtITSMIoneSPDInAccTOFbc0wSDD; //! output hist
  TH1F         *fHistPtITSMI2InAccwoSDD; //! output hist
  TH1F         *fHistPtITSMI3InAccwoSDD; //! output hist
  TH1F         *fHistPtITSMI4InAccwoSDD; //! output hist
  TH1F         *fHistPtITSMI5InAccwoSDD; //! output hist
  TH1F         *fHistPtITSMI6InAccwoSDD; //! output hist
  TH1F         *fHistPtITSMISPDInAccwoSDD; //! output hist
  TH1F         *fHistPtITSMIoneSPDInAccwoSDD; //! output hist
  TH1F         *fHistPtITSMI2InAccTOFbc0woSDD; //! output hist
  TH1F         *fHistPtITSMI3InAccTOFbc0woSDD; //! output hist
  TH1F         *fHistPtITSMI4InAccTOFbc0woSDD; //! output hist
  TH1F         *fHistPtITSMI5InAccTOFbc0woSDD; //! output hist
  TH1F         *fHistPtITSMI6InAccTOFbc0woSDD; //! output hist
  TH1F         *fHistPtITSMISPDInAccTOFbc0woSDD; //! output hist
  TH1F         *fHistPtITSMIoneSPDInAccTOFbc0woSDD; //! output hist
  TH1F         *fHistPtITSMIoneSPDInAccShared; //! output hist
  TH1F         *fHistPtITSMIoneSPDInAccSharedSPD; //! output hist
  TH1F         *fHistPtITSMISPD1InAccShared; //! output hist
  TH1F         *fHistPtITSMISPD2InAccShared; //! output hist
  TH1F         *fHistPtITSMIoneSPDInAccSharedFake; //! output hist
  TH1F         *fHistPtITSMIoneSPDInAccSharedSPDFake; //! output hist
  TH1F         *fHistPhiITSMI2InAcc; //! output hist
  TH1F         *fHistPhiITSMI3InAcc; //! output hist
  TH1F         *fHistPhiITSMI4InAcc; //! output hist
  TH1F         *fHistPhiITSMI5InAcc; //! output hist
  TH1F         *fHistPhiITSMI6InAcc; //! output hist
  TH1F         *fHistPhiITSMI7InAcc; //! output hist
  TH1F         *fHistPhiITSMISPDInAcc; //! output hist
  TH1F         *fHistPhiITSMIoneSPDInAcc; //! output hist
  TH1F         *fHistPhiITSMItwoSPDInAcc; //! output hist
  TH1F         *fHistEtaITSMI2InAcc; //! output hist
  TH1F         *fHistEtaITSMI3InAcc; //! output hist
  TH1F         *fHistEtaITSMI4InAcc; //! output hist
  TH1F         *fHistEtaITSMI5InAcc; //! output hist
  TH1F         *fHistEtaITSMI6InAcc; //! output hist
  TH1F         *fHistEtaITSMI7InAcc; //! output hist
  TH1F         *fHistEtaITSMISPDInAcc; //! output hist
  TH1F         *fHistEtaITSMIoneSPDInAcc; //! output hist
  TH1F         *fHistEtaITSMItwoSPDInAcc; //! output hist
  TH1F         *fHistPtITSMI2InAccFake; //! output hist
  TH1F         *fHistPtITSMI3InAccFake; //! output hist
  TH1F         *fHistPtITSMI4InAccFake; //! output hist
  TH1F         *fHistPtITSMI5InAccFake; //! output hist
  TH1F         *fHistPtITSMI6InAccFake; //! output hist
  TH1F         *fHistPtITSMI7InAccFake; //! output hist
  TH1F         *fHistPtITSMISPDInAccFake; //! output hist
  TH1F         *fHistPtITSMIoneSPDInAccFake; //! output hist
  TH1F         *fHistPtITSMItwoSPDInAccFake; //! output hist
  TH1F         *fHistPtITSMIoneSPDthreeSDDSSDInAcc; //! output hist
  TH1F         *fHistPtITSTPCsel; //! output hist
  TH1F         *fHistPtITSTPCselTOFbc0; //! output hist
  TH1F         *fHistPtITSTPCselwSDD; //! output hist
  TH1F         *fHistPtITSTPCselTOFbc0wSDD; //! output hist
  TH1F         *fHistPtITSTPCselwoSDD; //! output hist
  TH1F         *fHistPtITSTPCselTOFbc0woSDD; //! output hist
  TH1F         *fHistPtITSTPCselP; //! output hist
  TH1F         *fHistPtITSTPCselS; //! output hist
  TH1F         *fHistPtITSTPCselFake; //! output hist
  TH1F         *fHistPtITSTPCselPfromStrange; //! output hist
  TH1F         *fHistPtITSTPCselSfromStrange; //! output hist
  TH1F         *fHistPtITSTPCselSfromMat; //! output hist
  TH1F         *fHistPtITSMI2InAccP; //! output hist
  TH1F         *fHistPtITSMI3InAccP; //! output hist
  TH1F         *fHistPtITSMI4InAccP; //! output hist
  TH1F         *fHistPtITSMI5InAccP; //! output hist
  TH1F         *fHistPtITSMI6InAccP; //! output hist
  TH1F         *fHistPtITSMI7InAccP; //! output hist
  TH1F         *fHistPtITSMISPDInAccP; //! output hist
  TH1F         *fHistPtITSMIoneSPDInAccP; //! output hist
  TH1F         *fHistPtITSMItwoSPDInAccP; //! output hist
  TH1F         *fHistPtITSMI2InAccS; //! output hist
  TH1F         *fHistPtITSMI3InAccS; //! output hist
  TH1F         *fHistPtITSMI4InAccS; //! output hist
  TH1F         *fHistPtITSMI5InAccS; //! output hist
  TH1F         *fHistPtITSMI6InAccS; //! output hist
  TH1F         *fHistPtITSMI7InAccS; //! output hist
  TH1F         *fHistPtITSMISPDInAccS; //! output hist
  TH1F         *fHistPtITSMIoneSPDInAccS; //! output hist
  TH1F         *fHistPtITSMItwoSPDInAccS; //! output hist
  TH1F         *fHistPtITSMIokbadoutinz6; //! output hist
  TH1F         *fHistPtITSMIokbadoutinz4InAcc; //! output hist
  TH1F         *fHistPtITSMIokbadoutinz5InAcc; //! output hist
  TH1F         *fHistPtITSMIokbadoutinz6InAcc; //! output hist
  TH1F         *fHistPhiITSMIokbadoutinz6InAcc; //! output hist
  TH1F         *fHistRProdVtxInAccP; //! output hist
  TH1F         *fHistRProdVtxInAccS; //! output hist
  TH1F     *fHistd0rphiTPCInAccP150200; //! output hist
  TH1F     *fHistd0rphiTPCInAccP500700; //! output hist
  TH1F     *fHistd0rphiTPCInAccP10001500; //! output hist
  TH1F     *fHistd0rphiTPCInAccS150200; //! output hist
  TH1F     *fHistd0rphiTPCInAccS500700; //! output hist
  TH1F     *fHistd0rphiTPCInAccS10001500; //! output hist
  TH1F     *fHistd0rphiITSMISPDInAccP150200; //! output hist
  TH1F     *fHistd0rphiITSMISPDInAccP500700; //! output hist
  TH1F     *fHistd0rphiITSMISPDInAccP10001500; //! output hist
  TH1F     *fHistd0rphiITSMISPDInAccS150200; //! output hist
  TH1F     *fHistd0rphiITSMISPDInAccS500700; //! output hist
  TH1F     *fHistd0rphiITSMISPDInAccS10001500; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccP150200; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccP350450; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccP500700; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccP10001500; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccP25004000; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccP40008000; //! output hist
  TH1F     *fHistd0zITSMIoneSPDInAccP150200; //! output hist
  TH1F     *fHistd0zITSMIoneSPDInAccP500700; //! output hist
  TH1F     *fHistd0zITSMIoneSPDInAccP10001500; //! output hist
  TH2F     *fHistd0zVSetaTPCInAccP10001500; //! output hist
  TH2F     *fHistd0rphiVSphiITSMIoneSPDInAccP10001500; //! output hist
  TH2F     *fHistd0rphiVSetaITSMIoneSPDInAccP10001500; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS150200; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS350450; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS500700; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS500700from22; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS500700from211; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS500700from310; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS500700from321; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS500700from3122; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS10001500; //! output hist  
  TH1F     *fHistd0rphiITSMIoneSPDInAccS25004000; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS40008000; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS150200fromStrange; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS150200fromMat; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS350450fromStrange; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS350450fromMat; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS500700fromStrange; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS500700fromMat; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS10001500fromStrange; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS10001500fromMat; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS25004000fromStrange; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS25004000fromMat; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS40008000fromStrange; //! output hist
  TH1F     *fHistd0rphiITSMIoneSPDInAccS40008000fromMat; //! output hist
  TH1F     *fHistd0zITSMIoneSPDInAccS150200; //! output hist
  TH1F     *fHistd0zITSMIoneSPDInAccS500700; //! output hist
  TH1F     *fHistd0zITSMIoneSPDInAccS10001500; //! output hist
  TH1F     *fHistPDGMoth; //! output hist
  TH1F     *fHistPDGMoth150200; //! output hist
  TH1F     *fHistPDGMoth500700; //! output hist
  TH1F     *fHistPDGMoth10001500; //! output hist
  TH1F     *fHistPDGTrk; //! output hist
  TH1F     *fHistITSRedChi2NonFakePt02; //! output hist
  TH1F     *fHistITSRedChi2FakePt02; //! output hist
  TH1F     *fHistITSRedChi2NonFakePt05; //! output hist
  TH1F     *fHistITSRedChi2FakePt05; //! output hist
  TH1F     *fHistITSRedChi2NonFakePt1; //! output hist
  TH1F     *fHistITSRedChi2FakePt1; //! output hist
  TNtuple      *fNtupleESDTracks; //! output ntuple
  TNtuple      *fNtupleITSAlignExtra; //! output ntuple
  TNtuple      *fNtupleITSAlignSPDTracklets; //! output ntuple
  Int_t         fCountsPerPtBin[11]; // track per pt bin
  AliESDtrackCuts *fESDtrackCutsTPC; // cuts for TPC track
  AliESDtrackCuts *fESDtrackCutsITSTPC; // cuts for TPC+ITS track

 private:    

  AliAnalysisTaskITSTrackingCheck(const AliAnalysisTaskITSTrackingCheck&); // not implemented
  AliAnalysisTaskITSTrackingCheck& operator=(const AliAnalysisTaskITSTrackingCheck&); // not implemented
  

  Int_t NumberOfITSClustersMC(Int_t label,Int_t nModules=2198) const;
  Int_t NumberOfITSClusters(Int_t idet,Float_t &xloc) const;
  Double_t ParticleImpParMC(TParticle *part,AliESDVertex *vert,Double_t bzT) const;
  Bool_t SelectPt(Double_t pt);
  Int_t MakeITSflag(AliESDtrack *track) const;
  Bool_t IsSelectedCentrality() const;
  void FillNClustersModuleMap();
  Int_t NPointsInnerBarrel(AliESDtrack *track) const;
  Bool_t ConditionSPD(AliESDtrack *track) const;
  Bool_t ConditionSPDone(AliESDtrack *track) const;
  Bool_t ConditionSPDtwo(AliESDtrack *track) const;



  ClassDef(AliAnalysisTaskITSTrackingCheck,16); // ITS tracks analysis
};

#endif
