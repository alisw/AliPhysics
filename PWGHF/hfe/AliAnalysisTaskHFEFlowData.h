/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// Flow task class for the ALICE HFE group
//
//
#ifndef ALIANALYSISTASKHFEFLOWDATA_H
#define ALIANALYSISTASKHFEFLOWDATA_H




#include <AliAnalysisTaskSE.h>

class TList;
class AliFlowTrackCuts;
class AliFlowCandidateTrack;
class AliHFEcuts;
class AliHFEpid;
class TH1D;
class TH2D;
class TProfile;
class TProfile2D;
class THnSparse;
class AliHFEpidQAmanager;
class AliFlowEvent;
class AliHFEVZEROEventPlane;

class AliAnalysisTaskHFEFlowData: public AliAnalysisTaskSE {
  
public:
  AliAnalysisTaskHFEFlowData();
  AliAnalysisTaskHFEFlowData(const char *name);
  AliAnalysisTaskHFEFlowData(const AliAnalysisTaskHFEFlowData &ref);
  AliAnalysisTaskHFEFlowData& operator=(const AliAnalysisTaskHFEFlowData &ref);
  virtual void Copy(TObject &o) const;
  virtual ~AliAnalysisTaskHFEFlowData();
  
  virtual void  UserExec(Option_t */*option*/);
  virtual void  UserCreateOutputObjects();

  void SetAODAnalysis(Bool_t aodAnalysis) { fAODAnalysis = aodAnalysis; };
  void SetUseFlagAOD(Bool_t useFlagAOD) { fUseFlagAOD = useFlagAOD; }
  void SetApplyCut(Bool_t applyCut) { fApplyCut = applyCut; }
  void SetFlags(ULong_t flags)          { fFlags = flags; }
  
  AliHFEpid *GetPID() const { return fPID; }
  AliHFEpidQAmanager *GetPIDQAManager() const { return fPIDqa; }

  void SetHFECuts(AliHFEcuts * const cuts) { fHFECuts = cuts; };
  void SetSubEtaGapTPC(Bool_t  subEtaGapTPC) { fSubEtaGapTPC = subEtaGapTPC; };
  void SetEtaGap(Double_t  etaGap) { fEtaGap = etaGap; };
  void SetVZEROEventPlane(Bool_t vzeroEventPlane) { fVZEROEventPlane = vzeroEventPlane; };
  void SetVZEROEventPlaneA(Bool_t vzeroEventPlaneA) { fVZEROEventPlaneA = vzeroEventPlaneA; };
  void SetVZEROEventPlaneC(Bool_t vzeroEventPlaneC) { fVZEROEventPlaneC = vzeroEventPlaneC; };
  
  void SetDebugLevel(Int_t debugLevel) { fDebugLevel = debugLevel;};
  
private:
  TList     *fListHist;         //! TH list
  Bool_t    fAODAnalysis;       // AOD analysis
  Bool_t    fUseFlagAOD;        // Use the preselected AOD track
  Bool_t    fApplyCut;       // Apply the analysis cut for AOD tracks
  ULong_t   fFlags;             // reconstruction AOD status flags 
  
  Bool_t    fVZEROEventPlane;  // Use Event Planes from VZERO
  Bool_t    fVZEROEventPlaneA; // Use Event Planes from VZERO A
  Bool_t    fVZEROEventPlaneC; // Use Event Planes from VZERO C

  Bool_t    fSubEtaGapTPC;    // bool to fill with eta gap
  Double_t  fEtaGap;          // Value of the eta gap
 
  Int_t     fDebugLevel; // Debug Level  

  // Cuts for HFE
  AliHFEcuts *fHFECuts;           // HFE cuts
  AliHFEpid  *fPID;               // PID cuts 
  AliHFEpidQAmanager *fPIDqa;     // QA Manager
  
  // Histos
  TH2D *fHistEV;               //! Number of events
  
  // A Event plane as function of phiepa, phiepb, phiepc, phiepd centrality 
  // a V0A, b V0C, c TPC, d V0
  THnSparseF *fEventPlane;     //! Event plane
  
  // Fbis Resolution as function of cosres, cosres, cosres, centrality for three subevents (V0)
  // a V0A, b V0C, c TPC
  THnSparseF *fCosResabc; //! Res
  
  // F Resolution as function of cosres, centrality for two subevents (TPC)
  THnSparseF *fCosRes; //! Res
  
  // G Maps delta phi as function of deltaphi, centrality, pt
  THnSparseF *fDeltaPhiMaps; //! Delta phi
  
  // H Maps cos phi : cos, centrality, pt
  THnSparseF *fCosPhiMaps;         //! Cos
    
  
  ClassDef(AliAnalysisTaskHFEFlowData, 1); // analysisclass
};

#endif
