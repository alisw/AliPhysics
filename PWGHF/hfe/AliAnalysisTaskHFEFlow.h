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
#ifndef ALIANALYSISTASKHFEFLOW_H
#define ALIANALYSISTASKHFEFLOW_H




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

class AliAnalysisTaskHFEFlow: public AliAnalysisTaskSE {
  
public:
  AliAnalysisTaskHFEFlow();
  AliAnalysisTaskHFEFlow(const char *name);
  AliAnalysisTaskHFEFlow(const AliAnalysisTaskHFEFlow &ref);
  AliAnalysisTaskHFEFlow& operator=(const AliAnalysisTaskHFEFlow &ref);
  virtual void Copy(TObject &o) const;
  virtual ~AliAnalysisTaskHFEFlow();
  
  virtual void  UserExec(Option_t */*option*/);
  virtual void  UserCreateOutputObjects();

  AliHFEpid *GetPID() const { return fPID; }
  AliHFEpidQAmanager *GetPIDQAManager() const { return fPIDqa; }

  void SetHFECuts(AliHFEcuts * const cuts) { fHFECuts = cuts; };
  void SetSubEtaGapTPC(Bool_t  subEtaGapTPC) { fSubEtaGapTPC = subEtaGapTPC; };
  void SetEtaGap(Double_t  etaGap) { fEtaGap = etaGap; };
  void SetVZEROEventPlane(Bool_t vzeroEventPlane) { fVZEROEventPlane = vzeroEventPlane; };
  void SetVZEROEventPlaneA(Bool_t vzeroEventPlaneA) { fVZEROEventPlaneA = vzeroEventPlaneA; };
  void SetVZEROEventPlaneC(Bool_t vzeroEventPlaneC) { fVZEROEventPlaneC = vzeroEventPlaneC; };
  void SetHFEVZEROEventPlane(AliHFEVZEROEventPlane *hfeVZEROEventPlane) { fHFEVZEROEventPlane = hfeVZEROEventPlane; };

  void SetNbBinsCentralityQCumulant(Int_t nbBinsCentralityQCumulant) { fNbBinsCentralityQCumulant = nbBinsCentralityQCumulant; };
  void SetBinCentralityLess(Int_t k, Float_t value)  { fBinCentralityLess[k] = value; };
  void SetNbBinsPtQCumulant(Int_t nbBinsPtQCumulant) { fNbBinsPtQCumulant = nbBinsPtQCumulant; };
  void SetMinPtQCumulant(Double_t minPtQCumulant) { fMinPtQCumulant = minPtQCumulant; };
  void SetMaxPtQCumulant(Double_t maxPtQCumulant) { fMaxPtQCumulant = maxPtQCumulant; };

  void SetAfterBurnerOn(Bool_t afterBurnerOn)     { fAfterBurnerOn = afterBurnerOn; };
  void SetNonFlowNumberOfTrackClones(Int_t nonFlowNumberOfTrackClones) { fNonFlowNumberOfTrackClones = nonFlowNumberOfTrackClones; };
  void SetV1V2V3V4V5(Double_t v1,Double_t v2,Double_t v3,Double_t v4,Double_t v5) {fV1 = v1; fV2 = v2; fV3 = v3; fV4 = v4; fV5 = v5; };
  void SetMaxNumberOfIterations(Int_t maxNumberOfIterations) { fMaxNumberOfIterations = maxNumberOfIterations; };
  void SetPrecisionPhi(Double_t precisionPhi) { fPrecisionPhi = precisionPhi;};
  void SetUseMCReactionPlane(Bool_t useMCReactionPlane) { fUseMCReactionPlane = useMCReactionPlane;};
  void SetMCPID(Bool_t mcPID) { fMCPID = mcPID;};
  void SetNoPID(Bool_t noPID) { fNoPID = noPID;};

  void SetDebugLevel(Int_t debugLevel) { fDebugLevel = debugLevel;};

  Int_t GetNbBinsCentralityQCumulant() const { return  fNbBinsCentralityQCumulant; };
  Double_t GetBinCentralityLess(Int_t k) const { return fBinCentralityLess[k]; };
  
  AliFlowCandidateTrack *MakeTrack( Double_t mass, Double_t pt, Double_t phi, Double_t eta) ;
  Double_t GetPhiAfterAddV2(Double_t phi,Double_t reactionPlaneAngle) const;
  
private:
  TList       *fListHist;       //! TH list

  Bool_t    fVZEROEventPlane;  // Use Event Planes from VZERO
  Bool_t    fVZEROEventPlaneA; // Use Event Planes from VZERO A
  Bool_t    fVZEROEventPlaneC; // Use Event Planes from VZERO C

  Bool_t    fSubEtaGapTPC;    // bool to fill with eta gap
  Double_t  fEtaGap;          // Value of the eta gap

  Int_t     fNbBinsCentralityQCumulant;  // Number of Bins Q Cumulant
  Double_t  fBinCentralityLess[10];      // Centrality Bin lower value
  Int_t     fNbBinsPtQCumulant;          // Nbbinspt QCumulant method
  Double_t  fMinPtQCumulant;             // Min pt QCumulant method
  Double_t  fMaxPtQCumulant;             // Max pt QCumulant method
  Bool_t    fAfterBurnerOn;              // Add flow to all tracks
  Int_t     fNonFlowNumberOfTrackClones; // number of times to clone the particles (nonflow) 
  Double_t  fV1;        // Add Flow. Must be in range [0,0.5].
  Double_t  fV2;        // Add Flow. Must be in range [0,0.5].
  Double_t  fV3;        // Add Flow. Must be in range [0,0.5].
  Double_t  fV4;        // Add Flow. Must be in range [0,0.5].
  Double_t  fV5;        // Add Flow. Must be in range [0,0.5].
  Int_t     fMaxNumberOfIterations; // Max number of iteration for adding v2
  Double_t  fPrecisionPhi;  // precision phi for adding v2
  Bool_t    fUseMCReactionPlane; // use MC reaction plane

  Bool_t    fMCPID; // MC PID for electrons
  Bool_t    fNoPID; // No PID for checks

  Int_t     fDebugLevel; // Debug Level  

  // Cuts for FLOW PWG2
  AliFlowTrackCuts* fcutsRP;  // Reference particle cut
  AliFlowTrackCuts* fcutsPOI; // Particle Of Interest cut
  
  // Cuts for HFE
  AliHFEcuts *fHFECuts;           // HFE cuts
  AliHFEpid  *fPID;               // PID cuts 
  AliHFEpidQAmanager *fPIDqa;     // QA Manager
  AliFlowEvent *fflowEvent;       //! Flow event   

  // VZERO Event plane after calibration 2010
  AliHFEVZEROEventPlane *fHFEVZEROEventPlane; // VZERO event plane calibrated
  
  // Histos
  TH2D *fHistEV;               //! Number of events
  
  // A Event plane as function of phiepa, phiepb, phiepc, phiepd centrality 
  // a V0A, b V0C, c TPC, d V0
  THnSparseF *fEventPlane;     //! Event plane
  
  // B Event Plane after subtraction as function of phiep, centrality 
  THnSparseF *fEventPlaneaftersubtraction; //! Event plane

  // Monitoring Event plane: cos2phi, sin2phi, centrality
  THnSparseF *fCosSin2phiep;        //! Cos(2phi), Sin(2phi)
  
  // E Monitoring Event plane after subtraction of the track: cos, centrality, pt, eta
  THnSparseF *fCos2phie;  //! Monitoring
  THnSparseF *fSin2phie;  //! Monitoring
  THnSparseF *fCos2phiep;  //! Monitoring
  THnSparseF *fSin2phiep;  //! Monitoring
  THnSparseF *fSin2phiephiep;  //! Monitoring

  // Fbis Resolution as function of cosres, cosres, cosres, centrality for three subevents (V0)
  // a V0A, b V0C, c TPC
  THnSparseF *fCosResabc; //! Res
  THnSparseF *fSinResabc; //! Res
  TProfile   *fProfileCosResab; //! Profile Res_a_b
  TProfile   *fProfileCosResac; //! Profile Res_a_c
  TProfile   *fProfileCosResbc; //! Profile Res_b_c
  
  // F Resolution as function of cosres, centrality for two subevents (TPC)
  THnSparseF *fCosRes; //! Res
  THnSparseF *fSinRes; //! Res
  TProfile   *fProfileCosRes; //! Profile Res
  
  // G Maps delta phi as function of deltaphi, centrality, pt
  THnSparseF *fDeltaPhiMaps; //! Delta phi
  
  // H Maps cos phi : cos, centrality, pt
  THnSparseF *fCosPhiMaps;         //! Cos
  TProfile2D *fProfileCosPhiMaps;  //! Profile Cos
  
  
  ClassDef(AliAnalysisTaskHFEFlow, 1); // analysisclass
};

#endif
