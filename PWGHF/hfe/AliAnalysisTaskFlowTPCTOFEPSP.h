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
#ifndef ALIANALYSISTASKFLOWTPCTOFEPSP_H
#define ALIANALYSISTASKFLOWTPCTOFEPSP_H




#include <AliAnalysisTaskSE.h>

class TList;
class AliVTrack;
class AliVEvent;
class AliESDtrack;
class AliESDEvent;
class AliMCEvent;
class AliFlowTrackCuts;
class AliFlowCandidateTrack;
class AliHFEcuts;
class AliHFEpid;
class TH1D;
class TH2D;
class TF1;
class TProfile;
class TProfile2D;
class THnSparse;
class AliHFEpidQAmanager;
class AliFlowEvent;
class AliESDtrackCuts;
class AliHFEVZEROEventPlane;
class TArrayI;
class AliAODMCHeader;
class TClonesArray;
class AliHFENonPhotonicElectron;
class TTreeSRedirector;

class AliAnalysisTaskFlowTPCTOFEPSP: public AliAnalysisTaskSE {
public:

  typedef enum{
    kElectronfromconversion = 0,
    kElectronfromconversionboth = 1,
    kElectronfrompi0 = 2,
    kElectronfrompi0both = 3,
    kElectronfrometa = 4,
    kElectronfrometaboth = 5,
    kElectronfromC = 6,
    kElectronfromB = 7,
    kElectronfromother = 8,
    kNoElectron = 9
  } FlowSource_t;
  
  typedef enum{
    kS = 0,
    kOp = 1
  } FlowSign_t;




  AliAnalysisTaskFlowTPCTOFEPSP();
  AliAnalysisTaskFlowTPCTOFEPSP(const char *name);
  AliAnalysisTaskFlowTPCTOFEPSP(const AliAnalysisTaskFlowTPCTOFEPSP &ref);
  AliAnalysisTaskFlowTPCTOFEPSP& operator=(const AliAnalysisTaskFlowTPCTOFEPSP &ref);
  virtual void Copy(TObject &o) const;
  virtual ~AliAnalysisTaskFlowTPCTOFEPSP();
  
  virtual void  UserExec(Option_t */*option*/);
  virtual void  UserCreateOutputObjects();

  void SetAODAnalysis(Bool_t aodAnalysis)   { fAODAnalysis = aodAnalysis; };
  void SetUseFilterAOD(Bool_t useFilterAOD) { fUseFilterAOD = useFilterAOD; }
  void SetApplyCut(Bool_t applyCut)         { fApplyCut = applyCut; }
  void SetFilter(ULong_t filter)            { fFilter = filter; }
  
  AliHFEpid *GetPID() const { return fPID; }
  AliHFEpid *GetPIDTOFOnly() const { return fPIDTOFOnly; }
  AliHFEpidQAmanager *GetPIDQAManager() const { return fPIDqa; }
  AliHFEpid *GetPIDBackground() const { return fPIDBackground; }
  AliHFEpidQAmanager *GetPIDBackgroundQAManager() const { return fPIDBackgroundqa; }
  AliHFENonPhotonicElectron *GetHFEBackgroundSubtraction() const { return fBackgroundSubtraction; }


  void SetContamination(TF1 * const function,Int_t k) { fContamination[k] = function; };
  void SetV2Contamination(TF1 * const function,Int_t k) { fv2contamination[k] = function; };
  void SetHFECuts(AliHFEcuts * const cuts) { fHFECuts = cuts; };
  void SetRejectKinkMother(Bool_t rejectKinkMother = kFALSE) { fRejectKinkMother = rejectKinkMother; };
  void SetHFEBackgroundSubtraction(AliHFENonPhotonicElectron * const backgroundSubtraction) { fBackgroundSubtraction = backgroundSubtraction; };
  void SetHFEBackgroundCuts(AliESDtrackCuts * const cuts) { fHFEBackgroundCuts = cuts; };
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
  void SetUseSP(Bool_t useSP) { fSP = useSP;}
  void SetMCPID(Bool_t mcPID) { fMCPID = mcPID;};
  void SetNoPID(Bool_t noPID) { fNoPID = noPID;};

  void SetDebugLevel(Int_t debugLevel) { fDebugLevel = debugLevel;};
  void SetMonitorEventPlane(Bool_t monitorEventPlane) { fMonitorEventPlane = monitorEventPlane;};
  void SetMonitorContamination(Bool_t monitorContamination) { fMonitorContamination = monitorContamination;};
  void SetMonitorPhotonic(Bool_t monitorPhotonic) { fMonitorPhotonic = monitorPhotonic;};
  void SetMonitorWithoutPID(Bool_t monitorWithoutPID) { fMonitorWithoutPID = monitorWithoutPID;};
  void SetMonitorTrackCuts(Bool_t monitorTrackCuts) { fMonitorTrackCuts = monitorTrackCuts;};
  void SetMonitorQCumulant(Bool_t monitorQCumulant) { fMonitorQCumulant = monitorQCumulant;};

  Int_t GetNbBinsCentralityQCumulant() const { return  fNbBinsCentralityQCumulant; };
  Double_t GetBinCentralityLess(Int_t k) const { return fBinCentralityLess[k]; };
  
  AliFlowCandidateTrack *MakeTrack( Double_t mass, Double_t pt, Double_t phi, Double_t eta) ;
  Double_t GetPhiAfterAddV2(Double_t phi,Double_t reactionPlaneAngle) const;

  void  SetMaxInvmass(Double_t maxInvmass) { fMaxInvmass = maxInvmass; };
  void  SetMaxopening3D(Double_t maxOpening3D) { fMaxopening3D = maxOpening3D; };
  void  SetMaxopeningtheta(Double_t maxOpeningtheta) { fMaxopeningtheta = maxOpeningtheta; };
  void  SetMaxopeningphi(Double_t maxOpeningphi) { fMaxopeningphi = maxOpeningphi; };
  void  SetAlgorithmMA(Bool_t algorithmMA) { fAlgorithmMA = algorithmMA; };
  void  SetMassConstraint(Bool_t massConstraint) { fSetMassConstraint = massConstraint; };
  void  SetPileUpCut(Bool_t cut=kTRUE) { fPileUpCut=cut; }

  Int_t    LookAtNonHFE(Int_t iTrack1, AliVTrack *track1, AliVEvent *fESD, AliMCEvent *mcEvent,Int_t binct,Double_t deltaphi,Int_t source,Int_t indexmother);
  
private:
  TList     *fListHist;         //! TH list
  Bool_t    fAODAnalysis;       // AOD analysis
  Bool_t    fUseFilterAOD;     // Use the preselected AOD track
  Bool_t    fApplyCut;       // Apply the analysis cut for AOD tracks
  ULong_t   fFilter;             // reconstruction AOD status flags 
  AliAODMCHeader *fAODMCHeader;         // ! MC info AOD
  TClonesArray *fAODArrayMCInfo;        // ! MC info particle AOD
  AliHFENonPhotonicElectron *fBackgroundSubtraction; // Background subtraction
  
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
  Bool_t    fSP;        // calculate using scalar product method (instead of event plane method)

  Bool_t    fMCPID; // MC PID for electrons
  Bool_t    fNoPID; // No PID for checks

  Double_t  fChi2OverNDFCut;   // Limit chi2
  Double_t  fMaxdca;           // Limit dca
  Double_t  fMaxopeningtheta;  // Limit opening angle in theta
  Double_t  fMaxopeningphi;    // Limit opening angle in phi
  Double_t  fMaxopening3D;     // Limit opening 3D
  Double_t  fMaxInvmass;       // Limit invariant mass
  Bool_t    fSetMassConstraint; // Set mass constraint
  

  Int_t     fDebugLevel; // Debug Level  
  Bool_t    fMonitorEventPlane; // Monitor event plane
  Bool_t    fMonitorContamination; // Monitor contamination
  Bool_t    fMonitorPhotonic;// Monitor photonic
  Bool_t    fMonitorWithoutPID;// Monitor without PID
  Bool_t    fMonitorTrackCuts;// Monitor track cuts
  Bool_t    fMonitorQCumulant;// Monitor Q cumulant
  
  // Cuts for FLOW PWG2
  AliFlowTrackCuts* fcutsRP;  //! Reference particle cut
  AliFlowTrackCuts* fcutsPOI; //! Particle Of Interest cut
  
  // Cuts for HFE
  AliHFEcuts *fHFECuts;           // HFE cuts
  Bool_t fRejectKinkMother;       // Reject Kink Mother 
  AliHFEpid  *fPID;               // PID cuts 
  AliHFEpid  *fPIDTOFOnly;        // PID cuts TOF only
  AliHFEpidQAmanager *fPIDqa;     // QA Manager
  AliFlowEvent *fflowEvent;       //! Flow event 

  // Hadron Contamination
  TF1 *fContamination[11];        // Parametrization of the contamination (0-5,5-10,10-20,20-30,30-40,40-50,50-60,60-70,70-80,80-90,90-100)
  TF1 *fv2contamination[11];      // Parametrization of the v2 of charged pions (0-5,5-10,10-20,20-30,30-40,40-50,50-60,60-70,70-80,80-90,90-100)

  // Cuts for background study
  AliESDtrackCuts *fHFEBackgroundCuts;    // HFE background cuts
  AliHFEpid  *fPIDBackground;             // PID background cuts 
  AliHFEpidQAmanager *fPIDBackgroundqa;   // QA Manager Background  
  Bool_t fAlgorithmMA;                    // algorithm MA

  // List of tracks
  TArrayI *fArraytrack;                    //! list of tracks
  Int_t fCounterPoolBackground;            // number of tracks

  // VZERO Event plane after calibration 2010
  AliHFEVZEROEventPlane *fHFEVZEROEventPlane; // VZERO event plane calibrated
  
  // Histos
  TH2D *fHistEV;               //! Number of events
  THnSparseF *fHistPileUp;     //! Pile up histogram
  Bool_t fPileUpCut;

  // A Event plane as function of phiepa, phiepb, phiepc, phiepd centrality 
  // a V0A, b V0C, c TPC,
  THnSparseF *fEventPlane;     //! Event plane
  
  // B Event Plane after subtraction as function of phiep, centrality 
  THnSparseF *fEventPlaneaftersubtraction; //! Event plane

  // Contamination
  THnSparseF *fFractionContamination;    //! Fraction of contamination as function of pt
  TProfile2D *fContaminationv2;          //! v2 of contamination

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
  
  // Debuging Cuts step by step all centrality together: pt, step (6)
  THnSparseF *fTrackingCuts; //! Tracking Cuts

  // Before PID cut
  // G Maps delta phi as function of deltaphi, centrality, pt
  THnSparseF *fDeltaPhiMapsBeforePID; //! Delta phi
  // H Maps cos phi : cos, centrality, pt
  THnSparseF *fCosPhiMapsBeforePID; //! Cos

  // G Maps delta phi as function of deltaphi, centrality, pt
  THnSparseF *fDeltaPhiMaps; //! Delta phi
  THnSparseF *fDeltaPhiMapsContamination; //! Delta phi for contamination substraction
  // H Maps cos phi : cos, centrality, pt
  THnSparseF *fCosPhiMaps;         //! Cos
  TProfile2D *fProfileCosPhiMaps;  //! Profile Cos

  // Background study: not statistic but tagged 
  THnSparseF *fDeltaPhiMapsTaggedPhotonic; //! Delta phi
  //THnSparseF *fCosPhiMapsTaggedPhotonic; //! Cos
  THnSparseF *fDeltaPhiMapsTaggedNonPhotonic; //! Delta phi
  //THnSparseF *fCosPhiMapsTaggedNonPhotonic; //! Cos
  THnSparseF *fDeltaPhiMapsTaggedPhotonicLS; //! Delta phi
  //THnSparseF *fCosPhiMapsTaggedPhotonicLS; //! Cos

  // Background study: centrality, pt, source
  THnSparseF *fMCSourceDeltaPhiMaps; //! Source MC
  // Background study: deltaphi, centrality, pt, minv, source
  THnSparseF *fOppSignDeltaPhiMaps;  //! Delta phi
  THnSparseF *fSameSignDeltaPhiMaps; //! Delta phi
  // Background study: angle, centrality, source
  THnSparseF *fOppSignAngle;         // ! Opening Angles
  THnSparseF *fSameSignAngle;        // ! Opening Angles

  TTreeSRedirector  *fDebugStreamer;               //!Debug streamer

 Int_t FindMother(Int_t tr, AliMCEvent *mcEvent, Int_t &indexmother);
  Int_t CheckPdg(Int_t tr, AliMCEvent* mcEvent);
  Int_t IsMotherGamma(Int_t tr, AliMCEvent* mcEvent);
  Int_t IsMotherPi0(Int_t tr, AliMCEvent* mcEvent);
  Int_t IsMotherC(Int_t tr, AliMCEvent* mcEvent);
  Int_t IsMotherB(Int_t tr, AliMCEvent* mcEvent);
  Int_t IsMotherEta(Int_t tr, AliMCEvent* mcEvent);
    
  
  ClassDef(AliAnalysisTaskFlowTPCTOFEPSP, 1); // analysisclass
};

#endif
