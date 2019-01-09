/**************************************************************************
 * Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
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

//=========================================================================//
//             AliEbyE Analysis for Net-Particle Higher Moment study       //
//                           Nirbhay K. Behera                             //
//                           nbehera@cern.ch                               //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                                                                         //
//                        (Last Modified 2018/08/27)                       //
//                 Dealing with Wide pT Window Modified to ESDs            //
//Some parts of the code are taken from J. Thaeder/ M. Weber NetParticle   //
//analysis task.                                                           //
//=========================================================================//
//ver:2018/08/27 : tested only for proton

#ifndef AliEbyEPidEfficiencyContamination_H
#define AliEbyEPidEfficiencyContamination_H

class TList;
class TH2F;
class TH3F;
class TH1D;
class THnSparse;

class AliVEventHandler;
class AliMCEventHandler;
class AliVEvent;
class AliESDtrackCuts;
class AliMCEvent;
class AliStack;
class AliEventCuts;
class AliVParticle;
class AliVTrack;
class TClonesArray;
class AliPID;
class AliPIDResponse;


class AliEbyEPidEfficiencyContamination: public AliAnalysisTaskSE {
public:
AliEbyEPidEfficiencyContamination();
AliEbyEPidEfficiencyContamination( const char *name );
  virtual ~AliEbyEPidEfficiencyContamination();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetRunPeriod( TString runperiod) { fRun = runperiod; }
  void SetIsAOD(Bool_t IsAOD) {fIsAOD = IsAOD;}
  void SetAODtrackCutBit(Int_t bit) {fAODtrackCutBit = bit; }
  void SetIsMC(Bool_t Ismc) {fIsMC = Ismc;}
  
  void SetAnalysisCutObject(AliESDtrackCuts *const trackCuts) {
    fESDtrackCuts = trackCuts;}
  void SetDca(Double_t dcaxy,Double_t dcaz) { fDcaXy = dcaxy; fDcaZ = dcaz;}
  void SetVertexDiamond(Double_t vx, Double_t vy, Double_t vz) {
    fVxMax = vx;fVyMax = vy; fVzMax = vz;
  }
  
  void SetIsRapidityCut(Bool_t IsRapCut){ fIsRapCut = IsRapCut; }
  void SetUseTotalMomentumCut( Bool_t IsTotalMom){ fTotP = IsTotalMom; }
  void SetKinematicsCuts(Double_t ptl, Double_t pth, Double_t eta) {
    fPtMin = ptl; fPtMax = pth; fEtaMin = -eta; fEtaMax = eta; 
  }

  void SetNumberOfPtBins(Int_t nPtBins){ fNptBins = nPtBins;}
  
  void SetTPCTrackQualityCuts( Int_t NcrossRows, Double_t Chi2NDF){
    fNcrossRows = NcrossRows; fChi2NDF = Chi2NDF;
    
  }
  
  void SetCentralityEstimator(const char* centralityEstimator) { fCentralityEstimator = centralityEstimator;}
  
  void SetPidType(Int_t i);
  void SetPidStrategy(Int_t i)                       {fPidStrategy         = i;}
  void SetNSigmaMaxITS(Float_t f)                    {fNSigmaMaxITS        = f;}
  void SetNSigmaMaxTPC(Float_t f)                    {fNSigmaMaxTPC        = f;}
  void SetNSigmaMaxTOF(Float_t f)                    {fNSigmaMaxTOF        = f;}
  
 private:
  TList               *fThnList;       //
  AliVEventHandler    *fInputHandler;  // for Event handler
  AliMCEventHandler   *fMCEventHandler;         // for MC event handler---
  AliVEvent 	      *fVevent;        // V event
  TClonesArray        *fArrayMC;       // AOD MC stack
  AliESDtrackCuts     *fESDtrackCuts;  // ESD Track Cuts
  AliMCEvent          *fMCEvent;       // Current MC Event
  AliStack            *fMCStack;       // Stak tree
  AliEventCuts        *fEventCuts;      //using for pileup check
  TString             fRun;            //Run  production name 
  TString             fCentralityEstimator;   // "V0M","TRK","CL1"
 
  Int_t        fAODtrackCutBit; //
    
  Double_t   fVxMax;                        // X vertex  Range
  Double_t   fVyMax;                        // Y vertex Range
  Double_t   fVzMax;                        // Z vertex Range
  Double_t   fPtMin;                        // pT Minimum
  Double_t   fPtMax;                        // Pt Maximum
  Double_t   fEtaMin;                       // Eta Minimum
  Double_t   fEtaMax;                       // Eta Maximum
  Int_t      fNptBins;                     //pt bin numbers  
  Double_t   fDcaXy;                        // DCA Xy
  Double_t   fDcaZ;                         // DCA Z
  Int_t      fNcrossRows;                   //TPC nCross rows
  Double_t   fChi2NDF;                      //TPC track fit chi2/NDF 
  Bool_t     fIsMC;                         // Is MC event - Auto set by Add Task
  Bool_t     fIsAOD;                        // analysis mode: 0 = ESDs  | 1 = AODs
  Bool_t     fIsRapCut;                     // Use rapidity cut 1= yes, 0= no
  Bool_t     fTotP;                        // Swith to use total momentum cut

  Int_t   fNTracks;            // Number of Tracks of Current Events
  Float_t fCentrality;         //
 
  TH1D         *fEventCounter;  //!
  TH1F         *fHistCent;      //!

  TH2F  *fHitCentRec[2];
  TH2F  *fHitCentGen[2];
  

  TH3F *fCentPtEtaPhiThnGen[2];           //!
  TH3F *fCentPtEtaPhiThnRec[2];           //!
  TH3F *fCentPtEtaPhiThnRecPrim[2];       //!
  TH3F *fCentPtEtaPhiThnSec[2];           //!
  TH3F *fCentPtEtaPhiThnMat[2];           //!
  TH3F *fCentPtEtaPhiThnMisId[2];         //!
  
  AliPIDResponse   *fPIDResponse;              //! Ptr to PID response Object
  AliPIDCombined   *fPIDCombined;              //

  Int_t            fPidType;                  // 0=charge, 1=pion, 2=kaon, 3=proton
  Int_t            fMcPid;                    // MC PDG PID
  Int_t            fPidStrategy;              // 0: default, 1: ITS+TPC, 2: TPC+TOF
  Float_t          fNSigmaMaxITS;             //  N Sigma for ITS PID
  Float_t          fNSigmaMaxTPC;             //  N Sigma for TPC PID
  Float_t          fNSigmaMaxTOF;             //  N Sigma for TOF PID

  AliPID::EParticleType fParticleSpecies;     //  Particle species on basis of AliPID


  THnSparse *fPtBinNplusNminusCh;
  THnSparse *fPtBinNplusNminusChTruth;

  enum EbyEPIDType_t { kNSigmaTPC = 0, kNSigmaTOF, kNSigmaTPCTOF, kBayes };
  enum EbyEDetectorType_t { kITS = 0, kTPC,  kTOF,  kNDetectors };
  enum EbyECharge_t { kChPos = 0, kChNeg, kNCharge};

  /****** Some Private Functions ******/

  Bool_t   IsPidPassed(AliVTrack * track);
  Double_t TOFBetaCalc(AliVTrack *track) const;
  Double_t GetMass(AliPID::EParticleType id) const;
  
  Bool_t   AcceptTrackLMC(AliVParticle *particle) const; //
  Bool_t   AcceptTrackL(AliVTrack *track) const; //
  Double_t GetRapidity(Float_t pt, Float_t pz) const; //to get rapidity
  Int_t    GetPtBin(Double_t pt); //to get the bin number of the pt
  Int_t    GetEtaBin(Float_t eta);   //To get the eta bin
  void     LocalPost(); //
  void     CreateEffCont();

  //________________________________
  AliEbyEPidEfficiencyContamination(const AliEbyEPidEfficiencyContamination&);
  AliEbyEPidEfficiencyContamination& operator = (const AliEbyEPidEfficiencyContamination&);
  ClassDef(AliEbyEPidEfficiencyContamination, 6);
};

#endif

 
