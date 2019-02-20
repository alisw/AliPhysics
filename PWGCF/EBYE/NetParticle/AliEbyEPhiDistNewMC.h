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

//=========================================================================//
//             AliEbyE Analysis for Net-Particle study                     //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                            Surya Prakash Pathak                         //
//                       surya.prakash.pathak@cern.ch                      //
//                         (Last Modified 2019/01/23)                      //
//                                                                         //
//=========================================================================//


#ifndef AliEbyEPhiDistNewMC_H
#define AliEbyEPhiDistNewMC_H

class TList;
class TH2F;
class TH1D;
class THnSparse;

class AliMCEvent;
class AliStack;
class AliESDtrackCuts;
class AliVParticle;
class AliVTrack;
class TClonesArray;
class AliPID;
class AliPIDResponse;
class AliPIDCombined;
class AliAnalysisUtils;


class AliEbyEPhiDistNewMC: public AliAnalysisTaskSE {
public:
AliEbyEPhiDistNewMC();
AliEbyEPhiDistNewMC( const char *name );
  virtual ~AliEbyEPhiDistNewMC();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetRunPeriod( TString runperiod) { fRun = runperiod; }
  void SetIsAOD(Bool_t IsAOD) {fIsAOD = IsAOD;}
  void SetAODtrackCutBit(Int_t bit) {fAODtrackCutBit = bit; }
  void SetIsMC(Bool_t Ismc) {fIsMC = Ismc;}
  void RunQA(Bool_t IsQA) {fIsQA = IsQA;}
  
  void SetAnalysisCutObject(AliESDtrackCuts *const trackCuts) {
    fESDtrackCuts = trackCuts;}
  void SetDca(Double_t dcaxy,Double_t dcaz) { fDcaXy = dcaxy; fDcaZ = dcaz;}
  void SetVertexDiamond(Double_t vx, Double_t vy, Double_t vz) {
    fVxMax = vx;fVyMax = vy; fVzMax = vz;
  }
  
  void SetKinematicsCuts(Double_t ptl, Double_t pth, Double_t eta) {
    fPtMin = ptl; fPtMax = pth; fEtaMin = -eta; fEtaMax = eta; 
  }

  void SetNumberOfPtBins(Int_t nPtBins){ fNptBins = nPtBins;}
    
    void SetPhi(Double_t phil, Double_t phih) // i added
    {
        fPhiMin = phil ; fPhiMax = phih ;
    }

    void SetNumberOfPhiBins(Int_t nPhiBins){fNphiBins = nPhiBins;}
  
  void SetTPCTrackQualityCuts( Int_t NcrossRows, Double_t Chi2NDF){
    fNcrossRows = NcrossRows; fChi2NDF = Chi2NDF;
    
  }
  
  void SetCentralityEstimator(const char* centralityEstimator) { fCentralityEstimator = centralityEstimator;}
  void SetIsKMb() { fIsTrig = kTRUE; } 
  static const Int_t kTrack = 15000;
  
  void SetPidType(Int_t i);
  void SetPidStrategy(Int_t i)                       {fPidStrategy         = i;}
  void SetNSigmaMaxITS(Float_t f)                    {fNSigmaMaxITS        = f;}
  void SetNSigmaMaxTPC(Float_t f)                    {fNSigmaMaxTPC        = f;}
  void SetNSigmaMaxTOF(Float_t f)                    {fNSigmaMaxTOF        = f;}
  void SetNSigmaMaxTPClow(Float_t f)                 {fNSigmaMaxTPClow     = f;}
  void SetMinPtForTOFRequired(Float_t f)             {fMinPtForTOFRequired = f;}
  void SetMaxPtForTPClow(Float_t f)                  {fMaxPtForTPClow      = f;}
  
 private:
  TList           *fThnList;       //
  TClonesArray    *fArrayMC;       // AOD MC stack
  AliESDtrackCuts *fESDtrackCuts;  // ESD Track Cuts
  AliMCEvent      *fMCEvent;       // Current MC Event
  AliStack        *fMCStack;       // Stak tree
  AliAnalysisUtils *fanaUtils;     //using for pileup check
    AliEventCuts        *fEventCuts;      //!using for pileup check
  TString          fRun;           //Run  production name 
  TString          fCentralityEstimator;   // "V0M","TRK","CL1"
 
  Int_t        fAODtrackCutBit; //
    
  Double_t   fVxMax;                        // X vertex  Range
  Double_t   fVyMax;                        // Y vertex Range
  Double_t   fVzMax;                        // Z vertex Range
  Double_t   fPtMin;                        // pT Minimum
  Double_t   fPtMax;                        // Pt Maximum
  Double_t   fEtaMin;                       // Eta Minimum
  Double_t   fEtaMax;                       // Eta Maximum
  Int_t      fNptBins;                     //pt bin numbers
    Double_t   fPhiMin;                       //phi minimum
    Double_t   fPhiMax;                      //phi maximum
    Int_t      fNphiBins;                    //phi bin numbers

  Double_t   fDcaXy;                        // DCA Xy
  Double_t   fDcaZ;                         // DCA Z
  Int_t      fNcrossRows;                   //TPC nCross rows
  Double_t   fChi2NDF;                      //TPC track fit chi2/NDF 
  Bool_t     fIsMC;                         // Is MC event - Auto set by Add Task
  Bool_t     fIsAOD;                        // analysis mode: 0 = ESDs  | 1 = AODs
  Bool_t     fIsQA;                         // Check for QA
  Bool_t     fIsTrig;           //
  Bool_t     fIsThn;            //

  Int_t   fNTracks;            // Number of Tracks of Current Events
  Float_t fCentrality;         //
 
  TH1D         *fEventCounter;  //!
  TH1F         *fHistCent;      //!

  TH2F  *fHitCentRec[2];
  TH2F  *fHitCentGen[2];
  
  TH2F  *fHistERec[2][3];         //! Reconstructed Plus
  TH2F  *fHistERecPri[2][3];      //! Reconstructed Primary Plus
  TH2F  *fHistEGen[2][3];         //!
  
  TH2F  *fHistCSec[2][3];         //!
  TH2F  *fHistCMat[2][3];         //!
  TH2F  *fHistCMisId[2][3];       //!
  
  AliPIDResponse   *fPIDResponse;              //! Ptr to PID response Object
  AliPIDCombined   *fPIDCombined;              //

  Int_t            fPidType;                  // 0=charge, 1=pion, 2=kaon, 3=proton
  Int_t            fMcPid;                    // MC PDG PID
  Int_t            fPidStrategy;              // 0: default, 1: ITS+TPC, 2: TPC+TOF
  Float_t          fNSigmaMaxITS;             //  N Sigma for ITS PID
  Float_t          fNSigmaMaxTPC;             //  N Sigma for TPC PID
  Float_t          fNSigmaMaxTPClow;          //  N Sigma for TPC PID lower part
  Float_t          fNSigmaMaxTOF;             //  N Sigma for TOF PID
  Float_t          fMinPtForTOFRequired;      //  Min pt from where TOF is required
  Float_t          fMaxPtForTPClow;           //  Max pt until TPClow is used
  Double_t nPidRec[2];
  Double_t nPidRecP[2];
  Double_t nPidRecMid[2];
  Double_t nPidRecSec[2];
  Double_t nPidRecWD[2];
  Double_t nPidWoPID[2];

  AliPID::EParticleType fParticleSpecies;     //  Particle species on basis of AliPID


  TH2F *fHistTPC; //! 
  TH2F *fHistTOF; //!
  TH2F *fHistITS; //!

  TH2F *fHistTPCc; //!
  TH2F *fHistTOFc; //!
  TH2F *fHistITSc; //!

  TH2F *fHistTPCTOF;  //!
  TH2F *fHistTPCTOFc; //!
  

  TH2F *fHistNsTPC; //!
  TH2F *fHistNsTOF; //!
  TH2F *fHistNsITS; //!

  TH2F *fHistNsTPCc; //!
  TH2F *fHistNsTOFc; //!
  TH2F *fHistNsITSc; //!

  THnSparse *fPhiBinNplusNminusCh;
  THnSparse *fPhiBinNplusNminusChTruth;
  THnSparse *fTHnCentNplusNminusCh;

  enum EbyEPIDType_t { kNSigmaTPC = 0, kNSigmaTOF, kNSigmaTPCTOF, kBayes };
  enum EbyEDetectorType_t { kITS = 0, kTPC,  kTOF,  kNDetectors };
  enum EbyECharge_t { kChPos = 0, kChNeg, kNCharge};

  /****** Some Private Functions ******/

  Bool_t   IsPidPassed(AliVTrack * track);
  Double_t TOFBetaCalc(AliVTrack *track) const;
  Double_t GetMass(AliPID::EParticleType id) const;
  
  Bool_t   AcceptTrackLMC(AliVParticle *particle) const; //
  Bool_t   AcceptTrackL(AliVTrack *track) const; //
  Int_t    GetPtBin(Double_t pt); //to get the bin number of the pt
    Int_t    GetPhiBin(Double_t Phi);
  void     LocalPost(); //
  void     CreatePhiHistMC();

  //________________________________
  AliEbyEPhiDistNewMC(const AliEbyEPhiDistNewMC&);
  AliEbyEPhiDistNewMC& operator = (const AliEbyEPhiDistNewMC&);
  ClassDef(AliEbyEPhiDistNewMC, 2);
};

#endif

 
