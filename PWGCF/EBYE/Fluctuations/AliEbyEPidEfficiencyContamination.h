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

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//                                                                         //
//             AliEbyE Analysis  Net-Particle Higher Moment study          //
//              Author:   Deepika Rathee  || Satyajit Jena                 //
//                        drathee@cern.ch || sjena@cern.ch                 //
//                                Nirbhay K. Behera                        //              
//                          (Last modified: 2017/04/07)                    //
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#ifndef AliEbyEPidEfficiencyContamination_H
#define AliEbyEPidEfficiencyContamination_H

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


class AliEbyEPidEfficiencyContamination: public AliAnalysisTaskSE {
public:
AliEbyEPidEfficiencyContamination();
AliEbyEPidEfficiencyContamination( const char *name );
  virtual ~AliEbyEPidEfficiencyContamination();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

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
  
  void SetTPCTrackQualityCuts( Int_t NcrossRows, Double_t Chi2NDF){
    fNcrossRows = NcrossRows; fChi2NDF = Chi2NDF;
    
  }
  
  void SetCentralityEstimator(const char* centralityEstimator) { fCentralityEstimator = centralityEstimator;}
  void SetIsKMb() { fIsTrig = kTRUE; } 
  static const Int_t kTrack = 15000;
  
  void SetPidType(Int_t i);
  void SetNSigmaMaxITS(Float_t f)                    {fNSigmaMaxITS        = f;}
  void SetNSigmaMaxTPC(Float_t f)                    {fNSigmaMaxTPC        = f;}
  void SetNSigmaMaxTOF(Float_t f)                    {fNSigmaMaxTOF        = f;}
  void SetNSigmaMaxTPClow(Float_t f)                 {fNSigmaMaxTPClow     = f;}
  void SetMinPtForTOFRequired(Float_t f)             {fMinPtForTOFRequired = f;}
  void SetMaxPtForTPClow(Float_t f)                  {fMaxPtForTPClow      = f;}
  
 private:
  TList           *fThnList;       //!
  Double_t        *fPtArray;       //pt array
  TClonesArray    *fArrayMC;       // AOD MC stack
  AliESDtrackCuts *fESDtrackCuts;  // ESD Track Cuts
  AliMCEvent      *fMCEvent;       // Current MC Event
  AliStack        *fMCStack;       // Stak tree
  AliAnalysisUtils *fanaUtils;     //using for pileup check
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

  Int_t fPidType; //
  Int_t fMcPid; //
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

  THnSparse *fPtBinNplusNminusCh;
  THnSparse *fPtBinNplusNminusChTruth;
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
  void     LocalPost(); //
  void     CreateEffCont();

  //________________________________
  AliEbyEPidEfficiencyContamination(const AliEbyEPidEfficiencyContamination&);
  AliEbyEPidEfficiencyContamination& operator = (const AliEbyEPidEfficiencyContamination&);
  ClassDef(AliEbyEPidEfficiencyContamination, 1);
};

#endif

 
