#ifndef ALIANALYSISTASKCLQA_H
#define ALIANALYSISTASKCLQA_H

class TClonesArray;
class TString;
class TH1;
class TH1F;
class TH2F;
class TH3F;
class TNtuple;
class TNtupleD;
class TTree;

#include "AliAnalysisTaskEmcal.h"
#include "Cumulants.h"

class AliNtupHetInfo;
class AliNtupCumInfo;
class AliNtupZdcInfo;

class AliAnalysisTaskCLQA : public AliAnalysisTaskEmcal {
 public:
  AliAnalysisTaskCLQA();
  AliAnalysisTaskCLQA(const char *name);
  virtual ~AliAnalysisTaskCLQA();

  void                        SetCentCL1In(TH1F *h)                 { fCentCL1In       = h; }
  void                        SetCentV0AIn(TH1F *h)                 { fCentV0AIn       = h; }
  void                        SetCumParams(Double_t Mmin, Double_t ptmin, Double_t ptmax, Double_t etamin, Double_t etamax);
  void                        SetDoCumulants(Bool_t b, Bool_t bn=0) { fDoCumulants     = b; fDoCumNtuple = bn; }
  void                        SetDoMuonTracking(Bool_t b)           { fDoMuonTracking  = b; }
  void                        SetDoTrackProp(Bool_t b)              { fDoProp          = b; }
  void                        SetDoTracking(Bool_t b)               { fDoTracking      = b; }
  void                        SetDoVertexCut(Bool_t b)              { fDoVertexCut     = b; }
  void                        SetHetParams(Double_t Etmin);
  void                        SetDoHet(Bool_t b)                    { fDoHet           = b; }
  void                        SetQCEtaGap(Double_t e)               { fQC4EG           = e; }

  void                        UserCreateOutputObjects();

 protected:
  Bool_t                      FillHistograms();
  Bool_t                      RetrieveEventObjects();
  Bool_t                      Run();
  void                        RunCumulants();
  void                        RunCumulants(Double_t Mmin, Double_t ptmin, Double_t ptmax, Double_t etamin, Double_t etamax);
  void                        RunHet(Double_t Etmin);

  Bool_t                      fDoVertexCut;      // if true then use special (pileup) vertex checks
  Bool_t                      fDoTracking;       // if true run tracking analysis
  Bool_t                      fDoMuonTracking;   // if true run muon tracking analysis
  Bool_t                      fDoCumulants;      // if true run cumulant analysis
  Bool_t                      fDoCumNtuple;      // if true save cumulant ntuple
  Bool_t                      fDoProp;           // if true do track propagation
  Bool_t                      fDoHet;            // if true run het analysis
  Double_t                    fCumPtMin;         // minimum pt for cumulants
  Double_t                    fCumPtMax;         // maximum pt for cumulants
  Double_t                    fCumEtaMin;        // minimum eta for cumulants
  Double_t                    fCumEtaMax;        // maximum eta for cumulants
  Double_t                    fCumMmin;          // minimum number of tracks for cumulants 
  Int_t                       fCumMbins;         // number of bins for M
  Double_t                    fQC4EG;            // value for etagap (+-fQC4EG)
  Double_t                    fHetEtmin;         // minimum et cut for het
  TH1F                       *fCentCL1In;        // input for MC based CL1 centrality
  TH1F                       *fCentV0AIn;        // input for MC based V0A centrality
  TTree                      *fNtupCum;          //!ntuple for cumulant analysis
  AliNtupCumInfo             *fNtupCumInfo;      //!object holding cumulant results
  AliNtupZdcInfo             *fNtupZdcInfo;      //!object holding zdc info
  TTree                      *fNtupHet;          //!ntuple for het analysis
  AliNtupHetInfo             *fNtupHetInfo;      //!object holding het info
  TH1                        *fHists[1000];      //!pointers to histograms
  Cumulants                  *fCum;              //!pointer to cumulant class
 private:
  Double_t                    DeltaPhi(Double_t phia, Double_t phib,
                                       Double_t rangeMin = -TMath::Pi()/2, 
                                       Double_t rangeMax = 3*TMath::Pi()/2) const;
  AliAnalysisTaskCLQA(const AliAnalysisTaskCLQA&);            // not implemented
  AliAnalysisTaskCLQA &operator=(const AliAnalysisTaskCLQA&); // not implemented

  ClassDef(AliAnalysisTaskCLQA, 9) // Constantin's Task
};

class AliNtupHetInfo {
 public:
  AliNtupHetInfo() : fTrig(0), fRun(0), fVz(0), fIsFEC(0), fIsVSel(0), fIsP(0)
                     {;}
  virtual ~AliNtupHetInfo() {;}

 public:
  UInt_t        fTrig;         // trigger bits
  Int_t         fRun;          // run number
  Double32_t    fVz;           //[-32,32,8] vertex z
  Bool_t        fIsFEC;        // is first event in chunk
  Bool_t        fIsVSel;       // is vertex selected
  Bool_t        fIsP;          // is SPD pileup

  ClassDef(AliNtupHetInfo,3) // High energy cluster info
};

class AliNtupCumInfo {
 public:
  AliNtupCumInfo() : fTrig(0), fRun(0), fVz(0), fIsFEC(0), fIsVSel(0), fIsP(0),
                     fMall(0), fMall2(0), fPtMaxall(0), fMPtall(0), 
                     fMPt2all(0), fMPtall2(0), fTSall(0),
                     fM(0), fQ2abs(0), fQ4abs(0), fQ42re(0), fCos2phi(0), fSin2phi(0),
                     fPtMax(0), fMPt(0), fMPt2(0), fTS(0), fMV0M(0), 
                     fCl1(0), fV0M(0), fV0MEq(0), fV0A(0), fV0AEq(0), fZNA(0) {;}
  virtual ~AliNtupCumInfo() {;}

 public:
  UInt_t        fTrig;         // trigger selection
  Int_t         fRun;          // run number
  Double_t      fVz;           // vertex z
  Bool_t        fIsFEC;        // is first event in chunk
  Bool_t        fIsVSel;       // is vertex selected
  Bool_t        fIsP;          // is SPD pileup
  Int_t         fMall;         // multiplicity (tracks in eta range)
  Int_t         fMall2;        // multiplicity (tracks above 1 GeV/c in eta range)
  Double32_t    fPtMaxall;     //[0,0,16] maximum pT
  Double32_t    fMPtall;       //[0,0,16] mean pT
  Double32_t    fMPt2all;      //[0,0,16] mean pT2
  Double32_t    fMPtall2;      //[0,0,16] mean pT truncated above 1 GeV/c
  Double32_t    fTSall;        //[0,0,16] transverse sphericity
  Int_t         fM;            // multiplicity (tracks in pT range)
  Double32_t    fQ2abs;        // Q2 absolute
  Double32_t    fQ4abs;        // Q4 absolute
  Double32_t    fQ42re;        // Re(Q2Q*Q*)
  Double32_t    fCos2phi;      // Cos(2phi)
  Double32_t    fSin2phi;      // Sin(2phi)
  Double32_t    fPtMax;        //[0,0,16] maximum pT
  Double32_t    fMPt;          //[0,0,16] mean pT
  Double32_t    fMPt2;         //[0,0,16] mean pT2
  Double32_t    fTS;           //[0,0,16] transverse sphericity
  Double32_t    fMV0M;         // V0M amplitude
  Double32_t    fCl1;          //[0,0,16] class CL1
  Double32_t    fV0M;          //[0,0,16] class V0M
  Double32_t    fV0MEq;        //[0,0,16] class V0M Eq
  Double32_t    fV0A;          //[0,0,16] class V0A
  Double32_t    fV0AEq;        //[0,0,16] class V0A Eq
  Double32_t    fZNA;          //[0,0,16] class ZNA

  ClassDef(AliNtupCumInfo,3) // Cumulant storage class
};

class AliNtupZdcInfo {
 public:
   AliNtupZdcInfo() : fZna0(0), fZna1(0), fZna2(0), fZna3(0), fZna4(0) {;}
   virtual ~AliNtupZdcInfo() {;}

 public:
  Double32_t    fZna0;         // ZNA energy 0 
  Double32_t    fZna1;         // ZNA energy 1 
  Double32_t    fZna2;         // ZNA energy 2 
  Double32_t    fZna3;         // ZNA energy 3 
  Double32_t    fZna4;         // ZNA energy 4 

  ClassDef(AliNtupZdcInfo,1) // ZDC storage class
};
#endif
