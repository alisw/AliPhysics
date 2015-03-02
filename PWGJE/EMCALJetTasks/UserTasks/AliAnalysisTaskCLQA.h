#ifndef ALIANALYSISTASKCLQA_H
#define ALIANALYSISTASKCLQA_H

// $Id $

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

class AliNtupCumInfo;
class AliNtupZdcInfo;

class AliAnalysisTaskCLQA : public AliAnalysisTaskEmcal {
 public:
  AliAnalysisTaskCLQA();
  AliAnalysisTaskCLQA(const char *name);
  virtual ~AliAnalysisTaskCLQA();

  void                        SetCentCL1In(TH1F *h)             { fCentCL1In       = h; }
  void                        SetCentV0AIn(TH1F *h)             { fCentV0AIn       = h; }
  void                        SetCumParams(Double_t Mmin, Double_t ptmin, Double_t ptmax, Double_t etamin, Double_t etamax);
  void                        SetDoCumulants(Bool_t b)          { fDoCumulants     = b; }
  void                        SetDoMuonTracking(Bool_t b)       { fDoMuonTracking  = b; }
  void                        SetDoTracking(Bool_t b)           { fDoTracking      = b; }
  void                        SetDo2013VertexCut(Bool_t b)      { fDo2013VertexCut = b; }

  void                        UserCreateOutputObjects();

 protected:
  Bool_t                      FillHistograms();
  Bool_t                      RetrieveEventObjects();
  Bool_t                      Run();
  void                        RunCumulants(Double_t Mmin, Double_t ptmin, Double_t ptmax, Double_t etamin, Double_t etamax);

  Bool_t                      fDo2013VertexCut;  // if true then use 2013 pA vertex check (only if 2013 pPb run)
  Bool_t                      fDoTracking;       // if true run tracking analysis
  Bool_t                      fDoMuonTracking;   // if true run muon tracking analysis
  Bool_t                      fDoCumulants;      // if true run cumulant analysis
  Bool_t                      fDoCumNtuple;      // if true then save cumulant ntuple
  Double_t                    fCumPtMin;         // minimum pt for cumulants
  Double_t                    fCumPtMax;         // maximum pt for cumulants
  Double_t                    fCumEtaMin;        // minimum eta for cumulants
  Double_t                    fCumEtaMax;        // maximum eta for cumulants
  Double_t                    fCumMmin;          // minimum number of tracks for cumulants 
  Int_t                       fCumMbins;         // number of bins for M
  TH1F                       *fCentCL1In;        // input for MC based CL1 centrality
  TH1F                       *fCentV0AIn;        // input for MC based V0A centrality
  TTree                      *fNtupCum;          //!ntuple for cumulant analysis
  AliNtupCumInfo             *fNtupCumInfo;      //!object holding cumulant results
  AliNtupZdcInfo             *fNtupZdcInfo;      //!object holding zdc info
  TH1                        *fHists[1000];      //!pointers to histograms

 private:
  Double_t                    DeltaPhi(Double_t phia, Double_t phib,
                                       Double_t rangeMin = -TMath::Pi()/2, 
                                       Double_t rangeMax = 3*TMath::Pi()/2) const;
  AliAnalysisTaskCLQA(const AliAnalysisTaskCLQA&);            // not implemented
  AliAnalysisTaskCLQA &operator=(const AliAnalysisTaskCLQA&); // not implemented

  ClassDef(AliAnalysisTaskCLQA, 6) // Constantin's Task
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
