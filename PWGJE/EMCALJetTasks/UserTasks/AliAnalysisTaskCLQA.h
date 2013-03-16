#ifndef ALIANALYSISTASKCLQA_H
#define ALIANALYSISTASKCLQA_H

// $Id: AliAnalysisTaskCLQA.h 58847 2012-09-30 18:11:49Z loizides $

class TClonesArray;
class TString;
class TH1F;
class TH2F;
class TH3F;
class TNtuple;
class TNtupleD;
class TTree;

#include "AliAnalysisTaskEmcalJet.h"

class AliNtupCumInfo;

class AliAnalysisTaskCLQA : public AliAnalysisTaskEmcalJet {
 public:
  AliAnalysisTaskCLQA();
  AliAnalysisTaskCLQA(const char *name);
  virtual ~AliAnalysisTaskCLQA();

  void                        UserCreateOutputObjects();
  void                        SetDoCumulants(Bool_t b)          { fDoCumulants = b; }
  void                        SetCumParams(Double_t Mmin, Double_t ptmin, Double_t ptmax, Double_t etamin, Double_t etamax);

 protected:
  Bool_t                      FillHistograms();
  Bool_t                      RetrieveEventObjects();
  Bool_t                      Run();
  void                        RunCumulants(Double_t Mmin, Double_t ptmin, Double_t ptmax, Double_t etamin, Double_t etamax);

  Bool_t                      fDoCumulants;      // if true run cumulant analysis
  Double_t                    fCumPtMin;         // minimum pt for cumulants
  Double_t                    fCumPtMax;         // maximum pt for cumulants
  Double_t                    fCumEtaMin;        // minimum eta for cumulants
  Double_t                    fCumEtaMax;        // maximum eta for cumulants
  Double_t                    fCumMmin;          // minimum number of tracks for cumulants 
  TTree                      *fNtupCum;          //!ntuple for cumulant analysis
  AliNtupCumInfo             *fNtupCumInfo;      //!object holding cumulant results

 private:
  AliAnalysisTaskCLQA(const AliAnalysisTaskCLQA&);            // not implemented
  AliAnalysisTaskCLQA &operator=(const AliAnalysisTaskCLQA&); // not implemented

  ClassDef(AliAnalysisTaskCLQA, 1) // Constantin's Task
};

class AliNtupCumInfo {
 public:
    AliNtupCumInfo() : fTrig(0), fRun(0), fVz(0), fIsFEC(0), fIsVSel(0), fIsP(0),
                       fMall(0), fMall2(0), fPtMaxall(0), fMPtall(0), 
                       fMPt2all(0), fMPtall2(0), fTSall(0),
                       fM(0), fQ2abs(0), fQ4abs(0), fQ42re(0),
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
  Double_t      fQ2abs;        // Q2 absolute
  Double_t      fQ4abs;        // Q4 absolute
  Double_t      fQ42re;        // Re(Q2Q*Q*)
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
 
  ClassDef(AliNtupCumInfo,2) // Cumulant storage class
};

#endif
