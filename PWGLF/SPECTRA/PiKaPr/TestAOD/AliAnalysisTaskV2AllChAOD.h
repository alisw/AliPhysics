#ifndef ALIANALYSISTASKV2ALLCHAOD_H
#define ALIANALYSISTASKV2ALLCHAOD_H

/*  See cxx source for full Copyright notice */

//-------------------------------------------------------------------------
//                      AliAnalysisTaskV2AllChAOD
//
//
//
//
// Author: Leonardo Milano, CERN
//-------------------------------------------------------------------------

class AliAODEvent;
class AliSpectraAODTrackCuts;
class AliSpectraAODEventCuts;

#include "AliAnalysisTaskSE.h"
#include "TFile.h"
#include "TKey.h"
#include <TProfile.h>

class AliAnalysisTaskV2AllChAOD : public AliAnalysisTaskSE
{
public:
  // constructors
  AliAnalysisTaskV2AllChAOD() : AliAnalysisTaskSE(),
  fAOD(0x0),
  fTrackCuts(0x0),
  fEventCuts(0x0),
  fIsMC(0),
  fCharge(0),
  fVZEROside(0),
  fOutput(0x0),
  fOutput_lq(0x0),
  fOutput_sq(0x0),
  fnCentBins(20),
  fnQvecBins(40),
  fQvecUpperLim(100),
  fCutLargeQperc(9.),
  fCutSmallQperc(10.),
  fEtaGapMin(-0.5),
  fEtaGapMax(0.5),
  fTrkBit(128),
  fEtaCut(0.8),
  fMinPt(0.2),
  fMaxPt(20.0),
  fMinTPCNcls(70),
  fFillTHn(kTRUE),
  fCentrality(0),
  fQvector(0),
  fQvector_lq(0),
  fQvector_sq(0),
  fResSP(0),
  fResSP_vs_Cent(0),
  fEta_vs_Phi_bef(0),
  fEta_vs_PhiA(0),
  fEta_vs_PhiB(0),
  fResSP_lq(0),
  fResSP_vs_Cent_lq(0),
  fResSP_sq(0),
  fResSP_vs_Cent_sq(0),
  fResSP_inclusive(0),
  fv2SPGap1A_inclusive_mb(0),
  fv2SPGap1B_inclusive_mb(0),
  fv2SPGap1A_inclusive_lq(0),
  fv2SPGap1B_inclusive_lq(0),
  fv2SPGap1A_inclusive_sq(0),
  fv2SPGap1B_inclusive_sq(0),
  fResSPmc_inclusive(0),
  fv2SPGap1Amc_inclusive_mb(0),
  fv2SPGap1Bmc_inclusive_mb(0),
  fv2SPGap1Amc_inclusive_lq(0),
  fv2SPGap1Bmc_inclusive_lq(0),
  fv2SPGap1Amc_inclusive_sq(0),
  fv2SPGap1Bmc_inclusive_sq(0),
  fResGap1w(0),
  fV2IntGap1w(0),
  fResSP_qbin(0),
  fIsRecoEff(0),
  fRecoEffList(0),
  fQvecGen(0),
  fQgenType(0),
  fnNchBins(400),
  fDoCentrSystCentrality(0)
{
    for(Int_t j=0; j<9; j++){
      fv2SPGap1A[j]=0x0;
      fv2SPGap1B[j]=0x0;
      fSinGap1Aq[j]=0x0;
      fCosGap1Aq[j]=0x0;
      fSinGap1Bq[j]=0x0;
      fCosGap1Bq[j]=0x0;
      fSinGap1A[j]=0x0;
      fCosGap1A[j]=0x0;
      fSinGap1B[j]=0x0;
      fCosGap1B[j]=0x0;

      fv2SPGap1A_lq[j]=0x0;
      fv2SPGap1B_lq[j]=0x0;
      fSinGap1Aq_lq[j]=0x0;
      fCosGap1Aq_lq[j]=0x0;
      fSinGap1Bq_lq[j]=0x0;
      fCosGap1Bq_lq[j]=0x0;
      fSinGap1A_lq[j]=0x0;
      fCosGap1A_lq[j]=0x0;
      fSinGap1B_lq[j]=0x0;
      fCosGap1B_lq[j]=0x0;

      fv2SPGap1A_sq[j]=0x0;
      fv2SPGap1B_sq[j]=0x0;
      fSinGap1Aq_sq[j]=0x0;
      fCosGap1Aq_sq[j]=0x0;
      fSinGap1Bq_sq[j]=0x0;
      fCosGap1Bq_sq[j]=0x0;
      fSinGap1A_sq[j]=0x0;
      fCosGap1A_sq[j]=0x0;
      fSinGap1B_sq[j]=0x0;
      fCosGap1B_sq[j]=0x0;

      fResSP_vs_Qvec[j]=0x0;
      fV2IntGap1wq[j]=0x0;
    }
    
    for(Int_t j=0; j<10; j++){
      fv2SPGap1A_qbin[j]=0x0;
      fv2SPGap1B_qbin[j]=0x0;
    }
}
  AliAnalysisTaskV2AllChAOD(const char *name);
  virtual ~AliAnalysisTaskV2AllChAOD() {
    Printf("calling detructor of AliAnalysisTaskV2AllChAOD - To be implemented");
  }

  void SetIsMC(Bool_t isMC = kFALSE)    {fIsMC = isMC; };
  Bool_t GetIsMC()           const           { return fIsMC;};

  void SetCharge(Int_t charge = 0)    {fCharge = charge; };
  Int_t GetCharge()           const           { return fCharge;};

  void SetVZEROside(Int_t side = 0)    {fVZEROside = side; };
  Int_t GetVZEROside()           const           { return fVZEROside;};

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  AliSpectraAODTrackCuts      * GetTrackCuts()         {  return fTrackCuts; }
  AliSpectraAODEventCuts      * GetEventCuts()         {  return fEventCuts; }
  TList                          * GetOutputList()         { return fOutput; }

  void SetTrackCuts(AliSpectraAODTrackCuts * tc)       { fTrackCuts = tc; }
  void SetEventCuts(AliSpectraAODEventCuts * vc)       { fEventCuts = vc; }
  void SetnCentBins(Int_t val)                             { fnCentBins = val; }
  void SetnQvecBins(Int_t val)                             { fnQvecBins = val; }
  void SetQvecUpperLimit(Double_t val)                { fQvecUpperLim = val; }

  void SetTrackBits(UInt_t TrackBits) {fTrkBit=TrackBits;}
  void SetEtaCut(Double_t val) {fEtaCut=val;}
  void SetMinPt(Double_t val) {fMinPt=val;}
  void SetMaxPt(Double_t val) {fMaxPt=val;}
  void SetMinTPCNcls(Double_t val) {fMinTPCNcls=val;}

  Bool_t GetDCA(const AliAODTrack* trk, Double_t * p);
  void MCclosure(Double_t qvec);

  void EnableRecoEff (Bool_t val) { fIsRecoEff = val; }
  Double_t GetRecoEff(Double_t pt, Int_t iC);

  void SetRecoEffFile(TFile *f)    {
    TIter next(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
      TH1D * h=(TH1D*)key->ReadObj();
      fRecoEffList->Add(h);
    }
  };

  void     SetEtaGap(Float_t etamin,Float_t etamax)   { fEtaGapMin = etamin; fEtaGapMax = etamax; }
  void     SetQvecCut(Float_t qmin,Float_t qmax)      { fCutSmallQperc = qmin; fCutLargeQperc = qmax; }
  void     SetFillTHn (Bool_t val) { fFillTHn = val; }

  void SetQvecGen(Bool_t val) { fQvecGen = val; } //enable Qvec from generated
  void SetQgenType(Int_t val) { fQgenType = val ; } // type==0 qgen from tracks - type==1 qgen from vzero

  void SetnNchBins(Int_t val) { fnNchBins = val; }

  void SetDoCentrSystCentrality(Bool_t val) { fDoCentrSystCentrality = val; } //enable systematic for centrality

private:

  AliAODEvent                   * fAOD;                         //! AOD object
  AliSpectraAODTrackCuts      * fTrackCuts;                   // Track Cuts
  AliSpectraAODEventCuts      * fEventCuts;                   // Event Cuts
  Bool_t                          fIsMC;                         // true if processing MC
  Int_t                            fCharge;                      // charge to be selected
  Int_t                            fVZEROside;                  // 0: VZERO-A 1: VZERO-C
  TList                          * fOutput;                     // output list
  TList                          * fOutput_lq;                  // output list large Q
  TList                          * fOutput_sq;                  // output list small Q
  Int_t                            fnCentBins;                  // number of bins for the centrality axis
  Int_t                            fnQvecBins;                 // number of bins for the q vector axis
  Double_t                         fQvecUpperLim;             //Upper limit for Qvector

  Int_t                            fCutLargeQperc; // cut on 10% large Q-vec events
  Int_t                            fCutSmallQperc; // cut on 10% small Q-vec events

  Double_t fEtaGapMin;  // TBD
  Double_t fEtaGapMax;   // TBD

  UInt_t    fTrkBit;   // TBD
  Double_t  fEtaCut;   // TBD
  Double_t  fMinPt;   // TBD
  Double_t  fMaxPt;   // TBD
  Double_t  fMinTPCNcls;   // TBD

  Bool_t fFillTHn;   // TBD

  TH1D * fCentrality;   //! TBD
  TH1D * fQvector;   //! TBD
  TH1D * fQvector_lq;   //! TBD
  TH1D * fQvector_sq;   //! TBD

  //output object
  TProfile*     fResSP;             //! resolution
  TProfile*     fResSP_vs_Cent;   //! TBD
  TH2D*         fEta_vs_Phi_bef;        //! eta vs phi distribution before sub events cut 
  TH2D*         fEta_vs_PhiA;            //! eta vs phi distribution after sub events cut 
  TH2D*         fEta_vs_PhiB;            //! eta vs phi distribution after sub events cut 
  TProfile*     fv2SPGap1A[9];         //! v2{2} eta gap 1 for all events
  TProfile*     fv2SPGap1B[9];         //! v2{2} eta gap 1 for all events

  TProfile*     fSinGap1Aq[9];      //! <sin> vs pT gap 1
  TProfile*     fCosGap1Aq[9];      //! <cos> vs pT gap 1
  TProfile*     fSinGap1Bq[9];      //! <sin> vs pT gap 1
  TProfile*     fCosGap1Bq[9];      //! <cos> vs pT gap 1

  TProfile*     fSinGap1A[9];      //! <sin> vs pT gap 1
  TProfile*     fCosGap1A[9];      //! <cos> vs pT gap 1
  TProfile*     fSinGap1B[9];      //! <sin> vs pT gap 1
  TProfile*     fCosGap1B[9];      //! <cos> vs pT gap 1

  //large q
  TProfile*     fResSP_lq;             //! resolution
  TProfile*     fResSP_vs_Cent_lq;   //! TBD
  TProfile*     fv2SPGap1A_lq[9];         //! v2{2} eta gap 1 for all events
  TProfile*     fv2SPGap1B_lq[9];         //! v2{2} eta gap 1 for all events
  TProfile*     fSinGap1Aq_lq[9];      //! <sin> vs pT gap 1
  TProfile*     fCosGap1Aq_lq[9];      //! <cos> vs pT gap 1
  TProfile*     fSinGap1Bq_lq[9];      //! <sin> vs pT gap 1
  TProfile*     fCosGap1Bq_lq[9];      //! <cos> vs pT gap 1
  TProfile*     fSinGap1A_lq[9];      //! <sin> vs pT gap 1
  TProfile*     fCosGap1A_lq[9];      //! <cos> vs pT gap 1
  TProfile*     fSinGap1B_lq[9];      //! <sin> vs pT gap 1
  TProfile*     fCosGap1B_lq[9];      //! <cos> vs pT gap 1

  //small q
  TProfile*     fResSP_sq;             //! resolution
  TProfile*     fResSP_vs_Cent_sq;   //! TBD
  TProfile*     fv2SPGap1A_sq[9];         //! v2{2} eta gap 1 for all events
  TProfile*     fv2SPGap1B_sq[9];         //! v2{2} eta gap 1 for all events
  TProfile*     fSinGap1Aq_sq[9];      //! <sin> vs pT gap 1
  TProfile*     fCosGap1Aq_sq[9];      //! <cos> vs pT gap 1
  TProfile*     fSinGap1Bq_sq[9];      //! <sin> vs pT gap 1
  TProfile*     fCosGap1Bq_sq[9];      //! <cos> vs pT gap 1
  TProfile*     fSinGap1A_sq[9];      //! <sin> vs pT gap 1
  TProfile*     fCosGap1A_sq[9];      //! <cos> vs pT gap 1
  TProfile*     fSinGap1B_sq[9];      //! <sin> vs pT gap 1
  TProfile*     fCosGap1B_sq[9];      //! <cos> vs pT gap 1

  // MC closure test

  TProfile* fResSP_inclusive;   //! TBD
  TProfile* fv2SPGap1A_inclusive_mb;   //! TBD
  TProfile* fv2SPGap1B_inclusive_mb;   //! TBD
  TProfile* fv2SPGap1A_inclusive_lq;   //! TBD
  TProfile* fv2SPGap1B_inclusive_lq;   //! TBD
  TProfile* fv2SPGap1A_inclusive_sq;   //! TBD
  TProfile* fv2SPGap1B_inclusive_sq;   //! TBD

  TProfile* fResSPmc_inclusive;   //! TBD
  TProfile* fv2SPGap1Amc_inclusive_mb;   //! TBD
  TProfile* fv2SPGap1Bmc_inclusive_mb;   //! TBD
  TProfile* fv2SPGap1Amc_inclusive_lq;   //! TBD
  TProfile* fv2SPGap1Bmc_inclusive_lq;   //! TBD
  TProfile* fv2SPGap1Amc_inclusive_sq;   //! TBD
  TProfile* fv2SPGap1Bmc_inclusive_sq;   //! TBD

  // v2 vs qvec...

  TProfile*     fResGap1w;           //!
  TProfile*     fV2IntGap1w;         //! integrated v2 for gap 0.8 w
  TProfile*     fResSP_vs_Qvec[9];   //! TBD
  TProfile*     fV2IntGap1wq[9];     //!
  
  // v2 vs pt in q-vec bins
  
  TProfile*     fResSP_qbin;                //! resolution
  TProfile*     fv2SPGap1A_qbin[10];         //! v2{2} eta gap 1 for all events
  TProfile*     fv2SPGap1B_qbin[10];         //! v2{2} eta gap 1 for all events

  Bool_t fIsRecoEff;   // TBD
  TList * fRecoEffList; // reconstruction efficiency file

  Bool_t fQvecGen;  //enable Qvec from generated
  Int_t  fQgenType; // type==0 qgen from tracks - type==1 qgen from vzero
  Int_t  fnNchBins; //Ncharged

  Bool_t fDoCentrSystCentrality; //systematic check on centrality estimation


  AliAnalysisTaskV2AllChAOD(const AliAnalysisTaskV2AllChAOD&);
  AliAnalysisTaskV2AllChAOD& operator=(const AliAnalysisTaskV2AllChAOD&);

  ClassDef(AliAnalysisTaskV2AllChAOD, 18);
};

#endif
