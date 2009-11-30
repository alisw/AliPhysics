#ifndef ALIRESONANCEKINKLIKESIGN_H
#define ALIRESONANCEKINKLIKESIGN_H

/*  See cxx source for full Copyright notice */

//--------------------------------------------------------------------------------
//                   class AliResonanceKinkLikeSign
//         This task is an example of an analysis task
//        for producing a like-sign background for resonances having at least one 
//        kaon-kink in their decay products. 
//        Author: Paraskevi Ganoti, University of Athens (pganoti@phys.uoa.gr)
//---------------------------------------------------------------------------------

#include "AliAnalysisTaskSE.h"
#include "TVector3.h"

class TF1;
class TH1D;
class TH2D;
class AliESDtrack;
class AliESDVertex;
class AliESDEvent;

class AliResonanceKinkLikeSign : public AliAnalysisTaskSE {
 public:
  AliResonanceKinkLikeSign(const char *name = "AliResonanceKinkLikeSign");
  virtual ~AliResonanceKinkLikeSign() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void SetPDGCodes(Int_t d1, Int_t d2) {fdaughter1pdg=d1; fdaughter2pdg=d2;}
  void SetHistoSettings(Int_t nbins, Float_t nlowx, Float_t nhighx, Int_t nptbins, Float_t nlowpt, Float_t nupperpt) {fnbins=nbins; fnlowx=nlowx; fnhighx=nhighx; fptbins=nptbins; flowpt=nlowpt; fupperpt=nupperpt;}
  void SetEtaLimits(Float_t nloweta, Float_t nuppereta) {floweta=nloweta; fuppereta=nuppereta;} 
  Float_t GetSigmaToVertex(AliESDtrack* esdTrack) const ; 
  const AliESDVertex *GetEventVertex(const AliESDEvent* esd) const;
  void SetDebugLevel(Int_t level) {fDebug = level;}
  Bool_t IsAcceptedForKink(AliESDEvent *localesd, const AliESDVertex *localvertex, AliESDtrack *localtrack);
  Bool_t IsAcceptedForTrack(AliESDEvent *localesd, const AliESDVertex *localvertex, AliESDtrack *localtrack);
  Bool_t IsKink(AliESDEvent *localesd, Int_t kinkIndex, TVector3 trackMom); 
  
  void SetMaxNsigmaToVertex(Float_t maxNSigmaToVertex) {
   fMaxNSigmaToVertex=maxNSigmaToVertex;
  }
  Float_t GetMaxNsigmaToVertex() const {return fMaxNSigmaToVertex;} 
  
  void SetPtTrackCut(Double_t minPtTrackCut) {
   fMinPtTrackCut=minPtTrackCut;
  }
  Double_t GetPtTrackCut() const {return fMinPtTrackCut;} 
  
  void SetMaxDCAxy(Double_t maxDCAxy) {
   fMaxDCAxy=maxDCAxy;
  }
  Double_t GetMaxDCAxy() const {return fMaxDCAxy;}
  
  void SetMaxDCAzaxis(Double_t maxDCAzaxis) {
   fMaxDCAzaxis=maxDCAzaxis;
  }
  Double_t GetMaxDCAzaxis() const {return fMaxDCAzaxis;}   
  
  void SetMinTPCclusters(Int_t minTPCclusters) {
   fMinTPCclusters=minTPCclusters;
  }
  Int_t GetMinTPCclusters() const {return fMinTPCclusters;}    
  
  void SetMaxChi2PerTPCcluster(Double_t maxChi2PerTPCcluster) {
   fMaxChi2PerTPCcluster=maxChi2PerTPCcluster;
  }
  Double_t GetMaxChi2PerTPCcluster() const {return fMaxChi2PerTPCcluster;} 
  
  void SetMaxCov0(Double_t maxCov0) {
   fMaxCov0=maxCov0;
  }
  Double_t GetMaxCov0() const {return fMaxCov0;}     
  
  void SetMaxCov2(Double_t maxCov2) {
   fMaxCov2=maxCov2;
  }
  Double_t GetMaxCov2() const {return fMaxCov2;}   
  
  void SetMaxCov5(Double_t maxCov5) {
   fMaxCov5=maxCov5;
  }
  Double_t GetMaxCov5() const {return fMaxCov5;}   
  
  void SetMaxCov9(Double_t maxCov9) {
   fMaxCov9=maxCov9;
  }
  Double_t GetMaxCov9() const {return fMaxCov9;}   
  
  void SetMaxCov14(Double_t maxCov14) {
   fMaxCov14=maxCov14;
  }
  Double_t GetMaxCov14() const {return fMaxCov14;}
   
  void SetMinKinkRadius(Float_t minKinkRadius) {
   fminKinkRadius=minKinkRadius;
  }
  Float_t GetMinKinkRadius() const {return fminKinkRadius;}

  void SetMaxKinkRadius(Float_t maxKinkRadius) {
   fmaxKinkRadius=maxKinkRadius;
  }
  Float_t GetMaxKinkRadius() const {return fmaxKinkRadius;} 
  
  void SetQtLimits(Float_t minQt, Float_t maxQt) {fminQt=minQt; fmaxQt=maxQt;}     
  
 private:
  Int_t       fDebug;        //  Debug flag
 // AliESDEvent *fESD;    // ESD object
  TList       *fListOfHistos; // List 
  TF1         *f1; //upper limit curve for the decay K->mu
  TF1         *f2;  //upper limit curve for the decay pi->mu
  TH1D        *fPosKaonLikeSign; // negative spectrum
  TH2D        *fLikeSignInvmassPt; // negative spectrum
  Float_t     fMaxNSigmaToVertex; // standard cut to select primary tracks
  Double_t    fMinPtTrackCut; // lower pt cut for the tracks
  Double_t    fMaxDCAxy; // impact parameter in the xy plane
  Double_t    fMaxDCAzaxis; // impact parameter in the z axis
  Int_t       fMinTPCclusters; // standard cut for the TPC clusters
  Double_t    fMaxChi2PerTPCcluster; // standard cut for the chi2 of the TPC clusters
  Double_t    fMaxCov0;  // standard cut
  Double_t    fMaxCov2; // standard cut
  Double_t    fMaxCov5; // standard cut
  Double_t    fMaxCov9; // standard cut
  Double_t    fMaxCov14; // standard cut
  Int_t       fdaughter1pdg; // pdg code of the resonance's first daughter
  Int_t       fdaughter2pdg;  // pdg code of the resonance's second daughter  
  Int_t       fnbins; // Inv mass histo number of bins
  Float_t     fnlowx; // Inv mass histo lower limit
  Float_t     fnhighx; // Inv mass histo upper limit
  Float_t     floweta; // lower eta limit
  Float_t     fuppereta; // upper eta limit
  Float_t     fminKinkRadius; // min accepted radius for the kink vertex
  Float_t     fmaxKinkRadius; // max accepted radius for the kink vertex
  Float_t     fminQt; //min Qt cut
  Float_t     fmaxQt; //max Qt cut  
  Int_t       fptbins; // number of bins in pt
  Float_t     flowpt; // pt lower limit
  Float_t     fupperpt; // pt upper limit 
    
  AliResonanceKinkLikeSign(const AliResonanceKinkLikeSign&); // not implemented
  AliResonanceKinkLikeSign& operator=(const AliResonanceKinkLikeSign&); // not implemented

  ClassDef(AliResonanceKinkLikeSign, 1); // example of analysis
};

#endif
