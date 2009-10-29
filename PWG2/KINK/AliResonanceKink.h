#ifndef ALIRESONANCEKINK_H
#define ALIRESONANCEKINK_H

/*  See cxx source for full Copyright notice */

//------------------------------------------------------------------------------
//                   class AliResonanceKink
//         This task is an example of an analysis task
//        for analysing resonances having one kaon kink
//Author: Paraskevi Ganoti, University of Athens (pganoti@phys.uoa.gr)
//------------------------------------------------------------------------------
#include "TVector3.h"
class TF1;
class TH1D;
class TH2D;
class AliESDEvent;
class AliESDtrack;
class AliESDVertex;
class AliMCEvent;
class TList;
class TString;

class AliResonanceKink : public TObject {
 public:
 
  enum ResonanceType {kPhi=333, kKstar0=313, kLambda1520=3124};
  enum DaughterType {kdaughterPion=211, kdaughterKaon=321, kdaughterProton=2212};
  
  AliResonanceKink();
  AliResonanceKink(Int_t nbins, Float_t nlowx, Float_t nhighx, Int_t daughter1, Int_t daughter2, Int_t resonancePDG);
  virtual ~AliResonanceKink();
  
  TList* GetHistogramList();
  void Analyse(AliESDEvent* esd, AliMCEvent* mcEvent);  
  Float_t GetSigmaToVertex(AliESDtrack* esdTrack) const ; 
  const AliESDVertex *GetEventVertex(const AliESDEvent* esd) const;
  void SetDebugLevel(Int_t level) {fDebug = level;}
  void SetAnalysisType(TString type) {fAnalysisType=type;}
  void SetPDGCodes(Int_t d1, Int_t d2, Int_t res) {fdaughter1pdg=d1; fdaughter2pdg=d2; fresonancePDGcode=res;}
  void InitOutputHistograms(Int_t nbins, Float_t nlowx, Float_t nhighx);
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
  
  //void SetTPCrefit() {Int_t fTPCrefitFlag=kTRUE;}
     
 private:
 
  Int_t       fDebug;        //  Debug flag
  TList       *fListOfHistos; // List 
  TH1D        *fOpeningAngle; // Opening  
  TH1D        *fInvariantMass; // invMass spectrum   
  TH1D        *fInvMassTrue; // invMassTrue spectrum  
  TH1D        *fPhiBothKinks; // bothKaonsKinks   
  TH1D        *fRecPt; // pT spectrum  
  TH1D        *fRecEta; // Eta spectrum
  TH2D        *fRecEtaPt; // Eta pT spectrum  
  TH1D        *fSimPt; // pT Sim spectrum  
  TH1D        *fSimEta; // Eta Sim spectrum
  TH2D        *fSimEtaPt; // Eta pT Sim spectrum 
  TH1D        *fSimPtKink; // pT Sim one kaon kink spectrum  
  TH1D        *fSimEtaKink; // Eta Sim one kaon kink spectrum spectrum
  TH2D        *fSimEtaPtKink; // Eta pT Sim one kaon kink spectrum   
  TH1D        *fhdr ; // radial impact  
  TH1D        *fhdz ; // z impact
  TF1         *f1; //upper limit curve for the decay K->mu
  TF1         *f2;  //upper limit curve for the decay pi->mu
  TString     fAnalysisType;//"ESD" or "MC"
  TH1D        *fvtxz ; // vtx z component
  Int_t       fNbins; // bins
  Float_t     fLowX;// lowx
  Float_t     fHighX; // high x
  Int_t       fdaughter1pdg; // pdg code of the resonance's first daughter
  Int_t       fdaughter2pdg;  // pdg code of the resonance's second daughter
  Int_t       fresonancePDGcode; // pdg code of the resonance
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
  TH2D        *fInvmassPt;  // 2D histo for invariant mass calculation in pt bins (all pairs, ESD) 
  TH2D        *fInvmassPtTrue;  // 2D histo for invariant mass calculation in pt bins (true pairs, ESD) 
  TH2D        *fMCInvmassPt; // 2D histo for invariant mass calculation in pt bins (all pairs, MC)
  TH2D        *fMCInvmassPtTrue;  // 2D histo for invariant mass calculation in pt bins (true pairs, MC)       
//  Bool_t      fTPCrefitFlag;
  
  AliResonanceKink(const AliResonanceKink&); // not implemented
  AliResonanceKink& operator=(const AliResonanceKink&); // not implemented

  ClassDef(AliResonanceKink, 1); // example of analysis
};

#endif
