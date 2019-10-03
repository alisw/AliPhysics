#ifndef ALIJETBKG_H
#define ALIJETBKG_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//--------------------------------------------------
// Method implementation for background studies and background subtraction with UA1 algorithms
//
// Author: magali.estienne@subatech.in2p3.fr
//-------------------------------------------------

class TH1F;
class TH2F;
class TList;
class AliAODJetEventBackground;

class AliJetBkg : public TObject
{
 public:
  AliJetBkg();
  AliJetBkg(const AliJetBkg &input);
  ~AliJetBkg();
  void    SetHeader(AliJetHeader *header)  {fHeader=header;}
  void    SetCalTrkEvent(AliJetCalTrkEvent *evt)  {fEvent=evt;}
  Bool_t  PtCutPass(Int_t id, Int_t nTracks);
  Bool_t  SignalCutPass(Int_t id, Int_t nTracks);
  Float_t CalcJetAreaEtaCut(Float_t radius, Float_t etaJet);
  void    CalcJetAndBckgAreaEtaCut(Bool_t calcOutsideArea, Float_t rc, Int_t nJ, const Float_t* etaJet, Float_t* &areaJet, Float_t &areaOut);
 
  void    SubtractBackg(const Int_t& nIn, const Int_t&nJ, Float_t&EtbgTotalN, Float_t&sigmaN, 
			const Float_t* ptT, const Float_t* etaT, const Float_t* phiT, 
			Float_t* etJet, const Float_t* etaJet, const Float_t* phiJet, 
			Float_t* etsigJet, Int_t* multJetT, Int_t* multJetC, Int_t* multJet, 
			Int_t* injet, Float_t* &areaJet);
  
  void    SubtractBackgCone(const Int_t& nIn, const Int_t&nJ,Float_t& EtbgTotalN, Float_t&sigmaN,
			    const Float_t* ptT, const Float_t* etaT, const Float_t* phiT, Float_t* etJet, 
			    const Float_t* etaJet, const Float_t* phiJet, Float_t* etsigJet, 
			    Int_t* multJetT, Int_t* multJetC, Int_t* multJet, Int_t* injet, Float_t* &/*areaJet*/);

  void    SubtractBackgRatio(const Int_t& nIn, const Int_t&nJ,Float_t& EtbgTotalN, Float_t&sigmaN,
			     const Float_t* ptT, const Float_t* etaT, const Float_t* phiT, Float_t* etJet, 
			     const Float_t* etaJet, const Float_t* phiJet, Float_t* etsigJet, 
			     Int_t* multJetT, Int_t* multJetC, Int_t* multJet, Int_t* injet, Float_t* &/*areaJet*/);

  void    SubtractBackgStat(const Int_t& nIn, const Int_t&nJ,Float_t&EtbgTotalN, Float_t&sigmaN,
			    const Float_t* ptT, const Float_t* etaT, const Float_t* phiT, Float_t* etJet, 
			    const Float_t* etaJet, const Float_t* phiJet, Float_t* etsigJet, 
			    Int_t* multJetT, Int_t* multJetC, Int_t* multJet, Int_t* injet, Float_t* &areaJet);
  void    SetDebug(Int_t debug){fDebug = debug;}

  enum {kMaxJets = 60};

 private:
  AliJetBkg& operator=(const AliJetBkg& source); // not implemented
  //    Double_t CalcRho(vector<fastjet::PseudoJet> input_particles,Double_t RparamBkg,TString method);

  AliJetCalTrkEvent* fEvent;    //! reader
  AliJetHeader*      fHeader;   //! header
  Int_t              fDebug;    //  Debug option

  // temporary histos for background, reset for each event, no need to stream
  TH1F*  fhEtJet[kMaxJets];   //! histogram for background subtraction
  TH1F*  fhAreaJet[kMaxJets]; //! histogram for background subtraction (store global not to create it with every event
  TH1F*  fhEtBackg;           //! histogram for background subtraction
  TH1F*  fhAreaBackg;         //! histogram for background subtraction

  ClassDef(AliJetBkg, 2)      // background jet analysis

};
 
#endif
