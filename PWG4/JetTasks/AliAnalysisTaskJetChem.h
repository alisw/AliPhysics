/*************************************************************************
 *                                                                       *
 * Task for Jet Chemistry Analysis in PWG4 Jet Task Force Train          *
 *                                                                       *
 *                                                                       *
 * contact: Oliver Busch                                                 *
 * o.busch@gsi.de                                                        *
 *                                                                       *
 *************************************************************************/

#ifndef ALIANALYSISTASKJETCHEM_H
#define ALIANALYSISTASKJETCHEM_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliAnalysisTaskFragmentationFunction.h"

class AliAnalysisTaskJetChem : public AliAnalysisTaskFragmentationFunction {

 public:
  
  //----------------------------------------
  class AliFragFuncHistosInvMass : public TObject
  {
    
    public:
    
    AliFragFuncHistosInvMass(const char* name = "FFIMhistos", 
			     Int_t nJetPt = 0, Float_t jetPtMin = 0, Float_t jetPtMax = 0,
			     Int_t nInvMass = 0, Float_t invMassMin=0, Float_t invMassMax=0,  
			     Int_t nPt = 0, Float_t ptMin = 0, Float_t ptMax = 0,
			     Int_t nXi = 0, Float_t xiMin = 0, Float_t xiMax = 0,
			     Int_t nZ  = 0, Float_t zMin  = 0, Float_t zMax  = 0);
    AliFragFuncHistosInvMass(const AliFragFuncHistosInvMass& copy);
    AliFragFuncHistosInvMass& operator=(const AliFragFuncHistosInvMass &o);
    virtual ~AliFragFuncHistosInvMass();
    
    virtual void DefineHistos();
    virtual void FillFF(Float_t trackPt, Float_t invM, Float_t jetPt,Bool_t incrementJetPt);
    virtual void AddToOutput(TList* list) const;

  private:

    Int_t   fNBinsJetPt;    // FF histos bins
    Float_t fJetPtMin;      // FF histos limits
    Float_t fJetPtMax;      // FF histos limits
    Int_t   fNBinsInvMass;  // FF histos bins
    Float_t fInvMassMin;    // FF histos limits
    Float_t fInvMassMax;    // FF histos limits
    Int_t   fNBinsPt;       // FF histos bins
    Float_t fPtMin;         // FF histos limits
    Float_t fPtMax;         // FF histos limits
    Int_t   fNBinsXi;       // FF histos bins
    Float_t fXiMin;         // FF histos limits
    Float_t fXiMax;         // FF histos limits
    Int_t   fNBinsZ;        // FF histos bins
    Float_t fZMin;          // FF histos limits
    Float_t fZMax;          // FF histos limits
  
    TH3F*   fh3TrackPt;     //! FF: track transverse momentum 
    TH3F*   fh3Xi;          //! FF: xi 
    TH3F*   fh3Z;           //! FF: z  
    TH1F*   fh1JetPt;       //! jet pt 

    TString fNameFF;        // histo names prefix
    
    ClassDef(AliFragFuncHistosInvMass, 1);
  };
  

 //----------------------------------------
  class AliFragFuncHistosPhiCorrInvMass : public TObject
  {
				   
    public:
    
    AliFragFuncHistosPhiCorrInvMass(const char* name = "FFPhiCorrIMhistos", 
				    Int_t nPt = 0, Float_t ptMin = 0, Float_t ptMax = 0,
				    Int_t nPhi = 0, Float_t phiMin = 0, Float_t phiMax = 0,
				    Int_t nInvMass = 0, Float_t invMassMin=0, Float_t invMassMax=0);

    AliFragFuncHistosPhiCorrInvMass(const AliFragFuncHistosPhiCorrInvMass& copy);
    AliFragFuncHistosPhiCorrInvMass& operator=(const AliFragFuncHistosPhiCorrInvMass &o);
    virtual ~AliFragFuncHistosPhiCorrInvMass();
    
    virtual void DefineHistos();
    virtual void FillPhiCorr(Float_t pt, Float_t phi, Float_t invM);
    virtual void AddToOutput(TList* list) const;

  private:

    Int_t   fNBinsPt;       // FF histos bins
    Float_t fPtMin;         // FF histos limits
    Float_t fPtMax;         // FF histos limits

    Int_t   fNBinsPhi;      // FF histos bins
    Float_t fPhiMin;        // FF histos limits
    Float_t fPhiMax;        // FF histos limits
    
    Int_t   fNBinsInvMass;  // FF histos bins
    Float_t fInvMassMin;    // FF histos limits
    Float_t fInvMassMax;    // FF histos limits
  
    TH3F*   fh3PhiCorr;     //! FF: phi correlation histo 

    TString fNamePhiCorr;   // histo names prefix
    
    ClassDef(AliFragFuncHistosPhiCorrInvMass, 1);
  };
  
  //----------------------------------------

  AliAnalysisTaskJetChem(); 
  AliAnalysisTaskJetChem(const char *name);
  AliAnalysisTaskJetChem(const  AliAnalysisTaskJetChem &copy);
  AliAnalysisTaskJetChem& operator=(const  AliAnalysisTaskJetChem &o);
  virtual ~AliAnalysisTaskJetChem();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);

  static  void   SetProperties(TH3F* h,const char* x, const char* y,const char* z);

  Bool_t IsAccepteddEdx(const Double_t mom,const Double_t signal, AliPID::EParticleType n, const Double_t cutnSig) const;
  Bool_t IsK0InvMass(const Double_t mass) const; 
  Int_t  GetListOfK0s(TList *list, const Int_t type);
  virtual void SetK0Type(Int_t i){ fK0Type = i; }
  virtual void SetFilterMaskK0(UInt_t i) {fFilterMaskK0 = i;}


  void   SetFFInvMassHistoBins(Int_t nJetPt = 19, Float_t jetPtMin = 5, Float_t jetPtMax = 100, 
			       Int_t nInvM = 50, Float_t invMMin = 0.450,  Float_t invMMax = 0.550,
			       Int_t nPt = 20, Float_t ptMin = 0., Float_t ptMax = 20., 
			       Int_t nXi = 35, Float_t xiMin = 0., Float_t xiMax = 7.,
			       Int_t nZ = 11,  Float_t zMin = 0.,  Float_t zMax = 1.1)
  { fFFIMNBinsJetPt = nJetPt; fFFIMJetPtMin = jetPtMin; fFFIMJetPtMax = jetPtMax; 
    fFFIMNBinsInvM = nInvM; fFFIMInvMMin = invMMin; fFFIMInvMMax = invMMax; fFFIMNBinsPt = nPt; fFFIMPtMin = ptMin; fFFIMPtMax = ptMax; 
    fFFIMNBinsXi = nXi; fFFIMXiMin = xiMin; fFFIMXiMax = xiMax; fFFIMNBinsZ  = nZ;  fFFIMZMin  = zMin;  fFFIMZMax  = zMax; }

  void   SetPhiCorrInvMassHistoBins(Int_t nPt = 40, Float_t ptMin = 0., Float_t ptMax = 20., 
				    Int_t nPhi = 20, Float_t phiMin = 0., Float_t phiMax = 2*TMath::Pi(),
				    Int_t nInvM = 50, Float_t invMMin = 0.450,  Float_t invMMax = 0.550)
				    
  { fPhiCorrIMNBinsPt = nPt; fPhiCorrIMPtMin = ptMin; fPhiCorrIMPtMax = ptMax;
    fPhiCorrIMNBinsPhi = nPhi; fPhiCorrIMPhiMin = phiMin; fPhiCorrIMPhiMax = phiMax;
    fPhiCorrIMNBinsInvM = nInvM; fPhiCorrIMInvMMin = invMMin; fPhiCorrIMInvMMax = invMMax;
  }
  
  
  // consts
  enum { kTrackUndef =0, kOnFly, kOnFlyPID, kOnFlydEdx, kOnFlyPrim, kOffl, kOfflPID, kOffldEdx, kOfflPrim };  
  
 private:
  
  Int_t fK0Type;                                           //! K0 cuts
  UInt_t fFilterMaskK0;                                    //! K0 legs cuts
  TList* fListK0s;                                         //! K0 list 

  AliFragFuncQATrackHistos*  fV0QAK0;                      //! track QA: V0s in K0 inv mass range
  AliFragFuncHistos*         fFFHistosRecCutsK0Evt;        //! inclusive FF for K0 evt
  AliFragFuncHistosInvMass*  fFFHistosIMK0AllEvt;          //! K0 pt spec for all events
  AliFragFuncHistosInvMass*  fFFHistosIMK0Jet;             //! K0 FF all dPhi   
  AliFragFuncHistosInvMass*  fFFHistosIMK0Cone;            //! K0 FF jet cone   
  AliFragFuncHistosPhiCorrInvMass*  fFFHistosPhiCorrIMK0;  //! K0 corelation to jet axis 

  // histogram bins  

  Int_t   fFFIMNBinsJetPt;    // FF histos bins
  Float_t fFFIMJetPtMin;      // FF histos limits
  Float_t fFFIMJetPtMax;      // FF histos limits

  Int_t   fFFIMNBinsInvM;     // FF histos bins
  Float_t fFFIMInvMMin;       // FF histos bins
  Float_t fFFIMInvMMax;       // FF histos bins

  Int_t   fFFIMNBinsPt;       // FF histos bins
  Float_t fFFIMPtMin;         // FF histos limits
  Float_t fFFIMPtMax;         // FF histos limits

  Int_t   fFFIMNBinsXi;       // FF histos bins
  Float_t fFFIMXiMin;         // FF histos limits
  Float_t fFFIMXiMax;         // FF histos limits

  Int_t   fFFIMNBinsZ;        // FF histos bins
  Float_t fFFIMZMin;          // FF histos limits
  Float_t fFFIMZMax;          // FF histos limits


  Int_t fPhiCorrIMNBinsPt;    // FF histos bins
  Float_t fPhiCorrIMPtMin;    // FF histos limits
  Float_t fPhiCorrIMPtMax;    // FF histos limits

  Int_t fPhiCorrIMNBinsPhi;   // FF histos bins
  Float_t fPhiCorrIMPhiMin;   // FF histos limits
  Float_t fPhiCorrIMPhiMax;   // FF histos limits
		
  Int_t fPhiCorrIMNBinsInvM;  // FF histos bins
  Float_t fPhiCorrIMInvMMin;  // FF histos limits
  Float_t fPhiCorrIMInvMMax;  // FF histos limits
  


  // Histograms

  TH1F* fh1K0Mult;                  //!
  TH1F* fh1dPhiJetK0;               //!


  ClassDef(AliAnalysisTaskJetChem, 3);
};

#endif
