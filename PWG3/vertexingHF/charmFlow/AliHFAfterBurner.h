#ifndef ALIANALYSISTASKAFTERBURNER_H
#define ALIANALYSISTASKAFTERBURNER_H

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// AliAnalysisTaskSEHFv2 gives the needed tools for the D 
// mesons v2 analysis 
// Authors: Giacomo Ortona, ortona@to.infn.it,
//          
//*************************************************************************

/* $Id$ */

class TH1F;
class TH2D;
class AliRDHFCuts;
class TVector2;

class AliHFAfterBurner : public TObject
{
 public:

  AliHFAfterBurner();
  AliHFAfterBurner(Int_t decChannel);

  virtual ~AliHFAfterBurner();

  Double_t GetNewAngle(AliAODRecoDecayHF *d,TClonesArray *arrayMC);
  Double_t NewtonMethodv2(Double_t phi,Double_t v2,Double_t phi0);
  Double_t NewtonMethodv2(Double_t phi,Double_t v2){
    return NewtonMethodv2(phi,v2,phi);}
  void SetMCv2(Float_t v2sig,Float_t v2bkg);
  Int_t CheckOrigin(const AliAODMCParticle* mcPart,TClonesArray *arrayMC) const;
  Double_t GetPhi(Double_t phi,Float_t v2);
  Float_t GetPhi02Pi(Float_t phi);
  void SetPrecision(Bool_t prec){fPrecisionNewton=prec;}
  void SetUseNewton(Bool_t newt){fUseNewton=newt;}
  void SetDecChannel(Int_t decch);
  void SetEventPlane(Double_t ep){fEventPlane=ep;}
  void SetEventPlane();
  void SetEventPlaneMethod(Int_t method);

  Bool_t GetIsSignal(){return fSignal;}
  Double_t GetEventPlane(){return fEventPlane;}
 private:
  AliHFAfterBurner(const AliHFAfterBurner &source);
  AliHFAfterBurner& operator=(const AliHFAfterBurner& source); 

  Float_t fSigv2;       //v2 for signal
  Float_t fBkgv2;       //v2 for background
  Bool_t fUseNewton;    //Switch between Newton and fast methods
  Double_t fPrecisionNewton;//precision to stop Newton method
  Int_t fDecChannel;   //D2H decay channel
  Bool_t fSignal;      //kTRUE if candidate is signal
  Double_t fEventPlane; //event plane
  Int_t fMethodEP;     //switch between different EP method calculation

  ClassDef(AliHFAfterBurner,1); // AliFlowAB class
};
#endif
