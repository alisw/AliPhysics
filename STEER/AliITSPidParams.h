#ifndef ALIITSPIDPARAMS_H
#define ALIITSPIDPARAMS_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to store parameters of ITS response funcions            //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include <TFormula.h>
#include <TNamed.h>

class AliITSPidParams : public TNamed {

 public:
  AliITSPidParams();
  AliITSPidParams(Char_t * name);
  ~AliITSPidParams();

  void InitMC();
  Double_t GetLandauGausNormPdgCode(Double_t dedx, Int_t pdgCode, Double_t mom, Int_t lay) const;
  Double_t GetLandauGausNorm(Double_t dedx, Int_t partType, Double_t mom, Int_t lay) const;

  // pion setters
  void SetSDDPionMPV(const TFormula* form){
    if(fSDDPionMPV) delete fSDDPionMPV;
    fSDDPionMPV=new TFormula(*form);
  }
  void SetSDDPionLandauWidth(const TFormula* form){
    if(fSDDPionLandauWidth) delete fSDDPionLandauWidth;
    fSDDPionLandauWidth=new TFormula(*form);
  }
  void SetSDDPionGaussWidth(const TFormula* form){
    if(fSDDPionGaussWidth) delete fSDDPionGaussWidth;
    fSDDPionGaussWidth=new TFormula(*form);
  }
  void SetSSDPionMPV(const TFormula* form){
    if(fSSDPionMPV) delete fSSDPionMPV;
    fSSDPionMPV=new TFormula(*form);
  }
  void SetSSDPionLandauWidth(const TFormula* form){
    if(fSSDPionLandauWidth) delete fSSDPionLandauWidth;
    fSSDPionLandauWidth=new TFormula(*form);
  }
  void SetSSDPionGaussWidth(const TFormula* form){
    if(fSSDPionGaussWidth) delete fSSDPionGaussWidth;
    fSSDPionGaussWidth=new TFormula(*form);
  }

  // kaon setters
  void SetSDDKaonMPV(const TFormula* form){
    if(fSDDKaonMPV) delete fSDDKaonMPV;
    fSDDKaonMPV=new TFormula(*form);
  }
  void SetSDDKaonLandauWidth(const TFormula* form){
    if(fSDDKaonLandauWidth) delete fSDDKaonLandauWidth;
    fSDDKaonLandauWidth=new TFormula(*form);
  }
  void SetSDDKaonGaussWidth(const TFormula* form){
    if(fSDDKaonGaussWidth) delete fSDDKaonGaussWidth;
    fSDDKaonGaussWidth=new TFormula(*form);
  }
  void SetSSDKaonMPV(const TFormula* form){
    if(fSSDKaonMPV) delete fSSDKaonMPV;
    fSSDKaonMPV=new TFormula(*form);
  }
  void SetSSDKaonLandauWidth(const TFormula* form){
    if(fSSDKaonLandauWidth) delete fSSDKaonLandauWidth;
    fSSDKaonLandauWidth=new TFormula(*form);
  }
  void SetSSDKaonGaussWidth(const TFormula* form){
    if(fSSDKaonGaussWidth) delete fSSDKaonGaussWidth;
    fSSDKaonGaussWidth=new TFormula(*form);
  }


  // proton setters
  void SetSDDProtMPV(const TFormula* form){
    if(fSDDProtMPV) delete fSDDProtMPV;
    fSDDProtMPV=new TFormula(*form);
  }
  void SetSDDProtLandauWidth(const TFormula* form){
    if(fSDDProtLandauWidth) delete fSDDProtLandauWidth;
    fSDDProtLandauWidth=new TFormula(*form);
  }
  void SetSDDProtGaussWidth(const TFormula* form){
    if(fSDDProtGaussWidth) delete fSDDProtGaussWidth;
    fSDDProtGaussWidth=new TFormula(*form);
  }
  void SetSSDProtMPV(const TFormula* form){
    if(fSSDProtMPV) delete fSSDProtMPV;
    fSSDProtMPV=new TFormula(*form);
  }
  void SetSSDProtLandauWidth(const TFormula* form){
    if(fSSDProtLandauWidth) delete fSSDProtLandauWidth;
    fSSDProtLandauWidth=new TFormula(*form);
  }
  void SetSSDProtGaussWidth(const TFormula* form){
    if(fSSDProtGaussWidth) delete fSSDProtGaussWidth;
    fSSDProtGaussWidth=new TFormula(*form);
  }


  // pion getters
  Double_t GetSDDPionMPV(Double_t mom) const {
    return fSDDPionMPV->Eval(mom);
  }
  Double_t GetSDDPionLandauWidth(Double_t mom) const {
    return fSDDPionLandauWidth->Eval(mom);
  }
  Double_t GetSDDPionGaussWidth(Double_t mom) const {
    return fSDDPionGaussWidth->Eval(mom);
  }
  Double_t GetSSDPionMPV(Double_t mom) const {
    return fSSDPionMPV->Eval(mom);
  }
  Double_t GetSSDPionLandauWidth(Double_t mom) const {
    return fSSDPionLandauWidth->Eval(mom);
  }
  Double_t GetSSDPionGaussWidth(Double_t mom) const {
    return fSSDPionGaussWidth->Eval(mom);
  }

  // kaon getters
  Double_t GetSDDKaonMPV(Double_t mom) const {
    return fSDDKaonMPV->Eval(mom);
  }
  Double_t GetSDDKaonLandauWidth(Double_t mom) const {
    return fSDDKaonLandauWidth->Eval(mom);
  }
  Double_t GetSDDKaonGaussWidth(Double_t mom) const {
    return fSDDKaonGaussWidth->Eval(mom);
  }
  Double_t GetSSDKaonMPV(Double_t mom) const {
    return fSSDKaonMPV->Eval(mom);
  }
  Double_t GetSSDKaonLandauWidth(Double_t mom) const {
    return fSSDKaonLandauWidth->Eval(mom);
  }
  Double_t GetSSDKaonGaussWidth(Double_t mom) const {
    return fSSDKaonGaussWidth->Eval(mom);
  }

  // proton getters
  Double_t GetSDDProtMPV(Double_t mom) const {
    return fSDDProtMPV->Eval(mom);
  }
  Double_t GetSDDProtLandauWidth(Double_t mom) const {
    return fSDDProtLandauWidth->Eval(mom);
  }
  Double_t GetSDDProtGaussWidth(Double_t mom) const {
    return fSDDProtGaussWidth->Eval(mom);
  }
  Double_t GetSSDProtMPV(Double_t mom) const {
    return fSSDProtMPV->Eval(mom);
  }
  Double_t GetSSDProtLandauWidth(Double_t mom) const {
    return fSSDProtLandauWidth->Eval(mom);
  }
  Double_t GetSSDProtGaussWidth(Double_t mom) const {
    return fSSDProtGaussWidth->Eval(mom);
  }

 private:

  AliITSPidParams(const AliITSPidParams& rec);
  AliITSPidParams& operator=(const AliITSPidParams &source);

  TFormula* fSDDPionMPV;         // pion dE/dx MPV vs. p in SDD
  TFormula* fSDDPionLandauWidth; // pion dE/dx Landau width vs. p in SDD
  TFormula* fSDDPionGaussWidth;  // pion dE/dx Gaussian width vs. p in SDD

  TFormula* fSSDPionMPV;         // pion dE/dx MPV vs. p in SSD
  TFormula* fSSDPionLandauWidth; // pion dE/dx Landau width vs. p in SSD
  TFormula* fSSDPionGaussWidth;  // pion dE/dx Gaussian width vs. p in SSD

  TFormula* fSDDKaonMPV;         // kaon dE/dx MPV vs. p in SDD
  TFormula* fSDDKaonLandauWidth; // kaon dE/dx Landau width vs. p in SDD
  TFormula* fSDDKaonGaussWidth;  // kaon dE/dx Gaussian width vs. p in SDD

  TFormula* fSSDKaonMPV;         // kaon dE/dx MPV vs. p in SSD
  TFormula* fSSDKaonLandauWidth; // kaon dE/dx Landau width vs. p in SSD
  TFormula* fSSDKaonGaussWidth;  // kaon dE/dx Gaussian width vs. p in SSD

  TFormula* fSDDProtMPV;         // kaon dE/dx MPV vs. p in SDD
  TFormula* fSDDProtLandauWidth; // kaon dE/dx Landau width vs. p in SDD
  TFormula* fSDDProtGaussWidth;  // kaon dE/dx Gaussian width vs. p in SDD

  TFormula* fSSDProtMPV;         // kaon dE/dx MPV vs. p in SSD
  TFormula* fSSDProtLandauWidth; // kaon dE/dx Landau width vs. p in SSD
  TFormula* fSSDProtGaussWidth;  // kaon dE/dx Gaussian width vs. p in SSD

  ClassDef(AliITSPidParams,1);
};
#endif
