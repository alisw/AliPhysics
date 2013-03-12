#ifndef ALIITSPIDPARAMS_H
#define ALIITSPIDPARAMS_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to store parameters of ITS response functions           //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include <TFormula.h>
#include <TNamed.h>
#include "AliPID.h"

class AliITSPidParams : public TNamed {

 public:
  AliITSPidParams(Bool_t isMC = kFALSE);
  AliITSPidParams(Char_t * name, Bool_t isMC = kFALSE);
  ~AliITSPidParams();
  
  void InitMC();
  void InitData();
  Double_t GetLandauGausNormPdgCode(Double_t dedx, Int_t pdgCode, Double_t mom, Int_t lay) const;
  Double_t GetLandauGausNorm(Double_t dedx, Int_t partType, Double_t mom, Int_t lay) const;
  
  // electron setters
  void SetSDDElecMPVparams(const Double_t* params){
    for(Int_t ipar=0; ipar<5; ipar++) fSDDElecMPVBetheParams[ipar]=params[ipar];
  }  
  void SetSSDElecMPVparams(const Double_t* params){
    for(Int_t ipar=0; ipar<5; ipar++) fSSDElecMPVBetheParams[ipar]=params[ipar];
  }
  void SetSDDElecLandauWidth(const TFormula* form){
    if(fSDDElecLandauWidth) delete fSDDElecLandauWidth;
    fSDDElecLandauWidth=new TFormula(*form);
  }
  void SetSDDElecGaussWidth(const TFormula* form){
    if(fSDDElecGaussWidth) delete fSDDElecGaussWidth;
    fSDDElecGaussWidth=new TFormula(*form);
  } 
  void SetSSDElecLandauWidth(const TFormula* form){
    if(fSSDElecLandauWidth) delete fSSDElecLandauWidth;
    fSSDElecLandauWidth=new TFormula(*form);
  }
  void SetSSDElecGaussWidth(const TFormula* form){
    if(fSSDElecGaussWidth) delete fSSDElecGaussWidth;
    fSSDElecGaussWidth=new TFormula(*form);
  }

  //Hadron setters
  void SetSDDHadronMPVparams(const Double_t* params){
    for(Int_t ipar=0; ipar<5; ipar++) fSDDHadronMPVBetheParams[ipar]=params[ipar];
  }  
  void SetSSDHadronMPVparams(const Double_t* params){
    for(Int_t ipar=0; ipar<5; ipar++) fSSDHadronMPVBetheParams[ipar]=params[ipar];
  }
  
  //Pion setters
  void SetSDDPionLandauWidth(const TFormula* form){
    if(fSDDPionLandauWidth) delete fSDDPionLandauWidth;
    fSDDPionLandauWidth=new TFormula(*form);
  }
  void SetSDDPionGaussWidth(const TFormula* form){
    if(fSDDPionGaussWidth) delete fSDDPionGaussWidth;
    fSDDPionGaussWidth=new TFormula(*form);
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
  void SetSDDKaonLandauWidth(const TFormula* form){
    if(fSDDKaonLandauWidth) delete fSDDKaonLandauWidth;
    fSDDKaonLandauWidth=new TFormula(*form);
  }
  void SetSDDKaonGaussWidth(const TFormula* form){
    if(fSDDKaonGaussWidth) delete fSDDKaonGaussWidth;
    fSDDKaonGaussWidth=new TFormula(*form);
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
  void SetSDDProtLandauWidth(const TFormula* form){
    if(fSDDProtLandauWidth) delete fSDDProtLandauWidth;
    fSDDProtLandauWidth=new TFormula(*form);
  }
  void SetSDDProtGaussWidth(const TFormula* form){
    if(fSDDProtGaussWidth) delete fSDDProtGaussWidth;
    fSDDProtGaussWidth=new TFormula(*form);
  }
  void SetSSDProtLandauWidth(const TFormula* form){
    if(fSSDProtLandauWidth) delete fSSDProtLandauWidth;
    fSSDProtLandauWidth=new TFormula(*form);
  }
  void SetSSDProtGaussWidth(const TFormula* form){
    if(fSSDProtGaussWidth) delete fSSDProtGaussWidth;
    fSSDProtGaussWidth=new TFormula(*form);
  }

  // electron getters
  Double_t GetSDDElecMPV(Double_t mom) const {
    Double_t mass = AliPID::ParticleMass(AliPID::kElectron);
    return BetheBloch(mom,mass,fSDDElecMPVBetheParams);
  }
  Double_t GetSDDElecLandauWidth(Double_t mom) const {
    return fSDDElecLandauWidth->Eval(mom);
  }
  Double_t GetSDDElecGaussWidth(Double_t mom) const {
    return fSDDElecGaussWidth->Eval(mom);
  }
  Double_t GetSSDElecMPV(Double_t mom) const {
    Double_t mass = AliPID::ParticleMass(AliPID::kElectron);
    return BetheBloch(mom,mass,fSSDElecMPVBetheParams);
  }
  Double_t GetSSDElecLandauWidth(Double_t mom) const {
    return fSSDElecLandauWidth->Eval(mom);
  }
  Double_t GetSSDElecGaussWidth(Double_t mom) const {
    return fSSDElecGaussWidth->Eval(mom);
  }

  // pion getters
  Double_t GetSDDPionMPV(Double_t mom) const {
    Double_t mass = AliPID::ParticleMass(AliPID::kPion);
    return BetheBloch(mom,mass,fSDDHadronMPVBetheParams);
  }
  Double_t GetSDDPionLandauWidth(Double_t mom) const {
    return fSDDPionLandauWidth->Eval(mom);
  }
  Double_t GetSDDPionGaussWidth(Double_t mom) const {
    return fSDDPionGaussWidth->Eval(mom);
  }
  Double_t GetSSDPionMPV(Double_t mom) const {
    Double_t mass = AliPID::ParticleMass(AliPID::kPion);
    return BetheBloch(mom,mass,fSSDHadronMPVBetheParams);
  }
  Double_t GetSSDPionLandauWidth(Double_t mom) const {
    return fSSDPionLandauWidth->Eval(mom);
  }
  Double_t GetSSDPionGaussWidth(Double_t mom) const {
    return fSSDPionGaussWidth->Eval(mom);
  }

  // kaon getters
  Double_t GetSDDKaonMPV(Double_t mom) const {
    Double_t mass = AliPID::ParticleMass(AliPID::kKaon);
    return BetheBloch(mom,mass,fSDDHadronMPVBetheParams);
  }
  Double_t GetSDDKaonLandauWidth(Double_t mom) const {
    Double_t xmin = 0.2;
    mom = (mom > xmin) ? mom : xmin;
    return fSDDKaonLandauWidth->Eval(mom);
  }
  Double_t GetSDDKaonGaussWidth(Double_t mom) const {
    Double_t xmin = 0.2;
    mom = (mom > xmin) ? mom : xmin;
    return fSDDKaonGaussWidth->Eval(mom);
  }
  Double_t GetSSDKaonMPV(Double_t mom) const {
    Double_t mass = AliPID::ParticleMass(AliPID::kKaon);
    return BetheBloch(mom,mass,fSSDHadronMPVBetheParams);
  }
  Double_t GetSSDKaonLandauWidth(Double_t mom) const {
    Double_t xmin = 0.2;
    mom = (mom > xmin) ? mom : xmin;
    return fSSDKaonLandauWidth->Eval(mom);
  }
  Double_t GetSSDKaonGaussWidth(Double_t mom) const {
    Double_t xmin = 0.2;
    mom = (mom > xmin) ? mom : xmin;
    return fSSDKaonGaussWidth->Eval(mom);
  }

  // proton getters
  Double_t GetSDDProtMPV(Double_t mom) const {
    Double_t mass = AliPID::ParticleMass(AliPID::kProton);
    return BetheBloch(mom,mass,fSDDHadronMPVBetheParams);
  }
  Double_t GetSDDProtLandauWidth(Double_t mom) const {
    Double_t xmin = 0.4;
    mom = (mom > xmin) ? mom : xmin;
    return fSDDProtLandauWidth->Eval(mom);
  }
  Double_t GetSDDProtGaussWidth(Double_t mom) const {
    Double_t xmin = 0.4;
    mom = (mom > xmin) ? mom : xmin;
    return fSDDProtGaussWidth->Eval(mom);
  }
  Double_t GetSSDProtMPV(Double_t mom) const {
    Double_t mass = AliPID::ParticleMass(AliPID::kProton);
    return BetheBloch(mom,mass,fSSDHadronMPVBetheParams);
  }
  Double_t GetSSDProtLandauWidth(Double_t mom) const {
    Double_t xmin = 0.4;
    mom = (mom > xmin) ? mom : xmin;
    return fSSDProtLandauWidth->Eval(mom);
  }
  Double_t GetSSDProtGaussWidth(Double_t mom) const {
    Double_t xmin = 0.4;
    mom = (mom > xmin) ? mom : xmin;
    return fSSDProtGaussWidth->Eval(mom);
  }

 private:
  
  AliITSPidParams(const AliITSPidParams& rec);
  AliITSPidParams& operator=(const AliITSPidParams &source);

  Double_t BetheBloch(Double_t mom, Double_t mass, const Double_t* p) const;
  Double_t ExtrapolateWidth(Double_t mom, Double_t x1, Double_t y1, Double_t x2, Double_t y2) const;
  
  // Electron parameterizations
  Double_t fSDDElecMPVBetheParams[5];         // Electron bethe block parameters in SDD
  Double_t fSSDElecMPVBetheParams[5];         // Electron bethe block parameters in SSD
 
  TFormula* fSDDElecLandauWidth; // Electron dE/dx Landau width vs. p in SDD
  TFormula* fSDDElecGaussWidth;  // Electron dE/dx Gaussian width vs. p in SDD

  TFormula* fSSDElecLandauWidth; // Electron dE/dx Landau width vs. p in SSD
  TFormula* fSSDElecGaussWidth;  // Electron dE/dx Gaussian width vs. p in SSD
  
  // Hadron parameterizations
  Double_t fSDDHadronMPVBetheParams[5];         // hadrons bethe block parameters in SDD
  Double_t fSSDHadronMPVBetheParams[5];         // hadrons bethe block parameters in SSD
 
  TFormula* fSDDPionLandauWidth; // pion dE/dx Landau width vs. p in SDD
  TFormula* fSDDPionGaussWidth;  // pion dE/dx Gaussian width vs. p in SDD

  TFormula* fSSDPionLandauWidth; // pion dE/dx Landau width vs. p in SSD
  TFormula* fSSDPionGaussWidth;  // pion dE/dx Gaussian width vs. p in SSD

  TFormula* fSDDKaonLandauWidth; // kaon dE/dx Landau width vs. p in SDD
  TFormula* fSDDKaonGaussWidth;  // kaon dE/dx Gaussian width vs. p in SDD

  TFormula* fSSDKaonLandauWidth; // kaon dE/dx Landau width vs. p in SSD
  TFormula* fSSDKaonGaussWidth;  // kaon dE/dx Gaussian width vs. p in SSD

  TFormula* fSDDProtLandauWidth; // Proton dE/dx Landau width vs. p in SDD
  TFormula* fSDDProtGaussWidth;  // Proton dE/dx Gaussian width vs. p in SDD

  TFormula* fSSDProtLandauWidth; // Proton dE/dx Landau width vs. p in SSD
  TFormula* fSSDProtGaussWidth;  // Proton dE/dx Gaussian width vs. p in SSD


  
  ClassDef(AliITSPidParams,2);
};
#endif

