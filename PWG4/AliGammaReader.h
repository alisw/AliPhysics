#ifndef ALIGAMMAREADER_H
#define ALIGAMMAREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.2  2007/08/17 12:40:04  schutz
 * New analysis classes by Gustavo Conesa
 *
 * Revision 1.1.2.1  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
 *
 */

//_________________________________________________________________________
// Base class for reading data in order to do prompt gamma correlations 
//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include <TParticle.h> 
#include <TFormula.h>
#include <TClonesArray.h> 
#include "AliStack.h"
#include "TObject.h" 
 
class AliESD ; 

class TH2F ; 

class AliGammaReader : public TObject {

public: 

  AliGammaReader() ; // ctor
  AliGammaReader(const AliGammaReader & g) ; // cpy ctor
  AliGammaReader & operator = (const AliGammaReader & g) ;//cpy assignment
  virtual ~AliGammaReader() {;} //virtual dtor

  enum PidType {
    kPhoton = 22,
    kPi0 = 111,
    kEta = 221, 
    kElectron = 11, 
    kEleCon = 13, 
    kNeutralHadron = 2112, 
    kChargedHadron = 211, 
    kNeutralUnknown = 130, 
    kChargedUnknown=321
  };

  enum datatype_t {kData, kMC, kMCData};
  
  void InitParameters();

  Int_t GetDataType(){ return fDataType ; }
  void SetDataType(Int_t data ){fDataType = data ; }

  virtual Float_t   GetCTSEtaCut() const {return fCTSEtaCut ; }
  virtual Float_t   GetEMCALEtaCut() const {return fEMCALEtaCut ; }
  virtual Float_t   GetPHOSEtaCut() const {return fPHOSEtaCut ; }
  virtual Float_t  GetPhiEMCALCut(Int_t i) { return  fPhiEMCALCut[i]  ; }
  virtual Float_t  GetPhiPHOSCut(Int_t i) { return  fPhiPHOSCut[i]  ; }
  virtual Float_t  GetNeutralPtCut()    {  return fNeutralPtCut  ; }
  virtual Float_t  GetChargedPtCut()  {  return fChargedPtCut  ; }

  virtual Float_t  GetEMCALIPDistance()  {  return fEMCALIPDistance ; }
  virtual Float_t  GetPHOSIPDistance()  {  return fPHOSIPDistance ; }
  virtual Float_t  GetEMCALMinAngle()  {  return fEMCALMinAngle ; }
  virtual Float_t  GetPHOSMinAngle()  {  return fPHOSMinAngle ; }

  virtual Bool_t   IsEMCALPIDOn() const {return fEMCALPID ; }
  virtual Bool_t   IsPHOSPIDOn() const {return fPHOSPID ; }
  virtual Float_t  GetEMCALPhotonWeight() { return  fEMCALPhotonWeight  ; }
  virtual Float_t  GetEMCALPi0Weight()    {  return fEMCALPi0Weight  ; }
  virtual Float_t  GetEMCALElectronWeight() { return  fEMCALElectronWeight  ; }
  virtual Float_t  GetEMCALChargeWeight()    {  return fEMCALChargeWeight  ; }
  virtual Float_t  GetEMCALNeutralWeight()    {  return fEMCALNeutralWeight  ; }
  virtual Float_t  GetPHOSPhotonWeight()  {  return fPHOSPhotonWeight  ; }
  virtual Float_t  GetPHOSPi0Weight()  {  return fPHOSPi0Weight  ; }
  virtual Float_t  GetPHOSElectronWeight()  {  return fPHOSElectronWeight  ; }
  virtual Float_t  GetPHOSChargeWeight()  {  return fPHOSChargeWeight  ; }
  virtual Float_t  GetPHOSNeutralWeight()  {  return fPHOSNeutralWeight  ; }

  virtual Bool_t  IsPHOSPIDWeightFormulaOn()  {  return fPHOSWeightFormula  ; } 
  virtual TFormula * GetPHOSPhotonWeightFormula()    {  return fPHOSPhotonWeightFormula  ; } 
  virtual TFormula * GetPHOSPi0WeightFormula()   {  return fPHOSPi0WeightFormula  ; }

  virtual void Print(const Option_t * opt)const;
  
  virtual void SetCTSEtaCut(Float_t eta){ fCTSEtaCut= eta ; }
  virtual void SetEMCALEtaCut(Float_t eta){ fEMCALEtaCut= eta ; }
  virtual void SetPHOSEtaCut(Float_t eta){ fPHOSEtaCut= eta ; }
  virtual void SetPhiEMCALCut(Float_t  phi0, Float_t  phi1)
  { fPhiEMCALCut[0]= phi0 ; fPhiEMCALCut[1]= phi1 ;}
  virtual void SetPhiPHOSCut(Float_t  phi0, Float_t  phi1)
  { fPhiPHOSCut[0]= phi0 ; fPhiPHOSCut[1]= phi1 ;}
  virtual void SetNeutralPtCut(Float_t  pt){  fNeutralPtCut = pt ; }
  virtual void SetChargedPtCut(Float_t  pt){  fChargedPtCut = pt ; }

 virtual void SetEMCALIPDistance(Float_t  d){  fEMCALIPDistance = d ; }
  virtual void SetPHOSIPDistance(Float_t  d){  fPHOSIPDistance = d ; }
  virtual void SetEMCALMinAngle(Float_t  d){  fEMCALMinAngle = d ; }
  virtual void SetPHOSMinAngle(Float_t  d){  fPHOSMinAngle = d ; }
  
  virtual void SetEMCALPIDOn(Bool_t pid){ fEMCALPID= pid ; }
  virtual void SetPHOSPIDOn(Bool_t pid){ fPHOSPID= pid ; }
  virtual void SetEMCALPhotonWeight(Float_t  w){  fEMCALPhotonWeight = w ; }
  virtual void SetEMCALPi0Weight(Float_t  w){  fEMCALPi0Weight = w ; }
  virtual void SetEMCALElectronWeight(Float_t  w){  fEMCALElectronWeight = w ; }
  virtual void SetEMCALChargeWeight(Float_t  w){  fEMCALChargeWeight = w ; }
  virtual void SetEMCALNeutralWeight(Float_t  w){  fEMCALNeutralWeight = w ; }
  virtual void SetPHOSPhotonWeight(Float_t  w){  fPHOSPhotonWeight = w ; }
  virtual void SetPHOSPi0Weight(Float_t  w){  fPHOSPi0Weight = w ; }
  virtual void SetPHOSElectronWeight(Float_t  w){  fPHOSElectronWeight = w ; }
  virtual void SetPHOSChargeWeight(Float_t  w){  fPHOSChargeWeight = w ; }
  virtual void SetPHOSNeutralWeight(Float_t  w){  fPHOSNeutralWeight = w ; }

  virtual void UsePHOSPIDWeightFormula(Bool_t par)  { fPHOSWeightFormula  = par; } 
  virtual void SetPHOSPhotonWeightFormula(TFormula * photon)    {  fPHOSPhotonWeightFormula  = photon; } 
  virtual void SetPHOSPi0WeightFormula(TFormula * pi0)   {  fPHOSPi0WeightFormula  = pi0; }

  virtual void CreateParticleList(TObject* data, TObject * data2, 
				  TClonesArray * plCh, TClonesArray * plEMCAL, TClonesArray * plPHOS, 	  
				  TClonesArray * parton, TClonesArray * plPrimEMCAL, TClonesArray * plPrimPHOS) {;}
 protected:
  Int_t        fDataType ;
  Float_t      fCTSEtaCut ;//CTS  pseudorapidity acceptance
  Float_t      fEMCALEtaCut ;//EMCAL pseudorapidity acceptance
  Float_t      fPHOSEtaCut ;//PHOS pseudorapidity acceptance
  Float_t      fPhiEMCALCut[2]; //EMCAL phi acceptance 
  Float_t      fPhiPHOSCut[2];  //PHOS phi acceptance
  Float_t      fNeutralPtCut; //
  Float_t      fChargedPtCut;  // 

  Float_t      fEMCALIPDistance; //Calorimeter IP distance.
  Float_t      fPHOSIPDistance; //Calorimeter IP distance
  Float_t      fEMCALMinAngle; //Gamma decay minimum aperture angle for overlapping.
  Float_t      fPHOSMinAngle; //Gamma decay minimum aperture angle for overlapping.

  Bool_t       fEMCALPID ;//Fill EMCAL particle lists with particles with corresponding pid
  Bool_t       fPHOSPID;  //Fill PHOS particle lists with particles with corresponding pid
  Float_t      fEMCALPhotonWeight; //Bayesian PID weight for photons in EMCAL 
  Float_t      fEMCALPi0Weight;  //Bayesian PID weight for pi0 in EMCAL 
  Float_t      fEMCALElectronWeight; //Bayesian PID weight for electrons in EMCAL 
  Float_t      fEMCALChargeWeight;  //Bayesian PID weight for charged hadrons in EMCAL 
  Float_t      fEMCALNeutralWeight;  //Bayesian PID weight for neutral hadrons in EMCAL 
  Float_t      fPHOSPhotonWeight; //Bayesian PID weight for photons in PHOS 
  Float_t      fPHOSPi0Weight; //Bayesian PID weight for pi0 in PHOS 
  Float_t      fPHOSElectronWeight; //Bayesian PID weight for electrons in PHOS 
  Float_t      fPHOSChargeWeight; //Bayesian PID weight for charged hadrons in PHOS 
  Float_t      fPHOSNeutralWeight; //Bayesian PID weight for neutral hadrons in PHOS 

  Bool_t  fPHOSWeightFormula ; //Use parametrized weight threshold, function of energy
  TFormula * fPHOSPhotonWeightFormula ; //Formula for photon weight
  TFormula * fPHOSPi0WeightFormula ; //Formula for pi0 weight
  
  ClassDef(AliGammaReader,1)
} ;


#endif //ALIGAMMAREADER_H



