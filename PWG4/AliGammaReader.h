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
// Base class for reading data: MonteCarlo, ESD or AOD, of PHOS EMCAL and 
// Central Barrel Tracking detectors.
// Not all MC particles/tracks/clusters are kept, some kinematical restrictions are done.
// Mother class of : AliGammaDataReader: Fills ESD data in 3 TClonesArrays (PHOS, EMCAL, CTS)
//                 : AliGammaMCReader: Fills Kinematics data in 3 TClonesArrays (PHOS, EMCAL, CTS)
//                 : AliGammaMCDataReader: Fills ESD data in 3 TClonesArrays (PHOS, EMCAL, CTS) 
//                             and MC data in other 3 TClonesArray
//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include "TObject.h" 

class TClonesArray ; 
class TFormula ;
class TParticle ; 
class Riostream ;
//--- AliRoot system ---

class AliStack ;
class AliESDEvent ; 
class AliLog ;

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

  enum PhotonStatusType {
    kPromptPhoton=2,
    kFragmentPhoton=3,
    kPi0DecayPhoton=4, 
    kEtaDecayPhoton=5, 
    kOtherDecayPhoton=6,
    kUnknown=7
  };

  enum Datatype {kData, kMC, kMCData};
  
  void InitParameters();

  Int_t GetDataType() const { return fDataType ; }
  void SetDataType(Int_t data ){fDataType = data ; }

  virtual Float_t  GetCTSEtaCut() const {return fCTSEtaCut ; }
  virtual Float_t  GetEMCALEtaCut() const {return fEMCALEtaCut ; }
  virtual Float_t  GetPHOSEtaCut() const {return fPHOSEtaCut ; }
  virtual Float_t  GetPhiEMCALCut(Int_t i) { return  fPhiEMCALCut[i]  ; }
  virtual Float_t  GetPhiPHOSCut(Int_t i) { return  fPhiPHOSCut[i]  ; }
  virtual Float_t  GetNeutralPtCut()    {  return fNeutralPtCut  ; }
  virtual Float_t  GetChargedPtCut()  {  return fChargedPtCut  ; }

  virtual Float_t  GetEMCALIPDistance()  const {  return fEMCALIPDistance ; }
  virtual Float_t  GetPHOSIPDistance()  const {  return fPHOSIPDistance ; }
  virtual Float_t  GetEMCALMinAngle()  const {  return fEMCALMinAngle ; }
  virtual Float_t  GetPHOSMinAngle()  const {  return fPHOSMinAngle ; }

  virtual Bool_t   IsEMCALPIDOn() const {return fEMCALPID ; }
  virtual Bool_t   IsPHOSPIDOn() const {return fPHOSPID ; }
  virtual Float_t  GetEMCALPhotonWeight() const  { return  fEMCALPhotonWeight  ; }
  virtual Float_t  GetEMCALPi0Weight() const     {  return fEMCALPi0Weight  ; }
  virtual Float_t  GetEMCALElectronWeight() const  { return  fEMCALElectronWeight  ; }
  virtual Float_t  GetEMCALChargeWeight() const     {  return fEMCALChargeWeight  ; }
  virtual Float_t  GetEMCALNeutralWeight() const     {  return fEMCALNeutralWeight  ; }
  virtual Float_t  GetPHOSPhotonWeight() const   {  return fPHOSPhotonWeight  ; }
  virtual Float_t  GetPHOSPi0Weight() const   {  return fPHOSPi0Weight  ; }
  virtual Float_t  GetPHOSElectronWeight() const   {  return fPHOSElectronWeight  ; }
  virtual Float_t  GetPHOSChargeWeight() const   {  return fPHOSChargeWeight  ; }
  virtual Float_t  GetPHOSNeutralWeight() const   {  return fPHOSNeutralWeight  ; }

  virtual Bool_t  IsPHOSPIDWeightFormulaOn() const   {  return fPHOSWeightFormula  ; } 
  virtual TFormula * GetPHOSPhotonWeightFormula() const     {  return fPHOSPhotonWeightFormula  ; } 
  virtual TFormula * GetPHOSPi0WeightFormula() const    {  return fPHOSPi0WeightFormula  ; }

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

  virtual Bool_t IsEMCALOn() const  { return fSwitchOnEMCAL ; }
  virtual Bool_t IsPHOSOn() const  { return fSwitchOnPHOS ; }
  virtual Bool_t IsCTSOn() const  { return fSwitchOnCTS ; }
  
  virtual void SwitchOnEMCAL(Bool_t sw) {fSwitchOnEMCAL = sw ; }
  virtual void SwitchOnPHOS(Bool_t sw) {fSwitchOnPHOS = sw ; }
  virtual void SwitchOnCTS(Bool_t sw) {fSwitchOnCTS = sw ; }

  virtual void CreateParticleList(TObject* , TObject * , 
				  TClonesArray * , TClonesArray * , TClonesArray *,  
				  TClonesArray * , TClonesArray * , TClonesArray * ) {;}
  
 protected:
  
  Int_t        fDataType ; //Select MC:Kinematics, Data:ESD/AOD, MCData:Both
  
  Bool_t     fSwitchOnEMCAL ; //Do not consider EMCAL neutral particles/clusters
  Bool_t     fSwitchOnPHOS ;//Do not consider PHOS neutral particles/clusters
  Bool_t     fSwitchOnCTS ;//Do not consider tracks/ charged particles

  //Kinematical cuts
  Float_t      fCTSEtaCut ;//CTS  pseudorapidity acceptance
  Float_t      fEMCALEtaCut ;//EMCAL pseudorapidity acceptance
  Float_t      fPHOSEtaCut ;//PHOS pseudorapidity acceptance
  Float_t      fPhiEMCALCut[2]; //EMCAL phi acceptance 
  Float_t      fPhiPHOSCut[2];  //PHOS phi acceptance
  Float_t      fNeutralPtCut; //pT Threshold on neutral particles
  Float_t      fChargedPtCut;  // pT  Threshold on charged particles

  //Overlapping 
  Float_t      fEMCALIPDistance; //Calorimeter IP distance.
  Float_t      fPHOSIPDistance; //Calorimeter IP distance
  Float_t      fEMCALMinAngle; //Gamma decay minimum aperture angle for overlapping.
  Float_t      fPHOSMinAngle; //Gamma decay minimum aperture angle for overlapping.

  //PID
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



