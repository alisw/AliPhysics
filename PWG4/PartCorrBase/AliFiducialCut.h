#ifndef ALIFIDUCIALCUT_H
#define ALIFIDUCIALCUT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//_________________________________________________________________________
// Class for track/cluster acceptance selection
// Selection in Central barrel, EMCAL and PHOS
//  
// Several selection regions possible for the different
// detectors
//
//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include <TObject.h> 
#include <TArrayF.h> 

class TString ;
class TLorentzVector ;

//--- AliRoot system ---


class AliFiducialCut : public TObject {

 public: 
  AliFiducialCut() ; // ctor
  virtual ~AliFiducialCut() ;//virtual dtor

 private:
  AliFiducialCut(const AliFiducialCut & g) ; // cpy ctor
  AliFiducialCut & operator = (const AliFiducialCut & g) ;//cpy assignment
  
 public:

  void InitParameters();
  
  void Print(const Option_t * opt)const;
  Bool_t CheckFiducialRegion(const TLorentzVector lv, const TArrayF* minphi, const TArrayF* maxphi, const TArrayF* mineta, const TArrayF* maxeta) const ;
  Bool_t IsInFiducialCut    (const TLorentzVector lv, const TString det) const ;
  
  void DoCTSFiducialCut(Bool_t b) {fCTSFiducialCut = b; }
  void DoEMCALFiducialCut(Bool_t b) {fEMCALFiducialCut = b; }
  void DoPHOSFiducialCut(Bool_t b) {fPHOSFiducialCut = b; }

  Bool_t GetCTSFiducialCutStatus() const {return fCTSFiducialCut ; }
  Bool_t GetEMCALFiducialCut() const {return fEMCALFiducialCut ; }
  Bool_t GetPHOSFiducialCutStatus() const {return fPHOSFiducialCut ; }

  void SetSimpleCTSFiducialCut(const Float_t abseta, const Float_t phimin, const Float_t phimax) ;
  void SetSimpleEMCALFiducialCut(const Float_t abseta, const Float_t phimin, const Float_t phimax) ;
  void SetSimplePHOSFiducialCut(const Float_t abseta, const Float_t phimin, const Float_t phimax) ;

  void AddCTSFidCutMaxEtaArray(TArrayF & array)  
  {  fCTSFidCutMaxEta  = new TArrayF(array) ; } 
  TArrayF * GetCTSFidCutMaxEtaArray() const   {return   fCTSFidCutMaxEta;}
  
  void AddCTSFidCutMaxPhiArray(TArrayF & array)  
  {  fCTSFidCutMaxPhi  = new TArrayF(array) ; }
  TArrayF * GetCTSFidCutMaxPhiArray() const   {return   fCTSFidCutMaxPhi;}
  
  void AddCTSFidCutMinEtaArray(TArrayF & array)  
  {  fCTSFidCutMinEta  = new TArrayF(array) ; } 
  TArrayF * GetCTSFidCutMinEtaArray() const   {return   fCTSFidCutMinEta;}
 
  void AddCTSFidCutMinPhiArray(TArrayF & array)  
  {  fCTSFidCutMinPhi  = new TArrayF(array) ; }
  TArrayF * GetCTSFidCutMinPhiArray() const   {return   fCTSFidCutMinPhi;}
  
  void AddEMCALFidCutMaxEtaArray(TArrayF & array)  
  {  fEMCALFidCutMaxEta  = new TArrayF(array) ; } 
  TArrayF * GetEMCALFidCutMaxEtaArray() const   {return   fEMCALFidCutMaxEta;}
  
  void AddEMCALFidCutMaxPhiArray(TArrayF & array)  
  {  fEMCALFidCutMaxPhi  = new TArrayF(array) ; }
  TArrayF * GetEMCALFidCutMaxPhiArray() const   {return   fEMCALFidCutMaxPhi;}
 
  void AddEMCALFidCutMinEtaArray(TArrayF & array)  
  {  fEMCALFidCutMinEta  = new TArrayF(array) ; } 
  TArrayF * GetEMCALFidCutMinEtaArray() const   {return   fEMCALFidCutMinEta;}
  
  void AddEMCALFidCutMinPhiArray(TArrayF & array)  
  {  fEMCALFidCutMinPhi  = new TArrayF(array) ; }
  TArrayF * GetEMCALFidCutMinPhiArray() const   {return   fEMCALFidCutMinPhi;}
  
  void AddPHOSFidCutMaxEtaArray(TArrayF & array)  
  {  fPHOSFidCutMaxEta  = new TArrayF(array) ; } 
  TArrayF * GetPHOSFidCutMaxEtaArray() const   {return   fPHOSFidCutMaxEta;}
  
   void AddPHOSFidCutMaxPhiArray(TArrayF & array)  
   {  fPHOSFidCutMaxPhi  = new TArrayF(array) ; }
   TArrayF * GetPHOSFidCutMaxPhiArray() const   {return   fPHOSFidCutMaxPhi;}
   
   void AddPHOSFidCutMinEtaArray(TArrayF & array)  
   {  fPHOSFidCutMinEta  = new TArrayF(array) ; } 
   TArrayF * GetPHOSFidCutMinEtaArray() const   {return   fPHOSFidCutMinEta;}
 
   void AddPHOSFidCutMinPhiArray(TArrayF & array)  
   {  fPHOSFidCutMinPhi  = new TArrayF(array) ; }
   TArrayF * GetPHOSFidCutMinPhiArray() const   {return   fPHOSFidCutMinPhi;}
   
 
   
 protected:
   
   //Detector acceptance cuts
   Bool_t     fEMCALFiducialCut ; // Apply fiducial cuts to EMCAL clusters
   Bool_t     fPHOSFiducialCut ;// Apply fiducial cuts to PHOS clusters
   Bool_t     fCTSFiducialCut ;//Apply fiducial cuts to  CTS tracks
   
   TArrayF * fCTSFidCutMinEta ; //Take particles in CTS with eta > fCTSFidCutMinEta
   TArrayF * fCTSFidCutMinPhi ; //Take particles in CTS with phi > fCTSFidCutMinPhi
   TArrayF * fCTSFidCutMaxEta ; //Take particles in CTS with eta < fCTSFidCutMaxEta
   TArrayF * fCTSFidCutMaxPhi ; //Take particles in CTS with phi > fCTSFidCutMaxPhi
   
   TArrayF * fEMCALFidCutMinEta ; //Take particles in EMCAL with eta > fEMCALFidCutMinEta
   TArrayF * fEMCALFidCutMinPhi ; //Take particles in EMCAL with phi > fEMCALFidCutMinPhi
   TArrayF * fEMCALFidCutMaxEta ; //Take particles in EMCAL with eta < fEMCALFidCutMaxEta
   TArrayF * fEMCALFidCutMaxPhi ; //Take particles in EMCAL with phi > fEMCALFidCutMaxPhi
   
   TArrayF * fPHOSFidCutMinEta ; //Take particles in PHOS with eta > fPHOSFidCutMinEta
   TArrayF * fPHOSFidCutMinPhi ; //Take particles in PHOS with phi > fPHOSFidCutMinPhi
   TArrayF * fPHOSFidCutMaxEta ; //Take particles in PHOS with eta < fPHOSFidCutMaxEta
   TArrayF * fPHOSFidCutMaxPhi ; //Take particles in PHOS with phi > fPHOSFidCutMaxPhi
  
  ClassDef(AliFiducialCut,1)
    } ;


#endif //ALIFIDUCIALCUT_H



