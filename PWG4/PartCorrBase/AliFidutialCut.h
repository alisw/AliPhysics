#ifndef ALIFIDUTIALCUT_H
#define ALIFIDUTIALCUT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//_________________________________________________________________________
// Class for track/cluster acceptance selection
// Selection in Central barrel, EMCAL and PHOS
//
//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
#include <TObject.h> 
#include <TArrayF.h> 

class TString ;
class TLorentzVector ;

//--- AliRoot system ---


class AliFidutialCut : public TObject {

public: 

  AliFidutialCut() ; // ctor
  AliFidutialCut(const AliFidutialCut & g) ; // cpy ctor
  AliFidutialCut & operator = (const AliFidutialCut & g) ;//cpy assignment
  virtual ~AliFidutialCut() ;//virtual dtor

  void InitParameters();
  
  void Print(const Option_t * opt)const;
  
  Bool_t IsInFidutialCut(TLorentzVector l, TString det) const ;
  
  void DoCTSFidutialCut(Bool_t b) {fCTSFidutialCut = b; }
  void DoEMCALFidutialCut(Bool_t b) {fEMCALFidutialCut = b; }
  void DoPHOSFidutialCut(Bool_t b) {fPHOSFidutialCut = b; }

  Bool_t GetCTSFidutialCutStatus() const {return fCTSFidutialCut ; }
  Bool_t GetEMCALFidutialCut() const {return fEMCALFidutialCut ; }
  Bool_t GetPHOSFidutialCutStatus() const {return fPHOSFidutialCut ; }

  void SetSimpleCTSFidutialCut(const Float_t abseta, const Float_t phimin, const Float_t phimax) ;
  void SetSimpleEMCALFidutialCut(const Float_t abseta, const Float_t phimin, const Float_t phimax) ;
  void SetSimplePHOSFidutialCut(const Float_t abseta, const Float_t phimin, const Float_t phimax) ;

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
   Bool_t     fEMCALFidutialCut ; // Apply fidutial cuts to EMCAL clusters
   Bool_t     fPHOSFidutialCut ;// Apply fidutial cuts to PHOS clusters
   Bool_t     fCTSFidutialCut ;//Apply fidutial cuts to  CTS tracks
   
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
  
  ClassDef(AliFidutialCut,1)
    } ;


#endif //ALIFIDUTIALCUT_H



