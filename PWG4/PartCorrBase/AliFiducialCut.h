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
  virtual  ~AliFiducialCut() ;//virtual dtor
  
  void      InitParameters();
    
  Bool_t    CheckFiducialRegion(const TLorentzVector lv, 
                                const TArrayF* minphi, const TArrayF* maxphi, 
                                const TArrayF* mineta, const TArrayF* maxeta) const ;
  
  Bool_t    IsInFiducialCut    (const TLorentzVector lv, const TString det) const ;
  
  void      DoCTSFiducialCut  (Bool_t b) {fCTSFiducialCut   = b ; }
  void      DoEMCALFiducialCut(Bool_t b) {fEMCALFiducialCut = b ; }
  void      DoPHOSFiducialCut (Bool_t b) {fPHOSFiducialCut  = b ; }
  
  Bool_t    GetCTSFiducialCutStatus()  const {return fCTSFiducialCut   ; }
  Bool_t    GetEMCALFiducialCut()      const {return fEMCALFiducialCut ; }
  Bool_t    GetPHOSFiducialCutStatus() const {return fPHOSFiducialCut  ; }
  
  void      SetSimpleCTSFiducialCut  (const Float_t abseta, const Float_t phimin, const Float_t phimax) ;
  void      SetSimpleEMCALFiducialCut(const Float_t abseta, const Float_t phimin, const Float_t phimax) ;
  void      SetSimplePHOSFiducialCut (const Float_t abseta, const Float_t phimin, const Float_t phimax) ;
  
  void      Print(const Option_t * opt)const;

  
  void      AddCTSFidCutMaxEtaArray(Int_t size, Float_t* array)  
                                              { fCTSFidCutMaxEta->Set(size,array)   ; } 
  TArrayF * GetCTSFidCutMaxEtaArray() const   { return fCTSFidCutMaxEta             ; }
  
  void      AddCTSFidCutMaxPhiArray(Int_t size, Float_t* array)  
                                              { fCTSFidCutMaxPhi->Set(size,array)   ; }
  TArrayF * GetCTSFidCutMaxPhiArray() const   { return fCTSFidCutMaxPhi             ; }
  
  void      AddCTSFidCutMinEtaArray(Int_t size, Float_t* array)  
                                              { fCTSFidCutMinEta->Set(size,array)   ; } 
  TArrayF * GetCTSFidCutMinEtaArray() const   { return fCTSFidCutMinEta             ; }
  
  void      AddCTSFidCutMinPhiArray(Int_t size, Float_t* array)  
                                              { fCTSFidCutMinPhi->Set(size,array)   ; }
  TArrayF * GetCTSFidCutMinPhiArray() const   { return fCTSFidCutMinPhi             ; }
  
  void      AddEMCALFidCutMaxEtaArray(Int_t size, Float_t* array)  
                                              { fEMCALFidCutMaxEta->Set(size,array) ; } 
  TArrayF * GetEMCALFidCutMaxEtaArray() const { return fEMCALFidCutMaxEta           ; }
  
  void      AddEMCALFidCutMaxPhiArray(Int_t size, Float_t* array)  
                                              { fEMCALFidCutMaxPhi->Set(size,array) ; }
  TArrayF * GetEMCALFidCutMaxPhiArray() const { return fEMCALFidCutMaxPhi           ; }
  
  void      AddEMCALFidCutMinEtaArray(Int_t size, Float_t* array)  
                                              { fEMCALFidCutMinEta->Set(size,array) ; } 
  TArrayF * GetEMCALFidCutMinEtaArray() const { return fEMCALFidCutMinEta           ; }
  
  void      AddEMCALFidCutMinPhiArray(Int_t size, Float_t* array)  
                                              { fEMCALFidCutMinPhi->Set(size,array) ; }
  TArrayF * GetEMCALFidCutMinPhiArray() const { return fEMCALFidCutMinPhi           ; }
  
  void      AddPHOSFidCutMaxEtaArray(Int_t size, Float_t* array)  
                                              { fPHOSFidCutMaxEta->Set(size,array)  ; } 
  TArrayF * GetPHOSFidCutMaxEtaArray() const  { return fPHOSFidCutMaxEta            ; }
  
  void      AddPHOSFidCutMaxPhiArray(Int_t size, Float_t* array)  
                                              { fPHOSFidCutMaxPhi->Set(size,array)  ; }
  TArrayF * GetPHOSFidCutMaxPhiArray() const  { return fPHOSFidCutMaxPhi            ; }
  
  void      AddPHOSFidCutMinEtaArray(Int_t size, Float_t* array)  
                                              { fPHOSFidCutMinEta->Set(size,array)  ; } 
  TArrayF * GetPHOSFidCutMinEtaArray() const  { return fPHOSFidCutMinEta            ; }

  void      AddPHOSFidCutMinPhiArray(Int_t size, Float_t* array)  
                                              { fPHOSFidCutMinPhi->Set(size,array)  ; }
  TArrayF * GetPHOSFidCutMinPhiArray() const  { return fPHOSFidCutMinPhi            ; }
  
protected:
  
  //Detector acceptance cuts
  Bool_t    fEMCALFiducialCut ;  // Apply fiducial cuts to EMCAL clusters
  Bool_t    fPHOSFiducialCut ;   // Apply fiducial cuts to PHOS clusters
  Bool_t    fCTSFiducialCut ;    // Apply fiducial cuts to  CTS tracks
  
  TArrayF * fCTSFidCutMinEta ;   // Take particles in CTS with eta > fCTSFidCutMinEta
  TArrayF * fCTSFidCutMinPhi ;   // Take particles in CTS with phi > fCTSFidCutMinPhi
  TArrayF * fCTSFidCutMaxEta ;   // Take particles in CTS with eta < fCTSFidCutMaxEta
  TArrayF * fCTSFidCutMaxPhi ;   // Take particles in CTS with phi > fCTSFidCutMaxPhi
  
  TArrayF * fEMCALFidCutMinEta ; // Take particles in EMCAL with eta > fEMCALFidCutMinEta
  TArrayF * fEMCALFidCutMinPhi ; // Take particles in EMCAL with phi > fEMCALFidCutMinPhi
  TArrayF * fEMCALFidCutMaxEta ; // Take particles in EMCAL with eta < fEMCALFidCutMaxEta
  TArrayF * fEMCALFidCutMaxPhi ; // Take particles in EMCAL with phi > fEMCALFidCutMaxPhi
  
  TArrayF * fPHOSFidCutMinEta ;  // Take particles in PHOS with eta > fPHOSFidCutMinEta
  TArrayF * fPHOSFidCutMinPhi ;  // Take particles in PHOS with phi > fPHOSFidCutMinPhi
  TArrayF * fPHOSFidCutMaxEta ;  // Take particles in PHOS with eta < fPHOSFidCutMaxEta
  TArrayF * fPHOSFidCutMaxPhi ;  // Take particles in PHOS with phi > fPHOSFidCutMaxPhi
  
  AliFiducialCut(const AliFiducialCut & g) ;              // cpy ctor
  AliFiducialCut & operator = (const AliFiducialCut & g) ;// cpy assignment
  
  ClassDef(AliFiducialCut,1)
  
} ;


#endif //ALIFIDUCIALCUT_H



