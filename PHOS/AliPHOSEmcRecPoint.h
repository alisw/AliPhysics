#ifndef ALIPHOSEMCRECPOINT_H
#define ALIPHOSEMCRECPOINT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  RecPoint implementation for PHOS-EMC 
//  An EmcRecPoint is a cluster of digits   
//           
//*-- Author: Dmitri Peressounko (RRC KI & SUBATECH)

// --- ROOT system ---

#include "TObject.h"
#include "TArrayI.h"
 
// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSDigit.h"
#include "AliPHOSRecPoint.h"

class AliPHOSEmcRecPoint : public AliPHOSRecPoint  {

public:

  AliPHOSEmcRecPoint(){
    // default ctor
  } ;                    
  AliPHOSEmcRecPoint(Float_t W0, Float_t LocMaxCut) ;
  virtual ~AliPHOSEmcRecPoint() ;  

  virtual void  AddDigit(AliPHOSDigit & digit, Float_t Energy) ;  // add a digit to the digits list  
  Int_t       Compare(TObject * obj) ;                         // method for sorting  
  
  Float_t     GetDelta (){ 
    // gets the fDelta data member 
    return fDelta ; }    
  Float_t     GetDispersion() ;                               // computes the dispersion of the shower
  void        GetElipsAxis(Float_t * lambda) ;                // computes the axis of shower ellipsoide
  Float_t *   GetEnergiesList(){
    // gets the list of energies makink this recpoint
    return fEnergyList ;} 
  virtual void ExecuteEvent(Int_t event, Int_t px, Int_t py) ; 
  Float_t     GetLocMaxCut () {
    // gets the cut of the local maximum search 
    return fLocMaxCut ; }
  Float_t     GetLogWeightCut (){
    // gets the logarythmic weight for the center of gravity calculation
    return fW0 ; }
  Float_t     GetMaximalEnergy(void) ;                        // get the highest energy in the cluster
  Int_t       GetMaximumMultiplicity() { 
    // gets the maximum number of digits allowed
    return   fMaxDigit ; } 
  Int_t       GetMultiplicity(void) const { 
    // gets the number of digits making this recpoint
    return fMulDigit ; } 
  Int_t       GetMultiplicityAtLevel(const Float_t level) ;   // computes multiplicity of digits with energy above relative level
  Int_t       GetNumberOfLocalMax(Int_t *  maxAt, Float_t * maxAtEnergy) ; // searches for the local maxima 
 
  Float_t     GetTotalEnergy(void) const { 
    // gets the total amplitude of this recpoint (in EMC RecPoint Amp = Energy)
    return fAmp ; }    
  void        GetLocalPosition(TVector3 &Lpos) ;  // computes the position in the PHOS module 
  Bool_t      IsEmc(void) {
    // true if the recpoint is in EMC
    return kTRUE ; } 
  Bool_t      IsSortable() const { 
    // says that emcrecpoints are sortable objects 
    return kTRUE ; } 
  void        Print(Option_t * opt = "void") ; 

private:

  Bool_t AreNeighbours(AliPHOSDigit * digit1, AliPHOSDigit * digit2 ) ;

  Float_t  fDelta ;          // parameter used to sort the clusters    
  Float_t  *fEnergyList ;    //[fMulDigit] energy of digits
  Float_t  fLocMaxCut ;      // minimum energy difference to distinguish two maxima 
  Float_t  fW0 ;             // logarithmic weight factor for center of gravity calculation
  
  ClassDef(AliPHOSEmcRecPoint,1)  // EMC RecPoint (cluster)

};

#endif // AliPHOSEMCRECPOINT_H
