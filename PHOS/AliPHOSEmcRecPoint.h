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
    fEnergyList = 0;  
  } ;                    
  AliPHOSEmcRecPoint(Float_t W0, Float_t LocMaxCut) ;
  AliPHOSEmcRecPoint(const AliPHOSEmcRecPoint & rp) {
    // cpy ctor requested by Coding Convention 
    // but not yet needed
    assert(0==1) ; 
  } 
 
  virtual ~AliPHOSEmcRecPoint() ;  

  virtual void  AddDigit(AliPHOSDigit & digit, Float_t Energy) ;  // add a digit to the digits list  
  Int_t       Compare(const TObject * obj) const;                         // method for sorting  
  Float_t     CoreEnergy() ;
  void        EvalAll() ;
  void        EvalLocalPosition() ;                                // computes the position in the PHOS module 
  Float_t     GetDelta () const {     return fDelta ; }    
  Float_t     GetDispersion() const ;                              // computes the dispersion of the shower
  void        GetElipsAxis(Float_t * lambda) ;                     // computes the axis of shower ellipsoide
  Float_t *   GetEnergiesList() const {    return fEnergyList ;}   // gets the list of energies making this recpoint
  virtual void ExecuteEvent(Int_t event, Int_t px, Int_t py) ; 
  Float_t     GetLocMaxCut ()const  {    return fLocMaxCut ; }     // gets the cut of the local maximum search  
  Float_t     GetLogWeightCut ()const { return fW0 ; }             // gets the logarythmic weight for the 
                                                                   // center of gravity calculation
  Float_t     GetMaximalEnergy(void) const ;                       // get the highest energy in the cluster
  Int_t       GetMaximumMultiplicity() const {return fMaxDigit ;}  // gets the maximum number of digits allowed
  Int_t       GetMultiplicity(void) const { return fMulDigit ; }   // gets the number of digits making this recpoint
  Int_t       GetMultiplicityAtLevel(const Float_t level) const ;  // computes multiplicity of digits with 
                                                                   // energy above relative level
  Int_t       GetNumberOfLocalMax(Int_t *  maxAt, Float_t * maxAtEnergy) const ; // searches for the local maxima 
 
  Bool_t      IsEmc(void) const { return kTRUE ; }                 // true if the recpoint is in EMC
  Bool_t      IsSortable() const {return kTRUE ; }                 // says that emcrecpoints are sortable objects 
  void        Print(Option_t * opt = "void") ; 

  AliPHOSEmcRecPoint & operator = (const AliPHOSEmcRecPoint & rvalue)  {
    // assignement operator requested by coding convention
    // but not needed
    assert(0==1) ;
    return *this ; 
  }

 private:

  Bool_t AreNeighbours(AliPHOSDigit * digit1, AliPHOSDigit * digit2 ) const ;

  Float_t  fDelta ;          // parameter used to sort the clusters    
  Float_t  *fEnergyList ;    //[fMulDigit] energy of digits
  Float_t  fLocMaxCut ;      // minimum energy difference to distinguish two maxima 
  Float_t  fW0 ;             // logarithmic weight factor for center of gravity calculation
  
  ClassDef(AliPHOSEmcRecPoint,1)  // EMC RecPoint (cluster)

};

#endif // AliPHOSEMCRECPOINT_H
