#ifndef ALIPHOSCPVRECPOINT_H
#define ALIPHOSCPVRECPOINT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  RecPoint implementation for PHOS-CPV
//  An CpvRecPoint is a cluster of digits   
//*-- Author: Yuri Kharlov
//  (after Dmitri Peressounko (RRC KI & SUBATECH))
//  30 October 2000 
// --- ROOT system ---

#include "TObject.h"
#include "TArrayI.h"
 
// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSDigit.h"
#include "AliPHOSRecPoint.h"

class AliPHOSCpvRecPoint : public AliPHOSRecPoint  {

public:

  AliPHOSCpvRecPoint(){
   // default ctor
    fEnergyList = 0;  
  } ;                    
  AliPHOSCpvRecPoint(Float_t W0, Float_t LocMaxCut) ;
  AliPHOSCpvRecPoint(const AliPHOSCpvRecPoint & rp) {
    // cpy ctor requested by Coding Convention 
    // but not yet needed
    assert(0==1) ; 
  } 
 
  virtual ~AliPHOSCpvRecPoint() ;  

  virtual void  AddDigit(AliPHOSDigit & digit, Float_t Energy) ;  // add a digit to the digits list  
  Int_t       Compare(const TObject * obj) const;                 // method for sorting  
  void        EvalAll( void ) ;
  void        EvalLocalPosition(void ) ;  // computes the position in the PHOS module 
  Float_t     GetDelta () const {     return fDelta ; }           // gets the fDelta data member 
  Float_t     GetDispersion() const ;                             // computes the dispersion of the shower
  void        GetElipsAxis(Float_t * lambda) const;               // computes the axis of shower ellipsoide
  Float_t *   GetEnergiesList()const {return fEnergyList ;}  // gets the list of energies making this recpoint
  virtual void ExecuteEvent(Int_t event, Int_t px, Int_t py) ; 
  Float_t     GetLocMaxCut () const {return fLocMaxCut ; }        // gets the cut of the local maximum search 
  Float_t     GetLogWeightCut ()const {return fW0 ; }             // gets the logarythmic weight for the 
                                                                  // center of gravity calculation
  Float_t     GetMaximalEnergy(void) const ;                      // get the highest energy in the cluster
  Int_t       GetMaximumMultiplicity() const {return fMaxDigit ;} // gets the maximum number of digits allowed
  Int_t       GetMultiplicity(void) const {return fMulDigit ; }   // gets the number of digits making this recpoint
  Int_t       GetMultiplicityAtLevel(const Float_t level) const;  // computes multiplicity of digits with energy 
                                                                  // above relative level
  Int_t       GetNumberOfLocalMax(Int_t *  maxAt, Float_t * maxAtEnergy) const ; // searches for the local maxima 
 
  void        GetClusterLengths(Int_t &lengX, Int_t &lengZ);      // cluster lengths along x and z
  Bool_t      IsEmc(void) const {return kFALSE ;   }              // tells that this is not a EMC
  Bool_t      IsCPV(void) const {return (fPHOSMod <= ((AliPHOSGeometry*) fGeom)->GetNCPVModules()) ; }     
                                                                  // true if the recpoint is in CPV
  Bool_t      IsSortable() const {  return kTRUE ; }              // says that emcrecpoints are sortable objects
  void        Print(Option_t * opt = "void") ; 

  AliPHOSCpvRecPoint & operator = (const AliPHOSCpvRecPoint & rvalue)  {
    // assignement operator requested by coding convention but not needed
    assert(0==1) ;
    return *this ; 
  }

 private:

  Bool_t AreNeighbours(AliPHOSDigit * digit1, AliPHOSDigit * digit2 ) const ;

  Float_t  fDelta ;          // parameter used to sort the clusters    
  Float_t  *fEnergyList ;    //[fMulDigit] energy of digits
  Float_t  fLocMaxCut ;      // minimum energy difference to distinguish two maxima 
  Float_t  fW0 ;             // logarithmic weight factor for center of gravity calculation
  Int_t    fLengX ;          // cluster length along x
  Int_t    fLengZ ;          // cluster length along z
  
  ClassDef(AliPHOSCpvRecPoint,1)  // CPV RecPoint (cluster)

};

#endif // AliPHOSCPVRECPOINT_H
