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

//#include "TObject.h"
#include "TArrayI.h"
 
// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSDigit.h"
#include "AliPHOSRecPoint.h"

class AliPHOSEmcRecPoint : public AliPHOSRecPoint  {

public:

  AliPHOSEmcRecPoint() ;
  AliPHOSEmcRecPoint(const AliPHOSEmcRecPoint & rp) {
    // cpy ctor requested by Coding Convention 
    // but not yet needed
    assert(0==1) ; 
  } 
 
  virtual ~AliPHOSEmcRecPoint() ;  

  virtual void  AddDigit(AliPHOSDigit & digit, Float_t Energy) ;          // add a digit to the digits list  
  Int_t       Compare(const TObject * obj) const;                         // method for sorting  

  virtual void  EvalAll(Float_t logWeight,TClonesArray * digits) ;
          void  EvalCoreEnergy(TClonesArray * digits) ;             
  virtual void  EvalLocalPosition(Float_t logWeight,TClonesArray * digits) ;// computes the position in the PHOS module 
  virtual void  EvalDispersion(Float_t logWeight,TClonesArray * digits) ;   // computes the dispersion of the shower
  virtual void  EvalElipsAxis(Float_t logWeight, TClonesArray * digits );   // computes the axis of shower ellipsoide
  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py) ; 

  Float_t         GetCoreEnergy()const {return fCoreEnergy ;}
  virtual Float_t GetDispersion()const {return fDispersion ;}
  virtual void    GetElipsAxis(Float_t * lambda)const { lambda[0] = fLambda[0] ;
                                                        lambda[1] = fLambda[1] ; }
  Float_t *   GetEnergiesList() const {return fEnergyList ;}       // gets the list of energies making this recpoint
  Float_t     GetMaximalEnergy(void) const ;                       // get the highest energy in the cluster
  Int_t       GetMaximumMultiplicity() const {return fMaxDigit ;}  // gets the maximum number of digits allowed
  Int_t       GetMultiplicity(void) const { return fMulDigit ; }   // gets the number of digits making this recpoint
  Int_t       GetMultiplicityAtLevel(const Float_t level) const ;  // computes multiplicity of digits with 
                                                                   // energy above relative level
  virtual Int_t GetNumberOfLocalMax(Int_t *  maxAt, Float_t * maxAtEnergy,
                                    Float_t locMaxCut,TClonesArray * digits ) const ; 
                                                                   // searches for the local maxima 
  Bool_t      IsEmc(void) const { return kTRUE ; }                 // true if the recpoint is in EMC
  Bool_t      IsSortable() const {return kTRUE ; }                 // says that emcrecpoints are sortable objects 
  void        Print(Option_t * opt = "void") ; 

  AliPHOSEmcRecPoint & operator = (const AliPHOSEmcRecPoint & rvalue)  {
    // assignement operator requested by coding convention
    // but not needed
    assert(0==1) ;
    return *this ; 
  }

 protected:

  virtual Bool_t AreNeighbours(AliPHOSDigit * digit1, AliPHOSDigit * digit2 ) const ;

  Float_t fCoreEnergy ;
  Float_t fLambda[2] ;        //
  Float_t fDispersion ;
  Float_t *fEnergyList ;    //[fMulDigit] energy of digits
  
  ClassDef(AliPHOSEmcRecPoint,1)  // EMC RecPoint (cluster)

};

#endif // AliPHOSEMCRECPOINT_H
