#ifndef ALIEMCALTOWERRECPOINT_H
#define ALIEMCALTOWERRECPOINT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  RecPoint implementation for EMCAL-EMC 
//  An TowerRecPoint is a cluster of digits   
//           
//*-- Author: Dmitri Peressounko (RRC KI & SUBATECH)
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---

#include "TArrayI.h"
class TVector3 ;  

// --- Standard library ---

// --- AliRoot header files ---

#include "AliEMCALDigit.h"
#include "AliEMCALRecPoint.h"

class AliEMCALTowerRecPoint : public AliEMCALRecPoint  {

public:

  AliEMCALTowerRecPoint() ;
  AliEMCALTowerRecPoint(const char * opt) ;
  AliEMCALTowerRecPoint(const AliEMCALTowerRecPoint & rp):AliEMCALRecPoint(rp) {
    // cpy ctor requested by Coding Convention 
    // but not yet needed
    assert(0==1) ; 
  } 
 
  virtual ~AliEMCALTowerRecPoint() ;  

  virtual void  AddDigit(AliEMCALDigit & digit, Float_t Energy) ;          // add a digit to the digits list  
  Int_t       Compare(const TObject * obj) const;                         // method for sorting  

  virtual void  EvalAll(Float_t logWeight,TClonesArray * digits) ;
  virtual void  EvalGlobalPosition(Float_t logWeight, TClonesArray * digits) ;

  virtual void  ExecuteEvent(Int_t /*event*/, Int_t, Int_t) const; 

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
  virtual Int_t GetNumberOfLocalMax(AliEMCALDigit **  maxAt, Float_t * maxAtEnergy,
                                    Float_t locMaxCut,TClonesArray * digits ) const ; 
                                                                   // searches for the local maxima 
  Float_t     GetTime(void) const{return  fTime ; } 
  Bool_t      IsSortable() const {return kTRUE ; }                 // says that emcrecpoints are sortable objects 
  void        Print(Option_t * /*opt = "void"*/) ; 
  const TVector3 XYZInAlice(Float_t r = 9999., Float_t theta = 9999., Float_t phi = 9999.) const ;  

  AliEMCALTowerRecPoint & operator = (const AliEMCALTowerRecPoint & /*rvalue*/)  {
    // assignement operator requested by coding convention but not needed
    assert(0==1) ;
    return *this ; 
  }

 protected:
          void  EvalCoreEnergy(Float_t logWeight,TClonesArray * digits) ;             
  virtual void  EvalLocalPosition(Float_t /*logWeight*/,TClonesArray * /*digits*/) {;}// computes the position in the EMCAL module 
  virtual void  EvalDispersion(Float_t logWeight,TClonesArray * digits) ;   // computes the dispersion of the shower
  virtual void  EvalElipsAxis(Float_t logWeight, TClonesArray * digits );   // computes the axis of shower ellipsoide
          void  EvalTime( TClonesArray * digits );
  virtual Bool_t AreNeighbours(AliEMCALDigit * digit1, AliEMCALDigit * digit2 ) const ;

  Float_t fCoreEnergy ;       // energy in a shower core 
  Float_t fLambda[2] ;        // shower ellipse axes
  Float_t fDispersion ;       // shower dispersion
  Float_t *fEnergyList ;      //[fMulDigit] energy of digits
  Float_t fTime ;             // Time of the digit with maximal energy deposition
  
  ClassDef(AliEMCALTowerRecPoint,2)  // Tower RecPoint (cluster)

};

#endif // AliEMCALTOWERRECPOINT_H
