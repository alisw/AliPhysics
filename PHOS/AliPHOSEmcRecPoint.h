#ifndef ALIPHOSEMCRECPOINT_H
#define ALIPHOSEMCRECPOINT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////
//  Rec Point in the EM calorimeter of PHOS     //
//                                              //
//  Author Dmitri Peressounko RRC KI            //
//   comment: contains list of AliPHOSDigit's * //  
//     and evaluates a few average values       //
//////////////////////////////////////////////////

// --- ROOT system ---

#include "TObject.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSDigit.h"
#include "AliPHOSRecPoint.h"

class AliPHOSEmcRecPoint : public AliPHOSRecPoint  {

public:

  AliPHOSEmcRecPoint() ;                    
  AliPHOSEmcRecPoint(Float_t W0, Float_t LocMaxCut) ;
  //  virtual ~AliPHOSEmcRecPoint() ; 
  void        AddDigit(AliDigitNew & digit, Float_t Energy) ;  // add a digit to the digits list  
  Int_t       Compare(TObject * obj) ;                      // method for sorting  
  
  Float_t     GetDelta (){ return fDelta ; }    
  Float_t     GetDispersion() ;                             // computes the dispersion of the shower
  void        GetElipsAxis(Float_t * lambda) ;              // computes the axis of shower ellipsoide
  Float_t *   GetEnergiesList(){return fEnergyList ;} 
  Float_t     GetLocMaxCut () {return fLocMaxCut ; }
  Float_t     GetLogWeightCut (){return fW0 ; }
  Float_t     GetMaximalEnergy(void) ;                      // get the highest energy in the cluster
  Int_t       GetMaximumMultiplicity() { return   fMaxDigit ; } 
  Int_t       GetMultiplicity(void) const { return fMulDigit ; } 
  Int_t       GetMultiplicityAtLevel(const Float_t level) ; // computes multiplicity of digits with energy above relative level
  Int_t       GetNumberOfLocalMax(int *  maxAt, Float_t * maxAtEnergy) ;     // searches for the local maxima 
 
  Float_t     GetTotalEnergy(void) const { return fAmp ; } // in EMC RecPoint Amp = Energy
  void        GetLocalPosition(TVector3 &Lpos) ;            // computes the position in the PHOS module 
  Bool_t      IsEmc(void) {return kTRUE ; } 
  Bool_t IsSortable() const { return kTRUE ; } 
  void Print(Option_t * opt = "void") ; 

private:

  Bool_t AreNeighbours(AliPHOSDigit * digit1, AliPHOSDigit * digit2 ) ;

public:

  //  AliPHOSEmcRecPoint& operator = (AliPHOSEmcRecPoint clu) ;  
  
private:
  Float_t        fDelta ;        // parameter used to sort the clusters   
  Float_t        fLocMaxCut ;    // parameter used for local maximum searc
  Float_t    *   fEnergyList ;   //energy of digits
  Float_t        fW0 ;           // logarithmic weight factor for center of gravity calculation

public: 

ClassDef(AliPHOSEmcRecPoint,1)  // EMC cluster, version 1

};

#endif // AliPHOSEMCRECPOINT_H
