#ifndef ALIPHOSPPSDRECPOINT_H
#define ALIPHOSPPSDRECPOINT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Cluster in the PPSD of PHOS               //
//  Version SUBATECH                          //
//  Author Yves Schutz     SUBATECH           //
//      comment: its a list of AliPHOSDigit's //  
//                                            //
////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSDigit.h"
#include "AliPHOSRecPoint.h"

class AliPHOSPpsdRecPoint : public AliPHOSRecPoint {

public:

  AliPHOSPpsdRecPoint() ;          // ctor   
  virtual ~AliPHOSPpsdRecPoint() ; // dtor
 void AddDigit(AliDigitNew & digit, Float_t Energy) ;
  Int_t   Compare(TObject * obj) ;                    // method to sort clusters

  Float_t GetDelta(void) {return fDelta ;}
  Int_t   GetMultiplicity(void) const { return fMulDigit ; } 
  Int_t   GetMaximumMultiplicity() { return   fMaxDigit ; } 
  void    GetLocalPosition(TVector3 &LPos) ; // computes the position in the module of the cluster center 
  Float_t GetTotalEnergy(void) const { return fAmp ; }    // in Ppsd EMC RecPoint Amp = Energy                                        //projection of ALICE axes on PHOS Module, y = 0 .
  Bool_t  GetUp() ;               // true if cluster is in upper ppsd 
  Bool_t  IsEmc(void) {return kFALSE ; } 
  Bool_t  IsSortable() const { return kTRUE ; }
  virtual void  Paint(Option_t * option="");
  void    Print(Option_t * opt = "void") ; 

  //  AliPHOSPpsdRecPoint&  operator = (AliPHOSPpsdRecPoint Clu) ;

private:

  Float_t        fDelta ;         // parameter used for sorting


public: 
  
  ClassDef(AliPHOSPpsdRecPoint,1)  // PPSD RecPoint, version 1

};

#endif // AliPHOSPPSDRECPOINT_H
