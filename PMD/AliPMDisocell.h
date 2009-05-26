#ifndef ALIPMDISOCELL_H
#define ALIPMDISOCELL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//  Date   : May 22, 2009                              //
//                                                     //
//  Store isolated single cell information             //
//  to be used in offline calibartion                  //
//                                                     //
//-----------------------------------------------------//
// Author -  B.K. Nandi
//
//#include "Rtypes.h"
#include "TObject.h"

class TClonesArray;

class AliPMDisocell : public TObject
{
 public:
  AliPMDisocell();
  AliPMDisocell( Int_t idet, Int_t ismn, Int_t irow, Int_t icol, Float_t cadc);
  AliPMDisocell(AliPMDisocell *pmdisocell);
  AliPMDisocell (const AliPMDisocell &pmdisocell);  //copy constructor
  AliPMDisocell &operator=(const AliPMDisocell &pmdisocell); //assignment op
  
  virtual ~AliPMDisocell();

  Int_t   GetDetector() const;
  Int_t   GetSmn() const;
  Int_t   GetRow() const;
  Int_t   GetCol() const;
  Float_t GetADC() const;

  
 protected:

  Int_t fDet;          // Pre = 0, CPV =1 plane 
  Int_t fSmn;
  Int_t fRow;
  Int_t fCol;
  Float_t fAdc;

  ClassDef(AliPMDisocell,0) // Keep isolated single cell information
};
#endif
