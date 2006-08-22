#ifndef ALIPMDCELL_H
#define ALIPMDCELL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//  Date   : August 05 2003                            //
//                                                     //
//  Store cell/track info which is used to assign      //
//  the correct track number to a multiple hit cell    //
//                                                     //
//-----------------------------------------------------//
// Author - B.K. Nandi
//
#include "TObject.h"
class TClonesArray;

class AliPMDcell : public TObject
{
 public:
  AliPMDcell();
  AliPMDcell(Int_t trnumber, Int_t smnumber,
	      Int_t xpos, Int_t ypos, Float_t edep);
  AliPMDcell(AliPMDcell *pmdcell);
  AliPMDcell (const AliPMDcell &alipmdcell);  // copy constructor
  AliPMDcell &operator=(const AliPMDcell &alipmdcell); // assignment op

  virtual ~AliPMDcell();

  Int_t   GetTrackNumber() const;
  Int_t   GetSMNumber() const;
  Int_t   GetX() const;
  Int_t   GetY() const;
  Float_t GetEdep() const;
  
 protected:
  Int_t   fTrNumber;     // Track Number
  Int_t   fSMNumber;     // Serial Module Number
  Int_t   fXpos;         // x-position of the cell
  Int_t   fYpos;         // y-position of the cell
  Float_t fEdep;         // Energy deposition in a cell
  
  ClassDef(AliPMDcell,4) // To keep cell information
};

#endif
