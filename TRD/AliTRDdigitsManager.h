#ifndef ALITRDDIGITSMANAGER_H
#define ALITRDDIGITSMANAGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDdigitsManager.h,v */

/////////////////////////////////////////////////////////////
//  Manages the TRD digits                                 //
/////////////////////////////////////////////////////////////

#include "TObject.h"

#include "AliTRDsegmentArray.h"
#include "AliTRDdataArrayI.h"
#include "AliTRDdigit.h"

const Int_t  kNDict = 3;

class AliTRDdigitsManager : public TObject {

 public:

  AliTRDdigitsManager();
  AliTRDdigitsManager(AliTRDdigitsManager &m);
  virtual ~AliTRDdigitsManager();

  virtual void                Copy(AliTRDdigitsManager &m);
  virtual Bool_t              MakeBranch();
  virtual Bool_t              ReadDigits();
  virtual Bool_t              WriteDigits();

  virtual void                SetRaw();

  virtual Bool_t              IsRaw()                { return fIsRaw;         };
  virtual AliTRDsegmentArray *GetDigits()            { return fDigits;        };
  virtual AliTRDsegmentArray *GetDictionary(Int_t i) { return fDictionary[i]; };

          AliTRDdigit        *GetDigit(Int_t row, Int_t col, Int_t time, Int_t det);
          Int_t               GetTrack(Int_t track, Int_t row, Int_t col, Int_t time, Int_t det);

  inline  AliTRDdataArrayI   *GetDigits(Int_t det);
  inline  AliTRDdataArrayI   *GetDictionary(Int_t det, Int_t i);
  inline  Int_t               GetTrack(Int_t track, AliTRDdigit *Digit);

  inline  AliTRDdigitsManager &operator=(AliTRDdigitsManager &m);

 protected:

  AliTRDsegmentArray *fDigits;             //! Digits data Array
  AliTRDsegmentArray *fDictionary[kNDict]; //! Track dictionary data array

  Bool_t              fIsRaw;              //  Flag indicating raw digits

  ClassDef(AliTRDdigitsManager,1)          //  Manages the TRD digits

};

//_____________________________________________________________________________
AliTRDdataArrayI *AliTRDdigitsManager::GetDigits(Int_t det) 
{
  //
  // Returns the digits array for one detector
  //

  return (AliTRDdataArrayI *) fDigits->At(det);

}

//_____________________________________________________________________________
AliTRDdataArrayI *AliTRDdigitsManager::GetDictionary(Int_t det, Int_t i) 
{
  //
  // Returns the dictionary for one detector
  //

  return (AliTRDdataArrayI *) fDictionary[i]->At(det);

}

//_____________________________________________________________________________
Int_t AliTRDdigitsManager::GetTrack(Int_t track, AliTRDdigit *Digit)
{
  // 
  // Returns the MC-track numbers from the dictionary for a given digit
  //

  Int_t row  = Digit->GetRow();
  Int_t col  = Digit->GetCol();
  Int_t time = Digit->GetTime();
  Int_t det  = Digit->GetDetector();

  return GetTrack(track,row,col,time,det);

}

//_____________________________________________________________________________
AliTRDdigitsManager &AliTRDdigitsManager::operator=(AliTRDdigitsManager &m)
{
  //
  // Assignment operator
  //

  if (this != &m) m.Copy(*this);
  return *this;

}

#endif
