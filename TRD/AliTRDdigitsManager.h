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
  AliTRDdigitsManager(const AliTRDdigitsManager &m);
  virtual ~AliTRDdigitsManager();
  AliTRDdigitsManager &operator=(const AliTRDdigitsManager &m);

  virtual void                Copy(TObject &m);
  virtual Bool_t              MakeBranch();
  virtual Bool_t              ReadDigits();
  virtual Bool_t              WriteDigits();

  virtual void                SetRaw();

  virtual Bool_t              IsRaw()                { return fIsRaw;         };
  virtual AliTRDsegmentArray *GetDigits()            { return fDigits;        };
  virtual AliTRDsegmentArray *GetDictionary(Int_t i) { return fDictionary[i]; };

          AliTRDdigit        *GetDigit(Int_t row, Int_t col, Int_t time, Int_t det);
          Int_t               GetTrack(Int_t track, Int_t row, Int_t col, Int_t time, Int_t det);

          AliTRDdataArrayI   *GetDigits(Int_t det);
          AliTRDdataArrayI   *GetDictionary(Int_t det, Int_t i);
          Int_t               GetTrack(Int_t track, AliTRDdigit *Digit);

 protected:

  AliTRDsegmentArray *fDigits;             //! Digits data Array
  AliTRDsegmentArray *fDictionary[kNDict]; //! Track dictionary data array

  Bool_t              fIsRaw;              //  Flag indicating raw digits

  ClassDef(AliTRDdigitsManager,1)          //  Manages the TRD digits

};

#endif
