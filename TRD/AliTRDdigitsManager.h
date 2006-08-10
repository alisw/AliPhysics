#ifndef ALITRDDIGITSMANAGER_H
#define ALITRDDIGITSMANAGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDdigitsManager.h,v */

/////////////////////////////////////////////////////////////
//  Manages the TRD digits                                 //
/////////////////////////////////////////////////////////////

#include <TObject.h>

class TFile;
class TTree;

class AliTRDsegmentArray;
class AliTRDdataArrayI;
class AliTRDdigit;

class AliTRDdigitsManager : public TObject {

 public:

  enum { kNDict = 3 };

  AliTRDdigitsManager();
  AliTRDdigitsManager(const AliTRDdigitsManager &m);
  virtual ~AliTRDdigitsManager();
  AliTRDdigitsManager &operator=(const AliTRDdigitsManager &m);

  virtual void                CreateArrays();
  virtual void                ResetArrays();
  virtual void                Copy(TObject &m) const;

  virtual Bool_t              MakeBranch(TTree *tree);
          
  virtual Bool_t              ReadDigits(TTree *tree);
  virtual Bool_t              WriteDigits();

  virtual void                SetRaw();
  virtual void                SetEvent(Int_t evt)          { fEvent   = evt; };
  virtual void                SetSDigits(Int_t v = 1)      { fSDigits = v;   };

  virtual Bool_t              IsRaw() const                { return fIsRaw;         };
  static  Int_t               NDict()                      { return fgkNDict;       }; 

  virtual AliTRDsegmentArray *GetDigits() const            { return fDigits;        };
  virtual AliTRDsegmentArray *GetDictionary(Int_t i) const { return fDictionary[i]; };

          AliTRDdigit        *GetDigit(Int_t row, Int_t col, Int_t time, Int_t det) const;
          Int_t               GetTrack(Int_t track, Int_t row, Int_t col
                                     , Int_t time, Int_t det) const;

          AliTRDdataArrayI   *GetDigits(Int_t det) const;
          AliTRDdataArrayI   *GetDictionary(Int_t det, Int_t i) const;
          Int_t               GetTrack(Int_t track, AliTRDdigit *Digit) const;
          Short_t             GetDigitAmp(Int_t row, Int_t col, Int_t time, Int_t det) const;

 protected:

  static const Int_t  fgkNDict;            //  Number of track dictionary arrays

  Int_t               fEvent;              //  Event number

  TTree              *fTree;               //! Tree for the digits arrays

  AliTRDsegmentArray *fDigits;             //! Digits data array
  AliTRDsegmentArray *fDictionary[kNDict]; //! Track dictionary data array

  Bool_t              fIsRaw;              //  Flag indicating raw digits
  Bool_t              fSDigits;            //  Switch for the summable digits

  ClassDef(AliTRDdigitsManager,5)          //  Manages the TRD digits

};

#endif
