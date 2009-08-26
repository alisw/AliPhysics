#ifndef ALITRDDIGITSMANAGER_H
#define ALITRDDIGITSMANAGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////////
//  Manages the TRD digits                                 //
/////////////////////////////////////////////////////////////

#include <TObject.h>

class TFile;
class TTree;
class TBranch;  
class AliTRDdigit;
class AliTRDSignalIndex;
class AliTRDarrayADC;  
class AliTRDarraySignal; 
class AliTRDarrayDictionary;
class AliTRDdigitsParam;

class AliTRDdigitsManager : public TObject {

 public:

  enum { kNDict = 3 };

  AliTRDdigitsManager(Bool_t rawRec = kFALSE);  //if true digitsmanager uses only one entry in the TObjectArrays
  AliTRDdigitsManager(const AliTRDdigitsManager &m);
  virtual ~AliTRDdigitsManager();
  AliTRDdigitsManager &operator=(const AliTRDdigitsManager &m);

  virtual void                Copy(TObject &m) const;

  virtual void                CreateArrays();
  virtual void                ResetArrays();
  virtual void                ResetArrays(Int_t det);
  virtual Bool_t              BuildIndexes(Int_t det);

  virtual Bool_t              MakeBranch(TTree * const tree);
  virtual Bool_t              ReadDigits(TTree * consttree);
  virtual Bool_t              WriteDigits();

  virtual void                SetEvent(Int_t evt)             { fEvent           = evt;  };
  virtual void                SetSDigits(Int_t v = 1)         { fHasSDigits      = v;    };
  virtual void                SetUseDictionaries(Bool_t kval) { fUseDictionaries = kval; };

  virtual Bool_t              UsesDictionaries() const        { return fUseDictionaries; };
  virtual Bool_t              HasSDigits() const              { return fHasSDigits;      };
  static  Int_t               NDict()                         { return fgkNDict;         }; 

  virtual TObjArray          *GetDigits() const               { return fDigits;          };  
  virtual TObjArray          *GetDictionary(Int_t i) const    { return fDict[i];         }; 

  AliTRDdigit                *GetDigit(Int_t row, Int_t col, Int_t time, Int_t det) const;
  Int_t                       GetTrack(Int_t track, Int_t row, Int_t col, Int_t time, Int_t det) const;
  
  AliTRDarrayADC             *GetDigits(Int_t det)  const;
  AliTRDarraySignal          *GetSDigits(Int_t det) const;    
  AliTRDarrayDictionary      *GetDictionary(Int_t det, Int_t i) const;  
  AliTRDdigitsParam          *GetDigitsParam() const          { return fDigitsParam;     };
  AliTRDSignalIndex          *GetIndexes(Int_t det);
  TObjArray                  *GetIndexes() const              { return fSignalIndexes;   };

  void                        RemoveDigits(Int_t det);
  void                        RemoveDictionaries(Int_t det);
  void                        RemoveIndexes(Int_t det);
  void                        ClearIndexes(Int_t det);
  
  Int_t                       GetTrack(Int_t track, AliTRDdigit * const digit) const;
  Short_t                     GetDigitAmp(Int_t row, Int_t col, Int_t time, Int_t det) const;
  UChar_t                     GetPadStatus(Int_t row, Int_t col, Int_t time, Int_t det) const;

  Bool_t                      LoadArrayDigits();
  Bool_t                      LoadArrayDict();
  Bool_t                      LoadDigitsParam();
  Bool_t                      StoreArrayDigits();
  Bool_t                      StoreArrayDict();
  Bool_t                      StoreDigitsParam();

 protected:
  
  static const Int_t  fgkNDict;            //  Number of track dictionary arrays
  Int_t               fEvent;              //  Event number
  TTree              *fTree;               //! Tree for the digits arrays
  TObjArray          *fDigits;             //  Digits data array               
  TObjArray          *fDict[kNDict];       //  Track dictionary data array   
  Bool_t              fHasSDigits;         //  Switch for the summable digits
  TObjArray          *fSignalIndexes;      //  Provides access to the active pads and tbins
  Bool_t              fUseDictionaries;    //  Use dictionaries or not (case of real data)
  Int_t               fDets;               //  No of Detectors
  Bool_t              fRawRec;             //  Reconstruct from raw files? If its kTRUE then the TObjArrays have only one entry.
  AliTRDdigitsParam  *fDigitsParam;        //  Parameters of the digits

  ClassDef(AliTRDdigitsManager,8)          //  Manages the TRD digits

};
#endif
