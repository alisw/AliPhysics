#ifndef AliOADBContainer_H
#define AliOADBContainer_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     Offline Analysis Database Container and Service Class 
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include <TNamed.h>
#include <TList.h>
#include <TArrayI.h>
#include <TObjArray.h>

class TObjArray;
class TArrayI;
class TFile;
class AliOADBContainer : public TNamed {

 public :
  AliOADBContainer();
  AliOADBContainer(const char* name);
  virtual ~AliOADBContainer();
  AliOADBContainer(const AliOADBContainer& cont); 
  AliOADBContainer& operator=(const AliOADBContainer& cont);
// Object adding and removal
  void   AppendObject(TObject* obj, Int_t lower, Int_t upper, TString passName="");
  void   UpdateObject(Int_t index, TObject* obj, Int_t lower, Int_t upper, TString passName="");
  void   RemoveObject(Int_t index);
  void   AddDefaultObject(TObject* obj);
  void   CleanDefaultList();
  TList* GetDefaultList() const {return fDefaultList;}
// I/O  
  void  WriteToFile(const char* fname)  const;
  Int_t InitFromFile(const char* fname, const char* key);
// Getters
  Int_t GetNumberOfEntries()    const {return fEntries;}
  Int_t LowerLimit(Int_t idx)   const {return fLowerLimits[idx];}
  Int_t UpperLimit(Int_t idx)   const {return fUpperLimits[idx];}
  TObjArray* GetObjArray() {return fArray;}
  void SetToZeroObjArray() {fArray=0;}
  TObject* GetObject(Int_t run, const char* def = "", TString passName="") const;
  TObject* GetObjectFromFile(TFile* file, Int_t run, const char* def = "", TString passName="") const;
  TObject* GetObjectByIndex(Int_t run) const;
  TObject* GetPassNameByIndex(Int_t idx) const;
  TObject* GetDefaultObject(const char* key) 
  {return(fDefaultList->FindObject(key));}
// Debugging  
  void List();
// Browsable
  virtual Bool_t	IsFolder() const { return kTRUE; }
  void Browse(TBrowser *b);
  Int_t GetIndexForRun(Int_t run, TString passName="") const;
//
  static const char*   GetOADBPath();
 private:
  Int_t HasOverlap(Int_t lower, Int_t upper, TString passName) const;
 private :
  TObjArray*               fArray;         //Array with objects corresponding to run ranges
  TList*                   fDefaultList;   // List with default arrays
  TObjArray*               fPassNames;     // Pass names
  TArrayI                  fLowerLimits;   // lower limit of run range
  TArrayI                  fUpperLimits;   // upper limit of run range
  Int_t                    fEntries;       // Number of entries
  ClassDef(AliOADBContainer, 2);
};

#endif
