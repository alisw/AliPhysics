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

class AliOADBContainer : public TNamed {

 public :
  AliOADBContainer();
  AliOADBContainer(char* name);
  virtual ~AliOADBContainer();
  AliOADBContainer(const AliOADBContainer& cont); 
  AliOADBContainer& operator=(const AliOADBContainer& cont);
// Object adding and removal
  void   AppendObject(TObject* obj, Int_t lower, Int_t upper);
  void   UpdateObject(Int_t index, TObject* obj, Int_t lower, Int_t upper);
  void   RemoveObject(Int_t index);
  void   AddDefaultObject(TObject* obj);
  void   CleanDefaultList();
  TList* GetDefaultList() const {return fDefaultList;}
// I/O  
  void  WriteToFile(char* fname)  const;
  Int_t InitFromFile(char* fname, char* key);
// Getters
  Int_t GetNumberOfEntries()    const {return fEntries;}
  Int_t LowerLimit(Int_t idx)   const {return fLowerLimits[idx];}
  Int_t UpperLimit(Int_t idx)   const {return fUpperLimits[idx];}
  TObject* GetObject(Int_t run, char* def = "") const;
  TObject* GetObjectByIndex(Int_t run) const;
  TObject* GetDefaultObject(char* key) 
  {return(fDefaultList->FindObject(key));}
// Debugging  
  void List();
// Browsable
  virtual Bool_t	IsFolder() const { return kTRUE; }
  void Browse(TBrowser *b);

 private:
  Int_t HasOverlap(Int_t lower, Int_t upper) const;
  Int_t GetIndexForRun(Int_t run) const;
 private :
  TObjArray*               fArray;         // Array with objects corresponding to run ranges
  TList*                   fDefaultList;   // List with default arrays
  TArrayI                  fLowerLimits;   // lower limit of run range
  TArrayI                  fUpperLimits;   // upper limit of run range
  Int_t                    fEntries;       // Number of entries
//  TString                  fRelPath;       // Relative path to object
  
  ClassDef(AliOADBContainer, 1);
};

#endif
