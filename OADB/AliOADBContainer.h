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
#include <TArrayI.h>


class TObjArray;
class TArrayI;

class AliOADBContainer : public TNamed {

 public :
  AliOADBContainer();
  AliOADBContainer(char* name);
  virtual ~AliOADBContainer();
  AliOADBContainer(const AliOADBContainer& cont); 
  AliOADBContainer& operator=(const AliOADBContainer& cont);
  void AppendObject(TObject* obj, Int_t lower, Int_t upper);
  void UpdateObject(Int_t index, TObject* obj, Int_t lower, Int_t upper);
  void RemoveObject(Int_t index);
  Int_t GetIndexForRun(Int_t run);
  void WriteToFile(char* fname);
  void List();
 private :
  TObjArray*               fArray;         // Array with objects
  TArrayI                  fLowerLimits;   // lower limit of run range
  TArrayI                  fUpperLimits;   // upper limit of run range
  Int_t                    fEntries;       // Number of entries
  ClassDef(AliOADBContainer, 1);
};

#endif
