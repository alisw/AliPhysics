#ifndef ALITRDSEGMENTARRAY_H
#define ALITRDSEGMENTARRAY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliTRDsegmentArrayBase.h"

////////////////////////////////////////////////////////
//  Array for TRD detector segments containing digits //
////////////////////////////////////////////////////////

class AliTRDdataArray;

//_____________________________________________________________________________
class AliTRDsegmentArray : public AliTRDsegmentArrayBase {

 public:

  AliTRDsegmentArray();
  AliTRDsegmentArray(Text_t *classname, Int_t n);
  AliTRDsegmentArray(AliTRDsegmentArray &a);
  virtual ~AliTRDsegmentArray();

  virtual void             Copy(TObject &a);
  virtual void             Delete();
  virtual void             Delete(const char *) { Delete(); };

  virtual Bool_t           LoadArray(const Char_t *branchname);
  virtual Bool_t           StoreArray(const Char_t *branchname);

  virtual AliTRDdataArray *GetDataArray(Int_t det) const;
  virtual AliTRDdataArray *GetDataArray(Int_t sec, Int_t cha, Int_t pla) const;

 protected:

  ClassDef(AliTRDsegmentArray,1)        // TRD detector segment array 

};

#endif
