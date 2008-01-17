#ifndef ALITRDSEGMENTARRAY_H
#define ALITRDSEGMENTARRAY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliTRDsegmentArrayBase.h"

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Array for TRD detector segments containing digits                     //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>

class TTree;
class TBranch;
class TObjArray;

class AliTRDarrayI;
class AliTRDdataArray;
class AliTRDsegmentID;

class AliTRDsegmentArray : public TNamed {

 public:

  AliTRDsegmentArray();
  AliTRDsegmentArray(const char *classname, Int_t n);
  AliTRDsegmentArray(AliTRDsegmentArray &a);
  virtual ~AliTRDsegmentArray();
  AliTRDsegmentArray &operator=(const AliTRDsegmentArray &a);

  const   AliTRDsegmentID *At(Int_t i) const; 
  const   AliTRDsegmentID *operator[](Int_t i) const; 

          Bool_t           AddSegment(AliTRDsegmentID *segment);
          AliTRDsegmentID *AddSegment(Int_t index);  
  virtual AliTRDsegmentID *NewSegment(); 
  virtual AliTRDsegmentID *LoadSegment(Int_t index);
  virtual AliTRDsegmentID *LoadEntry(Int_t index); 
  virtual void             StoreSegment(Int_t index);
          void             ClearSegment(Int_t index);
 
  virtual void             Copy(TObject &a) const;
  virtual Bool_t           ConnectTree(const char *treeName);
  virtual void             Delete();
  virtual void             Delete(const char *) { Delete(); };

          Bool_t           MakeArray(Int_t n);    
  virtual void             MakeTree(char *file = 0);           
          Bool_t           MakeDictionary(Int_t size);

          Bool_t           SetClass(const char *classname);
 
          TClass          *GetClass() const { return fClass; };
          TTree           *GetTree() const  { return fTree;  };   

  virtual Bool_t           LoadArray(const Char_t *branchname, TTree *tree = 0);
  virtual Bool_t           StoreArray(const Char_t *branchname, TTree *tree = 0);

  virtual AliTRDdataArray *GetDataArray(Int_t det) const;
  virtual AliTRDdataArray *GetDataArray(Int_t sec, Int_t cha, Int_t pla) const;

 protected:

  TObjArray    *fSegment;            //! Pointer to an array of pointers to a segment
  AliTRDarrayI *fTreeIndex;          //! Pointers(index) table
  Int_t         fNSegment;           //  Number of segments 
  TTree        *fTree;               //! Tree with segment objects
  TBranch      *fBranch;             //! Branchaddress
  TClass       *fClass;              //! Class type of included objects 

  ClassDef(AliTRDsegmentArray,2)     //  TRD detector segment array 

};

#endif
