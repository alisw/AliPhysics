#ifndef ALITRDSEGMENTARRAYBASE_H
#define ALITRDSEGMENTARRAYBASE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDsegmentArrayBase.h,v */

////////////////////////////////////////////////
//  Manager class for a general Alice segment // 
////////////////////////////////////////////////

#include "TNamed.h"
#include "TError.h"
#include "TObjArray.h"

class TTree;
class TBranch;
class AliTRDarrayI;
class AliTRDsegmentID;
class TObjArray;
 
class AliTRDsegmentArrayBase: public TNamed {

 public:

  AliTRDsegmentArrayBase();
  AliTRDsegmentArrayBase(Text_t *classname, Int_t n); 
  AliTRDsegmentArrayBase(AliTRDsegmentArrayBase &a);
  virtual ~AliTRDsegmentArrayBase();
 
  const AliTRDsegmentID *At(Int_t i); 
  const AliTRDsegmentID *operator[](Int_t i); 

          Bool_t           AddSegment(AliTRDsegmentID *segment);
          AliTRDsegmentID *AddSegment(Int_t index);  
          void             ClearSegment(Int_t index); 
  virtual void             Copy(AliTRDsegmentArrayBase &a);
  virtual Bool_t           ConnectTree(const char *treeName);
          Bool_t           MakeArray(Int_t n);    
  virtual AliTRDsegmentID *NewSegment(); 
  virtual void             MakeTree();           
  virtual AliTRDsegmentID *LoadSegment(Int_t index);
  virtual AliTRDsegmentID *LoadEntry(Int_t index); 
  virtual void             StoreSegment(Int_t index);
          Bool_t           MakeDictionary(Int_t size);

          Bool_t           SetClass(Text_t *classname);
 
          TClass          *GetClass() { return fClass; };
          TTree           *GetTree()  { return fTree;  };   

  inline  AliTRDsegmentArrayBase &operator=(AliTRDsegmentArrayBase &a);

 protected:

  TObjArray    *fSegment;            //! Pointer to an array of pointers to a segment
  AliTRDarrayI *fTreeIndex;          //! Pointers(index) table
  Int_t         fNSegment;           //  Number of segments 
  TTree        *fTree;               //! Tree with segment objects
  TBranch      *fBranch;             //! Branchaddress
  TClass       *fClass;              //! Class type of included objects 

  ClassDef(AliTRDsegmentArrayBase,1) // TRD detextor segment array base class

};

//_____________________________________________________________________________
AliTRDsegmentArrayBase &AliTRDsegmentArrayBase
                        ::operator=(AliTRDsegmentArrayBase &a)
{
  //
  // Assignment operator
  //

  if (this != &a) a.Copy(*this);
  return *this;

}

#endif 
