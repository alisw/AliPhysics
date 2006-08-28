#ifndef ALITRDSEGMENTARRAYBASE_H
#define ALITRDSEGMENTARRAYBASE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Manager class for a general Alice segment                             // 
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>

class TTree;
class TBranch;
class TObjArray;

class AliTRDarrayI;
class AliTRDsegmentID;
 
class AliTRDsegmentArrayBase: public TNamed {

 public:

  AliTRDsegmentArrayBase();
  AliTRDsegmentArrayBase(const char *classname, Int_t n); 
  AliTRDsegmentArrayBase(const AliTRDsegmentArrayBase &a);
  virtual ~AliTRDsegmentArrayBase();
  AliTRDsegmentArrayBase &operator=(const AliTRDsegmentArrayBase &a);
 
  const AliTRDsegmentID *At(Int_t i) const; 
  const AliTRDsegmentID *operator[](Int_t i) const; 

          Bool_t           AddSegment(AliTRDsegmentID *segment);
          AliTRDsegmentID *AddSegment(Int_t index);  
  virtual AliTRDsegmentID *NewSegment(); 
  virtual AliTRDsegmentID *LoadSegment(Int_t index);
  virtual AliTRDsegmentID *LoadEntry(Int_t index); 
  virtual void             StoreSegment(Int_t index);
          void             ClearSegment(Int_t index);
 
  virtual void             Copy(TObject &a) const;
  virtual Bool_t           ConnectTree(const char *treeName);

          Bool_t           MakeArray(Int_t n);    
  virtual void             MakeTree(char *file = 0);           
          Bool_t           MakeDictionary(Int_t size);

          Bool_t           SetClass(const char *classname);
 
          TClass          *GetClass() const { return fClass; };
          TTree           *GetTree() const  { return fTree;  };   

 protected:

  TObjArray    *fSegment;            //! Pointer to an array of pointers to a segment
  AliTRDarrayI *fTreeIndex;          //! Pointers(index) table
  Int_t         fNSegment;           //  Number of segments 
  TTree        *fTree;               //! Tree with segment objects
  TBranch      *fBranch;             //! Branchaddress
  TClass       *fClass;              //! Class type of included objects 

  ClassDef(AliTRDsegmentArrayBase,1) // TRD detextor segment array base class

};

#endif 
