#ifndef ALISEGARRAYBASE_H
#define ALISEGARRAYBASE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDsegmentArrayBase.h,v */

////////////////////////////////////////////////
//  Manager class generaol Alice segment 
//  segment is for example one pad row in TPC //
////////////////////////////////////////////////

#include "TNamed.h"
#include "TError.h"
#include "TObjArray.h"

//#include "AliTRDsegmentID.h"

class TTree;
class TBranch;
class AliTRDarrayI;
class AliTRDsegmentID;
class TObjArray;
 
class AliTRDsegmentArrayBase: public TNamed{
public:
  AliTRDsegmentArrayBase();
  AliTRDsegmentArrayBase(Text_t *classname, Int_t n);  //
  Bool_t  SetClass(Text_t *classname);  //set class of stored object
  ~AliTRDsegmentArrayBase();
  const AliTRDsegmentID * At(Int_t i); //return pointer to segment with index i 
  const AliTRDsegmentID * operator[](Int_t i); //return pointer to segment with index i

  Bool_t AddSegment(AliTRDsegmentID *segment); // add segment to array
  AliTRDsegmentID * AddSegment(Int_t index);   //create objet and set index
  Bool_t   MakeArray(Int_t n);       //make array of pointers to Segments
  void ClearSegment(Int_t index); //remove segment from active   
  virtual AliTRDsegmentID * NewSegment(); //dynamicaly create new segment 
  //input output functions
  TTree * GetTree(){return fTree;}      //return pointer to connected tree
  
  virtual void MakeTree();              //Make tree with the name
  virtual Bool_t ConnectTree(const char * treeName); //connect tree from current directory 
  virtual AliTRDsegmentID * LoadSegment(Int_t index);//load segment with index to the memory
  virtual AliTRDsegmentID * LoadEntry(Int_t index); //load segment entry from position index in tree
  virtual void StoreSegment(Int_t index);//write segmen persistent  
  Bool_t  MakeDictionary(Int_t size);//create index table for tree
  TClass * GetClass() {return fClass;}
public:
  TObjArray  * fSegment;  //!pointer to array of pointers to segment
  AliTRDarrayI    * fTreeIndex; //!pointers(index) table in tree
  Int_t      fNSegment;   
  TTree    * fTree;   //!tree with segment objects
  TBranch  * fBranch; //!total branch 
private: 
  TClass  *   fClass;    //!class type of included objects 
  ClassDef(AliTRDsegmentArrayBase,1) 
};

#endif //ALISEGARRAY_H
