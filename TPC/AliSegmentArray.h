#ifndef ALISEGARRAY_H
#define ALISEGARRAY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class general Alice segment 
//  segment is for example one pad row in TPC //
////////////////////////////////////////////////

#include "TNamed.h"
#include "TError.h"
//#include "AliSegmentID.h"

class TTree;
class TBranch;
class AliArrayI;
class AliSegmentID;
class TObjArray;
 
class AliSegmentArray: public TNamed{
public:
  AliSegmentArray();
  AliSegmentArray(Text_t *classname, Int_t n);  //
  Bool_t  SetClass(Text_t *classname);  //set class of stored object
  ~AliSegmentArray();
  inline const AliSegmentID * At(Int_t i); //return pointer to segment with index i 
  inline const AliSegmentID * operator[](Int_t i); //return pointer to segment with index i

  Bool_t AddSegment(AliSegmentID *segment); // add segment to array
  AliSegmentID * AddSegment(Int_t index);   //create objet and set index
  Bool_t   MakeArray(Int_t n);       //make array of pointers to Segments
  void ClearSegment(Int_t index); //remove segment from active   
  virtual AliSegmentID * NewSegment(); //dynamicaly create new segment 
  //input output functions
  TTree * GetTree(){return fTree;}      //return pointer to connected tree
  
  virtual void MakeTree();              //Make tree with the name
  virtual Bool_t ConnectTree(const char * treeName); //connect tree from current directory 
  virtual AliSegmentID * LoadSegment(Int_t index);//load segment with index to the memory
  virtual AliSegmentID * LoadEntry(Int_t index); //load segment entry from position index in tree
  virtual void StoreSegment(Int_t index);//write segmen persistent  
  Bool_t  MakeDictionary(Int_t size);//create index table for tree
  TClass * GetClass() {return fClass;}
  
public:
  TObjArray  * fSegment;  //!pointer to array of pointers to segment
  AliArrayI    * fTreeIndex; //!pointers(index) table in tree
  Int_t      fNSegment;   
  TTree    * fTree;   //!tree with segment objects
  TBranch  * fBranch; //!total branch 
private: 
  TClass  *   fClass;    //!class type of included objects 
  ClassDef(AliSegmentArray,1) 
};



const AliSegmentID*  AliSegmentArray::operator[](Int_t i)
{
  //
  //return segment with given index
  //
  if ( (i<0) || (i>=fNSegment)) return 0; 
  return (AliSegmentID *)((*fSegment)[i]);
}
const AliSegmentID*  AliSegmentArray::At(Int_t i)
{
  //
  //return segment with given index
  //
  if ( (i<0) || (i>=fNSegment)) return 0; 
  return (AliSegmentID *)((*fSegment)[i]);
}

#endif //ALISEGARRAY_H
