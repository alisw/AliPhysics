#ifndef ALIARRAYBRANCH_H
#define ALIARRAYBRANCH_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliArrayBranch                                                        //
//                                                                      //
// A Branch for the case of an array of clone objects.                  //
//////////////////////////////////////////////////////////////////////////


#include "TBranch.h"
#include "TTree.h"
#include "TBranchObject.h"
class AliObjectArray;
 
class AliArraySubBranch : public TBranch {
public: 
  AliArraySubBranch(){;}
  AliArraySubBranch(const char* name, void* address, const char* leaflist, Int_t basketsize = 32000, 
    Int_t compress = -1):TBranch(name, address, leaflist, basketsize, compress){;}
    virtual Int_t  GetEntryExport(Int_t entry, Int_t getall, AliObjectArray* list, Int_t n);
  virtual void ReadBasketExport(TBuffer &b, TLeaf *leaf, AliObjectArray *list, Int_t n); 
  ClassDef(AliArraySubBranch,1)  //Branch in case of an array of clone objects
};

class AliArrayBranch : public TBranch {

private: 
  void Import(TLeaf * leaf, Int_t n);  //integer fill leef buffer 
protected:
    AliObjectArray     *fList;           //Pointer to the clonesarray
    Int_t            fRead;            //flag = 1 if clonesarray has been read
    Int_t            fN;               //Number of elements in ClonesArray
    Int_t            fNdataMax;        //Maximum value of fN
    TString          fClassName;       //name of the class of the objets in the ClonesArray
    TBranch          *fBranchCount;    //Branch with clones count

public:
    AliArrayBranch();
    AliArrayBranch(const Text_t *name, void *clonesaddress, TTree * tree, Int_t basketsize=32000,Int_t compress=-1);
    virtual ~AliArrayBranch();

    virtual void    Browse(TBrowser *b);
    virtual Int_t   Fill();
    virtual Int_t   GetEntry(Int_t entry=0, Int_t getall = 0);
    virtual Int_t   GetN() {return fN;}
    AliObjectArray    *GetList() {return fList;}
    Bool_t          IsFolder() {return kTRUE;}
    virtual void    Print(Option_t *option="");
    virtual void    Reset(Option_t *option="");
    virtual void    SetAddress(void *add);
    virtual void    SetBasketSize(Int_t buffsize);
    virtual Bool_t          IsFolder() const {return kTRUE;}
    ClassDef(AliArrayBranch,1)  //Branch in case of an array of clone objects
};


class AliObjectBranch: public TBranchObject{
public:
  AliObjectBranch():TBranchObject(){;}
  AliObjectBranch(const Text_t *name, const Text_t *classname, void *addobj, TTree * tree, 
		  Int_t basketsize=32000, Int_t splitlevel = 0, Int_t compress=-1);
  void SetAddress(void *add);
  ClassDef(AliObjectBranch,1) 
};

class AliTree : public TTree {
public:
  AliTree():TTree(){;}
  AliTree(const char *name,const char *title, Int_t maxvirtualsize=0);
  TBranch* AliBranch(const char *name, void *clonesaddress, Int_t bufsize =32000, 
		     Int_t splitlevel=1,Int_t compres=1);
  TBranch* AliBranch(const char *name, const char *classname, void *addobj, 
  		     Int_t bufsize=32000, Int_t splitlevel=1);
  ClassDef(AliTree,1)  
};

#endif
