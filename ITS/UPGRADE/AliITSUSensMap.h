#ifndef ALIITSUSENSMAP_H
#define ALIITSUSENSMAP_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//***********************************************************************
//
// It consist of a TClonesArray of objects
// and B-tree for fast sorted access
// This array can be accessed via 2 indexed
// it is used at digitization level by 
// all the ITS subdetectors
//
// The items should be added to the map like this:
// map->RegisterItem( new(map->GetFree()) ItemConstructor(...) );
//
// ***********************************************************************
#include <TClonesArray.h>
#include <TBtree.h>

#define _ROWWISE_SORT_

class AliITSUSensMap: public TObject 
{

 public:
  enum {kDisableBit=BIT(14)};
  //
  AliITSUSensMap();
  AliITSUSensMap(const char* className, UInt_t dimCol,UInt_t dimRow,UInt_t dimCycle=1);
  virtual ~AliITSUSensMap();
  AliITSUSensMap(const AliITSUSensMap &source);
  AliITSUSensMap& operator=(const AliITSUSensMap &source);
  void Clear(Option_t* option = "");
  void DeleteItem(UInt_t col,UInt_t row);
  void DeleteItem(TObject* obj);
  //
  void  SetDimensions(UInt_t dimCol,UInt_t dimRow,UInt_t dimCycle=1);
  void  GetMaxIndex(UInt_t &ni,UInt_t &nj,UInt_t &cycle) const {ni=fDimCol; nj=fDimRow; cycle=fDimCycle;}
  Int_t GetMaxIndex()                      const {return fDimCol*fDimRow;}
  Int_t GetEntries()                       const {return fBTree->GetEntries();}
  Int_t GetEntriesUnsorted()               const {return fItems->GetEntriesFast();}
  void  GetMapIndex(UInt_t index,UInt_t &col,UInt_t &row,UInt_t &cycle) const {return GetCell(index, col, row, cycle);}
  TObject* GetItem(UInt_t col,UInt_t row,UInt_t cycle)  {SetUniqueID(GetIndex(col,row,cycle)); return fBTree->FindObject(this);}
  TObject* GetItem(UInt_t index)                 {SetUniqueID(index);         return fBTree->FindObject(this);}
  TObject* GetItem(const TObject* obj)           {return fBTree->FindObject(obj);}
  TObject* At(Int_t i)                     const {return fBTree->At(i);}             //!!! Access in sorted order !!!
  TObject* AtUnsorted(Int_t i)             const {return fItems->At(i);}             //!!! Access in unsorted order !!!
  TObject* RegisterItem(TObject* obj)            {fBTree->Add(obj); return obj;}
  TObject* GetFree()                             {return (*fItems)[fItems->GetEntriesFast()];}
  //
  void   GetCell(UInt_t index,UInt_t &col,UInt_t &row, UInt_t &cycle)   const;
  UInt_t GetIndex(UInt_t col,UInt_t row,UInt_t cycle=0)  const;
  //
  TClonesArray* GetItems()                 const {return fItems;}
  TBtree*       GetItemsBTree()            const {return fBTree;}
  //
  Bool_t        IsSortable()                const {return kTRUE;}
  Bool_t        IsEqual(const TObject* obj) const {return GetUniqueID()==obj->GetUniqueID();}
  Int_t         Compare(const TObject* obj) const {return (GetUniqueID()<obj->GetUniqueID()) ? -1 : ((GetUniqueID()>obj->GetUniqueID()) ? 1 : 0 );}
  //
  static Bool_t IsDisabled(TObject* obj)         {return obj ? obj->TestBit(kDisableBit) : kFALSE;}
  static void   Disable(TObject* obj)            {if (obj) obj->SetBit(kDisableBit);}
  static void   Enable(TObject* obj)             {if (obj) obj->ResetBit(kDisableBit);}
  //
 protected:
  //
  UInt_t fDimCol;              // 1st dimension of the matrix
  UInt_t fDimRow;              // 2nd dimention of the matrix
  UInt_t fDimCycle;          // readout cycle
  TClonesArray*    fItems;   // pListItems array
  TBtree*          fBTree;   // tree for ordered access
  //
  ClassDef(AliITSUSensMap,1) // list of sensor signals (should be sortable objects)
};	

inline UInt_t AliITSUSensMap::GetIndex(UInt_t col,UInt_t row, UInt_t cycle) const  
{
  // linearized ID of digit
#ifdef _ROWWISE_SORT_
  return fDimCol*(cycle*fDimRow+row)+col; // sorted in row, then in column
#else
  return fDimRow*(cycle*fDimCol+col)+row; // sorted in column, then in row
#endif
}

#endif
