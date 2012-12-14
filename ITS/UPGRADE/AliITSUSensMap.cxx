/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id: AliITSUSensMap.cxx 39464 2010-03-09 14:59:19Z masera $ */

//***********************************************************************
//
// It consist of a TClonesArray of 
// AliITSUSensMapItem objects
// This array can be accessed via 2 indexed
// it is used at digitization level by 
// all the 3 ITS subdetectors
//
//
// The items should be added to the map like this:
// map->RegisterItem( new(map->GetFree()) ItemConstructor(...) );
//
// The items must be sortable with the same sorting algorithm like 
// for AliITSUSensMap::IsSortable,IsEqual,Compare
//
// ***********************************************************************

#include "AliITSUSensMap.h"
#include "AliLog.h"
//______________________________________________________________________

ClassImp(AliITSUSensMap)
//______________________________________________________________________
AliITSUSensMap::AliITSUSensMap() 
:  fDim0(0)
  ,fDim1(0)
  ,fItems(0)
  ,fBTree(0)
{
  // Default constructor
}

//______________________________________________________________________
AliITSUSensMap::AliITSUSensMap(const char* className, UInt_t dim0,UInt_t dim1)
  :fDim0(dim0)
  ,fDim1(dim1)
  ,fItems(new TClonesArray(className,100))
  ,fBTree(new TBtree())
{
  // Standard constructor
}

//______________________________________________________________________
AliITSUSensMap::~AliITSUSensMap() 
{
  // Default destructor
  delete fItems;
  delete fBTree;
}


//______________________________________________________________________
AliITSUSensMap::AliITSUSensMap(const AliITSUSensMap &source)
  :TObject(source)
  ,fDim0(source.fDim0)
  ,fDim1(source.fDim1)
  ,fItems( source.fItems ? new TClonesArray(*source.fItems) : 0)
  ,fBTree( 0 )
{
  if (source.fBTree) {
    fBTree = new TBtree();
    if (fItems) {
      for (int i=fItems->GetEntriesFast();i--;) {
	TObject* obj = fItems->At(i);
	if (obj && ! IsDisabled(obj)) continue;
	RegisterItem(obj);
      }
    }
  }
}

//______________________________________________________________________
AliITSUSensMap& AliITSUSensMap::operator=(const AliITSUSensMap &source)
{
  // = operator
  if (this!=&source) {
    this->~AliITSUSensMap();
    new(this) AliITSUSensMap(source);
  }
  return *this;
}

//______________________________________________________________________
void AliITSUSensMap::Clear(Option_t*) 
{
  // clean everything
  if (fItems) fItems->Clear();
  if (fBTree) fBTree->Clear();
}

//______________________________________________________________________
void AliITSUSensMap::DeleteItem(UInt_t i,UInt_t j)
{
  // Delete a particular AliITSUSensMapItems.
  SetUniqueID( GetIndex(i,j) );
  TObject* fnd = fBTree->FindObject(this);
  if (!fnd) return;
  Disable(fnd);
  fBTree->Remove(fnd);
}

//______________________________________________________________________
void AliITSUSensMap::DeleteItem(TObject* obj)
{
  // Delete a particular AliITSUSensMapItems.
  TObject* fnd = fBTree->FindObject(obj);
  if (!fnd) return;
  Disable(fnd);
  fBTree->Remove(fnd);
}

//______________________________________________________________________
void AliITSUSensMap::GetCell(UInt_t index,UInt_t &i,UInt_t &j) const 
{
  // returns the i,j index numbers from the linearized index computed
  // with GetIndex
  if(index>=fDim0*fDim1){
    Warning("GetCell","Index out of range 0<=index=%d<%d",index,fDim0*fDim1);
    i=-1;j=-1;
    return;
  } // end if
#ifdef _ROWWISE_SORT_
  i = index%fDim0;   // sorted in row, then in column
  j = index/fDim0;
#else
  i = index/fDim1;   // sorted in column, then in row
  j = index%fDim1;
#endif  
  //
  return;
}
