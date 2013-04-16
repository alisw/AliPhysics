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
:  fDimCol(0)
  ,fDimRow(0)
  ,fDimCycle(0)
  ,fItems(0)
  ,fBTree(0)
{
  // Default constructor
}

//______________________________________________________________________
AliITSUSensMap::AliITSUSensMap(const char* className, UInt_t dimCol,UInt_t dimRow,UInt_t dimCycle)
  :fDimCol(dimCol)
  ,fDimRow(dimRow)
  ,fDimCycle(dimCycle)
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
  ,fDimCol(source.fDimCol)
  ,fDimRow(source.fDimRow)
  ,fDimCycle(source.fDimCycle)
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
void AliITSUSensMap::DeleteItem(UInt_t col,UInt_t row,Int_t cycle)
{
  // Delete a particular AliITSUSensMapItems.
  SetUniqueID( GetIndex(col,row,cycle) );
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
void  AliITSUSensMap::SetDimensions(UInt_t dimCol,UInt_t dimRow,UInt_t dimCycle) 
{
  // set dimensions for current sensor
  const UInt_t kMaxPackDim = 0xffffffff;
  fDimCol = dimCol; 
  fDimRow = dimRow; 
  fDimCycle=dimCycle;
  if ((fDimCol*fDimRow*fDimCycle*2)>kMaxPackDim) AliFatal(Form("Dimension %dx%dx%d*2 cannot be packed to UInt_t",fDimCol,fDimRow,fDimCycle));
}

