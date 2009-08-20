//-*- Mode: C++ -*-

// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveHOMERSourceList.h"
#include "AliEveHOMERSourceMap.h"
#include "AliEveHOMERManager.h"

//______________________________________________________________________________
// AliEveHOMERSourceList
//

ClassImp(AliEveHOMERSourceList)

AliEveHOMERSourceList::AliEveHOMERSourceList(const Text_t* n, const Text_t* t) :
  TEveElementList(n, t),
  fManager (0),
  fSrcMap  (0)
{

}

AliEveHOMERSourceList::~AliEveHOMERSourceList()
{
  // !!!!! delete maps
}

/******************************************************************************/

void AliEveHOMERSourceList::CreateByDet()
{
  delete fSrcMap;
  fSrcMap = AliEveHOMERSourceMap::Create(AliEveHOMERSourceMap::kSG_ByDet);
  RebuildSourceReps();
}

void AliEveHOMERSourceList::CreateByType()
{
  delete fSrcMap;
  fSrcMap = AliEveHOMERSourceMap::Create(AliEveHOMERSourceMap::kSG_ByType);
  RebuildSourceReps();
}

void AliEveHOMERSourceList::RebuildSourceReps()
{
  DestroyElements();
  TList* srcList = fManager->GetSourceList();
  fSrcMap->FillMap(srcList, 1);

  List_t parentStack;
  parentStack.push_back(this);
  Int_t parentLvl = 1;
  for (AliEveHOMERSourceMap::iterator i=fSrcMap->begin(); i!=fSrcMap->end(); ++i)
  {
    while (parentLvl > i.level()) { parentStack.pop_back(); --parentLvl; }

    AliEveHOMERSource* src = new AliEveHOMERSource(i.description());
    src->SetSource(&i.id(), &i.state());

    parentStack.back()->AddElement(src);

    parentStack.push_back(src); ++parentLvl;
    
    printf("%*s%s [state=%d, handle=0x%lx] {ssdet='%s'}\n", 4*i.level(), "",
	   i.description().Data(), i.state().fState,
	   (ULong_t) i.state().fHandle,
	   i.id().fSSDet.Data());

    
  }
}


Bool_t AliEveHOMERSourceList::GetSelectedSources() {
  // Set selected source in HOMER sources list, of HOMERManager

  if ( ! fManager ) {
    printf ( "Error : no ptr to HomerManager!");
    return kFALSE;
  }
    

  Bool_t bResult = kFALSE;

  for ( AliEveHOMERSourceMap::iterator iter=fSrcMap->begin(); iter!=fSrcMap->end(); ++iter ) {

    if ( ! iter.state().fHandle ) 
      continue;
    
    fManager->SetSourceState( (AliHLTHOMERSourceDesc*) iter.state().fHandle,iter.state().fState );
    bResult = kTRUE;

#if 0 // EVE_DEBUG   
    printf("%*s%s [state=%d, handle=0x%lx] {ssdet='%s'}\n", 4*iter.level(), "",
	   iter.description().Data(), iter.state().fState,
	   (ULong_t) iter.state().fHandle,
	   iter.id().fSSDet.Data());
    
#endif
    


  }


  return bResult;
}


/******************************************************************************/
/*
void AliEveHOMERSourceList::SelectAll()
{
  EnableListElements(kTRUE, kTRUE);
}

void AliEveHOMERSourceList::DeselectAll()
{
  DisableListElements (kFALSE, kFALSE);
}
*/
