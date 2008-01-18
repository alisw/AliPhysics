// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveHOMERSourceList.h"

//______________________________________________________________________________
// AliEveHOMERSourceList
//

ClassImp(AliEveHOMERSourceList)

AliEveHOMERSourceList::AliEveHOMERSourceList(const Text_t* n, const Text_t* t) :
  TEveElementList(n, t)
{

}

void AliEveHOMERSourceList::SelectAll()
{
  EnableListElements(kTRUE, kTRUE);
}

void AliEveHOMERSourceList::DeselectAll()
{
  DisableListElements (kFALSE, kFALSE);
}
