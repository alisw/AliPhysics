// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveHOMERSourceMap.h"
#include "AliEveHOMERSourceMapByDet.h"
#include "AliEveHOMERSourceMapByType.h"

//______________________________________________________________________
//
// AliEveHOMERSourceMap is an abstract container for HLT HOMER sources,
// see AliHLTHOMERSourceDesc.
//
// The concrete implementations AliEveHOMERSourceMapByDet and
// AliEveHOMERSourceMapByType allow retrieval of HOMER sources in proper
// order as required for their display in EVE object browser.
// 

ClassImp(AliEveHOMERSourceMap)

AliEveHOMERSourceMap::AliEveHOMERSourceMap(ESourceGrouping_e grouping) :
  fGrouping(grouping)
{
  // Constructor.
}

AliEveHOMERSourceMap* AliEveHOMERSourceMap::Create(ESourceGrouping_e grouping)
{
  // Static constructor - instantiates appropriate sub-class.

  switch (grouping)
  {
    case kSG_ByDet:  return new AliEveHOMERSourceMapByDet(grouping);
    case kSG_ByType: return new AliEveHOMERSourceMapByType(grouping);
  }
  return 0;
}

Int_t AliEveHOMERSourceMap::iterator::level()
{
  // Returns the depth in iteration:
  // Det / Sub-Det / Sub-Sub-Det / Data-Type.

  const AliEveHOMERSource::SourceId& sid = id();

  Int_t lvl = 0;
  if ( ! sid.fDet.IsNull())   ++lvl;
  if ( ! sid.fSDet.IsNull())  ++lvl;
  if ( ! sid.fSSDet.IsNull()) ++lvl;
  if ( ! sid.fType.IsNull())  ++lvl;
  return lvl;
}

void AliEveHOMERSourceMap::PrintXXX()
{
  // Print entries in the map.

  for (iterator i = begin(); i != end(); ++i)
  {
    printf("%*s%s [state=%d, handle=0x%lx] {ssdet='%s'}\n", 4*i.level(), "",
	   i.description().Data(), i.state().fState,
	   (ULong_t) i.state().fHandle,
	   i.id().fSSDet.Data());
  }
}
