// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveHOMERSourceMap.h"
#include <AliHLTHOMERSourceDesc.h>

#include <TList.h>

//______________________________________________________________________
// AliEveHOMERSourceMap
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

/******************************************************************************/
// AliEveHOMERSourceMapByDet
/******************************************************************************/

AliEveHOMERSourceMapByDet::AliEveHOMERSourceMapByDet(ESourceGrouping_e grouping) :
  AliEveHOMERSourceMap(grouping),
  fMap()
{}

TString AliEveHOMERSourceMapByDet::iterator_imp::description() const
{
  // Return identifier string for current entry.

  const AliEveHOMERSource::SourceId& sid = id();

  if ( ! sid.fType.IsNull())  return sid.fType;
  if ( ! sid.fSSDet.IsNull()) return sid.fSSDet;
  if ( ! sid.fSDet.IsNull())  return sid.fSDet;
  if ( ! sid.fDet.IsNull())   return sid.fDet;
  return "<empty>";
}

void AliEveHOMERSourceMapByDet::insert(AliEveHOMERSource::SourceId& sid,
				       AliEveHOMERSource::SourceState& sst,
				       Bool_t def_state)
{
  // Insert source-state for given source-id.
  // Does nothing if the entry already exists.

  Map_i i = fMap.find(sid);
  if (i == fMap.end())
  {
    // Check wildcard, else
    sst.fState = def_state;
    fMap.insert(std::make_pair(sid, sst));
  }
}

void AliEveHOMERSourceMapByDet::FillMap(const TList* handles, Bool_t def_state)
{
  // Fill the map from the list of HOMER source handles.

  TIter next(handles);
  AliHLTHOMERSourceDesc* h;
  while ((h = (AliHLTHOMERSourceDesc*) next()))
  {
    AliEveHOMERSource::SourceId    srcid;
    AliEveHOMERSource::SourceState srcst;

    srcid.fDet = h->GetDetector();
    insert(srcid, srcst, def_state);

    srcid.fSDet = h->GetSubDetector();
    if ( ! srcid.fSDet.IsNull())  insert(srcid, srcst, def_state);

    srcid.fSSDet = h->GetSubSubDetector();
    if ( ! srcid.fSSDet.IsNull()) insert(srcid, srcst, def_state);

    srcid.fType   = h->GetDataType();
    srcst.fHandle = h;
    insert(srcid, srcst, def_state);
  }
}


/******************************************************************************/
// AliEveHOMERSourceMapByType
/******************************************************************************/

AliEveHOMERSourceMapByType::AliEveHOMERSourceMapByType(ESourceGrouping_e grouping) :
  AliEveHOMERSourceMap(grouping),
  fMap()
{
  // Constructor.
}

TString AliEveHOMERSourceMapByType::iterator_imp::description() const
{
  // Return identifier string for current entry.

  const AliEveHOMERSource::SourceId& sid = id();

  if ( ! sid.fSSDet.IsNull()) return sid.fSSDet;
  if ( ! sid.fSDet.IsNull())  return sid.fSDet;
  if ( ! sid.fDet.IsNull())   return sid.fDet;
  if ( ! sid.fType.IsNull())  return sid.fType;
  return "<empty>";
}

void AliEveHOMERSourceMapByType::insert(AliEveHOMERSource::SourceId& sid,
					AliEveHOMERSource::SourceState& sst,
					Bool_t def_state)
{
  // Insert source-state for given source-id.
  // Does nothing if the entry already exists.

  Map_i i = fMap.find(sid);
  if (i == fMap.end())
  {
    // Check wildcard, else
    sst.fState = def_state;
    fMap.insert(std::make_pair(sid, sst));
  }
}

void AliEveHOMERSourceMapByType::FillMap(const TList* handles, Bool_t def_state)
{
  // Fill the map from the list of HOMER source handles.

  TIter next(handles);
  AliHLTHOMERSourceDesc* h;
  while ((h = (AliHLTHOMERSourceDesc*) next()))
  {
    AliEveHOMERSource::SourceId    srcid;
    AliEveHOMERSource::SourceState srcst;

    srcid.fType   = h->GetDataType();
    insert(srcid, srcst, def_state);

    srcid.fDet = h->GetDetector();
    if (h->GetSubDetector() == 0)
    {
      srcst.fHandle = h;
      insert(srcid, srcst, def_state);
      continue;
    }
    insert(srcid, srcst, def_state);

    srcid.fSDet = h->GetSubDetector();
    if (h->GetSubSubDetector() == 0)
    {
      srcst.fHandle = h;
      insert(srcid, srcst, def_state);
      continue;
    }
    insert(srcid, srcst, def_state);

    srcid.fSSDet  = h->GetSubSubDetector();
    srcst.fHandle = h;
    insert(srcid, srcst, def_state);
  }
}
