// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveHOMERSourceMapByDet.h"
#include <AliHLTHOMERSourceDesc.h>

#include <TList.h>

//______________________________________________________________________
//
// AliEveHOMERSourceMap is an abstract container for HLT HOMER sources,
// see AliHLTHOMERSourceDesc.
//
// The concrete implementations AliEveHOMERSourceMapByDet and
// AliEveHOMERSourceMapByType allow retrieval of HOMER sources in proper
// order as required for their display in EVE object browser.
// 

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
