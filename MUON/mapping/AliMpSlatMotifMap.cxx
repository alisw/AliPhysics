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

/// 
/// \class AliMpSlatMotifMap
//
/// Basically this class provide a garbage collection of AliMpMotif and
/// AliMpMotifType objects.
///
///
/// \author Laurent Aphecetche


// $Id$

#include "AliMpSlatMotifMap.h"

#include "AliMpVMotif.h"
#include "AliMpMotifType.h"
#include "AliLog.h"
#include "TList.h"
#include "TObjString.h"
#include "TString.h"
#include "Riostream.h"

ClassImp(AliMpSlatMotifMap)

//_____________________________________________________________________________
AliMpSlatMotifMap::AliMpSlatMotifMap()
: TObject(),
fMotifs(),
fMotifTypes()
{
  /// ctor
  fMotifs.SetOwner(kTRUE);
  fMotifTypes.SetOwner(kTRUE);
}

//_____________________________________________________________________________
AliMpSlatMotifMap::~AliMpSlatMotifMap()
{
  /// dtor
  Reset();
}

//_____________________________________________________________________________
void
AliMpSlatMotifMap::Reset()
{
  /// Clear
  fMotifs.DeleteAll();
  fMotifTypes.DeleteAll();
}

//_____________________________________________________________________________
Bool_t 
AliMpSlatMotifMap::AddMotif(AliMpVMotif* motif, Bool_t warn)
{
  /// Add a motif to the map
  AliDebug(1,Form("Adding motif %s",motif->GetID().Data()));
  
  AliMpVMotif* found = FindMotif(motif->GetID());
  if (found) {    
    if (warn && found == motif) 
    {
      AliWarning(Form("The motif %s is already in map",motif->GetID().Data()));
    }
    if (warn && found != motif) 
    {
      AliError(Form("Another motif with the same ID=%s is already in map",
                    motif->GetID().Data()));
    }	     
    return false;
  }  
  
  fMotifs.Add(new TObjString(motif->GetID()),motif);
  return true;
}

//_____________________________________________________________________________
Bool_t 
AliMpSlatMotifMap::AddMotifType(AliMpMotifType* motifType, Bool_t warn)
{
  /// Add a motif to the map

  AliDebug(1,Form("Adding motifType %s",motifType->GetID().Data()));

  AliMpMotifType* found = FindMotifType(motifType->GetID());
  if (found) {    
    if (warn && found == motifType) 
    {
      AliWarning(Form("The motifType %s is already in map",
                      motifType->GetID().Data()));
    }
    if (warn && found != motifType)   
    {
      AliError(Form("Another motifType with the same ID=%s is already in map",
                    motifType->GetID().Data()));
    }	     
    return false;
  }  
  
  fMotifTypes.Add(new TObjString(motifType->GetID()),motifType);
  return true;
  
}

//_____________________________________________________________________________
AliMpVMotif* 
AliMpSlatMotifMap::FindMotif(const TString& id) const
{
  /// Search a given motif in the map and returns it if it's there.
  
  AliDebug(1,Form("Looking for motif %s",id.Data()));
  
  TObject* object = fMotifs.GetValue(id.Data());
  
  if (object)
  {
    AliMpVMotif* motif = static_cast<AliMpVMotif*>(object);
    AliDebug(1,Form("Found : %p id=%s",motif,motif->GetID().Data()));
    return motif;
  }
  AliDebug(1,"Not found");
  return 0x0;
}

//_____________________________________________________________________________
AliMpMotifType* 
AliMpSlatMotifMap::FindMotifType(const TString& id) const
{
  /// Search a given motifType in the map and returns it if it's there.
  AliDebug(1,Form("Looking for motifType %s",id.Data()));
  
  TObject* object = fMotifTypes.GetValue(id.Data());
  
  if (object)
  {
    AliMpMotifType* motifType = static_cast<AliMpMotifType*>(object);
    AliDebug(1,Form("Found : %p id=%s",motifType,motifType->GetID().Data()));
    return motifType;
  }
  AliDebug(1,"Not found");
  return 0x0;
  
}

//_____________________________________________________________________________
void
AliMpSlatMotifMap::Print(Option_t*) const
{
  /// printout
  cout << "Motifs=" << endl;
  TObject* key;
  TIter next(&fMotifs);
  while ( ( key = next() ) ) 
  {
    AliMpVMotif* motif = dynamic_cast<AliMpVMotif*>(fMotifs.GetValue(key));
    if (motif) cout << motif->GetID() << endl;
  }

  cout << "MotifTypes=" << endl;
  TIter tnext(&fMotifTypes);
  while ( ( key = tnext() ) ) 
  {
    AliMpMotifType* motifType = dynamic_cast<AliMpMotifType*>(fMotifTypes.GetValue(key));
    if (motifType) cout << motifType->GetID() << endl;
  }
  
}
