// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE Project            * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Sedat Altinpinar <Sedat.Altinpinar@cern.ch>           *
//*                  Hege Erdal       <hege.erdal@gmail.com>               *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliDxHFEParticleSelection.cxx
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-03-19
/// @brief  Base class for particle selection
///

#include "AliDxHFEParticleSelection.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "TObjArray.h"
#include <cerrno>

/// ROOT macro for the implementation of ROOT specific class methods
ClassImp(AliDxHFEParticleSelection)

AliDxHFEParticleSelection::AliDxHFEParticleSelection(const char* opt)
  : TObject()
  , fOption(opt)
  , fSelectedTracks(NULL)
{
  // constructor
  // 
  // 
  // 
  // 
}

AliDxHFEParticleSelection::~AliDxHFEParticleSelection()
{
  // destructor
  if (fSelectedTracks) delete fSelectedTracks;
  fSelectedTracks=NULL;
}

TObjArray* AliDxHFEParticleSelection::Select(const AliVEvent* pEvent)
{
  /// create selection, array contains only pointers but does not own the objects
  /// object array needs to be deleted by caller
  if (!pEvent) return NULL;
  TObjArray* selectedTracks=new TObjArray;
  if (!selectedTracks) return NULL;
  int nofTracks=pEvent->GetNumberOfTracks();
  for (int itrack=0; itrack<nofTracks; itrack++) {
    AliVParticle* track=pEvent->GetTrack(itrack);
    if (!IsSelected(track)) continue;
    selectedTracks->Add(track);
  }
  return selectedTracks;
}

TObjArray* AliDxHFEParticleSelection::Select(TObjArray* pTracks)
{
  /// create selection, array contains only pointers but does not own the objects
  /// object array needs to be deleted by caller
  if (!pTracks) return NULL;
  TObjArray* selectedTracks=new TObjArray;
  if (!selectedTracks) return NULL;
  TIter itrack(pTracks);
  TObject* pObj=NULL;
  while ((pObj=itrack())!=NULL) {
    AliVParticle* track=dynamic_cast<AliVParticle*>(pObj);
    if (!track) continue;
    if (!IsSelected(track)) continue;
    selectedTracks->Add(track);
  }
  return selectedTracks;
}

int AliDxHFEParticleSelection::CheckAndAdd(AliVParticle* /*p*/)
{
  /// check and add track to internal array
  /// TODO: check if needed
  return -ENOSYS;
}

bool AliDxHFEParticleSelection::IsSelected(AliVParticle* /*p*/)
{
  /// check particle if it passes the selection criteria
  /// childs can overload, by default all tracks are selected
  return true;
}

void AliDxHFEParticleSelection::AliDxHFEParticleSelection::Clear(Option_t * /*option*/)
{
  /// inherited from TObject: cleanup
}

void AliDxHFEParticleSelection::Print(Option_t */*option*/) const
{
  /// inherited from TObject: print info
}
 
void AliDxHFEParticleSelection::SaveAs(const char */*filename*/,Option_t */*option*/) const
{
  /// inherited from TObject: safe selection criteria
}
