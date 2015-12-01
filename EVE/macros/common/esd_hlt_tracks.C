// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TObject.h>
#include <TClass.h>
#include <TFile.h>
#include <TTree.h>
#include <TEveManager.h>
#include <TEveTrack.h>
#include <TEveUtil.h>

#include <AliESDEvent.h>
#include <AliEveEventManager.h>

#include <EVE/macros/esd_tracks.C>
#endif

/**
 * Display ESD Tracks from the HLTesdTree in AliEVE.
 *
 * Usage:
 * <pre>
 *   alieve $ALICE_ROOT/EVE/macros/event_next.C \
 *          $ALICE_ROOT/EVE/macros/alieve_init.C \
 *          $ALICE_ROOT/EVE/macros/geom_simple.C \
 *          $ALICE_ROOT/EVE/macros/esd_hlt_tracks.C
 * </pre>
 * Display is changed to next event by executing event_next()
 * from the root prompt.
 * <pre>
 *   event_next(); esd_tracks(); esd_hlt_tracks();
 * </pre>
 *
 * @ingroup alihlt_tpc
 * @author Matthias.Richter@ift.uib.no
 * @date   2008-11-22
 */
TEveTrackList* esd_hlt_tracks()
{
  if (!TClass::GetClass("AliEveEventManager")) {
    Error("hlt_tpc_clusters.C", "EVE library not loaded, please start alieve correctly");
    return NULL;
  }

  AliEveEventManager* eveManager=AliEveEventManager::GetMaster();
  if (!eveManager) {
    Error("esd_hlt_tracks.C", "EVE manager not initialized");
    return NULL;
  }

  TEveUtil::LoadMacro("esd_tracks.C");
   
  int eventId=eveManager->GetEventId();
  TFile* esdFile=eveManager->GetESDFile();
  if (!esdFile) {
    Warning("esd_hlt_tracks.C", "can not get esd file from EVE manager");
    return NULL;
  }

  TObject* pObj=NULL;
  TTree* pHLTTree=NULL;
  esdFile->GetObject("HLTesdTree", pObj);
  if (!pObj || (pHLTTree=dynamic_cast<TTree*>(pObj))==NULL) {
    Info("esd_hlt_tracks.C", "no HLT ESD tree in ESD file");
    return NULL;
  }
  if (pHLTTree->GetEntries()<=eventId) {
    Warning("esd_hlt_tracks.C", "skiping event %d: out of range %lld", eventId, pHLTTree->GetEntries());
    return NULL;
  }

  AliESDEvent* esd=new AliESDEvent;
  esd->ReadFromTree(pHLTTree);
  pHLTTree->GetEntry(eventId);

  TEveTrackList* cont = new TEveTrackList("HLT ESD Tracks");
  cont->SetMainColor(kCyan+3);
  esd_track_propagator_setup(cont->GetPropagator(),
			     0.1*esd->GetMagneticField(), 520);

  eveManager->AddElement(cont);

  Int_t count = 0;
  for (Int_t n = 0; n < esd->GetNumberOfTracks(); ++n)
  {
    ++count;
    TEveTrack* track = esd_make_track(esd->GetTrack(n), cont);

    cont->AddElement(track);
  }
  cont->SetTitle(Form("N=%d", count));
  cont->MakeTracks();

  gEve->Redraw3D();

  return cont;
}
