// $Id$
// Author: Matevz Tadel 2009

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTrack.h"

#include <TROOT.h>

//______________________________________________________________________________
// Full description of AliEveTrack
//

ClassImp(AliEveTrack)

//______________________________________________________________________________
AliEveTrack::AliEveTrack() :
  TEveTrack()
{
  // Constructor.
}

//______________________________________________________________________________
AliEveTrack::AliEveTrack(TParticle* t, Int_t label, TEveTrackPropagator* rs) :
  TEveTrack(t, label, rs)
{
  // Constructor.
}

//______________________________________________________________________________
AliEveTrack::AliEveTrack(TEveMCTrack*  t, TEveTrackPropagator* rs) :
  TEveTrack(t, rs)
{
  // Constructor.
}

//______________________________________________________________________________
AliEveTrack::AliEveTrack(TEveRecTrack* t, TEveTrackPropagator* rs) :
  TEveTrack(t, rs)
{
  // Constructor.
}

//______________________________________________________________________________
AliEveTrack::AliEveTrack(const AliEveTrack& t) :
  TEveTrack(t)
{
  // Copy constructor.
}

//______________________________________________________________________________
AliEveTrack::~AliEveTrack()
{
  // Destructor.
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveTrack::ImportHits()
{
  // Import hits with same label as the track.
  // Uses macro "hits_from_label.C".

  TEveUtil::LoadMacro("hits_from_label.C");
  gROOT->ProcessLine(Form("hits_from_label(%d, (TEveElement*)%p);",
                          fLabel, this));
}

//______________________________________________________________________________
void AliEveTrack::ImportClusters()
{
  // Import clusters with same label as the track.
  // Uses macro "clusters_from_label.C".

  TEveUtil::LoadMacro("clusters_from_label.C");
  gROOT->ProcessLine(Form("clusters_from_label(%d, (TEveElement*)%p);",
                          fLabel, this));
}

//______________________________________________________________________________
void AliEveTrack::ImportClustersFromIndex()
{
  // Import clusters marked with same reconstructed track index as the track.
  // Uses macro "clusters_from_index.C".

  static const TEveException kEH("AliEveTrack::ImportClustersFromIndex ");

  if (fIndex == kMinInt)
    throw(kEH + "index not set.");

  TEveUtil::LoadMacro("clusters_from_index.C");
  gROOT->ProcessLine(Form("clusters_from_index(%d, (TEveElement*)%p);",
                          fIndex, this));
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveTrack::ImportKine()
{
   // Import kinematics of the track's label recursively.
   // Uses macro "kine_tracks.C".

   static const TEveException kEH("AliEveTrack::ImportKine ");

   if (fLabel == kMinInt)
      throw(kEH + "label not set.");

   Int_t label;
   if (fLabel < 0) {
      Warning(kEH, "label negative, taking absolute value.");
      label = -fLabel;
   } else {
      label = fLabel;
   }

   TEveUtil::LoadMacro("kine_tracks.C");
   gROOT->ProcessLine(Form("kine_track(%d, kTRUE, kTRUE, kTRUE, kTRUE, (TEveElement*)%p);",
                           label, this));

}

//______________________________________________________________________________
void AliEveTrack::ImportKineWithArgs(Bool_t importMother, Bool_t importDaugters,
                                     Bool_t colorPdg,     Bool_t recurse)
{
  // Import kinematics of the track's label. Arguments steer the
  // import process:
  //   importMother     import particle with track's label
  //   importDaugters   import direct daughters of label
  //   colorPdg         color kinematics by PDG code
  //   recurse          recursive import of daughters' daughters
  // Uses macro "kine_tracks.C".

  static const TEveException kEH("AliEveTrack::ImportKineWithArgs ");

  if (fLabel == kMinInt)
    throw(kEH + "label not set.");

  Int_t label;
  if (fLabel < 0) {
    Warning(kEH, "label negative, taking absolute value.");
    label = -fLabel;
  } else {
    label = fLabel;
  }

  TEveUtil::LoadMacro("kine_tracks.C");
  gROOT->ProcessLine(Form("kine_track(%d, %d, %d, %d, %d, (TEveElement*)%p);",
                          label, importMother, importDaugters, colorPdg, recurse, this));
}

//______________________________________________________________________________
void AliEveTrack::PrintKineStack()
{
  // Print kinematics pertaining to track's label.
  // Uses macro "print_kine_from_label.C".

  static const TEveException kEH("AliEveTrack::PrintKineStack ");

  if (fLabel == kMinInt)
    throw(kEH + "label not set.");

  Int_t label;
  if (fLabel < 0) {
    Warning(kEH, "label negative, taking absolute value.");
    label = -fLabel;
  } else {
    label = fLabel;
  }

  TEveUtil::LoadMacro("print_kine_from_label.C");
  gROOT->ProcessLine(Form("print_kine_from_label(%d);", label));
}
