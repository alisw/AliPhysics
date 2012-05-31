// $Id$
// Author: Matevz Tadel 2009

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTrack.h"

#include "AliESDtrack.h"
#include "AliAODTrack.h"

#include <TROOT.h>
#include <TMath.h>

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
AliEveTrack::AliEveTrack(TParticle* t, Int_t label, TEveTrackPropagator* prop) :
  TEveTrack(t, label, prop)
{
  // Constructor.
}

//______________________________________________________________________________
AliEveTrack::AliEveTrack(TEveMCTrack*  t, TEveTrackPropagator* prop) :
  TEveTrack(t, prop)
{
  // Constructor.
}

//______________________________________________________________________________
AliEveTrack::AliEveTrack(TEveRecTrack* t, TEveTrackPropagator* prop) :
  TEveTrack(t, prop)
{
  // Constructor.
}

//______________________________________________________________________________
AliEveTrack::AliEveTrack(AliESDtrack* t, TEveTrackPropagator* prop) :
  TEveTrack()
{
  // Constructor.

  Double_t buf[3];
  t->GetXYZ(buf);    fV.Set(buf);
  t->GetPxPyPz(buf); fP.Set(buf);

  Double_t ep = t->GetP(), mc = t->GetMass();
  fBeta = ep/TMath::Sqrt(ep*ep + mc*mc);
  // fPdg = 0; // ??? Use PID ?
  fCharge = TMath::Nint(t->GetSign());
  
  fLabel = t->GetLabel();
  fIndex = t->GetID();
  // fStatus = (Int_t) t->GetStatus(); // RRRR Uncomment for root-5.26.

  SetPropagator(prop);
}

//______________________________________________________________________________
AliEveTrack::AliEveTrack(AliAODTrack* t, TEveTrackPropagator* prop) :
  TEveTrack()
{
  // Constructor.

  Double_t buf[3];

  t->GetXYZ(buf); fV.Set(buf);
  t->PxPyPz(buf); fP.Set(buf);

  // fBeta = 0; // Unknown, no mass function
  // fPdg = 0;  // ??? Use PID ?
  fCharge= t->Charge();
  
  fLabel = t->GetLabel();
  fIndex = t->GetID();
  // fStatus = (Int_t) t->GetStatus(); // RRRR Uncomment for root-5.26.

  SetPropagator(prop);
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

//______________________________________________________________________________
void AliEveTrack::SetStartParams(const AliExternalTrackParam* tp)
{
  // Set the initial vertex / momentum of eve track from 'tp'.

  Double_t buf[3];

  tp->GetXYZ(buf); fV.Set(buf);
  tp->PxPyPz(buf); fP.Set(buf);
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
void AliEveTrack::ImportClustersFromLabel()
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
    throw kEH + "index not set.";

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
      throw kEH + "label not set.";

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
    throw kEH + "label not set.";

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
    throw kEH + "label not set.";

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

//______________________________________________________________________________
void AliEveTrack::SecSelected(TEveTrack* track)
{
  // Emits "SecSelected(TEveTrack*)" signal.
  // Called from TEveTrackGL on secondary-selection.

  Emit("SecSelected(TEveTrack*)", (Long_t)track);
  SecSelectedTrack((AliEveTrack*) track);
}

//______________________________________________________________________________
void AliEveTrack::SecSelectedTrack(AliEveTrack* track)
{
  // Emits "SecSelectedTrack(AliEveTrack*)" signal.

  Emit("SecSelectedTrack(AliEveTrack*)", (Long_t)track);
}

//______________________________________________________________________________
AliESDtrack* AliEveTrack::GetESDTrack() const
{
  // Return source object dyn-casted to AliESDtrack.

  return dynamic_cast<AliESDtrack*>(GetSourceObject());
}

//______________________________________________________________________________
AliAODTrack* AliEveTrack::GetAODTrack() const
{
  // Return source object dyn-casted to AliAODTrack.

  return dynamic_cast<AliAODTrack*>(GetSourceObject());
}
