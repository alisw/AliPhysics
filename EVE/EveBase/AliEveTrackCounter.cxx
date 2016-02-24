// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTrackCounter.h"

#include "TEveManager.h"
#include "AliEveEventManager.h"
#include "AliEveTrack.h"
#include "AliEveTracklet.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliMultiplicity.h"

#include <TEveGedEditor.h>

#include <TFile.h>

//==============================================================================
// AliEveTrackCounter
//==============================================================================

//______________________________________________________________________________
//
// Provides event-based method for tagging of good / bad (or primary /
// secondary) tracks. A report can be written into a text file.
//
// AliEveTrack status is toggled by using secondary-selection / ctrl-click
// functionality of the GL viewer.
//
// Some of the functionality is implemented in AliEveTrackCounterEditor
// class.

ClassImp(AliEveTrackCounter)

//______________________________________________________________________________
AliEveTrackCounter* AliEveTrackCounter::fgInstance = 0;

Bool_t AliEveTrackCounter::IsActive()
{
  // Check if instance exists and is active.

  return fgInstance && fgInstance->fActive;
}

//______________________________________________________________________________
AliEveTrackCounter::AliEveTrackCounter(const Text_t* name, const Text_t* title) :
  TEveElement(),
  TNamed(name, title),

  fBadLineStyle (6),
  fClickAction  (kCA_ToggleTrack),
  fEventId      (-1),
  fAllTracks    (0), fGoodTracks   (0),
  fAllTracklets (0), fGoodTracklets(0),
  fTrackLists   (),  fTrackletLists(),
  fActive       (kFALSE)
{
  // Constructor.
  // Connects to global signal "AliEveTrack", "SecSelected(AliEveTrack*)".

  if (fgInstance == 0) fgInstance = this;

  TQObject::Connect("AliEveTrack", "SecSelectedTrack(AliEveTrack*)",
                    "AliEveTrackCounter", this, "DoTrackAction(AliEveTrack*)");
  TQObject::Connect("AliEveTracklet", "SecSelectedTracklet(AliEveTracklet*)",
                    "AliEveTrackCounter", this, "DoTrackletAction(AliEveTracklet*)");

  AliEveEventManager::Instance()->Connect("NewEventDataLoaded()", "AliEveTrackCounter", this, "Reset()");
}

//______________________________________________________________________________
AliEveTrackCounter::~AliEveTrackCounter()
{
  // Destructor.
  // Disconnect from the global track signals.

  AliEveEventManager::Instance()->Disconnect("NewEventDataLoaded()", this);

  TQObject::Disconnect("AliEveTrack", "DoTrackAction(AliEveTrack*)");
  TQObject::Disconnect("AliEveTracklet", "DoTrackletAction(AliEveTracklet*)");
  if (fgInstance == this) fgInstance = 0;
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveTrackCounter::Reset()
{
  // Reset internal track-counters and track-list.

  fAllTracks     = fGoodTracks    = 0;
  fAllTracklets  = fGoodTracklets = 0;
  {
    TIter next(&fTrackLists);
    TEveTrackList* tlist;
    while ((tlist = dynamic_cast<TEveTrackList*>(next())))
      tlist->DecDenyDestroy();
    //fTrackLists.Clear("nodelete");
        fTrackLists.Clear();
  }
  {
    TIter next(&fTrackletLists);
    TEveTrackList* tlist;
    while ((tlist = dynamic_cast<TEveTrackList*>(next())))
      tlist->DecDenyDestroy();
    //fTrackletLists.Clear("nodelete");
        fTrackletLists.Clear();
  }

  fEventId = AliEveEventManager::Instance()->GetEventId();
}

//______________________________________________________________________________
void AliEveTrackCounter::RegisterTracks(TEveTrackList* tlist, Bool_t goodTracks)
{
  // Register tracks from tlist and tlist itself.
  // If goodTracks is true, they are considered as primary/good
  // tracks.

  tlist->IncDenyDestroy();
  fTrackLists.Add(tlist);

  List_i i = tlist->BeginChildren();
  while (i != tlist->EndChildren())
  {
    AliEveTrack* t = dynamic_cast<AliEveTrack*>(*i);
    if (t != 0)
    {
      if (goodTracks)
      {
        ++fGoodTracks;
	t->GetESDTrack()->SetLabel(3);
      } else {
        t->SetLineStyle(fBadLineStyle);
	t->GetESDTrack()->SetLabel(0);
      }
      ++fAllTracks;
    }
    ++i;
  }
}

//______________________________________________________________________________
void AliEveTrackCounter::RegisterTracklets(TEveTrackList* tlist, Bool_t goodTracks)
{
  // Register tracklets from tlist and tlist itself.
  // If goodTracks is true, they are considered as primary/good
  // tracks.

  AliESDEvent     *esd = AliEveEventManager::Instance()->AssertESD();
  AliMultiplicity *mul = const_cast<AliMultiplicity*>(esd->GetMultiplicity());

  tlist->IncDenyDestroy();
  fTrackletLists.Add(tlist);

  List_i i = tlist->BeginChildren();
  while (i != tlist->EndChildren())
  {
    AliEveTracklet* t = dynamic_cast<AliEveTracklet*>(*i);
    if (t != 0)
    {
      if (goodTracks)
      {
	mul->SetLabel(t->GetIndex(), 0, 3);
        ++fGoodTracklets;
      } else {
	mul->SetLabel(t->GetIndex(), 0, 0);
        t->SetLineStyle(fBadLineStyle);
      }
      ++fAllTracklets;
    }
    ++i;
  }
}

//______________________________________________________________________________
void AliEveTrackCounter::DoTrackAction(AliEveTrack* track)
{
   // Slot called when track is secondary selected.
   //
   // No check is done if track actually belongs to one of the
   // registered track-lists.
   //
   // Probably it would be safer to copy good/bad tracks into special
   // sub-containers.
   // In this case one should also override RemoveElementLocal.

   static const TEveException eh("AliEveTrackCounter::DoTrackAction ");

   if (!fActive)
     return;

   switch (fClickAction)
   {

      case kCA_PrintTrackInfo:
      {
         printf("AliEveTrack '%s'\n", track->GetObject(eh)->GetName());
         const TEveVector &v = track->GetVertex();
         const TEveVector &p = track->GetMomentum();;
         printf("  Vx=%f, Vy=%f, Vz=%f; Pt=%f, Pz=%f, phi=%f)\n",
                v.fX, v.fY, v.fZ, p.Perp(), p.fZ, TMath::RadToDeg()*p.Phi());
         break;
      }

      case kCA_ToggleTrack:
      {
	 AliESDtrack *esdt = track->GetESDTrack();
         if (track->GetLineStyle() == 1)
         {
            track->SetLineStyle(fBadLineStyle);
	    esdt->SetLabel(esdt->GetLabel() & ~1);
            --fGoodTracks;
         } else {
            track->SetLineStyle(1);
	    esdt->SetLabel(esdt->GetLabel() | 1);
            ++fGoodTracks;
         }
         track->ElementChanged();
         gEve->Redraw3D();

         //printf("AliEveTrackCounter::DoTrackAction All=%d, Good=%d, Bad=%d\n",
         //       fAllTracks, fGoodTracks, fAllTracks-fGoodTracks);

         if (gEve->GetEditor()->GetModel() == GetObject(eh))
            gEve->EditElement(this);

         break;
      }

   } // end switch fClickAction
}

//______________________________________________________________________________
void AliEveTrackCounter::DoTrackletAction(AliEveTracklet* track)
{
   // Slot called when tracklet is secondary selected.
   //
   // No check is done if track actually belongs to one of the
   // registered track-lists.
   //
   // Probably it would be safer to copy good/bad tracks into special
   // sub-containers.
   // In this case one should also override RemoveElementLocal.

   static const TEveException eh("AliEveTrackCounter::DoTrackletAction ");

   if (!fActive)
     return;

   switch (fClickAction)
   {

      case kCA_PrintTrackInfo:
      {
         printf("AliEveTracklet '%s'\n", track->GetObject(eh)->GetName());
         const TEveVector &v = track->GetVertex();
         const TEveVector &p = track->GetMomentum();;
         printf("  Vx=%f, Vy=%f, Vz=%f; Pt=%f, Pz=%f, phi=%f)\n",
                v.fX, v.fY, v.fZ, p.Perp(), p.fZ, TMath::RadToDeg()*p.Phi());
         break;
      }

      case kCA_ToggleTrack:
      {
         AliESDEvent     *esd = AliEveEventManager::Instance()->AssertESD();
	 AliMultiplicity *mul = const_cast<AliMultiplicity*>(esd->GetMultiplicity());

         if (track->GetLineStyle() == 1)
         {
            track->SetLineStyle(fBadLineStyle);
	    mul->SetLabel(track->GetIndex(), 0, mul->GetLabel(track->GetIndex(), 0) & ~1);
            --fGoodTracklets;
         } else {
            track->SetLineStyle(1);
	    mul->SetLabel(track->GetIndex(), 0, mul->GetLabel(track->GetIndex(), 0) | 1);
            ++fGoodTracklets;
         }
         track->ElementChanged();
         gEve->Redraw3D();

         // printf("AliEveTrackCounter::DoTrackletAction All=%d, Good=%d, Bad=%d\n",
         //        fAllTracklets, fGoodTracklets, fAllTracklets-fGoodTracklets);

         if (gEve->GetEditor()->GetModel() == GetObject(eh))
            gEve->EditElement(this);

         break;
      }

   } // end switch fClickAction
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveTrackCounter::OutputEventTracks()
{
  // Print good-track summary into a plain-text file by iteration
  // through all registered track-lists.
  // State of each track is determined by its line-style, it is
  // considered a good track if it's line style is solid.

  {
    TFile *f = TFile::Open("scan_results.root", "UPDATE");

    AliESDEvent     *esd = AliEveEventManager::Instance()->AssertESD();
    TClonesArray    *trk = static_cast<TClonesArray*>   (esd->GetList()->FindObject("Tracks"));
    AliMultiplicity *mul = const_cast <AliMultiplicity*>(esd->GetMultiplicity());

    trk->Write(TString::Format("Tracks_%04d",    fEventId), kWriteDelete | kSingleKey);
    mul->Write(TString::Format("Tracklets_%04d", fEventId), kWriteDelete);

    esd->GetPrimaryVertexTracks()->Write(TString::Format("PrimVertTracks_%04d", fEventId), kWriteDelete);
    esd->GetPrimaryVertexTPC()   ->Write(TString::Format("PrimVertTPC_%04d",    fEventId), kWriteDelete);
    esd->GetPrimaryVertexSPD()   ->Write(TString::Format("PrimVertSPD_%04d",    fEventId), kWriteDelete);

    f->Close();
  }
}

//______________________________________________________________________________
void AliEveTrackCounter::PrintEventTracks()
{
  // Print good-track summary to stdout by iteration
  // through all registered track-lists.
  // State of each track is determined by its line-style, it is
  // considered a good track if it's line style is solid.

  FILE* out = stdout;

  fprintf(out, "AliEveTrackCounter::PrintEventTracks()\n");

  fprintf(out, "Event=%d\n", fEventId);
  fprintf(out, "GoodTracks=%d  AllTracks=%d\n", fGoodTracks, fAllTracks);

  {
    TIter tlists(&fTrackLists);
    TEveTrackList* tlist;
    Int_t cnt = 0;
    while ((tlist = (TEveTrackList*) tlists()) != 0)
    {
      List_i i = tlist->BeginChildren();
      while (i != tlist->EndChildren())
      {
        AliEveTrack* t = dynamic_cast<AliEveTrack*>(*i);
        if (t != 0 && t->GetLineStyle() == 1)
        {
          ++cnt;
          fprintf(out, " %2d: chg=%+2d  pt=%8.5f  eta=%+8.5f  phi=%+8.5f\n",
                  cnt, t->GetCharge(), t->GetMomentum().Perp(),
                  t->GetMomentum().Eta(), t->GetMomentum().Phi());
        }
        ++i;
      }
    }
  }

  fprintf(out, "GoodTracklets=%d  AllTracklets=%d\n", fGoodTracklets, fAllTracklets);
  {
    TIter tlists(&fTrackletLists);
    TEveTrackList* tlist;
    Int_t cnt = 0;
    while ((tlist = (TEveTrackList*) tlists()) != 0)
    {
      List_i i = tlist->BeginChildren();
      while (i != tlist->EndChildren())
      {
        AliEveTracklet* t = dynamic_cast<AliEveTracklet*>(*i);
        if (t != 0 && t->GetLineStyle() == 1)
        {
          ++cnt;
          fprintf(out, " %2d: theta=%+8.5f  eta=%+8.5f  phi=%+8.5f\n",
                  cnt, t->GetMomentum().Theta(), t->GetMomentum().Eta(), t->GetMomentum().Phi());
        }
        ++i;
      }
    }
  }
}
