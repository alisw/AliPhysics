// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveTrackCounter.h"

#include "TEveManager.h"
#include "AliEveTrack.h"
#include "TEveGedEditor.h"

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

//______________________________________________________________________________
AliEveTrackCounter::AliEveTrackCounter(const Text_t* name, const Text_t* title) :
   TEveElement(),
   TNamed(name, title),

   fBadLineStyle (6),
   fClickAction  (kCA_ToggleTrack),
   fEventId      (-1),
   fAllTracks    (0),
   fGoodTracks   (0),
   fTrackLists   ()
{
   // Constructor.
   // Connects to global signal "AliEveTrack", "SecSelected(AliEveTrack*)".

   if (fgInstance == 0) fgInstance = this;
   TQObject::Connect("AliEveTrack", "SecSelected(AliEveTrack*)",
                     "AliEveTrackCounter", this, "DoTrackAction(AliEveTrack*)");
}

//______________________________________________________________________________
AliEveTrackCounter::~AliEveTrackCounter()
{
   // Destructor.
   // Disconnect from the global track signals.

   TQObject::Disconnect("AliEveTrack", "DoTrackAction(AliEveTrack*)");
   if (fgInstance == this) fgInstance = 0;
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveTrackCounter::Reset()
{
   // Reset internal track-counters and track-list.

   fAllTracks  = 0;
   fGoodTracks = 0;
   TIter next(&fTrackLists);
   TEveTrackList* tlist;
   while ((tlist = dynamic_cast<TEveTrackList*>(next())))
      tlist->DecDenyDestroy();
   fTrackLists.Clear("nodelete");
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
         } else {
            t->SetLineStyle(fBadLineStyle);
         }
         ++fAllTracks;
      }
      ++i;
   }
}

//______________________________________________________________________________
void AliEveTrackCounter::DoTrackAction(AliEveTrack* track)
{
   // Slot called when track is ctrl-clicked.
   //
   // No check is done if track actually belongs to one of the
   // registered track-lists.
   //
   // Probably it would be safer to copy good/bad tracks into special
   // sub-containers.
   // In this case one should also override RemoveElementLocal.

   static const TEveException eh("AliEveTrackCounter::DoTrackAction ");

   switch (fClickAction)
   {

      case kCA_PrintTrackInfo:
      {
         printf("AliEveTrack '%s'\n", track->GetObject(eh)->GetName());
         const TEveVector &v = track->GetVertex();
         const TEveVector &p = track->GetMomentum();;
         printf("  Vx=%f, Vy=%f, Vz=%f; Pt=%f, Pz=%f, phi=%f)\n",
                v.fX, v.fY, v.fZ, p.Perp(), p.fZ, TMath::RadToDeg()*p.Phi());
         printf("  <other information should be printed ... full AliESDtrack>\n");
         break;
      }

      case kCA_ToggleTrack:
      {
         if (track->GetLineStyle() == 1)
         {
            track->SetLineStyle(fBadLineStyle);
            --fGoodTracks;
         } else {
            track->SetLineStyle(1);
            ++fGoodTracks;
         }
         track->ElementChanged();
         gEve->Redraw3D();

         printf("AliEveTrackCounter::CountTrack All=%d, Good=%d, Bad=%d\n",
                fAllTracks, fGoodTracks, fAllTracks-fGoodTracks);

         if (gEve->GetEditor()->GetModel() == GetObject(eh))
            gEve->EditElement(this);

         break;
      }

   } // end switch fClickAction
}

/******************************************************************************/

//______________________________________________________________________________
void AliEveTrackCounter::OutputEventTracks(FILE* out)
{
   // Print good-track summary into a plain-text file by iteration
   // through all registered track-lists.
   // State of each track is determined by its line-style, it is
   // considered a good track if it's line style is solid.

   if (out == 0)
   {
      out = stdout;
      fprintf(out, "AliEveTrackCounter::FinalizeEvent()\n");
   }

   fprintf(out, "Event = %d  Ntracks = %d\n", fEventId, fGoodTracks);

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
            fprintf(out, " %2d: chg=%+2d  pt=%8.5f  eta=%+8.5f\n",
                    cnt, t->GetCharge(), t->GetMomentum().Perp(), t->GetMomentum().Eta());
         }
         ++i;
      }
   }
}
