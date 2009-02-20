// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveTrackCounter_H
#define AliEveTrackCounter_H

#include "TEveElement.h"
#include "TList.h"

class AliEveTrack;
class TEveTrackList;

//______________________________________________________________________________
// Short description of AliEveTrackCounter
//

class AliEveTrackCounter :  public TEveElement, public TNamed
{
  friend class AliEveTrackCounterEditor;

public:
  enum EClickAction_e { kCA_PrintTrackInfo, kCA_ToggleTrack };

  AliEveTrackCounter(const Text_t* name="AliEveTrackCounter", const Text_t* title="");
  virtual ~AliEveTrackCounter();

  Int_t GetEventId() const { return fEventId; }
  void  SetEventId(Int_t id) { fEventId = id; }

  void Reset();

  void RegisterTracks(TEveTrackList* tlist, Bool_t goodTracks);

  void DoTrackAction(AliEveTrack* track);

  Int_t GetClickAction() const  { return fClickAction; }
  void  SetClickAction(Int_t a) { fClickAction = a; }

  void OutputEventTracks(FILE* out=0);

  static AliEveTrackCounter* fgInstance;

protected:
  Int_t fBadLineStyle;  // TEveLine-style used for secondary/bad tracks.
  Int_t fClickAction;   // Action to take when a track is ctrl-clicked.

  Int_t fEventId;       // Current event-id.

  Int_t fAllTracks;     // Counter of all tracks.
  Int_t fGoodTracks;    // Counter of good tracks.

  TList fTrackLists;    // List of TrackLists registered for management.

private:
  AliEveTrackCounter(const AliEveTrackCounter&);            // Not implemented
  AliEveTrackCounter& operator=(const AliEveTrackCounter&); // Not implemented

  ClassDef(AliEveTrackCounter, 0); // Class for selection of good/primary tracks with basic processing functionality.
};

#endif
