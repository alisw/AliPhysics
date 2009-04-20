/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */


//-------------------------------------------------------------------------
//                          Class AliMixedEvent
// VEvent which is the container of several VEvents 
// Use Case: Event Mixing     
// Origin: Andreas Morsch, CERN, Andreas.Morsch@cern.ch 
//-------------------------------------------------------------------------


#include "AliMixedEvent.h"
#include <TMath.h>

ClassImp(AliMixedEvent)


AliMixedEvent::AliMixedEvent() :
    AliVEvent(),
    fEventList(), 
    fNEvents(0),       
    fNumberOfTracks(0),
    fNTracksCumul(0)
{
    // Default constructor
}


AliMixedEvent::AliMixedEvent(const AliMixedEvent& Evnt) :
    AliVEvent(Evnt),
    fEventList(), 
    fNEvents(0),       
    fNumberOfTracks(0),
    fNTracksCumul(0)
{ } // Copy constructor

AliMixedEvent& AliMixedEvent::operator=(const AliMixedEvent& vEvnt)
{ if (this!=&vEvnt) { 
    AliVEvent::operator=(vEvnt); 
  }
  
  return *this; 
}


void AliMixedEvent::AddEvent(AliVEvent* evt)
{
    // Add a new event to the list
    fEventList.Add(evt);
}


void AliMixedEvent::Init()
{
    // Initialize meta information
    fNEvents = fEventList.GetEntries();
    fNTracksCumul = new Int_t[fNEvents];
    fNumberOfTracks = 0;
    TIter next(&fEventList);
    AliVEvent* event;
    Int_t iev = 0;
    
    while((event = (AliVEvent*)next())) {
	fNTracksCumul[iev++] = fNumberOfTracks;
	fNumberOfTracks += (event->GetNumberOfTracks());
    }
}


AliVParticle* AliMixedEvent::GetTrack(Int_t i) const
{
    // Return track # i
    Int_t iEv  = TMath::BinarySearch(fNEvents, fNTracksCumul, i);
    while((iEv < (fNEvents - 1)) && (fNTracksCumul[iEv] == fNTracksCumul[iEv+1])) {iEv++;}
    
    Int_t irel = i - fNTracksCumul[iEv];
    AliVEvent* evt = (AliVEvent*) (fEventList.At(iEv));
    return (evt->GetTrack(irel));
}


void AliMixedEvent::Reset()
{
    // Reset the event
    fEventList.Clear();
    fNEvents = 0;
    fNumberOfTracks = 0;
    if (fNTracksCumul) {
	delete[]  fNTracksCumul;
	fNTracksCumul = 0;
    }
}

Int_t AliMixedEvent::EventIndex(Int_t itrack)
{
  // Return the event index for track #itrack
  return  TMath::BinarySearch(fNEvents, fNTracksCumul, itrack);
}
