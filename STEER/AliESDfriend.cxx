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

//-------------------------------------------------------------------------
//               Implementation of the AliESDfriend class
//  This class contains some additional to the ESD information like
//  the clusters associated to tracks.
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-------------------------------------------------------------------------

#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"
#include "AliESD.h"

ClassImp(AliESDfriend)

AliESDfriend::AliESDfriend(): TObject(), fTracks("AliESDfriendTrack",15000)
{
 //
 // Default constructor
 //
}

AliESDfriend::AliESDfriend(const AliESDfriend &f):TObject(f),fTracks(f.fTracks)
{
 //
 // Copy constructor
 //
}

AliESDfriend::AliESDfriend(const AliESD &event): TObject(event),
fTracks("AliESDfriendTrack",event.GetNumberOfTracks()) {
  //
  // Extracts the additional info from the ESD
  //
  Int_t ntrk=event.GetNumberOfTracks();

  for (Int_t i=0; i<ntrk; i++) {
    const AliESDtrack *t=event.GetTrack(i);
    new (fTracks[fTracks.GetEntriesFast()]) AliESDfriendTrack(*t); 
  }
}

AliESDfriend::~AliESDfriend() {
  //
  // Destructor
  //
  fTracks.Delete();
}
