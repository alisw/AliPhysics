/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notifce   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

////////////////////////////////////////////////////
//
//  \file AliPhJTrackList.cxx
//  \brief
//  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
//  \email: djkim@jyu.fi
//  \version $Revision: 1.4 $
//  \date $Date: 2008/05/08 13:44:46 $
//
//  Class contains a list (TClonesArray) of tracks  
/////////////////////////////////////////////////////

// $Id: AliPhJTrackList.cxx,v 1.4 2008/05/08 13:44:46 djkim Exp $

#include "AliJTrack.h"
#include "AliPhJTrackList.h"

using namespace std;

ClassImp(AliPhJTrackList)

#define kNumTracks 1500

//______________________________________________________________________________

AliPhJTrackList::AliPhJTrackList(expName exp):fTrackList(0x0),fTracks(0){
  //constructor
  switch (exp){
   case kPHENIX:
       break;
   case kALICE:
       fTrackList = new TClonesArray("AliJTrack",kNumTracks);
       break;
   }
}

//______________________________________________________________________________

AliPhJTrackList::AliPhJTrackList():fTrackList(0x0),fTracks(0){
  //constructor
  fTrackList = new TClonesArray("AliJTrack",kNumTracks);
}

//______________________________________________________________________________

AliPhJTrackList::AliPhJTrackList(const AliPhJTrackList& a):
  TObject(a),
  fTrackList(new TClonesArray(*a.fTrackList)),
  fTracks(a.fTracks)
{
  //copy constructor
}

//______________________________________________________________________________


AliPhJTrackList::~AliPhJTrackList(){
  //destructor
  fTrackList->Clear();
  delete fTrackList;
}

//______________________________________________________________________________

void AliPhJTrackList::Reset(){
  //reset list
  fTrackList->Clear();
  if(fTracks>kNumTracks){
    fTrackList->Expand(kNumTracks);
  }
  fTracks = 0;
}

//______________________________________________________________________________

int AliPhJTrackList::SetTClonesArraySize(const unsigned int ntrk){
  //set list size
  if(ntrk>kNumTracks){
    fTrackList->Expand(kNumTracks);
  }
  return ntrk;
}

//______________________________________________________________________________

void AliPhJTrackList::AddAliJTrack(const unsigned int itrk){
  //add a new track to the list
  new((*fTrackList)[itrk]) AliJTrack();
}

//______________________________________________________________________________

AliPhJBaseTrack* AliPhJTrackList::GetTrack(const unsigned int itrk){
  //retrieve a track from the list
  return (AliPhJBaseTrack*)fTrackList->UncheckedAt(itrk);
}

//______________________________________________________________________________

AliJTrack* AliPhJTrackList::GetAliJTrack(const unsigned int itrk){
   // ALICE getter
   return (AliJTrack*)fTrackList->UncheckedAt(itrk);
}


//______________________________________________________________________________

AliPhJTrackList& AliPhJTrackList::operator=(const AliPhJTrackList& list){
  //operator=
  if(this != &list){
    TObject::operator=(list);
    fTrackList = list.fTrackList;
    fTracks    = list.fTracks;
  }

  return *this;
}

