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

// $Id: AliPhJMCTrackList.cxx,v 1.4 2008/05/08 13:44:46 djkim Exp $

////////////////////////////////////////////////////
//
//  \file AliPhJMCTrackList.cxx
//  \brief
//  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
//  \email: djkim@jyu.fi
//  \version $Revision: 1.4 $
//  \date $Date: 2008/05/08 13:44:46 $
//
//  Class containing a list (TClonesArray) of MC tracks
////////////////////////////////////////////////////

#include "AliPhJMCTrackList.h"

using namespace std;

ClassImp(AliPhJMCTrackList)

#define kNumTracks 1500

//______________________________________________________________________________

AliPhJMCTrackList::AliPhJMCTrackList(expName exp):fMcTrackList(0x0),fTracks(0){
  //constructor
  switch (exp){
   case kPHENIX:
       break;
   case kALICE:
       fMcTrackList = new TClonesArray("AliJMCTrack",kNumTracks);
       break;
   }
}

//______________________________________________________________________________

AliPhJMCTrackList::AliPhJMCTrackList():fMcTrackList(0x0),fTracks(0){
  //constructor
  fMcTrackList = new TClonesArray("AliJMCTrack",kNumTracks);
}

//______________________________________________________________________________

AliPhJMCTrackList::AliPhJMCTrackList(const AliPhJMCTrackList& a):
  TObject(a),
  fMcTrackList(new TClonesArray(*a.fMcTrackList)),
  fTracks(a.fTracks)
{
//copy constructor
}

//______________________________________________________________________________

AliPhJMCTrackList::~AliPhJMCTrackList(){
  //destructor
  fMcTrackList->Clear();
  delete fMcTrackList;
}

//______________________________________________________________________________

void AliPhJMCTrackList::Reset(){
  //reset list
  fMcTrackList->Clear();
  if(fTracks>kNumTracks){
    fMcTrackList->Expand(kNumTracks);
  }
  fTracks = 0;
}

//______________________________________________________________________________

int AliPhJMCTrackList::SetTClonesArraySize(const unsigned int ntrk){
  //set size of the list
  if(ntrk>kNumTracks){
    fMcTrackList->Expand(kNumTracks);
  }
  return ntrk;
}

//______________________________________________________________________________

void AliPhJMCTrackList::AddJMCTrack(const unsigned int itrk){
  //add Mc track to the list
  new((*fMcTrackList)[itrk]) AliJMCTrack();
}

//______________________________________________________________________________


AliJMCTrack* AliPhJMCTrackList::GetTrack(const unsigned int itrk){
  //get track from the list
  return (AliJMCTrack*)fMcTrackList->UncheckedAt(itrk);
}

//______________________________________________________________________________

AliPhJMCTrackList& AliPhJMCTrackList::operator=(const AliPhJMCTrackList&  list){
  //operator=
  if(this != &list){
    TObject::operator=(list);
    fMcTrackList=list.fMcTrackList;
    fTracks=list.fTracks;
  }

  return *this;
}



