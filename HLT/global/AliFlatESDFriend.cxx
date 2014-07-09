/* $Id$ */

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

/**
 * >> Flat structure representing an ESD friend <<
 *
 * To be used in the online and offline calibration schema.
 *
 * Class provides interface methods for 
 *   - Filling from AliESDfriend, but also from HLT 
 *   - Getter methods
 *
 * In the online case, the structure can be directly written into a shared 
 * memory, in the offline case, the size has to be estimated first.
 *
 * 
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli
 *
 * ************************************************************************
 **/

#include "AliFlatESDFriend.h"
#include "AliFlatESDFriendTrack.h"
#include "Riostream.h"

// _______________________________________________________________________________________________________
AliFlatESDFriend::AliFlatESDFriend() :
  fContentSize(0),
  fBitFlags(0),
  fNTracks(0),
  fNTrackEntries(0),
  fTrackTablePointer(0),
  fTracksPointer(0)
{
  // Default constructor
  Reset();
}

void AliFlatESDFriend::Reset() 
{
  fBitFlags = 0;
  fNTracks = 0;
  fNTrackEntries = 0; 
  fTrackTablePointer = 0;
  fTracksPointer = 0;
  for( int i=0; i<72; i++ ){
    fNclustersTPC[i]=0;
    fNclustersTPCused[i]=0;
  }
  // We set size of the fContent array such, that it reaches the end of the AliFlatESDFriend structure. 
  // To be sure that actual GetSize() is always >= size of the structure. 
  // First, it makes the memory alignment better. Second, just for a case..
  fContentSize = sizeof(AliFlatESDFriend) - (fContent - reinterpret_cast<const Byte_t*>(this));
  for( UInt_t i=0; i<fContentSize; i++ ) fContent[i]=0;
}
 
// _______________________________________________________________________________________________________
AliFlatESDFriend::AliFlatESDFriend(AliFlatESDSpecialConstructorFlag f) :
  fContentSize(),
  fBitFlags(),
  fNTracks(),
  fNTrackEntries(),
  fTrackTablePointer(),
  fTracksPointer()
{
  //special constructor, used to restore the vtable pointer
  //uses the special dummy constructors of contained objects

  // the vtable pointer for this AliFlatESDFriend object is already reset when this constructor is called
  // we should only initialise vtable pointers for all contained objects

  if(f == AliFlatESDReinitialize){   
    for( int i=0; i<fNTracks; i++ ){
      AliFlatESDFriendTrack  *tr = GetFlatTrackNonConst(i);
      if( tr ) tr->Reinitialize();
    }
  }
  else{
    AliFlatESDFriend();
  }
}

void AliFlatESDFriend::Ls() const
{
  cout<<"Flat ESD friend: "<<endl;
  cout<<"  N tracks: "<<fNTracks<<endl;
  cout<<"  N track entries: "<<fNTrackEntries<<endl;
}
