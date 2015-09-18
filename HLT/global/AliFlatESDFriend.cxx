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
#include "AliESDfriend.h"

// _______________________________________________________________________________________________________


using std::cout;
using std::endl;



void AliFlatESDFriend::Ls() const
{
  cout<<"Flat ESD friend: "<<endl;
  cout<<"  N tracks: "<<fNTracks<<endl;
  cout<<"  N track entries: "<<fNTrackEntries<<endl;
}



// _______________________________________________________________________________________________________
  ULong64_t AliFlatESDFriend::EstimateSize(AliESDfriend* esdFriend) 
{
  // Estimate upper limit of the object size
  // -> Added objects have to be added here as well
  if(esdFriend == NULL) return 0;
  ULong64_t size = sizeof(AliFlatESDFriend);
  size+= AliFlatESDVZEROFriend::GetSize();

  // one Long64_t per track for tracks table
  size += esdFriend->GetNumberOfTracks() *  (AliFlatESDFriendTrack::EstimateSize() + sizeof(Long64_t) );
  return size;
}

Int_t AliFlatESDFriend::SetVZEROFriend( const AliESDVZEROfriend *vzero, size_t allocatedVZEROMemory )
{
  // fill VZERO info
  fVZEROFriendPointer = -1;
  if( !vzero ) return 0;
  if( allocatedVZEROMemory < sizeof(AliFlatESDVZEROFriend) ) return -1;
  fVZEROFriendPointer = fContentSize;
  AliFlatESDVZEROFriend *flatVZERO = reinterpret_cast<AliFlatESDVZEROFriend*> (fContent + fContentSize);
  new (flatVZERO) AliFlatESDVZEROFriend ;
  flatVZERO->SetFromESDVZEROfriend( *vzero );
  fContentSize += flatVZERO->GetSize();
  return 0;
}

// _______________________________________________________________________________________________________

Int_t AliFlatESDFriend::SetFromESDfriend( const size_t allocatedMemorySize, const AliESDfriend *esdFriend )
{
  // Fill flat ESD friend from ALiESDfriend
 
  if( allocatedMemorySize < sizeof(AliFlatESDFriend ) ) return -1;

  Reset();
  
  if( !esdFriend ) return 0;
  
  Int_t err = 0;
  size_t freeSpace = allocatedMemorySize - GetSize();

  // fill event info
  {
    SetSkipBit( esdFriend->TestSkipBit() );
    for( int iSector=0; iSector<72; iSector++ ){
      SetNclustersTPC( iSector, esdFriend->GetNclustersTPC(iSector) );
      SetNclustersTPCused( iSector, esdFriend->GetNclustersTPCused(iSector) );
    }
  }

  // fill VZERO info
  {   
    const AliESDVZEROfriend *esdVZEROfriend = esdFriend->GetVZEROfriendConst();
    if (esdVZEROfriend) {
      err = SetVZEROFriend( esdVZEROfriend, freeSpace );      
      freeSpace = allocatedMemorySize - GetSize();
    }
  }
   
  if( err!=0 ) return err;

  // fill track friends
  {
   size_t trackSize = 0;
   int nTracks = 0;
   int nTrackEntries = 0;
   Long64_t *table = NULL;
   AliFlatESDFriendTrack *flatTrack = NULL;
   err = SetTracksStart( flatTrack, table, esdFriend->GetNumberOfTracks(), freeSpace );
   if( err!=0 ) return err;
   freeSpace = allocatedMemorySize - GetSize();
   
   for (Int_t idxTrack = 0; idxTrack < esdFriend->GetNumberOfTracks(); ++idxTrack) {
     const AliESDfriendTrack *esdTrack = esdFriend->GetTrack(idxTrack);
     table[idxTrack] = -1;
     if (esdTrack) {
       table[idxTrack] = trackSize;
       if( freeSpace<flatTrack->EstimateSize() ) return -1;
       new (flatTrack) AliFlatESDFriendTrack;       
       if( flatTrack->SetFromESDfriendTrack( esdTrack, freeSpace ) ) return -1;

       trackSize += flatTrack->GetSize();
       freeSpace -= flatTrack->GetSize();
       nTrackEntries++;
       flatTrack = flatTrack->GetNextTrackNonConst();
     }
     nTracks++;
    }
   
   SetTracksEnd( nTracks, nTrackEntries, trackSize );
  }

  return 0;
}

