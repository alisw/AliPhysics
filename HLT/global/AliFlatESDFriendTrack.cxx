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
 * >> Flat structure representing an ESDTrack <<
 *
 * To be used in the online and offline calibration schema.
 *
 * Class provides interface methods for 
 *   - Filling from AliESDtrack and AliExternalTrackParam, as well 
 *     as clusters from ESD friends (if requested)
 *   - HLT Filling to be added
 * 
 *
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli
 *
 **************************************************************************/


#include "AliFlatESDFriendTrack.h"
#include "AliExternalTrackParam.h"
#include "Riostream.h"

// _______________________________________________________________________________________________________
AliFlatESDFriendTrack::AliFlatESDFriendTrack() :
  AliVVfriendTrack()
{
  // Default constructor
}


// _______________________________________________________________________________________________________
AliFlatESDFriendTrack::AliFlatESDFriendTrack(AliFlatESDSpecialConstructorFlag f) :
  AliVVfriendTrack()
{
  //special constructor, used to restore the vtable pointer
  //uses the special dummy constructors of contained objects

  // the vtable pointer for this AliFlatESDFriend object is already reset when this constructor is called
  // we should only initialise vtable pointers for all contained objects

  if(f == AliFlatESDReinitialize){   
  }
  else{
    AliFlatESDFriendTrack();
  }
}
