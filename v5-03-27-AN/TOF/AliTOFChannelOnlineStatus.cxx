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

/*
$Log$
*/  

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for TOF Online calibration - defining channel status                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTOFChannelOnlineStatus.h"

ClassImp(AliTOFChannelOnlineStatus)

//________________________________________________________________
AliTOFChannelOnlineStatus::AliTOFChannelOnlineStatus():
fStatus(kTOFOnlineUnknown)
{
  //default constructor
}

//________________________________________________________________
AliTOFChannelOnlineStatus::AliTOFChannelOnlineStatus(UChar_t status):
fStatus(status)
{
  // constructor with status 
}

//________________________________________________________________
AliTOFChannelOnlineStatus::AliTOFChannelOnlineStatus(const AliTOFChannelOnlineStatus& channel) :
  TObject(channel),
  fStatus(kTOFOnlineUnknown)
{
// copy constructor

}


//________________________________________________________________
AliTOFChannelOnlineStatus &AliTOFChannelOnlineStatus::operator =(const AliTOFChannelOnlineStatus& channel)
{
// assignment operator

  if (this == &channel)
    return *this;

  TObject::operator=(channel);
  fStatus= channel.GetStatus();
  return *this;
}

//

