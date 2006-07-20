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

//
// This class stores the status of the Shuttle processing for a given run and a given detector
//
// This class stores the status of the processing, the number of retries and the timestamp of the last action
// The detector and run number are stored using the CDB framework
//
//

#include "AliShuttleStatus.h"

ClassImp(AliShuttleStatus)

AliShuttleStatus::AliShuttleStatus() : TObject(),
  fTimeStamp(0),
  fStatus(kInvalid),
  fCount(0)
{
  // default constructor
}

AliShuttleStatus::AliShuttleStatus(Status status) : TObject(),
  fTimeStamp(0),
  fStatus(status),
  fCount(1)
{
  // constructor

  fTimeStamp = time(0);
}

AliShuttleStatus::AliShuttleStatus(const AliShuttleStatus& c) : TObject(c)
{
  // copy constructor

  ((AliShuttleStatus &)c).Copy(*this);
}

AliShuttleStatus::~AliShuttleStatus()
{
  // destructor
}

AliShuttleStatus &AliShuttleStatus::operator=(const AliShuttleStatus &c)
{
  // assigment operator

  if (this != &c) 
    ((AliShuttleStatus &) c).Copy(*this);

  return *this;
}

void AliShuttleStatus::Copy(TObject& c) const
{
  // copy function

  AliShuttleStatus& target = (AliShuttleStatus &) c;

  target.fTimeStamp = fTimeStamp;
  target.fStatus = fStatus;
  target.fCount = fCount;
}

void AliShuttleStatus::SetStatus(Status status)
{
  // sets a new status, add the same time the timestamp is set to now

  fStatus = status;
  fTimeStamp = time(0);
}

const char* AliShuttleStatus::GetStatusName(Status status)
{
  // returns a name (string) of the status

  switch (status)
  {
    case kInvalid: return "Invalid";
    case kStarted: return "Started";
    case kDCSStarted: return "DCSStarted";
    case kDCSError: return "DCSError";
    case kPPStarted: return "PPStarted";
    case kPPError: return "PPError";
    case kPPDone: return "PPDone";
    case kDone: return "Done";
    case kFailed: return "Failed";
  }

  return 0;
}
