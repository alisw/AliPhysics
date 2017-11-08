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

////////////////////////////////////////////////////////////////////////////
///                                                                       //
///                                                                       //
/// Class container for identified charged hadron spectra: TOF            //
///                                                                       //
///                                                                       //
/// Authors:                                                              //
/// N. Jacazio,  nicolo.jacazio[AROBASe]bo.infn.it                        //
////////////////////////////////////////////////////////////////////////////

#define LOG_NO_INFO
#define LOG_NO_DEBUG
#include "AliAnTOFevent.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliLog.h"
#include <vector>

ClassImp(AliAnTOFevent);

//________________________________________________________________________
AliAnTOFevent::AliAnTOFevent()
    : fEvtMultBin(-1)
    , fAliAnTOFtracks()
{

  //
  // standard constructur which should be used
  //
  AliInfo("**** CONSTRUCTOR CALLED ****");
  Reset();
  AliInfo("**** END OF CONSTRUCTOR ****");
}

//________________________________________________________________________
AliAnTOFevent::~AliAnTOFevent()
{ //Destructor
  AliInfo("**** DESTRUCTOR CALLED ****");

  AliInfo("**** END OF DESTRUCTOR ****");
}

//________________________________________________________________________
AliAnTOFtrack* AliAnTOFevent::GetTrack(const Int_t i)
{
  if (i >= 0 && i < GetNtracks()) {
    AliInfoF("Returning track at potition i = %i", i);
    return &fAliAnTOFtracks.at(i);
  } else if (i >= 0) {
    AliInfoF("Creating new track track at potition i = %i", fNTracks);
    fAliAnTOFtracks.push_back(AliAnTOFtrack());
  }

  return &fAliAnTOFtracks.back();
}

//________________________________________________________________________
Int_t AliAnTOFevent::GetNtracks()
{
  return fAliAnTOFtracks.size();
}

//________________________________________________________________________
void AliAnTOFevent::Reset()
{
  fEvtMultBin = -1;
  fAliAnTOFtracks.clear();
}
