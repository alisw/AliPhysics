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
#include "TError.h"
#include <iostream>
#include <vector>

//________________________________________________________________________
AliAnTOFevent::AliAnTOFevent()
    : fEvtMultBin(-1)
    , fVtxX(-999)
    , fVtxY(-999)
    , fVtxZ(-999)
    , fAliAnTOFtracks()
{ // standard constructor which should be used

#ifndef LOG_NO_INFO
  ::Info("AliAnTOFevent::AliAnTOFevent", "**** CONSTRUCTOR CALLED ****");
#endif

  Reset();

#ifndef LOG_NO_INFO
  ::Info("AliAnTOFevent::AliAnTOFevent", "**** END OF CONSTRUCTOR ****");
#endif
}

//________________________________________________________________________
AliAnTOFevent::~AliAnTOFevent()
{ //Destructor
#ifndef LOG_NO_INFO
  ::Info("AliAnTOFevent::~AliAnTOFevent", "**** DESTRUCTOR CALLED ****");
#endif

#ifndef LOG_NO_INFO
  ::Info("AliAnTOFevent::~AliAnTOFevent", "**** END OF DESTRUCTOR ****");
#endif
}

//________________________________________________________________________
AliAnTOFtrack* AliAnTOFevent::GetTrack(const Int_t i)
{
  if (i >= 0 && i < GetNtracks()) {
#ifndef LOG_NO_INFO
    ::Info("AliAnTOFevent::GetTrack", "Returning track at potition i = %i", i);
#endif
    return &fAliAnTOFtracks.at(i);
  } else if (i >= 0) {
#ifndef LOG_NO_INFO
    ::Info("AliAnTOFevent::GetTrack", "Creating new track track at potition i = %i", GetNtracks());
#endif
    fAliAnTOFtracks.push_back(AliAnTOFtrack());
  }
  return &fAliAnTOFtracks.back();
}

//________________________________________________________________________
void AliAnTOFevent::AdoptVertex(const AliESDVertex* vtx)
{
  fVtxX = vtx->GetX();
  fVtxY = vtx->GetY();
  fVtxZ = vtx->GetZ();
}

//________________________________________________________________________
void AliAnTOFevent::Reset()
{
  fEvtMultBin = -1;
  fAliAnTOFtracks.clear();
  fVtxX = -999;
  fVtxY = -999;
  fVtxZ = -999;
}
