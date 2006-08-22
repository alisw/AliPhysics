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
Revision 1.2  2006/02/13 17:22:26  arcelli
just Fixing Log info

Revision 1.1  2006/02/13 16:10:48  arcelli
Add classes for TOF Calibration (C.Zampolli)

author: Chiara Zampolli, zampolli@bo.infn.it
*/  

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for TOF calibration                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTOFChSim.h"

ClassImp(AliTOFChSim)

//________________________________________________________________
AliTOFChSim::AliTOFChSim():
fSlewedStatus(kFALSE),
fSpectrum(-1)
{}
//________________________________________________________________
AliTOFChSim::AliTOFChSim(const AliTOFChSim& channel) :
  TObject(channel),
  fSlewedStatus(kFALSE),
  fSpectrum(-1)
{
// copy constructor
  fSlewedStatus = channel.fSlewedStatus;
  fSpectrum= channel.fSpectrum;
}


//________________________________________________________________
AliTOFChSim &AliTOFChSim::operator =(const AliTOFChSim& channel)
{
// assignment operator
  fSlewedStatus = channel.IsSlewed();
  fSpectrum= channel.GetSpectrum();
  return *this;
}

