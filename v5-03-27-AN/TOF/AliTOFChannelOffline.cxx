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
Revision 1.1  2007/06/27 16:55:24  zampolli
Task for TOF Calibration, optimized for pp runs

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

#include "AliTOFChannelOffline.h"
#include "AliLog.h"

ClassImp(AliTOFChannelOffline)

//________________________________________________________________
AliTOFChannelOffline::AliTOFChannelOffline()
{
  for(Int_t i=0;i<6;i++) fSlewPar[i]=0.;
}

//________________________________________________________________
AliTOFChannelOffline::AliTOFChannelOffline(Float_t* slewingPar)
{
  for(Int_t i = 0; i<6;i++) fSlewPar[i]=slewingPar[i];
}

//________________________________________________________________
AliTOFChannelOffline::~AliTOFChannelOffline()
{
  //dtor
  
}

//________________________________________________________________
AliTOFChannelOffline::AliTOFChannelOffline(const AliTOFChannelOffline& channel) :
  TObject(channel)
{
// copy constructor
  for(Int_t i = 0; i<6;i++) fSlewPar[i]=channel.GetSlewPar(i);
}


//________________________________________________________________
AliTOFChannelOffline &AliTOFChannelOffline::operator =(const AliTOFChannelOffline& channel)
{
// assignment operator

  if (this == &channel)
    return *this;

  TObject::operator=(channel);
  for(Int_t i = 0; i<6;i++) fSlewPar[i]=channel.GetSlewPar(i);
  return *this;
}

//
//________________________________________________________________
//virtual AliTOFChannelOffline::~AliTOFChannelOffline()/
//{

//}
//*/

//________________________________________________________________
void AliTOFChannelOffline::SetSlewPar(Float_t* slewingPar)
{
  if(slewingPar) for(Int_t i = 0; i<6;i++) fSlewPar[i]=slewingPar[i];
  else for(int t=0; t<6; t++) fSlewPar[t] = 0.;
}

