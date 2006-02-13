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

/*$Log$
author: Chiara Zampolli, zampolli@bo.infn.it
 */  

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for TOF calibration                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTOFChannel.h"

ClassImp(AliTOFChannel)

//________________________________________________________________
AliTOFChannel::AliTOFChannel():
fStatus(kFALSE),
fDelay(0)
{
  for(Int_t i=0;i<6;i++) fSlewPar[i]=0.;
}

//________________________________________________________________
AliTOFChannel::AliTOFChannel(Bool_t status, Float_t delay, Float_t* slewingPar):
fStatus(status),
fDelay(delay)
{
  for(Int_t i = 0; i<6;i++) fSlewPar[i]=slewingPar[i];
}

//________________________________________________________________
AliTOFChannel::AliTOFChannel(const AliTOFChannel& channel) :
  TObject(channel)
{
// copy constructor

}


//________________________________________________________________
AliTOFChannel &AliTOFChannel::operator =(const AliTOFChannel& channel)
{
// assignment operator
  fStatus= channel.GetStatus();
  fDelay=channel.GetDelay();
  return *this;
}

//
//________________________________________________________________
//virtual AliTOFChannel::~AliTOFChannel()/
//{

//}
//*/

//________________________________________________________________
void AliTOFChannel::SetSlewPar(Float_t* slewingPar)
{
  if(slewingPar) for(Int_t i = 0; i<6;i++) fSlewPar[i]=slewingPar[i];
  else for(int t=0; t<6; t++) fSlewPar[t] = 0.;
}



