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

#include "AliAlgSensTPC.h"
#include "AliAlgAux.h"
#include "AliLog.h"
ClassImp(AliAlgSensTPC)

using namespace AliAlgAux;
using namespace TMath;

//_________________________________________________________
AliAlgSensTPC::AliAlgSensTPC(const char* name,Int_t vid, Int_t iid, Int_t isec) 
  :AliAlgSens(name,vid,iid)
  ,fSector(isec)
{
  // def c-tor
}

//_________________________________________________________
AliAlgSensTPC::~AliAlgSensTPC()
{
  // d-tor
}

//__________________________________________________________________
void AliAlgSensTPC::SetTrackingFrame()
{
  // define tracking frame of the sensor: just rotation by sector angle
  fAlp = Sector2Alpha(fSector);
  fX = 0;
}

