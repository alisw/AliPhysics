/***************************************************************************
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
//-----------------------------------------------------//
//                                                     //
//  Date   : August 05 2003                            //
//                                                     //
//  Store cell/track info which is used to assign      //
//  the correct track number to a multiple hit cell    //
//                                                     //
//-----------------------------------------------------//

#include "Riostream.h"
#include "Rtypes.h"
#include "AliPMDcell.h"

ClassImp(AliPMDcell)

AliPMDcell::AliPMDcell()
{
  // Standard constructor
  fTrNumber = 0;
  fSMNumber = 0;
  fXpos     = 0;
  fYpos     = 0;
  fEdep     = 0.;
}

AliPMDcell::AliPMDcell(Int_t trnumber, Int_t smnumber, 
			 Int_t xpos, Int_t ypos, Float_t edep)
{
  // Constructor
  fTrNumber = trnumber;
  fSMNumber = smnumber;
  fXpos     = xpos;
  fYpos     = ypos;
  fEdep     = edep;

}

AliPMDcell::AliPMDcell(const AliPMDcell& source):TObject(source) {
  //Copy Constructor 
  if(&source == this) return;
  this->fTrNumber = source.fTrNumber;
  this->fSMNumber = source.fSMNumber;
  this->fXpos = source.fXpos;
  this->fYpos = source.fYpos;
  this->fEdep = source.fEdep;
  return;
}

AliPMDcell& AliPMDcell::operator=(const AliPMDcell& source) {
  //Copy Constructor 
  if(&source == this) return *this;
  this->fTrNumber = source.fTrNumber;
  this->fSMNumber = source.fSMNumber;
  this->fXpos = source.fXpos;
  this->fYpos = source.fYpos;
  this->fEdep = source.fEdep;
  return *this;
}

AliPMDcell::~AliPMDcell()
{
  // Default destructor
}

Int_t AliPMDcell::GetTrackNumber() const
{
  return fTrNumber;
}
Int_t AliPMDcell::GetSMNumber() const
{
  return fSMNumber;
}
Int_t AliPMDcell::GetX() const
{
  return fXpos;
}
Int_t AliPMDcell::GetY() const
{
  return fYpos;
}

Float_t AliPMDcell::GetEdep() const
{
  return fEdep;
}

