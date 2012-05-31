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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Output of T0 DAQ		                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliT0Dqclass.h"

ClassImp(AliT0Dqclass)

//________________________________________________________________
AliT0Dqclass::AliT0Dqclass()
{
//  fHistMeanPed=0;
  Reset();
}

//________________________________________________________________
AliT0Dqclass::AliT0Dqclass(const char* name)
{
  // Constructor
  TString namst = "Calib_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
//  fHistMeanPed=0;
  Reset();
}


//________________________________________________________________
AliT0Dqclass::~AliT0Dqclass()
{
//  CleanHistos();
}

//________________________________________________________________
void AliT0Dqclass::Reset()
{
  // Reset
  memset(fTime,0,24*sizeof(Float_t));
  memset(fAmplitude,0,24*sizeof(Float_t));
}                                                                                       


//________________________________________________________________
void  AliT0Dqclass::Print(Option_t *) const
{

} 

