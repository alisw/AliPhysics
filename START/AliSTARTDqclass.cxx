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

#include "AliSTARTDqclass.h"

ClassImp(AliSTARTDqclass)

//________________________________________________________________
AliSTARTDqclass::AliSTARTDqclass()
{
//  fHistMeanPed=0;
  Reset();
}

//________________________________________________________________
AliSTARTDqclass::AliSTARTDqclass(const char* name)
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
AliSTARTDqclass::~AliSTARTDqclass()
{
//  CleanHistos();
}

//________________________________________________________________
void AliSTARTDqclass::Reset()
{
  // Reset
  memset(fTime,0,24*sizeof(Float_t));
  memset(fAmplitude,0,24*sizeof(Float_t));
}                                                                                       


//________________________________________________________________
void  AliSTARTDqclass::Print(Option_t *) const
{

} 

