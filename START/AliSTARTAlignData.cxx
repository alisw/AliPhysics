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
// class for START calibrationalignment                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliSTARTAlignData.h"

ClassImp(AliSTARTAlignData)

//________________________________________________________________
AliSTARTAlignData::AliSTARTAlignData()
{
  // Default constructor
  Reset();
}

//________________________________________________________________
AliSTARTAlignData::AliSTARTAlignData(const char* name)
{
  // Constructor
  TString namst = "Align_";
  namst += name;
  SetName(namst.Data());
  SetTitle(namst.Data());
  Reset();
}

//________________________________________________________________
AliSTARTAlignData::AliSTARTAlignData(const AliSTARTAlignData& alignda) :
  TNamed(alignda)
{
  // copy constructor
  SetName(alignda.GetName());
  SetTitle(alignda.GetName());
  Reset();
  fSTARTzPosition[0] = alignda.GetZposition(0);
  fSTARTzPosition[1] = alignda.GetZposition(1);
}

//________________________________________________________________
AliSTARTAlignData &AliSTARTAlignData::operator =(const AliSTARTAlignData& alignda)
{
  // assignment operator
  SetName(alignda.GetName());
  SetTitle(alignda.GetName());
  Reset();
  fSTARTzPosition[0] = alignda.GetZposition(0);
  fSTARTzPosition[1] = alignda.GetZposition(1);
  return *this;
}

//________________________________________________________________
AliSTARTAlignData::~AliSTARTAlignData()
{
  // Destructor
}

//________________________________________________________________
void AliSTARTAlignData::Reset()
{
  // Set all pedestals to 0 and all ADC channels to 1
  memset(fSTARTzPosition,0,2*sizeof(Float_t));

}

//________________________________________________________________
void  AliSTARTAlignData::Print(Option_t *option) const
{
  // Print tables of pedestals and ADC channels

  printf("START aignment data:\n");
  printf("Z(A) = %f.2 cm, Z(C) = %f.2 cm\n",
	 fSTARTzPosition[0],fSTARTzPosition[1]);

}
