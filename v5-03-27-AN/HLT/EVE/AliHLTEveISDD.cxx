/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Svein Lindal <slindal@fys.uio.no   >                  *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   AliHLTEvePhos.cxx
/// @author Svein Lindal <slindal@fys.uio.no>
/// @brief  ISDD class for the HLT EVE display

#include "AliHLTEveISDD.h"
#include "TEvePointSet.h"

ClassImp(AliHLTEveISDD);

AliHLTEveISDD::AliHLTEveISDD() : 
AliHLTEveITS("ISDD")
{
  // Constructor.
}

AliHLTEveISDD::~AliHLTEveISDD()
{
  //Destructor
}

void AliHLTEveISDD::SetUpPointSet(TEvePointSet * ps) {
  //See header file for documentation
  ps->SetMainColor(kBlack);
  ps->SetMarkerStyle((Style_t)kFullDotMedium);
}
