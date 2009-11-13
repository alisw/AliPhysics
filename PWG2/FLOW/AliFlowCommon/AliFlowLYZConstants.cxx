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
*/ 

#include <TNamed.h> 
#include "AliFlowLYZConstants.h" 

ClassImp(AliFlowLYZConstants)
// AliFlowLYZConstant:
// Description: constants for the LYZ flow analysis
// Author: Naomi van der Kolk (kolk@nikhef.nl)
// modified: Mikolaj Krzewicki, Nikhef (mikolaj.krzewicki@cern.ch)

AliFlowLYZConstants* AliFlowLYZConstants::fgPMasterConfig = NULL;

//______________________________________________________________________________
AliFlowLYZConstants::AliFlowLYZConstants():
  TNamed(),
  fNtheta(5),
  fNbins(1200),
  fMaxSUM(120.),
  fMaxPROD(1.)
{
  //def ctor
}

//______________________________________________________________________________
AliFlowLYZConstants::~AliFlowLYZConstants()
{
  //dtor
}

//______________________________________________________________________________
AliFlowLYZConstants* AliFlowLYZConstants::GetMaster()
{
  //get global master config
  if (!fgPMasterConfig) fgPMasterConfig = new AliFlowLYZConstants();
  return fgPMasterConfig;
}
