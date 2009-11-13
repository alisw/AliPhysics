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

#include <TNamed.h>
#include "AliFlowCumuConstants.h"  

ClassImp(AliFlowCumuConstants)
// Description: constants for the Cumulant flow analysis

AliFlowCumuConstants* AliFlowCumuConstants::fgPMasterConfig = NULL;

//______________________________________________________________________________
AliFlowCumuConstants::AliFlowCumuConstants():
  TNamed(),
  fQmax(11),
  fPmax(5),
  fQmax4(5),
  fPmax4(2),     
  fQmax6(7),
  fPmax6(3),   
  fQmax8(9),
  fPmax8(4), 
  fQmax16(17),
  fPmax16(8),     
  fFlow(2),  
  fMltpl(1),
  fBinWidth(0.1),   
  fR0(2.2),
  fPtMax(3.1),
  fPtMin(0.0),
  fOtherEquations(kFALSE)
{
  //def ctor
}

//______________________________________________________________________________
AliFlowCumuConstants::~AliFlowCumuConstants()
{
  //dtor
}

//______________________________________________________________________________
AliFlowCumuConstants* AliFlowCumuConstants::GetMaster()
{
  //get global master config
  if (!fgPMasterConfig) fgPMasterConfig = new AliFlowCumuConstants();
  return fgPMasterConfig;
}

