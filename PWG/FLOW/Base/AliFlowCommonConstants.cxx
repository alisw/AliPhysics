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
// AliFlowCommonConstants:
//
// Constants for the common histograms in the flow analysis
// semi "singleton" mode: apart from any objects instantiated via
// the constructor there may be a global Master object. get it with
// AliFlowCommonConstants::GetMaster()->GetMaster() static method.
//
// Author: Naomi van der Kolk (kolk@nikhef.nl)
// mod: Mikolaj Krzewicki, Nikhef (mikolaj.krzewicki@cern.ch)

/*
$Log$
*/ 

#include <TMath.h>
#include "AliFlowCommonConstants.h" 

ClassImp(AliFlowCommonConstants)

AliFlowCommonConstants* AliFlowCommonConstants::fgPMasterConfig = NULL;

//______________________________________________________________________________
AliFlowCommonConstants::AliFlowCommonConstants():
  TNamed(),
  fNbinsMult(10000),
  fNbinsPt(100),   
  fNbinsPhi(72),
  fNbinsEta(200),
  fNbinsQ(500),
  fMultMin(0.),            
  fMultMax(10000.),
  fPtMin(0.),	     
  fPtMax(10.),
  fPhiMin(0.),	     
  fPhiMax(TMath::TwoPi()),
  fEtaMin(-5.),	     
  fEtaMax(5.),	     
  fQMin(0.),	     
  fQMax(3.),
  fHistWeightvsPhiMin(0.),
  fHistWeightvsPhiMax(3.)
{
  //def ctor
}

//______________________________________________________________________________
AliFlowCommonConstants::~AliFlowCommonConstants()
{
  //dtor
}

//______________________________________________________________________________
AliFlowCommonConstants* AliFlowCommonConstants::GetMaster()
{
  //return pointer to master config object
  if (!fgPMasterConfig) fgPMasterConfig = new AliFlowCommonConstants();
  return fgPMasterConfig;
}
