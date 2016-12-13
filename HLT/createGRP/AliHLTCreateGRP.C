// **************************************************************************
// This file is property of and copyright by the ALICE HLT Project          *
// ALICE Experiment at CERN, All rights reserved.                           *
//                                                                          *
// Primary Authors: David Rohr <drohr@kip.uni-heidelberg.de>                *
//                  for The ALICE HLT Project.                              *
//                                                                          *
// Permission to use, copy, modify and distribute this software and its     *
// documentation strictly for non-commercial purposes is hereby granted     *
// without fee, provided that the above copyright notice appears in all     *
// copies and that both the copyright notice and this permission notice     *
// appear in the supporting documentation. The authors make no claims       *
// about the suitability of this software for any purpose. It is            *
// provided "as is" without express or implied warranty.                    *
//                                                                          *
//***************************************************************************

#include <stdio.h>

Int_t AliHLTCreateGRP(Int_t runNumber, TString detectorList, TString beamType, TString runType, UInt_t startShift, UInt_t endShift, Int_t defaults)
{
	printf("Running ALICE HLT GRP Creation\n");
	gSystem->Load("libAliHLTCreateGRP");
	
	AliHLTCreateGRP* GrpBuilder = new AliHLTCreateGRP;
	
	int retVal = GrpBuilder->CreateGRP(runNumber, detectorList, beamType, runType, startShift, endShift, defaults);
	
	delete GrpBuilder;
	return(retVal == 0 ? 2 : 100 + retVal);
}
