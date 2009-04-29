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
 

// Realisation of an AliVEventPool via
// on the flight (OTF) generation of the bin using AliTagAnalysis.
// Author Andreas Morsch
// andreas.morsch@cern.ch

#include "AliEventPoolOTF.h"

#include "AliRunTagCuts.h"
#include "AliLHCTagCuts.h"
#include "AliDetectorTagCuts.h"
#include "AliEventTagCuts.h"
#include "AliTagAnalysis.h"

#include <TMath.h>

ClassImp(AliEventPoolOTF)


////////////////////////////////////////////////////////////////////////

AliEventPoolOTF::AliEventPoolOTF():
    AliVEventPool(),
    fTagAnalysis(0),
    fRunCuts(0),
    fLHCCuts(0),
    fDetectorCuts(0),
    fEventCuts(0),
    fTagDirectory(0),
    fValueMin(),
    fValueMax(),
    fValueStep(),
    fValue(),    
    fBinNumber(0),
    fNoMore(0)

{
  // Default constructor
    InitArrays();
}

AliEventPoolOTF::AliEventPoolOTF(const char* name, const char* title):
    AliVEventPool(name, title),
    fTagAnalysis(new AliTagAnalysis(title)),
    fRunCuts(new AliRunTagCuts()),
    fLHCCuts(new AliLHCTagCuts()),
    fDetectorCuts(new AliDetectorTagCuts()),
    fEventCuts(new AliEventTagCuts()),
    fTagDirectory("."),
    fValueMin(),
    fValueMax(),
    fValueStep(),
    fValue(),    
    fBinNumber(0),
    fNoMore(0)

{
  // Constructor
    InitArrays();
}


AliEventPoolOTF::AliEventPoolOTF(const AliEventPoolOTF& obj):
    AliVEventPool(obj),
    fTagAnalysis(0),
    fRunCuts(0),
    fLHCCuts(0),
    fDetectorCuts(0),
    fEventCuts(0),
    fTagDirectory(0),
    fValueMin(),
    fValueMax(),
    fValueStep(),
    fValue(),    
    fBinNumber(0),
    fNoMore(0)
{
    // Copy constructor
    InitArrays();
}

AliEventPoolOTF& AliEventPoolOTF::operator=(const AliEventPoolOTF& other)
{
// Assignment operator
    AliVEventPool::operator=(other);
    return *this;
}


void AliEventPoolOTF::Init()
{
    //
    fTagAnalysis->ChainLocalTags(fTagDirectory);
    for (Int_t i = 0; i < 4; i++) fValue[i] = fValueMin[i];    
}

TChain* AliEventPoolOTF::GetNextChain()
{
    //
    TChain* chain = 0;
    fBinNumber++;
    if (fNoMore) {
	return 0;
    } else {
	printf("Current bin (lower) %13.3f %13.3f %13.3f \n", fValue[kMultiplicity], fValue[kZVertex], fValue[kEventPlane]);
	printf("Current bin (upper) %13.3f %13.3f %13.3f \n", fValue[kMultiplicity] + fValueStep[kMultiplicity], 
	                                                      fValue[kZVertex]      + fValueStep[kZVertex], 
	                                                      fValue[kEventPlane]   + fValueStep[kEventPlane]);
	
	fEventCuts->SetMultiplicityRange(Int_t(fValue[kMultiplicity]) , Int_t(fValue[kMultiplicity] + fValueStep[kMultiplicity]));
	fEventCuts->SetPrimaryVertexZRange(fValue[kZVertex] , fValue[kZVertex] + fValueStep[kZVertex]);
	fEventCuts->SetEventPlaneAngleRange(fValue[kEventPlane] , fValue[kEventPlane] + fValueStep[kEventPlane]);
	chain = fTagAnalysis->QueryTags(fRunCuts, fLHCCuts, fDetectorCuts, fEventCuts);
//
//      Next bin 
//
	for (Int_t i = 2; i >= 0; i--) 
	{
	    fValue[i] += fValueStep[i];
	    if (i > 0  && fValue[i] >= fValueMax[i]) {
		fValue[i] = fValueMin[i];
	    } else if (i == 0 && fValue[i] >= fValueMax[i]) {
		fNoMore = kTRUE;
	    } else {
		break;
	    }
	}
	return chain;
    }
}

void  AliEventPoolOTF::GetCurrentBin(Float_t* /*bin*/)
{
    //
}

Int_t AliEventPoolOTF::GetDimension()
{
    //
    return (3);
}

void AliEventPoolOTF::InitArrays()
{
    // Initializes the pool axis
    
    SetMultiplicityBinning(0, 20000, 20000);
    SetZVertexBinning(-1000., 1000., 2000.);
    SetEventPlaneBinning(-1000., 1000., 2000.);
    SetLeadingParticleEtaBinning(-1.0, 1.0, 2.);    
    for (Int_t i = 0; i < 4; i++) fValue[i] = fValueMin[i];
}


