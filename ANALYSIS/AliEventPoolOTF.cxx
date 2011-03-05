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
#include <TChain.h>
#include <TGridResult.h>

ClassImp(AliEventPoolOTF)


////////////////////////////////////////////////////////////////////////

AliEventPoolOTF::AliEventPoolOTF():
    AliVEventPool(),
    fTagAnalysis(0),
    fRunCuts(0),
    fLHCCuts(0),
    fDetectorCuts(0),
    fEventCuts(0),
    fGridTags(0),
    fChain(0),
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
    fGridTags(0),
    fChain(0),
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
    fGridTags(0),
    fChain(0),
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


AliEventPoolOTF::~AliEventPoolOTF()
{
    // Destructor
    delete fTagAnalysis;
    delete fRunCuts;
    delete fEventCuts;
    delete fLHCCuts;
    delete fDetectorCuts;
    delete fChain;
}

AliEventPoolOTF& AliEventPoolOTF::operator=(const AliEventPoolOTF& other)
{
// Assignment operator
    AliVEventPool::operator=(other);
    return *this;
}


void AliEventPoolOTF::Init()
{
    // Initialisation
    if (!fGridTags) {
	fTagAnalysis->ChainLocalTags(fTagDirectory);
    } else {
	fTagAnalysis->ChainGridTags(fGridTags);
    }
    
    
    for (Int_t i = 0; i < 5; i++) fValue[i] = fValueMin[i];    
}

TChain* AliEventPoolOTF::GetNextChain()
{
    // Get Next Chain
    if (fChain) {
	delete fChain;
	fChain = 0;
    }

    fBinNumber++;
    if (fNoMore) {
 	return 0;
    } else {
    printf("Current bin (lower) %13.3f %13.3f %13.3f %13.3f %13.3f \n", fValue[kMultiplicity], fValue[kZVertex], fValue[kEventPlane],fValue[kLeadingParticleEta],fValue[kLeadingParticlePhi]);
    printf("Current bin (upper) %13.3f %13.3f %13.3f %13.3f %13.3f \n", fValue[kMultiplicity] + fValueStep[kMultiplicity] - 1, 
	   fValue[kZVertex]      + fValueStep[kZVertex], 
	   fValue[kEventPlane]   + fValueStep[kEventPlane],
	   fValue[kLeadingParticleEta]   + fValueStep[kLeadingParticleEta],
           fValue[kLeadingParticlePhi]   + fValueStep[kLeadingParticlePhi]
    
	   );

	fEventCuts->SetMultiplicityRange(Int_t(fValue[kMultiplicity]) , Int_t(fValue[kMultiplicity] + fValueStep[kMultiplicity] - 1));
	fEventCuts->SetPrimaryVertexZRange(fValue[kZVertex] , fValue[kZVertex] + fValueStep[kZVertex]);
        fEventCuts->SetEtaLeadingParticleRange(fValue[kLeadingParticleEta] , fValue[kLeadingParticleEta] + fValueStep[kLeadingParticleEta]);
        fEventCuts->SetPhiLeadingParticleRange(fValue[kLeadingParticlePhi] , fValue[kLeadingParticlePhi] + fValueStep[kLeadingParticlePhi]);
        fEventCuts->SetEventPlaneAngleRange(fValue[kEventPlane] , fValue[kEventPlane] + fValueStep[kEventPlane]);
    
	fChain = fTagAnalysis->QueryTags(fRunCuts, fLHCCuts, fDetectorCuts, fEventCuts);
//
//      Next bin 
//
	for (Int_t i = 5; i >= 0; i--) 
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
	return fChain;
    }
}

void  AliEventPoolOTF::GetCurrentBin(Float_t* /*bin*/)
{
    //
}

Int_t AliEventPoolOTF::GetDimension()
{
    //
    return (5);
}

void AliEventPoolOTF::InitArrays()
{
    // Initializes the pool axis
    
    SetMultiplicityBinning(0, 20000, 20000);
    SetZVertexBinning(-1000., 1000., 2000.);
    SetEventPlaneBinning(-1000., 1000., 2000.);
    SetLeadingParticleEtaBinning(-13.0, 13.0, 27.);
    SetLeadingParticlePhiBinning(0., 2*(TMath::Pi()),2*(TMath::Pi()));
    for (Int_t i = 0; i < 5; i++) fValue[i] = fValueMin[i];
}


