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
    fMultMin(0),
    fMultMax(0),
    fMultStep(0),
    fMultiplicity(),
    fBinNumber(0)
{
  // Default constructor
}

AliEventPoolOTF::AliEventPoolOTF(const char* name, const char* title):
    AliVEventPool(name, title),
    fTagAnalysis(new AliTagAnalysis("AOD")),
    fRunCuts(new AliRunTagCuts()),
    fLHCCuts(new AliLHCTagCuts()),
    fDetectorCuts(new AliDetectorTagCuts()),
    fEventCuts(new AliEventTagCuts()),
    fTagDirectory("."),
    fMultMin(0),
    fMultMax(0),
    fMultStep(0),
    fMultiplicity(),
    fBinNumber(0)
{
  // Constructor
}


AliEventPoolOTF::AliEventPoolOTF(const AliEventPoolOTF& obj):
    AliVEventPool(obj),
    fTagAnalysis(0),
    fRunCuts(0),
    fLHCCuts(0),
    fDetectorCuts(0),
    fEventCuts(0),
    fTagDirectory(0),
    fMultMin(0),
    fMultMax(0),
    fMultStep(0),
    fMultiplicity(),
    fBinNumber(0)
{
    // Copy constructor
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
    fMultiplicity = fMultMin;
}

TChain* AliEventPoolOTF::GetNextChain()
{
    //
    TChain* chain = 0;
    fBinNumber++;
    Int_t mmax = fMultiplicity + fMultStep - 1;
    if (mmax > fMultMax) {
	return 0;
    } else {
	fEventCuts->SetMultiplicityRange(fMultiplicity, mmax);
	chain = fTagAnalysis->QueryTags(fRunCuts, fLHCCuts, fDetectorCuts, fEventCuts);
	fMultiplicity += fMultStep;
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
    return (1);
}

