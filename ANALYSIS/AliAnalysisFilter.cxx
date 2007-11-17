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

//
// Manager class for filter decisions based on cuts
// The filter contains a list of sets of cuts.
// A bit field is filled in order to store the decision of each cut-set. 
// Author: Andreas Morsch
// andreas.morsch@cern.ch

#include <TObject.h>
#include <TList.h>
#include "AliAnalysisFilter.h"
#include "AliAnalysisCuts.h"


ClassImp(AliAnalysisFilter)


////////////////////////////////////////////////////////////////////////

AliAnalysisFilter::AliAnalysisFilter():
    TNamed(),
    fCuts(0)
{
  // Default constructor
}

AliAnalysisFilter::AliAnalysisFilter(const char* name, const char* title):
    TNamed(name, title),
    fCuts(new TList())
{
  // Constructor
}

AliAnalysisFilter::AliAnalysisFilter(const AliAnalysisFilter& obj):
    TNamed(obj)
{
// Copy constructor
}


UInt_t AliAnalysisFilter::IsSelected(TObject* obj)
{
    //
    // Loop over all set of cuts
    // and store the decision
    UInt_t result = 0;
    TIter next(fCuts);
    AliAnalysisCuts *cuts;
    Int_t iCutB = 1;
	
    while((cuts = (AliAnalysisCuts*)next())) {
	Bool_t acc = cuts->IsSelected(obj);
	if (acc) {result |= iCutB & 0x00ffffff;}
	iCutB *= 2;
    }  

    return result;
}

void AliAnalysisFilter::Init()
{
    //
    // Loop over all set of cuts and call Init
    TIter next(fCuts);
    AliAnalysisCuts *cuts;
    while((cuts = (AliAnalysisCuts*)next())) cuts->Init();
}

void AliAnalysisFilter::AddCuts(AliAnalysisCuts* cuts)
{
    // Add a set of cuts
    fCuts->Add(cuts);
}
