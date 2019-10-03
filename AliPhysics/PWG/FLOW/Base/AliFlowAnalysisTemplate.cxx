/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

#define AliFlowAnalysisTemplate_CXX

// root includes
#include "TFile.h"      
#include "TList.h"
#include "TMath.h"
#include "TProfile.h"
// aliroot includes
#include "AliFlowCommonConstants.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliFlowAnalysisTemplate.h"

// Description: Template maker to serve as a starting point for flow analysis
// Author:      Redmer Alexander Bertens, Utrecht University, 2013
//              rbertens@cern.ch, rbertens@nikhef.nl, r.a.bertens@uu.nl 


ClassImp(AliFlowAnalysisTemplate)

AliFlowAnalysisTemplate::AliFlowAnalysisTemplate() :
    fDebug                      (0),
    fUsePhiWeights              (kFALSE),
    fApplyCorrectionForNUA      (kFALSE),
    fHarmonic                   (2),
    fWeightsList                (0x0),
    fHistList                   (0x0),
    fCommonHists                (0x0),
    fCommonHistsRes             (0x0)
{ /* constructor */ }
//_____________________________________________________________________________
AliFlowAnalysisTemplate::~AliFlowAnalysisTemplate() 
{
  // destructor
}
//_____________________________________________________________________________
void AliFlowAnalysisTemplate::Init() 
{
    // initialize the histograms
    fHistList = new TList();
    fHistList->SetName("TemplateList");
    fHistList->SetOwner();

    // common histogram container
    fCommonHists = new AliFlowCommonHist("AliFlowCommonHist_Template","AliFlowCommonHist");
    (fCommonHists->GetHarmonic())->Fill(0.5,fHarmonic); // store harmonic 
    fHistList->Add(fCommonHists);
    
    // common results container
    fCommonHistsRes = new AliFlowCommonHistResults("AliFlowCommonHistResults_Template","",fHarmonic);
    fHistList->Add(fCommonHistsRes);
}
//_____________________________________________________________________________
void AliFlowAnalysisTemplate::Make(AliFlowEventSimple* anEvent) 
{
  // core method, called for each event
  if (!anEvent) return;
  // test statement
  printf("Numer of POIs %i", anEvent->NumberOfTracks());
}
//_____________________________________________________________________________
void AliFlowAnalysisTemplate::GetOutputHistograms(TList *outputListHistos)
{
    //get pointers to all output histograms (called before Finish())
    fHistList = outputListHistos;
    fCommonHists = (AliFlowCommonHist*) fHistList->FindObject("AliFlowCommonHist_SP");
    fCommonHistsRes = (AliFlowCommonHistResults*) fHistList->FindObject("AliFlowCommonHistResults_SP");
}   
//_____________________________________________________________________________
void AliFlowAnalysisTemplate::Finish() {
    // calculate flow and fill the AliFlowCommonHistResults
}
//_____________________________________________________________________________
void AliFlowAnalysisTemplate::WriteHistograms(TDirectoryFile *outputFileName) const {
    // store the final results in output .root file
    outputFileName->Add(fHistList);
    outputFileName->Write(outputFileName->GetName(), TObject::kSingleKey);
}
//_____________________________________________________________________________
