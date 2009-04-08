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
 
#include <TChain.h>
#include <TTree.h>
#include <TString.h>
#include <TFile.h>

#include "AliAnalysisTaskTagCreator.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODHandler.h"
#include "AliRunTag.h"
#include "AliEventTag.h"
#include "AliAODTagCreator.h"
#include "AliLog.h"

ClassImp(AliAnalysisTaskTagCreator)

////////////////////////////////////////////////////////////////////////

AliAnalysisTaskTagCreator::AliAnalysisTaskTagCreator():
    AliAnalysisTaskSE(),
    fCreateTags(kFALSE),
    fFirstFile(kTRUE),
    fRunTag(0), 
    fTreeT(0),
    fTagCreator(0)
{
  // Default constructor
}

AliAnalysisTaskTagCreator::AliAnalysisTaskTagCreator(const char* name):
    AliAnalysisTaskSE(name),
    fCreateTags(kFALSE),
    fFirstFile(kTRUE),
    fRunTag(0), 
    fTreeT(0),
    fTagCreator(0)
{
  // Constructor
    DefineOutput(1, TTree::Class()); 	
}

void AliAnalysisTaskTagCreator::UserCreateOutputObjects()
{
// Create the output container
    OpenFile(1);
    fTreeT  = new TTree("T", "AOD Tags");
    fRunTag = new AliRunTag();
    TBranch * btag = fTreeT->Branch("AliTAG", "AliRunTag", &fRunTag);
    btag->SetCompressionLevel(9);
    fTagCreator = new AliAODTagCreator();
}

void AliAnalysisTaskTagCreator::Init()
{
    // Initialization
    if (fDebug > 1) AliInfo("Init() \n");
    // Call configuration file
}


void AliAnalysisTaskTagCreator::UserExec(Option_t */*option*/)
{
    // Create Tags for the current event
    AliEventTag* evtTag = new AliEventTag();
    fTagCreator->FillEventTag(AODEvent(), evtTag);
    // Reference to the input file
    TString fturl, fturltemp, fguid;
    
    TString opt(fInputHandler->GetAnalysisType());
    opt.ToLower();
    
    TFile *file = OutputTree()->GetCurrentFile();
    const TUrl *url = file->GetEndpointUrl();
    fguid = file->GetUUID().AsString();
    if (opt.Contains("grid")) {
	fturltemp = "alien://"; fturltemp += url->GetFile();
	fturl = fturltemp(0,fturltemp.Index(".root",5,0,TString::kExact)+5);
    } else {
	fturl = url->GetFile();
    }
    evtTag->SetGUID(fguid);
    if(opt.Contains("grid")) {
	evtTag->SetMD5(0);
	evtTag->SetTURL(fturl);
	evtTag->SetSize(0);
    }
    else evtTag->SetPath(fturl);
    //
    // Add the event tag
    fRunTag->AddEventTag(*evtTag);
    PostData(1, fTreeT);
}


void AliAnalysisTaskTagCreator::FinishTaskOutput()
{
// Terminate analysis
//
    if (fInputHandler->GetRunTag()) fRunTag->CopyStandardContent(fInputHandler->GetRunTag());	    
    fTreeT->Fill();
}

Bool_t AliAnalysisTaskTagCreator::Notify()
{
    // Notify file change
    if (!fFirstFile) {
	if (fInputHandler->GetRunTag()) fRunTag->CopyStandardContent(fInputHandler->GetRunTag());	    
	fTreeT->Fill();
	fRunTag->Clear();
    } else {
	fFirstFile = kFALSE;
    }
    return kTRUE;
}


void AliAnalysisTaskTagCreator::Terminate(Option_t */*option*/)
{
// Terminate analysis
//
    if (fDebug > 1) printf("AnalysisTagCreator: Terminate() \n");
}

