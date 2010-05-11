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

#include <Riostream.h>

#include <TChain.h>
#include <TTree.h>
#include <TString.h>
#include <TFile.h>
#include <TSystem.h>

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
    fTagCreator(0),
    fAODFileName("")
{
  // Default constructor
}

AliAnalysisTaskTagCreator::AliAnalysisTaskTagCreator(const char* name):
    AliAnalysisTaskSE(name),
    fCreateTags(kFALSE),
    fFirstFile(kTRUE),
    fRunTag(0), 
    fTreeT(0),
    fTagCreator(0),
    fAODFileName("")
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

}

void AliAnalysisTaskTagCreator::ConnectInputData(Option_t * /*option*/)
{
    // Initialization
    const char* turl = gSystem->Getenv("ALIEN_JDL_OUTPUTDIR");
    TString sturl = turl;
    
    if (sturl.Length() != 0) {
      fAODFileName = "alien://";
      fAODFileName += turl;
      fAODFileName += "/AliAOD.root";
    }  
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
    if (fAODFileName.Length() != 0) {
	fturl = fAODFileName;
	GetGUID(fguid);
    } else {
	fturl = url->GetFile();
    }

    evtTag->SetGUID(fguid);
    if(fAODFileName.Length() != 0) {
	evtTag->SetMD5("");
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
    fInputHandler = (AliInputEventHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());

    if (!fFirstFile) {
	if (fInputHandler->GetRunTag()) fRunTag->CopyStandardContent(fInputHandler->GetRunTag());	    
	fTreeT->Fill();
	fRunTag->Clear();
    } else {
	fFirstFile = kFALSE;
    }
    return kTRUE;
}


void AliAnalysisTaskTagCreator::GetGUID(TString &guid) {
    // Get the guid of the AliAOD.root file
    ofstream myfile ("guid.txt");
    if (myfile.is_open()) {
	TFile *f = TFile::Open("AliAOD.root","read");
	if(f->IsOpen()) {
	    guid = f->GetUUID().AsString();
	    myfile << "AliAOD.root \t"<<f->GetUUID().AsString();
	    cout<<guid.Data()<<endl;
	    myfile.close();
	}
	else cout<<"Input file not found"<<endl;
    }
    else cout<<"Output file can't be created..."<<endl;
}



void AliAnalysisTaskTagCreator::Terminate(Option_t */*option*/)
{
// Terminate analysis
//
    if (fDebug > 1) printf("AnalysisTagCreator: Terminate() \n");
}


