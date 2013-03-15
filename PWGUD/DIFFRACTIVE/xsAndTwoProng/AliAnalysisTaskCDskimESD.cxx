/*************************************************************************
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
//
// $Id$
//
// Task to skim ESD files, basic idea found at PWGGA/EMCALTasks as
// AliEsdSkimTask. This version is modified to produce a exact copy of the ESD
// stream containing only events which are of interest for the central
// diffractive analysis (PWGUD) and a control sample. In order to be able to
// process MC data as well, the usual concept using an output slot is for
// output data is violated and an output file is directly created without
// a corresponding output slot.
//
// Author:
//  Felix Reidt <Felix.Reidt@cern.ch>

// ROOT headers
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TAxis.h>
#include <TRandom3.h>
//#include <TKey.h>
//#include <TROOT.h>

// AliRoot headers
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDFMD.h"
#include "AliEventplane.h"
#include "AliInputEventHandler.h"
#include "AliMultiplicity.h"
#include "AliVEventHandler.h"
#include "AliMCEventHandler.h"
//#include "AliMCEvent.h"
//#include "AliStack.h"
//#include "AliGenEventHeader.h"
//#include "AliRunLoader.h"

// own classes
#include "AliCDMesonBase.h"
#include "AliCDMesonUtils.h"

// own header
#include "AliAnalysisTaskCDskimESD.h"

ClassImp(AliAnalysisTaskCDskimESD)

//______________________________________________________________________________
AliAnalysisTaskCDskimESD::AliAnalysisTaskCDskimESD(const char *opt,
                                                   Bool_t reduceGapEvents)
	: AliAnalysisTaskSE(opt)
	, fInputTree(0x0)
	, fInputTreeCopy(0x0)
	, fOutputTree(0x0)
	, fDoMC(kFALSE)
	, fWrkDir(0x0)
	, fKinematicsFile(0x0)
	, fTrackRefsFile(0x0)
	, fTEinput(0x0)
	, fTE(0x0)
	, fgaliceFile(0x0)
	, fDoGapCond(kTRUE)
	, fCurrentGapCondition(0)
	, fRequestedGapCondition(AliCDMesonBase::kBitV0A + AliCDMesonBase::kBitV0C)
	                        // only events containing a V0 double gap or
	, fNoGapFraction(0.006) // 0.6% of the minimum bias events are stored
	, fReduceGapEvents(reduceGapEvents)
	, fRefinedGapCondition(AliCDMesonBase::kBitV0A + AliCDMesonBase::kBitV0C +
	                       AliCDMesonBase::kBitFMDA + AliCDMesonBase::kBitFMDC)
	, fReducedGapFraction(0.02) // 2% of the events with a loose gap are stored
	, fStatsFlow(0x0)
	, fSkimmingList(0x0)
	, fFileName()
	, fEventNumberInFile(-1)
	, fRunNumber(-1)
	, fEventTime(0)
{
	//
	// Constructor.
	//

	if (!opt)
		return;

	//fBranchNames = "ESD:AliESDHeader.,AliESDRun."; // TODO don't we need all?

	// check whether we are running on MC data
	AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
	AliMCEventHandler* mcH =
		dynamic_cast<AliMCEventHandler*>(am->GetMCtruthEventHandler());
	if(mcH) {
		fDoMC = kTRUE;

		fWrkDir = gDirectory; // save the old directory
		fKinematicsFile = new TFile("Kinematics.root", "recreate");
		fKinematicsFile->Close();
		delete fKinematicsFile;
		fKinematicsFile = 0x0;
		gDirectory = fWrkDir; // restore the old directory

		fTrackRefsFile = new TFile("TrackRefs.root", "recreate");
		fTrackRefsFile->Close();
		delete fTrackRefsFile;
		fTrackRefsFile = 0x0;
		gDirectory = fWrkDir; // restore the old directory
	}

	DefineOutput(1, TTree::Class());
	DefineOutput(2, TH1::Class());
	DefineOutput(3, TTree::Class());
	if (fDoMC) DefineOutput(4, TTree::Class());
}


//______________________________________________________________________________
AliAnalysisTaskCDskimESD::AliAnalysisTaskCDskimESD()
	: AliAnalysisTaskSE()
	, fInputTree(0x0)
	, fInputTreeCopy(0x0)
	, fOutputTree(0x0)
	, fDoMC(kFALSE)
	, fWrkDir(0x0)
	, fKinematicsFile(0x0)
	, fTrackRefsFile(0x0)
	, fTEinput(0x0)
	, fTE(0x0)
	, fgaliceFile(0x0)
	, fDoGapCond(kTRUE)
	, fCurrentGapCondition(0)
	, fRequestedGapCondition(AliCDMesonBase::kBitV0A + AliCDMesonBase::kBitV0C)
	                        // only events containing a V0 double gap or
	, fNoGapFraction(0.006) // 0.6% of the minimum bias events are stored
	, fReduceGapEvents(kFALSE)
	, fRefinedGapCondition(AliCDMesonBase::kBitV0A + AliCDMesonBase::kBitV0C +
	                       AliCDMesonBase::kBitFMDA + AliCDMesonBase::kBitFMDC)
	, fReducedGapFraction(0.02) // 2% of the events with a loose gap are stored
	, fStatsFlow(0x0)
	, fSkimmingList(0x0)
	, fFileName()
	, fEventNumberInFile(-1)
	, fRunNumber(-1)
	, fEventTime(0)
{
	//
	// Default Constructor
	//
}


//______________________________________________________________________________
AliAnalysisTaskCDskimESD::~AliAnalysisTaskCDskimESD()
{
	//
	// deconstructor
	//

	/* // this delete lead to chrases caused by the garbage collection ...
	if (fOutputTree) {
		delete fOutputTree;
		fOutputTree = 0x0;
		}
	if (fInputTreeCopy) {
		delete fInputTreeCopy;
		fInputTreeCopy = 0x0;
	}
	if (fStatsFlow) {
		delete fStatsFlow;
		fStatsFlow = 0x0;
	} */
	if (fDoMC && fKinematicsFile) {
		fKinematicsFile->Close();
		delete fKinematicsFile;
		fKinematicsFile = 0x0;
	}
	if (fDoMC && fTrackRefsFile) {
		fTrackRefsFile->Close();
		delete fTrackRefsFile;
		fTrackRefsFile = 0x0;
	}
}


//______________________________________________________________________________
void AliAnalysisTaskCDskimESD::UserExec(Option_t */*opt*/)
{
	//
  // Process event.
	//

	fStatsFlow->Fill(0.); // UserExec(...) <= bin name

	if (fDoMC) {
		// replace dummy tree by the correct tree
		AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
		if(!fTEinput) {
			fTEinput =
				((AliMCEventHandler*)am->GetMCtruthEventHandler())->GetTree();
			if (fTEinput) {
				delete fTE;
				fTE = fTEinput->CloneTree(0);
				fTE->SetDirectory(fgaliceFile);
				fTE->SetAutoFlush(-10*1024*1024);
				//CopygALICE(am);
			}
		}
		else if (fTEinput != am->GetMCtruthEventHandler()->GetTree()) {
			fTEinput =
				((AliMCEventHandler*)am->GetMCtruthEventHandler())->GetTree();
			fTEinput->CopyAddresses(fTE);
		}
	}

	// check whether the ESD input tree is still the same - NOT NECESSARY
	//TTree* currentInput =
	//	AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetTree();
	//if (fInputTree != currentInput) { // ESD input tree has changed
	//	puts("CHANGED INPUT TREE");
	//	fInputTree = currentInput;
	//	fInputTreeCopy = (TTree*)fInputTree->Clone("esdTree");
	//	fInputTreeCopy->CopyAddresses(fOutputTree);
	//}

	AliESDEvent *esdin = dynamic_cast<AliESDEvent*>(InputEvent());
	if (!esdin || !fInputTree || !fInputTreeCopy || !fOutputTree || !fTEinput
	    || !fTE) {
		// check whether everything is ready
		PostData(1, fOutputTree);
		PostData(2, fStatsFlow);
		PostData(3, fSkimmingList);
		if (fDoMC) PostData(4, fTE);
		return;
	}

	fStatsFlow->Fill(1.); // Ready <= bin name

	//============================================================================
	// event selection, every event passing the next lines is written to the new
	// ESD stream
	//
	// modify this lines in order to adjust the ESD selection to your needs

	// determine current gap condition
	fCurrentGapCondition = AliCDMesonUtils::GetGapConfig(esdin, 0x0, 0x0, 0x0,
	                                                     0x0, 0x0, 0x0, 0x0, 0x0);

	TRandom3 rnd(0);
	if (!(AliCDMesonBase::kBitBaseLine & fCurrentGapCondition) ||
	     (fRequestedGapCondition & fCurrentGapCondition) || fReduceGapEvents) {
		// check whether the event satifies the required gap condition and whether
		// all double-gap events will be stored

		if (fReduceGapEvents &&
		    (AliCDMesonBase::kBitBaseLine & fCurrentGapCondition) &&
		    !(fRequestedGapCondition & fCurrentGapCondition)) {
			// not all double gap events can be stored, but the gap determination
			// was successful and at least a according to the less stringent condition
			// a double gap was found
			if (!(fRefinedGapCondition & fCurrentGapCondition)) {
				// event fulfilled the more stringent gap and hence is will be stored
				fStatsFlow->Fill(5.);
			}
			else if (rnd.Rndm() > (1.-fReducedGapFraction)) {
				// randomly select a fraction of
				fStatsFlow->Fill(4.); // Double Gap Event <= bin name
			}
			else {
				// reject event (no refined gap and not randomly selected)
				fStatsFlow->Fill(3.); // Event Rejected <= bin name
				PostData(1, fOutputTree);
				PostData(2, fStatsFlow);
				PostData(3, fSkimmingList);
				if (fDoMC) PostData(4, fTE);
				return;
			}
		}
		else if (rnd.Rndm() > (1.-fNoGapFraction)) { // randomly selected
			fStatsFlow->Fill(2.); // Control Sample Event <= bin name
		}
		else {
			fStatsFlow->Fill(3.); // Event Rejected <= bin name
			PostData(1, fOutputTree);
			PostData(2, fStatsFlow);
			PostData(3, fSkimmingList);
			if (fDoMC) PostData(4, fTE);
			return;
		}
	}
	else {
		fStatsFlow->Fill(4.); // Double Gap Event <= bin name
	}
	// end of event selection
	//============================================================================

	// load the current event in the cloned input tree
	//
	// unfortunately the Entry() gives the entry number with in the current input
	// tree, not within the chain for MC data
	// => use GetReadEntry() of the ESD tree (value return by Entry() for real
	//data)
	fInputTreeCopy->GetEvent(fInputTree->GetReadEntry());
	fOutputTree->Fill(); // fill the current event into the new ESD stream
	//printf("Entry()=%lld\nEntries=%lld\nChainOffset=%lld\nReadEntry=%lld\n",
	//       Entry(), fInputTreeCopy->GetEntries(),
	//       fInputTreeCopy->GetChainOffset(), fInputTree->GetReadEntry());

	// MC specific stuff ---------------------------------------------------------
	AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();

	if (fDoMC) { // do the MC related part of the skimming
		Long64_t iEntry = fOutputTree->GetEntries() - 1;
		AliMCEventHandler* mcHandler =
			(AliMCEventHandler*)am->GetMCtruthEventHandler();
		TTree* treeK = mcHandler->TreeK();
		TTree* treeTR = mcHandler->TreeTR();
		if (!treeK) AliFatal("TreeK not found!");
		if (!treeTR) AliFatal("TreeTR not found!");

		fWrkDir = gDirectory;
		gDirectory = fKinematicsFile;
		TString dirName = TString(Form("Event%ld", (Long_t)iEntry));
		gDirectory->mkdir(dirName);
		gDirectory->cd(dirName);
		TTree* outputTreeK = treeK->CloneTree(); // copy the whole tree
		outputTreeK->Write();
		treeK->CopyAddresses(outputTreeK, kTRUE); // separate clone again
		treeK->GetListOfClones()->Remove((TObject*)outputTreeK);
		outputTreeK->Delete();
		outputTreeK = 0x0;

		gDirectory = fTrackRefsFile;
		gDirectory->mkdir(dirName);
		gDirectory->cd(dirName);
		TTree* outputTreeTR = treeTR->CloneTree(); // copy the whole tree
		outputTreeTR->Write();
		treeTR->CopyAddresses(outputTreeTR, kTRUE); // separate clone again
		treeTR->GetListOfClones()->Remove((TObject*)outputTreeTR);
		outputTreeTR->Delete();
		outputTreeTR = 0x0;
		gDirectory = fWrkDir;

		TTree* currentTEinput =
			((AliMCEventHandler*)am->GetMCtruthEventHandler())->GetTree();
		if (fTEinput != currentTEinput) {
			// when the input file is changed, the kinematics tree address changes
			// hence our output tree fTE has to be linked to the new one!
			fTEinput = currentTEinput;
			fTEinput->CopyAddresses(fTE);
		}
		fTE->Fill(); // fill MC event headers
	}

	fStatsFlow->Fill(6.); // Tree Filled <= bin name

	// remember which event we stored
	fEventNumberInFile = esdin->GetEventNumberInFile();
	fRunNumber = esdin->GetRunNumber();
	fFileName = fInputTree->GetCurrentFile()->GetName();
	fEventTime = esdin->GetTimeStamp();
	fSkimmingList->Fill();

	PostData(1, fOutputTree);
	PostData(2, fStatsFlow);
	PostData(3, fSkimmingList);
	if (fDoMC) PostData(4, fTE);
}

//______________________________________________________________________________
void AliAnalysisTaskCDskimESD::UserCreateOutputObjects()
{
	// Create output objects.
	fStatsFlow = (TH1*)(new TH1F("skimmingStatsFlow", "", 7, 0., 7.));
	TAxis* x = fStatsFlow->GetXaxis();
	x->SetBinLabel(1, "UserExec(...)");
	x->SetBinLabel(2, "Ready");
	x->SetBinLabel(3, "Control Sample Event");
	x->SetBinLabel(4, "Rejected Event");
	x->SetBinLabel(5, "Double Gap Event");
	x->SetBinLabel(6, "Refined Double Gap Event");
	x->SetBinLabel(7, "Tree Filled");

	PostData(2, fStatsFlow);

	// TODO implement fSkimmingList
	fSkimmingList = new TTree("SkimmingList", "SkimmingList");
	fSkimmingList->Branch("FileName", &fFileName);
	fSkimmingList->Branch("EventNumberInFile", &fEventNumberInFile);
	fSkimmingList->Branch("RunNumber", &fRunNumber);
	fSkimmingList->Branch("TimeStamp", &fEventTime);
	PostData(3, fSkimmingList);

	// Get input information
	AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
	fInputTree = am->GetInputEventHandler()->GetTree();
	if (!fInputTree) {
		puts("AliAnalysisTaskCDskimESD: input tree not found\n");
		return;
	}
	fInputTreeCopy = (TTree*)fInputTree->Clone("esdTree");


	// prevent the task from being run on proof
	if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() ==
	    AliAnalysisManager::kProofAnalysis) {
		AliFatal("AliAnalysisTaskCDskimESD: cannot be run on PROOF!");
	}

	TFile *file = 0x0;
	file = OpenFile(1); // open file for the first output slot
	if (!file) {
		puts("AliAnalysisTaskCDskimESD: file not opened\n");
			return;
	}

	fOutputTree = fInputTreeCopy->CloneTree(0);
	if (!fOutputTree) {
		puts("AliAnalysisTaskCDskimESD: cloning tree not successful\n");
		return;
	}

	file->SetCompressionLevel(5); // caused a crash
	fOutputTree->SetDirectory(file);
	fOutputTree->SetAutoFlush(-10*1024*1024);

	if (fDoGapCond) {
		if(!fOutputTree->GetBranch("gapCondition")) {
			fOutputTree->Branch("gapCondition", &fCurrentGapCondition);
		}
	}

	PostData(1, fOutputTree);

	// process the mc information
	if(fDoMC) {
		fWrkDir = gDirectory; // save the old directory
		fKinematicsFile = new TFile("Kinematics.root", "update");
		gDirectory = fWrkDir; // restore the old directory

		fTrackRefsFile = new TFile("TrackRefs.root", "update");
		gDirectory = fWrkDir; // restore the old directory

		fTEinput =
			((AliMCEventHandler*)am->GetMCtruthEventHandler())->GetTree();
		if (!fTEinput) {
			// this trick is done in order to post an output tree, although the tree
			// which is cloned is not available during the UserCreateOutputObjects()
			fTE = new TTree("TE", "TE");
		}
		else {
			fTE = fTEinput->CloneTree(0);
			//CopygALICE(am); // not needed, not properly working/tested
		}
		fgaliceFile = OpenFile(4);
		if (!file) {
			puts("AliAnalysisTaskCDskimESD: galice.root file not opened\n");
			return;
		}
		fgaliceFile->SetCompressionLevel(5);
		fTE->SetDirectory(fgaliceFile);
		fTE->SetAutoFlush(-10*1024*1024);


		PostData(4, fTE);
	}
}


//______________________________________________________________________________
/* // not needed => not yet working / properly tested
void AliAnalysisTaskCDskimESD::CopygALICE(AliAnalysisManager* am)
{
	//
	// copy all objects contained in the galice root excluding the TE tree
	//

	TString galiceInput =
		*(((AliMCEventHandler*)am->GetMCtruthEventHandler())->GetInputPath());
	galiceInput += "galice.root";
	printf("galiceInput=%s\n", galiceInput.Data());
	fWrkDir = gDirectory;
	TFile* input = new TFile(galiceInput.Data(), "READ");

	input->cd();
	gDirectory->ls();

	//loop on all entries of of the input directory
	TKey *key;
	TIter nextkey(input->GetListOfKeys());
	input->GetListOfKeys()->ls();
	while ((key = (TKey*)nextkey())) {
		TString name = key->GetName();
		printf("name=%s\n", name.Data());
		const TString classname = key->GetClassName();
		printf("classname=%s\n", classname.Data());
		TClass *cl = gROOT->GetClass(classname);
		if (!cl) continue;
		if (!name.Contains("TE") && !classname.Contains("AliRun")) {
			TObject *obj = key->ReadObj();
			fgaliceFile->cd();
			obj->Write();
			delete obj;
			input->cd();
		}
		else if (cl->InheritsFrom("AliRunLoader")) {
			// TODO write the runloader stuff
			//AliRunLoader* rl = AliRunLoader::Instance();
			//printf("FILENAME=%s\n", rl->GetFileName().Data());
		}
		fgaliceFile->ls();
	}
	fgaliceFile->SaveSelf(kTRUE);
	fWrkDir->cd();
}
*/
