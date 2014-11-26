/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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
/*
 * A small task dumping all EMCal trigger related information into a TTree
 *      Author: Markus Fasel
 */
#include <iostream>
#include <string>
#include <TMath.h>
#include <TTree.h>

#include "AliInputEventHandler.h"
#include "AliVTrack.h"
#include "AliVCluster.h"

#include "AliAnalysisTaskEmcalTriggerTreeWriter.h"

AliAnalysisTaskEmcalTriggerTreeWriter::AliAnalysisTaskEmcalTriggerTreeWriter():
	AliAnalysisTaskSE(),
	fOutputTree(NULL),
	fOutputInfo()
{
	/*
	 * Dummy constructor
	 */
}

AliAnalysisTaskEmcalTriggerTreeWriter::AliAnalysisTaskEmcalTriggerTreeWriter(const char *name):
	AliAnalysisTaskSE(name),
	fOutputTree(NULL),
	fOutputInfo()
{
	/*
	 * Constructor
	 */
	DefineOutput(1,TTree::Class());
}

AliAnalysisTaskEmcalTriggerTreeWriter::~AliAnalysisTaskEmcalTriggerTreeWriter() {
	/*
	 * Destructor
	 */
	if(fOutputTree) delete fOutputTree;
}

void AliAnalysisTaskEmcalTriggerTreeWriter::UserCreateOutputObjects() {
	/*
	 * Create output tree, with two branches, one for the tracks matched and one for the clusters
	 */

	// Build the tree
	OpenFile(1);
	fOutputTree = new TTree("EMCalTree", "A tree with emcal information");
	fOutputTree->Branch("run", &fOutputInfo.fRun, "pdg/I");
	fOutputTree->Branch("col", &fOutputInfo.fCol, "col/I");
	fOutputTree->Branch("row", &fOutputInfo.fRow, "isUnique/i");
	fOutputTree->Branch("NL0Times", &fOutputInfo.fNL0Times, "NL0Times/I");
	fOutputTree->Branch("Level0Times", fOutputInfo.fLevel0Times, "Level0Times[10]/I");
	fOutputTree->Branch("ADC", &fOutputInfo.fADC, "ADC/I");
	fOutputTree->Branch("Amplitude", &fOutputInfo.fAmplitude, "Amplitude/F");
	fOutputTree->Branch("Time", &fOutputInfo.fTime, "Time/F");
	fOutputTree->Branch("TriggerBits", &fOutputInfo.fTriggerBits, "TriggerBits/I");
	fOutputTree->Branch("L1Threshold", &fOutputInfo.fL1Threshold, "L1Threshold/I");
	fOutputTree->Branch("L1V0", &fOutputInfo.fL1V0, "L1V0/I");
	PostData(1, fOutputTree);
}

void AliAnalysisTaskEmcalTriggerTreeWriter::UserExec(Option_t *) {
	/*
	 * Build the tree
	 */

	AliVCaloTrigger *emctrigger = fInputEvent->GetCaloTrigger(strcmp(fInputHandler->GetDataType(), "ESD" ) == 0 ? "EMCALTrigger" : "EMCALTrigger");
	emctrigger->Reset();
	while(emctrigger->Next()){
		fOutputInfo.Reset();
		fOutputInfo.fRun = fInputEvent->GetRunNumber();
		emctrigger->GetPosition(fOutputInfo.fCol, fOutputInfo.fRow);
		emctrigger->GetNL0Times(fOutputInfo.fNL0Times);
		if(fOutputInfo.fNL0Times > 0 && fOutputInfo.fNL0Times < 10)
			emctrigger->GetL0Times(fOutputInfo.fLevel0Times);
		emctrigger->GetL1TimeSum(fOutputInfo.fADC);
		emctrigger->GetAmplitude(fOutputInfo.fAmplitude);
		emctrigger->GetL1V0(fOutputInfo.fL1V0);
		emctrigger->GetTriggerBits(fOutputInfo.fTriggerBits);
		emctrigger->GetTime(fOutputInfo.fTime);
		emctrigger->GetL1Threshold(fOutputInfo.fL1Threshold);
		fOutputTree->Fill();
	}

	PostData(1, fOutputTree);
}

