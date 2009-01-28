/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Artur Szostak <artursz@iafrica.com>                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// \ingroup macros
/// \file RunLinkToMCAnalysisExample.C
/// \brief Example macro for running the AliAnalysisTaskLinkToMC analysis task.
/// \author Artur Szostak <artursz@iafrica.com>
///
/// This macro shows an example of how to run the AliAnalysisTaskLinkToMC analysis
/// task in local mode. It can be used as a quick check of the analysis task.
/// It will generate AliAOD.root and hists.root files. The hists.root can then
/// be used in the PlotEfficiency.C macro.
/// Run this macro as follows:
///
///  $ aliroot RunLinkToMCAnalysisExample.C\(\"AliESDs.root\"\)
///
/// where AliESDs.root should be the correct path to a root file containing ESD
/// objects created by the offline reconstruction.


#if !defined(__CINT__) || defined(__MAKECINT__)
#error This macro must be run in interpreted mode.
#endif


void RunLinkToMCAnalysisExample(const char* esdFile = "./AliESDs.root")
{
	// Load needed libraries
	gSystem->Load("libTree.so");
	gSystem->Load("libGeom.so");
	gSystem->Load("libVMC.so");
	gSystem->Load("libPhysics.so");
	gSystem->Load("libSTEERBase.so");
	gSystem->Load("libESD.so");
	gSystem->Load("libAOD.so");
	gSystem->Load("libANALYSIS.so");
	gSystem->Load("libANALYSISalice.so");
	gSystem->Load("libPWG3base.so");
	gSystem->Load("libPWG3muon.so");
	
	// Create the TChain for esdTrees in the AliESDs.root file.
	TChain* chain = new TChain("esdTree");
	chain->Add(esdFile);
	if (!chain) return;
	
	// Create the analysis manager and event handlers.
	AliAnalysisManager* mgr = new AliAnalysisManager("Analysis Train", "An example analysis train setup for AliAnalysisTaskLinkToMC.");
	AliESDInputHandler* esdHandler = new AliESDInputHandler();
	mgr->SetInputEventHandler(esdHandler);
	AliMCEventHandler* mcHandler = new AliMCEventHandler();
	mgr->SetMCtruthEventHandler(mcHandler);
	mcHandler->SetReadTR(kTRUE); 
	AliAODHandler* aodHandler = new AliAODHandler();
	mgr->SetOutputEventHandler(aodHandler);
	aodHandler->SetOutputFileName("AliAOD.root");
	
	// Create the analysis task and setup the parameters.
	AliAnalysisTaskLinkToMC* linktask = new AliAnalysisTaskLinkToMC("Task to link ESD tracks to corresponding MC tracks.");
	linktask->MinClusters(6);
	linktask->HardCutLimitX(4);
	linktask->HardCutLimitY(4);
	linktask->SigmaCut(5.);
	linktask->MinClustersInSt45(3);
	linktask->StationMustMatch(1, true);  // At least one cluster in station 1 must match.
	linktask->StationMustMatch(2, true);  // At least one cluster in station 2 must match.
	linktask->StationMustMatch(3, true);  // At least one cluster in station 3 must match.
	linktask->GenerateHistograms(true);
	mgr->AddTask(linktask);
	
	// Create the input and output containers and connect them up to the analysis task.
	AliAnalysisDataContainer* cinEsd = mgr->CreateContainer("cESD", TChain::Class(), AliAnalysisManager::kInputContainer);
	AliAnalysisDataContainer* coutAod = mgr->CreateContainer("cAOD", TTree::Class(), AliAnalysisManager::kOutputContainer, "default");
	AliAnalysisDataContainer* coutHists = mgr->CreateContainer("cHists", TList::Class(), AliAnalysisManager::kOutputContainer, "hists.root");
	mgr->ConnectInput(linktask, 0, cinEsd);
	mgr->ConnectOutput(linktask, 0, coutAod);
	mgr->ConnectOutput(linktask, 1, coutHists);
	
	if (mgr->InitAnalysis())
	{
		mgr->PrintStatus();
		mgr->StartAnalysis("local", chain);
	}
}

