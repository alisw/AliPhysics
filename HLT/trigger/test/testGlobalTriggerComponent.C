/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Artur Szostak <artursz@iafrica.com>                   *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/**
 * @file   testGlobalTriggerComponent.C
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   19 Dec 2008
 *
 * This macro is used to test the AliHLTGlobalTriggerComponent class.
 * A number of tests are run with the AliHLTSystem framework to check that
 * the automatically generated global trigger logic is generated correctly.
 */

#if defined(__CINT__) && (! defined(__MAKECINT__))
#error This macro must be compiled. Try running as testGlobalTriggerComponent.C++
#endif

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TSystem.h"
#include "TClassTable.h"
#include "TFile.h"
#include "AliLog.h"
#include "AliHLTReadoutList.h"
#include "AliHLTTriggerDomain.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTGlobalTriggerDecision.h"
#include "AliHLTSystem.h"
#include "AliHLTConfiguration.h"
#include "AliHLTConfigurationHandler.h"
#include "Riostream.h"
#endif

/**
 * Generates some sample input data and writes it into 8 files named
 * testInputFile1.root ... testInputFile8.root
 */
void GenerateInputData()
{
	if (gClassTable->GetID("AliHLTGlobalTriggerComponent") < 0)
	{
		gSystem->Load("libAliHLTUtil.so");
		gSystem->Load("libAliHLTTRD.so");
		gSystem->Load("libAliHLTTrigger.so");
	}

	AliHLTReadoutList readoutList1("TPC");
	AliHLTTriggerDomain triggerDomain1;
	triggerDomain1.Add("CLUSTERS", "TPC ");
	triggerDomain1.Add("TRACKS", "TPC ");
	triggerDomain1.Add(readoutList1);
	AliHLTTriggerDecision decision1(true, "triggerTPC", triggerDomain1, "TPC has data");
	
	AliHLTReadoutList readoutList2a("MUONTRK");
	AliHLTReadoutList readoutList2b("MUONTRG");
	AliHLTTriggerDomain triggerDomain2;
	triggerDomain2.Add("TRACKS", "MUON");
	triggerDomain2.Add(readoutList2a);
	triggerDomain2.Add(readoutList2b);
	AliHLTTriggerDecision decision2(true, "triggerMUON", triggerDomain2, "MUON has data");
	
	AliHLTReadoutList readoutList3("ITSSSD");
	AliHLTTriggerDomain triggerDomain3;
	triggerDomain3.Add("*******", "SSD ");
	triggerDomain3.Add(readoutList3);
	AliHLTTriggerDecision decision3(true, "triggerSSD", triggerDomain3, "SSD has data");
	
	TFile* file = new TFile("testInputFile1.root", "RECREATE");
	decision1.Write("triggerTPC");
	decision2.Write("triggerMUON");
	decision3.Write("triggerSSD");
	delete file;
	
	file = new TFile("testInputFile2.root", "RECREATE");
	decision1.Write("triggerTPC");
	delete file;
	
	file = new TFile("testInputFile3.root", "RECREATE");
	decision2.Write("triggerMUON");
	delete file;
	
	file = new TFile("testInputFile4.root", "RECREATE");
	decision3.Write("triggerSSD");
	delete file;
	
	file = new TFile("testInputFile5.root", "RECREATE");
	decision1.Write("triggerTPC");
	decision2.Write("triggerMUON");
	delete file;
	
	file = new TFile("testInputFile6.root", "RECREATE");
	decision1.Write("triggerTPC");
	decision3.Write("triggerSSD");
	delete file;
	
	file = new TFile("testInputFile7.root", "RECREATE");
	decision2.Write("triggerMUON");
	decision3.Write("triggerSSD");
	delete file;
	
	file = new TFile("testInputFile8.root", "RECREATE");
	delete file;
}

/**
 * Runs a small global trigger test chain with the different configuration as specified
 * in TriggerConfig.C.
 * \param config  The configuration version to pass to TriggerConfig.C
 * \param usecint  If true then the global trigger component uses CINT to interpret
 *     the code rather than compiling it.
 * \param debug  If true then the global trigger component generates extra debug
 *     statements in the on the fly AliHLTGlobalTriggerImp_*.cxx file.
 * \param numOfEvents  The number of events to run the chain for.
 * \param customClass  Names the custom class that should be loaded from the file
 *     <i>\<customClass\>.cxx</i>. This is useful for debugging only. i.e. you can
 *     edit a generated logic file and test it by hand.
 */
void RunTrigger(int config = 0, bool usecint = false, bool debug = false, int numOfEvents = 8, const char* customClass = NULL)
{
	AliHLTSystem sys;
	sys.fpConfigurationHandler->SetLocalLoggingLevel(kHLTLogAll);
	sys.LoadComponentLibraries("libAliHLTUtil.so");
	sys.LoadComponentLibraries("libAliHLTTRD.so");
	sys.LoadComponentLibraries("libAliHLTTrigger.so");
	if (debug)
	{
		AliLog::SetGlobalLogLevel(AliLog::kMaxType);
		sys.SetGlobalLoggingLevel(kHLTLogAll);
	}
	
	TString cmdline = "-datatype ROOTTOBJ 'HLT ' ";
	for (int i = 1; i <= 8; i++)
	{
		if (i > 1) cmdline += " -nextevent";
		cmdline += Form(" -datafile testInputFile%d.root", i);
	}
	AliHLTConfiguration pub("pub", "ROOTFilePublisher", NULL, cmdline.Data());
	
	cmdline = Form("-config TriggerConfig.C(%d) -includepath $ALICE_ROOT/include -includepath $ALICE_ROOT/HLT/BASE"
		" -includepath $ALICE_ROOT/HLT/trigger -include AliHLTEventSummary.h",
		config
		);
	if (customClass != NULL) cmdline += Form(" -usecode %s.cxx %s", customClass, customClass);
	if (usecint) cmdline += " -cint";
	if (debug) cmdline += " -debug";
	AliHLTConfiguration proc("proc", "HLTGlobalTrigger", "pub", cmdline.Data());
	
	AliHLTConfiguration sink("sink", "ROOTFileWriter", "proc", "-datafile testOutputFile.root -concatenate-events");
	
	sys.BuildTaskList("sink");
	sys.Run(numOfEvents);
}

/**
 * Checks that a particular decision is as expected and prints error messages
 * if it is not.
 * \param testName  The name of the test being run.
 * \param eventNum  The number of the event being checked.
 * \param decision  The global trigger decision being checked.
 * \param expectedResult  The expected global trigger result.
 * \param expectedDomain  The expected resulting global trigger domain.
 * \param expectedDescription  The expected resulting trigger description.
 * \returns true if the decision is as expected.
 */
bool Check(
		const char* testName,
		int eventNum,
		AliHLTGlobalTriggerDecision* decision,
		bool expectedResult,
		AliHLTTriggerDomain expectedDomain,
		TString expectedDescription
	)
{
	if (decision == NULL)
	{
		cerr << "ERROR (Test: " << testName
		     << "): No decision found where expected for event "
		     << eventNum << "." << endl;
		return false;
	}
	if (decision->Result() != expectedResult)
	{
		cerr << "ERROR (Test: " << testName
		     << "): The result does not match the expected value for event "
		     << eventNum << ". Got " << decision->Result() << " but expected "
		     << expectedResult << "." << endl;
		return false;
	}
	if (decision->TriggerDomain() != expectedDomain)
	{
		cerr << "ERROR (Test: " << testName
		     << "): The domain does not match the expected value for event "
		     << eventNum << ". Got:" << endl;
		decision->TriggerDomain().Print();
		cout << "but expected:" << endl;
		expectedDomain.Print();
		return false;
	}
	if (decision->Description() != expectedDescription)
	{
		cerr << "ERROR (Test: " << testName
		     << "): The description does not match the expected value for event "
		     << eventNum << ". Got '" << decision->Description() << "' but expected '"
		     << expectedDescription << "'." << endl;
		return false;
	}
	return true;
}


/// Routine for checking the result of the PriorityGroupTestConfig() config in TriggerConfig.C
bool CheckPriorityGroupTestConfig(const char* testName = "Priority group config")
{
	AliHLTGlobalTriggerDecision* decision = NULL;
	bool result = false;
	
	AliHLTTriggerDomain domainTPC("CLUSTERS:TPC ,TRACKS:TPC ");
	domainTPC.Add(AliHLTReadoutList("TPC"));
	AliHLTTriggerDomain domainMUON("TRACKS:MUON");
	domainMUON.Add(AliHLTReadoutList("MUONTRK"));
	domainMUON.Add(AliHLTReadoutList("MUONTRG"));
	AliHLTTriggerDomain domainSSD("*******:SSD ");
	domainSSD.Add(AliHLTReadoutList("ITSSSD"));

	TFile* file = new TFile("testOutputFile.root", "READ");
	
	// Triggers in events (i.e. input triggers):
	// event 1: TPC MUON SSD
	// event 2: TPC
	// event 3: MUON
	// event 4: SSD
	// event 5: TPC MUON
	// event 6: TPC SSD
	// event 7: MUON SSD
	// event 8: (empty)
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;1"));
	result = Check(testName, 1, decision, true, domainSSD, "Fast SSD trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;2"));
	result = Check(testName, 2, decision, true, domainTPC, "TPC trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;3"));
	result = Check(testName, 3, decision, true, domainMUON, "MUON trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;4"));
	result = Check(testName, 4, decision, true, domainSSD, "MUON trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;5"));
	result = Check(testName, 5, decision, true, domainMUON, "MUON trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;6"));
	result = Check(testName, 6, decision, true, domainSSD, "Fast SSD trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;7"));
	result = Check(testName, 7, decision, true, domainMUON | domainSSD, "MUON trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;8"));
	result = Check(testName, 8, decision, false, AliHLTTriggerDomain(), "No trigger");
	if (! result) goto cleanup;
	
	delete file;
	return true;
	
cleanup:
	if (decision != NULL)
	{
		cout << "========== Dumping incorrect decision ========== " << endl;
		decision->Print();
	}
	delete file;
	return false;
}


/// Routine for checking the result of the SingleGroupTestConfig() config in TriggerConfig.C
bool CheckSingleGroupTestConfig(const char* testName = "Single group config")
{
	AliHLTGlobalTriggerDecision* decision = NULL;
	bool result = false;
	
	AliHLTTriggerDomain domainTPC("CLUSTERS:TPC ,TRACKS:TPC ");
	domainTPC.Add(AliHLTReadoutList("TPC"));
	AliHLTTriggerDomain domainMUON("TRACKS:MUON");
	domainMUON.Add(AliHLTReadoutList("MUONTRK"));
	domainMUON.Add(AliHLTReadoutList("MUONTRG"));
	AliHLTTriggerDomain domainSSD("*******:SSD ");
	domainSSD.Add(AliHLTReadoutList("ITSSSD"));

	TFile* file = new TFile("testOutputFile.root", "READ");
	
	// Triggers in events (i.e. input triggers):
	// event 1: TPC MUON SSD
	// event 2: TPC
	// event 3: MUON
	// event 4: SSD
	// event 5: TPC MUON
	// event 6: TPC SSD
	// event 7: MUON SSD
	// event 8: (empty)
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;1"));
	result = Check(testName, 1, decision, true, domainTPC | domainMUON | domainSSD, "TPC trigger,MUON trigger,SSD trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;2"));
	result = Check(testName, 2, decision, true, domainTPC, "TPC trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;3"));
	result = Check(testName, 3, decision, true, domainMUON, "MUON trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;4"));
	result = Check(testName, 4, decision, true, domainSSD, "SSD trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;5"));
	result = Check(testName, 5, decision, true, domainTPC | domainMUON, "TPC trigger,MUON trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;6"));
	result = Check(testName, 6, decision, true, domainTPC | domainSSD, "TPC trigger,SSD trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;7"));
	result = Check(testName, 7, decision, true, domainMUON | domainSSD, "MUON trigger,SSD trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;8"));
	result = Check(testName, 8, decision, false, AliHLTTriggerDomain(), "");
	if (! result) goto cleanup;
	
	delete file;
	return true;
	
cleanup:
	if (decision != NULL)
	{
		cout << "========== Dumping incorrect decision ========== " << endl;
		decision->Print();
	}
	delete file;
	return false;
}


typedef bool (*CheckFunctionType)(const char* testName);

/**
 * This will check the results of a global trigger run with a particular checking
 * routine and for different combinations of -cint and -debug flags passed to the
 * global trigger component.
 * It is important to check these flag combinations to make sure that everything
 * is functioning the same under all the different combinations as it should.
 * \param function  The checking routine to use.
 * \param version  The trigger menu configuration version to use in <i>RunTrigger</i>.
 * \param testName  The name of the test being run.
 * \param numOfEvents  The number of events to run the chain for.
 * \param customClass  Name of the custom class as passed to <i>RunTrigger</i>.
 * \returns true if the different checks succeeded and false otherwise.
 */
bool CheckDifferentModes(CheckFunctionType function, int version, const char* testName, int numOfEvents = 8, const char* customClass = NULL)
{
	TString name = testName;
	RunTrigger(version, false, false, numOfEvents, customClass);
	if (! function(testName)) return false;
	gSystem->Exec("rm -f testOutputFile.root");  // Cleanup output file for next test.
	
	name = testName;
	name += " interpreted with CINT";
	RunTrigger(version, true, false, numOfEvents, customClass);
	if (! function(testName)) return false;
	gSystem->Exec("rm -f testOutputFile.root");  // Cleanup output file for next test.
	
	name = testName;
	name += " in debug mode";
	RunTrigger(version, false, true, numOfEvents, customClass);
	if (! function(testName)) return false;
	gSystem->Exec("rm -f testOutputFile.root");  // Cleanup output file for next test.
	
	name = testName;
	name += " interpreted with CINT in debug mode";
	RunTrigger(version, true, true, numOfEvents, customClass);
	if (! function(testName)) return false;
	gSystem->Exec("rm -f testOutputFile.root");  // Cleanup output file for next test.
	
	return true;
}

/**
 * Runs several tests for the AliHLTGlobalTriggerComponent class.
 * We specifically test if the global trigger menu configuration is interpreted
 * correctly and the trigger logic generated correctly on the fly.
 * \param numOfEvents  The number of events to run the chain for.
 * \param customClass  Name of the custom class as passed to <i>CheckDifferentModes</i>.
 * \returns true if the different checks succeeded and false otherwise.
 */
bool testGlobalTriggerComponent(int numOfEvents = 8, const char* customClass = NULL)
{
	GenerateInputData();
	
	if (! CheckDifferentModes(CheckPriorityGroupTestConfig, 0, "Priority group config", numOfEvents, customClass)) return false;
	if (! CheckDifferentModes(CheckSingleGroupTestConfig, 1, "Single group config", numOfEvents, customClass)) return false;
	
	// Cleanup all temporary files generated.
	gSystem->Exec("rm -f testOutputFile.root testInputFile*.root AliHLTGlobalTriggerImpl*");
	return true;
}



#ifndef __MAKECINT__

int main(int /*argc*/, const char** /*argv*/)
{
	gSystem->Exec("if test ! -f TriggerConfig.C ; then cp $ALICE_ROOT/HLT/trigger/test/TriggerConfig.C ./; fi");
	bool resultOk = testGlobalTriggerComponent();
	if (not resultOk) return 1;
	return 0;
}

#endif // __MAKECINT__
