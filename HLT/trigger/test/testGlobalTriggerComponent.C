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
 *
 * This macro can also be used to debug the configuration that fails.
 * If a configuration test fails then the global trigger logic implementation
 * is left in a generated file with a name of the form \<AliHLTGlobalTriggerImpl_*.cxx\>.
 * For example a file named like: AliHLTGlobalTriggerImpl_08869e64_c54b_11de_9717_0101007fbeef.cxx
 * One can make manual modifications to this file and then rerun the test with
 * these manual modifications with the following command in aliroot:
 *  .x testGlobalTriggerComponent.C+(\<configVersion\>,"\<AliHLTGlobalTriggerImpl_*\>")
 * where \<configVersion\> is the appropriate config version number as passed to the
 * TriggerConfig.C file to initialise the configuration we want to test. Also take note
 * that we only specify the root of the file name without the .cxx file name extention.
 * For our example file name the command to execute in aliroot would be:
 *  .x testGlobalTriggerComponent.C+(2,"AliHLTGlobalTriggerImpl_08869e64_c54b_11de_9717_0101007fbeef")
 * where we are testing the 2nd trigger configuration in TriggerConfig.C.
 */

#if defined(__CINT__) && (! defined(__MAKECINT__))
#error This macro must be compiled. Try running as testGlobalTriggerComponent.C++, but remember to load the libAliHLTTrigger.so library first.
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
#include "AliHLTScalars.h"
#include "AliHLTSystem.h"
#include "AliHLTConfiguration.h"
#include "Riostream.h"
#endif

/**
 * Generates some sample input data and writes it into 8 files named
 * testInputFile1.root ... testInputFile8.root
 */
void GenerateInputData()
{
	bool loadedLibs = false;
	if (gClassTable->GetID("AliHLTGlobalTriggerComponent") < 0)
	{
		gSystem->Load("libAliHLTUtil.so");
		gSystem->Load("libAliHLTTRD.so");
		gSystem->Load("libAliHLTMUON.so");
		gSystem->Load("libAliHLTTrigger.so");
		loadedLibs = true;
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
	
	AliHLTScalars summary1;
	summary1.GetScalar("TrigClass").Value(0x1);
	
	AliHLTScalars summary2;
	summary2.GetScalar("TrigClass").Value(0x2);
	
	TFile* file = new TFile("testInputFile1.root", "RECREATE");
	decision1.Write("triggerTPC");
	decision2.Write("triggerMUON");
	decision3.Write("triggerSSD");
	summary2.Write("summary");
	delete file;
	
	file = new TFile("testInputFile2.root", "RECREATE");
	decision1.Write("triggerTPC");
	summary2.Write("summary");
	delete file;
	
	file = new TFile("testInputFile3.root", "RECREATE");
	decision2.Write("triggerMUON");
	summary1.Write("summary");
	delete file;
	
	file = new TFile("testInputFile4.root", "RECREATE");
	decision3.Write("triggerSSD");
	summary2.Write("summary");
	delete file;
	
	file = new TFile("testInputFile5.root", "RECREATE");
	decision1.Write("triggerTPC");
	decision2.Write("triggerMUON");
	summary1.Write("summary");
	delete file;
	
	file = new TFile("testInputFile6.root", "RECREATE");
	decision1.Write("triggerTPC");
	decision3.Write("triggerSSD");
	summary2.Write("summary");
	delete file;
	
	file = new TFile("testInputFile7.root", "RECREATE");
	decision2.Write("triggerMUON");
	decision3.Write("triggerSSD");
	summary1.Write("summary");
	delete file;
	
	file = new TFile("testInputFile8.root", "RECREATE");
	delete file;
	
	/*FIXME This has stopped working in AliRoot, it causes a segfault now within the AliHLTReadoutList vtable.
	if (loadedLibs)
	{
		gSystem->Unload("libAliHLTTrigger.so");
		gSystem->Unload("libAliHLTMUON.so");
		gSystem->Unload("libAliHLTTRD.so");
		gSystem->Unload("libAliHLTUtil.so");
	}
	*/
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
	sys.ScanOptions("ECS=CTP_TRIGGER_CLASS=00:TRIGGER-ALL:00-01-02-03-04-05-06-07-08-09-10-11-12-13-14-15-16-17");
	sys.LoadComponentLibraries("libAliHLTUtil.so");
	sys.LoadComponentLibraries("libAliHLTTRD.so");
	sys.LoadComponentLibraries("libAliHLTMUON.so");
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
	
	cmdline = Form("-config $ALICE_ROOT/HLT/trigger/test/TriggerConfig.C(%d)"
		" -includepath $ALICE_ROOT/include -includepath $ALICE_ROOT/HLT/BASE"
		" -includepath $ALICE_ROOT/HLT/trigger -include AliHLTScalars.h",
		config
		);
	if (customClass != NULL) cmdline += Form(" -usecode %s.cxx %s", customClass, customClass);
	if (usecint) cmdline += " -cint";
	if (debug) cmdline += " -debug";
	AliHLTConfiguration proc("proc", "HLTGlobalTrigger", "pub", cmdline.Data());
	
	AliHLTConfiguration sink("sink", "ROOTFileWriter", "proc", "-datafile testOutputFile.root -concatenate-events");
	
	sys.BuildTaskList("sink");
	sys.Run(
		numOfEvents,
		1,   // Stop chain at end of run.
		0x1, // Active CTP trigger mask.
		0,   // Time stamp.
		0    // Event type.
	);
}

/**
 * This method calls the RunTrigger method in an independant aliroot process.
 * This is necessary since we get memory corruption if we run too many instances of
 * AliHLTSystem in the same process.
 */
void CallRunTrigger(
		int config = 0, bool usecint = false, bool debug = false,
		int numOfEvents = 8, const char* customClass = NULL,
		bool showOutput = false
	)
{
	const char* redirection = "> /dev/null";
	const char* classNameString = "NULL";
	if (showOutput) redirection = "";
	if (customClass != NULL) classNameString = Form("\"%s\"", customClass);
	const char* command = Form(
			"aliroot %s <<EOF\n"
			"gSystem->Load(\"libAliHLTUtil.so\");\n"
			"gSystem->Load(\"libAliHLTTRD.so\");\n"
			"gSystem->Load(\"libAliHLTMUON.so\");\n"
			"gSystem->Load(\"libAliHLTTrigger.so\");\n"
			"gSystem->SetIncludePath(\"-I${ALICE_ROOT}/include"
			" -I${ALICE_ROOT}/HLT/BASE -I${ALICE_ROOT}/HLT/trigger\");\n"
			".L $ALICE_ROOT/HLT/trigger/test/testGlobalTriggerComponent.C+\n"
			"RunTrigger(%d,%d,%d,%d,%s);\n"
			"EOF\n",
			redirection,
			config,
			usecint,
			debug,
			numOfEvents,
			classNameString
		);
	gSystem->Exec(command);
}

/**
 * Runs a global trigger test chain to test L0 software triggers.
 * \param config  The configuration version to pass to TriggerConfig.C
 * \param usecint  If true then the global trigger component uses CINT to interpret
 *     the code rather than compiling it.
 * \param debug  If true then the global trigger component generates extra debug
 *     statements in the on the fly AliHLTGlobalTriggerImp_*.cxx file.
 * \param customClass  Names the custom class that should be loaded from the file
 *     <i>\<customClass\>.cxx</i>. This is useful for debugging only. i.e. you can
 *     edit a generated logic file and test it by hand.
 */
void RunTriggerSW(int config = 0, bool usecint = false, bool debug = false, const char* customClass = NULL)
{
	AliHLTSystem sys;
	sys.ScanOptions("ECS=CTP_TRIGGER_CLASS=00:TRIGGER-ALL:00-01-02-03-04-05-06-07-08-09-10-11-12-13-14-15-16-17");
	sys.LoadComponentLibraries("libAliHLTUtil.so");
	sys.LoadComponentLibraries("libAliHLTTRD.so");
	sys.LoadComponentLibraries("libAliHLTMUON.so");
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
	
	cmdline = Form("-config $ALICE_ROOT/HLT/trigger/test/TriggerConfig.C(%d)"
		" -includepath $ALICE_ROOT/include -includepath $ALICE_ROOT/HLT/BASE"
		" -includepath $ALICE_ROOT/HLT/trigger -include AliHLTScalars.h"
		" -process-all-events",
		config
		);
	if (customClass != NULL) cmdline += Form(" -usecode %s.cxx %s", customClass, customClass);
	if (usecint) cmdline += " -cint";
	if (debug) cmdline += " -debug";
	AliHLTConfiguration proc("proc", "HLTGlobalTrigger", "pub", cmdline.Data());
	
	AliHLTConfiguration sink("sink", "ROOTFileWriter", "proc", "-datafile testOutputFile.root -concatenate-events");
	
	sys.BuildTaskList("sink");
	sys.Run(
		1,   // Number of events to process.
		0,   // Stop chain at end of run.
		0x1, // Active CTP trigger mask.
		0,   // Time stamp.
		gkAliEventTypeSoftware  // Event type.
	);
	sys.Run(
		1,   // Number of events to process.
		0,   // Stop chain at end of run.
		0x1, // Active CTP trigger mask.
		0,   // Time stamp.
		gkAliEventTypeCalibration  // Event type.
	);
	sys.Run(
		1,   // Number of events to process.
		0,   // Stop chain at end of run.
		0x1, // Active CTP trigger mask.
		0,   // Time stamp.
		0    // Event type.
	);
	sys.Run(
		1,   // Number of events to process.
		1,   // Stop chain at end of run.
		0x1, // Active CTP trigger mask.
		0,   // Time stamp.
		gkAliEventTypeSoftware,  // Event type.
		AliHLTReadoutList::kTPC
	);
}

/**
 * This method calls the RunTriggerSW method in an independant aliroot process.
 * This is necessary since we get memory corruption if we run too many instances of
 * AliHLTSystem in the same process.
 */
void CallRunTriggerSW(
		int config = 0, bool usecint = false, bool debug = false,
		const char* customClass = NULL, bool showOutput = false
	)
{
	const char* redirection = "> /dev/null";
	const char* classNameString = "NULL";
	if (showOutput) redirection = "";
	if (customClass != NULL) classNameString = Form("\"%s\"", customClass);
	const char* command = Form(
			"aliroot %s <<EOF\n"
			"gSystem->Load(\"libAliHLTUtil.so\");\n"
			"gSystem->Load(\"libAliHLTTRD.so\");\n"
			"gSystem->Load(\"libAliHLTMUON.so\");\n"
			"gSystem->Load(\"libAliHLTTrigger.so\");\n"
			"gSystem->SetIncludePath(\"-I${ALICE_ROOT}/include"
			" -I${ALICE_ROOT}/HLT/BASE -I${ALICE_ROOT}/HLT/trigger\");\n"
			".L $ALICE_ROOT/HLT/trigger/test/testGlobalTriggerComponent.C+\n"
			"RunTriggerSW(%d,%d,%d,%s);\n"
			"EOF\n",
			redirection,
			config,
			usecint,
			debug,
			classNameString
		);
	gSystem->Exec(command);
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
		cerr << "but expected:" << endl;
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


/// Routine for checking the result of the PrescalarTestConfig() config in TriggerConfig.C
bool CheckPrescalarTestConfig(const char* testName = "Prescalar config")
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
	AliHLTTriggerDomain defaultDomain("*******:HLT ");

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
	result = Check(testName, 1, decision, true, domainTPC, "TPC trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;2"));
	result = Check(testName, 2, decision, false, defaultDomain, "No trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;3"));
	result = Check(testName, 3, decision, true, domainMUON, "MUON trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;4"));
	result = Check(testName, 4, decision, false, defaultDomain, "No trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;5"));
	result = Check(testName, 5, decision, true, domainMUON, "MUON trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;6"));
	result = Check(testName, 6, decision, true, domainTPC, "TPC trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;7"));
	result = Check(testName, 7, decision, true, domainMUON, "MUON trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;8"));
	result = Check(testName, 8, decision, false, defaultDomain, "No trigger");
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


/// Routine for checking the result of the SymbolTestConfig() config in TriggerConfig.C
bool CheckSymbolTestConfig(const char* testName = "Symbol config")
{
	AliHLTGlobalTriggerDecision* decision = NULL;
	bool result = false;
	
	AliHLTTriggerDomain domainAll("*******:***,-DAQRDOUT:TST");
	AliHLTTriggerDomain domainPHOS("CLUSTERS:PHOS,TRACKS:PHOS");
	AliHLTTriggerDomain domainTPC("CLUSTERS:TPC ,TRACKS:TPC ");
	domainTPC.Add(AliHLTReadoutList("TPC"));
	AliHLTTriggerDomain domainMUON("TRACKS:MUON");
	domainMUON.Add(AliHLTReadoutList("MUONTRK"));
	domainMUON.Add(AliHLTReadoutList("MUONTRG"));
	AliHLTTriggerDomain domainSSD("*******:SSD ");
	domainSSD.Add(AliHLTReadoutList("ITSSSD"));

	TFile* file = new TFile("testOutputFile.root", "READ");
	
	// Triggers in events (i.e. input triggers) and trigger classes in AliHLTScalars:
	// event 1: TPC MUON SSD, 0x2
	// event 2: TPC         , 0x2
	// event 3: MUON        , 0x1
	// event 4: SSD         , 0x2
	// event 5: TPC MUON    , 0x1
	// event 6: TPC SSD     , 0x2
	// event 7: MUON SSD    , 0x1
	// event 8: (empty)     , 0x0
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;1"));
	result = Check(testName, 1, decision, true, domainAll - AliHLTTriggerDomain("DAQRDOUT:EMC"), "Pass through");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;2"));
	result = Check(testName, 2, decision, true, domainTPC | domainPHOS, "Trigger class 2");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;3"));
	result = Check(testName, 3, decision, false, AliHLTTriggerDomain(), "");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;4"));
	result = Check(testName, 4, decision, true, domainPHOS, "Trigger class 2");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;5"));
	result = Check(testName, 5, decision, true, domainAll - AliHLTTriggerDomain("DAQRDOUT:EMC"), "Pass through");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;6"));
	result = Check(testName, 6, decision, true, domainTPC | domainPHOS, "Trigger class 2");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;7"));
	result = Check(testName, 7, decision, false, AliHLTTriggerDomain(), "");
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


/// Routine for checking the result of the ComplexTestConfig() config in TriggerConfig.C
bool CheckComplexTestConfig(const char* testName = "Complex config")
{
	AliHLTGlobalTriggerDecision* decision = NULL;
	bool result = false;
	
	AliHLTTriggerDomain domainAll("*******:***,-DAQRDOUT:TST,-DAQRDOUT:EMC");
	AliHLTTriggerDomain domainPHOS("CLUSTERS:PHOS,TRACKS:PHOS");
	AliHLTTriggerDomain domainTPC("CLUSTERS:TPC ,TRACKS:TPC ");
	domainTPC.Add(AliHLTReadoutList("TPC"));
	AliHLTTriggerDomain domainMUON("TRACKS:MUON");
	domainMUON.Add(AliHLTReadoutList("MUONTRK"));
	domainMUON.Add(AliHLTReadoutList("MUONTRG"));
	AliHLTTriggerDomain domainSSD("*******:SSD ");
	domainSSD.Add(AliHLTReadoutList("ITSSSD"));

	TFile* file = new TFile("testOutputFile.root", "READ");
	
	// Triggers in events (i.e. input triggers) and trigger classes in AliHLTScalars:
	// event 1: TPC MUON SSD, 0x2
	// event 2: TPC         , 0x2
	// event 3: MUON        , 0x1
	// event 4: SSD         , 0x2
	// event 5: TPC MUON    , 0x1
	// event 6: TPC SSD     , 0x2
	// event 7: MUON SSD    , 0x1
	// event 8: (empty)     , 0x0
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;1"));
	result = Check(testName, 1, decision, true, domainAll, "Pass through");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;2"));
	result = Check(testName, 2, decision, true, domainTPC, "Slow trigger");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;3"));
	result = Check(testName, 3, decision, true, domainMUON, "MUON trigger 2");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;4"));
	result = Check(testName, 4, decision, true, domainSSD | domainPHOS, "SSD trigger 2");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;5"));
	result = Check(testName, 5, decision, true, domainMUON, "MUON trigger 2");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;6"));
	result = Check(testName, 6, decision, true, domainSSD | domainPHOS, "SSD trigger 2");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;7"));
	result = Check(testName, 7, decision, true, domainSSD | domainMUON, "SSD trigger 1,MUON trigger 1");
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;8"));
	result = Check(testName, 8, decision, true, domainAll, "Pass through");
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
 * \param showOutput  If true then the output from the RunTrigger method is not suppressed.
 * \returns true if the different checks succeeded and false otherwise.
 */
bool CheckDifferentModes(
		CheckFunctionType function, int version, const char* testName,
		int numOfEvents = 8, const char* customClass = NULL,
		bool showOutput = false
	)
{
	TString name = testName;
	name += " in debug mode";
	cout << "#################### Running test: " << name.Data() << " ####################" << endl;
	CallRunTrigger(version, false, true, numOfEvents, customClass, showOutput);
	if (! function(testName)) return false;
	gSystem->Exec("rm -f testOutputFile.root");  // Cleanup output file for next test.
	
	name = testName;
	cout << "#################### Running test: " << name.Data() << " ####################" << endl;
	CallRunTrigger(version, false, false, numOfEvents, customClass, showOutput);
	if (! function(testName)) return false;
	gSystem->Exec("rm -f testOutputFile.root");  // Cleanup output file for next test.
	
	name = testName;
	name += " interpreted with CINT in debug mode";
	cout << "#################### Running test: " << name.Data() << " ####################" << endl;
	CallRunTrigger(version, true, true, numOfEvents, customClass, showOutput);
	if (! function(testName)) return false;
	gSystem->Exec("rm -f testOutputFile.root");  // Cleanup output file for next test.
	
	name = testName;
	name += " interpreted with CINT";
	cout << "#################### Running test: " << name.Data() << " ####################" << endl;
	CallRunTrigger(version, true, false, numOfEvents, customClass, showOutput);
	if (! function(testName)) return false;
	gSystem->Exec("rm -f testOutputFile.root");  // Cleanup output file for next test.
	
	return true;
}

/**
 * This routine is used to check if the global Trigger counters are correct.
 * \param testName  The name of the test being run.
 * \param eventNum  The number of the event being checked.
 * \param decision  The global trigger decision being checked.
 * \param expectedCounters  The expected counters.
 * \returns true if the decision is as expected.
 */
bool CheckCounters(
		const char* testName,
		int eventNum,
		AliHLTGlobalTriggerDecision* decision,
		const TArrayL64& expectedCounters
	)
{
	if (decision->Counters().GetSize() != expectedCounters.GetSize())
	{
		cerr << "ERROR (Test: " << testName
		     << "): The result does not have the required number of counters for event "
		     << eventNum << ". Got " << decision->Counters().GetSize() << " but expected "
		     << expectedCounters.GetSize() << "." << endl;
		return false;
	}
	for (Int_t i = 0; i < expectedCounters.GetSize(); ++i)
	{
		if (decision->Counters()[i] != expectedCounters[i])
		{
			cerr << "ERROR (Test: " << testName
			     << "): The result does not have the correct counter value for event "
			     << eventNum << ". Got a value " << decision->Counters()[i]
			     << " for counter " << i << ", but expected a value of "
			     << expectedCounters[i] << "." << endl;
			return false;
		}
	}
	return true;
}

/// Routine for checking the result of the SoftwareTriggersTestConfig() config in TriggerConfig.C
bool CheckSoftwareTriggerTestConfig(const char* testName = "Software trigger config")
{
	AliHLTGlobalTriggerDecision* decision = NULL;
	bool result = false;
	
	AliHLTTriggerDomain domainPHOS("*******:PHOS");
	AliHLTTriggerDomain domainSPD("*******:SPD");
	AliHLTTriggerDomain domainTPC("DAQRDOUT:TPC");
	AliHLTTriggerDomain domainMUON("TRACKS:MUON");
	domainMUON.Add(AliHLTReadoutList("MUONTRK MUONTRG"));

	TFile* file = new TFile("testOutputFile.root", "READ");
	TArrayL64 expectedCounters;
	expectedCounters.Set(6);
	
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;1"));
	result = Check(testName, 1, decision, true, AliHLTTriggerDomain(), "Start of data");
	if (! result) goto cleanup;
	expectedCounters[0] = 1; expectedCounters[5] = 1;
	result = CheckCounters(testName, 1, decision, expectedCounters);
	if (! result) goto cleanup;
	
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;2"));
	result = Check(testName, 2, decision, true, domainSPD, "Software trigger");
	if (! result) goto cleanup;
	expectedCounters[2] = 1; expectedCounters[5] = 2;
	result = CheckCounters(testName, 2, decision, expectedCounters);
	if (! result) goto cleanup;
	
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;3"));
	result = Check(testName, 3, decision, true, domainPHOS, "Calibration trigger");
	if (! result) goto cleanup;
	expectedCounters[3] = 1; expectedCounters[5] = 3;
	result = CheckCounters(testName, 3, decision, expectedCounters);
	if (! result) goto cleanup;
	
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;4"));
	result = Check(testName, 4, decision, true, domainMUON, "MUON trigger");
	if (! result) goto cleanup;
	expectedCounters[4] = 1; expectedCounters[5] = 4;
	result = CheckCounters(testName, 4, decision, expectedCounters);
	if (! result) goto cleanup;
	
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;5"));
	result = Check(testName, 5, decision, true, domainSPD|domainTPC, "Software trigger");
	if (! result) goto cleanup;
	expectedCounters[2] = 2; expectedCounters[5] = 5;
	result = CheckCounters(testName, 5, decision, expectedCounters);
	if (! result) goto cleanup;
	
	decision = dynamic_cast<AliHLTGlobalTriggerDecision*>(file->Get("HLTGlobalTrigger;6"));
	result = Check(testName, 6, decision, true, AliHLTTriggerDomain(), "End of data");
	if (! result) goto cleanup;
	expectedCounters[1] = 1; expectedCounters[5] = 6;
	result = CheckCounters(testName, 6, decision, expectedCounters);
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

/**
 * This method performs the same task as for CheckDifferentModes, but trying to
 * test the behaviour of the global HLT trigger component with L0 software triggers.
 * \param version  The trigger menu configuration version to use in <i>RunTrigger</i>.
 * \param testName  The name of the test being run.
 * \param customClass  Name of the custom class as passed to <i>RunTrigger</i>.
 * \param showOutput  If true then the output from the RunTriggerSW method is not suppressed.
 * \returns true if the different checks succeeded and false otherwise.
 */
bool CheckDifferentSWTestModes(
		int version, const char* testName,
		const char* customClass = NULL, bool showOutput = false
	)
{
	TString name = testName;
	name += " in debug mode";
	cout << "#################### Running test: " << name.Data() << " ####################" << endl;
	CallRunTriggerSW(version, false, true, customClass, showOutput);
	if (! CheckSoftwareTriggerTestConfig(testName)) return false;
	gSystem->Exec("rm -f testOutputFile.root");  // Cleanup output file for next test.
	
	name = testName;
	cout << "#################### Running test: " << name.Data() << " ####################" << endl;
	CallRunTriggerSW(version, false, false, customClass, showOutput);
	if (! CheckSoftwareTriggerTestConfig(testName)) return false;
	gSystem->Exec("rm -f testOutputFile.root");  // Cleanup output file for next test.
	
	name = testName;
	name += " interpreted with CINT in debug mode";
	cout << "#################### Running test: " << name.Data() << " ####################" << endl;
	CallRunTriggerSW(version, true, true, customClass, showOutput);
	if (! CheckSoftwareTriggerTestConfig(testName)) return false;
	gSystem->Exec("rm -f testOutputFile.root");  // Cleanup output file for next test.
	
	name = testName;
	name += " interpreted with CINT";
	cout << "#################### Running test: " << name.Data() << " ####################" << endl;
	CallRunTriggerSW(version, true, false, customClass, showOutput);
	if (! CheckSoftwareTriggerTestConfig(testName)) return false;
	gSystem->Exec("rm -f testOutputFile.root");  // Cleanup output file for next test.
	
	return true;
}

/**
 * Runs several tests for the AliHLTGlobalTriggerComponent class.
 * We specifically test if the global trigger menu configuration is interpreted
 * correctly and the trigger logic generated correctly on the fly.
 * \param configVersion  The appropriate version number of the config being tested
 *     which is passed to TriggerConfig.C.
 * \param customClass  Name of the custom class as passed to <i>CheckDifferentModes</i>.
 * \returns true if the different checks succeeded and false otherwise.
 * \param numOfEvents  The number of events to run the chain for.
 * \returns true if all the tests succeeded and false otherwise.
 */
bool testGlobalTriggerComponent(int configVersion = -1, const char* customClass = NULL, int numOfEvents = 8)
{
	GenerateInputData();
	
	if (configVersion != -1)
	{
		CheckFunctionType function = NULL;
		switch (configVersion)
		{
		case 0: function = CheckPriorityGroupTestConfig; break;
		case 1: function = CheckSingleGroupTestConfig; break;
		case 2: function = CheckPrescalarTestConfig; break;
		case 3: function = CheckSymbolTestConfig; break;
		case 4: function = CheckComplexTestConfig; break;
		case 5: break;
		default:
			cerr << "ERROR: Invalid value for configVersion specified." << endl;
			return false;
		}
		bool result = false;
		if (configVersion != 5)
		{
			result = CheckDifferentModes(
					function,
					configVersion,
					Form("Config version %d", configVersion),
					numOfEvents,
					customClass,
					true
				);
		}
		else
		{
			result = CheckDifferentSWTestModes(
					configVersion,
					Form("Config version %d", configVersion),
					customClass,
					true
				);
		}
		return result;
	}
	
	if (! CheckDifferentModes(CheckPriorityGroupTestConfig, 0, "Priority group config", numOfEvents, customClass)) return false;
	if (! CheckDifferentModes(CheckSingleGroupTestConfig, 1, "Single group config", numOfEvents, customClass)) return false;
	if (! CheckDifferentModes(CheckPrescalarTestConfig, 2, "Prescalar config", numOfEvents, customClass)) return false;
	if (! CheckDifferentModes(CheckSymbolTestConfig, 3, "Symbol config", numOfEvents, customClass)) return false;
	if (! CheckDifferentModes(CheckComplexTestConfig, 4, "Complex config", numOfEvents, customClass)) return false;
	if (! CheckDifferentSWTestModes(5, "Software trigger config", customClass)) return false;
	
	// Cleanup all temporary files generated.
	gSystem->Exec("rm -f testOutputFile.root testInputFile*.root AliHLTGlobalTriggerImpl*");
	return true;
}


#ifndef __MAKECINT__

int main(int /*argc*/, const char** /*argv*/)
{
	bool resultOk = testGlobalTriggerComponent();
	if (not resultOk) return 1;
	return 0;
}

#endif // __MAKECINT__
