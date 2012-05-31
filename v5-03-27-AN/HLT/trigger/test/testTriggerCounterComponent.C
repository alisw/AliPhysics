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
 * @file   testTriggerCounterComponent.C
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   5 Nov 2010
 *
 * This macro is used to test the AliHLTTriggerCounterComponent class.
 * Tests are run with component within the AliHLTSystem framework to check that
 * the trigger counters are generated correctly by the component.
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TSystem.h"
#include "TFile.h"
#include "AliLog.h"
#include "AliHLTTriggerCounters.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTGlobalTriggerDecision.h"
#include "AliHLTSystem.h"
#include "AliHLTConfiguration.h"
#include "AliHLTCTPData.h"
#include "Riostream.h"
#endif

/**
 * Generates some sample input data and writes it into 9 files named
 * testInputFileTriggerCounter1.root ... testInputFileTriggerCounter9.root
 */
void GenerateInputData()
{
	AliHLTTriggerDecision td1(true, "TRIGGER-A", AliHLTTriggerDomain("TRACKS:TPC"));
	AliHLTTriggerDecision td1f(false, "TRIGGER-A", AliHLTTriggerDomain());
	AliHLTTriggerDecision td2(true, "TRIGGER-B", AliHLTTriggerDomain("CLUSTERS:TPC"));
	AliHLTTriggerDecision td2f(false, "TRIGGER-B", AliHLTTriggerDomain());
	AliHLTTriggerDecision td3(true, "TRIGGER-C", AliHLTTriggerDomain("*******:***"));
	AliHLTTriggerDecision td3f(false, "TRIGGER-C", AliHLTTriggerDomain());
	AliHLTGlobalTriggerDecision gd1(true, AliHLTTriggerDomain("TRACKS:TPC"), "TRIGGER-A");
	AliHLTGlobalTriggerDecision gd2(true, AliHLTTriggerDomain("CLUSTERS:TPC"), "TRIGGER-B");
	AliHLTGlobalTriggerDecision gd3(true, AliHLTTriggerDomain("*******:***"), "TRIGGER-B,TRIGGER-C");
	
	AliHLTCTPData ctp("CTP_TRIGGER_CLASS=00:CINTA:17,01:CINTB:00-01-02-03-04-05-06-07-08-09-10-11-12-13-14-15-16-17,02:CMUB:00-10-11-17");
	gd1.AddInputObjectRef(&ctp);
	gd2.AddInputObjectRef(&ctp);
	gd3.AddInputObjectRef(&ctp);
	
	TFile* file = new TFile("testInputFileTriggerCounter1.root", "RECREATE");
	ctp.Increment(int(0));
	gd1.Write("HLTGlobalTrigger");
	td1.Write("Trigger");
	delete file;
	
	file = new TFile("testInputFileTriggerCounter2.root", "RECREATE");
	ctp.Increment(int(1));
	gd2.Write("HLTGlobalTrigger");
	td1f.Write("Trigger");
	td2.Write("Trigger");
	td3f.Write("Trigger");
	delete file;
	
	file = new TFile("testInputFileTriggerCounter3.root", "RECREATE");
	ctp.Increment(int(2));
	gd3.Write("HLTGlobalTrigger");
	td1f.Write("Trigger");
	td2.Write("Trigger");
	td3.Write("Trigger");
	delete file;
	
	gd1.AddTriggerInput(td1);
	gd2.AddTriggerInput(td1f);
	gd2.AddTriggerInput(td2);
	gd2.AddTriggerInput(td3f);
	gd3.AddTriggerInput(td1f);
	gd3.AddTriggerInput(td2);
	gd3.AddTriggerInput(td3);
	
	file = new TFile("testInputFileTriggerCounter4.root", "RECREATE");
	ctp.Increment(int(0));
	gd1.Write("HLTGlobalTrigger");
	td1.Write("Trigger");
	delete file;
	
	file = new TFile("testInputFileTriggerCounter5.root", "RECREATE");
	ctp.Increment(int(1));
	gd2.Write("HLTGlobalTrigger");
	td1f.Write("Trigger");
	td2.Write("Trigger");
	delete file;
	
	file = new TFile("testInputFileTriggerCounter6.root", "RECREATE");
	ctp.Increment(int(2));
	gd3.Write("HLTGlobalTrigger");
	td1f.Write("Trigger");
	td2.Write("Trigger");
	td3.Write("Trigger");
	delete file;
	
	file = new TFile("testInputFileTriggerCounter7.root", "RECREATE");
	ctp.Increment(int(0));
	gd1.Write("HLTGlobalTrigger");
	delete file;
	
	file = new TFile("testInputFileTriggerCounter8.root", "RECREATE");
	ctp.Increment(int(1));
	gd2.Write("HLTGlobalTrigger");
	delete file;
	
	file = new TFile("testInputFileTriggerCounter9.root", "RECREATE");
	ctp.Increment(int(2));
	gd3.Write("HLTGlobalTrigger");
	delete file;
}

/**
 * Runs a small chain with the HLTTriggerCounter component to generate output that
 * can be checked.
 * \param debug  If true then full debug logging is enabled.
 * \param numOfEvents  The number of events to run the chain for.
 */
void RunTriggerCounter(bool debug = false, int numOfEvents = 9)
{
	AliHLTSystem sys;
	sys.ScanOptions("ECS=CTP_TRIGGER_CLASS=00:CINTA:17,01:CINTB:00-01-02-03-04-05-06-07-08-09-10-11-12-13-14-15-16-17,02:CMUB:00-10-11-17");
	sys.LoadComponentLibraries("libAliHLTUtil.so");
	sys.LoadComponentLibraries("libAliHLTTRD.so");
	sys.LoadComponentLibraries("libAliHLTMUON.so");
	sys.LoadComponentLibraries("libAliHLTTrigger.so");
	if (debug)
	{
		AliLog::SetGlobalLogLevel(AliLog::kMaxType);
		sys.SetGlobalLoggingLevel(kHLTLogAll);
	}
	
	TString cmdline = "-datatype GLOBTRIG 'HLT ' -objectname HLTGlobalTrigger ";
	for (int i = 1; i <= 9; i++)
	{
		if (i > 1) cmdline += " -nextevent";
		cmdline += Form(" -datafile testInputFileTriggerCounter%d.root", i);
	}
	AliHLTConfiguration pub1("pub1", "ROOTFilePublisher", NULL, cmdline.Data());
	cmdline = "-datatype TRIG_DEC 'HLT ' -objectname Trigger ";
	for (int i = 1; i <= 9; i++)
	{
		if (i > 1) cmdline += " -nextevent";
		cmdline += Form(" -datafile testInputFileTriggerCounter%d.root", i);
	}
	AliHLTConfiguration pub2("pub2", "ROOTFilePublisher", NULL, cmdline.Data());
	AliHLTConfiguration proc("proc", "HLTTriggerCounter", "pub1 pub2", "-skipcdb -config InitialCounterConfig.C");
	AliHLTConfiguration sink("sink", "ROOTFileWriter", "proc", "-datafile testOutputFileTriggerCounter.root -concatenate-events");
	sys.BuildTaskList("sink");
	AliHLTUInt64_t trigMask[9] = {0x1, 0x2, 0x4, 0x1, 0x2, 0x4, 0x1, 0x2, 0x4};
	for (int i = 0; i < 9 && i < numOfEvents; i++)
	{
		sys.Run(
			1,  // number of events
			(i == 8 ? 1 : 0),  // Stop chain at end of run.
			trigMask[i], // Active CTP trigger mask.
			0,   // Time stamp.
			0    // Event type.
		);
	}
}

/**
 * Checks that the input and output counters have the expected values or not.
 * If they do not an appropriate error message is printed.
 * \param eventNum  The number of the event being checked.
 * \param inputCounters  The input counters to check.
 * \param outputCounters  The output counters to check.
 * \param expectedInputCounters  The expected input counter values to check against.
 * \param expectedOutputCounters  The expected output counter values to check against.
 * \param expectedDescription  The expected resulting trigger description.
 * \returns true if the counters are correct and false otherwise.
 */
bool CheckCounters(
		int eventNum,
		AliHLTTriggerCounters* inputCounters,
		AliHLTTriggerCounters* outputCounters,
		AliHLTTriggerCounters& expectedInputCounters,
		AliHLTTriggerCounters& expectedOutputCounters
	)
{
	if (inputCounters == NULL)
	{
		cerr << "ERROR: No trigger input counters were found for event "
		     << eventNum << "." << endl;
		return false;
	}
	if (outputCounters == NULL)
	{
		cerr << "ERROR: No trigger output counters were found for event "
		     << eventNum << "." << endl;
		return false;
	}
	// Copy over the Rates because these might not be identical since they are a
	// run time parameter. The counter value and name must be identical however.
	for (UInt_t i = 0; i < inputCounters->NumberOfCounters() && i < expectedInputCounters.NumberOfCounters(); ++i)
	{
		expectedInputCounters.GetCounterN(i).Rate( inputCounters->GetCounterN(i).Rate() );
	}
	for (UInt_t i = 0; i < outputCounters->NumberOfCounters() && i < expectedOutputCounters.NumberOfCounters(); ++i)
	{
		expectedOutputCounters.GetCounterN(i).Rate( outputCounters->GetCounterN(i).Rate() );
	}
	bool inputCountersAsExpected = (inputCounters->NumberOfCounters() == expectedInputCounters.NumberOfCounters());
	bool outputCountersAsExpected = (outputCounters->NumberOfCounters() == expectedOutputCounters.NumberOfCounters());
	for (UInt_t i = 0; i < inputCounters->NumberOfCounters() && i < expectedInputCounters.NumberOfCounters(); ++i)
	{
		if (strcmp(expectedInputCounters.GetCounterN(i).Name(), inputCounters->GetCounterN(i).Name()) != 0)
		{
			inputCountersAsExpected = false;
			break;
		}
		if (expectedInputCounters.GetCounterN(i).Counter() != inputCounters->GetCounterN(i).Counter())
		{
			inputCountersAsExpected = false;
			break;
		}
	}
	for (UInt_t i = 0; i < outputCounters->NumberOfCounters() && i < expectedOutputCounters.NumberOfCounters(); ++i)
	{
		if (strcmp(expectedOutputCounters.GetCounterN(i).Name(), outputCounters->GetCounterN(i).Name()) != 0)
		{
			outputCountersAsExpected = false;
			break;
		}
		if (expectedOutputCounters.GetCounterN(i).Counter() != outputCounters->GetCounterN(i).Counter())
		{
			outputCountersAsExpected = false;
			break;
		}
	}
	if (! inputCountersAsExpected)
	{
		cerr << "ERROR: The trigger input counters do not match the expected set of counters for event "
		     << eventNum << "." << endl;
		cout << "Found the following input counters: " << endl;  // using cout since Print() goes there.
		inputCounters->Print();
		cout << "But expected the following input counters: " << endl;
		expectedInputCounters.Print();
		return false;
	}
	if (! outputCountersAsExpected)
	{
		cerr << "ERROR: The trigger output counters do not match the expected set of counters for event "
		     << eventNum << "." << endl;
		cout << "Found the following output counters: " << endl;  // using cout since Print() goes there.
		outputCounters->Print();
		cout << "But expected the following output counters: " << endl;
		expectedOutputCounters.Print();
		return false;
	}
	return true;
}


/// Routine for checking the results file generated by the RunTriggerCounter routine.
bool CheckTriggerCounterOutput()
{
	AliHLTTriggerCounters expectedInputCounters[11];
	AliHLTTriggerCounters expectedOutputCounters[11];
	
	// Start-of-run event (should contain no values)
	expectedInputCounters[0].Add("CINTA", "", 0, 0);
	expectedInputCounters[0].Add("CINTB", "", 0, 0);
	expectedInputCounters[0].Add("CMUB", "", 0, 0);
	expectedOutputCounters[0].Add("TRIGGER-A", "", 0, 0);
	
	// Data events
	expectedInputCounters[1].Add("CINTA", "", 0, 1);
	expectedInputCounters[1].Add("CINTB", "", 0, 0);
	expectedInputCounters[1].Add("CMUB", "", 0, 0);
	expectedInputCounters[1].Add("TRIGGER-A", "", 0, 1);
	expectedOutputCounters[1].Add("TRIGGER-A", "", 0, 1);
	
	expectedInputCounters[2].Add("CINTA", "", 0, 1);
	expectedInputCounters[2].Add("CINTB", "", 0, 1);
	expectedInputCounters[2].Add("CMUB", "", 0, 0);
	expectedInputCounters[2].Add("TRIGGER-A", "", 0, 1);
	expectedInputCounters[2].Add("TRIGGER-B", "", 0, 1);
	expectedOutputCounters[2].Add("TRIGGER-A", "", 0, 1);
	expectedOutputCounters[2].Add("TRIGGER-B", "", 0, 1);
	
	expectedInputCounters[3].Add("CINTA", "", 0, 1);
	expectedInputCounters[3].Add("CINTB", "", 0, 1);
	expectedInputCounters[3].Add("CMUB", "", 0, 1);
	expectedInputCounters[3].Add("TRIGGER-A", "", 0, 1);
	expectedInputCounters[3].Add("TRIGGER-B", "", 0, 2);
	expectedInputCounters[3].Add("TRIGGER-C", "", 0, 1);
	expectedOutputCounters[3].Add("TRIGGER-A", "", 0, 1);
	expectedOutputCounters[3].Add("TRIGGER-B", "", 0, 2);
	expectedOutputCounters[3].Add("TRIGGER-C", "", 0, 1);
	
	expectedInputCounters[4].Add("CINTA", "", 0, 2);
	expectedInputCounters[4].Add("CINTB", "", 0, 1);
	expectedInputCounters[4].Add("CMUB", "", 0, 1);
	expectedInputCounters[4].Add("TRIGGER-A", "", 0, 2);
	expectedInputCounters[4].Add("TRIGGER-B", "", 0, 2);
	expectedInputCounters[4].Add("TRIGGER-C", "", 0, 1);
	expectedOutputCounters[4].Add("TRIGGER-A", "", 0, 2);
	expectedOutputCounters[4].Add("TRIGGER-B", "", 0, 2);
	expectedOutputCounters[4].Add("TRIGGER-C", "", 0, 1);
	
	expectedInputCounters[5].Add("CINTA", "", 0, 2);
	expectedInputCounters[5].Add("CINTB", "", 0, 2);
	expectedInputCounters[5].Add("CMUB", "", 0, 1);
	expectedInputCounters[5].Add("TRIGGER-A", "", 0, 2);
	expectedInputCounters[5].Add("TRIGGER-B", "", 0, 3);
	expectedInputCounters[5].Add("TRIGGER-C", "", 0, 1);
	expectedOutputCounters[5].Add("TRIGGER-A", "", 0, 2);
	expectedOutputCounters[5].Add("TRIGGER-B", "", 0, 3);
	expectedOutputCounters[5].Add("TRIGGER-C", "", 0, 1);
	
	expectedInputCounters[6].Add("CINTA", "", 0, 2);
	expectedInputCounters[6].Add("CINTB", "", 0, 2);
	expectedInputCounters[6].Add("CMUB", "", 0, 2);
	expectedInputCounters[6].Add("TRIGGER-A", "", 0, 2);
	expectedInputCounters[6].Add("TRIGGER-B", "", 0, 4);
	expectedInputCounters[6].Add("TRIGGER-C", "", 0, 2);
	expectedOutputCounters[6].Add("TRIGGER-A", "", 0, 2);
	expectedOutputCounters[6].Add("TRIGGER-B", "", 0, 4);
	expectedOutputCounters[6].Add("TRIGGER-C", "", 0, 2);
	
	expectedInputCounters[7].Add("CINTA", "", 0, 3);
	expectedInputCounters[7].Add("CINTB", "", 0, 2);
	expectedInputCounters[7].Add("CMUB", "", 0, 2);
	expectedInputCounters[7].Add("TRIGGER-A", "", 0, 3);
	expectedInputCounters[7].Add("TRIGGER-B", "", 0, 4);
	expectedInputCounters[7].Add("TRIGGER-C", "", 0, 2);
	expectedOutputCounters[7].Add("TRIGGER-A", "", 0, 3);
	expectedOutputCounters[7].Add("TRIGGER-B", "", 0, 4);
	expectedOutputCounters[7].Add("TRIGGER-C", "", 0, 2);
	
	expectedInputCounters[8].Add("CINTA", "", 0, 3);
	expectedInputCounters[8].Add("CINTB", "", 0, 3);
	expectedInputCounters[8].Add("CMUB", "", 0, 2);
	expectedInputCounters[8].Add("TRIGGER-A", "", 0, 3);
	expectedInputCounters[8].Add("TRIGGER-B", "", 0, 5);
	expectedInputCounters[8].Add("TRIGGER-C", "", 0, 2);
	expectedOutputCounters[8].Add("TRIGGER-A", "", 0, 3);
	expectedOutputCounters[8].Add("TRIGGER-B", "", 0, 5);
	expectedOutputCounters[8].Add("TRIGGER-C", "", 0, 2);
	
	expectedInputCounters[9].Add("CINTA", "", 0, 3);
	expectedInputCounters[9].Add("CINTB", "", 0, 3);
	expectedInputCounters[9].Add("CMUB", "", 0, 3);
	expectedInputCounters[9].Add("TRIGGER-A", "", 0, 3);
	expectedInputCounters[9].Add("TRIGGER-B", "", 0, 6);
	expectedInputCounters[9].Add("TRIGGER-C", "", 0, 3);
	expectedOutputCounters[9].Add("TRIGGER-A", "", 0, 3);
	expectedOutputCounters[9].Add("TRIGGER-B", "", 0, 6);
	expectedOutputCounters[9].Add("TRIGGER-C", "", 0, 3);
	
	// End-of-run event
	expectedInputCounters[10].Add("CINTA", "", 0, 3);
	expectedInputCounters[10].Add("CINTB", "", 0, 3);
	expectedInputCounters[10].Add("CMUB", "", 0, 3);
	expectedInputCounters[10].Add("TRIGGER-A", "", 0, 3);
	expectedInputCounters[10].Add("TRIGGER-B", "", 0, 6);
	expectedInputCounters[10].Add("TRIGGER-C", "", 0, 3);
	expectedOutputCounters[10].Add("TRIGGER-A", "", 0, 3);
	expectedOutputCounters[10].Add("TRIGGER-B", "", 0, 6);
	expectedOutputCounters[10].Add("TRIGGER-C", "", 0, 3);
	
	TFile* file = new TFile("testOutputFileTriggerCounter.root", "READ");
	for (int i = 0; i < 11; ++i)
	{
		AliHLTTriggerCounters* inputCounters = dynamic_cast<AliHLTTriggerCounters*>( file->Get(Form("AliHLTTriggerCounters;%d", i*2+1)) );
		AliHLTTriggerCounters* outputCounters = dynamic_cast<AliHLTTriggerCounters*>( file->Get(Form("AliHLTTriggerCounters;%d", i*2+2)) );
		if (! CheckCounters(i-1, inputCounters, outputCounters, expectedInputCounters[i], expectedOutputCounters[i]))
		{
			delete file;
			return false;
		}
	}
	delete file;
	return true;
}

/**
 * Runs several tests for the AliHLTTriggerCounterComponent class.
 * We specifically test if the generated counters are correctly generated.
 * \param debug  If true then full debug logging is enabled.
 * \param numOfEvents  The number of events to run the chain over.
 * \returns true if all the tests succeeded and false otherwise.
 */
bool testTriggerCounterComponent(bool debug = false, int numOfEvents = 9)
{
	GenerateInputData();
	RunTriggerCounter(debug, numOfEvents);
	if (! CheckTriggerCounterOutput()) return false;
	// Cleanup all temporary files generated.
	gSystem->Exec("rm -f testOutputFileTriggerCounter.root testInputFileTriggerCounter*.root");
	return true;
}


#ifndef __MAKECINT__

int main(int /*argc*/, const char** /*argv*/)
{
	bool resultOk = testTriggerCounterComponent();
	if (not resultOk) return 1;
	return 0;
}

#endif // __MAKECINT__
