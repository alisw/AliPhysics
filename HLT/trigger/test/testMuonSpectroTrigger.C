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
 * @file   testMuonSpectroTrigger.C
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   11 Nov 2008
 *
 * This macro is used to test the AliHLTMuonSpectroTriggerComponent class.
 * A basic test is run with the AliHLTSystem framework to check that the dHLT
 * decision is interpreted correctly and the summary statistics generated properly.
 */

#if defined(__CINT__) && (! defined(__MAKECINT__))
#error This macro must be compiled. Try running as testMuonSpectroTrigger.C++, but remember to load the libAliHLTTrigger.so and libAliHLTMUON.so libraries first.
#endif

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TSystem.h"
#include "TClassTable.h"
#include "TFile.h"
#include "AliLog.h"
#include "AliHLTReadoutList.h"
#include "AliHLTTriggerDomain.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTMuonSpectroScalars.h"
#include "AliHLTMUONDataBlockWriter.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTSystem.h"
#include "AliHLTConfiguration.h"
#include "Riostream.h"
#include <stdio.h>
#endif

/**
 * Writes the indicated buffer to the given file name.
 * \returns true if the file was written correctly and false otherwise.
 */
bool WriteFile(const char* filename, void* buffer, unsigned int size)
{
	FILE* file = fopen(filename, "w");
	if (file == NULL)
	{
		cerr << "ERROR: Could not create file: " << filename << endl;
		return false;
	}
	fwrite(buffer, size, 1, file);
	fclose(file);
	return true;
}

/**
 * Generates some sample input data and writes it into 6 files named
 * testSinglesDecisionInputFile1.dat ... testSinglesDecisionInputFile3.dat and
 * testPairsDecisionInputFile1.dat ... testPairsDecisionInputFile3.dat
 */
void GenerateInputData()
{
	if (gClassTable->GetID("AliHLTMuonSpectroTriggerComponent") < 0)
	{
		gSystem->Load("libAliHLTUtil.so");
		gSystem->Load("libAliHLTTRD.so");
		gSystem->Load("libAliHLTTrigger.so");
	}
	
	// Allocate two 1 MByte buffers, this will be more than enough space.
	unsigned int bufferSize = 1024*1024;
	void* buffer1 = new char[bufferSize];
	void* buffer2 = new char[bufferSize];
	
	AliHLTMUONSinglesDecisionBlockWriter singlesBlock(buffer1, bufferSize);
	singlesBlock.InitCommonHeader();
	AliHLTMUONSinglesDecisionBlockStruct& singlesHeader = singlesBlock.BlockHeader();
	singlesHeader.fNlowPt = 2;
	singlesHeader.fNhighPt = 1;
	singlesBlock.SetNumberOfEntries(3);
	AliHLTMUONTrackDecisionStruct& track1 = singlesBlock[0];
	track1.fTrackId = 1;
	track1.fTriggerBits = AliHLTMUONUtils::PackTrackDecisionBits(false, false);
	track1.fPt = 0.9;
	AliHLTMUONTrackDecisionStruct& track2 = singlesBlock[1];
	track2.fTrackId = 2;
	track2.fTriggerBits = AliHLTMUONUtils::PackTrackDecisionBits(false, true);
	track2.fPt = 1.8;
	AliHLTMUONTrackDecisionStruct& track3 = singlesBlock[2];
	track3.fTrackId = 3;
	track3.fTriggerBits = AliHLTMUONUtils::PackTrackDecisionBits(true, true);
	track3.fPt = 3.1;
 
	AliHLTMUONPairsDecisionBlockWriter pairsBlock(buffer2, bufferSize);
	pairsBlock.InitCommonHeader();
	AliHLTMUONPairsDecisionBlockStruct& pairsHeader = pairsBlock.BlockHeader();
	pairsHeader.fNunlikeAnyPt = 1;
	pairsHeader.fNunlikeLowPt = 1;
	pairsHeader.fNunlikeHighPt = 1;
	pairsHeader.fNlikeAnyPt = 2;
	pairsHeader.fNlikeLowPt = 2;
	pairsHeader.fNlikeHighPt = 1;
	pairsHeader.fNmassAny = 1;
	pairsHeader.fNmassLow = 1;
	pairsHeader.fNmassHigh = 0;
	pairsBlock.SetNumberOfEntries(3);
	AliHLTMUONPairDecisionStruct& pair1 = pairsBlock[0];
	pair1.fTrackAId = 1;
	pair1.fTrackBId = 2;
	pair1.fTriggerBits = AliHLTMUONUtils::PackPairDecisionBits(false, false, false, 0, 1);
	pair1.fInvMass = 0.1;
	AliHLTMUONPairDecisionStruct& pair2 = pairsBlock[1];
	pair2.fTrackAId = 1;
	pair2.fTrackBId = 3;
	pair2.fTriggerBits = AliHLTMUONUtils::PackPairDecisionBits(false, true, true, 1, 1);
	pair2.fInvMass = 3.1;
	AliHLTMUONPairDecisionStruct& pair3 = pairsBlock[2];
	pair3.fTrackAId = 2;
	pair3.fTrackBId = 3;
	pair3.fTriggerBits = AliHLTMUONUtils::PackPairDecisionBits(true, true, false, 1, 2);
	pair3.fInvMass = 10.1;
	
	if (! WriteFile("testSinglesDecisionInputFile1.dat", buffer1, singlesBlock.BytesUsed())) return;
	if (! WriteFile("testPairsDecisionInputFile1.dat", buffer2, pairsBlock.BytesUsed())) return;
	
	singlesHeader.fNlowPt = 2;
	singlesHeader.fNhighPt = 1;
	singlesBlock.SetNumberOfEntries(2);
	track1.fTrackId = 4;
	track1.fTriggerBits = AliHLTMUONUtils::PackTrackDecisionBits(false, true);
	track1.fPt = 1.2;
	track2.fTrackId = 5;
	track2.fTriggerBits = AliHLTMUONUtils::PackTrackDecisionBits(true, true);
	track2.fPt = 2.3;
 
	pairsHeader.fNunlikeAnyPt = 1;
	pairsHeader.fNunlikeLowPt = 1;
	pairsHeader.fNunlikeHighPt = 1;
	pairsHeader.fNlikeAnyPt = 0;
	pairsHeader.fNlikeLowPt = 0;
	pairsHeader.fNlikeHighPt = 0;
	pairsHeader.fNmassAny = 1;
	pairsHeader.fNmassLow = 1;
	pairsHeader.fNmassHigh = 1;
	pairsBlock.SetNumberOfEntries(1);
	pair1.fTrackAId = 4;
	pair1.fTrackBId = 5;
	pair1.fTriggerBits = AliHLTMUONUtils::PackPairDecisionBits(true, true, true, 1, 2);
	pair1.fInvMass = 9.7;
	
	if (! WriteFile("testSinglesDecisionInputFile2.dat", buffer1, singlesBlock.BytesUsed())) return;
	if (! WriteFile("testPairsDecisionInputFile2.dat", buffer2, pairsBlock.BytesUsed())) return;
	
	singlesHeader.fNlowPt = 0;
	singlesHeader.fNhighPt = 0;
	singlesBlock.SetNumberOfEntries(1);
	track1.fTrackId = 6;
	track1.fTriggerBits = AliHLTMUONUtils::PackTrackDecisionBits(false, false);
	track1.fPt = 0.6;
 
	pairsHeader.fNunlikeAnyPt = 0;
	pairsHeader.fNunlikeLowPt = 0;
	pairsHeader.fNunlikeHighPt = 0;
	pairsHeader.fNlikeAnyPt = 0;
	pairsHeader.fNlikeLowPt = 0;
	pairsHeader.fNlikeHighPt = 0;
	pairsHeader.fNmassAny = 0;
	pairsHeader.fNmassLow = 0;
	pairsHeader.fNmassHigh = 0;
	pairsBlock.SetNumberOfEntries(0);
	
	if (! WriteFile("testSinglesDecisionInputFile3.dat", buffer1, singlesBlock.BytesUsed())) return;
	if (! WriteFile("testPairsDecisionInputFile3.dat", buffer2, pairsBlock.BytesUsed())) return;
}

/**
 * Runs a small test chain for the muon spectrometer trigger component.
 * \param debug  If true then full logging is turned on.
 * \param numOfEvents  The number of events to run the chain for.
 */
void RunTrigger(bool debug = false, int numOfEvents = 3)
{
	AliHLTSystem sys;
	sys.LoadComponentLibraries("libAliHLTUtil.so");
	sys.LoadComponentLibraries("libAliHLTTRD.so");
	sys.LoadComponentLibraries("libAliHLTMUON.so");
	sys.LoadComponentLibraries("libAliHLTTrigger.so");
	if (debug)
	{
		AliLog::SetGlobalLogLevel(AliLog::kMaxType);
		sys.SetGlobalLoggingLevel(kHLTLogAll);
	}
	
	TString cmdline = "-datatype DECIDSIN MUON ";
	for (int i = 1; i <= 3; i++)
	{
		if (i > 1) cmdline += " -nextevent";
		cmdline += Form(" -datafile testSinglesDecisionInputFile%d.dat", i);
	}
	AliHLTConfiguration pub1("pub1", "FilePublisher", NULL, cmdline.Data());
	
	cmdline = "-datatype DECIDPAR MUON ";
	for (int i = 1; i <= 3; i++)
	{
		if (i > 1) cmdline += " -nextevent";
		cmdline += Form(" -datafile testPairsDecisionInputFile%d.dat", i);
	}
	AliHLTConfiguration pub2("pub2", "FilePublisher", NULL, cmdline.Data());
	
	AliHLTConfiguration proc("proc", "MuonSpectroTrigger", "pub1 pub2", "-makestats -triggerdimuons");
	AliHLTConfiguration sink("sink", "ROOTFileWriter", "proc", "-datafile testMuonTriggerOutputFile.root -concatenate-events");
	
	sys.BuildTaskList("sink");
	sys.Run(numOfEvents);
}


/**
 * Checks that a particular decision and summary object are as expected and prints
 * error messages if there are any problems with the results.
 * \param eventNum  The number of the event being checked.
 * \param decision  The trigger decision being checked.
 * \param scalars  The scalars object being checked.
 * \param expectedResult  The expected global trigger result.
 * \param expectedDomain  The expected resulting trigger domain.
 * \param expectedDescription  The expected resulting trigger description.
 * \param expectedScalars  The expected resulting scalars.
 * \returns true if the decision and scalars are as expected.
 */
bool Check(
		int eventNum,
		AliHLTTriggerDecision* decision,
		AliHLTMuonSpectroScalars* scalars,
		bool expectedResult,
		AliHLTTriggerDomain expectedDomain,
		TString expectedDescription,
		AliHLTMuonSpectroScalars& expectedScalars
	)
{
	if (decision == NULL)
	{
		cerr << "ERROR: No decision found where expected for event "
		     << eventNum << "." << endl;
		return false;
	}
	if (decision->Result() != expectedResult)
	{
		cerr << "ERROR: The result does not match the expected value for event "
		     << eventNum << ". Got " << decision->Result() << " but expected "
		     << expectedResult << "." << endl;
		return false;
	}
	if (decision->TriggerDomain() != expectedDomain)
	{
		cerr << "ERROR: The domain does not match the expected value for event "
		     << eventNum << ". Got:" << endl;
		decision->TriggerDomain().Print();
		cerr << "but expected:" << endl;
		expectedDomain.Print();
		return false;
	}
	if (decision->Description() != expectedDescription)
	{
		cerr << "ERROR: The description does not match the expected value for event "
		     << eventNum << ". Got '" << decision->Description() << "' but expected '"
		     << expectedDescription << "'." << endl;
		return false;
	}
	if (*scalars != expectedScalars)
	{
		cerr << "ERROR: The scalars summary object does not match the expected one for event "
		     << eventNum << ". Expected the following scalar values:" << endl;
		expectedScalars.Print();
		return false;
	}
	return true;
}


/// Routine for checking the results of the chain run in RunTrigger.
bool CheckResults()
{
	AliHLTTriggerDecision* decision = NULL;
	AliHLTMuonSpectroScalars* scalars = NULL;
	bool result = false;
	
	AliHLTTriggerDomain domainMUON("*******:MUON");
	//domainMUON.Add(AliHLTReadoutList("MUONTRK"));
	
	AliHLTMuonSpectroScalars emptyScalars;
	
	AliHLTMuonSpectroScalars scalarsEvent1;
	scalarsEvent1.Add("NTracks", "", 3);
	scalarsEvent1.Add("NLowPt", "", 2);
	scalarsEvent1.Add("NHighPt", "", 1);
	scalarsEvent1.Add("MinPt", "", 0.9f);
	scalarsEvent1.Add("MaxPt", "", 3.1f);
	scalarsEvent1.Add("NLikeAny", "", 2);
	scalarsEvent1.Add("NLikeLow", "", 2);
	scalarsEvent1.Add("NLikeHigh", "", 1);
	scalarsEvent1.Add("NUnlikeAny", "", 1);
	scalarsEvent1.Add("NUnlikeLow", "", 1);
	scalarsEvent1.Add("NUnlikeHigh", "", 1);
	scalarsEvent1.Add("NLowMass", "", 1);
	scalarsEvent1.Add("NHighMass", "", 0);
	scalarsEvent1.Add("MinMass", "", 0.1f);
	scalarsEvent1.Add("MaxMass", "", 10.1f);
	
	AliHLTMuonSpectroScalars scalarsEvent2;
	scalarsEvent2.Add("NTracks", "", 2);
	scalarsEvent2.Add("NLowPt", "", 2);
	scalarsEvent2.Add("NHighPt", "", 1);
	scalarsEvent2.Add("MinPt", "", 1.2f);
	scalarsEvent2.Add("MaxPt", "", 2.3f);
	scalarsEvent2.Add("NLikeAny", "", 0);
	scalarsEvent2.Add("NLikeLow", "", 0);
	scalarsEvent2.Add("NLikeHigh", "", 0);
	scalarsEvent2.Add("NUnlikeAny", "", 1);
	scalarsEvent2.Add("NUnlikeLow", "", 1);
	scalarsEvent2.Add("NUnlikeHigh", "", 1);
	scalarsEvent2.Add("NLowMass", "", 1);
	scalarsEvent2.Add("NHighMass", "", 1);
	scalarsEvent2.Add("MinMass", "", 9.7f);
	scalarsEvent2.Add("MaxMass", "", 9.7f);
	
	AliHLTMuonSpectroScalars scalarsEvent3;
	scalarsEvent3.Add("NTracks", "", 1);
	scalarsEvent3.Add("NLowPt", "", 0);
	scalarsEvent3.Add("NHighPt", "", 0);
	scalarsEvent3.Add("MinPt", "", 0.6f);
	scalarsEvent3.Add("MaxPt", "", 0.6f);
	scalarsEvent3.Add("NLikeAny", "", 0);
	scalarsEvent3.Add("NLikeLow", "", 0);
	scalarsEvent3.Add("NLikeHigh", "", 0);
	scalarsEvent3.Add("NUnlikeAny", "", 0);
	scalarsEvent3.Add("NUnlikeLow", "", 0);
	scalarsEvent3.Add("NUnlikeHigh", "", 0);
	scalarsEvent3.Add("NLowMass", "", 0);
	scalarsEvent3.Add("NHighMass", "", 0);
	scalarsEvent3.Add("MinMass", "", -1);
	scalarsEvent3.Add("MaxMass", "", -1);

	TFile* file = new TFile("testMuonTriggerOutputFile.root", "READ");
	
	// First event is the Start-Of-Run
	decision = dynamic_cast<AliHLTTriggerDecision*>(file->Get("MuonSpectroTrigger;1"));
	scalars = dynamic_cast<AliHLTMuonSpectroScalars*>(file->Get("AliHLTMuonSpectroScalars;1"));
	result = Check(0, decision, scalars, false, AliHLTTriggerDomain(), "Not triggered", emptyScalars);
	if (! result) goto cleanup;
	// Now we have the 3 data events:
	decision = dynamic_cast<AliHLTTriggerDecision*>(file->Get("MuonSpectroTrigger;2"));
	scalars = dynamic_cast<AliHLTMuonSpectroScalars*>(file->Get("AliHLTMuonSpectroScalars;2"));
	result = Check(1, decision, scalars, true, domainMUON, "Dimuon in muon spectrometer", scalarsEvent1);
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTTriggerDecision*>(file->Get("MuonSpectroTrigger;3"));
	scalars = dynamic_cast<AliHLTMuonSpectroScalars*>(file->Get("AliHLTMuonSpectroScalars;3"));
	result = Check(2, decision, scalars, true, domainMUON, "Dimuon in muon spectrometer", scalarsEvent2);
	if (! result) goto cleanup;
	decision = dynamic_cast<AliHLTTriggerDecision*>(file->Get("MuonSpectroTrigger;4"));
	scalars = dynamic_cast<AliHLTMuonSpectroScalars*>(file->Get("AliHLTMuonSpectroScalars;4"));
	result = Check(3, decision, scalars, false, AliHLTTriggerDomain(), "Not triggered", scalarsEvent3);
	if (! result) goto cleanup;
	// and finally the End-Of-Run event.
	decision = dynamic_cast<AliHLTTriggerDecision*>(file->Get("MuonSpectroTrigger;5"));
	scalars = dynamic_cast<AliHLTMuonSpectroScalars*>(file->Get("AliHLTMuonSpectroScalars;5"));
	result = Check(4, decision, scalars, false, AliHLTTriggerDomain(), "Not triggered", emptyScalars);
	if (! result) goto cleanup;
	
	delete file;
	return true;
	
cleanup:
	if (decision != NULL)
	{
		cout << "========== Dumping received decision result ========== " << endl;
		decision->Print();
		cout << "=========== Dumping received scalars result ========== " << endl;
		scalars->Print();
	}
	delete file;
	return false;
}


/**
 * Runs the unit test for the AliHLTMuonSpectroTriggerComponent class.
 * \param debug  If specified then the HLT chain is run with full logging enabled.
 * \returns true if the class passed the test and false otherwise.
 */
bool testMuonSpectroTrigger(bool debug = false)
{
	GenerateInputData();
	RunTrigger(debug);
	if (! CheckResults()) return false;
	
	// Cleanup all temporary files generated.
	gSystem->Exec("rm -f testMuonTriggerOutputFile.root test*DecisionInputFile*.dat");
	return true;
}


#ifndef __MAKECINT__

int main(int /*argc*/, const char** /*argv*/)
{
	bool resultOk = testMuonSpectroTrigger();
	if (not resultOk) return 1;
	return 0;
}

#endif // __MAKECINT__
