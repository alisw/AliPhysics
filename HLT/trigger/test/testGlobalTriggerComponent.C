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
 */
 
void testGlobalTriggerComponent()
{
	gSystem->Load("libAliHLTTrigger.so");

	TFile* file = new TFile("testInputFile.root", "RECREATE");
	
	AliHLTReadoutList readoutList1("TPC");
	AliHLTTriggerDomain triggerDomain1;
	AliHLTTriggerDecision decision1(true, "Trigger1", readoutList1, triggerDomain1);
	
	AliHLTReadoutList readoutList2;
	AliHLTTriggerDomain triggerDomain2;
	triggerDomain2.Add("aaaa", "bbbb");
	AliHLTTriggerDecision decision2(true, "Trigger2", readoutList2, triggerDomain2);
	
	decision1.Write("Trigger1");
	decision2.Write("Trigger2");
	delete file;

	AliHLTSystem sys;
	sys.LoadComponentLibraries("libAliHLTUtil.so");
	sys.LoadComponentLibraries("libAliHLTTrigger.so");
	
	AliHLTConfiguration pub("pub", "ROOTFilePublisher", NULL, " -datatype ROOTTOBJ 'HLT ' -datafile testInputFile.root");
	AliHLTConfiguration proc("proc", "HLTGlobalTrigger", "pub", "");
	AliHLTConfiguration sink("sink", "ROOTFileWriter", "proc", "-datafile testOutput.root");
	
	sys.BuildTaskList("sink");
	sys.Run(1);
}

