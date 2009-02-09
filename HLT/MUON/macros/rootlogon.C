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

/* $Id: $ */

/**
 * \ingroup macros
 * \file rootlogon.C
 * \brief This macro is run when AliRoot starts and sets up the dHLT specific environment.
 *
 * It loads the required dimuon HLT specific libraries and include paths needed
 * to compile and run dHLT specific macros easily.
 *
 * \author Artur Szostak <artursz@iafrica.com>
 */

{
	cout << "Initialising AliRoot environment for dHLT..." << endl;
	
	// Check that we are actually loading aliroot or something compatible.
	if (gClassTable->GetID("AliHLTSystem") < 0)
	{
		cerr << "ERROR: you must run in AliRoot to run dHLT macros." << endl;
		gApplication->Terminate();
		return;
	}
	
	TString includePath = "-I${ALICE_ROOT}/include ";
	includePath        += "-I${ALICE_ROOT}/RAW ";
	includePath        += "-I${ALICE_ROOT}/MUON ";
	includePath        += "-I${ALICE_ROOT}/MUON/mapping ";
	includePath        += "-I${ALICE_ROOT}/HLT/BASE ";
	includePath        += "-I${ALICE_ROOT}/HLT/BASE/HOMER ";
	includePath        += "-I${ALICE_ROOT}/HLT/MUON ";
	includePath        += "-I${ALICE_ROOT}/HLT/MUON/macros ";
	includePath        += "-I${ALICE_ROOT}/HLT/MUON/OfflineInterface ";
	includePath        += "-I${ALICE_ROOT}/HLT/MUON/OnlineAnalysis ";
	gSystem->SetIncludePath(includePath.Data());
	
	TString macroPath = gROOT->GetMacroPath();
	macroPath        += "${ALICE_ROOT}/HLT/MUON/macros:";
	gROOT->SetMacroPath(macroPath);
	
	gSystem->Load("libAliHLTMUON.so");
	gSystem->Load("libAliHLTHOMER.so");
	
	// Setup the CDB default storage and run number if nothing was set.
	AliCDBManager* cdbManager = AliCDBManager::Instance();
	if (cdbManager == NULL)
	{
		cerr << "ERROR: Global CDB manager object does not exist." << endl;
		return;
	}
	if (cdbManager->GetDefaultStorage() == NULL)
	{
		cdbManager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
	}
	if (cdbManager->GetRun() == -1)
	{
		cdbManager->SetRun(0);
	}
}
