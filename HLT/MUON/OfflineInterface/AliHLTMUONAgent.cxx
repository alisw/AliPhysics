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

// $Id$

///
/// @file   AliHLTMUONAgent.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   28 May 2007
/// @brief  Implementation of the AliHLTMUONAgent class.
///

#include "AliHLTMUONAgent.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONRecHitsSource.h"
#include "AliHLTMUONTriggerRecordsSource.h"
#include "AliHLTMUONDigitPublisherComponent.h"
#include "AliHLTMUONRootifierComponent.h"
#include "AliHLTMUONHitReconstructorComponent.h"
#include "AliHLTMUONTriggerReconstructorComponent.h"
#include "AliHLTMUONMansoTrackerFSMComponent.h"
#include "AliHLTMUONFullTrackerComponent.h"
#include "AliHLTMUONDecisionComponent.h"
#include "AliHLTMUONESDMaker.h"
#include "AliHLTMUONEmptyEventFilterComponent.h"
#include "AliHLTMUONDataCheckerComponent.h"
#include "AliHLTMUONClusterFinderComponent.h"
#include "AliHLTMUONRawDataHistoComponent.h"
#include "AliHLTMUONClusterHistoComponent.h"
#include "AliHLTOUTHandlerChain.h"
#include "AliHLTOUTHandlerIgnore.h"
#include "AliRawReader.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliMUON.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TString.h"

// The single global instance of the dimuon HLT agent.
AliHLTMUONAgent AliHLTMUONAgent::fgkInstance;

AliHLTOUTHandlerChain AliHLTMUONAgent::fgkESDMakerChain("libAliHLTMUON.so chains=dHLT-make-esd");
AliHLTOUTHandlerChain AliHLTMUONAgent::fgkRootifyDumpChain("libAliHLTMUON.so chains=dHLT-rootify-and-dump");
AliHLTOUTHandlerIgnore AliHLTMUONAgent::fgkDataIgnoreHandler;
Int_t AliHLTMUONAgent::fgMuonModuleLoaded = 0;
bool AliHLTMUONAgent::fgRunRootifyChain = false;


ClassImp(AliHLTMUONAgent);


bool AliHLTMUONAgent::IsMuonModuleLoaded()
{
	/// Checks to see if the MUON module is loaded or not.

	// If the check was already done then use the cached value.
	if (fgMuonModuleLoaded > 0) return true;
	if (fgMuonModuleLoaded < 0) return false;
	
	if (gAlice != NULL)
	{
		// Search for a module in gAlice deriving from AliMUON.
		TIter next(gAlice->Modules());
		TObject* mod = NULL;
		while ((mod = next()) != NULL)
		{
			if (mod->InheritsFrom(AliMUON::Class()))
			{
				fgMuonModuleLoaded = 1;
				return true;
			}
		}
	}
	
	fgMuonModuleLoaded = -1;
	return false;
}


AliHLTMUONAgent::AliHLTMUONAgent() : AliHLTModuleAgent("MUON")
{
	///
	/// Default constructor.
	///
}

AliHLTMUONAgent::~AliHLTMUONAgent()
{
	///
	/// Default destructor.
	///
}

const char* AliHLTMUONAgent::GetReconstructionChains(AliRawReader* rawReader,
						     AliRunLoader* runloader
	) const
{
	///
	/// Inherited from AliHLTModuleAgent.
	/// Returns the top processing chain configurations for local event
	/// reconstruction.
	/// @param rawReader  [in] AliRoot rawreader instance.
	/// @param runloader  [in] AliRoot runloader
	/// @return string containing the top configurations separated by blanks.
	///
	/// If rawReader is not NULL then the standard dHLT chain is run taking
	/// data from raw DDL data. Otherwise runloader is checked and if it is
	/// not NULL then a dHLT chain is run with input data from digits.
	
	if (rawReader != NULL)
	{
		// Check if there is any data from the tracker and trigger.
		bool dataFromTracker = false;
		bool dataFromTrigger = false;
		rawReader->Select("MUONTRK", 0, 19);
		for (Int_t i = 0; i < 20; i++)
		{
			if (rawReader->ReadHeader()) dataFromTracker = true;
		}
		rawReader->Select("MUONTRG", 0, 1);
		for (Int_t i = 0; i < 2; i++)
		{
			if (rawReader->ReadHeader()) dataFromTrigger = true;
		}
		rawReader->Reset();
		
		// If raw data was found for our detector then select the
		// appropriate chain.
		if (dataFromTracker and dataFromTrigger)
			return "dHLT-sim-fromRaw";
	}
	
	if (runloader != NULL)
	{
		// IsMuonModuleLoaded() is used to check if the muon module was loaded
		// If there is no AliMUON module in the simulation then do not run the
		// MUON HLT chain.
		if (IsMuonModuleLoaded() and runloader->GetLoader("MUONLoader") != NULL)
			return "dHLT-sim";
	}
	
	return "";
}

const char* AliHLTMUONAgent::GetRequiredComponentLibraries() const
{
	///
	/// Inherited from AliHLTModuleAgent.
	/// Returns a list of libraries which the configurations registered by
	/// this module agent depend on.
	/// @return list of component libraries as a blank-separated string.
	///
	
	// List of libraries that we depend on.
	static const char* libs[] =
	{
		"libCore.so",
		"libCint.so",
		"libGraf.so",
		"libRIO.so",
		"libNet.so",
		"libHist.so",
		"libMatrix.so",
		"libMathCore.so",
		"libMinuit.so",
		"libTree.so",
		"libGraf3d.so",
		"libGpad.so",
		"libPhysics.so",
		"libGui.so",
		"libProofPlayer.so",
		"libProof.so",
		"libThread.so",
		"libGeom.so",
		"libEG.so",
		"libTreePlayer.so",
		"libXMLIO.so",
		"libVMC.so",
		"libESD.so",
		"libCDB.so",
		"libSTEERBase.so",
		"libSTEER.so",
		"libGui.so",
		"libRAWDatasim.so",
		"libRAWDatarec.so",
		"libRAWDatabase.so",
		"libMUONcore.so",
		"libMUONraw.so",
		"libMUONbase.so",
		"libMUONgeometry.so",
		"libMUONmapping.so",
		"libMUONcalib.so",
		"libMUONsim.so",
		"libMUONtrigger.so",
		"libMUONevaluation.so",
		"libMUONrec.so",
		"libANALYSIS.so",
		"libANALYSISalice.so",
		"libHLTbase.so",
		"libAliHLTUtil.so",
		NULL
	};
	
	// First check if the library is not already loaded. If it is then we have
	// no reason to declare it as needed, so that we do not load it again.
	static TString result;
	for (const char** lib = libs; *lib != NULL; lib++)
	{
		const char* list = gSystem->GetLibraries(*lib, "", kFALSE);
		if (list == NULL) continue;
		// Check if not found, i.e. result was an empty string.
		// If so then add the library to the required list.
		if (list[0] == '\0')
		{
			result += *lib;
			result += " ";
		}
	}
	return result.Data();
}


int AliHLTMUONAgent::CreateConfigurations(
		AliHLTConfigurationHandler* handler,
		AliRawReader* rawReader,
		AliRunLoader* runloader
	) const
{
	/// Register all processing configurations belonging to the dimuon HLT
	/// library with the AliHLTConfigurationHandler.
	/// @param handler      the configuration handler
	/// @param rawReader  [in] AliRoot rawreader instance.
	/// @param runloader    AliRoot runloader
	/// @return Zero on success and error code if failed.
	///
	/// Chains available:
	/// dHLT-sim          - standard dHLT simulation chain.
	/// dHLT-sim-fromRaw  - standard dHLT chain taking raw DDL data as input.
	/// dHLT-sim-fromMC   - dHLT chain taking Monte Carlo hit data as input.
	///                     So hit reconstruction is not checked, just the tracker.
	
	if (handler == NULL) return 0;
	
	const char* trackerId = AliHLTMUONConstants::MansoTrackerFSMId();
	const char* fullTrackerId = AliHLTMUONConstants::FullTrackerId();
	const char* decCompId = AliHLTMUONConstants::DecisionComponentId();
	
	if (rawReader != NULL)
	{
		// Implement the dHLT-sim-fromRaw dHLT simulation chain reading
		// from raw data.
		const char* rawPubComp = "AliRawReaderPublisher";
		const char* cmd1 = "-minid 2560 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000001";
		const char* cmd2 = "-minid 2561 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000002";
		const char* cmd3 = "-minid 2562 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000004";
		const char* cmd4 = "-minid 2563 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000008";
		const char* cmd5 = "-minid 2564 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000010";
		const char* cmd6 = "-minid 2565 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000020";
		const char* cmd7 = "-minid 2566 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000040";
		const char* cmd8 = "-minid 2567 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000080";
		const char* cmd9 = "-minid 2568 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000100";
		const char* cmd10 = "-minid 2569 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000200";
		const char* cmd11 = "-minid 2570 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000400";
		const char* cmd12 = "-minid 2571 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x000800";
		const char* cmd13 = "-minid 2572 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x001000";
		const char* cmd14 = "-minid 2573 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x002000";
		const char* cmd15 = "-minid 2574 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x004000";
		const char* cmd16 = "-minid 2575 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x008000";
		const char* cmd17 = "-minid 2576 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x010000";
		const char* cmd18 = "-minid 2577 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x020000";
		const char* cmd19 = "-minid 2578 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x040000";
		const char* cmd20 = "-minid 2579 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x080000";
		const char* cmd21 = "-minid 2816 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x100000";
		const char* cmd22 = "-minid 2817 -datatype 'DDL_RAW ' 'MUON' -dataspec 0x200000";
		handler->CreateConfiguration("RawDDL1", rawPubComp, NULL, cmd1);
		handler->CreateConfiguration("RawDDL2", rawPubComp, NULL, cmd2);
		handler->CreateConfiguration("RawDDL3", rawPubComp, NULL, cmd3);
		handler->CreateConfiguration("RawDDL4", rawPubComp, NULL, cmd4);
		handler->CreateConfiguration("RawDDL5", rawPubComp, NULL, cmd5);
		handler->CreateConfiguration("RawDDL6", rawPubComp, NULL, cmd6);
		handler->CreateConfiguration("RawDDL7", rawPubComp, NULL, cmd7);
		handler->CreateConfiguration("RawDDL8", rawPubComp, NULL, cmd8);
		handler->CreateConfiguration("RawDDL9", rawPubComp, NULL, cmd9);
		handler->CreateConfiguration("RawDDL10", rawPubComp, NULL, cmd10);
		handler->CreateConfiguration("RawDDL11", rawPubComp, NULL, cmd11);
		handler->CreateConfiguration("RawDDL12", rawPubComp, NULL, cmd12);
		handler->CreateConfiguration("RawDDL13", rawPubComp, NULL, cmd13);
		handler->CreateConfiguration("RawDDL14", rawPubComp, NULL, cmd14);
		handler->CreateConfiguration("RawDDL15", rawPubComp, NULL, cmd15);
		handler->CreateConfiguration("RawDDL16", rawPubComp, NULL, cmd16);
		handler->CreateConfiguration("RawDDL17", rawPubComp, NULL, cmd17);
		handler->CreateConfiguration("RawDDL18", rawPubComp, NULL, cmd18);
		handler->CreateConfiguration("RawDDL19", rawPubComp, NULL, cmd19);
		handler->CreateConfiguration("RawDDL20", rawPubComp, NULL, cmd20);
		handler->CreateConfiguration("RawDDL21", rawPubComp, NULL, cmd21);
		handler->CreateConfiguration("RawDDL22", rawPubComp, NULL, cmd22);
		const char* hrId = AliHLTMUONConstants::HitReconstructorId();
		const char* trId = AliHLTMUONConstants::TriggerReconstructorId();
		handler->CreateConfiguration("RecoRawDDL1", hrId, "RawDDL1", "-ddl 1 -cdb");
		handler->CreateConfiguration("RecoRawDDL2", hrId, "RawDDL2", "-ddl 2 -cdb");
		handler->CreateConfiguration("RecoRawDDL3", hrId, "RawDDL3", "-ddl 3 -cdb");
		handler->CreateConfiguration("RecoRawDDL4", hrId, "RawDDL4", "-ddl 4 -cdb");
		handler->CreateConfiguration("RecoRawDDL5", hrId, "RawDDL5", "-ddl 5 -cdb");
		handler->CreateConfiguration("RecoRawDDL6", hrId, "RawDDL6", "-ddl 6 -cdb");
		handler->CreateConfiguration("RecoRawDDL7", hrId, "RawDDL7", "-ddl 7 -cdb");
		handler->CreateConfiguration("RecoRawDDL8", hrId, "RawDDL8", "-ddl 8 -cdb");
		handler->CreateConfiguration("RecoRawDDL9", hrId, "RawDDL9", "-ddl 9 -cdb");
		handler->CreateConfiguration("RecoRawDDL10", hrId, "RawDDL10", "-ddl 10 -cdb");
		handler->CreateConfiguration("RecoRawDDL11", hrId, "RawDDL11", "-ddl 11 -cdb");
		handler->CreateConfiguration("RecoRawDDL12", hrId, "RawDDL12", "-ddl 12 -cdb");
		handler->CreateConfiguration("RecoRawDDL13", hrId, "RawDDL13", "-ddl 13 -cdb");
		handler->CreateConfiguration("RecoRawDDL14", hrId, "RawDDL14", "-ddl 14 -cdb");
		handler->CreateConfiguration("RecoRawDDL15", hrId, "RawDDL15", "-ddl 15 -cdb");
		handler->CreateConfiguration("RecoRawDDL16", hrId, "RawDDL16", "-ddl 16 -cdb");
		handler->CreateConfiguration("RecoRawDDL17", hrId, "RawDDL17", "-ddl 17 -cdb");
		handler->CreateConfiguration("RecoRawDDL18", hrId, "RawDDL18", "-ddl 18 -cdb");
		handler->CreateConfiguration("RecoRawDDL19", hrId, "RawDDL19", "-ddl 19 -cdb");
		handler->CreateConfiguration("RecoRawDDL20", hrId, "RawDDL20", "-ddl 20 -cdb");
		handler->CreateConfiguration("RecoRawDDL21", trId, "RawDDL21", "-ddl 21 -cdb -suppress_partial_triggers");
		handler->CreateConfiguration("RecoRawDDL22", trId, "RawDDL22", "-ddl 22 -cdb -suppress_partial_triggers");
		
		const char* recoSrcs = "RecoRawDDL13 RecoRawDDL14 RecoRawDDL15 RecoRawDDL16 RecoRawDDL17"
			" RecoRawDDL18 RecoRawDDL19 RecoRawDDL20 RecoRawDDL21 RecoRawDDL22";
		handler->CreateConfiguration("MansoTrackerForRaw", trackerId, recoSrcs, "");
		
		handler->CreateConfiguration("DecisionForRaw", decCompId, "MansoTrackerForRaw", "");
		
		TString outputSrcs = "DecisionForRaw MansoTrackerForRaw ";
		outputSrcs += recoSrcs;
		handler->CreateConfiguration("dHLT-sim-fromRaw", "BlockFilter", outputSrcs, "");
		
		// For reconstruction using full tracker.
		const char* recoSrcsFull = " RecoRawDDL1 RecoRawDDL2 RecoRawDDL3 RecoRawDDL4 RecoRawDDL5 RecoRawDDL6 "
			" RecoRawDDL7 RecoRawDDL8 RecoRawDDL9 RecoRawDDL10 RecoRawDDL11 RecoRawDDL12 "
			" RecoRawDDL13 RecoRawDDL14 RecoRawDDL15 RecoRawDDL16 RecoRawDDL17 "
			" RecoRawDDL18 RecoRawDDL19 RecoRawDDL20 RecoRawDDL21 RecoRawDDL22 ";
		handler->CreateConfiguration("FullTrackerForRaw", fullTrackerId, recoSrcsFull, "-cdb");
		
		handler->CreateConfiguration("DecisionForRawFullTrk", decCompId, "FullTrackerForRaw", "");
		
		TString outputSrcsFull = "DecisionForRawFullTrk FullTrackerForRaw ";
		outputSrcsFull += recoSrcsFull;
		handler->CreateConfiguration("dHLT-sim-fromRaw-fullTracker", "BlockFilter", outputSrcsFull, "");
	}
	
	if (IsMuonModuleLoaded() and runloader != NULL)
	{
		// Implement the dHLT-sim dHLT simulation chain reading from
		// simulated digits.
		const char* digitPub = AliHLTMUONConstants::DigitPublisherId();
		handler->CreateConfiguration("DigitDDL13", digitPub, NULL, "-simdata -ddl 13");
		handler->CreateConfiguration("DigitDDL14", digitPub, NULL, "-simdata -ddl 14");
		handler->CreateConfiguration("DigitDDL15", digitPub, NULL, "-simdata -ddl 15");
		handler->CreateConfiguration("DigitDDL16", digitPub, NULL, "-simdata -ddl 16");
		handler->CreateConfiguration("DigitDDL17", digitPub, NULL, "-simdata -ddl 17");
		handler->CreateConfiguration("DigitDDL18", digitPub, NULL, "-simdata -ddl 18");
		handler->CreateConfiguration("DigitDDL19", digitPub, NULL, "-simdata -ddl 19");
		handler->CreateConfiguration("DigitDDL20", digitPub, NULL, "-simdata -ddl 20");
		handler->CreateConfiguration("DigitDDL21", digitPub, NULL, "-simdata -ddl 21");
		handler->CreateConfiguration("DigitDDL22", digitPub, NULL, "-simdata -ddl 22");
		const char* hrId = AliHLTMUONConstants::HitReconstructorId();
		const char* trId = AliHLTMUONConstants::TriggerReconstructorId();
		handler->CreateConfiguration("RecoDDL13", hrId, "DigitDDL13", "-ddl 13 -cdb");
		handler->CreateConfiguration("RecoDDL14", hrId, "DigitDDL14", "-ddl 14 -cdb");
		handler->CreateConfiguration("RecoDDL15", hrId, "DigitDDL15", "-ddl 15 -cdb");
		handler->CreateConfiguration("RecoDDL16", hrId, "DigitDDL16", "-ddl 16 -cdb");
		handler->CreateConfiguration("RecoDDL17", hrId, "DigitDDL17", "-ddl 17 -cdb");
		handler->CreateConfiguration("RecoDDL18", hrId, "DigitDDL18", "-ddl 18 -cdb");
		handler->CreateConfiguration("RecoDDL19", hrId, "DigitDDL19", "-ddl 19 -cdb");
		handler->CreateConfiguration("RecoDDL20", hrId, "DigitDDL20", "-ddl 20 -cdb");
		handler->CreateConfiguration("RecoDDL21", trId, "DigitDDL21", "-ddl 21 -cdb -suppress_partial_triggers");
		handler->CreateConfiguration("RecoDDL22", trId, "DigitDDL22", "-ddl 22 -cdb -suppress_partial_triggers");
		
		const char* recoSrcs = "RecoDDL13 RecoDDL14 RecoDDL15 RecoDDL16 RecoDDL17"
			" RecoDDL18 RecoDDL19 RecoDDL20 RecoDDL21 RecoDDL22";
		handler->CreateConfiguration("MansoTracker", trackerId, recoSrcs, "");
		
		handler->CreateConfiguration("Decision", decCompId, "MansoTracker", "");
		
		TString outputSrcs = "Decision MansoTracker ";
		outputSrcs += recoSrcs;
		handler->CreateConfiguration("dHLT-sim", "BlockFilter", outputSrcs.Data(), "");
	
		// Implement the dHLT-sim-fromMC dHLT simulation chain reading
		// Monte Carlo geant hits and putting those into a tracker component.
		const char* rhsId = AliHLTMUONConstants::RecHitsSourceId();
		const char* trsId = AliHLTMUONConstants::TriggerRecordsSourceId();
		handler->CreateConfiguration("HitsDDL13", rhsId, NULL, "-simdata -plane left  -chamber 7");
		handler->CreateConfiguration("HitsDDL14", rhsId, NULL, "-simdata -plane right -chamber 7");
		handler->CreateConfiguration("HitsDDL15", rhsId, NULL, "-simdata -plane left  -chamber 8");
		handler->CreateConfiguration("HitsDDL16", rhsId, NULL, "-simdata -plane right -chamber 8");
		handler->CreateConfiguration("HitsDDL17", rhsId, NULL, "-simdata -plane left  -chamber 9");
		handler->CreateConfiguration("HitsDDL18", rhsId, NULL, "-simdata -plane right -chamber 9");
		handler->CreateConfiguration("HitsDDL19", rhsId, NULL, "-simdata -plane left  -chamber 10");
		handler->CreateConfiguration("HitsDDL20", rhsId, NULL, "-simdata -plane right -chamber 10");
		handler->CreateConfiguration("TrigRecsDDL21", trsId, NULL, "-hitdata -plane left");
		handler->CreateConfiguration("TrigRecsDDL22", trsId, NULL, "-hitdata -plane right");
		
		const char* dataSrcs = "HitsDDL13 HitsDDL14 HitsDDL15 HitsDDL16 HitsDDL17"
			" HitsDDL18 HitsDDL19 HitsDDL20 TrigRecsDDL21 TrigRecsDDL22";
		handler->CreateConfiguration("MansoTrackerForMC", trackerId, dataSrcs, "");
		
		handler->CreateConfiguration("DecisionForMC", decCompId, "MansoTrackerForMC", "");
		
		outputSrcs = "DecisionForMC MansoTrackerForMC ";
		outputSrcs += dataSrcs;
		handler->CreateConfiguration("dHLT-sim-fromMC", "BlockFilter", outputSrcs.Data(), "");
	}
	
	// Create a chain for generating AliESDEvent objects from dHLT raw reconstructed data.
	handler->CreateConfiguration("HLTOUTPubTrigRecs", "AliHLTOUTPublisher", NULL, "-datatype 'TRIGRECS' 'MUON'");
	handler->CreateConfiguration("HLTOUTPubMansoTracks", "AliHLTOUTPublisher", NULL, "-datatype 'MANTRACK' 'MUON'");
	handler->CreateConfiguration("HLTOUTPubTracks", "AliHLTOUTPublisher", NULL, "-datatype 'TRACKS  ' 'MUON'");
	handler->CreateConfiguration(
			"dHLT-make-esd",
			AliHLTMUONConstants::ESDMakerId(),
			"HLTOUTPubTrigRecs HLTOUTPubMansoTracks HLTOUTPubTracks",
			"-make_minimal_esd"
		);
	
	// Create a chain for rootifying the raw dHLT data and dumping to file.
	// This is used during AliRoot reconstruction.
	handler->CreateConfiguration("HLTOUTPubTrigDbg", "AliHLTOUTPublisher", NULL, "-datatype 'TRIGRDBG' 'MUON'");
	handler->CreateConfiguration("HLTOUTPubHits", "AliHLTOUTPublisher", NULL, "-datatype 'RECHITS ' 'MUON'");
	handler->CreateConfiguration("HLTOUTPubClusters", "AliHLTOUTPublisher", NULL, "-datatype 'CLUSTERS' 'MUON'");
	handler->CreateConfiguration("HLTOUTPubChannels", "AliHLTOUTPublisher", NULL, "-datatype 'CHANNELS' 'MUON'");
	handler->CreateConfiguration("HLTOUTPubCandidates", "AliHLTOUTPublisher", NULL, "-datatype 'MNCANDID' 'MUON'");
	handler->CreateConfiguration("HLTOUTPubSingles", "AliHLTOUTPublisher", NULL, "-datatype 'DECIDSIN' 'MUON'");
	handler->CreateConfiguration("HLTOUTPubPairs", "AliHLTOUTPublisher", NULL, "-datatype 'DECIDPAR' 'MUON'");
	handler->CreateConfiguration(
			"HLTOUTConverter",
			AliHLTMUONConstants::RootifierComponentId(),
			"HLTOUTPubTrigRecs HLTOUTPubTrigDbg HLTOUTPubHits HLTOUTPubClusters"
			 " HLTOUTPubChannels HLTOUTPubMansoTracks HLTOUTPubCandidates"
			 " HLTOUTPubTracks HLTOUTPubSingles HLTOUTPubPairs",
			""
		);
	handler->CreateConfiguration(
			"dHLT-rootify-and-dump",
			"ROOTFileWriter",
			"HLTOUTConverter",
			"-concatenate-events -datafile dHLTRawData.root -specfmt"
		);
	
	return 0;
}


int AliHLTMUONAgent::RegisterComponents(AliHLTComponentHandler* pHandler) const
{
	///
	/// Registers all available components of this module.
	/// @param pHandler  [in] instance of the component handler.
	///
	
	if (pHandler == NULL) return -EINVAL;
	pHandler->AddComponent(new AliHLTMUONRecHitsSource);
	pHandler->AddComponent(new AliHLTMUONTriggerRecordsSource);
	pHandler->AddComponent(new AliHLTMUONDigitPublisherComponent);
	pHandler->AddComponent(new AliHLTMUONRootifierComponent);
	pHandler->AddComponent(new AliHLTMUONHitReconstructorComponent);
	pHandler->AddComponent(new AliHLTMUONTriggerReconstructorComponent);
	pHandler->AddComponent(new AliHLTMUONMansoTrackerFSMComponent);
	pHandler->AddComponent(new AliHLTMUONFullTrackerComponent);
	pHandler->AddComponent(new AliHLTMUONDecisionComponent);
	pHandler->AddComponent(new AliHLTMUONESDMaker);
	pHandler->AddComponent(new AliHLTMUONEmptyEventFilterComponent);
	pHandler->AddComponent(new AliHLTMUONDataCheckerComponent);
	pHandler->AddComponent(new AliHLTMUONClusterFinderComponent);
	pHandler->AddComponent(new AliHLTMUONRawDataHistoComponent);
	pHandler->AddComponent(new AliHLTMUONClusterHistoComponent);
	return 0;
}


int AliHLTMUONAgent::GetHandlerDescription(
		AliHLTComponentDataType dt,
#ifdef __DEBUG
		AliHLTUInt32_t spec,
#else
		AliHLTUInt32_t /*spec*/,
#endif
		AliHLTOUTHandlerDesc& desc
	) const
{
	/// Get handler decription for MUON data in the HLTOUT data stream.
	
	if (dt == AliHLTMUONConstants::TriggerRecordsBlockDataType() or
	    dt == AliHLTMUONConstants::MansoTracksBlockDataType() or
	    dt == AliHLTMUONConstants::TracksBlockDataType()
	   )
	{
		HLTDebug("Indicating we can handle data type = %s and specification"
			" = 0x%8.8X with dHLT-make-esd chain",
			AliHLTComponent::DataType2Text(dt).c_str(),
			spec
		);
		desc = AliHLTOUTHandlerDesc(kChain, dt, "dHLT-make-esd");
		return 1;
	}
	
	if (dt == AliHLTMUONConstants::TriggerRecordsBlockDataType() or
	    dt == AliHLTMUONConstants::TrigRecsDebugBlockDataType() or
	    dt == AliHLTMUONConstants::RecHitsBlockDataType() or
	    dt == AliHLTMUONConstants::ClusterBlockDataType() or
	    dt == AliHLTMUONConstants::ChannelBlockDataType() or
	    dt == AliHLTMUONConstants::MansoTracksBlockDataType() or
	    dt == AliHLTMUONConstants::MansoCandidatesBlockDataType() or
	    dt == AliHLTMUONConstants::TracksBlockDataType() or
	    dt == AliHLTMUONConstants::SinglesDecisionBlockDataType() or
	    dt == AliHLTMUONConstants::PairsDecisionBlockDataType()
	   )
	{
		if (fgRunRootifyChain)
		{
			HLTDebug("Indicating we can handle data type = %s and specification"
				" = 0x%8.8X with dHLT-rootify-and-dump chain",
				AliHLTComponent::DataType2Text(dt).c_str(),
				spec
			);
			desc = AliHLTOUTHandlerDesc(kChain, dt, "dHLT-rootify-and-dump");
			return 1;
		}
		else
		{
			HLTDebug("Indicating we will ignore data type = %s and specification = 0x%8.8X",
				AliHLTComponent::DataType2Text(dt).c_str(),
				spec
			);
			desc = AliHLTOUTHandlerDesc(kProprietary, dt, "AliHLTOUTHandlerIgnore");
			return 1;
		}
	}

	if (dt == kAliHLTDataTypeEvent)
	{
		HLTDebug("Indicating we will ignore data type = %s and specification = 0x%8.8X",
			AliHLTComponent::DataType2Text(dt).c_str(),
			spec
		);
		desc = AliHLTOUTHandlerDesc(kProprietary, dt, "AliHLTOUTHandlerIgnore");
		return 1;
	}
	
	return 0;
}


AliHLTOUTHandler* AliHLTMUONAgent::GetOutputHandler(
		AliHLTComponentDataType dt,
#ifdef __DEBUG
		AliHLTUInt32_t spec
#else
		AliHLTUInt32_t /*spec*/
#endif
	)
{
	/// Get specific handler for MUON data in the HLTOUT data stream.
	
	HLTDebug("Trying to create HLTOUT handler for data type = %s and"
		" specification = 0x%8.8X",
		AliHLTComponent::DataType2Text(dt).c_str(),
		spec
	);
	
	if (dt == AliHLTMUONConstants::TriggerRecordsBlockDataType() or
	    dt == AliHLTMUONConstants::MansoTracksBlockDataType() or
	    dt == AliHLTMUONConstants::TracksBlockDataType()
	   )
	{
		return &fgkESDMakerChain;
	}

	if (dt == AliHLTMUONConstants::TriggerRecordsBlockDataType() or
	    dt == AliHLTMUONConstants::TrigRecsDebugBlockDataType() or
	    dt == AliHLTMUONConstants::RecHitsBlockDataType() or
	    dt == AliHLTMUONConstants::ClusterBlockDataType() or
	    dt == AliHLTMUONConstants::ChannelBlockDataType() or
	    dt == AliHLTMUONConstants::MansoTracksBlockDataType() or
	    dt == AliHLTMUONConstants::MansoCandidatesBlockDataType() or
	    dt == AliHLTMUONConstants::TracksBlockDataType() or
	    dt == AliHLTMUONConstants::SinglesDecisionBlockDataType() or
	    dt == AliHLTMUONConstants::PairsDecisionBlockDataType()
	   )
	{
		if (fgRunRootifyChain)
		{
			return &fgkRootifyDumpChain;
		}
		else
		{
			return &fgkDataIgnoreHandler;
		}
	}

	if (dt == kAliHLTDataTypeEvent)
	{
		return &fgkDataIgnoreHandler;
	}
	
	return NULL;
}


int AliHLTMUONAgent::DeleteOutputHandler(AliHLTOUTHandler* pInstance)
{
	/// Deletes the HLTOUT handlers. In this case since the handlers are
	/// allocated statically, we just check that the right pointer was
	/// given and exit.
	
	HLTDebug("Trying to delete HLTOUT handler: %p", pInstance);
	
	if (pInstance != &fgkESDMakerChain or pInstance != &fgkRootifyDumpChain
	    or pInstance != &fgkDataIgnoreHandler
	   )
	{
		return -EINVAL;
	}
	
	return 0;
}
