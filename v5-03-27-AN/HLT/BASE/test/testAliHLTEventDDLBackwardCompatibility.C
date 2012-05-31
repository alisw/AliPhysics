// $Id: $

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

/// @file   testAliHLTEventDDLBackwardCompatibility.C
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   25 June 2010
/// @brief  Test program for backward compatibility of the AliHLTEventDDL structure.
///

#if defined(__CINT__) && (! defined(__MAKECINT__))
#error This macro must be compiled. Try running as testAliHLTEventDDLBackwardCompatibility.C++.
#endif

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliHLTDataTypes.h"
#include "AliHLTReadoutList.h"
#include "AliHLTDAQ.h"
#include "AliHLTComponent.h"
#include "AliRawDataHeader.h"
#include "TRandom3.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "Riostream.h"
#include <vector>
#include <algorithm>
#endif

/**
 * Tests to see if the AliHLTReadoutList class handles both AliHLTEventDDLV0
 * and AliHLTEventDDLV1 correctly.
 */
bool CheckReadoutListConvertedCorrectly()
{
	// Initialise the different versions of the structures so that every detector
	// has only its first DDL set.
	union
	{
		AliHLTEventDDL eventddlV0;
		AliHLTEventDDLV0 bitsV0;
	};
	bitsV0.fCount = gkAliHLTDDLListSizeV0;
	if (gkAliHLTDDLListSizeV0 != 30)
	{
		cerr << "ERROR: gkAliHLTDDLListSizeV0 has a value of " << gkAliHLTDDLListSizeV0
			<< " but expected a value of 30." << endl;
		return false;
	}
	bitsV0.fList[0] = 0x00000001; // kITSSPD
	bitsV0.fList[1] = 0x00000001; // kITSSDD
	bitsV0.fList[2] = 0x00000001; // kITSSSD
	bitsV0.fList[3] = 0x00000001; // kTPC
	bitsV0.fList[4] = 0x00000000; // kTPC
	bitsV0.fList[5] = 0x00000000; // kTPC
	bitsV0.fList[6] = 0x00000000; // kTPC
	bitsV0.fList[7] = 0x00000000; // kTPC
	bitsV0.fList[8] = 0x00000000; // kTPC
	bitsV0.fList[9] = 0x00000000; // kTPC
	bitsV0.fList[10] = 0x00000000; // kTPC
	bitsV0.fList[11] = 0x00000001; // kTRD
	bitsV0.fList[12] = 0x00000001; // kTOF
	bitsV0.fList[13] = 0x00000000; // kTOF
	bitsV0.fList[14] = 0x00000000; // kTOF
	bitsV0.fList[15] = 0x00000001; // kHMPID
	bitsV0.fList[16] = 0x00000001; // kPHOS
	bitsV0.fList[17] = 0x00000001; // kCPV
	bitsV0.fList[18] = 0x00000001; // kPMD
	bitsV0.fList[19] = 0x00000001; // kMUONTRK
	bitsV0.fList[20] = 0x00000001; // kMUONTRG
	bitsV0.fList[21] = 0x00000001; // kFMD
	bitsV0.fList[22] = 0x00000001; // kT0
	bitsV0.fList[23] = 0x00000001; // kV0
	bitsV0.fList[24] = 0x00000001; // kZDC
	bitsV0.fList[25] = 0x00000001; // kACORDE
	bitsV0.fList[26] = 0x00000001; // kTRG
	bitsV0.fList[27] = 0x00000001; // kEMCAL
	bitsV0.fList[28] = 0x00000001; // kDAQTEST
	bitsV0.fList[29] = 0x00000001; // kHLT
	
	union
	{
		AliHLTEventDDL eventddlV1;
		AliHLTEventDDLV1 bitsV1;
	};
	bitsV1.fCount = gkAliHLTDDLListSizeV1;
	if (gkAliHLTDDLListSizeV1 != 31)
	{
		cerr << "ERROR: gkAliHLTDDLListSizeV1 has a value of " << gkAliHLTDDLListSizeV1
			<< " but expected a value of 31." << endl;
		return false;
	}
	bitsV1.fList[0] = 0x00000001; // kITSSPD
	bitsV1.fList[1] = 0x00000001; // kITSSDD
	bitsV1.fList[2] = 0x00000001; // kITSSSD
	bitsV1.fList[3] = 0x00000001; // kTPC
	bitsV1.fList[4] = 0x00000000; // kTPC
	bitsV1.fList[5] = 0x00000000; // kTPC
	bitsV1.fList[6] = 0x00000000; // kTPC
	bitsV1.fList[7] = 0x00000000; // kTPC
	bitsV1.fList[8] = 0x00000000; // kTPC
	bitsV1.fList[9] = 0x00000000; // kTPC
	bitsV1.fList[10] = 0x00000000; // kTPC
	bitsV1.fList[11] = 0x00000001; // kTRD
	bitsV1.fList[12] = 0x00000001; // kTOF
	bitsV1.fList[13] = 0x00000000; // kTOF
	bitsV1.fList[14] = 0x00000000; // kTOF
	bitsV1.fList[15] = 0x00000001; // kHMPID
	bitsV1.fList[16] = 0x00000001; // kPHOS
	bitsV1.fList[17] = 0x00000001; // kCPV
	bitsV1.fList[18] = 0x00000001; // kPMD
	bitsV1.fList[19] = 0x00000001; // kMUONTRK
	bitsV1.fList[20] = 0x00000001; // kMUONTRG
	bitsV1.fList[21] = 0x00000001; // kFMD
	bitsV1.fList[22] = 0x00000001; // kT0
	bitsV1.fList[23] = 0x00000001; // kV0
	bitsV1.fList[24] = 0x00000001; // kZDC
	bitsV1.fList[25] = 0x00000001; // kACORDE
	bitsV1.fList[26] = 0x00000001; // kTRG
	bitsV1.fList[27] = 0x00000001; // kEMCAL
	bitsV1.fList[28] = 0x00000000; // kEMCAL
	bitsV1.fList[29] = 0x00000001; // kDAQTEST
	bitsV1.fList[30] = 0x00000001; // kHLT
	
	AliHLTReadoutList rlV0(eventddlV0);
	AliHLTReadoutList rlV1(eventddlV1);
	
	// Check that for both readout list versions only the first DDLs are
	// enabled as expected.
	for (Int_t i = 0; i < AliHLTDAQ::NumberOfDetectors(); ++i)
	for (Int_t j = 0; j < AliHLTDAQ::NumberOfDdls(i); ++j)
	{
		Int_t ddlid = AliHLTDAQ::DdlIDOffset(i) | (j & 0xFF);
		if (j == 0)
		{
			if (rlV0.IsDDLDisabled(ddlid))
			{
				cerr << "ERROR: The first DDL for detector " << AliHLTDAQ::DetectorName(i)
					<< " was not enabled for readout list initialised from AliHLTEventDDLV0."
					<< endl;
				return false;
			}
			if (rlV1.IsDDLDisabled(ddlid))
			{
				cerr << "ERROR: The first DDL for detector " << AliHLTDAQ::DetectorName(i)
					<< " was not enabled for readout list initialised from AliHLTEventDDLV1."
					<< endl;
				return false;
			}
		}
		else
		{
			if (rlV0.IsDDLEnabled(ddlid))
			{
				cerr << "ERROR: DDL " << ddlid << " for detector " << AliHLTDAQ::DetectorName(i)
					<< " was marked enabled for readout list initialised from AliHLTEventDDLV0."
					<< endl;
				return false;
			}
			if (rlV1.IsDDLEnabled(ddlid))
			{
				cerr << "ERROR: DDL " << ddlid << " for detector " << AliHLTDAQ::DetectorName(i)
					<< " was marked enabled for readout list initialised from AliHLTEventDDLV1."
					<< endl;
				return false;
			}
		}
	}
	
	// Check that the internal buffers are identical.
	if (rlV0.BufferSize() != rlV1.BufferSize())
	{
		cerr << "ERROR: Buffer sizes for readout lists are different: rlV0.BufferSize() = "
			<< rlV0.BufferSize() << ", rlV1.BufferSize() = " << rlV1.BufferSize() << endl;
		return false;
	}
	if (memcmp(rlV0.Buffer(), rlV1.Buffer(), rlV0.BufferSize()) != 0)
	{
		cerr << "ERROR: Buffers for the two readout list versions are different." << endl;
		return false;
	}
	
	return true;
}

/**
 * Tests to see if AliHLTComponent::ExtractTriggerData recognises the old format
 * of AliHLTEventTriggerData and handles it correctly.
 */
bool CheckHandlingOfOldAliHLTEventTriggerData()
{
	// Prepare old structure format.
	struct AliHLTEventTriggerDataV0
	{
		AliHLTUInt8_t  fAttributes[gkAliHLTBlockDAttributeCount];
		AliHLTUInt64_t fHLTStatus;
		AliHLTUInt32_t fCommonHeaderWordCnt;
		AliHLTUInt32_t fCommonHeader[gkAliHLTCommonHeaderCount];
		AliHLTEventDDLV0 fReadoutListV0;
	};
	AliHLTEventTriggerDataV0 eventTrigData;
	memset(&eventTrigData, 0x0, sizeof(eventTrigData));
	eventTrigData.fCommonHeaderWordCnt = 8;
	eventTrigData.fHLTStatus = 0x123;
	eventTrigData.fReadoutListV0.fCount = gkAliHLTDDLListSizeV0;
	eventTrigData.fReadoutListV0.fList[0] = 0x00000001; // kITSSPD
	eventTrigData.fReadoutListV0.fList[1] = 0x00000001; // kITSSDD
	eventTrigData.fReadoutListV0.fList[2] = 0x00000001; // kITSSSD
	eventTrigData.fReadoutListV0.fList[3] = 0x00000001; // kTPC
	eventTrigData.fReadoutListV0.fList[4] = 0x00000000; // kTPC
	eventTrigData.fReadoutListV0.fList[5] = 0x00000000; // kTPC
	eventTrigData.fReadoutListV0.fList[6] = 0x00000000; // kTPC
	eventTrigData.fReadoutListV0.fList[7] = 0x00000000; // kTPC
	eventTrigData.fReadoutListV0.fList[8] = 0x00000000; // kTPC
	eventTrigData.fReadoutListV0.fList[9] = 0x00000000; // kTPC
	eventTrigData.fReadoutListV0.fList[10] = 0x00000000; // kTPC
	eventTrigData.fReadoutListV0.fList[11] = 0x00000001; // kTRD
	eventTrigData.fReadoutListV0.fList[12] = 0x00000001; // kTOF
	eventTrigData.fReadoutListV0.fList[13] = 0x00000000; // kTOF
	eventTrigData.fReadoutListV0.fList[14] = 0x00000000; // kTOF
	eventTrigData.fReadoutListV0.fList[15] = 0x00000001; // kHMPID
	eventTrigData.fReadoutListV0.fList[16] = 0x00000001; // kPHOS
	eventTrigData.fReadoutListV0.fList[17] = 0x00000001; // kCPV
	eventTrigData.fReadoutListV0.fList[18] = 0x00000001; // kPMD
	eventTrigData.fReadoutListV0.fList[19] = 0x00000001; // kMUONTRK
	eventTrigData.fReadoutListV0.fList[20] = 0x00000001; // kMUONTRG
	eventTrigData.fReadoutListV0.fList[21] = 0x00000001; // kFMD
	eventTrigData.fReadoutListV0.fList[22] = 0x00000001; // kT0
	eventTrigData.fReadoutListV0.fList[23] = 0x00000001; // kV0
	eventTrigData.fReadoutListV0.fList[24] = 0x00000001; // kZDC
	eventTrigData.fReadoutListV0.fList[25] = 0x00000001; // kACORDE
	eventTrigData.fReadoutListV0.fList[26] = 0x00000001; // kTRG
	eventTrigData.fReadoutListV0.fList[27] = 0x00000001; // kEMCAL
	eventTrigData.fReadoutListV0.fList[28] = 0x00000001; // kDAQTEST
	eventTrigData.fReadoutListV0.fList[29] = 0x00000001; // kHLT
	
	AliHLTComponentTriggerData trigData = {
		sizeof(AliHLTComponentTriggerData),
		sizeof(AliHLTEventTriggerDataV0),
		&eventTrigData
	};
	
	AliHLTReadoutList readoutlistExpected;
	for (Int_t i = 0; i < AliHLTDAQ::NumberOfDetectors(); ++i)
	{
		Int_t ddlid = AliHLTDAQ::DdlIDOffset(i);
		readoutlistExpected.EnableDDLBit(ddlid);
	}
	
	const AliHLTUInt8_t (*attribs)[gkAliHLTBlockDAttributeCount];
	AliHLTUInt64_t status = 0x0;
	const AliRawDataHeader* cdh = NULL;
	AliHLTReadoutList readoutlist;
	
	int result = AliHLTComponent::ExtractTriggerData(trigData, &attribs, &status, &cdh, &readoutlist, true);
	if (result != 0)
	{
		cerr << "ERROR: The method AliHLTComponent::ExtractTriggerData"
			" fails for the old structure format." << endl;
		return false;
	}
	if (attribs != &eventTrigData.fAttributes)
	{
		cerr << "ERROR: The method AliHLTComponent::ExtractTriggerData"
			" fails to locate the attributes structure correctly." << endl;
		return false;
	}
	if (status != 0x123)
	{
		cerr << "ERROR: The method AliHLTComponent::ExtractTriggerData"
			" fails to locate the HLT status word correctly." << endl;
		return false;
	}
	if ((const void*)cdh != (void*)&eventTrigData.fCommonHeader)
	{
		cerr << "ERROR: The method AliHLTComponent::ExtractTriggerData"
			" fails to locate the Common Data Header (CDH) structure correctly." << endl;
		return false;
	}
	if (memcmp(readoutlist.Buffer(), readoutlistExpected.Buffer(), readoutlist.BufferSize()) != 0)
	{
		cerr << "ERROR: The method AliHLTComponent::ExtractTriggerData"
			" fails to extract the readout list correctly." << endl;
		return false;
	}
	
	return true;
}

/**
 * Tests to see if reading old AliHLTReadoutList versions from root files works.
 * \param filename  The name of the file generated by GenerateReadoutListFile.C,
 *      which contains the readout list objects to test.
 */
bool CheckReadingOldFormat(const char* filename = "oldAliHLTReadoutListFormat.root")
{
	TFile file(filename, "READ");
	
	// Load the readout list objects.
	AliHLTReadoutList* rl[5] = {NULL, NULL, NULL, NULL, NULL};
	for (int i = 0; i < 5; ++i)
	{
		char name[1024];
		sprintf(name, "readoutlist%d", i+1);
		rl[i] = dynamic_cast<AliHLTReadoutList*>(file.Get(name));
		if (rl[i] == NULL)
		{
			cerr << "ERROR: Could not get object '" << name
				<< "' from file '" << filename << "'." << endl;
			return false;
		}
	}
	
	// Now load the tree and see that the objects are the same as the readout
	// list objects stored directly to the TFile.
	const char* treename = "rltree";
	TTree* tree = dynamic_cast<TTree*>(file.Get(treename));
	if (tree == NULL)
	{
		cerr << "ERROR: Could not get TTree '" << treename
			<< "' from file '" << filename << "'." << endl;
		return false;
	}
	AliHLTReadoutList* r = new AliHLTReadoutList;
	tree->SetBranchAddress("readoutlist", &r);
	for (int i = 0; i < 5; ++i)
	{
		tree->GetEvent(i);
		if (r->BufferSize() != rl[i]->BufferSize())
		{
			cerr << "ERROR: readoutlist" << i+1
				<< " and the one from the TTree have different sizes."
				<< endl;
			return false;
		}
#if 1 // ROOT_SVN_REVISION < 9999  //FIXME: after fixed https://savannah.cern.ch/bugs/?69241
		r->SetDDLBit(9999999, kTRUE);  // Triggers a reformating of the internal structure to the new version.
#endif
		if (memcmp(r->Buffer(), rl[i]->Buffer(), r->BufferSize()) != 0)
		{
			cerr << "ERROR: readoutlist" << i+1
				<< " and the one from the TTree are different."
				<< endl;
			return false;
		}
	}
	
	// Now check each readout list individually.
	typedef AliHLTReadoutList RL;
	Int_t alwaysoff = RL::kITSSPD
		| RL::kITSSDD
		| RL::kITSSSD
		| RL::kTPC
		| RL::kTRD
		| RL::kTOF
		| RL::kHMPID
		| RL::kPHOS
		| RL::kCPV
		| RL::kPMD
		| RL::kMUONTRK
		| RL::kMUONTRG
		| RL::kFMD
		| RL::kT0
		| RL::kV0
		| RL::kZDC
		| RL::kACORDE;
	
	// We will need to try and set the missing EMCAL DDL bits.
	// Otherwise the DetectorEnabled and DetectorDisabled methods will not
	// give the expected answers.
	for (int i = 24; i < AliHLTDAQ::NumberOfDdls("EMCAL"); ++i)
	{
		for (int j = 1; j <= 3; ++j)
		{
			Int_t ddlid = AliHLTDAQ::DdlIDOffset("EMCAL") | (i & 0xFF);
			rl[j]->EnableDDLBit(ddlid);
		}
	}
	
	int rlnum = 0;
	if (not (rl[rlnum]->DetectorEnabled(RL::kTRG) and
	         rl[rlnum]->DetectorDisabled(alwaysoff | RL::kEMCAL | RL::kDAQTEST | RL::kHLT)
	   ))
	{
		rl[0]->Print();
		cerr << "ERROR: readoutlist" << rlnum+1 << " does not have the correct bits set." << endl;
		rl[rlnum]->Print();
		return false;
	}
	rlnum = 1;
	if (not (rl[rlnum]->DetectorEnabled(RL::kTRG | RL::kEMCAL) and
	         rl[rlnum]->DetectorDisabled(alwaysoff | RL::kDAQTEST | RL::kHLT)
	   ))
	{
		cerr << "ERROR: readoutlist" << rlnum+1 << " does not have the correct bits set." << endl;
		rl[rlnum]->Print();
		return false;
	}
	rlnum = 2;
	if (not (rl[rlnum]->DetectorEnabled(RL::kTRG | RL::kEMCAL | RL::kDAQTEST) and
	         rl[rlnum]->DetectorDisabled(alwaysoff | RL::kHLT)
	   ))
	{
		cerr << "ERROR: readoutlist" << rlnum+1 << " does not have the correct bits set." << endl;
		rl[rlnum]->Print();
		return false;
	}
	rlnum = 3;
	if (not (rl[rlnum]->DetectorEnabled(RL::kTRG | RL::kEMCAL | RL::kDAQTEST | RL::kHLT) and
	         rl[rlnum]->DetectorDisabled(alwaysoff)
	   ))
	{
		cerr << "ERROR: readoutlist" << rlnum+1 << " does not have the correct bits set." << endl;
		rl[rlnum]->Print();
		return false;
	}
	rlnum = 4;
	if (not (rl[rlnum]->DetectorEnabled(RL::kTRG | RL::kDAQTEST | RL::kHLT) and
	         rl[rlnum]->DetectorDisabled(alwaysoff | RL::kEMCAL)
	   ))
	{
		cerr << "ERROR: readoutlist" << rlnum+1 << " does not have the correct bits set." << endl;
		rl[rlnum]->Print();
		return false;
	}
	
	return true;
}

/**
 * Runs the unit tests for backward compatibility of AliHLTEventDDL succeeded.
 * \returns true if the tests were passed and false otherwise.
 */
bool testAliHLTEventDDLBackwardCompatibility()
{
	if (not CheckReadoutListConvertedCorrectly()) return false;
	if (not CheckHandlingOfOldAliHLTEventTriggerData()) return false;
	if (not CheckReadingOldFormat()) return false;
	return true;
}

#ifndef __MAKECINT__

int main(int /*argc*/, const char** /*argv*/)
{
	bool resultOk = testAliHLTEventDDLBackwardCompatibility();
	if (not resultOk) return 1;
	return 0;
}

#endif // __MAKECINT__
