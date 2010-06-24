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

/// @file   testAliHLTReadoutList.C
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   8 June 2010
/// @brief  Test program for the AliHLTReadoutList class.
///

#if defined(__CINT__) && (! defined(__MAKECINT__))
#error This macro must be compiled. Try running as testAliHLTReadoutList.C++.
#endif

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliHLTDataTypes.h"
#include "AliHLTReadoutList.h"
#include "AliHLTDAQ.h"
#include "TRandom3.h"
#include "TString.h"
#include "Riostream.h"
#include <vector>
#include <algorithm>
#endif

// The detector codes as used by AliHLTReadoutList.
const int kgNumberOfCodes = 21;
const int kgDetCodes[kgNumberOfCodes] = {
	AliHLTReadoutList::kITSSPD,
	AliHLTReadoutList::kITSSDD,
	AliHLTReadoutList::kITSSSD,
	AliHLTReadoutList::kTPC,
	AliHLTReadoutList::kTRD,
	AliHLTReadoutList::kTOF,
	AliHLTReadoutList::kHMPID,
	AliHLTReadoutList::kPHOS,
	AliHLTReadoutList::kCPV,
	AliHLTReadoutList::kPMD,
	AliHLTReadoutList::kMUONTRK,
	AliHLTReadoutList::kMUONTRG,
	AliHLTReadoutList::kFMD,
	AliHLTReadoutList::kT0,
	AliHLTReadoutList::kV0,
	AliHLTReadoutList::kZDC,
	AliHLTReadoutList::kACORDE,
	AliHLTReadoutList::kTRG,
	AliHLTReadoutList::kEMCAL,
	AliHLTReadoutList::kDAQTEST,
	AliHLTReadoutList::kHLT
};
const char* kgDetCodeName[kgNumberOfCodes] = {
	"AliHLTReadoutList::kITSSPD",
	"AliHLTReadoutList::kITSSDD",
	"AliHLTReadoutList::kITSSSD",
	"AliHLTReadoutList::kTPC",
	"AliHLTReadoutList::kTRD",
	"AliHLTReadoutList::kTOF",
	"AliHLTReadoutList::kHMPID",
	"AliHLTReadoutList::kPHOS",
	"AliHLTReadoutList::kCPV",
	"AliHLTReadoutList::kPMD",
	"AliHLTReadoutList::kMUONTRK",
	"AliHLTReadoutList::kMUONTRG",
	"AliHLTReadoutList::kFMD",
	"AliHLTReadoutList::kT0",
	"AliHLTReadoutList::kV0",
	"AliHLTReadoutList::kZDC",
	"AliHLTReadoutList::kACORDE",
	"AliHLTReadoutList::kTRG",
	"AliHLTReadoutList::kEMCAL",
	"AliHLTReadoutList::kDAQTEST",
	"AliHLTReadoutList::kHLT"
};

/**
 * Converts a code to string.
 * \param code  The ID code of the detector. One of AliHLTReadoutList::EDetectorId
 * \returns the code name as a string given a the detector code.
 */
const char* CodeToString(int code)
{
	for (int i = 0; i < kgNumberOfCodes; ++i)
	{
		if (kgDetCodes[i] == code) return kgDetCodeName[i];
	}
	return "UNKNOWN";
}

/**
 * Checks if the basic empty and clear methods work.
 */
bool CheckEmptyAndClear()
{
	AliHLTReadoutList rl;
	if (rl.Empty() != true)
	{
		cerr << "ERROR: AliHLTReadoutList::Empty returns false for an empty readout list." << endl;
		return false;
	}
	
	// Enable all the detectors and check this operation.
	rl.Enable(AliHLTReadoutList::kALLDET);
	for (int i = 0; i < kgNumberOfCodes; ++i)
	{
		if (i == 19) continue; // This is the test DDL. off by default.
		if (not rl.DetectorEnabled(kgDetCodes[i]))
		{
			cerr << "ERROR: AliHLTReadoutList::Enable(AliHLTReadoutList::kALLDET) did not enable for "
				<< CodeToString(kgDetCodes[i]) << "." << endl;
			return false;
		}
	}
	if (rl.DetectorEnabled(AliHLTReadoutList::kDAQTEST))
	{
		cerr << "ERROR: AliHLTReadoutList::Enable(AliHLTReadoutList::kALLDET) enabled bits"
			" for AliHLTReadoutList::kDAQTEST but should not have." << endl;
		return false;
	}
	
	rl.Clear();
	// Fetch the raw bits for the readout list structure and check that they
	// are all zero, since we should have disabled everything in the loop above.
	AliHLTEventDDL bits = rl;
	if (bits.fCount != (unsigned int)gkAliHLTDDLListSize)
	{
		cerr << "ERROR: Typecast operator AliHLTEventDDL () is not"
			" setting the fCount of the structure correctly." << endl;
		return false;
	}
	for (int j = 0; j < gkAliHLTDDLListSize; ++j)
	{
		if (bits.fList[j] != 0x0)
		{
			cerr << "ERROR: Word " << j << " in internal AliHLTReadoutList"
				" bitfield structure is not zero as expected after a"
				" call to AliHLTReadoutList::Clear." << endl;
			return false;
		}
	}
	
	return true;
}

/**
 * Tests enabling and disabling of different detectors.
 */
bool CheckEnablingDisabling()
{
	for (int i = 0; i < 10000; ++i)
	{
		// Get 3 random detector codes.
		int detNum[3] = {
			gRandom->Integer(kgNumberOfCodes),
			gRandom->Integer(kgNumberOfCodes),
			gRandom->Integer(kgNumberOfCodes)
		};
		int code[3] = {
			kgDetCodes[detNum[0]],
			kgDetCodes[detNum[1]],
			kgDetCodes[detNum[2]]
		};
		// make sure the codes are not duplicated.
		while (code[1] == code[0])
		{
			detNum[1] = gRandom->Integer(kgNumberOfCodes);
			code[1] = kgDetCodes[detNum[1]];
		}
		while (code[2] == code[1] or code[2] == code[0])
		{
			detNum[2] = gRandom->Integer(kgNumberOfCodes);
			code[2] = kgDetCodes[detNum[2]];
		}
		
		// Choose the number of codes to use, from 1 to max 3.
		int codeCount = gRandom->Integer(3) + 1;
		
		// Build up the detector code list for the AliHLTReadoutList constructor.
		int totalCode = 0;
		for (int j = 0; j < codeCount; ++j) totalCode |= code[j];
		
		AliHLTReadoutList rl(totalCode);
		if (rl.Empty() == true)
		{
			cerr << "ERROR: AliHLTReadoutList::Empty returns true for a non empty readout list." << endl;
			return false;
		}
		
		// Check that the correct detectors have been enabled and
		// that we can disable a detector correctly.
		for (int j = 0; j < codeCount; ++j)
		{
			if (rl.DetectorEnabled(code[j]) == false)
			{
				cerr << "ERROR: Detector was not enabled for "
					<< CodeToString(code[j]) << " by constructor." << endl;
				return false;
			}
			
			// Also check each bit individualy according to AliHLTDAQ values.
			int det = detNum[j];
			int maxddls = AliHLTDAQ::NumberOfDdls(det);
			for (int ddlindex = 0; ddlindex < maxddls; ++ddlindex)
			{
				int ddlid = AliHLTDAQ::DdlIDOffset(det) | (ddlindex & 0xFF);
				if (rl.IsDDLDisabled(ddlid))
				{
					cerr << "ERROR: Bit not set for DDL " << ddlid
						<< ", even though detector "
						<< AliHLTDAQ::OnlineName(det)
						<< " was enabled." << endl;
					return false;
				}
			}
			
			rl.Disable(code[j]);
			if (rl.DetectorEnabled(code[j]) == true)
			{
				cerr << "ERROR: AliHLTReadoutList::Disable(x) is not working for x = "
					<< CodeToString(code[j]) << "." << endl;
				return false;
			}
		}
		
		// Fetch the raw bits for the readout list structure and check that they
		// are all zero, since we should have disabled everything in the loop above.
		AliHLTEventDDL bits = rl;
		for (int j = 0; j < gkAliHLTDDLListSize; ++j)
		{
			if (bits.fList[j] != 0x0)
			{
				cerr << "ERROR: Word " << j << " in internal AliHLTReadoutList"
					" bitfield structure is not zero as expected." << endl;
				return false;
			}
		}
	}
	return true;
}

/**
 * Tests enabling and disabling of different DDLs.
 */
bool CheckEnablingDisablingDDLs()
{
	for (int i = 0; i < 10000; ++i)
	{
		// Get random DDL IDs that are each unique.
		std::vector<int> ddls;
		int ddlCount = gRandom->Integer(100) + 1;
		for (int j = 0; j < ddlCount; ++j)
		{
			int ddlid = -1;
			do
			{
				int det = gRandom->Integer(AliHLTDAQ::NumberOfDetectors());
				int maxddls = AliHLTDAQ::NumberOfDdls(det);
				// The following is a special check since AliDAQ could be the newer version.
				// det = 18 is for EMCAL and gkAliHLTDDLListSize = 30 indicates old version
				// of AliHLTEventDDL before EMCAL expantion with DCAL.
				if (det == 18 and gkAliHLTDDLListSize == 30 and maxddls > 24) maxddls = 24;
				int ddlindex = gRandom->Integer(maxddls);
				ddlid = AliHLTDAQ::DdlID(det, ddlindex);
				if (std::find(ddls.begin(), ddls.end(), ddlid) != ddls.end()) ddlid = -1;
			}
			while (ddlid == -1);
			ddls.push_back(ddlid);
		}
		
		// Build up the enable list for the AliHLTReadoutList constructor.
		TString enableList;
		for (size_t j = 0; j < ddls.size(); ++j)
		{
			char num[32];
			sprintf(num, "%d", ddls[j]);
			enableList += " ";
			enableList += num;
		}
		
		AliHLTReadoutList rl(enableList.Data());
		if (rl.Empty() == true)
		{
			cerr << "ERROR: AliHLTReadoutList::Empty returns true for a non empty readout list." << endl;
			return false;
		}
		
		// Check that the correct DDLs have been enabled and
		// that we can disable a DDL correctly.
		for (size_t j = 0; j < ddls.size(); ++j)
		{
			if (rl.IsDDLEnabled(ddls[j]) == false)
			{
				cerr << "ERROR: DDL " << ddls[j] << " was not enabled by constructor." << endl;
				return false;
			}
			rl.DisableDDLBit(ddls[j]);
			if (rl.IsDDLDisabled(ddls[j]) == false)
			{
				cerr << "ERROR: AliHLTReadoutList::DisableDDLBit(x) is not working for x = "
					<< ddls[j] << "." << endl;
				return false;
			}
		}
		
		// Fetch the raw bits for the readout list structure and check that they
		// are all zero, since we should have disabled everything in the loop above.
		AliHLTEventDDL bits = rl;
		for (int j = 0; j < gkAliHLTDDLListSize; ++j)
		{
			if (bits.fList[j] != 0x0)
			{
				cerr << "ERROR: Word " << j << " in internal AliHLTReadoutList"
					" bitfield structure is not zero as expected." << endl;
				return false;
			}
		}
	}
	return true;
}

/**
 * Tests if using incorrect DDL IDs returns zero or is ignored as expected.
 */
bool CheckIncorrectIDs()
{
	for (int i = 0; i < 1000000; ++i)
	{
		// Get random DDL ID outside the valid range.
		int ddlid = -1;
		int det = gRandom->Integer(AliHLTDAQ::NumberOfDetectors()+1);
		if (det != AliHLTDAQ::NumberOfDetectors())
		{
			int maxddls = AliHLTDAQ::NumberOfDdls(det);
			int ddlindex = gRandom->Integer(0xFF - maxddls) + maxddls;
			ddlid = AliHLTDAQ::DdlIDOffset(det) | (ddlindex & 0xFF);
		}
		else
		{
			det = gRandom->Integer(11) + 20;
			if (det == 30) det = 31;
			int ddlindex = gRandom->Integer(0xFF);
			ddlid = (det << 8) | (ddlindex & 0xFF);
		}
		
		AliHLTReadoutList rl;
		if (rl.GetDDLBit(ddlid) != kFALSE)
		{
			cerr << "ERROR: Received a non zero result for invalid DDL " << ddlid << "." << endl;
			return false;
		}
		AliHLTEventDDL before = rl;
		rl.EnableDDLBit(ddlid);
		AliHLTEventDDL after = rl;
		if (memcmp(&before, &after, sizeof(AliHLTEventDDL)) != 0)
		{
			cerr << "ERROR: Modified AliHLTReadoutList structure using an invalid DDL "
				<< ddlid << "." << endl;
			cerr << "========== Dump of bits before modification ==========" << endl;
			for (unsigned int j = 0; j < before.fCount; ++j)
			{
				cerr << "[word " << dec << j << "] = " << hex << before.fList[j] << dec << endl;
			}
			cerr << "========== Dump of bits after modification ==========" << endl;
			for (unsigned int j = 0; j < after.fCount; ++j)
			{
				cerr << "[word " << dec << j << "] = " << hex << after.fList[j] << dec << endl;
			}
			return false;
		}
	}
	return true;
}

/**
 * Runs the unit test for the AliHLTReadoutList class.
 * \returns true if the class passed the test and false otherwise.
 */
bool testAliHLTReadoutList()
{
	gRandom->SetSeed(123);
	
	if (AliHLTDAQ::NumberOfDetectors() != kgNumberOfCodes)
	{
		cerr << "ERROR: The testAliHLTReadoutList.C macro needs to be updated"
			" or something went really wrong."
			" The number of detectors reported by AliHLTDAQ is not "
			<< kgNumberOfCodes << ", as expected, but "
			<< AliHLTDAQ::NumberOfDetectors() << "." << endl;
		return false;
	}
	
	if (not CheckEmptyAndClear()) return false;
	if (not CheckEnablingDisabling()) return false;
	if (not CheckEnablingDisablingDDLs()) return false;
	if (not CheckIncorrectIDs()) return false;
	
	return true;
}

#ifndef __MAKECINT__

int main(int /*argc*/, const char** /*argv*/)
{
	bool resultOk = testAliHLTReadoutList();
	if (not resultOk) return 1;
	return 0;
}

#endif // __MAKECINT__
