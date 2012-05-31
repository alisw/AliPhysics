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
 * @file   testAliHLTTPCDefinitions.C
 * @author Artur Szostak <artursz@iafrica.com>
 * @date   4 Aug 2010
 *
 * This macro is used to test the AliHLTTPCDefinitions class.
 * Specifically the mapping encoded in the different methods.
 */

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Riostream.h"
#include "TClassTable.h"
#include "TSystem.h"
#include "AliDAQ.h"
#include "AliHLTTPCDefinitions.h"
#endif


/**
 * Routine to check that the SlicePatchToDDLId and DDLIdToSlicePatch methods work.
 */
bool CheckDDLToSlicePatchConversion()
{
	// Check the conversion from slice + patch to DDL ID.
	Int_t minDDL = AliDAQ::DdlIDOffset("TPC");
	Int_t maxDDL = AliDAQ::DdlIDOffset("TPC") + AliDAQ::NumberOfDdls("TPC") - 1;
	for (AliHLTUInt16_t slice = 0; slice < 256; ++slice)
	for (AliHLTUInt16_t patch = 0; patch < 256; ++patch)
	{
		Int_t ddlid = AliHLTTPCDefinitions::SlicePatchToDDLId(AliHLTUInt8_t(slice), AliHLTUInt8_t(patch));
		if (slice < 36 && patch < 6)
		{
			// The slice and patch are valid so they should give a valid DDL ID.
			if (ddlid < minDDL || maxDDL < ddlid)
			{
				cerr << "ERROR: AliHLTTPCDefinitions::SlicePatchToDDLId("
					<< int(slice) << ", " << int(patch) << ") returned invalid DDL ID of "
					<< ddlid << " which is outside the valid range of ["
					<< minDDL << ".." << maxDDL << "]." << endl;
				return false;
			}
		}
		else
		{
			// The slice or patch are not valid so the result should be -1
			if (ddlid != -1)
			{
				cerr << "ERROR: AliHLTTPCDefinitions::SlicePatchToDDLId("
					<< int(slice) << ", " << int(patch) << ") returned invalid responce of "
					<< ddlid << " when is should have returned -1." << endl;
				return false;
			}
		}
	}
	
	// Check the conversion from DDL ID to slice + patch.
	for (AliHLTInt32_t ddlid2 = 0; ddlid2 < 8000; ++ddlid2)
	{
		AliHLTUInt8_t slice2 = 255;
		AliHLTUInt8_t patch2 = 255;
		bool result = AliHLTTPCDefinitions::DDLIdToSlicePatch(ddlid2, slice2, patch2);
		if (minDDL <= ddlid2 && ddlid2 <= maxDDL)
		{
			// The DDL ID was valid so the slice and patch should also be.
			if (result == false)
			{
				cerr << "ERROR: AliHLTTPCDefinitions::DDLIdToSlicePatch("
					<< ddlid2 << ") returned invalid result of 'false'."
					<< " But it should have been 'true'." << endl;
				return false;
			}
			if (slice2 > 35)
			{
				cerr << "ERROR: AliHLTTPCDefinitions::DDLIdToSlicePatch("
					<< ddlid2 << ") returned invalid slice value of "
					<< int(slice2) << ", but the valid range is [0..35]."
					<< endl;
				return false;
			}
			if (patch2 > 5)
			{
				cerr << "ERROR: AliHLTTPCDefinitions::DDLIdToSlicePatch("
					<< ddlid2 << ") returned invalid patch value of "
					<< int(patch2) << ", but the valid range is [0..5]."
					<< endl;
				return false;
			}
			if (AliHLTTPCDefinitions::SlicePatchToDDLId(slice2, patch2) != ddlid2)
			{
				cerr << "ERROR: AliHLTTPCDefinitions::DDLIdToSlicePatch("
					<< ddlid2 << ") generated a slice value of "
					<< int(slice2) << " and patch value of "
					<< int(patch2) << ", which gives a different DDL ID with the"
					<< " inverse map given by AliHLTTPCDefinitions::SlicePatchToDDLId("
					<< int(slice2) << ", " << int(patch2) << ")."
					<< endl;
				return false;
			}
		}
		else
		{
			// The DDL ID was not valid so the responce should be false.
			if (result == true)
			{
				cerr << "ERROR: AliHLTTPCDefinitions::DDLIdToSlicePatch("
					<< ddlid2 << ") returned invalid result of 'true'."
					<< " But it should have been 'false'." << endl;
				return false;
			}
		}
	}
	return true;
}

/**
 * This is the top level testing method which calls individual tests.
 * \returns true if all tests succeeded and false otherwise.
 */
bool testAliHLTTPCDefinitions()
{
	if (gClassTable->GetID("AliHLTTPCDefinitions") < 0)
	{
		gSystem->Load("libAliHLTUtil.so");
		gSystem->Load("libAliHLTTPC.so");
	}
	if (! CheckDDLToSlicePatchConversion()) return false;
	return true;
}

#ifndef __MAKECINT__

int main(int /*argc*/, const char** /*argv*/)
{
	bool resultOk = testAliHLTTPCDefinitions();
	if (not resultOk) return 1;
	return 0;
}

#endif // __MAKECINT__

