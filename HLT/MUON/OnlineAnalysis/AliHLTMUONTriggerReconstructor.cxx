/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Indranil Das <indra.das@saha.ac.in>                                  *
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

/* $Id$ */

/**********************************************************************
 Created on : 16/05/2007
 Purpose    : This class reads the tracker DDL files and gives the output
              as AliMUONTriggerRecordStruct structures.
 Author     : Indranil Das, HEP Division, SINP
 Email      : indra.das@saha.ac.in | indra.ehep@gmail.com

 Artur Szostak <artursz@iafrica.com>:
  Completely reimplemented the lookup table to a simplified format.
**********************************************************************/

///
///  @file   AliHLTMUONTriggerReconstructor.cxx
///  @author Indranil Das <indra.das@saha.ac.in>,
///          Artur Szostak <artursz@iafrica.com>
///  @date   16 May 2007
///  @brief  Implementation of the AliHLTMUONTriggerReconstructor class.
///
///  The trigger reconstructor class is designed to deal the rawdata inputfiles
///  to findout the the reconstructed hits at the trigger DDL. The output is send
///  to the output block for further processing.
///

#include "AliHLTMUONTriggerReconstructor.h"
#include "AliHLTMUONTriggerRecordsBlockStruct.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONCalculations.h"
#include "AliRawDataHeader.h"
#include <vector>
#include <cassert>


AliHLTMUONTriggerReconstructor::AliHLTMUONTriggerReconstructor() :
	AliHLTLogging(),
	fDecoder()
{
	/// Default constructor.
}


AliHLTMUONTriggerReconstructor::~AliHLTMUONTriggerReconstructor()
{
	/// Default destructor.
}


bool AliHLTMUONTriggerReconstructor::Run(
		const AliHLTUInt8_t* rawData,
		AliHLTUInt32_t rawDataSize,
		AliHLTMUONTriggerRecordStruct* trigRecord,
		AliHLTUInt32_t& nofTrigRec,
		bool suppressPartialTrigs
	)
{
	/// Runs the trigger reconstruction algorithm on the raw data.
	
	// Reset and initialise some variables in the decoder.
	fDecoder.GetHandler().MaxOutputTrigRecs(nofTrigRec);
	fDecoder.GetHandler().OutputTrigRecs(trigRecord);
	fDecoder.GetHandler().SuppressPartialTriggers(suppressPartialTrigs);
	
	fDecoder.Decode(rawData, rawDataSize);
	
	// nofTrigRec now becomes the output of how many trigger records were found.
	nofTrigRec = fDecoder.GetHandler().OutputTrigRecsCount();
	
	return not fDecoder.GetHandler().OverflowedOutputBuffer();
}


AliHLTMUONTriggerReconstructor::AliDecoderHandler::AliDecoderHandler() :
	AliMUONTriggerDDLDecoderEventHandler(),
	AliHLTLogging(),
	fLookupTable(),
	fBufferStart(NULL),
	fMaxOutputTrigRecs(0),
	fOutputTrigRecsCount(0),
	fOutputTrigRecs(NULL),
	fTrigRecId(0),
	fCurrentRegional(0),
	fCurrentLocal(0),
	fSuppressPartialTriggers(false),
	fOverflowed(false)
{
	/// Default constructor just resets the lookup table to zero.
	
	for (Int_t i = 0; i < 8; i++)
	for (Int_t j = 0; j < 16; j++)
	for (Int_t k = 0; k < 4; k++)
	for (Int_t n = 0; n < 2; n++)
	for (Int_t m = 0; m < 16; m++)
	{
		fLookupTable.fRow[i][j][k][n][m].fX = 0;
		fLookupTable.fRow[i][j][k][n][m].fY = 0;
		fLookupTable.fRow[i][j][k][n][m].fZ = 0;
	}
}


void AliHLTMUONTriggerReconstructor::AliDecoderHandler::OnNewBuffer(
		const void* buffer, UInt_t /*bufferSize*/
	)
{
	/// Called for each new buffer. Sets the buffer and resets the structure
	/// counters.
	
	assert( buffer != NULL );
	fBufferStart = buffer;
	
	// Start from -1 since we increment immediately in OnNewRegionalStruct.
	fCurrentRegional = fCurrentLocal = -1;
}


void AliHLTMUONTriggerReconstructor::AliDecoderHandler::OnLocalStruct(
		const AliMUONLocalInfoStruct* localStruct,
		const AliMUONLocalScalarsStruct* /*scalars*/
	)
{
	/// Converts a local trigger structure from the L0 into a trigger record.
	/// The dHLT trigger records is then used as a seed for tracking algorithms.
	/// \note fOutputTrigRecs must be set before calling the decoder to decode
	///    a new raw data buffer.
	/// \param localStruct  This is a pointer to the local L0 trigger structure data.

	assert(localStruct != NULL);
	assert(fOutputTrigRecs != NULL);

	fCurrentLocal++;
	AliHLTInt32_t iReg = fCurrentRegional;
	AliHLTInt32_t iLoc = fCurrentLocal;
	assert(iReg >= 0);
	assert(iLoc >= 0);

	// Check if there is anything in the trigger patterns at all.
	// If nothing then ignore this local L0 trigger.
	if (localStruct->fX2X1 == 0 and localStruct->fX4X3 == 0 and
	    localStruct->fY2Y1 == 0 and localStruct->fY4Y3 == 0
	   )
	{
		return;
	}
	
	// Check that we will not overflow the output buffer.
	if (fOutputTrigRecsCount >= fMaxOutputTrigRecs)
	{
		HLTError("Output buffer has overflowed maximum element count of %d.",
			fMaxOutputTrigRecs
		);
		fOverflowed = true;
		return;
	}
	
	UShort_t pattern[2][4]; // 2 stands for two cathode planes and the 4 stands for 4 chambers.
	pattern[0][0] = GetLocalX1(localStruct); // x-strip pattern for chamber 0
	pattern[0][1] = GetLocalX2(localStruct); // x-strip pattern for chamber 1
	pattern[0][2] = GetLocalX3(localStruct); // x-strip pattern for chamber 2
	pattern[0][3] = GetLocalX4(localStruct); // x-strip pattern for chamber 3
	pattern[1][0] = GetLocalY1(localStruct); // y-strip pattern for chamber 0
	pattern[1][1] = GetLocalY2(localStruct); // y-strip pattern for chamber 1
	pattern[1][2] = GetLocalY3(localStruct); // y-strip pattern for chamber 2
	pattern[1][3] = GetLocalY4(localStruct); // y-strip pattern for chamber 3

	bool setX[4] = {false, false, false, false};
	bool setY[4] = {false, false, false, false};

	for (int iChamber = 0; iChamber < 4; iChamber++) //4 chambers
	for (int iPlane = 0; iPlane < 2; iPlane++) // 2 cathode planes
	{
		for (Int_t ibitxy = 0; ibitxy < 16; ++ibitxy)
		{
			if (((pattern[iPlane][iChamber] >> ibitxy) & 0x1) != 0x1)
				continue;
			
			if (iPlane == 1)
			{
				fOutputTrigRecs[fOutputTrigRecsCount].fHit[iChamber].fX =
					fLookupTable.fRow[iReg][iLoc][iChamber][iPlane][ibitxy].fX;
				setX[iChamber] = true;
			}
			else
			{
				fOutputTrigRecs[fOutputTrigRecsCount].fHit[iChamber].fY =
					fLookupTable.fRow[iReg][iLoc][iChamber][iPlane][ibitxy].fY;
				fOutputTrigRecs[fOutputTrigRecsCount].fHit[iChamber].fZ =
					fLookupTable.fRow[iReg][iLoc][iChamber][iPlane][ibitxy].fZ;
				setY[iChamber] = true;
			}
		}
	}

	// hitset indicates which hits on chambers 7 to 10 have been found and filled.
	bool hitset[4] = {false, false, false, false};
	
	// Fill the hitset flags and make sure the hit structures that were not
	// filled (set) get set to a nil value.
	for (int i = 0; i < 4; i++)
	{
		hitset[i] = setX[i] and setY[i];
		
		if (not hitset[i])
		{
			fOutputTrigRecs[fOutputTrigRecsCount].fHit[i]
				= AliHLTMUONConstants::NilRecHitStruct();
		}
	}

	fOutputTrigRecs[fOutputTrigRecsCount].fId = fTrigRecId;

	// Increment trigger record Id and keep it positive.
	if (fTrigRecId < 0x7FFFFFFF)
		fTrigRecId++;
	else
		fTrigRecId = 0;
	
	AliHLTMUONRecHitStruct* hit1 = NULL;
	if (hitset[0])
		hit1 = &fOutputTrigRecs[fOutputTrigRecsCount].fHit[0];
	else if (hitset[1])
		hit1 = &fOutputTrigRecs[fOutputTrigRecsCount].fHit[1];
	AliHLTMUONRecHitStruct* hit2 = NULL;
	if (hitset[2])
		hit2 = &fOutputTrigRecs[fOutputTrigRecsCount].fHit[2];
	else if (hitset[3])
		hit2 = &fOutputTrigRecs[fOutputTrigRecsCount].fHit[3];
	
	if (hit1 != NULL and hit2 != NULL)
	{
		// Calculate the momentum and fill in the flags and momentum fields.
		AliHLTMUONCalculations::ComputeMomentum(
				hit1->fX,
				hit1->fY, hit2->fY,
				hit1->fZ, hit2->fZ
			);
		fOutputTrigRecs[fOutputTrigRecsCount].fPx = AliHLTMUONCalculations::Px();
		fOutputTrigRecs[fOutputTrigRecsCount].fPy = AliHLTMUONCalculations::Py();
		fOutputTrigRecs[fOutputTrigRecsCount].fPz = AliHLTMUONCalculations::Pz();

		fOutputTrigRecs[fOutputTrigRecsCount].fFlags =
			AliHLTMUONUtils::PackTriggerRecordFlags(
				AliHLTMUONCalculations::Sign(),
				hitset
			);
		
		fOutputTrigRecsCount++;
	}
	else if ((hit1 != NULL or hit2 != NULL) and not fSuppressPartialTriggers)
	{
		fOutputTrigRecs[fOutputTrigRecsCount].fPx = 0;
		fOutputTrigRecs[fOutputTrigRecsCount].fPy = 0;
		fOutputTrigRecs[fOutputTrigRecsCount].fPz = 0;

		fOutputTrigRecs[fOutputTrigRecsCount].fFlags =
			AliHLTMUONUtils::PackTriggerRecordFlags(
				kSignUnknown,
				hitset
			);
		
		fOutputTrigRecsCount++;
	}
}


void AliHLTMUONTriggerReconstructor::AliDecoderHandler::OnError(
		ErrorCode code, const void* location
	)
{
	/// Logs an error message if there was a decoding problem with the DDL payload.
	
	long bytepos = long(location) - long(fBufferStart) + sizeof(AliRawDataHeader);
	HLTError("There is a problem with decoding the raw data. %s (Error code: %d, at byte %d)",
		ErrorCodeToMessage(code), code, bytepos
	);
}

