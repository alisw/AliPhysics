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
		bool scalarEvent,
		AliHLTMUONTriggerRecordStruct* trigRecord,
		AliHLTUInt32_t& nofTrigRec
	)
{
	/// Runs the trigger reconstruction algorithm on the raw data.
	/// [in]  \param rawData  Pointer to the raw data DDL payload.
	/// [in]  \param rawDataSize  Size of the raw data DDL payload in bytes.
	/// [in]  \param scalarEvent  Indicates if the raw data should contain
	///      scalar data also.
	/// [out] \param trigRecord  Pointer to output buffer for reconstructed
	///      trigger records.
	/// [in/out] \param nofTrigRec  Initialy should indicate the number of
	///      elements that can be stored in the trigRecord array. It will
	///      contain the number of elements filled after this method has returned.
	/// \return true if raw data was decoded and false if there was a problem
	///     with the raw data or we overflowed the output buffer.
	///
	/// \note OverflowedOutputBuffer() can be used to check if the output
	/// buffer 'trigRecord' was overflowed during this method call.
	
	// Reset and initialise some variables in the decoder.
	fDecoder.GetHandler().MaxOutputTrigRecs(nofTrigRec);
	fDecoder.GetHandler().OutputTrigRecs(trigRecord);
	
	if (not fDecoder.Decode(rawData, rawDataSize, scalarEvent))
	{
		if (TryRecover())
		{
			HLTWarning("There was a problem with the raw data."
				" Recovered as much data as possible."
				" Will continue processing the next event."
			);
		}
		else
		{
			HLTError("Failed to decode the trigger DDL raw data.");
			return false;
		}
	}
	
	// nofTrigRec now becomes the output of how many trigger records were found.
	nofTrigRec = fDecoder.GetHandler().OutputTrigRecsCount();
	
	return not fDecoder.GetHandler().OverflowedOutputBuffer();
}


void AliHLTMUONTriggerReconstructor::TryRecover(bool value)
{
	/// Sets the flag indicating if the decoder should enable the error
	/// recovery logic.
	
	fDecoder.TryRecover(value);
	fDecoder.ExitOnError(not value);
	fDecoder.GetHandler().WarnOnly(value);
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
	fDDLBit(0),
	fSuppressPartialTriggers(false),
	fOverflowed(false),
	fWarnOnly(false),
	fUseLocalId(true),
	fUseCrateId(true),
	fCurrentCrateId(0),
	fCurrentRegional(0)
{
	/// Default constructor just resets the lookup table to zero and local
	/// structure marker pointers to NULL.
	
	for (AliHLTInt32_t i = 0; i < 16; i++)
	for (AliHLTInt32_t j = 0; j < 16; j++)
	for (AliHLTInt32_t k = 0; k < 4; k++)
	for (AliHLTInt32_t n = 0; n < 2; n++)
	for (AliHLTInt32_t m = 0; m < 16; m++)
	{
		fLookupTable.fRow[i][j][k][n][m].fIdFlags = 0x0;
		fLookupTable.fRow[i][j][k][n][m].fX = 0;
		fLookupTable.fRow[i][j][k][n][m].fY = 0;
		fLookupTable.fRow[i][j][k][n][m].fZ = 0;
	}
}


bool AliHLTMUONTriggerReconstructor::AliDecoderHandler::FindStripsOnMT1(
		const AliMUONLocalInfoStruct* localStruct,
		AliHLTInt32_t& xPos, AliHLTInt32_t& yPos
	)
{
	/// This method will find the X and Y strip positions on stations MT1
	/// of the trigger system which were fired for the corresponding
	/// L0 local trigger decision.
	/// [in]  \param localStruct  The local trigger structure as found in the DDL payload.
	/// [out] \param xPos  The X strip that was fired.
	/// [out] \param yPos  The Y strip that was fired.
	/// \return  true is returned if a strip was fired, otherwise a warning is
	///      generated and false is returned.

	// Try to identify the strips on MT1 (chambers 11 or 12) that fired
	// the trigger and set yPos and xPos to the correct values.
	// For the Y strips the yPos value might or might not have to be divided
	// by 2. This depends on the switches in the trigger electronics and how
	// they were configured. To avoid having to try to track this info we
	// just use a trial and error method.
	yPos = GetLocalYPos(localStruct);
	AliHLTUInt32_t yStrips = GetLocalY1(localStruct);
	if (((yStrips >> yPos) & 0x1) == 0x1)
	{
		// nothing to do, yPos is already correct.
	}
	else if (((yStrips >> (yPos + 1)) & 0x1) == 0x1)
	{
		yPos = yPos + 1;
	}
	else if (((yStrips >> (yPos / 2)) & 0x1) == 0x1)
	{
		yPos = yPos / 2;
	}
	else if (((yStrips >> (yPos / 2 + 1)) & 0x1) == 0x1)
	{
		yPos = yPos / 2 + 1;
	}
	else
	{
		// OK, move onto chamber 12.
		yStrips = GetLocalY2(localStruct);
		if (((yStrips >> (yPos)) & 0x1) == 0x1)
		{
			// nothing to do, yPos is already correct.
		}
		else if (((yStrips >> (yPos + 1)) & 0x1) == 0x1)
		{
			yPos = yPos + 1;
		}
		else if (((yStrips >> (yPos / 2)) & 0x1) == 0x1)
		{
			yPos = yPos / 2;
		}
		else if (((yStrips >> (yPos / 2 + 1)) & 0x1) == 0x1)
		{
			yPos = yPos / 2 + 1;
		}
		else
		{
			// At this point give up on the value of yPos and just
			// try find the first strip that was fired.
			yStrips = GetLocalY1(localStruct);
			for (AliHLTInt32_t i = 0; i < 16; i++)
			{
				if (((yStrips >> i) & 0x1) == 0x1)
				{
					yPos = i;
					goto foundYstrip;
				}
			}
			
			yStrips = GetLocalY2(localStruct);
			for (AliHLTInt32_t i = 0; i < 16; i++)
			{
				if (((yStrips >> i) & 0x1) == 0x1)
				{
					yPos = i;
					goto foundYstrip;
				}
			}
			
			// No y strip found in MT1 so this local trigger circuit
			// does not pass the 3/4 coincidence requirement,
			// so ignore it and continue.
			HLTWarning("Could not find fired Y strip for local trigger"
				" structure (regional structure = %d, crate ID = %d, ID = %d),"
				" which corresponds to triggered strip YPos = %d.",
				fCurrentRegional, fCurrentCrateId, GetLocalId(localStruct),
				GetLocalYPos(localStruct)
			);
			return false;
			
		foundYstrip: ;
		}
	}
	
	// Now find the X strip on MT1 that fired the trigger.
	xPos = GetLocalXPos(localStruct);
	AliHLTUInt32_t xStrips = GetLocalX1(localStruct);
	if (((xStrips >> (xPos / 2)) & 0x1) == 0x1)
	{
		xPos = xPos / 2;
	}
	else if (((xStrips >> (xPos / 2 + 1)) & 0x1) == 0x1)
	{
		xPos = xPos / 2 + 1;
	}
	else
	{
		// OK, move onto chamber 12.
		xStrips = GetLocalX2(localStruct);
		if (((xStrips >> (xPos / 2)) & 0x1) == 0x1)
		{
			xPos = xPos / 2;
		}
		else if (((xStrips >> (xPos / 2 + 1)) & 0x1) == 0x1)
		{
			xPos = xPos / 2 + 1;
		}
		else
		{
			// At this point give up on the value of xPos and just
			// try find the first strip that was fired.
			xStrips = GetLocalX1(localStruct);
			for (AliHLTInt32_t i = 0; i < 16; i++)
			{
				if (((xStrips >> i) & 0x1) == 0x1)
				{
					xPos = i;
					goto foundXstrip;
				}
			}
			
			xStrips = GetLocalX2(localStruct);
			for (AliHLTInt32_t i = 0; i < 16; i++)
			{
				if (((xStrips >> i) & 0x1) == 0x1)
				{
					xPos = i;
					goto foundXstrip;
				}
			}
			
			// No x strip found in MT1 so this local trigger circuit
			// does not pass the 3/4 coincidence requirement,
			// so ignore it and continue.
			HLTWarning("Could not find fired X strip for local trigger"
				" structure (regional structure = %d, crate ID = %d, ID = %d),"
				" which corresponds to triggered strip XPos = %d.",
				fCurrentRegional, fCurrentCrateId, GetLocalId(localStruct),
				GetLocalXPos(localStruct)
			);
			return false;
			
		foundXstrip: ;
		}
	}
	
	return true;
}


void AliHLTMUONTriggerReconstructor::AliDecoderHandler::FindXStrips(
		const AliMUONLocalInfoStruct* localStruct, AliHLTInt32_t startPos,
		AliHLTInt32_t pos[4]
	)
{
	/// Finds the X strips that were fired in the local trigger structure.
	/// [in] \param  localStruct  The local trigger structure as found in the DDL payload.
	/// [in] \param  startPos  The first X strip location to start looking from.
	/// [out] \param pos  Array of X strip positions on chambers 11 to 14. pos[0]
	///     is for chamber 11, pos[1] for chamber 12 and so on.
	///     The elements of the array will contain -1 if no valid strip position
	///     was found for that chamber.
	
	assert( startPos >= 0 );
	
	if (((GetLocalX1(localStruct) >> startPos) & 0x1) == 0x1)
	{
		pos[0] = startPos;
	}
	else
	{
		pos[0] = -1;
	}
	
	AliHLTUInt32_t xStrips = GetLocalX2(localStruct);
	if (GetLocalSXDev(localStruct)) // check the direction of the deviation.
	{
		// For hits on chamber 12 we have to look for fired strips
		// within 1 strip of startPos. Where startPos is the X position
		// as found by FindStripsOnMT1.
		if (((xStrips >> startPos) & 0x1) == 0x1)
		{
			pos[1] = startPos;
		}
		else if (((xStrips >> (startPos + 1)) & 0x1) == 0x1)
		{
			pos[1] = startPos + 1;
		}
		else
		{
			pos[1] = -1;
		}
		
		// Given the MT1 coordinate and the deviation information we can
		// identify the X strip on MT2 that corresponds to the L0 trigger.
		// For fired strips on MT2 we look for strips that are within 2
		// strips of the position endPos = (posX + deviation) / 2, where
		// posX = GetLocalXPos(localStruct);
		// deviation = GetLocalXDev(localStruct)
		// We use the 2 strip tollerance because there is a truncation of 1
		// bit when we apply integer divide by 2.
		// This procedure should thus apply the same constraints and
		// tollerance as the L0 electronics for the X strip 'mini-roads',
		// ref. section 3.4.2.2, "Technical Design Report of the Dimuon
		// Forward Spectrometer".
		AliHLTInt32_t endPos = (GetLocalXPos(localStruct) + GetLocalXDev(localStruct)) / 2;
		
		// Note the order of the checks are such that we choose the strip with
		// giving the smallest deviation.
		xStrips = GetLocalX3(localStruct);
		if (endPos >= 0 and ((xStrips >> endPos) & 0x1) == 0x1)
		{
			pos[2] = endPos;
		}
		else if (endPos - 1 >= 0 and ((xStrips >> (endPos - 1)) & 0x1) == 0x1)
		{
			pos[2] = endPos - 1;
		}
		else if (endPos + 1 >= 0 and ((xStrips >> (endPos + 1)) & 0x1) == 0x1)
		{
			pos[2] = endPos + 1;
		}
		else if (endPos - 2 >= 0 and ((xStrips >> (endPos - 2)) & 0x1) == 0x1)
		{
			pos[2] = endPos - 2;
		}
		else if (endPos + 2 >= 0 and ((xStrips >> (endPos + 2)) & 0x1) == 0x1)
		{
			pos[2] = endPos + 2;
		}
		else
		{
			pos[2] = -1;
		}
		
		xStrips = GetLocalX4(localStruct);
		if (endPos >= 0 and ((xStrips >> endPos) & 0x1) == 0x1)
		{
			pos[3] = endPos;
		}
		else if (endPos - 1 >= 0 and ((xStrips >> (endPos - 1)) & 0x1) == 0x1)
		{
			pos[3] = endPos - 1;
		}
		else if (endPos + 1 >= 0 and ((xStrips >> (endPos + 1)) & 0x1) == 0x1)
		{
			pos[3] = endPos + 1;
		}
		else if (endPos - 2 >= 0 and ((xStrips >> (endPos - 2)) & 0x1) == 0x1)
		{
			pos[3] = endPos - 2;
		}
		else if (endPos + 2 >= 0 and ((xStrips >> (endPos + 2)) & 0x1) == 0x1)
		{
			pos[3] = endPos + 2;
		}
		else
		{
			pos[3] = -1;
		}
	}
	else
	{
		// The following code is the same as for the
		// GetLocalSXDev(localStruct) == true case above, but with the
		// arithmetic changing sign.
		if (((xStrips >> startPos) & 0x1) == 0x1)
		{
			pos[1] = startPos;
		}
		else if (startPos - 1 >= 0 and ((xStrips >> (startPos - 1)) & 0x1) == 0x1)
		{
			pos[1] = startPos + 1;
		}
		else
		{
			pos[1] = -1;
		}
		
		AliHLTInt32_t endPos = (GetLocalXPos(localStruct) - GetLocalXDev(localStruct)) / 2;
		
		xStrips = GetLocalX3(localStruct);
		if (endPos >= 0 and ((xStrips >> endPos) & 0x1) == 0x1)
		{
			pos[2] = endPos;
		}
		else if (endPos + 1 >= 0 and ((xStrips >> (endPos + 1)) & 0x1) == 0x1)
		{
			pos[2] = endPos + 1;
		}
		else if (endPos - 1 >= 0 and ((xStrips >> (endPos - 1)) & 0x1) == 0x1)
		{
			pos[2] = endPos - 1;
		}
		else if (endPos + 2 >= 0 and ((xStrips >> (endPos + 2)) & 0x1) == 0x1)
		{
			pos[2] = endPos + 2;
		}
		else if (endPos - 2 >= 0 and ((xStrips >> (endPos - 2)) & 0x1) == 0x1)
		{
			pos[2] = endPos - 2;
		}
		else
		{
			pos[2] = -1;
		}
		
		xStrips = GetLocalX4(localStruct);
		if (endPos >= 0 and ((xStrips >> endPos) & 0x1) == 0x1)
		{
			pos[3] = endPos;
		}
		else if (endPos + 1 >= 0 and ((xStrips >> (endPos + 1)) & 0x1) == 0x1)
		{
			pos[3] = endPos + 1;
		}
		else if (endPos - 1 >= 0 and ((xStrips >> (endPos - 1)) & 0x1) == 0x1)
		{
			pos[3] = endPos - 1;
		}
		else if (endPos + 2 >= 0 and ((xStrips >> (endPos + 2)) & 0x1) == 0x1)
		{
			pos[3] = endPos + 2;
		}
		else if (endPos - 2 >= 0 and ((xStrips >> (endPos - 2)) & 0x1) == 0x1)
		{
			pos[3] = endPos - 2;
		}
		else
		{
			pos[3] = -1;
		}
	}
}


void AliHLTMUONTriggerReconstructor::AliDecoderHandler::FindYStrips(
		const AliMUONLocalInfoStruct* localStruct, AliHLTInt32_t startPos,
		AliHLTInt32_t pos[4]
	)
{
	/// Finds the Y strips that were fired in the local trigger structure.
	/// [in] \param  localStruct  The local trigger structure as found in the DDL payload.
	/// [in] \param  startPos  The first Y strip location to start looking from.
	/// [out] \param pos  Array of Y strip positions on chambers 11 to 14. pos[0]
	///     is for chamber 11, pos[1] for chamber 12 and so on.
	///     The elements of the array will contain -1 if no valid strip position
	///     was found for that chamber.
	
	assert( startPos >= 0 );
	
	// First we scan from the i'th = startPos strip upwards (i.e. i+1, i+2 etc..)
	// to find the first fired strip. Then we similarly scan downwards
	// (i.e. i-1, i-2 etc..) to find the first fired strip. We actually only need
	// to check in the range [i-1 .. i+1] due to the constraint that valid tracks
	// only have a +/- 1 Y strip deviation on consecutive chambers.
	// Ideally we should have all of posUp[i] == posDown[i] == startPos, but this
	// need not be the case due to multiple scattering.
	// This procedure should thus apply the same constraints and tollerance
	// as the L0 electronics for the Y strip 'roads',
	// ref. section 3.4.2.2, "Technical Design Report of the Dimuon Forward Spectrometer".
	AliHLTUInt32_t strips[4] = {
			GetLocalY1(localStruct), GetLocalY2(localStruct),
			GetLocalY3(localStruct), GetLocalY4(localStruct)
		};
	AliHLTUInt8_t posUpCount = 0, posDownCount = 0;
	AliHLTInt32_t posUp[4] = {-1, -1, -1, -1};
	AliHLTInt32_t posDown[4] = {-1, -1, -1, -1};
	for (AliHLTInt32_t n = 0; n < 4; n++)
	{
		for (AliHLTInt32_t i = startPos; i <= startPos+1; i++)
		{
			if (((strips[n] >> i) & 0x1) == 0x1)
			{
				posUp[n] = i;
				posUpCount++;
				break;
			}
		}
		for (AliHLTInt32_t i = startPos; i >= 0; i--)
		{
			if (((strips[n] >> i) & 0x1) == 0x1)
			{
				posDown[n] = i;
				posDownCount++;
				break;
			}
		}
	}
	
	// Now select either posUp or posDown, whichever has the most found strips.
	if (posUpCount > posDownCount)
	{
		for (AliHLTInt32_t n = 0; n < 4; n++)
			pos[n] = posUp[n];
	}
	else
	{
		for (AliHLTInt32_t n = 0; n < 4; n++)
			pos[n] = posDown[n];
	}
}


void AliHLTMUONTriggerReconstructor::AliDecoderHandler::ReconstructHit(
		AliHLTUInt32_t xStrips, AliHLTUInt32_t yStrips,
		AliHLTInt32_t xPos, AliHLTInt32_t yPos,
		AliHLTUInt8_t crateId, AliHLTUInt8_t locId, AliHLTUInt8_t chamber,
		AliHLTMUONRecHitStruct& hit
	)
{
	/// Reconstructs the hit coordinates for the given chamber from the
	/// strip and fired strip information provided.
	/// [in]  \param xStrips  The X strip pattern for the given chamber.
	/// [in]  \param yStrips  The Y strip pattern for the given chamber.
	/// [in]  \param xPos  The position of the X strip that was fired.
	/// [in]  \param yPos  The position of the Y strip that was fired.
	/// [in]  \param chamber  The chamber on which the strips were found.
	///      Valid range [0..3].
	/// [out] \param hit  This will be filled with the reconstructed hit.

	assert( 0 <= xPos and xPos < 16 );
	assert( 0 <= yPos and yPos < 16 );
	assert( ((xStrips >> xPos) & 0x1) == 0x1 );
	assert( ((yStrips >> yPos) & 0x1) == 0x1 );
	assert( chamber <= 3 );
	assert( crateId < 16 );
	assert( locId < 16 );
	
	// Decode the Y position of the hit from the strip position information.
	// If two neighbouring strips were fired then we decluster the strips like
	// the L0 electronics does by taking the middle position of the two strips.
	if (xPos > 0 and ((xStrips >> (xPos-1)) & 0x1) == 0x1)
	{
		if (((xStrips >> (xPos+1)) & 0x1) == 0x1)
		{
			// Strips fired on both sides of strip at xPos so just use the middle one.
			hit.fFlags = fLookupTable.fRow[crateId][locId][chamber][0][xPos].fIdFlags;
			hit.fY = fLookupTable.fRow[crateId][locId][chamber][0][xPos].fY;
			hit.fZ = fLookupTable.fRow[crateId][locId][chamber][0][xPos].fZ;
		}
		else
		{
			// Second strip fired below the one at xPos, so decluster.
			assert(xPos-1 < 16);
			hit.fFlags = fLookupTable.fRow[crateId][locId][chamber][0][xPos].fIdFlags;
			hit.fY = (fLookupTable.fRow[crateId][locId][chamber][0][xPos-1].fY
				+ fLookupTable.fRow[crateId][locId][chamber][0][xPos].fY) * 0.5;
			hit.fZ = (fLookupTable.fRow[crateId][locId][chamber][0][xPos-1].fZ
				+ fLookupTable.fRow[crateId][locId][chamber][0][xPos].fZ) * 0.5;
		}
	}
	else
	{
		if (((xStrips >> (xPos+1)) & 0x1) == 0x1)
		{
			// Second strip fired above the one at xPos, so decluster.
			assert(xPos+1 < 16);
			hit.fFlags = fLookupTable.fRow[crateId][locId][chamber][0][xPos].fIdFlags;
			hit.fY = (fLookupTable.fRow[crateId][locId][chamber][0][xPos+1].fY
				+ fLookupTable.fRow[crateId][locId][chamber][0][xPos].fY) * 0.5;
			hit.fZ = (fLookupTable.fRow[crateId][locId][chamber][0][xPos+1].fZ
				+ fLookupTable.fRow[crateId][locId][chamber][0][xPos].fZ) * 0.5;
		}
		else
		{
			// Only strip at xPos fired and neither of its two neighbours.
			hit.fFlags = fLookupTable.fRow[crateId][locId][chamber][0][xPos].fIdFlags;
			hit.fY = fLookupTable.fRow[crateId][locId][chamber][0][xPos].fY;
			hit.fZ = fLookupTable.fRow[crateId][locId][chamber][0][xPos].fZ;
		}
	}
	
	// Similarly decode the X position of the hit from the strip position information.
	// Also decluster like for the Y strips.
	if (yPos > 0 and ((yStrips >> (yPos-1)) & 0x1) == 0x1)
	{
		if (((yStrips >> (yPos+1)) & 0x1) == 0x1)
		{
			// Strips fired on both sides of strip at yPos so just use the middle one.
			hit.fX = fLookupTable.fRow[crateId][locId][chamber][1][yPos].fX;
		}
		else
		{
			// Second strip fired below the one at yPos, so decluster.
			assert(yPos-1 < 16);
			hit.fX = (fLookupTable.fRow[crateId][locId][chamber][1][yPos-1].fX
				+ fLookupTable.fRow[crateId][locId][chamber][1][yPos].fX) * 0.5;
		}
	}
	else
	{
		if (((yStrips >> (yPos+1)) & 0x1) == 0x1)
		{
			// Second strip fired above the one at yPos, so decluster.
			assert(yPos+1 < 16);
			hit.fX = (fLookupTable.fRow[crateId][locId][chamber][1][yPos+1].fX
				+ fLookupTable.fRow[crateId][locId][chamber][1][yPos].fX) * 0.5;
		}
		else
		{
			// Only strip at yPos fired and neither of its two neighbours.
			hit.fX = fLookupTable.fRow[crateId][locId][chamber][1][yPos].fX;
		}
	}
}


void AliHLTMUONTriggerReconstructor::AliDecoderHandler::OnLocalStructV2(
		UInt_t iloc,
		const AliMUONLocalInfoStruct* localStruct,
		const AliMUONLocalScalarsStruct* /*scalars*/
	)
{
	/// Converts a local trigger structure from the L0 into a trigger record.
	/// The dHLT trigger records is then used as a seed for tracking algorithms.
	/// \note fOutputTrigRecs must be set before calling the decoder to decode
	///    a new raw data buffer.
	/// \param localStruct  This is a pointer to the local L0 trigger structure data.

	assert(iloc < 16);
	assert(localStruct != NULL);
	assert(fOutputTrigRecs != NULL);
	
	// We must have at least one bit in each of the 4 strip words otherwise
	// if one of the words is zero it means we only have X or Y coordinate
	// information for a station or no coordinate on one of the stations at all.
	if (localStruct->fX2X1 == 0 or localStruct->fX4X3 == 0 or
	    localStruct->fY2Y1 == 0 or localStruct->fY4Y3 == 0
	   )
	{
		return;
	}
	
	// First try to identify the X and Y strips on MT1 that fired the trigger.
	// Note: X strips are for the Y axis in ALICE coordinate system,
	// i.e. bending plane. and Y strips for the X axis.
	AliHLTInt32_t xPos, yPos;
	if (not FindStripsOnMT1(localStruct, xPos, yPos)) return;
	
	// Check that we will not overflow the output buffer.
	if (fOutputTrigRecsCount >= fMaxOutputTrigRecs)
	{
		HLTError("Output buffer has overflowed maximum element count of %d.",
			fMaxOutputTrigRecs
		);
		fOverflowed = true;
		return;
	}
	
	AliHLTMUONTriggerRecordStruct& trigger = fOutputTrigRecs[fOutputTrigRecsCount];
	
	// Now try find all the X and Y strips on all 4 trigger chambers.
	AliHLTInt32_t stripPosX[4];
	FindXStrips(localStruct, xPos, stripPosX);
	AliHLTInt32_t stripPosY[4];
	FindYStrips(localStruct, yPos, stripPosY);
	
	// hitset indicates which hits on chambers 7 to 10 have been found and filled.
	bool hitset[4] = {false, false, false, false};
	
	// Reconstruct the hits from the found strips. Also, fill the hitset
	// flags and make sure the hits for which no valid strips were found get
	// set to a nil value.
	AliHLTUInt32_t xStrips[4] = {
			GetLocalX1(localStruct), GetLocalX2(localStruct),
			GetLocalX3(localStruct), GetLocalX4(localStruct)
		};
	AliHLTUInt32_t yStrips[4] = {
			GetLocalY1(localStruct), GetLocalY2(localStruct),
			GetLocalY3(localStruct), GetLocalY4(localStruct)
		};
	int hitCount = 0;
	AliHLTUInt8_t locId = GetLocalId(localStruct);
	for (int i = 0; i < 4; i++)
	{
		if (stripPosX[i] != -1 and stripPosY[i] != -1)
		{
			ReconstructHit(
					xStrips[i], yStrips[i],
					stripPosX[i], stripPosY[i],
					fCurrentCrateId, locId, i,
					trigger.fHit[i]
				);
			hitset[i] = true;
			hitCount++;
		}
		else
		{
			trigger.fHit[i] = AliHLTMUONConstants::NilRecHitStruct();
			hitset[i] = false;
		}
	}

	// Construct the ID from the running counter fTrigRecId and use the
	// regional counter, local counter and DDL id for the bottom 8 bits.
	trigger.fId = (fTrigRecId << 8) | fDDLBit | ((fCurrentRegional & 0x7) << 4) | (iloc & 0xF);

	// Increment the trigger record ID and warp it around at 0x7FFFFF since
	// the bottom 8 bits are filled with the regional + local counters and the
	// sign bit in fOutputTrigRecs[fOutputTrigRecsCount].fId must be positive.
	fTrigRecId = (fTrigRecId + 1) & 0x007FFFFF;
	
	// Set the ideal Z coordinate used in line fit for trigger record to
	// the same as the Z coordinate for the hits that were found, otherwise
	// use nominal coordinates.
	AliHLTFloat32_t chamberZ11 = -1603.5f;
	if (hitset[0]) chamberZ11 = trigger.fHit[0].fZ;
	if (hitset[1]) chamberZ11 = trigger.fHit[1].fZ;
	AliHLTFloat32_t chamberZ13 = -1703.5f;
	if (hitset[2]) chamberZ13 = trigger.fHit[2].fZ;
	if (hitset[3]) chamberZ13 = trigger.fHit[3].fZ;
	AliHLTMUONCalculations::IdealZ1(chamberZ11);
	AliHLTMUONCalculations::IdealZ2(chamberZ13);
	
	if (hitCount >= 3 and
	    AliHLTMUONCalculations::FitLineToTriggerRecord(trigger, hitset)
	   )
	{
		// Calculate the momentum and fill in the flags and momentum fields.
		AliHLTMUONCalculations::ComputeMomentum(
				AliHLTMUONCalculations::IdealX1(),
				AliHLTMUONCalculations::IdealY1(),
				AliHLTMUONCalculations::IdealY2(),
				AliHLTMUONCalculations::IdealZ1(),
				AliHLTMUONCalculations::IdealZ2()
			);
		
		trigger.fPx = AliHLTMUONCalculations::Px();
		trigger.fPy = AliHLTMUONCalculations::Py();
		trigger.fPz = AliHLTMUONCalculations::Pz();

		trigger.fFlags = AliHLTMUONUtils::PackTriggerRecordFlags(
				AliHLTMUONCalculations::Sign(),
				hitset
			);
		
		fOutputTrigRecsCount++;
	}
	else if ((hitset[0] or hitset[1] or hitset[2] or hitset[3])
	         and not fSuppressPartialTriggers
	        )
	{
		trigger.fPx = trigger.fPy = trigger.fPz = 0;
		
		trigger.fFlags = AliHLTMUONUtils::PackTriggerRecordFlags(
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
	if (fWarnOnly)
	{
		HLTWarning("There is a problem with decoding the raw data."
			" %s (Error code: %d, at byte %d). Trying to recover from corrupt data.",
			ErrorCodeToMessage(code), code, bytepos
		);
	}
	else
	{
		HLTError("There is a problem with decoding the raw data. %s (Error code: %d, at byte %d)",
			ErrorCodeToMessage(code), code, bytepos
		);
	}
}

