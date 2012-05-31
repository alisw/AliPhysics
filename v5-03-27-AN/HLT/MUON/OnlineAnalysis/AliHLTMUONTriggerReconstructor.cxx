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

// $Id$

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
#include "AliHLTMUONTrigRecsDebugBlockStruct.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONCalculations.h"
#include "AliMUONConstants.h"
#include "AliRawDataHeader.h"
#include <vector>
#include <cassert>


const AliMUONLocalInfoStruct AliHLTMUONTriggerReconstructor::AliDecoderHandler::fgkNullStruct =
{
	0x0, 0x0, 0x0, 0x0, 0x0
};


AliHLTMUONTriggerReconstructor::AliHLTMUONTriggerReconstructor() :
	AliHLTLogging(),
	fDecoder()
{
	/// Default constructor.
	
	fDecoder.MaxRegionals(8);
	fDecoder.MaxLocals(16);
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
	/// \param [in] rawData  Pointer to the raw data DDL payload.
	/// \param [in] rawDataSize  Size of the raw data DDL payload in bytes.
	/// \param [in] scalarEvent  Indicates if the raw data should contain
	///      scalar data also.
	/// \param [out] trigRecord  Pointer to output buffer for reconstructed
	///      trigger records.
	/// \param [in,out] nofTrigRec  Initialy should indicate the number of
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
			/// Fix as long as the DARC header problem is not fixed in hardware by trigger colleagues
			if (fDecoder.GetHandler().HadNonWrongEventTypeError() or
			    (fDecoder.GetHandler().HadWrongEventTypeError() and not fDecoder.GetHandler().DontPrintWrongEventError())
			   )
			{
				HLTWarning("There was a problem with the raw data."
					" Recovered as much data as possible."
					" Will continue processing the next event."
				);
			}
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
	fCurrentRegional(0),
	fNextLocalIndex(0),
	fPrevStruct(&fgkNullStruct),
	fCurrentStruct(&fgkNullStruct),
	fNextStruct(&fgkNullStruct),
	fStoreInfo(false),
	fInfoBufferSize(0),
	fInfoBufferCount(0),
	fInfoBuffer(NULL),
	fDontPrintWrongEventError(false),
	fHadWrongEventTypeError(false),
	fHadNonWrongEventTypeError(false)
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


AliHLTMUONTriggerReconstructor::AliDecoderHandler::~AliDecoderHandler()
{
	// Default destructor deletes allocated array.
	
	if (fInfoBuffer != NULL) delete [] fInfoBuffer;
}


bool AliHLTMUONTriggerReconstructor::AliDecoderHandler::FindStripsOnMT1(
		AliHLTInt32_t& xPos, AliHLTInt32_t& yPos
	)
{
	/// This method will find the X and Y strip positions on stations MT1 of the
	/// trigger system which were fired for the current L0 local trigger decision.
	/// \param [out] xPos  The X strip that was fired.
	/// \param [out] yPos  The Y strip that was fired.
	/// \return  true is returned if a strip was fired, otherwise a warning is
	///      generated and false is returned.
	/// \note Values for xPos and yPos are in the range [0..15].

	// Try to identify the strips on MT1 (chambers 11 or 12) that fired
	// the trigger and set yPos and xPos to the correct values.
	// For the Y strips the yPos value might or might not have to be divided
	// by 2. This depends on the switches in the trigger electronics and how
	// they were configured. To avoid having to try to track this info we
	// just use a trial and error method.
	yPos = GetLocalYPos(fCurrentStruct);
	AliHLTUInt32_t yStrips1 = GetLocalY1(fCurrentStruct);
	AliHLTUInt32_t yStrips2 = GetLocalY2(fCurrentStruct);
	if (((yStrips1 >> yPos) & 0x1) != 0x1 and ((yStrips2 >> yPos) & 0x1) != 0x1)
	{
		if (((yStrips1 >> (yPos / 2)) & 0x1) == 0x1 or ((yStrips2 >> (yPos / 2)) & 0x1) == 0x1)
		{
			yPos = yPos / 2;
		}
		else
		{
			for (AliHLTInt32_t i = 1; i < 16; ++i)
			{
				if (yPos + i < 16 and (((yStrips1 >> (yPos + i)) & 0x1) == 0x1 or
				                       ((yStrips2 >> (yPos + i)) & 0x1) == 0x1)
				   )
				{
					yPos = yPos + i;
					break;
				}
				else if (yPos / 2 + i < 16 and (((yStrips1 >> (yPos / 2 + i)) & 0x1) == 0x1 or
				                                ((yStrips2 >> (yPos / 2 + i)) & 0x1) == 0x1)
				        )
				{
					yPos = yPos / 2 + i;
					break;
				}
				else if (yPos - i >= 0 and (((yStrips1 >> (yPos - i)) & 0x1) == 0x1 or
				                            ((yStrips2 >> (yPos - i)) & 0x1) == 0x1)
				        )
				{
					yPos = yPos - i;
					break;
				}
				else if (yPos / 2 - i >= 0 and (((yStrips1 >> (yPos / 2 - i)) & 0x1) == 0x1 or
				                                ((yStrips2 >> (yPos / 2 - i)) & 0x1) == 0x1)
				        )
				{
					yPos = yPos / 2 - i;
					break;
				}
			}
			if (((yStrips1 >> yPos) & 0x1) != 0x1 and ((yStrips2 >> yPos) & 0x1) != 0x1)
			{
				// No y strip found in MT1 so this local trigger circuit does not
				// pass the 3/4 coincidence requirement, so ignore it and continue.
				HLTWarning("Could not find fired Y strip for local trigger"
					" structure (regional structure = %d, crate ID = %d, ID = %d),"
					" which corresponds to triggered strip YPos = %d.",
					fCurrentRegional, fCurrentCrateId, GetLocalId(fCurrentStruct),
					GetLocalYPos(fCurrentStruct)
				);
				return false;
			}
		}
	}
	
	// Now find the X strip on MT1 that fired the trigger.
	xPos = GetLocalXPos(fCurrentStruct) / 2;
	AliHLTUInt32_t xStrips1 = GetLocalX1(fCurrentStruct);
	AliHLTUInt32_t xStrips2 = GetLocalX2(fCurrentStruct);
	if (((xStrips1 >> xPos) & 0x1) != 0x1 and ((xStrips2 >> xPos) & 0x1) != 0x1)
	{
		for (AliHLTInt32_t i = 1; i < 16; ++i)
		{
			if (xPos + i < 16 and (((xStrips1 >> (xPos + i)) & 0x1) == 0x1 or
			                       ((xStrips2 >> (xPos + i)) & 0x1) == 0x1)
			   )
			{
				xPos = xPos + i;
				break;
			}
			else if (xPos - i >= 0 and (((xStrips1 >> (xPos - i)) & 0x1) == 0x1 or
			                            ((xStrips2 >> (xPos - i)) & 0x1) == 0x1)
			        )
			{
				xPos = xPos - i;
				break;
			}
		}
		if (((xStrips1 >> xPos) & 0x1) != 0x1 and ((xStrips2 >> xPos) & 0x1) != 0x1)
		{
			// No x strip found in MT1 so this local trigger circuit does not
			// pass the 3/4 coincidence requirement, so ignore it and continue.
			HLTWarning("Could not find fired X strip for local trigger"
				" structure (regional structure = %d, crate ID = %d, ID = %d),"
				" which corresponds to triggered strip XPos = %d.",
				fCurrentRegional, fCurrentCrateId, GetLocalId(fCurrentStruct),
				GetLocalXPos(fCurrentStruct)
			);
			return false;
		}
	}
	
	return true;
}


void AliHLTMUONTriggerReconstructor::AliDecoderHandler::SelectXPatterns(
		AliHLTUInt64_t strips[4]
	)
{
	// Select the correct X strips to use.
	
	assert( fCurrentStruct != NULL );
	
	strips[0] = AliHLTUInt64_t(GetLocalX1(fPrevStruct)) |
		(AliHLTUInt64_t(GetLocalX1(fCurrentStruct)) << 16) |
		(AliHLTUInt64_t(GetLocalX1(fNextStruct)) << 32);
	
	strips[1] = AliHLTUInt64_t(GetLocalX2(fPrevStruct)) |
		(AliHLTUInt64_t(GetLocalX2(fCurrentStruct)) << 16) |
		(AliHLTUInt64_t(GetLocalX2(fNextStruct)) << 32);
		
	strips[2] = AliHLTUInt64_t(GetLocalX3(fPrevStruct)) |
		(AliHLTUInt64_t(GetLocalX3(fCurrentStruct)) << 16) |
		(AliHLTUInt64_t(GetLocalX3(fNextStruct)) << 32);
		
	strips[3] = AliHLTUInt64_t(GetLocalX4(fPrevStruct)) |
		(AliHLTUInt64_t(GetLocalX4(fCurrentStruct)) << 16) |
		(AliHLTUInt64_t(GetLocalX4(fNextStruct)) << 32);
}


void AliHLTMUONTriggerReconstructor::AliDecoderHandler::SelectYPatterns(
		AliHLTInt32_t xpos[4], AliHLTUInt32_t strips[4], AliHLTUInt8_t locId[4]
	)
{
	// Select the correct Y strip patterns and local IDs based on the X strip positions found.
	
	AliHLTUInt8_t locIdnext = fUseLocalId ? GetLocalId(fNextStruct) : AliHLTUInt8_t(fNextLocalIndex);
	if (locIdnext >= 16) locIdnext = 0;
	AliHLTUInt8_t locIdcurr = fUseLocalId ? GetLocalId(fCurrentStruct) : AliHLTUInt8_t(fNextLocalIndex-1);
	if (locIdcurr >= 16) locIdcurr = 0;
	AliHLTUInt8_t locIdprev = fUseLocalId ? GetLocalId(fPrevStruct) : AliHLTUInt8_t(fNextLocalIndex-2);
	if (locIdprev >= 16) locIdprev = 0;
	
	UShort_t patterns[4][3] = {
		{GetLocalY1(fPrevStruct), GetLocalY1(fCurrentStruct), GetLocalY1(fNextStruct)},
		{GetLocalY2(fPrevStruct), GetLocalY2(fCurrentStruct), GetLocalY2(fNextStruct)},
		{GetLocalY3(fPrevStruct), GetLocalY3(fCurrentStruct), GetLocalY3(fNextStruct)},
		{GetLocalY4(fPrevStruct), GetLocalY4(fCurrentStruct), GetLocalY4(fNextStruct)}
	};
	
	for (int i = 0; i < 4; i++)
	{
		// Check if the Y strips for the different local structures are the
		// same physical strip. If they are then performs a bit or between the
		// patterns. This is necessary because the signal sometimes does not
		// propagate in time (in particular for cosmic runs). This has to do with
		// the calibration of the timings in the trigger electronics. The solution
		// here is to perform the bitwise or to form the correct strip pattern.
		UShort_t mergedPattern[3] = {patterns[i][0], patterns[i][1], patterns[i][2]};
		const AliHLTMUONTriggerRecoLutRow& lutnext = fLookupTable.fRow[fCurrentCrateId][locIdnext][i][1][0];
		const AliHLTMUONTriggerRecoLutRow& lutcurr = fLookupTable.fRow[fCurrentCrateId][locIdcurr][i][1][0];
		const AliHLTMUONTriggerRecoLutRow& lutprev = fLookupTable.fRow[fCurrentCrateId][locIdprev][i][1][0];
		if (lutprev.fX == lutcurr.fX and lutprev.fY == lutcurr.fY and lutprev.fZ == lutcurr.fZ)
		{
			mergedPattern[0] |= patterns[i][1];
			mergedPattern[1] |= patterns[i][0];
		}
		if (lutnext.fX == lutcurr.fX and lutnext.fY == lutcurr.fY and lutnext.fZ == lutcurr.fZ)
		{
			mergedPattern[1] |= patterns[i][2];
			mergedPattern[2] |= patterns[i][1];
		}
	
		if (xpos[i] >= 32)
		{
			strips[i] = mergedPattern[2];
			locId[i] = locIdnext;
		}
		else if (xpos[i] >= 16)
		{
			strips[i] = mergedPattern[1];
			locId[i] = locIdcurr;
		}
		else  if (xpos[i] >= 0)
		{
			strips[i] = mergedPattern[0];
			locId[i] = locIdprev;
		}
		else
		{
			// If the X strip could not be found then just look on the
			// current local board strips.
			strips[i] = mergedPattern[1];
			locId[i] = locIdcurr;
		}
	}
}


void AliHLTMUONTriggerReconstructor::AliDecoderHandler::FindXStrips(
		AliHLTInt32_t startPos, AliHLTUInt64_t strips[4], AliHLTInt32_t pos[4]
	)
{
	/// Finds the X strips that were fired in the local trigger structures.
	/// \param [in] startPos  The first X strip location to start looking from.
	/// \param [in] strips  The X strip patterns for chambers 11 to 14 to use.
	/// \param [out] pos  Array of X strip positions on chambers 11 to 14. pos[0]
	///     is for chamber 11, pos[1] for chamber 12 and so on.
	///     The elements of the array will contain -1 if no valid strip position
	///     was found for that chamber.
	/// \note Values for startPos and pos are in the range [0..47], where 0..15 is
	///     for strip positions in the fPrevStruct patterns, 16..31 for fCurrentStruct
	///     and 32..47 for fNextStruct.
	
	assert( startPos >= 0 );
	assert( fCurrentStruct != NULL );
	
	if (GetLocalSXDev(fCurrentStruct)) // check the direction of the deviation.
	{
		for (int i = 0; i < 2; ++i)
		{
			if (((strips[i] >> startPos) & 0x1) == 0x1)
			{
				pos[i] = startPos;
				continue;
			}
			for (AliHLTInt32_t j = 1; j < 16; ++j)
			{
				// We first check the straighter tracklet option, i.e higher momentum.
				if (startPos - j >= 0 and ((strips[i] >> (startPos - j)) & 0x1) == 0x1)
				{
					pos[i] = startPos - j;
					break;
				}
				else if (startPos + j < 48 and ((strips[i] >> (startPos + j)) & 0x1) == 0x1)
				{
					pos[i] = startPos + j;
					break;
				}
			}
			if (((strips[i] >> pos[i]) & 0x1) != 0x1) pos[i] = -1;
		}
		
		// Given the MT1 coordinate 'startPos' and the deviation information we can
		// identify the X strip on MT2 that corresponds to the L0 trigger.
		// We find fired strips on MT2 by looking for strips around the position
		// endPos = (posX + deviation) / 2, where posX = GetLocalXPos(fCurrentStruct);
		// deviation = GetLocalXDev(fCurrentStruct)
		AliHLTInt32_t endPos = (GetLocalXPos(fCurrentStruct) + GetLocalXDev(fCurrentStruct)) / 2;
		endPos += 16; // fCurrentStruct strips start at bit 16.
		
		for (int i = 2; i < 4; ++i)
		{
			if (((strips[i] >> endPos) & 0x1) == 0x1)
			{
				pos[i] = endPos;
				continue;
			}
			for (AliHLTInt32_t j = 1; j < 16; ++j)
			{
				if (endPos - j >= 0 and ((strips[i] >> (endPos - j)) & 0x1) == 0x1)
				{
					pos[i] = endPos - j;
					break;
				}
				else if (endPos + j < 48 and ((strips[i] >> (endPos + j)) & 0x1) == 0x1)
				{
					pos[i] = endPos + j;
					break;
				}
			}
			if (((strips[i] >> pos[i]) & 0x1) != 0x1) pos[i] = -1;
		}
	}
	else
	{
		// Similar logic to the positive deviation case above, but with the
		// arithmetic inversed.
		for (int i = 0; i < 2; ++i)
		{
			if (((strips[i] >> startPos) & 0x1) == 0x1)
			{
				pos[i] = startPos;
				continue;
			}
			for (AliHLTInt32_t j = 1; j < 16; ++j)
			{
				// We first check the straighter tracklet option, i.e higher momentum.
				if (startPos + j < 48 and ((strips[i] >> (startPos + j)) & 0x1) == 0x1)
				{
					pos[i] = startPos + j;
					break;
				}
				else if (startPos - j >= 0 and ((strips[i] >> (startPos - j)) & 0x1) == 0x1)
				{
					pos[i] = startPos - j;
					break;
				}
			}
			if (((strips[i] >> pos[i]) & 0x1) != 0x1) pos[i] = -1;
		}
		
		AliHLTInt32_t endPos = (GetLocalXPos(fCurrentStruct) - GetLocalXDev(fCurrentStruct)) / 2;
		endPos += 16; // fCurrentStruct strips start at bit 16.
		
		for (int i = 2; i < 4; ++i)
		{
			if (((strips[i] >> endPos) & 0x1) == 0x1)
			{
				pos[i] = endPos;
				continue;
			}
			for (AliHLTInt32_t j = 1; j < 16; ++j)
			{
				if (endPos + j < 48 and ((strips[i] >> (endPos + j)) & 0x1) == 0x1)
				{
					pos[i] = endPos + j;
					break;
				}
				else if (endPos - j >= 0 and ((strips[i] >> (endPos - j)) & 0x1) == 0x1)
				{
					pos[i] = endPos - j;
					break;
				}
			}
			if (((strips[i] >> pos[i]) & 0x1) != 0x1) pos[i] = -1;
		}
	}
}


void AliHLTMUONTriggerReconstructor::AliDecoderHandler::FindYStrips(
		AliHLTInt32_t startPos, AliHLTUInt32_t strips[4], AliHLTInt32_t pos[4]
	)
{
	/// Finds the Y strips that were fired in the local trigger structures.
	/// \param [in] startPos  The first Y strip location to start looking from.
	/// \param [in] strips  Array of Y strip patterns to look in for chamber 11 to 14.
	/// \param [out] pos  Array of Y strip positions on chambers 11 to 14. pos[0]
	///     is for chamber 11, pos[1] for chamber 12 and so on.
	///     The elements of the array will contain -1 if no valid strip position
	///     was found for that chamber.
	/// \note Values for startPos and pos are in the range [0..15].
	
	assert( startPos >= 0 );
	
	// First we scan from the i'th = startPos strip upwards (i.e. i-1, i, i+1, i+2 etc..)
	// to find the first fired strip. Then we similarly scan downwards
	// (i.e. i+1, i, i-1, i-2 etc..) to find the first fired strip.
	// Ideally we should have all of posUp[i] == posDown[i] == startPos, but this
	// need not be the case due to multiple scattering or if dealing with cosmic tracks.
	AliHLTUInt8_t posUpCount = 0, posDownCount = 0, posNearestCount = 0;
	AliHLTInt32_t posUp[4] = {-1, -1, -1, -1};
	AliHLTInt32_t posDown[4] = {-1, -1, -1, -1};
	AliHLTInt32_t posNearest[4] = {-1, -1, -1, -1};
	for (AliHLTInt32_t n = 0; n < 4; n++)
	{
		for (AliHLTInt32_t i = (startPos >= 1 ? startPos-1 : 0); i <= 15; i++)
		{
			if (((strips[n] >> i) & 0x1) == 0x1)
			{
				posUp[n] = i;
				posUpCount++;
				break;
			}
		}
		for (AliHLTInt32_t i = startPos+1; i >= 0; i--)
		{
			if (((strips[n] >> i) & 0x1) == 0x1)
			{
				posDown[n] = i;
				posDownCount++;
				break;
			}
		}
		
		// 20 Nov 2009: Added scanning on either side of startPos to find the
		// nearest strip to startPos for pathological cases, where there is
		// a larger angle or scatter in Y strips than +/- 1 strip, eg. cosmics.
		if (((strips[n] >> startPos) & 0x1) == 0x1)
		{
			posNearest[n] = startPos;
			posNearestCount++;
		}
		else
		{
			for (AliHLTInt32_t i = 1; i < 16; ++i)
			{
				if (((strips[n] >> (startPos + i)) & 0x1) == 0x1)
				{
					posNearest[n] = startPos + i;
					posNearestCount++;
					break;
				}
				else if (((strips[n] >> (startPos - i)) & 0x1) == 0x1)
				{
					posNearest[n] = startPos - i;
					posNearestCount++;
					break;
				}
			}
		}
	}
	
	// Now select either posUp or posDown, whichever has the most found strips.
	if (posUpCount >= posDownCount and posUpCount >= posNearestCount)
	{
		for (AliHLTInt32_t n = 0; n < 4; n++)
			pos[n] = posUp[n];
	}
	else if (posDownCount >= posUpCount and posDownCount >= posNearestCount)
	{
		for (AliHLTInt32_t n = 0; n < 4; n++)
			pos[n] = posDown[n];
	}
	else
	{
		for (AliHLTInt32_t n = 0; n < 4; n++)
			pos[n] = posNearest[n];
	}
}


const AliHLTMUONTriggerRecoLutRow& AliHLTMUONTriggerReconstructor::AliDecoderHandler::GetLutRowX(
		AliHLTInt32_t xPos, AliHLTUInt8_t chamber
	)
{
	// Fetches the appropriate LUT row for a given strip X and Y position.
	
	assert( chamber <= 3 );
	assert( fCurrentCrateId < 16 );
	
	int locId = 0;
	int pos = 0;
	if (xPos >= 32)
	{
		locId = fUseLocalId ? GetLocalId(fNextStruct) : fNextLocalIndex;
		pos = xPos - 32;
	}
	else if (xPos >= 16)
	{
		locId = fUseLocalId ? GetLocalId(fCurrentStruct) : fNextLocalIndex-1;
		pos = xPos - 16;
	}
	else if (xPos >= 0)
	{
		locId = fUseLocalId ? GetLocalId(fPrevStruct) : fNextLocalIndex-2;
		pos = xPos;
	}
	if (locId < 0 or locId >= 16) locId = 0;
	if (pos < 0 or pos >= 16) pos = 0;
	
	return fLookupTable.fRow[fCurrentCrateId][locId][chamber][0][pos];
}


void AliHLTMUONTriggerReconstructor::AliDecoderHandler::ReconstructHit(
		AliHLTUInt64_t xStrips, AliHLTUInt32_t yStrips,
		AliHLTInt32_t xPos, AliHLTInt32_t yPos, AliHLTUInt8_t yLocId,
		AliHLTUInt8_t chamber, AliHLTMUONRecHitStruct& hit
	)
{
	/// Reconstructs the hit coordinates for the given chamber from the
	/// strip and fired strip information provided.
	/// \param [in] xStrips  The X strip pattern for the given chamber.
	/// \param [in] yStrips  The Y strip pattern for the given chamber.
	/// \param [in] xPos  The position of the X strip that was fired.
	/// \param [in] yPos  The position of the Y strip that was fired.
	/// \param [in] chamber  The chamber on which the strips were found.
	///      Valid range [0..3].
	/// \param [out] hit  This will be filled with the reconstructed hit.

	assert( 0 <= xPos and xPos < 48 );
	assert( 0 <= yPos and yPos < 16 );
	assert( ((xStrips >> xPos) & 0x1) == 0x1 );
	assert( ((yStrips >> yPos) & 0x1) == 0x1 );
	assert( chamber <= 3 );
	assert( fCurrentCrateId < 16 );
	assert( yLocId < 16 );
	
	const AliHLTMUONTriggerRecoLutRow& lut = GetLutRowX(xPos, chamber);
	
	// Decode the Y position of the hit from the strip position information.
	// If two neighbouring strips were fired then we decluster the strips like
	// the L0 electronics does by taking the middle position of the two strips.
	if (xPos > 0 and ((xStrips >> (xPos-1)) & 0x1) == 0x1)
	{
		if (((xStrips >> (xPos+1)) & 0x1) == 0x1)
		{
			// Strips fired on both sides of strip at xPos so just use the middle one.
			hit.fFlags = lut.fIdFlags;
			hit.fY = lut.fY;
			hit.fZ = lut.fZ;
		}
		else
		{
			// Second strip fired below the one at xPos, so decluster.
			assert(xPos-1 < 48);
			const AliHLTMUONTriggerRecoLutRow& lut2 = GetLutRowX(xPos-1, chamber);
			hit.fFlags = lut.fIdFlags;
			hit.fY = (lut2.fY + lut.fY) * 0.5;
			hit.fZ = (lut2.fZ + lut.fZ) * 0.5;
		}
	}
	else
	{
		if (((xStrips >> (xPos+1)) & 0x1) == 0x1)
		{
			// Second strip fired above the one at xPos, so decluster.
			assert(xPos+1 < 48);
			const AliHLTMUONTriggerRecoLutRow& lut2 = GetLutRowX(xPos+1, chamber);
			hit.fFlags = lut.fIdFlags;
			hit.fY = (lut2.fY + lut.fY) * 0.5;
			hit.fZ = (lut2.fZ + lut.fZ) * 0.5;
		}
		else
		{
			// Only strip at xPos fired and neither of its two neighbours.
			hit.fFlags = lut.fIdFlags;
			hit.fY = lut.fY;
			hit.fZ = lut.fZ;
		}
	}
	
	// Similarly decode the X position of the hit from the strip position information.
	// Also decluster like for the Y strips.
	if (yPos > 0 and ((yStrips >> (yPos-1)) & 0x1) == 0x1)
	{
		if (((yStrips >> (yPos+1)) & 0x1) == 0x1)
		{
			// Strips fired on both sides of strip at yPos so just use the middle one.
			hit.fX = fLookupTable.fRow[fCurrentCrateId][yLocId][chamber][1][yPos].fX;
		}
		else
		{
			// Second strip fired below the one at yPos, so decluster.
			assert(yPos-1 < 16);
			hit.fX = (fLookupTable.fRow[fCurrentCrateId][yLocId][chamber][1][yPos-1].fX
				+ fLookupTable.fRow[fCurrentCrateId][yLocId][chamber][1][yPos].fX) * 0.5;
		}
	}
	else
	{
		if (((yStrips >> (yPos+1)) & 0x1) == 0x1)
		{
			// Second strip fired above the one at yPos, so decluster.
			assert(yPos+1 < 16);
			hit.fX = (fLookupTable.fRow[fCurrentCrateId][yLocId][chamber][1][yPos+1].fX
				+ fLookupTable.fRow[fCurrentCrateId][yLocId][chamber][1][yPos].fX) * 0.5;
		}
		else
		{
			// Only strip at yPos fired and neither of its two neighbours.
			hit.fX = fLookupTable.fRow[fCurrentCrateId][yLocId][chamber][1][yPos].fX;
		}
	}
}


void AliHLTMUONTriggerReconstructor::AliDecoderHandler::OnNewRegionalStructV2(
		UInt_t num,
		const AliMUONRegionalHeaderStruct* regionalStruct,
		const AliMUONRegionalScalarsStruct* /*scalars*/,
		const void* /*data*/
	)
{
	// Reset the local trigger structure pointers, and mark the current regional
	// structure number and Crate ID.
	
	fCurrentRegional = num;
	fCurrentCrateId = (fUseCrateId ? GetRegionalId(regionalStruct) : num);
	fPrevStruct = fCurrentStruct = fNextStruct = &fgkNullStruct;
}


void AliHLTMUONTriggerReconstructor::AliDecoderHandler::OnEndOfRegionalStructV2(
		UInt_t /*num*/,
		const AliMUONRegionalHeaderStruct* /*regionalStruct*/,
		const AliMUONRegionalScalarsStruct* /*scalars*/,
		const void* /*data*/
	)
{
	// Process the last local trigger structure.
	
	fPrevStruct = fCurrentStruct;
	fCurrentStruct = fNextStruct;
	fNextStruct = &fgkNullStruct;
	
	// The index numbers for fPrevStruct and fCurrentStruct are calculated from
	// fNextLocalIndex in ProcessLocalStruct so we need to increment it correctly.
	++fNextLocalIndex;
	
	ProcessLocalStruct();
}


void AliHLTMUONTriggerReconstructor::AliDecoderHandler::OnLocalStructV2(
		UInt_t iloc,
		const AliMUONLocalInfoStruct* localStruct,
		const AliMUONLocalScalarsStruct* /*scalars*/
	)
{
	// Update pointers and process the current local trigger structure.
	
	assert(iloc < 16);
	assert(localStruct != NULL);
	
	fPrevStruct = fCurrentStruct;
	fCurrentStruct = fNextStruct;
	fNextStruct = localStruct;
	fNextLocalIndex = iloc;
	ProcessLocalStruct();
}


void AliHLTMUONTriggerReconstructor::AliDecoderHandler::ProcessLocalStruct()
{
	/// Converts the fCurrentStruct local trigger structure from the L0 into a trigger record.
	/// The dHLT trigger records is then used as a seed for tracking algorithms.
	/// \note fOutputTrigRecs must be set before calling the decoder to decode
	///    a new raw data buffer.

	assert(fOutputTrigRecs != NULL);
	
	// If the current local trigger structure does not have a decision then skip it.
	if (GetLocalDec(fCurrentStruct) == 0) return;
	
	// First try to identify the X and Y strips on MT1 that fired the trigger.
	// Note: X strips are for the Y axis in ALICE coordinate system,
	// i.e. bending plane. and Y strips for the X axis.
	AliHLTInt32_t xPos, yPos;
	if (not FindStripsOnMT1(xPos, yPos)) return;
	
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
	
	// Now try find all the fired X and Y strips on all 4 trigger chambers.
	
	AliHLTUInt64_t xStrips[4];
	SelectXPatterns(xStrips);
	AliHLTInt32_t stripPosX[4];
	// Note: the +16 is because FindStripsOnMT1 returns value in the range [0..15]
	// for fCurrentStruct, but we need the value in the range [0..47].
	FindXStrips(xPos+16, xStrips, stripPosX);
	AliHLTUInt32_t yStrips[4]; AliHLTUInt8_t locId[4];
	SelectYPatterns(stripPosX, yStrips, locId);
	AliHLTInt32_t stripPosY[4];
	FindYStrips(yPos, yStrips, stripPosY);
	
	// hitset indicates which hits on chambers 7 to 10 have been found and filled.
	bool hitset[4] = {false, false, false, false};
	
	// Reconstruct the hits from the found strips. Also, fill the hitset
	// flags and make sure the hits for which no valid strips were found get
	// set to a nil value.
	int hitCount = 0;
	for (int i = 0; i < 4; i++)
	{
		if (stripPosX[i] != -1 and stripPosY[i] != -1)
		{
			ReconstructHit(
					xStrips[i], yStrips[i],
					stripPosX[i], stripPosY[i],
					locId[i], i, trigger.fHit[i]
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
	
	if (hitCount < 3)
	{
		// If we could not find at least 3 hits, but the trigger fired, then
		// maybe we have a pathalogical case where 3 X strips and 3 Y strips
		// fired but one chamber has an X, one a Y and only the other 2 have both
		// X and Y strips fired.
		// In such a case we need to try fit a line to X and Y independantly
		// and form the hits from the best line fit.
		
		AliHLTFloat32_t x[4], zx[4], y[4], zy[4];
		AliHLTUInt32_t nx = 0, ny = 0;
		for (int i = 0; i < 4; i++)
		{
			if (stripPosX[i] != -1)
			{
				const AliHLTMUONTriggerRecoLutRow& lut = GetLutRowX(stripPosX[i], i);
				y[ny] = lut.fY;
				zy[ny] = lut.fZ;
				++ny;
			}
			if (stripPosY[i] != -1)
			{
				const AliHLTMUONTriggerRecoLutRow& lut =
					fLookupTable.fRow[fCurrentCrateId][locId[i]][i][1][stripPosY[i]];
				x[nx] = lut.fX;
				zx[nx] = lut.fZ;
				++nx;
			}
		}
		
		AliHLTFloat32_t mx = 0, cx = 0, my = 0, cy = 0;
		bool xfitted = AliHLTMUONCalculations::FitLineToData(x, zx, nx);
		mx = AliHLTMUONCalculations::Mzx();
		cx = AliHLTMUONCalculations::Czx();
		bool yfitted = AliHLTMUONCalculations::FitLineToData(y, zy, ny);
		my = AliHLTMUONCalculations::Mzx();
		cy = AliHLTMUONCalculations::Czx();
		if (xfitted and yfitted)
		{
			for (int i = 0; i < 4; i++)
			{
				if (hitset[i]) continue;  // Leave the found hits alone.
				if (stripPosX[i] != -1)
				{
					// Got X strip but no hit, so Y strip is missing.
					// Thus we have a good Y coordinate but poor X.
					const AliHLTMUONTriggerRecoLutRow& lut = GetLutRowX(stripPosX[i], i);
					trigger.fHit[i].fFlags = lut.fIdFlags;
					trigger.fHit[i].fX = mx * lut.fZ + cx;
					trigger.fHit[i].fY = lut.fY;
					trigger.fHit[i].fZ = lut.fZ;
					hitset[i] = true;
					hitCount++;
				}
				else if (stripPosY[i] != -1)
				{
					// Got Y strip but no hit, so X strip is missing.
					// Thus we have a good X coordinate but poor Y.
					const AliHLTMUONTriggerRecoLutRow& lut =
						fLookupTable.fRow[fCurrentCrateId][locId[i]][i][1][stripPosY[i]];
					trigger.fHit[i].fFlags = lut.fIdFlags;
					trigger.fHit[i].fX = lut.fX;
					trigger.fHit[i].fY = my * lut.fZ + cy;
					trigger.fHit[i].fZ = lut.fZ;
					hitset[i] = true;
					hitCount++;
				}
			}
		}
	}
	
	// If 4 hits found then check if they are all good, otherwise find the 3
	// best fitting ones.
	if (hitCount > 3)
	{
		AliHLTFloat32_t dx = AliMUONConstants::TriggerNonBendingReso();
		AliHLTFloat32_t dy = AliMUONConstants::TriggerBendingReso();
		AliHLTMUONCalculations::SigmaX2(dx*dx);
		AliHLTMUONCalculations::SigmaY2(dy*dy);
		
		AliHLTFloat32_t chi2 = AliHLTMUONCalculations::ComputeChi2(trigger, hitset);
		if (chi2 != -1 and chi2 > 5.*4)  // check 5 sigma cut.
		{
			// Poor fit so look for best 3 points.
			int worstHit = -1;
			AliHLTFloat32_t bestchi2 = 1e38;
			for (int j = 0; j < 4; j++)
			{
				bool tmphitset[4] = {true, true, true, true};
				tmphitset[j] = false;
				AliHLTFloat32_t tmpchi2 = AliHLTMUONCalculations::ComputeChi2(trigger, tmphitset);
				if (tmpchi2 * 4 < chi2 * 3 and tmpchi2 < bestchi2)
				{
					bestchi2 = tmpchi2;
					worstHit = j;
				}
			}
			if (worstHit != -1)
			{
				for (int j = 0; j < 4; j++) hitset[j] = true;
				hitset[worstHit] = false;
				trigger.fHit[worstHit] = AliHLTMUONConstants::NilRecHitStruct();
			}
		}
	}

	// Construct the ID from the running counter fTrigRecId and use the
	// regional counter, local counter and DDL id for the bottom 8 bits.
	AliHLTUInt8_t iloc = fNextLocalIndex-1;
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
	
	bool trigAdded = false;
	
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
		trigAdded = true;
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
		trigAdded = true;
	}
	
	if (trigAdded and fStoreInfo)
	{
		// Allocate or reallocate buffer.
		if (fInfoBuffer == NULL)
		{
			try
			{
				fInfoBuffer = new AliHLTMUONTrigRecInfoStruct[256];
			}
			catch (...)
			{
				HLTError("Could not allocate buffer space for debug information.");
				return;
			}
			fInfoBufferSize = 256;
		}
		else if (fInfoBufferCount >= fInfoBufferSize)
		{
			AliHLTMUONTrigRecInfoStruct* newbuf = NULL;
			try
			{
				newbuf = new AliHLTMUONTrigRecInfoStruct[fInfoBufferSize*2];
			}
			catch (...)
			{
				HLTError("Could not allocate more buffer space for debug information.");
				return;
			}
			for (AliHLTUInt32_t i = 0; i < fInfoBufferSize; ++i) newbuf[i] = fInfoBuffer[i];
			delete [] fInfoBuffer;
			fInfoBuffer = newbuf;
			fInfoBufferSize = fInfoBufferSize*2;
		}
		
		fInfoBuffer[fInfoBufferCount].fTrigRecId = trigger.fId;
		for (int i = 0; i < 4; ++i)
		{
			if (trigger.fHit[i] != AliHLTMUONConstants::NilRecHitStruct())
			{
				fInfoBuffer[fInfoBufferCount].fDetElemId[i] =
					AliHLTMUONUtils::GetDetElemIdFromFlags(trigger.fHit[i].fFlags);
			}
			else
			{
				fInfoBuffer[fInfoBufferCount].fDetElemId[i] = -1;
			}
		}
		fInfoBuffer[fInfoBufferCount].fZmiddle = AliHLTMUONCalculations::Zf();
		fInfoBuffer[fInfoBufferCount].fBl = AliHLTMUONCalculations::QBL();
		fInfoBuffer[fInfoBufferCount].fL0Struct = *fCurrentStruct;
		fInfoBuffer[fInfoBufferCount].fL0StructPrev = *fPrevStruct;
		fInfoBuffer[fInfoBufferCount].fL0StructNext = *fNextStruct;
		++fInfoBufferCount;
	}
}


void AliHLTMUONTriggerReconstructor::AliDecoderHandler::OnError(
		ErrorCode code, const void* location
	)
{
	/// Logs an error message if there was a decoding problem with the DDL payload.
	
	long bytepos = long(location) - long(fBufferStart) + sizeof(AliRawDataHeader);
	if (code == kWrongEventType)
	{
		fHadWrongEventTypeError = true;
		
		/// Do not generate an error message if the fDontPrintWrongEventError option is set.
		if (fDontPrintWrongEventError) return;
	}
	else
	{
		fHadNonWrongEventTypeError = true;
	}
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

