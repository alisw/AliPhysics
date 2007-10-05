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

///*
//
//  The trigger reconstructor class is designed to deal the rawdata inputfiles
//  to findout the the reconstructed hits at the trigger DDL. The output is send
//  to the output block for further processing.
//
//  Author : Indranil Das ( indra.das@saha.ac.in || indra.ehep@gmail.com )
// 
//*/

#include "AliHLTMUONTriggerReconstructor.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONCalculations.h"
#include <vector>
#include <cassert>


AliHLTMUONTriggerReconstructor::AliHLTMUONTriggerReconstructor() :
	fMaxRecPointsCount(0),
	fTrigRecId(0)
{
	// ctor
	
	for (Int_t i = 0; i < 8; i++)
	for (Int_t j = 0; j < 16; j++)
	for (Int_t k = 0; k < 4; k++)
	for (Int_t n = 0; n < 2; n++)
	for (Int_t m = 0; m < 16; m++)
	{
		fLookupTable[i][j][k][n][m].fX = 0;
		fLookupTable[i][j][k][n][m].fY = 0;
		fLookupTable[i][j][k][n][m].fZ = 0;
	}
}


AliHLTMUONTriggerReconstructor::~AliHLTMUONTriggerReconstructor()
{
	// dtor
}


bool AliHLTMUONTriggerReconstructor::Run(
		const AliHLTUInt32_t* rawData,
		// TODO: if we are not checking rawDataSize then it means we are
		// not parsing the raw data safely or checking for corruption carefully.
		// This must be fixed at some point.
		AliHLTUInt32_t /*rawDataSize*/,
		AliHLTMUONTriggerRecordStruct* trigRecord,
		AliHLTUInt32_t& nofTrigRec,
		bool suppressPartialTrigs
	)
{
	fMaxRecPointsCount = nofTrigRec;
	
	// nofTrigRec now becomes the output of how many trigger records were found.
	nofTrigRec = 0;
	
	int index = 0;
	int reg_output, reg_phys_trig_occur;
	int iLocIndex,loc,locDec,triggY,sign,loDev,triggX;
	short pattern[2][4]; // 2 stands for two cathode planes and 4 stands for 4 chambers
	
	int phys_trig_occur = (rawData[index]>>30)&0x1; // 1 for physics trigger, 0 for software trigger
	
	if (not phys_trig_occur) // for software trigger
		index += 8 ;// corresponding to scalar words
	
	index += 1 ; // To skip the separator 0xDEADFACE
	
	index += 4 ; // corresponding to global input
	
	index += 1 ; // reaches to global output
	
	if (not phys_trig_occur) index += 10; // corresponds to scalar words
	
	index += 1; // separator 0xDEADBEEF 
	
	for (int iReg = 0; iReg < 8; iReg++)
	{
		index += 1; // DARC Status Word
		index += 1; // Regeional Word
		reg_output = rawData[index] & 0xFF;
		reg_phys_trig_occur = ( rawData[index] >> 31) & 0x1;
		
		index += 2; // 2 words for regional input
		
		index += 1; // L0 counter
		
		if (not reg_phys_trig_occur) index += 10;
		
		index += 1; // end of Regeonal header 0xBEEFFACE
		
		for(int iLoc = 0; iLoc < 16 ; iLoc++)
		{
			iLocIndex = index;
			
			loc = (rawData[index+5] >> 19) &  0xF ;
			
			locDec = (rawData[index+5] >> 15) & 0xF;
			triggY = (rawData[index+5] >> 14) & 0x1;
			sign = (rawData[index+5] >> 9) & 0x1;
			loDev = (rawData[index+5] >> 5) & 0xF ;
			triggX = (loDev >> 4 & 0x1 ) && !(loDev & 0xF);
			
			if( locDec != 0x9 )
			{ // check for Dec
			
				index += 1;
				pattern[0][0] = rawData[index] & 0xFFFF; // x-strip pattern for chamber 0 
				pattern[0][1] = (rawData[index] >> 16) & 0xFFFF; // x-strip pattern for chamber 1
				index += 1; 
				pattern[0][2] = rawData[index] & 0xFFFF; 
				pattern[0][3] = (rawData[index] >> 16) & 0xFFFF; 
				
				index += 1;
				pattern[1][0] = rawData[index] & 0xFFFF; // y-strip pattern for chamber 0
				pattern[1][1] = (rawData[index] >> 16) & 0xFFFF; // y-strip pattern for chamber 0 
				index += 1; 
				pattern[1][2] = rawData[index] & 0xFFFF; 
				pattern[1][3] = (rawData[index] >> 16) & 0xFFFF; 
			
				if (pattern[0][0] || pattern[0][1] || pattern[0][2] || pattern[0][3]
				   || pattern[1][0] || pattern[1][1] || pattern[1][2] || pattern[1][3]
				)
				{
					if (nofTrigRec == fMaxRecPointsCount)
					{
						HLTError("Output buffer is overflowed maximum assiged arraysize : %d, present array index : %d",
							fMaxRecPointsCount, nofTrigRec
						);
						return false;
					}
				
					bool Xset[4] = {false, false, false, false};
					bool Yset[4] = {false, false, false, false};
				
					for (int iChamber = 0; iChamber < 4 ; iChamber++) //4 chambers per DDL 
					for (int iPlane = 0; iPlane < 2 ; iPlane++) // 2 cathode plane
					{
						for (Int_t ibitxy = 0; ibitxy < 16; ++ibitxy)
						{
							if (((pattern[iPlane][iChamber] >> ibitxy) & 0x1) != 0x1)
								continue;
							
							if (iPlane == 1)
							{
								trigRecord[nofTrigRec].fHit[iChamber].fX =
									fLookupTable[iReg][iLoc][iChamber][iPlane][ibitxy].fX;
								Xset[iChamber] = true;
							}
							else
							{
								trigRecord[nofTrigRec].fHit[iChamber].fY =
									fLookupTable[iReg][iLoc][iChamber][iPlane][ibitxy].fY;
								trigRecord[nofTrigRec].fHit[iChamber].fZ =
									fLookupTable[iReg][iLoc][iChamber][iPlane][ibitxy].fZ;
								Yset[iChamber] = true;
							}
							
						}// loop of ibitxy
					}// ichamber, iplane
				
					// hitset indicates which hits on chambers 7 to 10 have been found and filled.
					bool hitset[4] = {false, false, false, false};
					
					// Fill the hitset flags and make sure the hit structures that were not
					// filled (set) get set to a nil value.
					for (int i = 0; i < 4; i++)
					{
						hitset[i] = Xset[i] and Yset[i];
						
						if (not hitset[i])
						{
							trigRecord[nofTrigRec].fHit[i]
								= AliHLTMUONConstants::NilRecHitStruct();
						}
					}
			
					trigRecord[nofTrigRec].fId = fTrigRecId;
				
					// Increment trigger record Id and keep it positive.
					//TODO: handle the wrapparound better.
					if (fTrigRecId < 0x7FFFFFFF)
						fTrigRecId++;
					else
						fTrigRecId = 0;
					
					AliHLTMUONRecHitStruct* hit1 = NULL;
					if (hitset[0])
						hit1 = &trigRecord[nofTrigRec].fHit[0];
					else if (hitset[1])
						hit1 = &trigRecord[nofTrigRec].fHit[1];
					AliHLTMUONRecHitStruct* hit2 = NULL;
					if (hitset[2])
						hit2 = &trigRecord[nofTrigRec].fHit[2];
					else if (hitset[3])
						hit2 = &trigRecord[nofTrigRec].fHit[3];
					
					if (hit1 != NULL and hit2 != NULL)
					{
						// Calculate the momentum and fill in the flags and momentum fields.
						AliHLTMUONCalculations::ComputeMomentum(
								hit1->fX,
								hit1->fY, hit2->fY,
								hit1->fZ, hit2->fZ
							);
						trigRecord[nofTrigRec].fPx = AliHLTMUONCalculations::Px();
						trigRecord[nofTrigRec].fPy = AliHLTMUONCalculations::Py();
						trigRecord[nofTrigRec].fPz = AliHLTMUONCalculations::Pz();
			
						trigRecord[nofTrigRec].fFlags =
							AliHLTMUONUtils::PackTriggerRecordFlags(
								AliHLTMUONCalculations::Sign(),
								hitset
							);
						
						nofTrigRec++;
					}
					else if ((hit1 != NULL or hit2 != NULL) and not suppressPartialTrigs)
					{
						trigRecord[nofTrigRec].fPx = 0;
						trigRecord[nofTrigRec].fPy = 0;
						trigRecord[nofTrigRec].fPz = 0;
			
						trigRecord[nofTrigRec].fFlags =
							AliHLTMUONUtils::PackTriggerRecordFlags(
								kSignUnknown,
								hitset
							);
						
						nofTrigRec++;
					}
				
				}// if any non zero pattern found
			
				index += 1 ; // the last word, important one
			}// Dec Condn
			
			if (not reg_phys_trig_occur)
				index += 45;
				
			index += 1; // end of local Data 0xCAFEFADE
			
			index = iLocIndex + 6; //important to reset the index counter for fake locids like 235 
		}// iLoc loop
	}// iReg Loop
	
	return true;
}
