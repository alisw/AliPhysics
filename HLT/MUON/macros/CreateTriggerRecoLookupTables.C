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

/**
 * This macro is used to generate the lookup tables for the trigger reconstructor
 * component. All alignment and geometry data is taken from the CDB.
 */

#include <iostream>
#include <fstream>

#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpDEIterator.h"
#include "AliMpVSegmentation.h"
#include "AliMpDEManager.h"
#include "AliMpLocalBoard.h"
#include "AliMpTriggerCrate.h"

using namespace std;


struct TriggerRecoLookupTableRow
{
	float fX, fY, fZ;
};

struct TriggerRecoLookupTable
{
	// [regional header index][local board ID][chamber][cathode - X/Y][bit set in bit pattern]
	TriggerRecoLookupTableRow fRow[8][16][4][2][16];
};


void CreateTriggerRecoLookupTables(const char* CDBpath = "local://$ALICE_ROOT")
{
	AliCDBManager* cdbManager = AliCDBManager::Instance();
	cdbManager->SetDefaultStorage(CDBpath);
	cdbManager->SetRun(0);
	AliGeomManager::LoadGeometry();
	
	AliMUONGeometryTransformer transformer;
	if (! transformer.LoadGeometryData())
	{
		cerr << "ERROR: Could not load geometry into transformer." << endl;
		return;
	}
	
	if (! AliMpCDB::LoadMpSegmentation())
	{
		cerr << "ERROR: Could not load segmentation mapping." << endl;
		return;
	}
	AliMpSegmentation* segmentation = AliMpSegmentation::Instance();
	if (segmentation == NULL)
	{
		cerr << "ERROR: AliMpSegmentation::Instance() was NULL." << endl;
		return;
	}
	if (! AliMpCDB::LoadDDLStore())
	{
		cerr << "ERROR: Could not load DDL mapping." << endl;
		return;
	}
	AliMpDDLStore* ddlStore = AliMpDDLStore::Instance();
	if (ddlStore == NULL)
	{
		cerr << "ERROR: AliMpDDLStore::Instance() was NULL." << endl;
		return;
	}
	
	cout << "Building LUTs..." << endl;
	
	TriggerRecoLookupTable* lookupTable21 = new TriggerRecoLookupTable;
	TriggerRecoLookupTable* lookupTable22 = new TriggerRecoLookupTable;
	
	for (Int_t i = 0; i < 8; i++)
	for (Int_t j = 0; j < 16; j++)
	for (Int_t k = 0; k < 4; k++)
	for (Int_t n = 0; n < 2; n++)
	for (Int_t m = 0; m < 16; m++)
	{
		lookupTable21->fRow[i][j][k][n][m].fX = 0;
		lookupTable21->fRow[i][j][k][n][m].fY = 0;
		lookupTable21->fRow[i][j][k][n][m].fZ = 0;
		lookupTable22->fRow[i][j][k][n][m].fX = 0;
		lookupTable22->fRow[i][j][k][n][m].fY = 0;
		lookupTable22->fRow[i][j][k][n][m].fZ = 0;
	}
	
	AliMpDEIterator detElemIter;
	for (Int_t iDDL = 20; iDDL <= 21; iDDL++)
	{
		for (Int_t iReg = 0; iReg < 8; iReg++)
		{
			AliMpTriggerCrate* crate = ddlStore->GetTriggerCrate(iDDL, iReg);
			if (crate == NULL)
			{
				cerr << "ERROR: Could not get crate for regional header = " << iReg
					<< ", and DDL ID = " << iDDL << endl;
				continue;
			}
			
			for (Int_t iLocBoard = 0; iLocBoard < 16; iLocBoard++)
			{
				Int_t boardId = crate->GetLocalBoardId(iLocBoard);
				if (boardId == 0) continue;
				
				AliMpLocalBoard* localBoard = ddlStore->GetLocalBoard(boardId);
				if (localBoard == NULL)
				{
					cerr << "ERROR: Could not get loacl board: " << boardId << endl;
					continue;
				}

				// skip copy cards
				if (! localBoard->IsNotified()) continue;
			
				for (Int_t iChamber = 0; iChamber < 4; iChamber++)
				{
					Int_t detElemId = ddlStore->GetDEfromLocalBoard(boardId, iChamber);
					
					const AliMUONGeometryDetElement* detElemTransform = transformer.GetDetElement(detElemId);
					if (detElemTransform == NULL)
					{
						cerr << "ERROR: Got NULL pointer for geometry transformer for detection element ID = "
							<< detElemId << endl;
						continue;
					}
					
					for (Int_t iCathode = 0; iCathode <= 1; iCathode++)
					{
						const AliMpVSegmentation* seg = segmentation->GetMpSegmentation(
								detElemId, AliMp::GetCathodType(iCathode)
							);
						
						for (Int_t bitxy = 0; bitxy < 16; bitxy++)
						{
							Int_t offset = 0;
							if (iCathode && localBoard->GetSwitch(6)) offset = -8;
							
							AliMpPad pad = seg->PadByLocation(AliMpIntPair(boardId, bitxy+offset), kFALSE);
						
							if (! pad.IsValid())
							{
								// There is no pad associated with the given local board and bit pattern.
								continue;
							}
							
							// Get the global coodinates of the pad.
							Float_t lx = pad.Position().X();
							Float_t ly = pad.Position().Y();
							Float_t gx, gy, gz;
							detElemTransform->Local2Global(lx, ly, 0, gx, gy, gz);
							
							// Fill the LUT
							if (crate->GetDdlId() == 20)
							{
								lookupTable21->fRow[iReg][iLocBoard][iChamber][iCathode][bitxy].fX = gx;
								lookupTable21->fRow[iReg][iLocBoard][iChamber][iCathode][bitxy].fY = gy;
								lookupTable21->fRow[iReg][iLocBoard][iChamber][iCathode][bitxy].fZ = gz;
							}
							else
							{
								lookupTable22->fRow[iReg][iLocBoard][iChamber][iCathode][bitxy].fX = gx;
								lookupTable22->fRow[iReg][iLocBoard][iChamber][iCathode][bitxy].fY = gy;
								lookupTable22->fRow[iReg][iLocBoard][iChamber][iCathode][bitxy].fZ = gz;
							}
						}
					}
				}
			}
		}
	}
	
	fstream file;
	file.open("Lut21.dat", fstream::out | fstream::binary | fstream::trunc);
	if (file)
	{
		file.write((char*)lookupTable21, sizeof(TriggerRecoLookupTable));
		if (! file)
		{
			cerr << "ERROR: There was a problem writing to file Lut21.dat" << endl;
		}
		file.close();
	}
	else
	{
		cerr << "ERROR: Could not open file Lut21.dat for writing." << endl;
	}
	
	file.open("Lut22.dat", fstream::out | fstream::binary | fstream::trunc);
	if (file)
	{
		file.write((char*)lookupTable22, sizeof(TriggerRecoLookupTable));
		if (! file)
		{
			cerr << "ERROR: There was a problem writing to file Lut22.dat" << endl;
		}
		file.close();
	}
	else
	{
		cerr << "ERROR: Could not open file Lut22.dat for writing." << endl;
	}
	
	delete lookupTable21;
	delete lookupTable22;
}
