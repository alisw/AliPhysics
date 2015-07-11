/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/*
 
 
EMCal trigger electronics manager L0/L1
can handle both simulated digits and raw data
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include "AliEMCALTriggerElectronics.h"
#include "AliEMCALTriggerTRU.h"
#include "AliEMCALTriggerSTU.h"
#include "AliEMCALGeometry.h"
#include "AliRunLoader.h"
#include "AliEMCAL.h" 
#include "AliRun.h" 
#include "AliEMCALTriggerDCSConfig.h"
#include "AliEMCALTriggerData.h"
#include "AliEMCALDigit.h"
#include "AliCaloRawStreamV3.h"
#include "AliEMCALTriggerSTURawStream.h"
#include "AliEMCALDigit.h"
#include "AliEMCALTriggerRawDigit.h"
#include "AliEMCALTriggerPatch.h"
#include "AliEMCALTriggerSTUDCSConfig.h"

#include <TVector2.h>

namespace
{
	const Int_t kNTRU = 30; // TODO: kNTRU should be set to / replaced by  fGeometry->GetNTotalTRU() (total number of TRU for a given geom)  after adding 1 STU for DCAL
}

ClassImp(AliEMCALTriggerElectronics)

//__________________
AliEMCALTriggerElectronics::AliEMCALTriggerElectronics(const AliEMCALTriggerDCSConfig *dcsConf) : TObject(),
fTRU(new TClonesArray("AliEMCALTriggerTRU",52)),
fSTU(0x0),
fGeometry(0)
{
	// Ctor
	
	TVector2 rSize;
	
	rSize.Set( 24.,  4. );
	
	AliRunLoader *rl = AliRunLoader::Instance();
	if (rl->GetAliRun() && rl->GetAliRun()->GetDetector("EMCAL")) {
		AliEMCAL *emcal = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"));
		if (emcal) fGeometry = emcal->GetGeometry();
	}
	
	if (!fGeometry) {		
		fGeometry =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
		AliError("Cannot access geometry, create a new default one!");
	}
	
	// 32 TRUs
	for (Int_t i = 0; i < kNTRU; i++) {
                if(i>=(fGeometry->GetNTotalTRU())) continue; //  i <= kNTRU < 62. Prevents fTRU to go out of bonds with EMCALFirstYEARV1 of EMCALCompleteV1 (NTRU<32) TODO: fix the logic
		AliEMCALTriggerTRUDCSConfig *truConf = dcsConf->GetTRUDCSConfig(fGeometry->GetOnlineIndexFromTRUIndex(i));
		if (truConf) new ((*fTRU)[i]) AliEMCALTriggerTRU(truConf, rSize, i % 2);
	}
	
	rSize.Set( 48., 64. );
	
	// 1 STU
	AliEMCALTriggerSTUDCSConfig* stuConf = dcsConf->GetSTUDCSConfig();
	fSTU = new AliEMCALTriggerSTU(stuConf, rSize);
	
	TString str = "map";
	for (Int_t i = 0; i < kNTRU; i++) {
		AliEMCALTriggerTRU *iTRU = static_cast<AliEMCALTriggerTRU*>(fTRU->At(i));
		if (!iTRU) continue;

		fSTU->Build(str,
					i,
					iTRU->Map(),
					iTRU->RegionSize() 
					);
	}
}

//________________
AliEMCALTriggerElectronics::~AliEMCALTriggerElectronics()
{
	// Dtor
	
	fTRU->Delete();
 	delete fSTU;
}

//__________________
void AliEMCALTriggerElectronics::Digits2Trigger(TClonesArray* digits, const Int_t V0M[], AliEMCALTriggerData* data)
{
	// Digits to trigger
	
	Int_t pos, px, py, id; 
	
	Int_t region[48][64], posMap[48][64];
	for (Int_t i = 0; i < 48; i++) for (Int_t j = 0; j < 64; j++) 
	{
		region[i][j] =  0;
		posMap[i][j] = -1;
	}
	
	for (Int_t i = 0; i < digits->GetEntriesFast(); i++)
	{
		AliEMCALTriggerRawDigit* digit = (AliEMCALTriggerRawDigit*)digits->At(i);
		
		id = digit->GetId();
		
		Int_t iTRU, iADC;
		
		Bool_t isOK1 = fGeometry->GetTRUFromAbsFastORIndex(id, iTRU, iADC);
		
		if (!isOK1 || iTRU > kNTRU-1) continue;

		for (Int_t j = 0; j < digit->GetNSamples(); j++)
		{
			Int_t time, amp;
			Bool_t isOK2 = digit->GetTimeSample(j, time, amp);
			
			if (isOK1 && isOK2 && amp) {
				AliDebug(999, Form("=== TRU# %2d ADC# %2d time# %2d signal %d ===", iTRU, iADC, time, amp));

				AliEMCALTriggerTRU * etr = (static_cast<AliEMCALTriggerTRU*>(fTRU->At(iTRU)));
				if (etr) {
				  if (data->GetMode())
				    etr->SetADC(iADC, time, 4 * amp);
				  else
				    etr->SetADC(iADC, time,     amp);
				}
			}
		}
		
		if (fGeometry->GetPositionInEMCALFromAbsFastORIndex(id, px, py)) posMap[px][py] = i;
	}

	Int_t iL0 = 0;

	Int_t timeL0[kNTRU] = {0}, timeL0min = 999;
	
	for (Int_t i = 0; i < kNTRU; i++) 
	{
		AliEMCALTriggerTRU *iTRU = static_cast<AliEMCALTriggerTRU*>(fTRU->At(i));
		if (!iTRU) continue;
		
		AliDebug(999, Form("===========< TRU %2d >============", i));
		
		if (iTRU->L0()) // L0 recomputation: *ALWAYS* done from FALTRO
		{
			iL0++;
			
			timeL0[i] = iTRU->GetL0Time();
			
			if (!timeL0[i]) AliWarning(Form("TRU# %d has 0 trigger time",i));
			
			if (timeL0[i] < timeL0min) timeL0min = timeL0[i];
			
			data->SetL0Trigger(0, i, 1); // TRU# i has issued a L0
		}
		else
			data->SetL0Trigger(0, i, 0);
	}

	AliDebug(999, Form("=== %2d TRU (out of %2d) has issued a L0 / Min L0 time: %d", iL0, fTRU->GetEntriesFast(), timeL0min));
	
	AliEMCALTriggerRawDigit* dig = 0x0;
	
	if (iL0 && (!data->GetMode() || !fSTU->GetDCSConfig()->GetRawData())) 
	{
		// Update digits after L0 calculation
		for (Int_t i = 0; i < kNTRU; i++)
		{
			AliEMCALTriggerTRU *iTRU = static_cast<AliEMCALTriggerTRU*>(fTRU->At(i));
			if (!iTRU) continue;
			
			Int_t reg[24][4];
			for (int j = 0; j < 24; j++) for (int k = 0; k < 4; k++) reg[j][k] = 0;
					
			iTRU->GetL0Region(timeL0min, reg);
			
			for (int j = 0; j < iTRU->RegionSize()->X(); j++)
			{
				for (int k = 0; k < iTRU->RegionSize()->Y(); k++)
				{
					if (reg[j][k]
						&& 
						fGeometry->GetAbsFastORIndexFromPositionInTRU(i, j, k, id)
						&&
						fGeometry->GetPositionInEMCALFromAbsFastORIndex(id, px, py))
					{
						pos = posMap[px][py];
									
						if (pos == -1)
						{
							// Add a new digit
							posMap[px][py] = digits->GetEntriesFast();

							new((*digits)[digits->GetEntriesFast()]) AliEMCALTriggerRawDigit(id, 0x0, 0);
										
							dig = (AliEMCALTriggerRawDigit*)digits->At(digits->GetEntriesFast() - 1);							
						}
						else
						{
							dig = (AliEMCALTriggerRawDigit*)digits->At(pos);
						}
						
						// 14b to 12b STU time sums
						reg[j][k] >>= 2; 
						
						dig->SetL1TimeSum(reg[j][k]);
					}
				}
			}
		}
	}

	if (iL0 && !data->GetMode())
	{
		for (Int_t i = 0; i < kNTRU; i++)
		{
			AliEMCALTriggerTRU *iTRU = static_cast<AliEMCALTriggerTRU*>(fTRU->At(i));
			if (!iTRU) continue;
			
			AliDebug(999, Form("=== TRU# %2d found %d patches", i, (iTRU->Patches()).GetEntriesFast()));
			
			TIter next(&iTRU->Patches());
			while (AliEMCALTriggerPatch* p = (AliEMCALTriggerPatch*)next())
			{
				p->Position(px, py);
			
				// Local 2 Global
				if (fGeometry->GetAbsFastORIndexFromPositionInTRU(i, px, py, id) 
					&& 
					fGeometry->GetPositionInEMCALFromAbsFastORIndex(id, px, py)) p->SetPosition(px, py);
				
				if (AliDebugLevel()) p->Print("");
				
				Int_t peaks = p->Peaks();
								
				Int_t sizeX = (Int_t) ((iTRU->PatchSize())->X() * (iTRU->SubRegionSize())->X());
				Int_t sizeY = (Int_t) ((iTRU->PatchSize())->Y() * (iTRU->SubRegionSize())->Y());
					
				for (Int_t j = 0; j < sizeX * sizeY; j++)
				{
					if (peaks & (1 << j))
					{
						pos = posMap[px + j % sizeX][py + j / sizeX];
							
						if (pos == -1) {
							// Add a new digit
							posMap[px + j % sizeX][py + j / sizeX] = digits->GetEntriesFast();

							new((*digits)[digits->GetEntriesFast()]) AliEMCALTriggerRawDigit(id, 0x0, 0);
								
							dig = (AliEMCALTriggerRawDigit*)digits->At(digits->GetEntriesFast() - 1);							
						} else {
							dig = (AliEMCALTriggerRawDigit*)digits->At(pos);
						}
							
						if (dig) dig->SetL0Time(p->Time());
					}
				}
					
				pos = posMap[px][py];
					
				if (pos == -1) {
					// Add a new digit
					posMap[px][py] = digits->GetEntriesFast();

					new((*digits)[digits->GetEntriesFast()]) AliEMCALTriggerRawDigit(id, 0x0, 0);
						
					dig = (AliEMCALTriggerRawDigit*)digits->At(digits->GetEntriesFast() - 1);
					
				}
				else
				{
					dig = (AliEMCALTriggerRawDigit*)digits->At(pos);
				}
					
				if (dig) dig->SetTriggerBit(kL0, 0);
			}
		}
	}
	
	// Prepare STU for L1 calculation
	for (int i = 0; i < (fSTU->RegionSize())->X(); i++) {
		for (int j = 0; j < (fSTU->RegionSize())->Y(); j++) {
			
			pos = posMap[i][j];
		
			if (pos >= 0) {
				
				AliEMCALTriggerRawDigit *digit = (AliEMCALTriggerRawDigit*)digits->At(pos);
		
				if (digit && digit->GetL1TimeSum() > -1) region[i][j] = digit->GetL1TimeSum();
			}
		}
	}
	
	AliDebug(999,"==================== STU  ====================");
	if (AliDebugLevel() >= 999) fSTU->Scan();
	AliDebug(999,"==============================================");
	
	fSTU->SetRegion(region);
	
	if (data->GetMode()) 
	{
		for (int ithr = 0; ithr < 2; ithr++) {
			AliDebug(999, Form(" THR %d / EGA %d / EJE %d", ithr, data->GetL1GammaThreshold(ithr), data->GetL1JetThreshold(ithr)));
							   
			fSTU->SetThreshold(kL1GammaHigh + ithr, data->GetL1GammaThreshold(ithr));
			fSTU->SetThreshold(kL1JetHigh + ithr, data->GetL1JetThreshold(  ithr));
		}
	}
	else
	{
		for (int ithr = 0; ithr < 2; ithr++) {
			//
			fSTU->ComputeThFromV0(kL1GammaHigh + ithr, V0M); 
			data->SetL1GammaThreshold(ithr, fSTU->GetThreshold(kL1GammaHigh + ithr));
			
			fSTU->ComputeThFromV0(kL1JetHigh + ithr,   V0M);
			data->SetL1JetThreshold(ithr, fSTU->GetThreshold(kL1JetHigh + ithr)  );
			
			AliDebug(999, Form("STU THR %d EGA %d EJE %d", ithr, fSTU->GetThreshold(kL1GammaHigh + ithr), fSTU->GetThreshold(kL1JetHigh + ithr)));
		}
	}

	for (int ithr = 0; ithr < 2; ithr++) {
		//
		fSTU->Reset();
		
		fSTU->L1(kL1GammaHigh + ithr);
		
		TIterator* nP = 0x0;
		
		nP = (fSTU->Patches()).MakeIterator();
		
		AliDebug(999, Form("=== STU found %d gamma patches", (fSTU->Patches()).GetEntriesFast()));
		
		while (AliEMCALTriggerPatch* p = (AliEMCALTriggerPatch*)nP->Next()) 
		{			
			p->Position(px, py);
			
			if (AliDebugLevel()) p->Print("");
			
			if (fGeometry->GetAbsFastORIndexFromPositionInEMCAL(px, py, id))
			{
				if (posMap[px][py] == -1)
				{
					posMap[px][py] = digits->GetEntriesFast();
					
					// Add a new digit
					new((*digits)[digits->GetEntriesFast()]) AliEMCALTriggerRawDigit(id, 0x0, 0);
					
					dig = (AliEMCALTriggerRawDigit*)digits->At(digits->GetEntriesFast() - 1);
				} 
				else
				{
					dig = (AliEMCALTriggerRawDigit*)digits->At(posMap[px][py]);								
				}
				
				if (AliDebugLevel()) dig->Print("");
				
				if (dig) dig->SetTriggerBit(kL1GammaHigh + ithr, 0);
			}
		}
		
		fSTU->Reset();
		
		fSTU->L1(kL1JetHigh + ithr);
		
		nP = (fSTU->Patches()).MakeIterator();
		
		AliDebug(999, Form("=== STU found %d jet patches", (fSTU->Patches()).GetEntriesFast()));
		
		while (AliEMCALTriggerPatch* p = (AliEMCALTriggerPatch*)nP->Next()) 
		{			
			p->Position(px, py);
			
			if (AliDebugLevel()) p->Print("");
			
			if (fGeometry->GetAbsFastORIndexFromPositionInEMCAL(px, py, id))
			{
				if (posMap[px][py] == -1)
				{
					posMap[px][py] = digits->GetEntriesFast();
					
					// Add a new digit
					new((*digits)[digits->GetEntriesFast()]) AliEMCALTriggerRawDigit(id, 0x0, 0);
					
					dig = (AliEMCALTriggerRawDigit*)digits->At(digits->GetEntriesFast() - 1);
				} 
				else
				{
					dig = (AliEMCALTriggerRawDigit*)digits->At(posMap[px][py]);
				}
				
				if (AliDebugLevel()) dig->Print("");
				
				if (dig) dig->SetTriggerBit(kL1JetHigh + ithr, 0);
			}
		}
	}
	
	// Now reset the electronics for a fresh start with next event
	Reset();
}

//__________________
void AliEMCALTriggerElectronics::Reset()
{
	// Reset
	
	TIter nextTRU(fTRU);
	while ( AliEMCALTriggerTRU *TRU = (AliEMCALTriggerTRU*)nextTRU() ) TRU->Reset();
	
	fSTU->Reset();
}
