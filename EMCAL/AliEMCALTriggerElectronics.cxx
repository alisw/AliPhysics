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
	const Int_t kNTRU = 30;
}

ClassImp(AliEMCALTriggerElectronics)

//__________________
AliEMCALTriggerElectronics::AliEMCALTriggerElectronics(const AliEMCALTriggerDCSConfig *dcsConf) : TObject(),
fTRU(new TClonesArray("AliEMCALTriggerTRU",32)),
fSTU(0x0)
{
	// Ctor
	
	TVector2 rSize;
	
	rSize.Set( 24.,  4. );

	// 32 TRUs
	for (Int_t i=0;i<kNTRU;i++) 
	{
		AliEMCALTriggerTRUDCSConfig* truConf = dcsConf->GetTRUDCSConfig(i);
		new ((*fTRU)[i]) AliEMCALTriggerTRU(truConf, rSize, i % 2);
	}
	
	rSize.Set( 48., 64. );
	
	// 1 STU
	AliEMCALTriggerSTUDCSConfig* stuConf = dcsConf->GetSTUDCSConfig();
	fSTU = new AliEMCALTriggerSTU(stuConf, rSize);
	
	TString str = "map";
	for (Int_t i=0;i<kNTRU;i++) fSTU->Build(str,
											i,
											(static_cast<AliEMCALTriggerTRU*>(fTRU->At(i)))->Map(),
											(static_cast<AliEMCALTriggerTRU*>(fTRU->At(i)))->RegionSize() 
										    );
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

	AliEMCALGeometry* geom = 0x0;
	
	AliRunLoader *rl = AliRunLoader::Instance();
	if (rl->GetAliRun() && rl->GetAliRun()->GetDetector("EMCAL")){
	  AliEMCAL* emcal = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"));
	  if(emcal)geom = emcal->GetGeometry();
	}
	
	if(!geom) geom =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
	
	if(!geom) AliError("Cannot access geometry!");
	
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
		
		Bool_t isOK1 = geom->GetTRUFromAbsFastORIndex(id, iTRU, iADC);
		
		if ((isOK1 && iTRU >= kNTRU) || !isOK1) continue;

		for (Int_t j = 0; j < digit->GetNSamples(); j++)
		{
			Int_t time, amp;
			Bool_t isOK2 = digit->GetTimeSample(j, time, amp);
			
			if (isOK1 && isOK2 && amp) (static_cast<AliEMCALTriggerTRU*>(fTRU->At(iTRU)))->SetADC(iADC, time, amp);
		}
		
		if (geom->GetPositionInEMCALFromAbsFastORIndex(id, px, py)) posMap[px][py] = i;
	}

	Int_t iL0 = 0;

	Int_t timeL0[kNTRU] = {0}, timeL0min = 999;
	
	for (Int_t i=0; i<kNTRU; i++) 
	{
		AliDebug(999, Form("===========< TRU %2d >============\n", i));
		
		AliEMCALTriggerTRU* iTRU = static_cast<AliEMCALTriggerTRU*>(fTRU->At(i));
		
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

	AliDebug(999, Form("=== %2d TRU (out of %2d) has issued a L0 / Min L0 time: %d\n", iL0, kNTRU, timeL0min));
	
	AliEMCALTriggerRawDigit* dig = 0x0;
	
	if (iL0 && (!data->GetMode() || !fSTU->GetDCSConfig()->GetRawData())) 
	{
		// Update digits after L0 calculation
		for (Int_t i = 0; i < kNTRU; i++)
		{
			AliEMCALTriggerTRU *iTRU = static_cast<AliEMCALTriggerTRU*>(fTRU->At(i));

			Int_t reg[24][4];
			for (int j = 0; j < 24; j++) for (int k = 0; k < 4; k++) reg[j][k] = 0;
					
			iTRU->GetL0Region(timeL0min, reg);
			
			for (int j = 0; j < iTRU->RegionSize()->X(); j++)
			{
				for (int k = 0; k < iTRU->RegionSize()->Y(); k++)
				{
					if (reg[j][k]
							&& 
							geom->GetAbsFastORIndexFromPositionInTRU(i, j, k, id)
							&&
							geom->GetPositionInEMCALFromAbsFastORIndex(id, px, py))
					{
						pos = posMap[px][py];
									
						if (pos == -1)
						{
							// Add a new digit
							new((*digits)[digits->GetEntriesFast()]) AliEMCALTriggerRawDigit(id, 0x0, 0);
										
							dig = (AliEMCALTriggerRawDigit*)digits->At(digits->GetEntriesFast() - 1);
						}
						else
						{
							dig = (AliEMCALTriggerRawDigit*)digits->At(pos);
						}
						
						dig->SetL1TimeSum(reg[j][k]);
					}
				}
			}
		}
	}

	if (iL0 && !data->GetMode())
	{
			// transform local to global 
		
		for (Int_t i = 0; i < kNTRU; i++)
		{
			AliEMCALTriggerTRU* iTRU = static_cast<AliEMCALTriggerTRU*>(fTRU->At(i));
			
			TIter next(&iTRU->Patches());
			while (AliEMCALTriggerPatch* p = (AliEMCALTriggerPatch*)next())
			{
				p->Position(px, py);
			
			// Local 2 Global	
				if (geom->GetAbsFastORIndexFromPositionInTRU(i, px, py, id) && 
								geom->GetPositionInEMCALFromAbsFastORIndex(id, px, py)) p->SetPosition(px, py);
				
				Int_t peaks = p->Peaks();
				
				Int_t sizeX = (Int_t) ((iTRU->PatchSize())->X() * (iTRU->SubRegionSize())->X());
				Int_t sizeY = (Int_t) ((iTRU->PatchSize())->Y() * (iTRU->SubRegionSize())->Y());
					
				for (Int_t j = 0; j < sizeX * sizeY; j++)
				{
					if (peaks & (1 << j))
					{
						pos = posMap[px + j % sizeX][py + j / sizeX];
							
						if (pos == -1)
						{
								// Add a new digit
							new((*digits)[digits->GetEntriesFast()]) AliEMCALTriggerRawDigit(id, 0x0, 0);
								
							dig = (AliEMCALTriggerRawDigit*)digits->At(digits->GetEntriesFast() - 1);
						}
						else
						{
							dig = (AliEMCALTriggerRawDigit*)digits->At(pos);
						}
							
						dig->SetL0Time(timeL0min);
					}
				}
					
				pos = posMap[px][py];
					
				if (pos == -1)
				{
						// Add a new digit
					new((*digits)[digits->GetEntriesFast()]) AliEMCALTriggerRawDigit(id, 0x0, 0);
						
					dig = (AliEMCALTriggerRawDigit*)digits->At(digits->GetEntriesFast() - 1);
				}
				else
				{
					dig = (AliEMCALTriggerRawDigit*)digits->At(pos);
				}
					
				dig->SetTriggerBit(kL0, 0);
			}
		}
	}
	
	//
	// Prepare STU for L1 calculation
	
	for (int i = 0; i < (fSTU->RegionSize())->X(); i++)
	{
		for (int j = 0; j < (fSTU->RegionSize())->Y(); j++)
		{
			pos = posMap[i][j];
		
			if (pos >= 0) 
			{
				AliEMCALTriggerRawDigit *digit = (AliEMCALTriggerRawDigit*)digits->At(pos);
		
				if (digit->GetL1TimeSum() > -1) region[i][j] = digit->GetL1TimeSum();
			}
		}
	}
	
	fSTU->SetRegion(region);
	
	if (data->GetMode()) 
	{
		fSTU->SetThreshold(kL1Gamma, data->GetL1GammaThreshold());
		fSTU->SetThreshold(kL1Jet,   data->GetL1JetThreshold()  );
	}
	else
	{
		fSTU->ComputeThFromV0(kL1Gamma, V0M); 
		data->SetL1GammaThreshold( fSTU->GetThreshold(kL1Gamma));
		fSTU->ComputeThFromV0(kL1Jet,   V0M);
		data->SetL1JetThreshold(   fSTU->GetThreshold(kL1Jet)  );
	}

	fSTU->L1(kL1Gamma);
		
	TIterator* nP = 0x0;
		
	nP = (fSTU->Patches()).MakeIterator();
	
	while (AliEMCALTriggerPatch* p = (AliEMCALTriggerPatch*)nP->Next()) 
	{			
		p->Position(px, py);
			
		if (geom->GetAbsFastORIndexFromPositionInEMCAL(px, py, id))
		{
			if (posMap[px][py] == -1)
			{
					// Add a new digit
				new((*digits)[digits->GetEntriesFast()]) AliEMCALTriggerRawDigit(id, 0x0, 0);
					
				dig = (AliEMCALTriggerRawDigit*)digits->At(digits->GetEntriesFast() - 1);
			} 
			else
			{
				dig = (AliEMCALTriggerRawDigit*)digits->At(posMap[px][py]);								
			}
				
			dig->SetTriggerBit(kL1Gamma,0);
		}
	}

	fSTU->Reset();

	fSTU->L1(kL1Jet);
		
	nP = (fSTU->Patches()).MakeIterator();
			
	while (AliEMCALTriggerPatch* p = (AliEMCALTriggerPatch*)nP->Next()) 
	{			
		p->Position(px, py);

		px *= (Int_t)((fSTU->SubRegionSize())->X());

		py *= (Int_t)((fSTU->SubRegionSize())->Y());
			
		if (geom->GetAbsFastORIndexFromPositionInEMCAL(px, py, id))
		{
			if (posMap[px][py] == -1)
			{
					// Add a new digit
				new((*digits)[digits->GetEntriesFast()]) AliEMCALTriggerRawDigit(id, 0x0, 0);
					
				dig = (AliEMCALTriggerRawDigit*)digits->At(digits->GetEntriesFast() - 1);
			} 
			else
			{
				dig = (AliEMCALTriggerRawDigit*)digits->At(posMap[px][py]);
			}
				
			dig->SetTriggerBit(kL1Jet, 0);
		}
	}

	if (AliDebugLevel() >= 999) data->Scan();
	
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
