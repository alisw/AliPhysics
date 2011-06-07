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

#include <TVector2.h>
#include <TClonesArray.h>

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
	//
	fTRU->Delete();
 	delete fSTU;
}

//__________________
void AliEMCALTriggerElectronics::Digits2Trigger(TClonesArray* digits, const Int_t V0M[], AliEMCALTriggerData* data)
{
	//
	AliEMCALGeometry* geom = 0x0;
	
	AliRunLoader *rl = AliRunLoader::Instance();
	if (rl->GetAliRun() && rl->GetAliRun()->GetDetector("EMCAL")){
	  AliEMCAL* emcal = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"));
	  if(emcal)geom = emcal->GetGeometry();
	}
	
	if(!geom) geom =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
	
	if(!geom) AliError("Cannot access geometry!");
	
	//	digits->Sort();
	
	Int_t region[48][64], posMap[48][64];
	for (Int_t i = 0; i < 48; i++) for (Int_t j = 0; j < 64; j++) 
	{
		region[i][j] =  0;
		posMap[i][j] = -1;
	}
	
	for (Int_t i = 0; i < digits->GetEntriesFast(); i++)
	{
		AliEMCALTriggerRawDigit* digit = (AliEMCALTriggerRawDigit*)digits->At(i);
		
		Int_t id = digit->GetId();
		
		Int_t iTRU, iADC;
		
		Bool_t isOK1 = geom->GetTRUFromAbsFastORIndex(id, iTRU, iADC);
		
		if ((isOK1 && iTRU >= kNTRU) || !isOK1) continue;

		for (Int_t j = 0; j < digit->GetNSamples(); j++)
		{
			Int_t time, amp;
			Bool_t isOK2 = digit->GetTimeSample(j, time, amp);
			
			if (isOK1 && isOK2 && amp) (static_cast<AliEMCALTriggerTRU*>(fTRU->At(iTRU)))->SetADC(iADC, time, amp);
		}
		
		Int_t px, py;
		if (geom->GetPositionInEMCALFromAbsFastORIndex(id, px, py))
		{
			posMap[px][py] = i;
			
			if (fSTU->GetRawData() && digit->GetL1TimeSum() >= 0) 
			{
				region[px][py] = digit->GetL1TimeSum();
			}
		}
	}

	Int_t iL0 = 0;

	for (Int_t i=0; i<kNTRU; i++) 
	{
		AliDebug(999, Form("===========< TRU %2d >============\n", i));
		
		AliEMCALTriggerTRU* iTRU = static_cast<AliEMCALTriggerTRU*>(fTRU->At(i));

		// L0 is always computed from F-ALTRO
		if (iTRU->L0()) 
		{
			iL0 += iTRU->L0();
			
			Int_t sizeX = (Int_t) ((iTRU->PatchSize())->X() * (iTRU->SubRegionSize())->X());
			
			Int_t sizeY = (Int_t) ((iTRU->PatchSize())->Y() * (iTRU->SubRegionSize())->Y());
			
			// transform local to global 
			TIter Next(&iTRU->Patches());
			while (AliEMCALTriggerPatch* p = (AliEMCALTriggerPatch*)Next())
			{
				Int_t px, py, id; p->Position(px, py);
				
				if (geom->GetAbsFastORIndexFromPositionInTRU(i, px, py, id) 
					&& 
					geom->GetPositionInEMCALFromAbsFastORIndex(id, px, py)) p->SetPosition(px, py);
				
				if (!data->GetMode()) // Simulation
				{
					Int_t peaks = p->Peaks();
					
					Int_t pos;
					AliEMCALTriggerRawDigit* dig = 0x0;
					
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
							
							dig->SetL0Time(p->Time());
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
					
					dig->SetTriggerBit(kL0,0);
				}
			}
			
			data->SetL0Trigger(0, i, 1);
		}
		else
			data->SetL0Trigger(0, i, 0);
	}

	// A L0 has been issued, run L1
	// Depending on raw data enabled or not in STU data: L1 computation 
	// should be done from F-ALTRO or directly on TRU time sums in STU raw data
	if (iL0) 
	{
		// Use L1 threshold from raw data when reconstructing raw data
		if (data->GetMode())
		{
			fSTU->SetThreshold(kL1Gamma, data->GetL1GammaThreshold());
			fSTU->SetThreshold(kL1Jet,   data->GetL1JetThreshold()  );			
		}
		else
		{
			fSTU->ComputeThFromV0(V0M); // C/A
			data->SetL1GammaThreshold(fSTU->GetThreshold(kL1Gamma));
			data->SetL1JetThreshold(  fSTU->GetThreshold(kL1Jet)  );
		}
		
		if (fSTU->GetRawData())
		{
			// Compute L1 from STU raw data
			fSTU->SetRegion(region);
		}
		else
		{
			// Build STU raw data from F-ALTRO
			TString str = "region";
			for (Int_t i = 0; i < kNTRU; i++) fSTU->Build(str,
														  i, 
														  (static_cast<AliEMCALTriggerTRU*>(fTRU->At(i)))->Region(), 
														  (static_cast<AliEMCALTriggerTRU*>(fTRU->At(i)))->RegionSize());
		}

		fSTU->L1(kL1Gamma);
		
		Int_t id, px, py;
		AliEMCALTriggerRawDigit* dig = 0x0;
		
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
				
				dig->SetTriggerBit(kL1Jet,0);
			}
		}
		
		Int_t** reg = fSTU->Region();
		
		if (!fSTU->GetRawData())
		{
			// Update digits w/ L1 time sum
			// Done in raw digit maker when raw data enabled
			for (Int_t i = 0; i < 48; i++)
			{
				for (Int_t j = 0; j < 64; j++)
				{
					if (reg[i][j])
					{
						if (geom->GetAbsFastORIndexFromPositionInEMCAL(i, j, id))
						{
							if (posMap[i][j] == -1)
							{
								// Add a new digit with L1 time sum
								new((*digits)[digits->GetEntriesFast()]) AliEMCALTriggerRawDigit(id, 0x0, 0);
								
								dig = (AliEMCALTriggerRawDigit*)digits->At(digits->GetEntriesFast() - 1);
							} 
							else
							{
								dig = (AliEMCALTriggerRawDigit*)digits->At(posMap[i][j]);								
							}
							
							dig->SetL1TimeSum(reg[i][j]);
						}
					}
				}
			}
		}
	}

	if (AliDebugLevel() >= 999) data->Scan();
	
	// Now reset the electronics for a fresh start with next event
	Reset();
}

//__________________
void AliEMCALTriggerElectronics::Reset()
{
	//
	TIter NextTRU(fTRU);
	while ( AliEMCALTriggerTRU *TRU = (AliEMCALTriggerTRU*)NextTRU() ) TRU->Reset();
	
	fSTU->Reset();
}
