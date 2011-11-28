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
 


 
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include "AliEMCALTriggerRawDigitMaker.h"
#include "AliEMCALTriggerSTURawStream.h"
#include "AliCaloRawAnalyzerFakeALTRO.h"
#include "AliEMCALTriggerRawDigit.h"
#include "AliCaloRawStreamV3.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliEMCAL.h"
#include "AliCaloBunchInfo.h"
#include "AliRawReader.h"
#include "AliEMCALTriggerDCSConfigDB.h"
#include "AliEMCALTriggerData.h"
#include "AliEMCALTriggerPatch.h"
#include "AliLog.h"

#include "AliRawDataHeader.h"
#include "AliRawVEvent.h"
#include "AliRawEventHeaderBase.h"
#include "AliRawEvent.h"
#include "AliRawVEquipment.h"
#include "AliRawEquipmentHeader.h"

#include "Riostream.h"

#include "AliCaloRawAnalyzerFactory.h"

namespace
{
	const Int_t kSTUEqId = 4652;
}

ClassImp(AliEMCALTriggerRawDigitMaker)

//_______________
AliEMCALTriggerRawDigitMaker::AliEMCALTriggerRawDigitMaker() : TObject(),
fGeometry(0x0),
fRawReader(0x0),
fCaloRawStream(0x0),
fSTURawStream(0x0),
fRawDigits(0x0),
fRawAnalyzer(0x0),
fDCSConfig(0x0),
fTriggerData(0x0)
{
  // def ctor
  
  AliRunLoader* rl = AliRunLoader::Instance();
  if (rl && rl->GetAliRun()){
    AliEMCAL * emcal = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"));
    if(emcal) fGeometry = emcal->GetGeometry();
  }
  
  if(!fGeometry)
    {
      AliDebug(1, Form("Using default geometry"));
      fGeometry =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
    }
  
  //  fRawAnalyzer = new AliCaloRawAnalyzerFakeALTRO ();
  
  fRawAnalyzer =  (AliCaloRawAnalyzerFakeALTRO*)AliCaloRawAnalyzerFactory::CreateAnalyzer(kFakeAltro);

  fDCSConfig = AliEMCALTriggerDCSConfigDB::Instance();
  
  for (Int_t i=0; i<3072; i++) fRawDigitIndex[i] = -1;
}	

//_______________
AliEMCALTriggerRawDigitMaker::~AliEMCALTriggerRawDigitMaker()
{
	// dtor
}

//_______________
void AliEMCALTriggerRawDigitMaker::SetIO(AliRawReader* reader, AliCaloRawStreamV3& in, AliEMCALTriggerSTURawStream& inSTU, TClonesArray* digits, AliEMCALTriggerData* data)
{
	//
	fRawReader     = reader;
	fCaloRawStream = &in;
	fRawDigits     = digits;
	fSTURawStream  = &inSTU;
	fTriggerData   = data;
}

//_______________
void AliEMCALTriggerRawDigitMaker::Add(const std::vector<AliCaloBunchInfo> &bunchlist)
{
	//
	Int_t    hwAdd   = fCaloRawStream->GetHWAddress();
	UShort_t iRCU    = fCaloRawStream->GetDDLNumber() % 2; // 0/1
	UShort_t iBranch = ( hwAdd >> 11 ) & 0x1;              // 0/1
	
	// TRU id	
	Int_t iTRU = ( (iRCU << 1) | iBranch ) - 1; // 0..2
	
	iTRU  = (fCaloRawStream->GetModule() % 2) ? 2 * (2 - iTRU) + 1 : 2 * iTRU;
	
	iTRU += 6 * int(fCaloRawStream->GetModule()/2);
	
	if (AliDebugLevel())
	{
		printf("===\n");
		printf("| Hw Adress: 0x%x => SM# %2d / RCU# %d / Branch# %d / TRU# %2d / ADC# %2d\n",
			   hwAdd, fCaloRawStream->GetModule(), iRCU, iBranch, iTRU, fCaloRawStream->GetColumn());
	}
	
	Int_t idx;
	
	AliEMCALTriggerRawDigit* dig = 0x0;	
	
	Int_t timeSamples[256]; for (Int_t j=0; j<256; j++) timeSamples[j] = 0;
	Int_t nSamples = 0;
	
	UInt_t iBin   = bunchlist.at(0).GetStartBin();
	 Int_t iBunch = 0;
	
	for (UInt_t i = 0; i < bunchlist.size(); i++)
	{
		AliCaloBunchInfo bunch = bunchlist.at(i);
		
		if (iBin > bunch.GetStartBin()) 
		{
			iBin   = bunch.GetStartBin();
			iBunch = i;
		}
		
		if (fCaloRawStream->GetColumn() < 96)
		{
			const UShort_t* sig = bunch.GetData();
			Int_t startBin = bunch.GetStartBin();
			
			for (Int_t iS = 0; iS < bunch.GetLength(); iS++) 
			{
				Int_t time = startBin--;
				Int_t amp  = sig[iS];
				
				if (amp) timeSamples[nSamples++] = ((time << 12) & 0xFF000) | (amp & 0xFFF);
				
				if (AliDebugLevel())
				{
					printf("ADC# %2d / time: %2d amplitude: %d\n", fCaloRawStream->GetColumn(), time, amp);
				}
			}
		}
	}
	
	if (fCaloRawStream->GetColumn() > 95 && fCaloRawStream->GetColumn() < 106)
	{
		Int_t nBits = (fCaloRawStream->GetColumn() == 105) ? 6 : 10;
		
		const UShort_t* sig = bunchlist.at(iBunch).GetData();
		
		if (AliDebugLevel()) printf("| L0 id in F-ALTRO => bunch length is: %d\n", bunchlist.at(iBunch).GetLength());
		
		for (Int_t i = 0; i < bunchlist.at(iBunch).GetLength(); i++) 
		{
			if (AliDebugLevel()) printf("| sig[%3d]: %x\n",i,sig[i]);
										
			for (Int_t j = 0; j < nBits; j++)
			{
				if (sig[i] & ( 1 << j ))
				{
					if (AliDebugLevel()) 
					{
						printf("| Add L0 patch index in TRU# %2d position %2d\n",iTRU,(fCaloRawStream->GetColumn() - 96) * 10 + j);
					}
					
					if (fGeometry->GetAbsFastORIndexFromTRU(iTRU, (fCaloRawStream->GetColumn() - 96) * 10 + j, idx))
					{
						if (fRawDigitIndex[idx] >= 0)
						{
							dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
						}
						else
						{
							AliDebug(100,"L0: Trying to update trigger info of a non-existent digit!");
							
							fRawDigitIndex[idx] = fRawDigits->GetEntriesFast();
							new((*fRawDigits)[fRawDigits->GetEntriesFast()]) AliEMCALTriggerRawDigit(idx, 0x0, 0);
							
							dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
						}
						
						dig->SetL0Time(iBin);
					}
				}
			}
			
			if (fCaloRawStream->GetColumn() == 105 && (sig[i] & (1 << 6))) 
			{
				fTriggerData->SetL0Trigger(1, iTRU, 1);
										   
				if (AliDebugLevel()) printf("=======TRU# %2d has issued a L0\n",iTRU);
			}
			
			iBin--;
		}
	} 
	else
	{
		if (nSamples && fGeometry->GetAbsFastORIndexFromTRU(iTRU, fCaloRawStream->GetColumn(), idx)) 
		{
			if (fRawDigitIndex[idx] < 0)
			{
				fRawDigitIndex[idx] = fRawDigits->GetEntriesFast();
				new((*fRawDigits)[fRawDigits->GetEntriesFast()]) AliEMCALTriggerRawDigit(idx, timeSamples, nSamples);
			}
			else
			{
				dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
				dig->SetTimeSamples(timeSamples, nSamples);
			}
			
			if (AliDebugLevel())
			{
				printf("| Add TRG digit of id# %4d from TRU# %2d ADC# %2d\n", idx, iTRU, fCaloRawStream->GetColumn());

				dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
				dig->Print("");

				Int_t iSm, iTru, iEta, iPhi, iD[4], iFor;
				if (fGeometry->GetPositionInTRUFromAbsFastORIndex(idx, iTru, iEta, iPhi))
				{
					printf("| Position => TRU: %2d Eta: %2d Phi: %2d\n", iTru, iEta, iPhi);
				}
				
				if (fGeometry->GetPositionInSMFromAbsFastORIndex(idx, iSm, iEta, iPhi))
				{
					printf("| Position =>  SM: %2d Eta: %2d Phi: %2d\n", iSm, iEta, iPhi);
				}
								
				if (fGeometry->GetCellIndexFromFastORIndex(idx, iD))
				{
					printf("| tower iDs: ");
					for (Int_t i = 0; i < 4; i++)
					{
						printf("%5d ",iD[i]); 
					}
					printf("\n");
					
					for (Int_t i = 0; i < 4; i++)
					{
						if (fGeometry->GetFastORIndexFromCellIndex(iD[i], iFor))
						{
							printf("| tower %d to F-OR %d\n",iD[i],iFor);
						}
					}
				}				
			}
		}
	}
}

//_______________
void AliEMCALTriggerRawDigitMaker::PostProcess()
{	
	//
        AliDebug(2,"Start post processing the raw digit maker");
	Int_t idx;
	
	AliEMCALTriggerRawDigit* dig = 0x0;
	
	Int_t sizeL1gsubr[2], sizeL1gpatch[2], sizeL1jsubr[2], sizeL1jpatch[2];
	
	fDCSConfig->GetSTUSegmentation(sizeL1gsubr, sizeL1gpatch, sizeL1jsubr, sizeL1jpatch);
	
	fRawReader->Reset();
	fRawReader->Select("EMCAL",44);	

	Bool_t STUin = kFALSE;
	
	Int_t nSubEv = fRawReader->GetEvent()->GetNSubEvents();
	
	for ( Int_t iSubEv=0; iSubEv<nSubEv; iSubEv++)
	{
		AliRawVEvent *SubEv = ((AliRawEvent*)fRawReader->GetEvent())->GetSubEvent(iSubEv);
		if ( !SubEv ) continue;
		
		for (Int_t iEquip = 0; iEquip < SubEv->GetNEquipments(); iEquip++)
		{
			Int_t eqId = SubEv->GetEquipment(iEquip)->GetEquipmentHeader()->GetId();
			
			if (eqId == kSTUEqId) STUin = kTRUE;
		}
	}
	
	fRawReader->Reset();
	
	if (STUin && fSTURawStream && fSTURawStream->ReadPayLoad())
	{
		fTriggerData->SetL1DataDecoded(1);
		
		fTriggerData->SetL1GammaThreshold(fSTURawStream->GetL1GammaThreshold());
		fTriggerData->SetL1JetThreshold(  fSTURawStream->GetL1JetThreshold()  );
		
		Int_t v0[2] = {fSTURawStream->GetV0A(), fSTURawStream->GetV0C()};
		
		Int_t type[8] = 
		{
			fSTURawStream->GetGA(),
			fSTURawStream->GetGB(),
			fSTURawStream->GetGC(),
			fSTURawStream->GetJA(),
			fSTURawStream->GetJB(),
			fSTURawStream->GetJC(),
			fSTURawStream->GetRegionEnable(), 
			fSTURawStream->GetFwVersion()
		};		
		
		fTriggerData->SetL1FrameMask(fSTURawStream->GetFrameReceived());
		fTriggerData->SetL1V0(v0);
		fTriggerData->SetL1TriggerType(type);
		
		Int_t iTRU, x, y;

		if (fSTURawStream->GetRawData())
		{
			if (AliDebugLevel()) printf("| STU => TRU raw data are there!\n");
			
			for (Int_t i = 0; i < 32; i++)
			{
				iTRU = fGeometry->GetTRUIndexFromSTUIndex(i);
				
				UInt_t adc[96]; for (Int_t j = 0; j < 96; j++) adc[j] = 0;
				
				fSTURawStream->GetADC(i, adc);
				/*
				ofstream outfile(Form("data_TRU%d.txt",i),ios_base::trunc);
				
				for (Int_t j = 0; j < 96; j++) 
				{
					outfile << adc[j] << endl;
				}
				
				outfile.close();
				*/
				for (Int_t j = 0; j < 96; j++)
				{
					//if (adc[j] < 5) continue;
					
					if (AliDebugLevel()) printf("| STU => TRU# %2d raw data: ADC# %2d: %d\n", iTRU, j, adc[j]);
					
					fGeometry->GetAbsFastORIndexFromTRU(iTRU, j, idx);
					
					if (fRawDigitIndex[idx] >= 0)
					{
						dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
						
						if (!dig->GetNSamples()) AliDebug(10,Form("TRG digit of id: %4d found in STU but has no time sample in F-ALTRO!",idx));
					}
					else
					{
						AliDebug(10,Form("TRG digit of id: %4d found in STU but not in F-ALTRO! Create a new digit!",idx));
						
						fRawDigitIndex[idx] = fRawDigits->GetEntriesFast();
						new((*fRawDigits)[fRawDigits->GetEntriesFast()]) AliEMCALTriggerRawDigit(idx, 0x0, 0);
						
						dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
					}
					
					dig->SetL1TimeSum(adc[j]);
				}
			}
		}
		
		// List of patches in EMCal coordinate system
		
		for (Int_t i = 0; i < fSTURawStream->GetNL0GammaPatch(); i++)
		{
			fSTURawStream->GetL0GammaPatch(i, iTRU, x);
			
			iTRU = fGeometry->GetTRUIndexFromSTUIndex(iTRU);
			
			if (AliDebugLevel()) printf("| STU => Found L0 patch id: %2d in TRU# %2d\n", x, iTRU);
			
			const Int_t sizePatchL0 = fDCSConfig->GetTRUSegmentation(iTRU) * fDCSConfig->GetTRUSegmentation(iTRU);
			
			Int_t idFastOR[4];
			for (Int_t j = 0; j < 4; j++) idFastOR[j] = -1;
			
			if (fGeometry->GetFastORIndexFromL0Index(iTRU, x, idFastOR, sizePatchL0))
			{
				idx = idFastOR[1];
				
				Int_t px, py;
				if (fGeometry->GetPositionInEMCALFromAbsFastORIndex(idx, px, py))
				{
					if (AliDebugLevel()) printf("| STU => Add L0 patch at (%2d , %2d)\n", px, py);
										
					if (fRawDigitIndex[idx] >= 0)
					{
						dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
					}
					else
					{
						fRawDigitIndex[idx] = fRawDigits->GetEntriesFast();
						new((*fRawDigits)[fRawDigits->GetEntriesFast()]) AliEMCALTriggerRawDigit(idx, 0x0, 0);
			
						dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
					}
		
					dig->SetTriggerBit(kL0,1);
				}
			}
		}
		
		for (Int_t i = 0; i < fSTURawStream->GetNL1GammaPatch(); i++)
		{
			if (fSTURawStream->GetL1GammaPatch(i, iTRU, x, y)) // col (0..23), row (0..3)
			{
				iTRU = fGeometry->GetTRUIndexFromSTUIndex(iTRU);
			
				if (AliDebugLevel()) printf("| STU => Found L1 gamma patch at (%2d , %2d) in TRU# %2d\n", x, y, iTRU);
				
				Int_t vx = 23 - x, vy = y + 4 * int(iTRU / 2); // Position in EMCal frame
				
				if (iTRU % 2) vx += 24; // C side
				
				vx = vx - sizeL1gsubr[0] * sizeL1gpatch[0] + 1;
				
				if (vx >= 0 && vy < 63) 
				{
					if (fGeometry->GetAbsFastORIndexFromPositionInEMCAL(vx, vy, idx))
					{
					if (AliDebugLevel()) printf("| STU => Add L1 gamma patch at (%2d , %2d)\n", vx, vy);
						
						if (fRawDigitIndex[idx] >= 0)
						{
							dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
				}
						else
						{
							fRawDigitIndex[idx] = fRawDigits->GetEntriesFast();
							new((*fRawDigits)[fRawDigits->GetEntriesFast()]) AliEMCALTriggerRawDigit(idx, 0x0, 0);

							dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
						}
		
						dig->SetTriggerBit(kL1Gamma,1);
					}
				}
			}
		}
		
		for (Int_t i = 0; i < fSTURawStream->GetNL1JetPatch(); i++)
		{
			if (fSTURawStream->GetL1JetPatch(i, x, y)) // col (0,15), row (0,11)
			{
				if (AliDebugLevel()) printf("| STU => Found L1 jet patch at (%2d , %2d)\n", x, y);
				
				Int_t ix = sizeL1jsubr[0] * (11 - y - 4 + 1);

				Int_t iy = sizeL1jsubr[1] * (15 - x - 4 + 1);
				
				// FIXME: x = 0 || y = 0 (Olivier's CS) patches a lost?
				
				if (ix >= 0 && iy >= 0)
				{	
					if (fGeometry->GetAbsFastORIndexFromPositionInEMCAL(ix, iy, idx))
					{
					if (AliDebugLevel()) printf("| STU => Add L1 jet patch at (%2d , %2d)\n", ix, iy);
		
						if (fRawDigitIndex[idx] >= 0)
						{
							dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
						}
						else
						{
							fRawDigitIndex[idx] = fRawDigits->GetEntriesFast();
							new((*fRawDigits)[fRawDigits->GetEntriesFast()]) AliEMCALTriggerRawDigit(idx, 0x0, 0);
		
							dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
						}
		
						dig->SetTriggerBit(kL1Jet,1);
					}
				}
			}
		}		
	}
}

//_______________
void AliEMCALTriggerRawDigitMaker::Reset()
{
	//	
	for (Int_t i = 0; i < 3072; i++) fRawDigitIndex[i] = -1;
}


