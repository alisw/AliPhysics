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
#include "AliEMCALTriggerDCSConfig.h"
#include "AliEMCALTriggerTRUDCSConfig.h"
#include "AliEMCALTriggerSTUDCSConfig.h"
#include "AliEMCALTriggerData.h"
//#include "AliEMCALTriggerPatch.h"
#include "AliLog.h"

#include "AliRawVEvent.h"
#include "AliRawEventHeaderBase.h"
#include "AliRawEvent.h"
#include "AliRawVEquipment.h"
#include "AliRawEquipmentHeader.h"
#include "AliDAQ.h"

#include "Riostream.h"

#include "AliCaloRawAnalyzerFactory.h"

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
//       AliDebug(1, Form("Using default geometry"));
      fGeometry =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
    }
  
  //  fRawAnalyzer = new AliCaloRawAnalyzerFakeALTRO ();
  
  fRawAnalyzer =  (AliCaloRawAnalyzerFakeALTRO*)AliCaloRawAnalyzerFactory::CreateAnalyzer(kFakeAltro);

  fDCSConfig = AliEMCALTriggerDCSConfigDB::Instance();

  //Int_t nRawDigits = 5952; //fGeometry->GetNTotalTRU() * 96;
  for (Int_t i=kMaxDigitIndex;i--;) fRawDigitIndex[i] = -1;
}	

//_______________
AliEMCALTriggerRawDigitMaker::~AliEMCALTriggerRawDigitMaker()
{
	// dtor
}

//_______________
void AliEMCALTriggerRawDigitMaker::SetIO(AliRawReader* reader, AliCaloRawStreamV3& in, AliEMCALTriggerSTURawStream& inSTU, TClonesArray* digits, TClonesArray *data)
{
	// Connect I/O
	
	fRawReader     = reader;
	fCaloRawStream = &in;
	fRawDigits     = digits;
	fSTURawStream  = &inSTU;
	fTriggerData   = data;
}

///
/// Add bunch list of type AliCaloBunchInfo
///
//_________________
void AliEMCALTriggerRawDigitMaker::Add(const std::vector<AliCaloBunchInfo> &bunchlist)
{	
  Int_t    hwAdd   = fCaloRawStream->GetHWAddress();
  UShort_t iRCU    = fCaloRawStream->GetDDLNumber() % 2; // 0/1
  Int_t    iSM     = fCaloRawStream->GetModule();
  
  Int_t iTRU = fGeometry->GetTriggerMapping()->GetTRUIndexFromOnlineHwAdd(hwAdd,iRCU,iSM);
      
  // if (AliDebugLevel())
  // 	{
  //   UShort_t iBranch = ( hwAdd >> 11 ) & 0x1; // 0/1
  // 		printf("===\n");
  // 		printf("| Hw Adress: 0x%x => SM# %2d / RCU# %d / Branch# %d / TRU# %2d / ADC# %2d\n",
  // 			   hwAdd, fCaloRawStream->GetModule(), iRCU, iBranch, iTRU, fCaloRawStream->GetColumn());
  // }
	
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
				
				if (amp) timeSamples[nSamples++] = ((time << 16) & 0xFF0000) | (amp & 0xFFFF);
				
// 				if (AliDebugLevel())
// 				{
// 					printf("ADC# %2d / time: %2d amplitude: %d\n", fCaloRawStream->GetColumn(), time, amp);
// 				}
	
			}
		}
	}
	
	if (fCaloRawStream->GetColumn() > 95 && fCaloRawStream->GetColumn() < 106)
	{
		Int_t nBits = (fCaloRawStream->GetColumn() == 105) ? 6 : 10;
		
		const UShort_t* sig = bunchlist.at(iBunch).GetData();
		
// 		if (AliDebugLevel()) printf("| L0 id in F-ALTRO => bunch length is: %d\n", bunchlist.at(iBunch).GetLength());
		
		for (Int_t i = 0; i < bunchlist.at(iBunch).GetLength(); i++) 
		{
// 			if (AliDebugLevel()) printf("| sig[%3d]: %x\n",i,sig[i]);
										
			for (Int_t j = 0; j < nBits; j++)
			{
				if (sig[i] & ( 1 << j ))
				{
// 					if (AliDebugLevel()) 
// 					{
// 						printf("| Add L0 patch index in TRU# %2d position %2d\n",iTRU,(fCaloRawStream->GetColumn() - 96) * 10 + j);
// 					}
					
					if (fGeometry->GetAbsFastORIndexFromTRU(iTRU, (fCaloRawStream->GetColumn() - 96) * 10 + j, idx))
					{
						if (fRawDigitIndex[idx] >= 0)
						{
							dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
						}
						else
						{
// 							AliDebug(100,"L0: Trying to update trigger info of a non-existent digit!");
							
							fRawDigitIndex[idx] = fRawDigits->GetEntriesFast();
							new((*fRawDigits)[fRawDigits->GetEntriesFast()]) AliEMCALTriggerRawDigit(idx, 0x0, 0);
							
							dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
						}
						
						dig->SetL0Time(iBin);
					}
				}
			}

// FIX ME			
// 			if (fCaloRawStream->GetColumn() == 105 && (sig[i] & (1 << 6))) 
// 			{
// 				fTriggerData->SetL0Trigger(1, iTRU, 1);
// 										   
// 				if (AliDebugLevel()) printf("=======TRU# %2d has issued a L0\n",iTRU);
// 			}
			
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
			
// 			if (AliDebugLevel())
// 			{
// 				printf("| Add TRG digit of id# %4d from TRU# %2d ADC# %2d\n", idx, iTRU, fCaloRawStream->GetColumn());
// 
// 				dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
// 				dig->Print("");
// 
// 				Int_t iSm, iTru, iEta, iPhi, iD[4], iFor;
// 				if (fGeometry->GetPositionInTRUFromAbsFastORIndex(idx, iTru, iEta, iPhi))
// 				{
// 					printf("| Position => TRU: %2d Eta: %2d Phi: %2d\n", iTru, iEta, iPhi);
// 				}
// 				
// 				if (fGeometry->GetPositionInSMFromAbsFastORIndex(idx, iSm, iEta, iPhi))
// 				{
// 					printf("| Position =>  SM: %2d Eta: %2d Phi: %2d\n", iSm, iEta, iPhi);
// 				}
// 								
// 				if (fGeometry->GetCellIndexFromFastORIndex(idx, iD))
// 				{
// 					printf("| tower iDs: ");
// 					for (Int_t i = 0; i < 4; i++)
// 					{
// 						printf("%5d ",iD[i]); 
// 					}
// 					printf("\n");
// 					
// 					for (Int_t i = 0; i < 4; i++)
// 					{
// 						if (fGeometry->GetFastORIndexFromCellIndex(iD[i], iFor))
// 						{
// 							printf("| tower %d to F-OR %d\n",iD[i],iFor);
// 						}
// 					}
// 				}				
// 			}
		}
	}
}

//_______________
void AliEMCALTriggerRawDigitMaker::PostProcess()
{	
  // Post process digits
  
//   AliDebug(2,"Start post processing the raw digit maker");
  Int_t idx;
  const Int_t offsetEMCAL = AliDAQ::DdlIDOffset("EMCAL");
  
  AliEMCALTriggerRawDigit* dig = 0x0;
  
  // Loop over both STU DDLs in order to read out both EMCAL and DCAL (if available)
  Int_t detectorID = 0;       // 0 - EMCAL, 1 - DCAL
  for(Int_t istu = AliDAQ::GetFirstSTUDDL(); istu <= AliDAQ::GetLastSTUDDL(); istu++){
    detectorID = istu - AliDAQ::GetFirstSTUDDL();
    
    AliEMCALTriggerData *trgData = (AliEMCALTriggerData*)fTriggerData->At(detectorID);
    
    fRawReader->Reset();
    fRawReader->Select("EMCAL",istu,istu);
    
    AliRawVEvent *event = (AliRawEvent*)fRawReader->GetEvent();
    if (!event) {
      AliError(Form("STU DDL# %d not available!",istu));
      continue;
    }
    
    Bool_t isSTUin = kFALSE;
    
    Int_t nSubEv = fRawReader->GetEvent()->GetNSubEvents();
    
    for ( Int_t iSubEv=0; iSubEv<nSubEv; iSubEv++)
    {
      AliRawVEvent *subEv = ((AliRawEvent*)fRawReader->GetEvent())->GetSubEvent(iSubEv);
      if ( !subEv ) continue;
      
      for (Int_t iEquip = 0; iEquip < subEv->GetNEquipments(); iEquip++)
      {
        Int_t eqId = subEv->GetEquipment(iEquip)->GetEquipmentHeader()->GetId();
        
        if (eqId == offsetEMCAL + istu) isSTUin = kTRUE;
      }
    }
    
    fRawReader->Reset();
    
    if (isSTUin && fSTURawStream && fSTURawStream->ReadPayLoad())
    {
      trgData->SetL1DataDecoded(1);
      
      for (int i = 0; i < 2; i++) {
        trgData->SetL1GammaThreshold(i, fSTURawStream->GetL1GammaThreshold(i));
        trgData->SetL1JetThreshold(  i, fSTURawStream->GetL1JetThreshold(i)  );
      }
      
      Int_t v0[2] = { static_cast<Int_t>(fSTURawStream->GetV0A()),  static_cast<Int_t>(fSTURawStream->GetV0C())};
      
      Int_t type[19] =
      {
        static_cast<Int_t>(fSTURawStream->GetG(0, 0)),
        static_cast<Int_t>(fSTURawStream->GetG(1, 0)),
        static_cast<Int_t>(fSTURawStream->GetG(2, 0)),
        static_cast<Int_t>(fSTURawStream->GetJ(0, 0)),
        static_cast<Int_t>(fSTURawStream->GetJ(1, 0)),
        static_cast<Int_t>(fSTURawStream->GetJ(2, 0)),
        static_cast<Int_t>(fSTURawStream->GetG(0, 1)),
        static_cast<Int_t>(fSTURawStream->GetG(1, 1)),
        static_cast<Int_t>(fSTURawStream->GetG(2, 1)),
        static_cast<Int_t>(fSTURawStream->GetJ(0, 1)),
        static_cast<Int_t>(fSTURawStream->GetJ(1, 1)),
        static_cast<Int_t>(fSTURawStream->GetJ(2, 1)),
        static_cast<Int_t>(fSTURawStream->GetRawData()),
        static_cast<Int_t>(fSTURawStream->GetRegionEnable()),
        static_cast<Int_t>(fSTURawStream->GetFwVersion())
      };
      
      // Modify DCS config from STU payload content
//       fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig()->SetG(0, 0, type[0]);
//       fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig()->SetG(1, 0, type[1]);
//       fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig()->SetG(2, 0, type[2]);
//       fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig()->SetJ(0, 0, type[3]);
//       fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig()->SetJ(1, 0, type[4]);
//       fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig()->SetJ(2, 0, type[5]);
//       fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig()->SetG(0, 1, type[6]);
//       fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig()->SetG(1, 1, type[7]);
//       fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig()->SetG(2, 1, type[8]);
//       fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig()->SetJ(0, 1, type[9]);
//       fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig()->SetJ(1, 1, type[10]);
//       fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig()->SetJ(2, 1, type[11]);
//       fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig()->SetRawData(type[12]);
//       fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig()->SetRegion(type[13]);
//       fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig()->SetFw(type[14]);
            
      if (detectorID) {
        for (int i=0;i<4;i++) {
          type[15+i]=static_cast<Int_t>(fSTURawStream->GetPHOSScale(i));
//           fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig(kTRUE)->SetPHOSScale(i,type[15+i]);
        }
        // Modify DCS config from STU payload content
//         fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig(kTRUE)->SetG(0, 0, type[0]);
//         fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig(kTRUE)->SetG(1, 0, type[1]);
//         fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig(kTRUE)->SetG(2, 0, type[2]);
//         fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig(kTRUE)->SetJ(0, 0, type[3]);
//         fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig(kTRUE)->SetJ(1, 0, type[4]);
//         fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig(kTRUE)->SetJ(2, 0, type[5]);
//         fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig(kTRUE)->SetG(0, 1, type[6]);
//         fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig(kTRUE)->SetG(1, 1, type[7]);
//         fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig(kTRUE)->SetG(2, 1, type[8]);
//         fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig(kTRUE)->SetJ(0, 1, type[9]);
//         fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig(kTRUE)->SetJ(1, 1, type[10]);
//         fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig(kTRUE)->SetJ(2, 1, type[11]);
//         fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig(kTRUE)->SetRawData(type[12]);
//         fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig(kTRUE)->SetRegion(type[13]);
//         fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig(kTRUE)->SetFw(type[14]);
      }
      
      trgData->SetL1FrameMask(fSTURawStream->GetFrameReceived());
      trgData->SetL1V0(v0);
      trgData->SetL1TriggerType(type);
      
      trgData->SetL1RawData(fSTURawStream->GetRawData());
      
      Int_t iTRU, jTRU, x, y, jADC;
      
      TVector2 sizeL1gsubr, sizeL1gpatch, sizeL1jsubr, sizeL1jpatch;
      fDCSConfig->GetTriggerDCSConfig()->GetSTUDCSConfig()->GetSegmentation(sizeL1gsubr, sizeL1gpatch, sizeL1jsubr, sizeL1jpatch);
      
      if (fSTURawStream->GetRawData())
      {
//         if (AliDebugLevel()) printf("| STU => TRU raw data are there!\n");
        
        Int_t nTRU = fSTURawStream->GetnTRU();
        for (Int_t i = 0; i < nTRU; i++)
        {
          iTRU = fGeometry->GetTRUIndexFromSTUIndex(i, detectorID);
          
          UInt_t adc[96]; for (Int_t j = 0; j < 96; j++) adc[j] = 0;
          
          fSTURawStream->GetADC(i, adc);
          
          for (Int_t j = 0; j < 96; j++)
          {
            if (adc[j] <= 0) continue;
            
//             AliDebug(10,Form("| STU => TRU# %2d raw data: ADC# %2d: %d\n", iTRU, j, adc[j]));
            
            if(fGeometry->GetTriggerMappingVersion() == 1){
              fGeometry->GetAbsFastORIndexFromTRU(iTRU, j, idx);
            } else {
              fGeometry->GetTRUFromSTU(i, j, jTRU, jADC, detectorID);
              fGeometry->GetAbsFastORIndexFromTRU(jTRU, jADC, idx);
            }
            
            if (fRawDigitIndex[idx] >= 0)
            {
              dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
              
//               if (!dig->GetNSamples()) AliDebug(10,Form("TRG digit of id: %4d found in STU but has no time sample in F-ALTRO!",idx));
            }
            else
            {
//               AliDebug(10,Form("TRG digit of id: %4d found in STU but not in F-ALTRO! Create a new digit!",idx));
              
              fRawDigitIndex[idx] = fRawDigits->GetEntriesFast();
              new((*fRawDigits)[fRawDigits->GetEntriesFast()]) AliEMCALTriggerRawDigit(idx, 0x0, 0);
              
              dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
            }
            
            dig->SetL1TimeSum(adc[j]);
          }
        }
      }
      
      // List of patches in EMCal coordinate system
      
//       for (Int_t i = 0; i < fSTURawStream->GetNL0GammaPatch(); i++)
//       {
//         fSTURawStream->GetL0GammaPatch(i, iTRU, x);
//         
//         iTRU = fGeometry->GetTRUIndexFromSTUIndex(iTRU, detectorID);
//         
//         const Int_t sizePatchL0 =
//         ((AliEMCALTriggerTRUDCSConfig*)fDCSConfig->GetTriggerDCSConfig()->GetTRUArr()->At(fGeometry->GetOnlineIndexFromTRUIndex(iTRU)))->GetSegmentation()
//         *
//         ((AliEMCALTriggerTRUDCSConfig*)fDCSConfig->GetTriggerDCSConfig()->GetTRUArr()->At(fGeometry->GetOnlineIndexFromTRUIndex(iTRU)))->GetSegmentation();
//         
//         if (AliDebugLevel()) printf("| STU => Found L0 patch id: %2d in TRU# %2d\n", x, iTRU);
//         
//         Int_t idFastOR[4];
//         for (Int_t j = 0; j < 4; j++) idFastOR[j] = -1;
//         
//         if (fGeometry->GetFastORIndexFromL0Index(iTRU, x, idFastOR, sizePatchL0))
//         {
//           idx = idFastOR[1];
//           
//           Int_t px, py;
//           if (fGeometry->GetPositionInEMCALFromAbsFastORIndex(idx, px, py))
//           {
//             if (AliDebugLevel()) printf("| STU => Add L0 patch at (%2d , %2d)\n", px, py);
//             
//             if (fRawDigitIndex[idx] >= 0)
//             {
//               dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
//             }
//             else
//             {
//               fRawDigitIndex[idx] = fRawDigits->GetEntriesFast();
//               new((*fRawDigits)[fRawDigits->GetEntriesFast()]) AliEMCALTriggerRawDigit(idx, 0x0, 0);
//               
//               dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
//             }
//             
//             dig->SetTriggerBit(kL0,1);
//           }
//         }
//       }
      
      for (int ithr = 0; ithr < 2; ithr++) {
        for (Int_t i = 0; i < fSTURawStream->GetNL1GammaPatch(ithr); i++) {
          if (fSTURawStream->GetL1GammaPatch(i, ithr, iTRU, x, y)) { // col (0..23), row (0..3)
//             if (AliDebugLevel()) printf("| STU => Found L1 gamma patch at (%2d , %2d) in TRU# %2d\n", x, y, iTRU);
            
            Int_t vx, vy, lphi;
            iTRU = fGeometry->GetTRUIndexFromSTUIndex(iTRU, detectorID);
            if (fGeometry->GetTriggerMappingVersion() == 1) {
             
              vx = 23 - x; 
              vy = y + 4 * int(iTRU / 2); // Position in EMCal frame
              if (iTRU % 2) vx += 24; // C side
              vx = vx - int(sizeL1gsubr.X()) * int(sizeL1gpatch.X()) + 1;
              
              fGeometry->GetAbsFastORIndexFromPositionInEMCAL(vx, vy, idx);
              lphi = 63;
            }
            else {
//               fGeometry->GetTRUFromSTU(iTRU, x, y, jTRU, vx, vy, detectorID);
              fGeometry->GetAbsFastORIndexFromPositionInTRU(iTRU, x, y, idx);
              lphi = 103;
            }
            
            if (vx >= 0 && vy < lphi) {
//               if (AliDebugLevel()) printf("| STU => Add L1 gamma [%d] patch at (%2d , %2d)\n", ithr, vx, vy);
              
              if (fRawDigitIndex[idx] >= 0) {
                dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
              }
              else {
                fRawDigitIndex[idx] = fRawDigits->GetEntriesFast();
                new((*fRawDigits)[fRawDigits->GetEntriesFast()]) AliEMCALTriggerRawDigit(idx, 0x0, 0);
                
                dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
              }
              
              dig->SetTriggerBit(kL1GammaHigh + ithr,1);
            }
          }
        }
        
        for (Int_t i = 0; i < fSTURawStream->GetNL1JetPatch(ithr); i++) {
          if (fSTURawStream->GetL1JetPatch(i, ithr, x, y)) { // col (0,15), row (0,11)
            Int_t vx, vy;
            if (fGeometry->GetTriggerMappingVersion() == 1) {
              vx = 11 - y - int(sizeL1jpatch.X()) + 1;              
              vy = 15 - x - int(sizeL1jpatch.Y()) + 1;
            } 
            else {
              vx = y;
              vy = x;
            }

            vx *= int(sizeL1jsubr.X());
            vy *= int(sizeL1jsubr.Y());
            
//             AliDebug(1, Form("| STU => Found L1 jet [%d] patch at (%2d , %2d)\n", ithr, x, y));
            
            if (x >= 0 && y >= 0) {
              if (fGeometry->GetAbsFastORIndexFromPositionInEMCAL(vx, vy, idx)) {
//                 if (AliDebugLevel()) printf("| STU => Add L1 jet patch at (%2d , %2d)\n", x, y);
                
                if (fRawDigitIndex[idx] >= 0) {
                  dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
                } 
                else {
                  fRawDigitIndex[idx] = fRawDigits->GetEntriesFast();
                  new((*fRawDigits)[fRawDigits->GetEntriesFast()]) AliEMCALTriggerRawDigit(idx, 0x0, 0);
                  
                  dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
                }
                
                dig->SetTriggerBit(kL1JetHigh + ithr,1);
              }
            }
          }
        }
      }
      
      if (detectorID) {
        UInt_t sregion[36] = {0};
        fSTURawStream->GetPHOSSubregion(sregion);
        for (int isr=0;isr<36;isr++) {
          if (fGeometry->GetAbsFastORIndexFromPHOSSubregion(isr, idx)) {
            if (fRawDigitIndex[idx] >= 0) {
              dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
            }
            else {
              fRawDigitIndex[idx] = fRawDigits->GetEntriesFast();
              new((*fRawDigits)[fRawDigits->GetEntriesFast()]) AliEMCALTriggerRawDigit(idx, 0x0, 0);
              
              dig = (AliEMCALTriggerRawDigit*)fRawDigits->At(fRawDigitIndex[idx]);
            }
            
            dig->SetL1SubRegion(sregion[isr]);
          }
        }
      }
    }
  }
}

//_______________
void AliEMCALTriggerRawDigitMaker::Reset()
{
	// Reset
	
	Int_t nRawDigits = 4992;
	for (Int_t i = 0; i < nRawDigits; i++) fRawDigitIndex[i] = -1;
}


