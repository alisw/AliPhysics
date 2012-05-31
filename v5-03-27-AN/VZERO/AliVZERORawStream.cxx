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

///////////////////////////////////////////////////////////////////////////////
///
/// This is a class for reading the VZERO DDL raw data
/// The format of the raw data corresponds to the one
/// implemented in AliVZEROBuffer class. 
///
///////////////////////////////////////////////////////////////////////////////

#include "AliVZERORawStream.h"
#include "AliRawReader.h"
#include "AliLog.h"
#include "AliDAQ.h"
#include "AliVZEROCalibData.h"
#include "AliVZEROTriggerData.h"
ClassImp(AliVZERORawStream)

//_____________________________________________________________________________
AliVZERORawStream::AliVZERORawStream(AliRawReader* rawReader) :
  fTrigger(0),
  fTriggerMask(0),
  fPosition(-1),
  fRawReader(rawReader),
  fData(NULL)
{
  // create an object to read VZERO raw data
  //
  // select the raw data corresponding to
  // the VZERO detector id
  fRawReader->Reset();
  AliDebug(1,Form("Selecting raw data for detector %d",AliDAQ::DetectorID("VZERO")));
  fRawReader->Select("VZERO");

  // Initalize the containers
  for(Int_t i = 0; i < kNChannels; i++) {
    fTime[i] = fWidth[i] = 0;
    for(Int_t j = 0; j < kNEvOfInt; j++) {
      fADC[i][j] = 0;
      fIsInt[i][j] = fIsBB[i][j] = fIsBG[i][j] = kFALSE;
    }
    fBBScalers[i] = fBGScalers[i] = 0;
    for(Int_t j = 0; j < kNBunches; j++) {
      fChargeMB[i][j] = 0;
      fIsIntMB[i][j] = fIsBBMB[i][j] = fIsBGMB[i][j] = kFALSE;
    }
  }
  for(Int_t i = 0; i < kNScalers; i++) fScalers[i] = 0;
  for(Int_t i = 0; i < kNBunches; i++) fBunchNumbers[i] = 0;
}

//_____________________________________________________________________________
AliVZERORawStream::~AliVZERORawStream()
{
  // destructor
}

//_____________________________________________________________________________
void AliVZERORawStream::Reset()
{
  // reset raw stream params

  // Reinitalize the containers
  for(Int_t i = 0; i < kNChannels; i++) {
    fTime[i] = fWidth[i] = 0;
    for(Int_t j = 0; j < kNEvOfInt; j++) {
      fADC[i][j] = 0;
      fIsInt[i][j] = fIsBB[i][j] = fIsBG[i][j] = kFALSE;
    }
    fBBScalers[i] = fBGScalers[i] = 0;
    for(Int_t j = 0; j < kNBunches; j++) {
      fChargeMB[i][j] = 0;
      fIsIntMB[i][j] = fIsBBMB[i][j] = fIsBGMB[i][j] = kFALSE;
    }
  }
  for(Int_t i = 0; i < kNScalers; i++) fScalers[i] = 0;
  for(Int_t i = 0; i < kNBunches; i++) fBunchNumbers[i] = 0;

  fTrigger = fTriggerMask = 0;
  fPosition = -1;
  fData = NULL;

  if (fRawReader) fRawReader->Reset();
}

//_____________________________________________________________________________
Bool_t AliVZERORawStream::Next()
{
  // read next digit from the VZERO raw data stream
  // return kFALSE in case of error or no digits left

  if (fPosition >= 0) return kFALSE;

  if (!fRawReader->ReadNextData(fData)) return kFALSE;
  if (fRawReader->GetDataSize() == 0) return kFALSE;
     
  if (fRawReader->GetDataSize() != 5936) {
     fRawReader->AddFatalErrorLog(kRawDataSizeErr,Form("size %d != 5936",fRawReader->GetDataSize()));
     AliWarning(Form("Wrong VZERO raw data size: %d, expected 5936 bytes!",fRawReader->GetDataSize()));
     return kFALSE;
  }

  fPosition = 0;

  fTrigger = GetNextWord() & 0xffff;
  fTriggerMask = GetNextWord() & 0xffff;

  for(Int_t iScaler = 0; iScaler < kNScalers; iScaler++)
     fScalers[iScaler] = GetNextWord();

  for(Int_t iBunch = 0; iBunch < kNBunches; iBunch++)
     fBunchNumbers[iBunch] = GetNextWord();
 
  for (Int_t  iCIU = 0; iCIU < 8; iCIU++) { 
 
  // decoding of one Channel Interface Unit numbered iCIU - there are 8 channels per CIU (and 8 CIUs) :
  
    for (Int_t iChannel_Offset = iCIU*8; iChannel_Offset < (iCIU*8)+8; iChannel_Offset=iChannel_Offset+4) { 
      for(Int_t iChannel = iChannel_Offset; iChannel < iChannel_Offset+4; iChannel++) {
        for(Int_t iEvOfInt = 0; iEvOfInt < kNEvOfInt; iEvOfInt++) {
          UShort_t data = GetNextShort();
          fADC[iChannel][iEvOfInt] = data & 0x3ff;
          fIsInt[iChannel][iEvOfInt] = (data >> 10) & 0x1;
        }
      }
      for(Int_t iEvOfInt = 0; iEvOfInt < kNEvOfInt; iEvOfInt=iEvOfInt+2) {
        UShort_t data = GetNextShort();
        for(Int_t iChannel = iChannel_Offset; iChannel < iChannel_Offset+4; iChannel++) {          
          fIsBB[iChannel][iEvOfInt] = (data >>  2*(iChannel-iChannel_Offset)) & 0x1;
          fIsBG[iChannel][iEvOfInt] = (data >> (2*(iChannel-iChannel_Offset)+1)) & 0x1; 
	  if(iEvOfInt < (kNEvOfInt - 1)) {      
             fIsBB[iChannel][iEvOfInt+1] = (data >> (8+ 2*(iChannel-iChannel_Offset))) & 0x1;
             fIsBG[iChannel][iEvOfInt+1] = (data >> (8+ 2*(iChannel-iChannel_Offset)+1)) & 0x1;
	  }
        }
      }

      GetNextShort();

      for(Int_t iChannel = iChannel_Offset; iChannel < iChannel_Offset+4; iChannel++) {
        for(Int_t iBunch = 0; iBunch < kNBunches; iBunch++) {
          UShort_t data = GetNextShort();
          fChargeMB[iChannel][iBunch] = data & 0x3ff;
          fIsIntMB[iChannel][iBunch] = (data >> 10) & 0x1;
        } 
      }
   
      for(Int_t iBunch = 0; iBunch < kNBunches; iBunch=iBunch+2) {
        UShort_t data = GetNextShort();
        for(Int_t iChannel = iChannel_Offset; iChannel < iChannel_Offset+4; iChannel++) {  
          fIsBBMB[iChannel][iBunch] = (data >>  2*iBunch) & 0x1;
          fIsBGMB[iChannel][iBunch] = (data >> (2*iBunch+1)) & 0x1;
	  if(iBunch < (kNBunches - 1)) {
             fIsBBMB[iChannel][iBunch+1] = (data >> (8+2*iBunch)) & 0x1;
             fIsBGMB[iChannel][iBunch+1] = (data >> (8+2*iBunch+1)) & 0x1;
	  }	  
        }
      }
  
      GetNextShort();
   
      for(Int_t iChannel = iChannel_Offset; iChannel < iChannel_Offset+4; iChannel++) {
        fBBScalers[iChannel] = ((ULong64_t)GetNextWord()) << 32;
        fBBScalers[iChannel] |= GetNextWord();
        fBGScalers[iChannel] = ((ULong64_t)GetNextWord()) << 32;
        fBGScalers[iChannel] |= GetNextWord();
      }

    } 

    for(Int_t iChannel = (iCIU*8) + 7; iChannel >= iCIU*8; iChannel--) { 
      UInt_t time = GetNextWord();
      fTime[iChannel]  = time & 0xfff;
      fWidth[iChannel] = ((time >> 12) & 0x7f); // HPTDC used in pairing mode
    }
    
    // End of decoding of one CIU card
    // printf("Number of bytes used at end of reading CIU card number %d %d \n\n", iCIU+1, fPosition); 
    
  } // end of decoding the eight CIUs
    
  return kTRUE;
}

//_____________________________________________________________________________
UInt_t AliVZERORawStream::GetNextWord()
{
  // This method returns the next 32 bit word
  // inside the raw data payload.
  // The method is supposed to be endian (platform)
  // independent.
  if (!fData || fPosition < 0) AliFatal("Raw data payload buffer is not yet initialized !");

  UInt_t word = 0;
  word |= fData[fPosition++];
  word |= fData[fPosition++] << 8;
  word |= fData[fPosition++] << 16;
  word |= fData[fPosition++] << 24;

  return word;
}

//_____________________________________________________________________________
UShort_t AliVZERORawStream::GetNextShort()
{
  // This method returns the next 16 bit word
  // inside the raw data payload.
  // The method is supposed to be endian (platform)
  // independent.
  if (!fData || fPosition < 0) AliFatal("Raw data payload buffer is not yet initialized !");

  UShort_t word = 0;
  word |= fData[fPosition++];
  word |= fData[fPosition++] << 8;

  return word;
}

//_____________________________________________________________________________
void AliVZERORawStream::CalculateChargeForCentrTriggers(AliVZEROTriggerData *triggerData,
							UShort_t &chargeA, UShort_t &chargeC) const
{
  // Use the raw-data payload
  // in order to calculate the total
  // charge (which is used in the
  // centrality triggers) on each side of V0
  if (!triggerData) {
    AliFatal("Trigger configuration data is not provided. Exiting...");
    return;
  }
  chargeA = chargeC = 0;
  for(Int_t iChannel=0; iChannel<64; iChannel++) {
    Int_t offlineCh = GetOfflineChannel(iChannel);
    Int_t board = AliVZEROCalibData::GetBoardNumber(offlineCh);
    Int_t feeChannel = AliVZEROCalibData::GetFEEChannelNumber(offlineCh);
    if (triggerData->GetEnableCharge(board,feeChannel)) {
      Bool_t integ10 = GetIntegratorFlag(iChannel,10);
      UShort_t ch10 = (UShort_t)GetPedestal(iChannel,10);
      UShort_t trPed = (integ10 == kFALSE) ? triggerData->GetPedestal(0,board,feeChannel) : triggerData->GetPedestal(1,board,feeChannel);
      UShort_t trPedCut = (integ10 == kFALSE) ? triggerData->GetPedestalCut(0,board,feeChannel) : triggerData->GetPedestalCut(1,board,feeChannel);
      if (!triggerData->GetPedestalSubtraction(board)) trPed = trPedCut = 0;
      if (ch10 > trPedCut) {
	if (offlineCh < 32) {
	  chargeC += (ch10 - trPed);
	}
	else {
	  chargeA += (ch10 - trPed);
	}
      }
    }
  }
}

//_____________________________________________________________________________
void AliVZERORawStream::CalculateBBandBGFlags(AliVZEROTriggerData *triggerData,
					      UChar_t &nBBA, UChar_t &nBBC,
					      UChar_t &nBGA, UChar_t &nBGC) const
{
  // Use the raw-data payload
  // in order to calculate the total
  // number of beam-beam and beam-gas flags
  // (which is used in the centrality and
  // multiplicity triggers) on each side of V0
  if (!triggerData) {
    AliFatal("Trigger configuration data is not provided. Exiting...");
    return;
  }
  nBBA = nBBC = nBGA = nBGC = 0;
  for(Int_t iChannel=0; iChannel<64; iChannel++) {
    Int_t offlineCh = GetOfflineChannel(iChannel);
    Int_t board = AliVZEROCalibData::GetBoardNumber(offlineCh);
    Int_t feeChannel = AliVZEROCalibData::GetFEEChannelNumber(offlineCh);
    if (triggerData->GetEnableTiming(board,feeChannel)) {
      if (offlineCh < 32) {
	if (GetBBFlag(iChannel,10)) nBBC++;
	if (GetBGFlag(iChannel,10)) nBGC++;
      }
      else {
	if (GetBBFlag(iChannel,10)) nBBA++;
	if (GetBGFlag(iChannel,10)) nBGA++;
      }
    }
  }
}

//_____________________________________________________________________________
void AliVZERORawStream::FillTriggerBits(AliVZEROTriggerData *triggerData)
{
  // Calculate the charge sums and
  // number of trigger flags and then
  // fill the V0 trigger bits word
  // following the trigger logic implemented
  // in the firmware
  UShort_t chargeA,chargeC;
  CalculateChargeForCentrTriggers(triggerData,chargeA,chargeC);
  UChar_t nBBA,nBBC,nBGA,nBGC;
  CalculateBBandBGFlags(triggerData,
			nBBA,nBBC,
			nBGA,nBGC);

  fTrigger = 0;
  //BBA and BBC
  if((nBBC >= triggerData->GetBBCThreshold()) && (nBBA >= triggerData->GetBBAThreshold())) fTrigger |= 1;
  //BBA or BBC
  if((nBBC >= triggerData->GetBBCThreshold()) || (nBBA >= triggerData->GetBBAThreshold())) fTrigger |= (1<<1);
  //BGA and BBC
  if((nBBC >= triggerData->GetBBCForBGThreshold()) && (nBGA >= triggerData->GetBGAThreshold())) fTrigger |= (1<<2);
  //BGA 
  if((nBGA >= triggerData->GetBGAThreshold())) fTrigger |= (1<<3);
  //BGC and BBA
  if((nBGC >= triggerData->GetBGCThreshold()) && (nBBA >= triggerData->GetBBAForBGThreshold())) fTrigger |= (1<<4);
  //BGC 
  if((nBGC >= triggerData->GetBGCThreshold())) fTrigger |= (1<<5);
  //CTA1 and CTC1
  if((chargeC >= triggerData->GetCentralityV0CThrLow()) && (chargeA >= triggerData->GetCentralityV0AThrLow())) fTrigger |= (1<<6);
  //CTA1 or CTC1
  if((chargeC >= triggerData->GetCentralityV0CThrLow()) || (chargeA >= triggerData->GetCentralityV0AThrLow())) fTrigger |= (1<<7);
  //CTA2 and CTC2
  if((chargeC >= triggerData->GetCentralityV0CThrHigh()) && (chargeA >= triggerData->GetCentralityV0AThrHigh())) fTrigger |= (1<<8);
  //CTA2 or CTC2
  if((chargeC >= triggerData->GetCentralityV0CThrHigh()) || (chargeA >= triggerData->GetCentralityV0AThrHigh())) fTrigger |= (1<<9);
  //MTA and MTC
  if(((nBBC >= triggerData->GetMultV0CThrLow()) && (nBBC <= triggerData->GetMultV0CThrHigh())) &&
     ((nBBA >= triggerData->GetMultV0AThrLow()) && (nBBA <= triggerData->GetMultV0AThrHigh()))) 
    fTrigger |= (1<<10);
  //MTA or MTC
  if(((nBBC >= triggerData->GetMultV0CThrLow()) && (nBBC <= triggerData->GetMultV0CThrHigh())) ||
     ((nBBA >= triggerData->GetMultV0AThrLow()) && (nBBA <= triggerData->GetMultV0AThrHigh()))) 
    fTrigger |= (1<<11);
  //BBA 
  if((nBBA >= triggerData->GetBBAThreshold())) fTrigger |= (1<<12);
  //BBC
  if((nBBC >= triggerData->GetBBCThreshold())) fTrigger |= (1<<13);
  //BGA or BGC 
  if((nBGC >= triggerData->GetBGCThreshold()) || (nBGA >= triggerData->GetBGAThreshold())) fTrigger |= (1<<14);
  //(BGA and BBC) or (BGC and BBA) 
  if(((nBBC >= triggerData->GetBBCForBGThreshold()) && (nBGA >= triggerData->GetBGAThreshold())) ||
     ((nBGC >= triggerData->GetBGCThreshold()) && (nBBA >= triggerData->GetBBAForBGThreshold()))) fTrigger |= (1<<15);

}
