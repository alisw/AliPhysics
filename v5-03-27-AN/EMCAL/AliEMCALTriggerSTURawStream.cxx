
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



This class provides access to STU DDL raw data.
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include "AliEMCALTriggerSTURawStream.h"
#include "AliRawReader.h"
#include "AliLog.h"

#include "Riostream.h"
#include "TBits.h"

#include <cstdlib>

namespace
{
	const Int_t kPayLoadSizeV0     = 236;
	const Int_t kPayLoadSizeV1     = 245;        
	const Int_t kPayLoadSizeV2     = 390;        
}

ClassImp(AliEMCALTriggerSTURawStream)

//_____________________________________________________________________________
AliEMCALTriggerSTURawStream::AliEMCALTriggerSTURawStream() : TObject(),
fRawReader(0x0),
fL1JetThreshold(),
fL1GammaThreshold(),
fL0GammaPatchIndex(),
fL1GammaPatchIndex(),
fL1JetPatchIndex(),
fNL0GammaPatch(0),
fNL1JetPatch(),
fNL1GammaPatch(),
fGetRawData(0),
fV0A(0),
fV0C(0),
fG(),
fJ(),
fRegionEnable(0),
fFrameReceived(0),
fFwVersion(0)
{
	//
	for (int i = 0; i < 2; i++) {
		//
		fL1JetThreshold[i] = fL1GammaThreshold[i] = 0;
		
		fNL1JetPatch[i] = fNL1GammaPatch[i] = 0;
	}
	
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 2; j++) {
			//
			fG[i][j] = fJ[i][j] = 0;
		}
	}
	
	for (int i = 0; i < 3100; i++) {
		//
		fL0GammaPatchIndex[i] = 0;
		
		for (int j = 0; j < 2; j++) {
			
			fL1GammaPatchIndex[i][j] = 0;
		}
	}
	
	for (int i = 0; i < 200; i++) {
		for (int j = 0; j < 2; j++) {
			
			fL1JetPatchIndex[i][j] = 0;
		}
	}
}

//_____________________________________________________________________________
AliEMCALTriggerSTURawStream::AliEMCALTriggerSTURawStream(AliRawReader* rawReader) : TObject(),
fRawReader(rawReader),
fL1JetThreshold(),
fL1GammaThreshold(),
fL0GammaPatchIndex(),
fL1GammaPatchIndex(),
fL1JetPatchIndex(),
fNL0GammaPatch(0),
fNL1JetPatch(),
fNL1GammaPatch(),
fGetRawData(0),
fV0A(0),
fV0C(0),
fG(),
fJ(),
fRegionEnable(0),
fFrameReceived(0),
fFwVersion(0)
{
	//
	fRawReader->Reset();
	fRawReader->Select("EMCAL",44);
	//
	for (int i = 0; i < 2; i++) {
		//
		fL1JetThreshold[i] = fL1GammaThreshold[i] = 0;
		
		fNL1JetPatch[i] = fNL1GammaPatch[i] = 0;
	}
	
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 2; j++) {
			//
			fG[i][j] = fJ[i][j] = 0;
		}
	}
	
	for (int i = 0; i < 3100; i++) {
		//
		fL0GammaPatchIndex[i] = 0;
		
		for (int j = 0; j < 2; j++) {
			
			fL1GammaPatchIndex[i][j] = 0;
		}
	}
	
	for (int i = 0; i < 200; i++) {
		for (int j = 0; j < 2; j++) {
			
			fL1JetPatchIndex[i][j] = 0;
		}
	}
}

//_____________________________________________________________________________
AliEMCALTriggerSTURawStream::~AliEMCALTriggerSTURawStream()
{
	// destructor
}

//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::Reset()
{
	// Reset
	
	if (fRawReader) fRawReader->Reset();
	
	fNL0GammaPatch = 0;
	fNL1GammaPatch[0] = fNL1GammaPatch[1] = 0;
	fNL1JetPatch[0] = fNL1JetPatch[1] = 0;
}

//_____________________________________________________________________________
Bool_t AliEMCALTriggerSTURawStream::ReadPayLoad()
{
	// STU data decoder from Olivier Bourrion LPSC CNRS-IN2P3
	// bourrion_at_lpsc_dot_in2p3_dot_fr
	
	UInt_t word32[kPayLoadSizeV2 + 1536]; // 32b words
	for (Int_t i = 0;i < kPayLoadSizeV2 + 1536; i++) word32[i] = 0;
	
	Int_t iword = 0;
	
	fNL0GammaPatch = 0;
	fNL1GammaPatch[0] = fNL1GammaPatch[1] = 0;
	fNL1JetPatch[0] = fNL1JetPatch[1] = 0;
	
	Int_t eqId = -1, eqSize = 0;
	
	UInt_t w32;
	
	while (fRawReader->ReadNextInt(w32)) 
	{
		if (!iword)
		{
			eqId   = fRawReader->GetEquipmentId();
			eqSize = fRawReader->GetEquipmentSize();
		}
		
		word32[iword++] = w32;
	}
 
	if (iword != kPayLoadSizeV0 && iword != kPayLoadSizeV1 && iword != kPayLoadSizeV2 
		&& 
		iword != (kPayLoadSizeV0 + 1536) && iword != (kPayLoadSizeV1 + 1536) && iword != (kPayLoadSizeV2 + 1536))
	{
		AliError(Form("STU payload (eqId: %d, eqSize: %d) doesn't match expected size! %d word32",
					  eqId, eqSize, iword));
		return kFALSE;
	}
	
	AliDebug(1, Form("STU (eqId: %d, eqSize: %d) payload size: %d word32",
					 eqId, eqSize, iword));
	
	int offset = 1;//, jetSize = 2;
	
	int nthres = 1;
	
	switch (iword) 
	{
		case kPayLoadSizeV0:
		case kPayLoadSizeV0 + 1536:
		{
			fL1JetThreshold[0]   = ((word32[0]>>16) & 0xFFFF);
			fL1GammaThreshold[0] =  (word32[0]      & 0xFFFF);
			
			break;
		}
		case kPayLoadSizeV1:
		case kPayLoadSizeV1 + 1536:
		{
			fV0A = ((word32[0]>>16) & 0xFFFF);
			fV0C =  (word32[0]      & 0xFFFF);
			
			fG[0][0]       = word32[1];
			fG[1][0]       = word32[2];
			fG[2][0]       = word32[3];
			fJ[0][0]       = word32[4];
			fJ[1][0]       = word32[5];
			fJ[2][0]       = word32[6];		
			fRegionEnable  = word32[7];
			fFrameReceived = word32[8];
			fFwVersion     = word32[9];
			
			fL1JetThreshold[0]   = GetThreshold(fJ[0][0], fJ[1][0], fJ[2][0], fV0A, fV0C);
			fL1GammaThreshold[0] = GetThreshold(fG[0][0], fG[1][0], fG[2][0], fV0A, fV0C); 
			
			offset = 10;
			
			break;
		}
		case kPayLoadSizeV2:
		case kPayLoadSizeV2 + 1536:
		{
			fV0A = ((word32[0]>>16) & 0xFFFF);
			fV0C =  (word32[0]      & 0xFFFF);
			
			fG[0][0]       = word32[1];
			fG[1][0]       = word32[2];
			fG[2][0]       = word32[3];
			fJ[0][0]       = word32[4];
			fJ[1][0]       = word32[5];
			fJ[2][0]       = word32[6];		
			
			fG[0][1]       = word32[7];
			fG[1][1]       = word32[8];
			fG[2][1]       = word32[9];
			fJ[0][1]       = word32[10];
			fJ[1][1]       = word32[11];
			fJ[2][1]       = word32[12];		
			
			fRegionEnable  = word32[13];
			fFrameReceived = word32[14];
			fFwVersion     = word32[15];
			
			for (int i = 0; i < 2; i++) {
				//
				fL1JetThreshold[i]   = GetThreshold(fJ[0][i], fJ[1][i], fJ[2][i], fV0A, fV0C);
				fL1GammaThreshold[i] = GetThreshold(fG[0][i], fG[1][i], fG[2][i], fV0A, fV0C); 
			}
			
			offset = 16;			
			
			nthres = 2;
			
			break;
		}
	}
	
//	jetSize += (fFwVersion >> 16);
	
	///////////
	// START DECODING
	//////////

	for (int i = 0; i < nthres; i++) {
		DecodeL1JetPatchIndexes(i, word32, offset);
		
		offset += 11;
	}
	
	DecodeL0GammaPatchIndexes(word32, offset);
	
	offset += 96;
	
	for (int i = 0; i < nthres; i++) {
		DecodeL1GammaPatchIndexes(i, word32, offset);	
		
		offset += 128;
	}
	
	if (iword == kPayLoadSizeV0 || iword == kPayLoadSizeV1 || iword == kPayLoadSizeV2) {
		fGetRawData = 0;
		return kTRUE;
	}
	
	fGetRawData = 1;
	
	DecodeTRUADC(word32, offset);
	
	return kTRUE;
}

//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::DecodeL0GammaPatchIndexes(UInt_t *word32, const int offset)
{
	//////////////////////////////////////////////////////////
	// index des L0                                         //
	//////////////////////////////////////////////////////////
	// FIXME: sounds like not valid data
	
	unsigned short truL0indexes[32][6];
	
	// extraction from stream
	for (Int_t index=0;index<6;index++)
	{
		for (Int_t tru_num=0;tru_num<16;tru_num++)
		{
			truL0indexes[2*tru_num  ][index] = ( word32[offset + index * 16 + tru_num]        & 0xFFFF);
			truL0indexes[2*tru_num+1][index] = ((word32[offset + index * 16 + tru_num] >> 16) & 0xFFFF);
		}
	}
	
	for (Int_t tru_num=0;tru_num<32;tru_num++) 
	{
		for (Int_t index=0;index<6;index++) 
		{
			for (Int_t bit_num=0;bit_num<12;bit_num++)
			{
				if ((truL0indexes[tru_num][index] & (1 << bit_num)))
				{
					Int_t idx = 12 * index + bit_num;
					
					fNL0GammaPatch++;
					
					fL0GammaPatchIndex[fNL0GammaPatch-1] = (((idx << 5) & 0x7E0) | (tru_num & 0x1F));
				}
			}
		}
	}
}

//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::DecodeL1JetPatchIndexes(const int i, UInt_t *word32, const int offset)
{
	//////////////////////////////////////////////////////////
	// index des L1 jet                                     //
	//////////////////////////////////////////////////////////
	
	int jetSize = 2 + (fFwVersion >> 16);
	
	for (Int_t jet_row = 0; jet_row < 12 - (jetSize - 1); jet_row++)
	{
		UInt_t currentrow = word32[offset + jet_row];
		
		for (Int_t jet_col = 0; jet_col < 15; jet_col++)
		{
			if (currentrow & (1 << jet_col))
			{
				fNL1JetPatch[i] = fNL1JetPatch[i] + 1;
				
				fL1JetPatchIndex[fNL1JetPatch[i] - 1][i] = ((jet_row << 8) & 0xFF00) | (jet_col & 0xFF);
			}
		}
	}
}

//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::DecodeL1GammaPatchIndexes(const int i, UInt_t *word32, const int offset)
{
	//////////////////////////////////////////////////////////
	// index des L1 gamma                                   //
	//////////////////////////////////////////////////////////
	
	unsigned short truL1indexes[32][8];
	
	// extraction from stream
	for (Int_t index=0;index<8;index++)
	{
		for (Int_t tru_num=0;tru_num<16;tru_num++)
		{
			truL1indexes[2*tru_num  ][index] = ( word32[offset + index * 16 + tru_num]        & 0xFFFF);
			truL1indexes[2*tru_num+1][index] = ((word32[offset + index * 16 + tru_num] >> 16) & 0xFFFF);
		}
	}	
	
	// interpretation
	int gammacolnum;
	short indexcopy;
	
	for (Int_t tru_num=0;tru_num<32;tru_num++)
	{
		for (Int_t index=0;index<8;index++)
		{
			for (Int_t bit_num=0; bit_num<12; bit_num++)
			{
				if ((truL1indexes[tru_num][index] & (1<<bit_num)) != 0)
				{
					if (index<4) // Even
					{
						gammacolnum = (2*bit_num  );
						indexcopy   = index;
					}
					else         // Odd
					{
						gammacolnum = (2*bit_num+1);
						indexcopy   = index-4;
					}						
					
					fNL1GammaPatch[i] = fNL1GammaPatch[i] + 1;
					
					fL1GammaPatchIndex[fNL1GammaPatch[i] - 1][i] = (((indexcopy << 10) & 0xC00) | ((gammacolnum << 5) & 0x3E0) | (tru_num & 0x1F));
				}
			}
		}
	}	
}

//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::DecodeTRUADC(UInt_t *word32, const int offset)
{
	// extraction from stream
	for (Int_t index=0;index<96;index++)
	{
		for (Int_t tru_num=0;tru_num<16;tru_num++)
		{
			fADC[2*tru_num  ][index] = ( word32[offset + index * 16 + tru_num]        & 0xFFFF);
			fADC[2*tru_num+1][index] = ((word32[offset + index * 16 + tru_num] >> 16) & 0xFFFF);
		}
	}	
	
	for (Int_t tru_num=16;tru_num<32;tru_num++) // A side
	{
		Int_t v[96];
		for (Int_t index=0;index<96;index++) v[index] = fADC[tru_num][95-index];
		
		for (Int_t index=0;index<96;index++) fADC[tru_num][index] = v[index];
	}
}


//_____________________________________________________________________________
Bool_t AliEMCALTriggerSTURawStream::GetL0GammaPatch(const Int_t i, Int_t& tru, Int_t& idx) const
{
	// L0 gamma patches sent to STU (original access to L0 patch indexes)
	
	if (i > fNL0GammaPatch) return kFALSE;
	
	tru =  fL0GammaPatchIndex[i] & 0x1F;
	idx = (fL0GammaPatchIndex[i] & 0x7E0) >> 5;
	
	return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliEMCALTriggerSTURawStream::GetL1GammaPatch(const Int_t i, const Int_t j, Int_t& tru, Int_t& col, Int_t& row) const
{
	// L1 gamma patch indexes
	
	if (j >= 2 || i > fNL1GammaPatch[j]) return kFALSE;
	
	tru =  fL1GammaPatchIndex[i][j] & 0x1F;
	col = (fL1GammaPatchIndex[i][j] & 0x3E0) >> 5;
	row = (fL1GammaPatchIndex[i][j] & 0xC00) >> 10;
	
	return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliEMCALTriggerSTURawStream::GetL1JetPatch(const Int_t i, const Int_t j, Int_t& col, Int_t& row) const
{
	// L1 jet patch indexes
	
	if (j >= 2 || i > fNL1JetPatch[j]) return kFALSE;
	
	col =  fL1JetPatchIndex[i][j] & 0xFF;
	row = (fL1JetPatchIndex[i][j] & 0xFF00) >> 8;
	
	return kTRUE;
}

//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::GetADC(Int_t iTRU, UInt_t ADC[])
{
	// Time sums
	
	for (Int_t i=0; i<96; i++) ADC[i] = fADC[iTRU][i];
}

//_____________________________________________________________________________
void AliEMCALTriggerSTURawStream::DumpPayLoad(const Option_t *option) const
{
	// Dump STU payload
	
	TString op = option;
	
	printf("V0A:             %d\n", fV0A);
	printf("V0C:             %d\n", fV0C);

	printf("RawData:         %d\n", fGetRawData);
	printf("RegionEnable:    %8x\n", fRegionEnable); 
	printf("FrameReceived:   %8x\n", fFrameReceived);
	printf("FwVersion:       %x\n", fFwVersion);
	printf("Number of L0:    %d\n", fNL0GammaPatch);

	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 3; j++) {
			printf("G[%d][%d]: %d\n", j, i, fG[j][i]);
			printf("J[%d][%d]: %d\n", j, i, fJ[j][i]);
		}
		
		printf("Gamma threshold[%d]: %d\n", i, fL1GammaThreshold[i]);
		printf("Jet Threshold[%d]:   %d\n", i, fL1JetThreshold[i]);

		printf("Number of L1-g[%d]: %d\n", i, fNL1GammaPatch[i]);
		printf("Number of L1-j[%d]: %d\n", i, fNL1JetPatch[i]);
	}
	
	Int_t itru, col, row;
	
	if (op.Contains("L0") || op.Contains("ALL"))
	{
		for (Int_t i = 0;i < fNL0GammaPatch; i++) {
			if (GetL0GammaPatch(i,itru,col)) 
				cout << "> Found L0 gamma in TRU #" << setw(2) << itru <<  " at idx: " << setw(2) << col << endl;
		}
	}
	
	if (op.Contains("L1") || op.Contains("ALL")) {
		for (int j = 0; j < 2; j++) {
			for (Int_t i = 0; i < fNL1GammaPatch[j]; i++) {
				if (GetL1GammaPatch(i, j, itru, col, row)) 
					cout << "> Found L1 gamma " << j << " in TRU #" << setw(2) << itru <<  " at: ( col: " << setw(2) << col << " , row: " << setw(2) << row << " )" << endl;
			}
			
			for (Int_t i = 0; i < fNL1JetPatch[j]; i++) {
				if (GetL1JetPatch(i, j, col, row)) cout << "> Found L1 jet " << j << " at: ( col: " << setw(2) << col << " , row: " << setw(2) << row << " )" << endl;
			}
		}
	}
	
	if ((op.Contains("ADC") || op.Contains("ALL")) && fGetRawData) {
		for (Int_t i = 0; i < 32; i++) {
			cout << "--------\n";
			cout << "TRU #" << setw(2) << i << ":";
			for (Int_t j = 0;j < 96; j++) { 
				//TBits xadc(12); xadc.Set(12,&fADC[i][j]); 
				if ((j % 4) == 0) cout << endl;
				//cout << setw(2) << j << ": " << xadc << " ";
				printf("%2d: %3x / ",j,fADC[i][j]); 
			}
			cout << "\n";
		}
	}
}

//_____________________________________________________________________________
UShort_t AliEMCALTriggerSTURawStream::GetThreshold(Short_t A, Short_t B, Short_t C, UShort_t V0A, UShort_t V0C) const
{
	// Get threshold 
  ULong64_t v0sum = V0A + V0C;
  
  ULong64_t sqrV0 = v0sum * v0sum;					
					
  sqrV0 *= A;
					
  sqrV0 >>= 32;
				
  v0sum *= B;
					
  v0sum >>= 16;
					
  return (UShort_t)(sqrV0 + v0sum + C);
}
