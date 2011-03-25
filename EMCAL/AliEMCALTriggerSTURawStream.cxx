
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
	const Int_t kPayLoadSize            =  944;
	const Int_t kPayLoadSizeWithRawData = 1772;
}

ClassImp(AliEMCALTriggerSTURawStream)

//_____________________________________________________________________________
AliEMCALTriggerSTURawStream::AliEMCALTriggerSTURawStream() : TObject(),
fRawReader(0x0),
fL1JetThreshold(0),
fL1GammaThreshold(0),
fL0GammaPatchIndex(),
fL1GammaPatchIndex(),
fL1JetPatchIndex(),
fNL0GammaPatch(0),
fNL1JetPatch(0),
fNL1GammaPatch(0),
fGetRawData(0)
{
	//
}

//_____________________________________________________________________________
AliEMCALTriggerSTURawStream::AliEMCALTriggerSTURawStream(AliRawReader* rawReader) : TObject(),
fRawReader(rawReader),
fL1JetThreshold(0),
fL1GammaThreshold(0),
fL0GammaPatchIndex(),
fL1GammaPatchIndex(),
fL1JetPatchIndex(),
fNL0GammaPatch(0),
fNL1JetPatch(0),
fNL1GammaPatch(0),
fGetRawData(0)
{
	//
	fRawReader->Reset();
	fRawReader->Select("EMCAL",44);
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
	fNL1GammaPatch = 0;
	fNL1JetPatch   = 0;	
}

//_____________________________________________________________________________
Bool_t AliEMCALTriggerSTURawStream::ReadPayLoad()
{
	// STU data decoder from Olivier Bourrion LPSC CNRS-IN2P3
	// bourrion_at_lpsc_dot_in2p3_dot_fr
	
	UInt_t word32[1772]; // 32b words
	for (Int_t i=0;i<1772;i++) word32[i] = 0;
	
	Int_t iword = 0;
	
	fNL0GammaPatch = 0;
	fNL1GammaPatch = 0;
	fNL1JetPatch   = 0;
	
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
	
	if (iword != kPayLoadSize && iword != kPayLoadSizeWithRawData)
	{
		AliError(Form("ERROR: STU payload (eqId: %d, eqSize: %d) doesn't match expected size! %d word32",
					  eqId, eqSize, iword));
		return kFALSE;
	}
	else if (AliDebugLevel())
	{
		AliInfo(Form("STU (eqId: %d, eqSize: %d) payload size: %d word32",
					 eqId, eqSize, iword));
	}
	
	  fL1JetThreshold = ((word32[0]>>16) & 0xFFF);
	fL1GammaThreshold =  (word32[0]      & 0xFFF);
	
	for (Int_t jet_row=0; jet_row<11; jet_row++)
	{
		UInt_t currentrow = word32[1+jet_row];
		
		for (Int_t jet_col=0; jet_col<15; jet_col++)
		{
			if (currentrow & (1 << jet_col))
			{
				fNL1JetPatch++;
				
				fL1JetPatchIndex[fNL1JetPatch-1] = ((jet_row << 8) & 0xFF00) | (jet_col & 0xFF);
			}
		}
	}
	
	//////////////////////////////////////////////////////////
	// index des L0                                         //
	//////////////////////////////////////////////////////////
	// FIXME: still not interpreted to be done with Jiri
	
	unsigned short truL0indexes[32][6];
	
	// extraction from stream
	for (Int_t index=0;index<6;index++)
	{
		for (Int_t tru_num=0;tru_num<16;tru_num++)
		{
			truL0indexes[2*tru_num  ][index] = ( word32[12+index*16+tru_num]      & 0xFFFF);
			truL0indexes[2*tru_num+1][index] = ((word32[12+index*16+tru_num]>>16) & 0xFFFF);
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

	//////////////////////////////////////////////////////////
	// index des L1 gamma                                   //
	//////////////////////////////////////////////////////////
	
	unsigned short truL1indexes[32][8];
	
	// extraction from stream
	for (Int_t index=0;index<8;index++)
	{
		for (Int_t tru_num=0;tru_num<16;tru_num++)
		{
			truL1indexes[2*tru_num  ][index] = ( word32[108+index*16+tru_num]      & 0xFFFF);
			truL1indexes[2*tru_num+1][index] = ((word32[108+index*16+tru_num]>>16) & 0xFFFF);
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
					
					fNL1GammaPatch++;
					
					fL1GammaPatchIndex[fNL1GammaPatch-1] = (((indexcopy << 10) & 0xC00) | ((gammacolnum << 5) & 0x3E0) | (tru_num & 0x1F));
				}
			}
		}
	}	

	//////////////////////////////////////////////////////////
	// raw output                                           //
	//////////////////////////////////////////////////////////
	
	if ( iword <= kPayLoadSize ) 
	{
		fGetRawData = 0;
		return kTRUE;
	}
	
	fGetRawData = 1;
	
	// extraction from stream
	for (Int_t index=0;index<96;index++)
	{
		for (Int_t tru_num=0;tru_num<16;tru_num++)
		{
			fADC[2*tru_num  ][index] = ( word32[236+index*16+tru_num]      & 0xFFFF);
			fADC[2*tru_num+1][index] = ((word32[236+index*16+tru_num]>>16) & 0xFFFF);
		}
	}	

	for (Int_t tru_num=16;tru_num<32;tru_num++) // A side
	{
		Int_t v[96];
		for (Int_t index=0;index<96;index++) v[index] = fADC[tru_num][95-index];
		
		for (Int_t index=0;index<96;index++) fADC[tru_num][index] = v[index];
	}
	
	return kTRUE;
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
Bool_t AliEMCALTriggerSTURawStream::GetL1GammaPatch(const Int_t i, Int_t& tru, Int_t& col, Int_t& row) const
{
	// L1 gamma patch indexes
	
	if (i > fNL1GammaPatch) return kFALSE;
	
	tru =  fL1GammaPatchIndex[i] & 0x1F;
	col = (fL1GammaPatchIndex[i] & 0x3E0) >> 5;
	row = (fL1GammaPatchIndex[i] & 0xC00) >> 10;
	
	return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliEMCALTriggerSTURawStream::GetL1JetPatch(const Int_t i, Int_t& col, Int_t& row) const
{
	// L1 jet patch indexes
	
	if (i > fNL1JetPatch) return kFALSE;
	
	col =  fL1JetPatchIndex[i] & 0xFF;
	row = (fL1JetPatchIndex[i] & 0xFF00) >> 8;
	
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
	
	cout << "Jet Threshold: " << fL1JetThreshold << " Gamma threshold: " << fL1GammaThreshold << endl;
	cout << "Number of L0:   " << fNL0GammaPatch << endl;
	cout << "Number of L1-g: " << fNL1GammaPatch << endl;
	cout << "Number of L1-j: " << fNL1JetPatch << endl;
	
	Int_t itru, col, row;
	
	if (op.Contains("L0") || op.Contains("ALL"))
	{
		for (Int_t i=0;i<fNL0GammaPatch;i++)
		{
			
			if (GetL0GammaPatch(i,itru,col)) 
				cout << "> Found L0 gamma in TRU #" << setw(2) << itru <<  " at idx: " << setw(2) << col << endl;
		}
	}
	
	if (op.Contains("L1") || op.Contains("ALL"))
	{
		for (Int_t i=0;i<fNL1GammaPatch;i++)
		{
			if (GetL1GammaPatch(i,itru,col,row)) 
				cout << "> Found L1 gamma in TRU #" << setw(2) << itru <<  " at: ( col: " << setw(2) << col << " , row: " << setw(2) << row << " )" << endl;
		}

		for (Int_t i=0;i<fNL1JetPatch;i++)
		{
			if (GetL1JetPatch(i,col,row)) cout << "> Found L1 jet at: ( col: " << setw(2) << col << " , row: " << setw(2) << row << " )" << endl;
		}
	}
	
	cout << "> RawData: " << fGetRawData << endl;
	
	if ( (op.Contains("ADC") || op.Contains("ALL")) && fGetRawData )
	{
		for (Int_t i=0;i<32;i++)
		{
			cout << "--------\n";
			cout << "TRU #" << setw(2) << i << ":";
			for (Int_t j=0;j<96;j++) 
			{ 
				TBits xadc(12); xadc.Set(12,&fADC[i][j]); 
				if ((j%4)==0) cout << endl;
				//cout << setw(2) << j << ": " << xadc << " ";
				printf("%2d: %3x / ",j,fADC[i][j]); 
			}
			cout << "\n";
		}
	}
}
