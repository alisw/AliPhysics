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


#include "AliEMCALTriggerTRU.h"
#include "AliEMCALTriggerPatch.h"
#include "AliEMCALTriggerTRUDCSConfig.h"
#include "AliLog.h"

#include <TClonesArray.h>
#include <TSystem.h>
#include <Riostream.h>
#include <TFile.h>
#include <TROOT.h>

namespace
{
	const Int_t kTimeBins       = 16; // number of sampling bins of the FastOR signal
	const Int_t kTimeWindowSize =  4; // 
}

using std::ofstream;
using std::endl;
using std::ios_base;
ClassImp(AliEMCALTriggerTRU)

//________________
AliEMCALTriggerTRU::AliEMCALTriggerTRU() : AliEMCALTriggerBoard(),
fDCSConfig(0x0),
fL0Time(0)
{
	// Ctor
	
	for (Int_t i=0;i<96;i++) for (Int_t j=0;j<256;j++) fADC[i][j] = 0;
}

//________________
AliEMCALTriggerTRU::AliEMCALTriggerTRU(AliEMCALTriggerTRUDCSConfig* dcsConf, const TVector2& rSize, Int_t mapType) : 
AliEMCALTriggerBoard(rSize),
fDCSConfig(dcsConf),
fL0Time(0)
{
	// Ctor
	
	for (Int_t i=0;i<96;i++) for (Int_t j=0;j<256;j++) fADC[i][j] = 0;
	
	TVector2 size;
	
	if (dcsConf->GetL0SEL() & 0x0001) // 4-by-4
	{
		size.Set( 1. , 1. );
		SetSubRegionSize( size );
		
		size.Set( 2. , 2. );
		SetPatchSize( size );
	}	
	else                              // 2-by-2
	{
		size.Set( 1. , 1. );
		SetSubRegionSize( size );
		
		size.Set( 1. , 1. );
		SetPatchSize( size );	
	}	
	
	for (Int_t ietam=0;ietam<24;ietam++)
	{
		for (Int_t iphim=0;iphim<4;iphim++)
		{
			// idx: 0..95 since iphim: 0..11 ietam: 0..23
			Int_t idx = ( !mapType ) ? ( 3 - iphim ) + ietam * 4 : iphim + (23 - ietam) * 4;	
	
			// Build a matrix used to get TRU digit id (ADC channel) from (eta,phi)|SM
			fMap[ietam][iphim] = idx; // [0..11][0..3] namely [eta][phi] in SM
		}
	}
}

//________________
AliEMCALTriggerTRU::~AliEMCALTriggerTRU()
{
	// Dtor
}

//________________
void AliEMCALTriggerTRU::ShowFastOR(Int_t iTimeWindow, Int_t iChannel)
{
	// Dump
	
	Int_t iChanF, iChanL;
	
	if (iChannel != -1) iChanF = iChanL = iChannel;
	else
	{
		iChanF =  0;
		iChanL = 95;
	}
	
	for (Int_t i=iChanF;i<iChanL+1;i++)
	{
		printf("\tChannel: %2d - ",i);
		for (Int_t j=0;j<60;j++) 
		{
			if (j == iTimeWindow)
				printf(" | %4d",fADC[i][j]);
			else if (j == iTimeWindow+kTimeWindowSize-1)
				printf(" %4d |",fADC[i][j]);
			else
				printf(" %4d",fADC[i][j]);
		}
		
		printf("\n");
	}
}

//________________
Int_t AliEMCALTriggerTRU::L0()
{
	// L0 algo depending on TRU fw version
	
	const Int_t xsize    = Int_t(fRegionSize->X());
	const Int_t ysize    = Int_t(fRegionSize->Y());

	Int_t asum = 0;
	for (Int_t j = 0; j < xsize; j++) {		
		for (Int_t k = 0; k < ysize; k++) {
			for (Int_t l = 0; l < kTimeBins; l++) {
				asum += fADC[fMap[j][k]][l];
			}
		}
	}	
	
	// TRU has no signal, return!
	if (!asum) {
		AliDebug(999,"=== TRU has no signal ===");
		return 0;
	}
	
	AliDebug(999,Form("=== TRU PF: %x",fDCSConfig->GetSELPF()));
 
	UInt_t ma = fDCSConfig->GetSELPF() & 0xffff;
	
  // Set default peak finder if null
  if (!ma) ma = 0x1e1f;
  
	int nb = ma & 0x7f;
	
	ma = (ma >> 8) & 0x7f;
	
	AliDebug(999,Form("=== TRU fw version %x ===",fDCSConfig->GetFw()));
    
  if (fDCSConfig->GetFw() < 0x4d) {
      return L0v0(nb, ma);
    } else {
      return L0v1(nb, ma);
    }
}

//________________
Int_t AliEMCALTriggerTRU::L0v0(int mask, int pattern)
{
	// L0 issuing condition is: (2x2 PF) AND (4x4 > thres)
	
	AliDebug(999,"=== Running TRU L0 algorithm version 0 ===");

	const Int_t xsize    = Int_t(fRegionSize->X());
	const Int_t ysize    = Int_t(fRegionSize->Y());

	Int_t **othr = new Int_t*[xsize];
	Int_t **patt = new Int_t*[xsize];	
	Int_t **buff = new Int_t*[xsize];	
	
	for (Int_t x = 0; x < xsize; x++) {
		othr[x] = new Int_t[ysize];
		patt[x] = new Int_t[ysize];
		buff[x] = new Int_t[ysize];
	}
	
	for (Int_t i = 0; i < xsize; i++) {
		for (Int_t j = 0; j < ysize; j++) {
			othr[i][j] = 0;	
			patt[i][j] = 0;
			buff[i][j] = 0;
		}
	}
			
	// Time sliding window algorithm
	for (int i = 0; i <= (kTimeBins - kTimeWindowSize); i++) 
	{
		AliDebug(999,Form("----------- Time window: %d\n",i));
		
		if (AliDebugLevel()) ShowFastOR(i, -1);
		
		for (int j = 0; j < xsize; j++) {		
			for (int k = 0; k < ysize; k++) {
				
//				if (
//					!(j % int(fSubRegionSize->X())) 
//					&& 
//					!(k % int(fSubRegionSize->Y())) 
//					&& 
//					(j + int(fPatchSize->X() * fSubRegionSize->X()) <= xsize)
//					&& 
//					(k + int(fPatchSize->Y() * fSubRegionSize->Y()) <= ysize)
//					) 
//				{
//					int sum = 0;
//					
//					for (int l = 0; l < int(fPatchSize->X() * fSubRegionSize->X()); l++) 
//						for (int m = 0; m < int(fPatchSize->Y() * fSubRegionSize->Y()); m++) sum += fRegion[j + l][k + m];
//					
//					if (sum > int(fDCSConfig->GetGTHRL0())) {
//						AliDebug(999,Form("----------- Patch (%2d,%2d) is over threshold\n", j, k));
//						othr[j][k] = sum;
//					}
//				}
				
				buff[j][k] = fRegion[j][k];
				
				fRegion[j][k] = 0;
				for (Int_t l = i; l < i + kTimeWindowSize; l++) fRegion[j][k] += fADC[fMap[j][k]][l];	
				
				if (fRegion[j][k] > buff[j][k]) {
					patt[j][k] |= 0x1;
				}
				
				if (patt[j][k]) AliDebug(999,Form("----------- (%2d,%2d) New: %d Old: %d patt: %x / pattern: %x / mask: %x", j, k, fRegion[j][k], buff[j][k], patt[j][k], pattern, mask));
			}
		}
		
		for (int j = 0; j <= int(fRegionSize->X() - fPatchSize->X() * fSubRegionSize->X()); j += int(fSubRegionSize->X())) {
			for (int k = 0; k <= int(fRegionSize->Y() - fPatchSize->Y() * fSubRegionSize->Y()); k += int(fSubRegionSize->Y())) {
				
//				if (!othr[j][k]) continue;
				int sizeX = int(fPatchSize->X() * fSubRegionSize->X());				
				int sizeY = int(fPatchSize->Y() * fSubRegionSize->Y());
				
				int foundPeak = 0;
				int sum       = 0;
				
				for (int l = 0; l < sizeX; l++) {
					for (int m = 0; m < sizeY; m++) {
						sum += fRegion[j + l][k + m];
						
						if ((patt[j + l][k + m] & mask) == pattern) foundPeak++;
					}
				}
				
				if (sum > int(fDCSConfig->GetGTHRL0())) othr[j][k] = sum;
		
				if (foundPeak && othr[j][k]) {
					
					new((*fPatches)[fPatches->GetEntriesFast()]) AliEMCALTriggerPatch(j, k, othr[j][k], i);
					
					AliEMCALTriggerPatch* p = (AliEMCALTriggerPatch*)fPatches->At(fPatches->GetEntriesFast() - 1);
					
					if (AliDebugLevel()) p->Print("");
					
					const Int_t psize =  sizeX * sizeY; // Number of FastOR in the patch
					
					Int_t* idx = new Int_t[psize];
					
					for (Int_t l = 0; l < sizeX; l++) 
					{
						for (Int_t m = 0; m < sizeY; m++) 
						{   
							Int_t index = l * sizeY + m;
							
							idx[index] = fMap[int(j * fSubRegionSize->X()) + l][int(k * fSubRegionSize->Y()) + m];
							
              if ((patt[j + l][k + m] & mask) == (pattern & mask)) {
//								cout << "setting peak at " << l << " " << m << endl;
								p->SetPeak(l, m, sizeX, sizeY);
							}
							
							if (AliDebugLevel() >= 999) ShowFastOR(i, idx[index]);
						}
					}
					
					delete [] idx;
				}
			}
		}
		
		if (fPatches->GetEntriesFast() && !fL0Time) {			
			// Stop the algo when at least one patch is found ( thres & max )
			fL0Time = i;

//			break;
		}
		
		for (int j = 0; j < xsize; j++) 	
			for (int k = 0; k < ysize; k++) patt[j][k] <<= 1;
	}
	
	for (Int_t x = 0; x < xsize; x++) {
		delete [] othr[x];
		delete [] patt[x];
		delete [] buff[x];
	}

	delete [] othr;
	delete [] patt;
	delete [] buff;

	return fPatches->GetEntriesFast();
}

//________________
Int_t AliEMCALTriggerTRU::L0v1(int mask, int pattern)
{
	// L0 issuing condition is: (4x4 PF) AND (4x4 > thres)
	
	AliDebug(999,"=== Running TRU L0 algorithm version 1 ===");
	
	const Int_t xsize    = Int_t(fRegionSize->X());
	const Int_t ysize    = Int_t(fRegionSize->Y());
		
	Int_t **othr = new Int_t*[xsize];
	Int_t **buff = new Int_t*[xsize];
	Int_t **patt = new Int_t*[xsize];
	
	for (Int_t i = 0; i < xsize; i++) {
		buff[i] = new Int_t[ysize];
		patt[i] = new Int_t[ysize];
		othr[i] = new Int_t[ysize];
	}
	
	for (Int_t i = 0; i < xsize; i++) for (Int_t j = 0; j < ysize; j++) {
		othr[i][j] = 0;
		patt[i][j] = 0;
		buff[i][j] = 0;
	}
	
	// Time sliding window algorithm
	for (Int_t i = 0; i <= (kTimeBins - kTimeWindowSize); i++) {
		
		AliDebug(999,Form("----------- Time window: %d\n",i));
		
		for (int j = 0; j < xsize; j++) {		
			for (int k = 0; k < ysize; k++) {
				
//				if (
//					!(j % int(fSubRegionSize->X())) 
//					&& 
//					!(k % int(fSubRegionSize->Y())) 
//					&& 
//					(j + int(fPatchSize->X() * fSubRegionSize->X()) <= xsize)
//					&& 
//					(k + int(fPatchSize->Y() * fSubRegionSize->Y()) <= ysize)
//					) 
//				{
//					int sum = 0;
//					
//					for (int l = 0; l < int(fPatchSize->X() * fSubRegionSize->X()); l++) 
//						for (int m = 0; m < int(fPatchSize->Y() * fSubRegionSize->Y()); m++) sum += fRegion[j + l][k + m];
//					
//					if (sum > buff[j][k]) patt[j][k] |= 0x1;
//					
//					AliDebug(999,Form("----------- Patch (%2d,%2d) has sum %d while its whole time pattern is %x\n", j, k, sum, patt[j][k]));
//					
//					buff[j][k] = sum;
//					
//					if (sum > int(fDCSConfig->GetGTHRL0())) {
//						AliDebug(999,Form("----------- Patch (%2d,%2d) is over threshold\n", j, k));
//						othr[j][k] = sum;
//					}
//				}
				
				fRegion[j][k] = 0;
				for (Int_t l = i; l < i + kTimeWindowSize; l++) fRegion[j][k] += fADC[fMap[j][k]][l];	
			}
		}
		
		for (int j = 0; j <= int(fRegionSize->X() - fPatchSize->X() * fSubRegionSize->X()); j += int(fSubRegionSize->X())) {
			for (int k = 0; k <= int(fRegionSize->Y() - fPatchSize->Y() * fSubRegionSize->Y()); k += int(fSubRegionSize->Y())) {
					
				int sum = 0;
				
				for (int l = 0; l < int(fPatchSize->X() * fSubRegionSize->X()); l++) 
					for (int m = 0; m < int(fPatchSize->Y() * fSubRegionSize->Y()); m++) sum += fRegion[j + l][k + m];
				
				if (sum > buff[j][k]) patt[j][k] |= 0x1;
				
				if (sum > int(fDCSConfig->GetGTHRL0())) {
					AliDebug(999,Form("----------- Patch (%2d,%2d) is over threshold\n", j, k));
					
					othr[j][k] = sum;
				}
				
				AliDebug(999,Form("----------- Patch (%2d,%2d) has sum %d while its whole time pattern is %x\n", j, k, sum, patt[j][k]));
				
				buff[j][k] = sum;
				
        if (othr[j][k] && ((patt[j][k] & mask) == (pattern & mask))) {
					
					new((*fPatches)[fPatches->GetEntriesFast()]) AliEMCALTriggerPatch(j, k, othr[j][k], i);
					
//					AliDebug(999,Form("=== New L0 patch at (%2d,%2d) time: %2d",j, k, i));
					
					int sizeX = int(fPatchSize->X() * fSubRegionSize->X());					
					int sizeY = int(fPatchSize->Y() * fSubRegionSize->Y());
					
					for (int xx = 0; xx < sizeX; xx++) {
						for (int yy = 0; yy < sizeY; yy++) {
							((AliEMCALTriggerPatch*)fPatches->At(fPatches->GetEntriesFast() - 1))->SetPeak(xx, yy, sizeX, sizeY);
						}
					}
					
					AliEMCALTriggerPatch* p = (AliEMCALTriggerPatch*)fPatches->At(fPatches->GetEntriesFast() - 1);
					
					if (AliDebugLevel()) p->Print("");
				}
			}
		}
		
		if (fPatches->GetEntriesFast() && !fL0Time) {
			fL0Time = i;
			
//			break;
		} 
		
		for (int j = 0; j < xsize; j++) 	
			for (int k = 0; k < ysize; k++) patt[j][k] <<= 1;
	}
	
	for (Int_t x = 0; x < xsize; x++) {
		delete [] othr[x];
		delete [] patt[x];
		delete [] buff[x];
	}
	
	delete [] othr;
	delete [] patt;
	delete [] buff;
	
	return fPatches->GetEntriesFast();
}

//________________
void AliEMCALTriggerTRU::SetADC( Int_t channel, Int_t bin, Int_t sig )
{
  //Set ADC value
  if (channel > 95 || bin > 255) {
    AliError("TRU has 96 ADC channels and 256 bins only!");
  }
  else{ 
	if (((fDCSConfig->GetMaskReg(int(channel / 16)) >> (channel % 16)) & 0x1) == 0) fADC[channel][bin] = sig;
  }
}

//________________
void AliEMCALTriggerTRU::GetL0Region(const int time, Int_t arr[][4])
{
	Int_t r0 = time - fDCSConfig->GetRLBKSTU();
	
	if (r0 < 0) 
	{
		AliError(Form("TRU buffer not accessible! time: %d rollback: %d", time, fDCSConfig->GetRLBKSTU()));
		return;
	}
	
	for (Int_t i = 0; i < fRegionSize->X(); i++) 
	{
		for (Int_t j = 0; j < fRegionSize->Y(); j++) 
		{
			for (Int_t k = r0; k < r0 + kTimeWindowSize; k++)
			{
				arr[i][j] += fADC[fMap[i][j]][k];
			}
		}
	}
}

//________________
void AliEMCALTriggerTRU::SaveRegionADC(Int_t iTRU, Int_t iEvent)
{
	// O for STU Hw
	//
	gSystem->Exec(Form("mkdir -p Event%d",iEvent));
	
	ofstream outfile(Form("Event%d/data_TRU%d.txt",iEvent,iTRU),ios_base::trunc);
	
	for (Int_t i=0;i<96;i++) 
	{
		Int_t ietam = 23 - i/4;
	
		Int_t iphim =  3 - i%4;
		
		outfile << fRegion[ietam][iphim] << endl;
	}

	outfile.close();
}

//________________
void AliEMCALTriggerTRU::Reset()
{
	// Reset
	
	fPatches->Delete();
	
	ZeroRegion();
	
	for (Int_t i=0;i<96;i++) for (Int_t j=0;j<256;j++) fADC[i][j] = 0;
}

