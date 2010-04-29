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
#include "AliEMCALDigit.h"
#include "AliEMCALTriggerSTU.h"
#include "AliEMCALCalibData.h"
#include "AliLog.h"

#include <TF1.h>
#include <TMath.h>
#include <TClonesArray.h>
#include <TSystem.h>
#include <Riostream.h>
  
namespace
{
	const Int_t kTimeBins       = 64; // number of sampling bins of the FastOR signal
	const Int_t kTimeWindowSize =  4; // 
	const Int_t kNup            =  2; // 
	const Int_t kNdown          =  1; // 
}

ClassImp(AliEMCALTriggerTRU)

//________________
AliEMCALTriggerTRU::AliEMCALTriggerTRU() : AliEMCALTriggerBoard()//,
//fDigits(    0x0 )
{
	//
	for (Int_t i=0;i<96;i++) for (Int_t j=0;j<256;j++) fADC[i][j] = 0;
}

//________________
AliEMCALTriggerTRU::AliEMCALTriggerTRU(AliEMCALCalibData *calibData, const TVector2& rSize, Int_t mapType) : 
AliEMCALTriggerBoard(calibData, rSize)
{
	//
	for (Int_t i=0;i<96;i++) for (Int_t j=0;j<256;j++) fADC[i][j] = 0;

	// FIXME: use of class AliEMCALTriggerParam to get size
	TVector2 size;
	
	size.Set( 1. , 1. );
	SetSubRegionSize( size ); // 1 by 1 FOR
	
	size.Set( 2. , 2. );
	SetPatchSize( size );     // 2 by 2 subregions
	
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
   // delete TRU digits only used as transient containers 
   // to compute FastOR from energy deposit

}

//________________
void AliEMCALTriggerTRU::Peaks(Int_t arr[96][2])
{
	// Return max time bin & max value for all channels
	for (Int_t i=0;i<96;i++)
	{
		arr[i][0] = arr[i][1] = 0;
		
		Int_t max = 0, pos = 0;
		for (Int_t j=0;j<256;j++)
		{
			if (fADC[i][j]>max) 
			{
				max = fADC[i][j];
				pos = j;
			}
		}
		
		arr[i][0] = max;
		arr[i][1] = pos;
	}
}

//________________
void AliEMCALTriggerTRU::ShowFastOR(Int_t iTimeWindow, Int_t iChannel)
{
	//
	Int_t iChanF, iChanL;
	
	if (iChannel != -1) iChanF = iChanL = iChannel;
	else
	{
		iChanF =  0;
		iChanL = 96;
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
Int_t AliEMCALTriggerTRU::L0v0()
{
	// Mimick the TRU L0 'virtual' since not yet released algo
	
	// L0 issuing condition is: (2 up & 1 down) AND (time sum > thres)
	// fill a matrix to support sliding window
	// compute the time sum for all the FastOR of a given TRU
	// and then move the space window
	
	Int_t sum[96][3] ;//= { 0 };
	for(Int_t i = 0; i < 96 ; i++)
		for(Int_t j = 0; j < 3 ; j++) sum[i][j] = 0;
	
	// Sliding window algorithm
	for (Int_t i=0; i<=(kTimeBins-kTimeWindowSize); i++) 
	{
		for(Int_t j=0; j<fRegionSize->X(); j++)
		{		
			for (Int_t k=0; k<fRegionSize->Y(); k++)
			{
				for (Int_t l=i; l<i+kTimeWindowSize; l++) 
				{
					// [eta][phi][time]
					fRegion[j][k] += fADC[fMap[j][k]][l];
				}
				
//				fRegion[j][k] = fRegion[j][k]>>2; // truncate time sum
			}
		}
				
		// Threshold 
		// FIXME: for now consider just one threshold for all patches, should consider one per patch?
		// return the list of patches above threshold
		// Should probably be checked out from OCDB
		
		Int_t vL0Threshold = 0;
		
		SlidingWindow( kGamma, vL0Threshold );
		
		Int_t nP = 0;
		
		for (Int_t j=0; j<fPatches->GetEntriesFast(); j++)
		{
			AliEMCALTriggerPatch *p = (AliEMCALTriggerPatch*)fPatches->At( j );

			TVector2 v;
			p->Position(v);
			
			Int_t idx = fMap[int(v.X())][int(v.Y())];
			
			if ( i>2 ) 
			{				
				// Now check the '2 up/1 down' on each patch
				if ( sum[idx][1]>sum[idx][0] && sum[idx][2]<=sum[idx][1] ) nP++;
				
				sum[idx][0] = sum[idx][1];
				sum[idx][1] = sum[idx][2];
				sum[idx][2] = p->Sum();
			}
			else
			{
				sum[idx][i] = p->Sum();
			}
		}
		
		if ( !nP ) 
			fPatches->Delete();
		else
			break;     // Stop the algo when at least one patch is found ( thres & max )
		
		ZeroRegion();  // Clear fRegion for this time window before computing the next one
		
	}
	
	return fPatches->GetEntriesFast();
}

//________________
Int_t AliEMCALTriggerTRU::L0v1()
{
	// Mimick the TRU L0 'virtual' since not yet released algo
	
	// L0 issuing condition is: (2 up & 1 down) AND (time sum > thres)
	// fill a matrix to support sliding window
	// compute the time sum for all the FastOR of a given TRU
	// and then move the space window

	AliDebug(1,"=== Running TRU L0 v1 version ===");
	
	// Time sliding window algorithm
	for (Int_t i=0; i<=(kTimeBins-kTimeWindowSize); i++) 
	{
		AliDebug(1,Form("----------- Time window: %d\n",i));
		
		for (Int_t j=0; j<fRegionSize->X(); j++)
		{		
			for (Int_t k=0; k<fRegionSize->Y(); k++)
			{
				for (Int_t l=i; l<i+kTimeWindowSize; l++) 
				{
					// [eta][phi][time]
					fRegion[j][k] += fADC[fMap[j][k]][l];
				}
				
//				if (kTimeWindowSize > 4) fRegion[j][k] = fRegion[j][k] >> 1; // truncate time sum to fit 14b
			}
		}
		
		// Threshold 
		// FIXME: for now consider just one threshold for all patches, should consider one per patch?
		// ANSWE: both solutions will be implemented in the TRU
		// return the list of patches above threshold
		// Should probably be checked out from OCDB
		
		Int_t vL0Threshold = 0;
		
		SlidingWindow( kGamma, vL0Threshold );
		
//		for(Int_t j=0; j<fRegionSize->X(); j++)
//			for (Int_t k=0; k<fRegionSize->Y(); k++) fRegion[j][k] = fRegion[j][k]>>2; // go to 12b before shipping to STU
		
		Int_t nP = 0;
		
		for (Int_t j=0; j<fPatches->GetEntriesFast(); j++)
		{
			AliEMCALTriggerPatch* p = (AliEMCALTriggerPatch*)fPatches->At( j );

			if ( AliDebugLevel() ) p->Print("");

			TVector2 v;
			p->Position(v);
			
			Int_t sizeX = (Int_t)(fPatchSize->X() * fSubRegionSize->X());
			Int_t sizeY = (Int_t)(fPatchSize->Y() * fSubRegionSize->Y());
			
			const Int_t psize =  sizeX * sizeY; // Number of FastOR in the patch
			
			Int_t *idx= new Int_t[psize];
			
			Int_t aPeaks = 0;
			
			for (Int_t xx=0;xx<sizeX;xx++) 
			{
				for (Int_t yy=0;yy<sizeY;yy++) 
				{   
					idx[xx*sizeY+yy] = fMap[int(v.X()*fSubRegionSize->X())+xx][int(v.Y()*fSubRegionSize->Y())+yy]; // Get current patch FastOR ADC channels 
					
					if (fRegion[int(v.X()*fSubRegionSize->X())+xx][int(v.Y()*fSubRegionSize->Y())+yy]) aPeaks++;
					
					if ( AliDebugLevel() ) ShowFastOR(i,idx[xx*sizeY+yy]);
				}
			}

			Int_t nPeaks = 0;
			
			for (Int_t k=i;k<=i+kTimeWindowSize-(kNup+kNdown);k++)
			{				
				// Now check the 'kNup up / kNdown down' on each FastOR of the patch
				PeakFinder( idx , psize , k , kNup , kNdown , nPeaks );
			}
			
			if (nPeaks == aPeaks) 
			{
				if ( AliDebugLevel() ) 
				{
					printf("\t----- Valid patch (all FastOR have crossed a maximum)\n");
				  for (Int_t xx=0;xx<sizeX;xx++) {
				    for (Int_t yy=0;yy<sizeY;yy++) {
				      Int_t index = xx*sizeY+yy;
				      ShowFastOR(i,idx[index]); 
				    }
				  }
				}
				
				nP++; // all FOR in the patch must have seen a max
			}
			
			delete [] idx;
		}
		
		if ( !nP ) 
			fPatches->Delete();
		else
		{
			AliDebug(1,Form("==========[ Found %4d valid patches out of %4d ]==========\n",nP,fPatches->GetEntriesFast()));
			break;     // Stop the algo when at least one patch is found ( thres & max )
		}

		ZeroRegion();  // Clear fRegion for this time window before computing the next one		
	}
	
	return fPatches->GetEntriesFast();
}


//________________
Int_t AliEMCALTriggerTRU::L0v2()
{
	// Activity trigger

	// Sliding window algorithm

	for(Int_t j=0; j<fRegionSize->X(); j++)
	{		
		for (Int_t k=0; k<fRegionSize->Y(); k++)
		{
			Int_t max = 0;
			for (Int_t l=0; l<kTimeBins; l++) 
			{
				if (fADC[fMap[j][k]][l] > max) max = fADC[fMap[j][k]][l];
			}	
			
			if (max>4) fRegion[j][k] = max;
		}
	}
		
	Int_t vL0Threshold = 0;
		
	SlidingWindow( kGamma, vL0Threshold );
	
	return fPatches->GetEntriesFast();
}


//________________
void AliEMCALTriggerTRU::SetADC( Int_t channel, Int_t bin, Int_t sig )
{
	//
	if (channel>95) AliError("TRU has 96 ADC channels only!");
	fADC[channel][bin] = sig;
}

//________________
void AliEMCALTriggerTRU::PeakFinder( const Int_t idx[], Int_t nfastor, Int_t start, Int_t nup, Int_t ndown, Int_t& nPeaks ) 
{
	//
	for (Int_t i=0;i<nfastor;i++)
	{ 
		Int_t foundU = 0;
		Int_t foundD = 0;
		
		for (Int_t j=start+  1;j<start+nup      ;j++) foundU = ( fADC[idx[i]][j]> fADC[idx[i]][j-1] && fADC[idx[i]][j-1] ) ? 1 : 0;
		for (Int_t j=start+nup;j<start+nup+ndown;j++) foundD = ( fADC[idx[i]][j]<=fADC[idx[i]][j-1] && fADC[idx[i]][j  ] ) ? 1 : 0; 
		
		if ( foundU && foundD ) nPeaks++;
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

/*
//________________
void AliEMCALTriggerTRU::Scan()
{
	//
	for (Int_t i=0;i<96;i++) 
	{
		Int_t ietam = 23 - i/4;
		
		Int_t iphim =  3 - i%4;
		
		printf("ADC: %2d fRegion[%2d][%2d]: %4d\n",i,ietam,iphim,fRegion[ietam][iphim]);
	}	
}	
*/
//________________
void AliEMCALTriggerTRU::Reset()
{
	//
	fPatches->Delete();
	
	ZeroRegion();
	
	for (Int_t i=0;i<96;i++) for (Int_t j=0;j<256;j++) fADC[i][j] = 0;
}

