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

#include "AliEMCALTriggerSTU.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliEMCALTriggerSTUDCSConfig.h"
#include "AliVZEROCalibData.h"
#include "AliVZEROdigit.h"
#include "AliEMCALTriggerPatch.h"
#include "AliESDVZERO.h"
#include "AliLog.h"

#include <TClonesArray.h>
#include <TSystem.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>

#include <fstream>
#include <Riostream.h>
#include <cstdlib>

ClassImp(AliEMCALTriggerSTU)

//_______________
AliEMCALTriggerSTU::AliEMCALTriggerSTU() : AliEMCALTriggerBoard()
,fDCSConfig(0x0)
{
	// Ctor
	fGammaTh[0] = fGammaTh[1] = 0;
	fJetTh[0] = fJetTh[1] = 0;
}

//_______________
AliEMCALTriggerSTU::AliEMCALTriggerSTU(AliEMCALTriggerSTUDCSConfig *dcsConf, const TVector2& RS) : AliEMCALTriggerBoard(RS)
,fDCSConfig(dcsConf)
{
	// Ctor
	fGammaTh[0] = fGammaTh[1] = 0;
	fJetTh[0] = fJetTh[1] = 0;	
}

//_______________
AliEMCALTriggerSTU::~AliEMCALTriggerSTU()
{
	// Dtor
}

//_______________
void AliEMCALTriggerSTU::Build( TString& str, Int_t iTRU, Int_t** M, const TVector2* rSize )
{
	// Build
	
	str.ToLower();
	
	Int_t ix = (iTRU % 2) ? 24 : 0;

	Int_t iy =  iTRU / 2;
	
	Int_t** v = 0x0;
	
	if (str.Contains("map"))
	{
		v = fMap;
	}
	else if (str.Contains("region"))
	{
		v = fRegion;
	}
	else
	{
		AliError("Operation not allowed: STU won't be configured properly!");
	}
  
  if(v){	
    for (Int_t i=0; i<rSize->X(); i++)
      for (Int_t j=0; j<rSize->Y(); j++) v[i + ix][j + iy * 4] = M[i][j];
  }
}

//_______________
void AliEMCALTriggerSTU::L1(int type)
{
	// L1
	
	TVector2 s1, s2, s3, s4;
	fDCSConfig->GetSegmentation(s1, s2, s3, s4);
	
	switch (type)
	{
		case kL1GammaHigh:
		case kL1GammaLow:
			SetSubRegionSize(s1); 
			SetPatchSize(s2);
			break;
		case kL1JetHigh:
		case kL1JetLow:
			SetSubRegionSize(s3);
			SetPatchSize(s4);
			break;
		default:
			AliError("Not supported L1 trigger type");
			return;
			break;
	}
	
	SlidingWindow(GetThreshold(type));	
	AliDebug(999, Form("STU type %d sliding window w/ thr %d found %d patches", type, GetThreshold(type), fPatches->GetEntriesFast()));
}

//___________
void AliEMCALTriggerSTU::ComputeThFromV0(int type, const Int_t M[])
{
	// Compute threshold from V0
	
	Short_t P[3] = {0};
	
	switch (type)
	{
		case kL1GammaHigh:
			P[0] = fDCSConfig->GetG(0, 0);
			P[1] = fDCSConfig->GetG(1, 0);
			P[2] = fDCSConfig->GetG(2, 0);			
			break;
		case kL1GammaLow:
			P[0] = fDCSConfig->GetG(0, 1);
			P[1] = fDCSConfig->GetG(1, 1);
			P[2] = fDCSConfig->GetG(2, 1);			
			break;
		case kL1JetHigh:
			P[0] = fDCSConfig->GetJ(0, 0);
			P[1] = fDCSConfig->GetJ(1, 0);
			P[2] = fDCSConfig->GetJ(2, 0);			
			break;
		case kL1JetLow:
			P[0] = fDCSConfig->GetJ(0, 1);
			P[1] = fDCSConfig->GetJ(1, 1);
			P[2] = fDCSConfig->GetJ(2, 1);			
			break;
		default:
			AliError("AliEMCALTriggerSTU::ComputeThFromV0(): Undefined trigger type, pls check!");
			return;
	}
	
	ULong64_t v0sum = M[0] + M[1];
	
	ULong64_t sqrV0 = v0sum * v0sum;
	
	sqrV0 *= P[0];	
	
	sqrV0 >>= 32;
	
	v0sum *= P[1];
	
	v0sum >>= 16;
	
	SetThreshold(type, (UShort_t)(sqrV0 + v0sum + P[2]));
}

//___________
void AliEMCALTriggerSTU::SetThreshold(int type, Int_t v)
{
	// Set threshold
	
	switch (type)
	{
		case kL1GammaHigh:
			fGammaTh[0] = v;
			break;
		case kL1GammaLow:
			fGammaTh[1] = v;
			break;
		case kL1JetHigh:
			fJetTh[0] = v;		
			break;
		case kL1JetLow:
			fJetTh[1] = v;		
			break;
		default:
			AliError("AliEMCALTriggerSTU::SetThreshold(): Undefined trigger type, pls check!");
	}
}

//___________
Int_t AliEMCALTriggerSTU::GetThreshold(int type)
{	
	// Compute threshold FIXME: need an access to the OCDB
	// to get f(V0) parameters depending on trigger type
	
	switch (type)
	{
		case kL1GammaHigh:
			return fGammaTh[0];
			break;
		case kL1GammaLow:
			return fGammaTh[1];
			break;
		case kL1JetHigh:
			return fJetTh[0];		
			break;
		case kL1JetLow:
			return fJetTh[1];		
			break;
		default:
			AliError("AliEMCALTriggerSTU::GetThreshold(): Undefined trigger type, pls check!");
	}
	
	return 0;
}

//__________
void AliEMCALTriggerSTU::Reset()
{
	// Reset
	
	fPatches->Delete();
}
