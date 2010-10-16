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
,fGammaTh(0)
,fJetTh(0)
,fDCSConfig(0x0)
{
	//
}

//_______________
AliEMCALTriggerSTU::AliEMCALTriggerSTU(AliEMCALTriggerSTUDCSConfig *dcsConf, const TVector2& RS) : AliEMCALTriggerBoard(RS)
,fGammaTh(0)
,fJetTh(0)
,fDCSConfig(dcsConf)
{
	//
}

//_______________
AliEMCALTriggerSTU::~AliEMCALTriggerSTU()
{
	//
}

//_______________
Int_t AliEMCALTriggerSTU::GetRawData() const
{
	//
	return fDCSConfig->GetRawData();
}

//_______________
void AliEMCALTriggerSTU::Build( TString& str, Int_t iTRU, Int_t** M, const TVector2* rSize )
{
	//
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
void AliEMCALTriggerSTU::L1(TriggerType_t type)
{
	//
	TVector2 s1, s2, s3, s4;
	fDCSConfig->GetSegmentation(s1, s2, s3, s4);
	
	switch (type)
	{
		case kL1Gamma:
			SetSubRegionSize(s1); 
			SetPatchSize(s2);
			break;
		case kL1Jet:
			SetSubRegionSize(s3);
			SetPatchSize(s4);
			break;
		default:
			AliError("Not supported L1 trigger type");
			return;
			break;
	}
	
	SlidingWindow(type, GetThreshold(type));	
}

//________________
void AliEMCALTriggerSTU::PrintADC( TriggerType_t type, TVector2& pos, TVector2& idx )
{
	//
	Int_t ix = (Int_t) (( pos.X() + fPatchSize->X() ) * fSubRegionSize->X());
	
	Int_t iy = (Int_t) (( pos.Y() + fPatchSize->Y() ) * fSubRegionSize->Y());
	
	TString subRegionADC[] = {"0->15", "16->31", "32->47", "48->63", "64->79", "80->95"};
	
	switch ( type )
	{
		case kL1Gamma:
		{
			Int_t iTRU = ( (ix-1) < 24 ) ? 31 - int(pos.Y() / 4) : 15 - int(pos.Y() / 4);
			
			printf("TRU #%d row #%d col #%d fastor: ",iTRU,int(idx.Y()),int(idx.X()));
			
			for (Int_t i=(Int_t)(pos.X() * fSubRegionSize->X());i<ix;i++) 
			{
			  for (Int_t j=(Int_t) (pos.Y() * fSubRegionSize->Y());j<iy;j++) 
				{
					Int_t jtru = ( i < 24 ) ? 31 - j / 4 : 15 - j / 4;
					printf("TRU#%d_%d ",jtru,fMap[i][j]);
				}
			}
			
			cout << endl;
		}	
		break;
		case kL1Jet:
		{
			//Int_t jTRU = ( (ix-1) < 24 ) ? 31 - (iy-1) / 4 : 15 - (iy-1) / 4;
			
			printf("jet found at row : %d and col : %d",int(idx.X()),int(idx.Y()));
			
			Char_t vPair[100];
			Int_t nSubRegion = 0;
			
			for (Int_t i=(Int_t)(pos.X() * fSubRegionSize->X());i<ix;i++) 
			{
			  for (Int_t j=(Int_t)(pos.Y() * fSubRegionSize->Y());j<iy;j++) 
				{
					Int_t itru = ( i < 24 ) ? 31 - j / 4 : 15 - j / 4;
					
					Int_t idSR = fMap[i][j]/16;
					
					Char_t value = ((itru << 3) & 0xF8) | (idSR & 0x7);
					
					Bool_t isFound = kFALSE;
					
					for (Int_t k=0;k<nSubRegion && nSubRegion;k++)
					{
						if (vPair[k] == value) isFound = kTRUE;
					}
					
					if (!isFound) 
					{	
						nSubRegion++;
						
						vPair[nSubRegion-1] = value;
					}
				}
			}
			
			cout << " fastor:"; 
			for (Int_t i=0;i<nSubRegion;i++)
			{
				cout << " TRU#" << ((vPair[i] & 0xF8) >> 3) << "_" << subRegionADC[(vPair[i] & 0x7)];
			}
		
			cout << endl;
		}
		break;
		default:
			AliError("AliEMCALTriggerSTU::PrintADC(): Undefined trigger type, pls check!");
	}
}

//___________
void AliEMCALTriggerSTU::PatchGenerator(const TClonesArray* lpos, Int_t val)
{
	ZeroRegion();
	
	Int_t vTRU[32][24][4] ;//= {0};
	
	for (Int_t i = 0; i < 32; i++) 
		for (Int_t j = 0; j < 24; j++) 
			for (Int_t k = 0; k < 4; k++) vTRU[i][j][k]=0;
	
	
	AliWarning("AliEMCALTriggerSTU::PatchGenerator(): STU region has been reset!");
	
	// Fill the patch FOR at 'pos' w/ value 'val'
	
	TIter NextPosition( lpos );
	
	while ( TVector2 *pos = (TVector2*)NextPosition() )
	{
		pos->Print();
		for (Int_t i=(Int_t)(pos->X()*fSubRegionSize->X()); i<int((pos->X()+fPatchSize->X())*fSubRegionSize->X()); i++)
		  for (Int_t j=(Int_t)(pos->Y()*fSubRegionSize->Y()); j<int((pos->Y()+fPatchSize->Y())*fSubRegionSize->Y()); j++) 
				fRegion[i][j] = val;
	}
	
	for (Int_t i=0; i<2; i++)
	{
		for (Int_t j=0; j<16; j++)
		{
			Int_t iTRU = ( !i ) ? 31 - j : 15 - j;
			
			for (Int_t i1=24*i; i1<24*(i+1); i1++)
			{
				for (Int_t j1=4*j; j1<4*(j+1); j1++)
				{
					vTRU[iTRU][i1%24][j1%4] = fRegion[i1][j1];
				}
			}
		}
	}
	
	gSystem->Exec(Form("mkdir -p GP"));

	for (Int_t i=0; i<32; i++)
	{
		std::ofstream outfile(Form("GP/data_TRU%d.txt",i),std::ios::trunc);
	
		for (Int_t j=0;j<96;j++) 
		{
			Int_t ietam, iphim;
			
			if ( int(i / 16) )
			{
				ietam =      j/4;
				iphim =  3 - j%4;
			}
			else
			{					
				ietam = 23 - j/4;
				iphim =      j%4;
			}
			
			outfile << vTRU[i][ietam][iphim] << endl;
		}
	
		outfile.close();
	}
}

//___________
void AliEMCALTriggerSTU::ComputeThFromV0(const Int_t M[])
{
	//
	if (!(M[0] + M[1])) AliWarning("V0A + V0C is null!"); // 0/1: V0C/V0A

	
	fGammaTh = fDCSConfig->GetGA()*(M[0] + M[1])*(M[0] + M[1]) + fDCSConfig->GetGB()*(M[0] + M[1]) + fDCSConfig->GetGC();
	
	fJetTh   = fDCSConfig->GetJA()*(M[0] + M[1])*(M[0] + M[1]) + fDCSConfig->GetJB()*(M[0] + M[1]) + fDCSConfig->GetJC();	
}

//___________
void AliEMCALTriggerSTU::SetThreshold(TriggerType_t type, Int_t v)
{
	switch (type)
	{
		case kL1Gamma:
			fGammaTh = v;
			break;
		case kL1Jet:
			fJetTh = v;		
			break;
		default:
			AliError("AliEMCALTriggerSTU::SetThreshold(): Undefined trigger type, pls check!");
	}
}

//___________
Int_t AliEMCALTriggerSTU::GetThreshold(TriggerType_t type)
{	
	// Compute threshold FIXME: need an access to the OCDB
	// to get f(V0) parameters depending on trigger type
	
	
	switch (type)
	{
		case kL1Gamma:
			return fGammaTh;
			break;
		case kL1Jet:
			return fJetTh;		
			break;
		default:
			AliError("AliEMCALTriggerSTU::GetThreshold(): Undefined trigger type, pls check!");
	}
	
	return 0;
}

//__________
void AliEMCALTriggerSTU::Reset()
{
	fPatches->Delete();
}
