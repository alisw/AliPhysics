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

/// \cond CLASSIMP
ClassImp(AliEMCALTriggerSTU) ;
/// \endcond

///
/// Default constructor
//_______________
AliEMCALTriggerSTU::AliEMCALTriggerSTU() : AliEMCALTriggerBoard()
,fDCSConfig(0x0)
{
  fGammaTh[0] = fGammaTh[1] = 0;
  fJetTh[0] = fJetTh[1] = 0;
  fBkgRho = 0;
}

///
/// Constructor
//_______________
AliEMCALTriggerSTU::AliEMCALTriggerSTU(AliEMCALTriggerSTUDCSConfig *dcsConf, const TVector2& RS) : AliEMCALTriggerBoard(RS)
,fDCSConfig(dcsConf)
{
  fGammaTh[0] = fGammaTh[1] = 0;
  fJetTh[0] = fJetTh[1] = 0;	
  fBkgRho = 0;
}

///
/// Destructor
//_______________
AliEMCALTriggerSTU::~AliEMCALTriggerSTU()
{ 
  fBkgRho = 0;
}

///
/// Build
//_______________
void AliEMCALTriggerSTU::Build( TString& str, Int_t iTRU, Int_t** M, const TVector2* rSize, Int_t triggerMapping = 1)
{	
  str.ToLower();
  
  Int_t ix,iy;

  switch (triggerMapping + 3*(iTRU >= 32)) {
    case 1:
      // Run 1 Geometry
      ix =  iTRU / 2; // Round down to even.  
      ix = 4 * ix; // Every two TRUs add 4 modules
      iy = (iTRU % 2) ? 24 : 0; // if odd, start at iEta = 24
      break;
    case 2: 
      // Run 2 EMCAL Geometry
      ix = iTRU / 6; // ix = number of rows  above.
      ix = 12 * ix; // Each row adds 12 modules
      iy = 8 * (iTRU % 6);  
      if ((iTRU == 30) || (iTRU == 31)) iy = 24 * (iTRU % 2); // EMCAL 1/3 SM TRUs
      break;
    case 4:	
      AliFatal("EMCAL STU object tried to use run 1 trigger mapping with a DCAL TRU");
  	case 5:
      iTRU -= 32;
      // Run2 DCAL Geometry
      // If using iTRU corresponding to virtual TRU index:
      ix = iTRU / 6; // ix = number of rows above
      ix = ix * 12; // Each row adds 12 modules
      iy = iTRU % 6; // iy = 0,1,2,3,4,5
      iy = iy * 8 ; // each TRU adds 8. 

      if ((iTRU == 18) || (iTRU == 19)) iy = (iTRU % 2) * 24; // DCAL 1/3 SM TRUs
      break;
    default:
      AliFatal("EMCAL STU object tried an unknown combination of Trigger Mapping and TRU index");
  }


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
      for (Int_t j=0; j<rSize->Y(); j++) { 
				v[i + ix][j + iy] = M[i][j];
			}
  }
}

///
/// L1
//_______________
void AliEMCALTriggerSTU::L1(int type)
{	
  TVector2 s1, s2, s3, s4;
  fDCSConfig->GetSegmentation(s1, s2, s3, s4);
  Int_t fBkg = 0;
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
      fBkg = fBkgRho * (PatchSize()->X() * PatchSize()->Y() / 4);
      break;
    default:
      AliError("Not supported L1 trigger type");
      return;
      break;
  }
  AliDebug(999,Form("STU L1 (type %d) using subregion size (%d,%d), patch size (%d,%d)\n",type,int(SubRegionSize()->X()),int(SubRegionSize()->Y()),int(PatchSize()->X()),int(PatchSize()->Y())));
  AliDebug(999,Form("STU L1 subtracting Bkg = %d\n",fBkg));

  if (GetThreshold(type) == 0) {
    AliDebug(999,Form("STU has 0 threshold for type %d.  Skipping this trigger calculation.  Check OCDB!!!",type));
    return; // Don't execute trigger with 0 threshold.
  }

  SlidingWindow(GetThreshold(type) + fBkg);
  AliDebug(999, Form("STU type %d sliding window w/ thr %d found %d patches", type, GetThreshold(type)+ fBkgRho * int(PatchSize()->X() * PatchSize()->Y() / 4.), fPatches->GetEntriesFast()));
}

///
/// Compute threshold from V0
//___________
void AliEMCALTriggerSTU::ComputeThFromV0(int type, const Int_t M[])
{	
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

///
/// Set threshold
//___________
void AliEMCALTriggerSTU::SetThreshold(int type, Int_t v)
{	
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

///
/// Compute threshold. 
/// FIXME: need an access to the OCDB to get f(V0) parameters 
/// depending on trigger type
//___________
Int_t AliEMCALTriggerSTU::GetThreshold(int type)
{	
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

///
/// Calculate median energy of patches
//___________
Int_t AliEMCALTriggerSTU::GetMedianEnergy(){

  // iterate over non-overlaping set of (8x8 tower) subregions 
	// 4x4 FastOR units
	const int kPatchXSize = 8; // How many FastOR units in X
	const int kPatchYSize = 8; // How many FastOR units in Y

  std::vector<int> fPatchEnergies ;
	// How many columns and rows in energy patches.
	int nX = (int) fRegionSize->X()/kPatchXSize; 
	int nY = (int) fRegionSize->Y()/kPatchYSize;

  Bool_t fZeroPHOS = (fDCSConfig->GetFw() & 0xffff) == 0xd007; // Whether to ignore PHOS patches

	for (int i = 0; i < nX; i++) {
		for (int j = 0; j < nY; j++) {
      if (fZeroPHOS && (j == 2 || j == 3))  continue; // Skip middle 10 patches
			int fLocalSum = 0;
			int fXAnchor = i * kPatchXSize;
			int fYAnchor = j * kPatchYSize;
			for (int k = 0; k < kPatchXSize; k++) {
				for (int l = 0; l < kPatchYSize; l++) {
					fLocalSum += fRegion[fXAnchor+k][fYAnchor+l];
				}
			}
			fPatchEnergies.push_back(fLocalSum);
		}	
	}
	sort(fPatchEnergies.begin(),fPatchEnergies.end());
	// If EMCAL, use 24th.
	// If DCAL, use 15th.  If DCAL without PHOS, use 10th.

  // Note: this is proper only for an even number of patches
  // Technically, this is the median excluding the highest patch
	return fPatchEnergies.at(fPatchEnergies.size()/2-1);
}

///
/// Reset
//__________
void AliEMCALTriggerSTU::Reset()
{	
  fPatches->Delete();
}
