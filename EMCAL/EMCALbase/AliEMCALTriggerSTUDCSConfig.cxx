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

#include "AliEMCALTriggerSTUDCSConfig.h"
#include "TClonesArray.h"
#include "TVector2.h"

#include <iostream>

ClassImp(AliEMCALTriggerSTUDCSConfig)
ClassImp(AliEMCALTriggerSTUDCSConfig::AliEMCALTriggerSTUTRUErrorCount)
  
//_____________________________________________________________________________
AliEMCALTriggerSTUDCSConfig::AliEMCALTriggerSTUDCSConfig() : TObject(),
fGetRawData(1),
fRegion(0xFFFFFFFF),
fFw(0x2A012)
{
	//
	// AliEMCALTriggerSTUDCSConfig default constructor
	//
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 2; j++) {
			fG[i][j] = 0;
			fJ[i][j] = 0;
		}
	}
	memset(fPHOSScale, 0, sizeof(Int_t) * 4);
	memset(fTRUErrorCounts, 0, sizeof(TClonesArray *) * 32);
}

//_____________________________________________________________________________
AliEMCALTriggerSTUDCSConfig::~AliEMCALTriggerSTUDCSConfig(){
  //
  // Destructor
  //
  for(int itru = 0; itru < 32; itru++){
    if(fTRUErrorCounts[itru]) delete fTRUErrorCounts[itru];
  }
}

//_____________________________________________________________________________
void AliEMCALTriggerSTUDCSConfig::GetSegmentation(TVector2& v1, TVector2& v2, TVector2& v3, TVector2& v4) const
{
	// Get Segmentation
	
	v1.Set(1., 1.);
	v2.Set(2., 2.);
	v3.Set(4., 4.);
	
	Double_t js = 2 + (fFw >> 16);
	v4.Set(js, js);
}

//_____________________________________________________________________________
void  AliEMCALTriggerSTUDCSConfig::SetTRUErrorCounts(Int_t itru, Int_t itime, ULong64_t errorcounts){
  //
  // Set TRU error counts
  //
  if(itru >= 32) return;
  if(!fTRUErrorCounts[itru])
    fTRUErrorCounts[itru] = new TClonesArray("AliEMCALTriggerSTUDCSConfig::AliEMCALTriggerSTUTRUErrorCount");
  AliEMCALTriggerSTUTRUErrorCount test(itime, errorcounts), *found(NULL);
  if((found = dynamic_cast<AliEMCALTriggerSTUTRUErrorCount *>(fTRUErrorCounts[itru]->FindObject(&test)))){
    found->SetValue(itime, errorcounts);
  } else {
    Int_t nErrorCountsTRU = fTRUErrorCounts[itru]->GetEntries();
    new((*(fTRUErrorCounts[itru]))[nErrorCountsTRU]) AliEMCALTriggerSTUTRUErrorCount(itime, errorcounts);
  }
}

//_____________________________________________________________________________
TClonesArray *AliEMCALTriggerSTUDCSConfig::GetErrorCountsForTRU(Int_t itru) const{
  //
  // Return time-dependent error counts for a given TRU
  //
  if(itru >= 32) return NULL;
  return fTRUErrorCounts[itru];
}

//_____________________________________________________________________________
Bool_t AliEMCALTriggerSTUDCSConfig::AliEMCALTriggerSTUTRUErrorCount::IsEqual(const TObject *o) const{
  //
  // Checks for equalness according to the time stamp
  //
  const AliEMCALTriggerSTUTRUErrorCount *test = dynamic_cast<const AliEMCALTriggerSTUTRUErrorCount *>(o);
  if(!test) return false;
  return test->fTime == fTime;
}

//_____________________________________________________________________________
Int_t AliEMCALTriggerSTUDCSConfig::AliEMCALTriggerSTUTRUErrorCount::Compare(const TObject *o) const{
  //
  // Compares time-dependent error counts based on the time information
  //
  const AliEMCALTriggerSTUTRUErrorCount *test = dynamic_cast<const AliEMCALTriggerSTUTRUErrorCount *>(o);
  if(!test) return 1;
  if(fTime > test->fTime) return 1;
  if(fTime < test->fTime) return -1;
  return 0;
}

