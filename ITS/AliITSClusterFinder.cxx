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

#include "AliITSClusterFinder.h"
#include "AliRun.h"
#include "AliITS.h"

//----------------------------------------------------------

//----------------------------------------------------------

ClassImp(AliITSClusterFinder)

AliITSClusterFinder::AliITSClusterFinder
(AliITSsegmentation *seg, AliITSresponse *response, TClonesArray *digits)   
{
  // cluster finder
    fSegmentation=seg;
    fResponse=response;
    fMap = 0;
    
    fDigits=digits;
    fNdigits = fDigits->GetEntriesFast();

    fNRawClusters=0;

    SetNperMax();
    SetClusterSize();
    SetDeclusterFlag();

    fNPeaks=-1;
}

//----------------------------------------------------------
AliITSClusterFinder::AliITSClusterFinder()
{
  // constructor
    fResponse=0;
    
    fDigits=0;
    fNdigits = 0;

    fNRawClusters=0;
    fMap = 0;


    SetNperMax();
    SetClusterSize();
    SetDeclusterFlag();

    fNPeaks=-1;
}
 
//__________________________________________________________________________
AliITSClusterFinder::AliITSClusterFinder(const AliITSClusterFinder &source){
  //     Copy Constructor 
  if(&source == this) return;
  this->fDigits = source.fDigits;
  this->fNdigits = source.fNdigits;
  this->fResponse = source.fResponse;
  this->fSegmentation = source.fSegmentation;
  this->fNRawClusters = source.fNRawClusters;
  this->fMap = source.fMap;
  this->fNperMax = source.fNperMax;
  this->fDeclusterFlag = source.fDeclusterFlag;
  this->fClusterSize = source.fClusterSize;
  this->fNPeaks = source.fNPeaks;
  return;
}

//_________________________________________________________________________
AliITSClusterFinder& 
  AliITSClusterFinder::operator=(const AliITSClusterFinder &source) {
  //    Assignment operator
  if(&source == this) return *this;
  this->fDigits = source.fDigits;
  this->fNdigits = source.fNdigits;
  this->fResponse = source.fResponse;
  this->fSegmentation = source.fSegmentation;
  this->fNRawClusters = source.fNRawClusters;
  this->fMap = source.fMap;
  this->fNperMax = source.fNperMax;
  this->fDeclusterFlag = source.fDeclusterFlag;
  this->fClusterSize = source.fClusterSize;
  this->fNPeaks = source.fNPeaks;
  return *this;
}

//----------------------------------------------------------
void AliITSClusterFinder::AddCluster(Int_t branch, AliITSRawCluster *c)
{
  //
  // Add a raw cluster copy to the list
  //
    AliITS *iTS=(AliITS*)gAlice->GetModule("ITS");
    iTS->AddCluster(branch,c); 
    fNRawClusters++;

}


//----------------------------------------------------------
void AliITSClusterFinder::AddCluster(Int_t branch, AliITSRawCluster *c, AliITSRecPoint &rp)
{
  //
  // Add a raw cluster copy to the list
  //
    AliITS *iTS=(AliITS*)gAlice->GetModule("ITS");
    iTS->AddCluster(branch,c); 
    fNRawClusters++;
    iTS->AddRecPoint(rp); 
}





