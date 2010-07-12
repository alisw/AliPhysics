// @(#) $Id$
// Original: AliHLTClustFinderNew.cxx,v 1.29 2005/06/14 10:55:21 cvetan Exp 

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Anders Vestbo, Constantin Loizides                    *
//* Developers:      Kenneth Aamodt <kenneth.aamodt@student.uib.no>        *
//*                  Kalliopi Kanaki                                       *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

//  @file   AliHLTTPCClusterFinder.cxx
//  @author Kenneth Aamodt, Kalliopi Kanaki
//  @date   
//  @brief  Cluster Finder for the TPC
//  @note 

#include "AliHLTTPCDigitReader.h"
#include "AliHLTTPCRootTypes.h"
#include "AliHLTTPCLogging.h"
#include "AliHLTTPCClusterFinder.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCMemHandler.h"
#include "AliHLTTPCPad.h"
#include <sys/time.h>
#include <algorithm>
#include <cmath>
#include "AliTPCcalibDB.h"
#include "AliTPCTransform.h"

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTPCClusterFinder)

AliHLTTPCClusterFinder::AliHLTTPCClusterFinder()
  :
  fClustersHWAddressVector(),
  fRowPadVector(),
  fSpacePointData(NULL),
  fDigitReader(NULL),
  fPtr(NULL),
  fSize(0),
  fDeconvTime(kFALSE),
  fDeconvPad(kFALSE),
  fStdout(kFALSE),
  fCalcerr(kTRUE),
  fRawSP(kFALSE),
  fFirstRow(0),
  fLastRow(0),
  fCurrentRow(0),
  fCurrentSlice(0),
  fCurrentPatch(0),
  fMatch(1),
  fThreshold(10),
  fNClusters(0),
  fMaxNClusters(0),
  fXYErr(0.2),
  fZErr(0.3),
  fOccupancyLimit(1.0),
  fUnsorted(0),
  fVectorInitialized(kFALSE),
  fClusters(),
  fClustersMCInfo(),
  fMCDigits(),
  fNumberOfPadsInRow(NULL),
  fNumberOfRows(0),
  fRowOfFirstCandidate(0),
  fDoPadSelection(kFALSE),
  fFirstTimeBin(0),
  fLastTimeBin(AliHLTTPCTransform::GetNTimeBins()),
  fTotalChargeOfPreviousClusterCandidate(0),
  fChargeOfCandidatesFalling(kFALSE),
  f32BitFormat(kFALSE),
  fDoMC(kFALSE),
  fClusterMCVector(),
  fOfflineTransform(NULL),
  fOfflineTPCRecoParam(),
  fTimeMeanDiff(2),
  fReleaseMemory(0)
{
  //constructor  
  fOfflineTransform = AliTPCcalibDB::Instance()->GetTransform(); 
  if(!fOfflineTransform){
    HLTError("AliHLTTPCClusterFinder():  Offline transform not in AliTPCcalibDB.");
  }
  else{
    fOfflineTPCRecoParam.SetUseExBCorrection(1);
    fOfflineTPCRecoParam.SetUseTOFCorrection(1);
    fOfflineTransform->SetCurrentRecoParam(&fOfflineTPCRecoParam);
  }
}

AliHLTTPCClusterFinder::~AliHLTTPCClusterFinder(){
  // see header file for class documentation
  
  //destructor
  if(fVectorInitialized){
    DeInitializePadArray();
  }
  if(fNumberOfPadsInRow){
    delete [] fNumberOfPadsInRow;
    fNumberOfPadsInRow=NULL;
  }
}

void AliHLTTPCClusterFinder::InitSlice(Int_t slice,Int_t patch,Int_t nmaxpoints){
  // see header file for class documentation

  //init slice
  fNClusters = 0;
  fMaxNClusters = nmaxpoints;
  fCurrentSlice = slice;
  fCurrentPatch = patch;
  fFirstRow=AliHLTTPCTransform::GetFirstRow(patch);
  fLastRow=AliHLTTPCTransform::GetLastRow(patch);

  fClusters.clear();
  fClustersMCInfo.clear();
  fMCDigits.clear();
  fClusterMCVector.clear();   
}

void AliHLTTPCClusterFinder::InitializePadArray(){
  // see header file for class documentation
  
  if(fCurrentPatch>5||fCurrentPatch<0){
    HLTFatal("Patch is not set");
    return;
  }

  HLTDebug("Patch number=%d",fCurrentPatch);

  fFirstRow = AliHLTTPCTransform::GetFirstRow(fCurrentPatch);
  fLastRow = AliHLTTPCTransform::GetLastRow(fCurrentPatch);

  fNumberOfRows=fLastRow-fFirstRow+1;
  fNumberOfPadsInRow= new UInt_t[fNumberOfRows];

  memset( fNumberOfPadsInRow, 0, sizeof(Int_t)*(fNumberOfRows));

  fRowPadVector.clear();

  for(UInt_t i=0;i<fNumberOfRows;i++){
    fNumberOfPadsInRow[i]=AliHLTTPCTransform::GetNPads(i+fFirstRow);
    AliHLTTPCPadVector tmpRow;
    for(UInt_t j=0;j<=fNumberOfPadsInRow[i];j++){
      AliHLTTPCPad *tmpPad = new AliHLTTPCPad(2);
      tmpPad->SetID(i,j);
      tmpRow.push_back(tmpPad);
    }
    fRowPadVector.push_back(tmpRow);
  }
  fVectorInitialized=kTRUE;
}

Int_t AliHLTTPCClusterFinder::DeInitializePadArray(){
  // see header file for class documentation

  if( fVectorInitialized ){
    for(UInt_t i=0;i<fNumberOfRows;i++){
      for(UInt_t j=0;j<=fNumberOfPadsInRow[i];j++){
	delete fRowPadVector[i][j];
	fRowPadVector[i][j]=NULL;
      }
      fRowPadVector[i].clear();
    }
    fRowPadVector.clear();
    delete[] fNumberOfPadsInRow;
    fNumberOfPadsInRow = 0;
  }
  fVectorInitialized=kFALSE;
  return 1;
}


void AliHLTTPCClusterFinder::SetOutputArray(AliHLTTPCSpacePointData *pt){
  // see header file for class documentation
  //set pointer to output
  fSpacePointData = pt;
}


void AliHLTTPCClusterFinder::ReadDataUnsorted(void* ptr,unsigned long size){
  // see header file for class documentation
  //set input pointer
  fPtr = (UChar_t*)ptr;
  fSize = size;

  if(!fVectorInitialized){
    InitializePadArray();
  }

  if (fDigitReader->InitBlock(fPtr,fSize,fFirstRow,fLastRow,fCurrentPatch,fCurrentSlice)<0) {
    HLTError("failed setting up digit reader (InitBlock)");
    return;
  }
  
  while(fDigitReader->NextChannel()){
    UInt_t row=fDigitReader->GetRow();
    UInt_t pad=fDigitReader->GetPad();

    if(row>=fRowPadVector.size()){
      HLTError("Row number is to large: %d, max is %d",row,fRowPadVector.size()-1);
      continue;
    }
    if(pad>=fRowPadVector[row].size()){
      HLTError("Pad number is to large: %d, max is %d",pad,fRowPadVector[row].size());
      continue;
    }

    while(fDigitReader->NextBunch()){
      if(fDigitReader->GetBunchSize()>1){//to remove single timebin values, this will have to change at some point
	UInt_t time = fDigitReader->GetTime();
	if((Int_t)time>=fFirstTimeBin && (Int_t)time+fDigitReader->GetBunchSize()<=fLastTimeBin){
	  // Kenneth: 20-04-09. The following if have been added because of inconsistency in the 40 bit decoder and the 32 bit decoder.
	  // GetSignals() in the 40 bit decoder returns an array of UInt_t while the 32 bit one returns UShort_t
	  // The same is true for the function ReadDataUnsortedDeconvoluteTime() below.
	  // In addition the signals are organized in the opposite direction
	  if(f32BitFormat){
	    const UShort_t *bunchData= fDigitReader->GetSignalsShort();
	    AliHLTTPCClusters candidate;
	    for(Int_t i=fDigitReader->GetBunchSize()-1;i>=0;i--){
	      candidate.fTotalCharge+=bunchData[i];	
	      candidate.fTime += time*bunchData[i];
	      candidate.fTime2 += time*time*bunchData[i];
	      if(bunchData[i]>candidate.fQMax){
		candidate.fQMax=bunchData[i];
	      }
	      time++;
	    }
	    if(candidate.fTotalCharge>0){
	      candidate.fMean=candidate.fTime/candidate.fTotalCharge;
	      candidate.fPad=candidate.fTotalCharge*pad;
	      candidate.fPad2=candidate.fPad*pad;
	      candidate.fLastMergedPad=pad;
	      candidate.fRowNumber=row+fDigitReader->GetRowOffset();
	    }
	    if(fRowPadVector[row][pad] != NULL){
	      fRowPadVector[row][pad]->AddClusterCandidate(candidate);
	    }
	  }
	  else{
	    const UInt_t *bunchData= fDigitReader->GetSignals();
	    AliHLTTPCClusters candidate;
	    const AliHLTTPCDigitData* digits = NULL;
	    if(fDoMC && (digits = fDigitReader->GetBunchDigits())!=NULL){
	      for(Int_t i=0;i<fDigitReader->GetBunchSize();i++){
		candidate.fTotalCharge+=bunchData[i];	
		candidate.fTime += time*bunchData[i];
		candidate.fTime2 += time*time*bunchData[i];
		if(bunchData[i]>candidate.fQMax){
		  candidate.fQMax=bunchData[i];
		}
		fMCDigits.push_back(digits[i]);
		time++;
	      }
	    }
	    else{
	      for(Int_t i=0;i<fDigitReader->GetBunchSize();i++){
		candidate.fTotalCharge+=bunchData[i];	
		candidate.fTime += time*bunchData[i];
		candidate.fTime2 += time*time*bunchData[i];
		if(bunchData[i]>candidate.fQMax){
		  candidate.fQMax=bunchData[i];
		}
		time++;
	      }
	    }
	    if(candidate.fTotalCharge>0){
	      candidate.fMean=candidate.fTime/candidate.fTotalCharge;
	      candidate.fPad=candidate.fTotalCharge*pad;
	      candidate.fPad2=candidate.fPad*pad;
	      candidate.fLastMergedPad=pad;
	      candidate.fRowNumber=row+fDigitReader->GetRowOffset();
	    }
	    if(fRowPadVector[row][pad] != NULL){
	      fRowPadVector[row][pad]->AddClusterCandidate(candidate);
	      if(fDoMC){
		fRowPadVector[row][pad]->AddCandidateDigits(fMCDigits);
		fMCDigits.clear();
	      }
	    }
	  }
	}
      }
    }
  }
}

void AliHLTTPCClusterFinder::ReadDataUnsortedDeconvoluteTime(void* ptr,unsigned long size){
  // see header file for class documentation

  //set input pointer
  fPtr = (UChar_t*)ptr;
  fSize = size;

  if(!fVectorInitialized){
    InitializePadArray();
  }

  if (fDigitReader->InitBlock(fPtr,fSize,fFirstRow,fLastRow,fCurrentPatch,fCurrentSlice)<0) {
    HLTError("failed setting up digit reader (InitBlock)");
    return;
  }
  
  while(fDigitReader->NextChannel()){
    UInt_t row=fDigitReader->GetRow();
    UInt_t pad=fDigitReader->GetPad();

    while(fDigitReader->NextBunch()){
      if(fDigitReader->GetBunchSize()>1){//to remove single timebin values, this will have to change at some point
	UInt_t time = fDigitReader->GetTime();
	if((Int_t)time>=fFirstTimeBin && (Int_t)time+fDigitReader->GetBunchSize()<=fLastTimeBin){
	  Int_t indexInBunchData=0;
	  Bool_t moreDataInBunch=kFALSE;
	  UInt_t prevSignal=0;
	  Bool_t signalFalling=kFALSE;

	  // Kenneth: 20-04-09. The following if have been added because of inconsistency in the 40 bit decoder and the 32 bit decoder.
	  // GetSignals() in the 40 bit decoder returns an array of UInt_t while the 32 bit one returns UShort_t
	  // The same is true for the function ReadDataUnsorted() above.
	  // In addition the signals are organized in the opposite direction
	  if(f32BitFormat){
	    indexInBunchData = fDigitReader->GetBunchSize();
	    const UShort_t *bunchData= fDigitReader->GetSignalsShort();
	    
	    do{
	      AliHLTTPCClusters candidate;
	      //for(Int_t i=indexInBunchData;i<fDigitReader->GetBunchSize();i++){
	      for(Int_t i=indexInBunchData;i>=0;i--){
		// Checks if one need to deconvolute the signals
		if(bunchData[i]>prevSignal && signalFalling==kTRUE){
		  if(i<fDigitReader->GetBunchSize()-1){ // means there are more than one signal left in the bunch
		    moreDataInBunch=kTRUE;
		    prevSignal=0;
		  }
		  break;
		}
		
		// Checks if the signal is 0, then quit processing the data.
		if(bunchData[i]==0 && i<fDigitReader->GetBunchSize()-1){//means we have 0 data fom the rcu, might happen depending on the configuration settings
		  moreDataInBunch=kTRUE;
		  prevSignal=0;
		  break;
		}
		
		if(prevSignal>bunchData[i]){//means the peak of the signal has been reached and deconvolution will happen if the signal rise again.
		  signalFalling=kTRUE;
		}
		candidate.fTotalCharge+=bunchData[i];	
		candidate.fTime += time*bunchData[i];
		candidate.fTime2 += time*time*bunchData[i];
		if(bunchData[i]>candidate.fQMax){
		  candidate.fQMax=bunchData[i];
		}
		
		prevSignal=bunchData[i];
		time++;
		indexInBunchData--;
	      }
	      if(candidate.fTotalCharge>0){
		candidate.fMean=candidate.fTime/candidate.fTotalCharge;
		candidate.fPad=candidate.fTotalCharge*pad;
		candidate.fPad2=candidate.fPad*pad;
		candidate.fLastMergedPad=pad;
		candidate.fRowNumber=row+fDigitReader->GetRowOffset();
	      }
	      fRowPadVector[row][pad]->AddClusterCandidate(candidate);
	      fRowPadVector[row][pad]->AddCandidateDigits(fMCDigits);
	      if(indexInBunchData<fDigitReader->GetBunchSize()-1){
		moreDataInBunch=kFALSE;
	      }
	    }while(moreDataInBunch);
	  }
	  else{
	    const UInt_t *bunchData= fDigitReader->GetSignals();
	    do{
	      AliHLTTPCClusters candidate;
	      for(Int_t i=indexInBunchData;i<fDigitReader->GetBunchSize();i++){
		// Checks if one need to deconvolute the signals
		if(bunchData[i]>prevSignal && signalFalling==kTRUE){
		  if(i<fDigitReader->GetBunchSize()-1){ // means there are more than one signal left in the bunch
		    moreDataInBunch=kTRUE;
		    prevSignal=0;
		  }
		  break;
		}
		
		// Checks if the signal is 0, then quit processing the data.
		if(bunchData[i]==0 && i<fDigitReader->GetBunchSize()-1){//means we have 0 data fom the rcu, might happen depending on the configuration settings
		  moreDataInBunch=kTRUE;
		  prevSignal=0;
		  break;
		}
		
		if(prevSignal>bunchData[i]){//means the peak of the signal has been reached and deconvolution will happen if the signal rise again.
		  signalFalling=kTRUE;
		}

		candidate.fTotalCharge+=bunchData[i];	
		candidate.fTime += time*bunchData[i];
		candidate.fTime2 += time*time*bunchData[i];
		if(bunchData[i]>candidate.fQMax){
		  candidate.fQMax=bunchData[i];
		}
		
		prevSignal=bunchData[i];
		time++;
		indexInBunchData++;
	      }
	      if(candidate.fTotalCharge>0){
		candidate.fMean=candidate.fTime/candidate.fTotalCharge;
		candidate.fPad=candidate.fTotalCharge*pad;
		candidate.fPad2=candidate.fPad*pad;
		candidate.fLastMergedPad=pad;
		candidate.fRowNumber=row+fDigitReader->GetRowOffset();
	      }
	      fRowPadVector[row][pad]->AddClusterCandidate(candidate);
	      if(indexInBunchData<fDigitReader->GetBunchSize()-1){
		moreDataInBunch=kFALSE;
	      }
	    }while(moreDataInBunch);
	  }
	}
      }
    }
  }
}

Bool_t AliHLTTPCClusterFinder::ComparePads(AliHLTTPCPad *nextPad,AliHLTTPCClusters* cluster,Int_t nextPadToRead){
  // see header file for class documentation

  //Checking if we have a match on the next pad
  for(UInt_t candidateNumber=0;candidateNumber<nextPad->fClusterCandidates.size();candidateNumber++){
    if(nextPad->fUsedClusterCandidates[candidateNumber] == 1){
      continue;
    }
    AliHLTTPCClusters *candidate =&nextPad->fClusterCandidates[candidateNumber]; 
    //    if(cluster->fMean-candidate->fMean==1 || candidate->fMean-cluster->fMean==1 || cluster->fMean-candidate->fMean==0){
    
    if( abs((Int_t)(cluster->fMean - candidate->fMean)) <= fTimeMeanDiff ){
      if(fDeconvPad){
	if(candidate->fTotalCharge<fTotalChargeOfPreviousClusterCandidate){//peak is reached
	  fChargeOfCandidatesFalling=kTRUE;
	}
	if(candidate->fTotalCharge>fTotalChargeOfPreviousClusterCandidate && fChargeOfCandidatesFalling==kTRUE){//we have deconvolution
	  return kFALSE;
	}
      }
      cluster->fMean=candidate->fMean;
      cluster->fTotalCharge+=candidate->fTotalCharge;
      cluster->fTime += candidate->fTime;
      cluster->fTime2 += candidate->fTime2;
      cluster->fPad+=candidate->fPad;
      cluster->fPad2+=candidate->fPad2;
      cluster->fLastMergedPad=candidate->fPad;
      if(candidate->fQMax>cluster->fQMax){
	cluster->fQMax=candidate->fQMax;
      }
      if(fDoMC){
	FillMCClusterVector(nextPad->GetCandidateDigits(candidateNumber));
      }

      if(fDoPadSelection){
	UInt_t rowNo = nextPad->GetRowNumber();
	UInt_t padNo = nextPad->GetPadNumber();
	if(padNo-1>0){
	  fRowPadVector[rowNo][padNo-2]->fSelectedPad=kTRUE;
	  fRowPadVector[rowNo][padNo-2]->fHWAddress=(AliHLTUInt16_t)fDigitReader->GetAltroBlockHWaddr(rowNo,padNo-2);
	}
	fRowPadVector[rowNo][padNo-1]->fSelectedPad=kTRUE;// quick solution to set the first pad to selected
	fRowPadVector[rowNo][padNo-1]->fHWAddress=(AliHLTUInt16_t)fDigitReader->GetAltroBlockHWaddr(rowNo,padNo-1);
	fRowPadVector[rowNo][padNo]->fSelectedPad=kTRUE;
	fRowPadVector[rowNo][padNo]->fHWAddress=(AliHLTUInt16_t)fDigitReader->GetAltroBlockHWaddr(rowNo,padNo);
      }

      //setting the matched pad to used
      nextPad->fUsedClusterCandidates[candidateNumber]=1;
      nextPadToRead++;
      if(nextPadToRead<(Int_t)fNumberOfPadsInRow[fRowOfFirstCandidate]){
	nextPad=fRowPadVector[fRowOfFirstCandidate][nextPadToRead];
	ComparePads(nextPad,cluster,nextPadToRead);
      }
      else{
	return kFALSE;
      }
    }
  }
  return kFALSE;
}

Int_t AliHLTTPCClusterFinder::FillHWAddressList(AliHLTUInt16_t *hwaddlist, Int_t maxHWadd){
  // see header file for class documentation

  Int_t counter=0;
  for(UInt_t row=0;row<fNumberOfRows;row++){
    for(UInt_t pad=0;pad<fNumberOfPadsInRow[row]-1;pad++){
      if(fRowPadVector[row][pad]->fSelectedPad){
       if(counter<maxHWadd){
	 hwaddlist[counter]=(AliHLTUInt16_t)fRowPadVector[row][pad]->fHWAddress;
	 counter++;
       }
       else{
	 HLTWarning("To many hardwareaddresses, skip adding");
       }
       
      }
    }
  }  
  return counter;
}
 

Int_t AliHLTTPCClusterFinder::FillOutputMCInfo(AliHLTTPCClusterFinder::ClusterMCInfo * outputMCInfo, Int_t maxNumberOfClusterMCInfo){
  // see header file for class documentation
  
  Int_t counter=0;
 
  for(UInt_t mc=0;mc<fClustersMCInfo.size();mc++){
    if(counter<maxNumberOfClusterMCInfo){
      outputMCInfo[counter] = fClustersMCInfo[mc];
      counter++;
    }
    else{
      HLTWarning("To much MCInfo has been added (no more space), skip adding");
    }
  }
  return counter;
}

void AliHLTTPCClusterFinder::FindClusters(){
  // see header file for function documentation

  AliHLTTPCClusters* tmpCandidate=NULL;
  for(UInt_t row=0;row<fNumberOfRows;row++){
    fRowOfFirstCandidate=row;
    for(UInt_t pad=0;pad<fNumberOfPadsInRow[row];pad++){
      AliHLTTPCPad *tmpPad=fRowPadVector[row][pad];
      for(size_t candidate=0;candidate<tmpPad->fClusterCandidates.size();candidate++){
	if(tmpPad->fUsedClusterCandidates[candidate]){
	  continue;
	}
	tmpCandidate=&tmpPad->fClusterCandidates[candidate];
	UInt_t tmpTotalCharge=tmpCandidate->fTotalCharge;

	if(fDoMC){
	  fClusterMCVector.clear();
	  FillMCClusterVector(tmpPad->GetCandidateDigits(candidate));
	}

	ComparePads(fRowPadVector[row][pad+1],tmpCandidate,pad+1);
	if(tmpCandidate->fTotalCharge>tmpTotalCharge){
	  //we have a cluster
	  fClusters.push_back(*tmpCandidate);
	  if(fDoMC){
	    //sort the vector (large->small) according to weight and remove elements above 2 (keep 0 1 and 2) 
	    sort(fClusterMCVector.begin(),fClusterMCVector.end(), MCWeight::CompareWeights );
	    ClusterMCInfo tmpClusterMCInfo;

	    MCWeight zeroMC;
	    zeroMC.fMCID=-1;
	    zeroMC.fWeight=0;

	    if(fClusterMCVector.size()>0){
	      tmpClusterMCInfo.fClusterID[0]=fClusterMCVector.at(0);
	    }
	    else{
	      tmpClusterMCInfo.fClusterID[0]=zeroMC;
	    }

	    if(fClusterMCVector.size()>1){
	    tmpClusterMCInfo.fClusterID[1]=fClusterMCVector.at(1);
	    }
	    else{
	      tmpClusterMCInfo.fClusterID[1]=zeroMC;
	    }

	    if(fClusterMCVector.size()>2){
	    tmpClusterMCInfo.fClusterID[2]=fClusterMCVector.at(2);
	    }
	    else{
	      tmpClusterMCInfo.fClusterID[2]=zeroMC;
	    }

	    fClustersMCInfo.push_back(tmpClusterMCInfo);
	  }
	  
	}
      }
      tmpPad->ClearCandidates();
    }
    fRowPadVector[row][fNumberOfPadsInRow[row]]->ClearCandidates();
  }

  HLTInfo("Found %d clusters.",fClusters.size());

  //TODO:  Change so it stores AliHLTTPCSpacePointData directly, instead of this copying
  
  AliClusterData * clusterlist = new AliClusterData[fClusters.size()]; //Clusterlist
  for(unsigned int i=0;i<fClusters.size();i++){
    clusterlist[i].fTotalCharge = fClusters[i].fTotalCharge;
    clusterlist[i].fPad = fClusters[i].fPad;
    clusterlist[i].fPad2 = fClusters[i].fPad2;
    clusterlist[i].fTime = fClusters[i].fTime;
    clusterlist[i].fTime2 = fClusters[i].fTime2;
    clusterlist[i].fMean = fClusters[i].fMean;
    clusterlist[i].fFlags = fClusters[i].fFlags;
    clusterlist[i].fChargeFalling = fClusters[i].fChargeFalling;
    clusterlist[i].fLastCharge = fClusters[i].fLastCharge;
    clusterlist[i].fLastMergedPad = fClusters[i].fLastMergedPad;
    clusterlist[i].fRow = fClusters[i].fRowNumber;
    clusterlist[i].fQMax = fClusters[i].fQMax;
  }

  WriteClusters(fClusters.size(),clusterlist);
  delete [] clusterlist;
  fClusters.clear();
  if( fReleaseMemory ) DeInitializePadArray();// call this when  the -releaseMemory flag is set
}


Bool_t AliHLTTPCClusterFinder::UpdateCalibDB(){
  
  //update the db
  AliTPCcalibDB::Instance()->Update();

  //uptate the transform class
  AliTPCTransform * tmp = AliTPCcalibDB::Instance()->GetTransform(); 
  if(!tmp){
    HLTError("AliHLTTPCClusterFinder::UpdateCAlibDB: Offline transform not in AliTPCcalibDB.");
    return 0;
  }
  fOfflineTransform = tmp;
  return 1;
}

//---------------------------------- Under this line the old sorted clusterfinder functions can be found --------------------------------


void AliHLTTPCClusterFinder::PrintClusters(){
  // see header file for class documentation

  for(size_t i=0;i<fClusters.size();i++){
    HLTInfo("Cluster number: %d",i);
    HLTInfo("Row: %d \t Pad: %d",fClusters[i].fRowNumber,fClusters[i].fPad/fClusters[i].fTotalCharge);
    HLTInfo("Total Charge:   %d",fClusters[i].fTotalCharge);
    HLTInfo("fPad:           %d",fClusters[i].fPad);
    HLTInfo("PadError:       %d",fClusters[i].fPad2);
    HLTInfo("TimeMean:       %d",fClusters[i].fTime/fClusters[i].fTotalCharge);
    HLTInfo("TimeError:      %d",fClusters[i].fTime2);
    HLTInfo("EndOfCluster:");
  }
}

void AliHLTTPCClusterFinder::FillMCClusterVector(vector<AliHLTTPCDigitData> digitData){
  // see header file for class documentation

  for(UInt_t d=0;d<digitData.size();d++){
    Int_t nIDsInDigit = (digitData.at(d).fTrackID[0]>=0) + (digitData.at(d).fTrackID[1]>=0) + (digitData.at(d).fTrackID[2]>=0);
    for(Int_t id=0; id<3; id++){
      if(digitData.at(d).fTrackID[id]>=0){
	Bool_t matchFound = kFALSE;
	MCWeight mc;
	mc.fMCID = digitData.at(d).fTrackID[id];
	mc.fWeight = ((Float_t)digitData.at(d).fCharge)/nIDsInDigit;
	for(UInt_t i=0;i<fClusterMCVector.size();i++){
	  if(mc.fMCID == fClusterMCVector.at(i).fMCID){
	    fClusterMCVector.at(i).fWeight += mc.fWeight;
	    matchFound = kTRUE;
	  }
	}
	if(matchFound == kFALSE){
	  fClusterMCVector.push_back(mc);
	}
      }
    }
  }
}


void AliHLTTPCClusterFinder::Read(void* ptr,unsigned long size){
  //set input pointer
  fPtr = (UChar_t*)ptr;
  fSize = size;
}

void AliHLTTPCClusterFinder::ProcessDigits(){
  // see header file for class documentation

  int iResult=0;
  bool readValue = true;
  Int_t newRow = 0;    
  Int_t rowOffset = 0;
  UShort_t time=0,newTime=0;
  UInt_t pad=0,newPad=0;
  AliHLTTPCSignal_t charge=0;

  fNClusters = 0;

  // initialize block for reading packed data
  iResult=fDigitReader->InitBlock(fPtr,fSize,fFirstRow,fLastRow,fCurrentPatch,fCurrentSlice);
  if (iResult<0) return;

  readValue = fDigitReader->Next();

  // Matthias 08.11.2006 the following return would cause termination without writing the
  // ClusterData and thus would block the component. I just want to have the commented line
  // here for information
  //if (!readValue)return;

  pad = fDigitReader->GetPad();
  time = fDigitReader->GetTime();
  fCurrentRow = fDigitReader->GetRow();

  if ( fCurrentPatch >= 2 ) // Outer sector, patches 2, 3, 4, 5
    rowOffset = AliHLTTPCTransform::GetFirstRow( 2 );

  fCurrentRow += rowOffset;

  UInt_t lastpad = 123456789;
  const UInt_t kPadArraySize=5000;
  const UInt_t kClusterListSize=10000;
  AliClusterData *pad1[kPadArraySize]; //2 lists for internal memory=2pads
  AliClusterData *pad2[kPadArraySize]; //2 lists for internal memory=2pads
  AliClusterData clusterlist[kClusterListSize]; //Clusterlist

  AliClusterData **currentPt;  //List of pointers to the current pad
  AliClusterData **previousPt; //List of pointers to the previous pad
  currentPt = pad2;
  previousPt = pad1;
  UInt_t nprevious=0,ncurrent=0,ntotal=0;

  /* quick implementation of baseline calculation and zero suppression
     open a pad object for each pad and delete it after processing.
     later a list of pad objects with base line history can be used
     The whole thing only works if we really get unprocessed raw data, if
     the data is already zero suppressed, there might be gaps in the time
     bins.
   */
  Int_t gatingGridOffset=50;
  if(fFirstTimeBin>0){
    gatingGridOffset=fFirstTimeBin;
  }
  AliHLTTPCPad baseline(gatingGridOffset, AliHLTTPCTransform::GetNTimeBins());
  // just to make later conversion to a list of objects easier
  AliHLTTPCPad* pCurrentPad=NULL;
  /*
    if (fSignalThreshold>=0) {
    pCurrentPad=&baseline;
    baseline.SetThreshold(fSignalThreshold);
  }
  */
  while ( readValue!=0 && iResult>=0){   // Reads through all digits in block
    iResult=0;

    if(pad != lastpad){
      //This is a new pad
      
      //Switch the lists:
      if(currentPt == pad2){
	currentPt = pad1;
	previousPt = pad2;
      }
      else {
	currentPt = pad2;
	previousPt = pad1;
      }
      nprevious = ncurrent;
      ncurrent = 0;
      if(pad != lastpad+1){
	//this happens if there is a pad with no signal.
	nprevious = ncurrent = 0;
      }
      lastpad = pad;
    }

    Bool_t newcluster = kTRUE;
    UInt_t seqcharge=0,seqaverage=0,seqerror=0;
    AliHLTTPCSignal_t lastcharge=0;
    UInt_t bLastWasFalling=0;
    Int_t newbin=-1;


    if(fDeconvTime){
      redo: //This is a goto.
      
      if(newbin > -1){
	//bin = newbin;
	newbin = -1;
      }
	  
      lastcharge=0;
      bLastWasFalling = 0;
    }

    while(iResult>=0){ //Loop over time bins of current pad
      // read all the values for one pad at once to calculate the base line
      if (pCurrentPad) {
	if (!pCurrentPad->IsStarted()) {
	  //HLTDebug("reading data for pad %d, padrow %d", fDigitReader->GetPad(), fDigitReader->GetRow()+rowOffset);
	  pCurrentPad->SetID(fDigitReader->GetRow()+rowOffset,fDigitReader->GetPad());
	  if ((pCurrentPad->StartEvent())>=0) {
	    do {
	      if ((fDigitReader->GetRow()+rowOffset)!=pCurrentPad->GetRowNumber()) break;
	      if (fDigitReader->GetPad()!=pCurrentPad->GetPadNumber()) break;
	      pCurrentPad->SetRawData(fDigitReader->GetTime(), fDigitReader->GetSignal());
	      //HLTDebug("set raw data to pad: bin %d charge %d", fDigitReader->GetTime(), fDigitReader->GetSignal());
	    } while ((readValue = fDigitReader->Next())!=0);
	  }
	  pCurrentPad->CalculateBaseLine(AliHLTTPCTransform::GetNTimeBins()/2);
	  if (pCurrentPad->Next(kTRUE/*do zero suppression*/)==0) {
	    //HLTDebug("no data available after zero suppression");
	    pCurrentPad->StopEvent();
	    pCurrentPad->ResetHistory();
	    break;
	  }
	  time=pCurrentPad->GetCurrentPosition();
	  if (time>pCurrentPad->GetSize()) {
	    HLTError("invalid time bin for pad");
	    break;
	  }
	}
      }
      if (pCurrentPad) {
	Float_t occupancy=pCurrentPad->GetOccupancy();
	//HLTDebug("pad %d occupancy level: %f", pCurrentPad->GetPadNumber(), occupancy);
	if ( occupancy < fOccupancyLimit ) {
	  charge = pCurrentPad->GetCorrectedData();
	} else {
	  charge = 0;
	  //HLTDebug("ignoring pad %d with occupancy level %f", pCurrentPad->GetPadNumber(), occupancy);
	}
      } else {
	charge = fDigitReader->GetSignal();
      }
      //HLTDebug("get next charge value: position %d charge %d", time, charge);


      // CHARGE DEBUG
      if (fDigitReader->GetRow() == 90){
/////	  LOG(AliHLTTPCLog::kFatal,"AliHLTTPCClusterFinder::Row","row90")  << "PAD=" <<  fDigitReader->GetPad() << "  TIME=" <<  fDigitReader->GetTime() 
	  //					   << "  SIGNAL=" <<  fDigitReader->GetSignal() << ENDLOG;

      }

      if(time >= AliHLTTPCTransform::GetNTimeBins()){
	HLTWarning("Pad %d: Timebin (%d) out of range (%d)", pad, time, AliHLTTPCTransform::GetNTimeBins());
	iResult=-ERANGE;
      }


      //Get the current ADC-value
      if(fDeconvTime){

	//Check if the last pixel in the sequence is smaller than this
	if(charge > lastcharge){
	  if(bLastWasFalling){
	    newbin = 1;
	    break;
	  }
	}
	else bLastWasFalling = 1; //last pixel was larger than this
	lastcharge = charge;
      }
	  
      //Sum the total charge of this sequence
      seqcharge += charge;
      seqaverage += time*charge;
      seqerror += time*time*charge;
      
      if (pCurrentPad) {
	
	if((pCurrentPad->Next(kTRUE/*do zero suppression*/))==0) {
	  pCurrentPad->StopEvent();
	  pCurrentPad->ResetHistory();
	  if(readValue) {
	    newPad = fDigitReader->GetPad();
	    newTime = fDigitReader->GetTime();
	    newRow = fDigitReader->GetRow() + rowOffset;
	  }
	  break;
	}

	newPad=pCurrentPad->GetPadNumber();
	newTime=pCurrentPad->GetCurrentPosition();
	newRow=pCurrentPad->GetRowNumber();
      } else {
      readValue = fDigitReader->Next();
      //Check where to stop:
      if(!readValue) break; //No more value

      newPad = fDigitReader->GetPad();
      newTime = fDigitReader->GetTime();
      newRow = fDigitReader->GetRow() + rowOffset;
      }

      if(newPad != pad)break; //new pad
      if(newTime != time+1) break; //end of sequence
      if(iResult<0) break;

      // pad = newpad;    is equal
      time = newTime;

    }//end loop over sequence

    //HLTDebug("ended time bin sequence loop: seqcharge=%d readValue=%d", seqcharge, readValue);
    //HLTDebug("pad=%d newpad=%d current row=%d newrow=%d", pad, newPad, fCurrentRow, newRow);
    if (seqcharge<=0) {
      // with active zero suppression zero values are possible
      continue;
    }

    //Calculate mean of sequence:
    Int_t seqmean=0;
    if(seqcharge)
      seqmean = seqaverage/seqcharge;
    else{
      LOG(AliHLTTPCLog::kFatal,"AliHLTTPCClusterFinder::ProcessRow","Data")
	<<"Error in data given to the cluster finder"<<ENDLOG;
      seqmean = 1;
      seqcharge = 1;
    }

    //Calculate mean in pad direction:
    Int_t padmean = seqcharge*pad;
    Int_t paderror = pad*padmean;

    //Compare with results on previous pad:
    for(UInt_t p=0; p<nprevious && p<kPadArraySize && ncurrent<kPadArraySize; p++){
      
      //dont merge sequences on the same pad twice
      if(previousPt[p]->fLastMergedPad==pad) continue;

      Int_t difference = seqmean - previousPt[p]->fMean;
      if(difference < -fMatch) break;

      if(difference <= fMatch){ //There is a match here!!
	AliClusterData *local = previousPt[p];
	
	if(fDeconvPad){
	  if(seqcharge > local->fLastCharge){
	    if(local->fChargeFalling){ //The previous pad was falling
	      break; //create a new cluster
	    }		    
	  }
	  else local->fChargeFalling = 1;
	  local->fLastCharge = seqcharge;
	}
	      
	//Don't create a new cluster, because we found a match
	newcluster = kFALSE;
	      
	//Update cluster on current pad with the matching one:
	local->fTotalCharge += seqcharge;
	local->fPad += padmean;
	local->fPad2 += paderror;
	local->fTime += seqaverage;
	local->fTime2 += seqerror;
	local->fMean = seqmean;
	local->fFlags++; //means we have more than one pad 
	local->fLastMergedPad = pad;

	currentPt[ncurrent] = local;
	ncurrent++;
	      
	break;
      } //Checking for match at previous pad
    } //Loop over results on previous pad.

    if(newcluster && ncurrent<kPadArraySize){
      //Start a new cluster. Add it to the clusterlist, and update
      //the list of pointers to clusters in current pad.
      //current pad will be previous pad on next pad.

      //Add to the clusterlist:
      AliClusterData *tmp = &clusterlist[ntotal];
      tmp->fTotalCharge = seqcharge;
      tmp->fPad = padmean;
      tmp->fPad2 = paderror;
      tmp->fTime = seqaverage;
      tmp->fTime2 = seqerror;
      tmp->fMean = seqmean;
      tmp->fFlags = 0;  //flags for single pad clusters
      tmp->fLastMergedPad = pad;

      if(fDeconvPad){
	tmp->fChargeFalling = 0;
	tmp->fLastCharge = seqcharge;
      }

      //Update list of pointers to previous pad:
      currentPt[ncurrent] = &clusterlist[ntotal];
      ntotal++;
      ncurrent++;
    }

    if(fDeconvTime)
      if(newbin >= 0) goto redo;
  
    // to prevent endless loop  
    if(time >= AliHLTTPCTransform::GetNTimeBins()){
      HLTWarning("Timebin (%d) out of range (%d)", time, AliHLTTPCTransform::GetNTimeBins());
      break;
    }


    if(!readValue) break; //No more value
    
    if (ntotal>=kClusterListSize || ncurrent>=kPadArraySize) {
      HLTWarning("pad array size exceeded ntotal=%d ncurrent=%d, skip rest of the data", ntotal, ncurrent);
      break;
    }

    if(fCurrentRow != newRow){
      WriteClusters(ntotal,clusterlist);

      lastpad = 123456789;

      currentPt = pad2;
      previousPt = pad1;
      nprevious=0;
      ncurrent=0;
      ntotal=0;
      
      fCurrentRow = newRow;
    }

    pad = newPad;
    time = newTime;

  } // END while(readValue)

  WriteClusters(ntotal,clusterlist);

  HLTInfo("ClusterFinder found %d clusters in slice %d patch %d", fNClusters, fCurrentSlice, fCurrentPatch);

} // ENDEND

void AliHLTTPCClusterFinder::WriteClusters(Int_t nclusters,AliClusterData *list){
  // see header file for class documentation

  //write cluster to output pointer
  Int_t thisrow=-1,thissector=-1;
  UInt_t counter = fNClusters;
  
  for(int j=0; j<nclusters; j++)
    {



      if(!list[j].fFlags){
	if(fDoMC){
	  if(j+(Int_t)fClustersMCInfo.size()-nclusters >=0 && j+fClustersMCInfo.size()-nclusters < fClustersMCInfo.size()){
	    fClustersMCInfo.erase(fClustersMCInfo.begin()+j+fClustersMCInfo.size()-nclusters); // remove the mc info for this cluster since it is not taken into account 
	  }
	}
	continue; //discard single pad clusters
      }
      if(list[j].fTotalCharge < fThreshold){
	if(fDoMC){
	  if(j+(Int_t)fClustersMCInfo.size()-nclusters >=0 && j+fClustersMCInfo.size()-nclusters < fClustersMCInfo.size()){
	    fClustersMCInfo.erase(fClustersMCInfo.begin()+j+fClustersMCInfo.size()-nclusters); // remove the mc info for this cluster since it is not taken into account 
	  }
	}
	continue; //noise cluster
      }
      Float_t xyz[3];      
      Float_t fpad =(Float_t)list[j].fPad / list[j].fTotalCharge;
      Float_t fpad2=fXYErr*fXYErr; //fixed given error
      Float_t ftime =(Float_t)list[j].fTime / list[j].fTotalCharge;
      Float_t ftime2=fZErr*fZErr;  //fixed given error



      if(fUnsorted){
	fCurrentRow=list[j].fRow;
      }

   
      if(fCalcerr) { //calc the errors, otherwice take the fixed error 
	Int_t patch = AliHLTTPCTransform::GetPatch(fCurrentRow);
	UInt_t q2=list[j].fTotalCharge*list[j].fTotalCharge;
	//	Float_t sy2=list[j].fPad2 * list[j].fTotalCharge - list[j].fPad * list[j].fPad;
	Float_t sy2=(Float_t)list[j].fPad2 * list[j].fTotalCharge - (Float_t)list[j].fPad * list[j].fPad;
	sy2/=q2;
	if(sy2 < 0) {
	    LOG(AliHLTTPCLog::kError,"AliHLTTPCClusterFinder::WriteClusters","Cluster width")
	      <<"SigmaY2 negative "<<sy2<<" on row "<<fCurrentRow<<" "<<fpad<<" "<<ftime<<ENDLOG;
	    continue;
	} else {
	  if(!fRawSP){
	    fpad2 = (sy2 + 1./12)*AliHLTTPCTransform::GetPadPitchWidth(patch)*AliHLTTPCTransform::GetPadPitchWidth(patch);
	    if(sy2 != 0){
	      fpad2*=0.108; //constants are from offline studies
	      if(patch<2)
		fpad2*=2.07;
	    }
	  } else fpad2=sy2; //take the width not the error
	}
	//	Float_t sz2=list[j].fTime2*list[j].fTotalCharge - list[j].fTime*list[j].fTime;
	Float_t sz2=(Float_t)list[j].fTime2*list[j].fTotalCharge - (Float_t)list[j].fTime*list[j].fTime;
	sz2/=q2;
	if(sz2 < 0){
	  LOG(AliHLTTPCLog::kError,"AliHLTTPCClusterFinder::WriteClusters","Cluster width")
	    <<"SigmaZ2 negative "<<sz2<<" on row "<<fCurrentRow<<" "<<fpad<<" "<<ftime<<ENDLOG;
	  continue;
	} else {
	  if(!fRawSP){
	    ftime2 = (sz2 + 1./12)*AliHLTTPCTransform::GetZWidth()*AliHLTTPCTransform::GetZWidth();
	    if(sz2 != 0) {
	      ftime2 *= 0.169; //constants are from offline studies
	      if(patch<2)
		ftime2 *= 1.77;
	    }
	  } else ftime2=sz2; //take the width, not the error
	}
      }
      if(fStdout==kTRUE)
	HLTInfo("WriteCluster: padrow %d pad %d +- %d time +- %d charge %d",fCurrentRow, fpad, fpad2, ftime, ftime2, list[j].fTotalCharge);
      
      if(!fRawSP){
	AliHLTTPCTransform::Slice2Sector(fCurrentSlice,fCurrentRow,thissector,thisrow);

	if(fOfflineTransform == NULL){
	  AliHLTTPCTransform::Raw2Local(xyz,thissector,thisrow,fpad,ftime);
	  
	  if(xyz[0]==0) LOG(AliHLTTPCLog::kError,"AliHLTTPCClustFinder","Cluster Finder")
			  <<AliHLTTPCLog::kDec<<"Zero cluster"<<ENDLOG;
	  if(fNClusters >= fMaxNClusters)
	    {
	      LOG(AliHLTTPCLog::kError,"AliHLTTPCClustFinder::WriteClusters","Cluster Finder")
		<<AliHLTTPCLog::kDec<<"Too many clusters "<<fNClusters<<ENDLOG;
	      return;
	    }  
	
	  fSpacePointData[counter].fX = xyz[0];
	  //	fSpacePointData[counter].fY = xyz[1];
	  if(fCurrentSlice<18){
	    fSpacePointData[counter].fY = xyz[1];
	  }
	  else{
	    fSpacePointData[counter].fY = -1*xyz[1];
	  }
	  fSpacePointData[counter].fZ = xyz[2];
	}
	else{
	  Double_t x[3]={thisrow,fpad+.5,ftime}; 
	  Int_t iSector[1]={thissector};
	  fOfflineTransform->Transform(x,iSector,0,1);
	  fSpacePointData[counter].fX = x[0];
	  fSpacePointData[counter].fY = x[1];
	  fSpacePointData[counter].fZ = x[2];
	}

      } 
      else {
	fSpacePointData[counter].fX = fCurrentRow;
	fSpacePointData[counter].fY = fpad;
	fSpacePointData[counter].fZ = ftime;
      }
      
      fSpacePointData[counter].fCharge = list[j].fTotalCharge;
      fSpacePointData[counter].fPadRow = fCurrentRow;
      fSpacePointData[counter].fSigmaY2 = fpad2;
      fSpacePointData[counter].fSigmaZ2  = ftime2;

      fSpacePointData[counter].fQMax = list[j].fQMax;

      fSpacePointData[counter].fUsed = kFALSE;         // only used / set in AliHLTTPCDisplay
      fSpacePointData[counter].fTrackN = -1;           // only used / set in AliHLTTPCDisplay

      Int_t patch=fCurrentPatch;
      if(patch==-1) patch=0; //never store negative patch number
      fSpacePointData[counter].fID = counter
	+((fCurrentSlice&0x7f)<<25)+((patch&0x7)<<22);//Uli

#ifdef do_mc
      Int_t trackID[3];
      GetTrackID((Int_t)rint(fpad),(Int_t)rint(ftime),trackID);

      fSpacePointData[counter].fTrackID[0] = trackID[0];
      fSpacePointData[counter].fTrackID[1] = trackID[1];
      fSpacePointData[counter].fTrackID[2] = trackID[2];

#endif
      
      fNClusters++;
      counter++;
    }
}

// STILL TO FIX  ----------------------------------------------------------------------------

#ifdef do_mc
void AliHLTTPCClusterFinder::GetTrackID(Int_t pad,Int_t time,Int_t *trackID) const {
  // see header file for class documentation

  //get mc id
  AliHLTTPCDigitRowData *rowPt = (AliHLTTPCDigitRowData*)fDigitRowData;
  
  trackID[0]=trackID[1]=trackID[2]=-2;
  for(Int_t i=fFirstRow; i<=fLastRow; i++){
    if(rowPt->fRow < (UInt_t)fCurrentRow){
      AliHLTTPCMemHandler::UpdateRowPointer(rowPt);
      continue;
    }
    AliHLTTPCDigitData *digPt = (AliHLTTPCDigitData*)rowPt->fDigitData;
    for(UInt_t j=0; j<rowPt->fNDigit; j++){
      Int_t cpad = digPt[j].fPad;
      Int_t ctime = digPt[j].fTime;
      if(cpad != pad) continue;
      if(ctime != time) continue;

      trackID[0] = digPt[j].fTrackID[0];
      trackID[1] = digPt[j].fTrackID[1];
      trackID[2] = digPt[j].fTrackID[2];
      
      break;
    }
    break;
  }
}
#endif



void AliHLTTPCClusterFinder::WriteClusters(Int_t nclusters,AliHLTTPCClusters *list){//This is used when using the AliHLTTPCClusters class for cluster data
  // see header file for class documentation

  //write cluster to output pointer
  Int_t thisrow,thissector;
  UInt_t counter = fNClusters;
  
  for(int j=0; j<nclusters; j++)
    {
      if(!list[j].fFlags) continue; //discard single pad clusters
      if(list[j].fTotalCharge < fThreshold) continue; //noise cluster

      Float_t xyz[3];      
      Float_t fpad =(Float_t)list[j].fPad / list[j].fTotalCharge;
      Float_t fpad2=fXYErr*fXYErr; //fixed given error
      Float_t ftime =(Float_t)list[j].fTime / list[j].fTotalCharge;
      Float_t ftime2=fZErr*fZErr;  //fixed given error


      if(fCalcerr) { //calc the errors, otherwice take the fixed error 
	Int_t patch = AliHLTTPCTransform::GetPatch(fCurrentRow);
	UInt_t q2=list[j].fTotalCharge*list[j].fTotalCharge;
	Float_t sy2=list[j].fPad2 * list[j].fTotalCharge - list[j].fPad * list[j].fPad;
	sy2/=q2;
	if(sy2 < 0) {
	    LOG(AliHLTTPCLog::kError,"AliHLTTPCClusterFinder::WriteClusters","Cluster width")
	      <<"SigmaY2 negative "<<sy2<<" on row "<<fCurrentRow<<" "<<fpad<<" "<<ftime<<ENDLOG;
	    continue;
	} else {
	  if(!fRawSP){
	    fpad2 = (sy2 + 1./12)*AliHLTTPCTransform::GetPadPitchWidth(patch)*AliHLTTPCTransform::GetPadPitchWidth(patch);
	    if(sy2 != 0){
	      fpad2*=0.108; //constants are from offline studies
	      if(patch<2)
		fpad2*=2.07;
	    }
	  } else fpad2=sy2; //take the width not the error
	}
	Float_t sz2=list[j].fTime2*list[j].fTotalCharge - list[j].fTime*list[j].fTime;
	sz2/=q2;
	if(sz2 < 0){
	  LOG(AliHLTTPCLog::kError,"AliHLTTPCClusterFinder::WriteClusters","Cluster width")
	    <<"SigmaZ2 negative "<<sz2<<" on row "<<fCurrentRow<<" "<<fpad<<" "<<ftime<<ENDLOG;
	  continue;
	} else {
	  if(!fRawSP){
	    ftime2 = (sz2 + 1./12)*AliHLTTPCTransform::GetZWidth()*AliHLTTPCTransform::GetZWidth();
	    if(sz2 != 0) {
	      ftime2 *= 0.169; //constants are from offline studies
	      if(patch<2)
		ftime2 *= 1.77;
	    }
	  } else ftime2=sz2; //take the width, not the error
	}
      }
      if(fStdout==kTRUE)
	HLTInfo("WriteCluster: padrow %d pad %d +- %d time +- %d charge %d",fCurrentRow, fpad, fpad2, ftime, ftime2, list[j].fTotalCharge);

      if(!fRawSP){
	AliHLTTPCTransform::Slice2Sector(fCurrentSlice,fCurrentRow,thissector,thisrow);
	AliHLTTPCTransform::Raw2Local(xyz,thissector,thisrow,fpad,ftime);
	
	if(xyz[0]==0) LOG(AliHLTTPCLog::kError,"AliHLTTPCClustFinder","Cluster Finder")
	  <<AliHLTTPCLog::kDec<<"Zero cluster"<<ENDLOG;
	if(fNClusters >= fMaxNClusters)
	  {
	    LOG(AliHLTTPCLog::kError,"AliHLTTPCClustFinder::WriteClusters","Cluster Finder")
	      <<AliHLTTPCLog::kDec<<"Too many clusters "<<fNClusters<<ENDLOG;
	    return;
	  }  
	
	fSpacePointData[counter].fX = xyz[0];
	//	fSpacePointData[counter].fY = xyz[1];
	if(fCurrentSlice<18){
	  fSpacePointData[counter].fY = xyz[1];
	}
	else{
	  fSpacePointData[counter].fY = -1*xyz[1];
	}
	fSpacePointData[counter].fZ = xyz[2];
	
      } else {
	fSpacePointData[counter].fX = fCurrentRow;
	fSpacePointData[counter].fY = fpad;
	fSpacePointData[counter].fZ = ftime;
      }
      
      fSpacePointData[counter].fCharge = list[j].fTotalCharge;
      fSpacePointData[counter].fPadRow = fCurrentRow;
      fSpacePointData[counter].fSigmaY2 = fpad2;
      fSpacePointData[counter].fSigmaZ2  = ftime2;

      fSpacePointData[counter].fQMax = list[j].fQMax;

      fSpacePointData[counter].fUsed = kFALSE;         // only used / set in AliHLTTPCDisplay
      fSpacePointData[counter].fTrackN = -1;           // only used / set in AliHLTTPCDisplay

      Int_t patch=fCurrentPatch;
      if(patch==-1) patch=0; //never store negative patch number
      fSpacePointData[counter].fID = counter
	+((fCurrentSlice&0x7f)<<25)+((patch&0x7)<<22);//Uli

#ifdef do_mc
      Int_t trackID[3];
      GetTrackID((Int_t)rint(fpad),(Int_t)rint(ftime),trackID);

      fSpacePointData[counter].fTrackID[0] = trackID[0];
      fSpacePointData[counter].fTrackID[1] = trackID[1];
      fSpacePointData[counter].fTrackID[2] = trackID[2];

#endif
      
      fNClusters++;
      counter++;
    }
}
