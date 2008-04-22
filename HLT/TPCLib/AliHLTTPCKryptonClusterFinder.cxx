// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Kenneth Aamodt <kenneth.aamodt@student.uib.no>        *
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

/** @file    AliHLTTPCKryptonClusterFinder.cxx
    @author  Kenneth Aamodt, Kalliopi Kanaki
    @date   
    @brief  Krypton Cluster Finder for the TPC
*/

#include "AliHLTTPCDigitReader.h"
#include "AliHLTTPCKryptonClusterFinder.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCPad.h"
#include "AliHLTTPCClusters.h"
#include "TFile.h"

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTPCKryptonClusterFinder)

AliHLTTPCKryptonClusterFinder::AliHLTTPCKryptonClusterFinder()
  : AliHLTTPCClusterFinder(),
    fTimebinsInBunch(),
    fIndexOfBunchStart(),
    fHKryptonSpectrumFullPatch(NULL),
    fHKryptonSpectrumSelection(NULL),
    fHNumberOfKryptonClusters(NULL),
    fHNumberOfKryptonClustersSelection(NULL),
    fHMaxQofKryptonClusterLast1000(NULL),
    fHMaxQofKryptonClusterSelection(NULL),
    fStartBinKryptonSpectrum(0),
    fEndBinKryptonSpectrum(3000),
    fStartBinMaxQ(0),
    fEndBinMaxQ(1000),
    fStartBinNumberOfKryptonClusters(0),
    fEndBinNumberOfKryptonClusters(1000),
    fSelectionMinRowNumber(0),
    fSelectionMaxRowNumber(0),
    fMaxQOfCluster(0),
    fMaxQOfClusterBin(0),
    fNumberOfKryptonClusters(0),
    fNumberOfKryptonClustersBin(0),
    fHistogramsInitialized(kFALSE)
{
  //constructor  
}

AliHLTTPCKryptonClusterFinder::~AliHLTTPCKryptonClusterFinder(){
  if(fHKryptonSpectrumFullPatch){
    delete fHKryptonSpectrumFullPatch;
    fHKryptonSpectrumFullPatch=NULL;
  }
  if(fHKryptonSpectrumSelection){
    delete fHKryptonSpectrumSelection;
    fHKryptonSpectrumSelection=NULL;
  }
  if(fHNumberOfKryptonClusters){
    delete fHNumberOfKryptonClusters;
    fHNumberOfKryptonClusters=NULL;
  }
  if(fHMaxQofKryptonClusterLast1000){
    delete fHMaxQofKryptonClusterLast1000;
    fHMaxQofKryptonClusterLast1000=NULL;
  }
  if(fHMaxQofKryptonClusterSelection){
    delete fHMaxQofKryptonClusterSelection;
    fHMaxQofKryptonClusterSelection=NULL;
  }
}

void AliHLTTPCKryptonClusterFinder::ReBunch(const UInt_t *bunchData , Int_t bunchSize){
  fTimebinsInBunch.clear();
  fIndexOfBunchStart.clear();
  Bool_t newBunch=kTRUE;
  Int_t currentBunchNumber=-1;
  for(Int_t i=0;i<bunchSize;i++){
    if(newBunch){
      if(bunchData[i]>5){
	fIndexOfBunchStart.push_back(i);
	fTimebinsInBunch.push_back(1);
	currentBunchNumber++;
	newBunch=kFALSE;
      }
    }
    else{
      if(bunchData[i]>0){
	fTimebinsInBunch[currentBunchNumber]++;
      }
      else{
	newBunch=kTRUE;
      }
    }
  }
}

void AliHLTTPCKryptonClusterFinder::ReadDataUnsorted(void* ptr,unsigned long size)
{
  //set input pointer
  fPtr = (UChar_t*)ptr;
  fSize = size;

  if(!fVectorInitialized){
    InitializePadArray();
   }

  fDigitReader->InitBlock(fPtr,fSize,fFirstRow,fLastRow,fCurrentPatch,fCurrentSlice);
  
  while(fDigitReader->NextChannel()){
    UInt_t row=fDigitReader->GetRow();
    UInt_t pad=fDigitReader->GetPad();

    fRowPadVector[row][pad]->ClearCandidates();

    while(fDigitReader->NextBunch()){
      if(fDigitReader->GetBunchSize()>1){
	const UInt_t *bunchData= fDigitReader->GetSignals();

	ReBunch(bunchData,fDigitReader->GetBunchSize());
	Int_t rebunchCount=fIndexOfBunchStart.size();
	
	for(Int_t r=0;r<rebunchCount;r++){
	  UInt_t time = fDigitReader->GetTime()+fIndexOfBunchStart[r];
	  AliHLTTPCClusters candidate;
	  candidate.fTotalCharge=0;
	  candidate.fChargeFalling=0;//PS This variable will store teh maximum charge of the candidate, used later. Keep this in mind 
	  if(fTimebinsInBunch[r]>2){
	    for(Int_t i=0;i<fTimebinsInBunch[r];i++){
	      candidate.fTotalCharge+=bunchData[i + fIndexOfBunchStart[r]];	
	      candidate.fTime += time*bunchData[i + fIndexOfBunchStart[r]];
	      candidate.fTime2 += time*time*bunchData[i + fIndexOfBunchStart[r]];
	      if(bunchData[i + fIndexOfBunchStart[r]]>candidate.fChargeFalling){
		candidate.fChargeFalling = bunchData[i + fIndexOfBunchStart[r]];
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
	    if(candidate.fTotalCharge>10 && candidate.fMean<924){
	      fRowPadVector[row][pad]->AddClusterCandidate(candidate);
	    }
	  }
	}
      }
    }
  }
}

Bool_t AliHLTTPCKryptonClusterFinder::ComparePads(AliHLTTPCPad *nextPad,AliHLTTPCClusters* cluster,Int_t nextPadToRead){
  //Checking if we have a match on the next pad
  for(UInt_t candidateNumber=0;candidateNumber<nextPad->fClusterCandidates.size();candidateNumber++){
    AliHLTTPCClusters *candidate =&nextPad->fClusterCandidates[candidateNumber]; 
    if(cluster->fMean-candidate->fMean==1 || candidate->fMean-cluster->fMean==1 || cluster->fMean-candidate->fMean==0){
      cluster->fMean=candidate->fMean;
      cluster->fTotalCharge+=candidate->fTotalCharge;
      cluster->fTime += candidate->fTime;
      cluster->fTime2 += candidate->fTime2;
      cluster->fPad+=candidate->fPad;
      cluster->fPad2=candidate->fPad2;
      cluster->fLastMergedPad=candidate->fPad;
      
      if(candidate->fChargeFalling>cluster->fChargeFalling){
	cluster->fChargeFalling = candidate->fChargeFalling;
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
    else{
      return kFALSE;
    }
  }
  return kFALSE;
}

void AliHLTTPCKryptonClusterFinder::FindRowClusters()
{
  // see header file for function documentation

  AliHLTTPCClusters* tmpCandidate=NULL;
  for(UInt_t row=5;row<fNumberOfRows-5;row++){
    fRowOfFirstCandidate=row;
    for(UInt_t pad=5;pad<fNumberOfPadsInRow[row]-1-5;pad++){
      AliHLTTPCPad *tmpPad=fRowPadVector[row][pad];
      for(size_t candidate=0;candidate<tmpPad->fClusterCandidates.size();candidate++){
	if(tmpPad->fUsedClusterCandidates[candidate]){
	  continue;
	}
	if((Int_t)tmpPad->fClusterCandidates[candidate].fMean<100 || (Int_t)tmpPad->fClusterCandidates[candidate].fMean>AliHLTTPCTransform::GetNTimeBins()-100){
	  continue;
	}
	tmpCandidate=&tmpPad->fClusterCandidates[candidate];
	UInt_t tmpTotalCharge=tmpCandidate->fTotalCharge;
	
	ComparePads(fRowPadVector[row][pad+1],tmpCandidate,pad+1);
	if(tmpCandidate->fTotalCharge>tmpTotalCharge){
	  //we have a cluster
	  tmpCandidate->fPad=tmpCandidate->fPad/tmpCandidate->fTotalCharge;
	  fClusters.push_back(*tmpCandidate);
	}
      }
    }
  }

  HLTDebug("Found %d normal clusters.",fClusters.size());
}

void AliHLTTPCKryptonClusterFinder::FindKryptonClusters()
{
  if(fClusters.size()<2){
    return;
  }

  for(UInt_t i=0; i<fClusters.size();i++){
    AliHLTTPCClusters * tmpCluster=&fClusters[i];
    if(tmpCluster->fFlags==99){//quickfix to check if a cluster is used already
      continue;
    }
    if((Int_t)tmpCluster->fChargeFalling>fMaxQOfCluster){
      fMaxQOfCluster = tmpCluster->fChargeFalling;
    }
    
    //adds "normal" clusters belonging to the krypton cluster
    for(UInt_t j=i+1;j<fClusters.size();j++){
      AliHLTTPCClusters * nextCluster=&fClusters[j];

      if(nextCluster->fFlags==99){//quickfix to check if a cluster is used already
	continue;
      }

      if(tmpCluster->fRowNumber==nextCluster->fRowNumber-1){//Checks if the row numbers are ok (next to eachother)
	if(abs((Int_t)(tmpCluster->fPad) - (Int_t)(nextCluster->fPad))<3){ // checks if the pad numbers are ok
	  if(abs((Int_t)tmpCluster->fMean-(Int_t)nextCluster->fMean)<2){
	    tmpCluster->fMean=nextCluster->fMean;
	    tmpCluster->fTotalCharge+=nextCluster->fTotalCharge;
	    tmpCluster->fPad=nextCluster->fPad;
	    
	    if((Int_t)nextCluster->fChargeFalling>fMaxQOfCluster){
	      fMaxQOfCluster = nextCluster->fChargeFalling;
	    }
	    
	    if(tmpCluster->fFlags!=99){//means that this is the first time normal clusters match
	      CheckForCandidateOnPreviousRow(tmpCluster);
	    }
	    tmpCluster->fRowNumber=nextCluster->fRowNumber;
	    nextCluster->fFlags=99;
	    tmpCluster->fFlags=99;
	    if(j!=fClusters.size()-1){
	      continue;
	    }
	    
	  }
	}
      }
      
      if(tmpCluster->fFlags==99){
	//adds a clustercandidate on next row if present TODO
	/*	if(tmpCluster->fFlags==99){
	  for(Int_t p=-1;p<2;p++){
	  AliHLTTPCPad *padAfter=fRowPadVector[tmpCluster->fRowNumber+1][tmpCluster->fPad+p];
	  for(UInt_t c=0;c<padAfter->fClusterCandidates.size();c++)
	  if(abs((Int_t)tmpCluster->fMean - (Int_t)padAfter->fClusterCandidates[c].fMean)<2){
	  tmpCluster->fTotalCharge+=padAfter->fClusterCandidates[c].fTotalCharge;
	  }
	  }//end add clustercandidate if present
	  }*/
	fHMaxQofKryptonClusterLast1000->Fill(fMaxQOfClusterBin,fMaxQOfCluster);
	fMaxQOfClusterBin++;
	if(fMaxQOfClusterBin>fEndBinMaxQ){
	  fMaxQOfClusterBin=0;
	}
	fHKryptonSpectrumFullPatch->Fill(tmpCluster->fTotalCharge);
	fHNumberOfKryptonClusters->Fill(fNumberOfKryptonClustersBin);
	if(fNumberOfKryptonClustersBin>fEndBinNumberOfKryptonClusters){
	  fNumberOfKryptonClustersBin=0;
	}
	/*
	fHKryptonSpectrumSelection
	fHNumberOfKryptonClustersSelection
	fHMaxQofKryptonClusterSelection
	*/
	HLTInfo("Krypton cluster found, charge: %d   in patch number: %d",tmpCluster->fTotalCharge,fCurrentPatch);
	break;
      }


    }
  }//end add "normal" clusters belonging to the krypton cluster
  fNumberOfKryptonClustersBin++;
}

void AliHLTTPCKryptonClusterFinder::CheckForCandidateOnPreviousRow(AliHLTTPCClusters *tmpCluster){
  if(tmpCluster->fRowNumber>1){
    for(Int_t p=-1;p<2;p++){
      if(tmpCluster->fPad+p>0 && tmpCluster->fPad+p<fNumberOfPadsInRow[tmpCluster->fRowNumber-1]){
	AliHLTTPCPad *prevPad=fRowPadVector[tmpCluster->fRowNumber-1-AliHLTTPCTransform::GetFirstRow(fCurrentPatch)][tmpCluster->fPad+p];
	for(UInt_t i=0;i<prevPad->fClusterCandidates.size();i++){
	  if(abs((Int_t)prevPad->fClusterCandidates[i].fMean - (Int_t)tmpCluster->fMean)<2 && prevPad->fUsedClusterCandidates[i]==0){
	    tmpCluster->fTotalCharge += prevPad->fClusterCandidates[i].fTotalCharge;
	  }
	}
      }
    }
  }
}

void AliHLTTPCKryptonClusterFinder::InitializeHistograms(){


  if(fHistogramsInitialized){
    return;
  }

  TString sliceStr("-Slice_");
  sliceStr+=fCurrentSlice;
  TString patchStr("-Patch_");
  patchStr+=fCurrentPatch;

  TString namefKryptonSpectrumFullPatch("KryptonSpectrumFullPatch"+sliceStr+patchStr);
  TString namefKryptonSpectrumSelection = "KryptonSpectrumSelection"+sliceStr+patchStr;
  TString namefNumberOfKryptonClusters = "NumberOfKryptonClusters"+sliceStr+patchStr;
  TString namefNumberOfKryptonClustersSelection = "NumberOfKryptonClustersSelection"+sliceStr+patchStr;
  TString namefMaxQLast1000 = "MaxQ"+sliceStr+patchStr;
  TString namefMaxQSelection = "MaxQSelection"+sliceStr+patchStr;

  fHKryptonSpectrumFullPatch = new TH1F(namefKryptonSpectrumFullPatch,namefKryptonSpectrumFullPatch,fEndBinKryptonSpectrum-fStartBinKryptonSpectrum,fStartBinKryptonSpectrum,fEndBinKryptonSpectrum);
 
 fHKryptonSpectrumSelection = new TH1F(namefKryptonSpectrumSelection,namefKryptonSpectrumSelection,fEndBinKryptonSpectrum-fStartBinKryptonSpectrum,fStartBinKryptonSpectrum,fEndBinKryptonSpectrum);

  fHNumberOfKryptonClusters = new TH1F(namefNumberOfKryptonClusters,namefNumberOfKryptonClusters,fEndBinNumberOfKryptonClusters-fStartBinNumberOfKryptonClusters,fStartBinNumberOfKryptonClusters,fEndBinNumberOfKryptonClusters);

  fHNumberOfKryptonClustersSelection = new TH1F(namefNumberOfKryptonClustersSelection,namefNumberOfKryptonClustersSelection,fEndBinNumberOfKryptonClusters-fStartBinNumberOfKryptonClusters,fStartBinNumberOfKryptonClusters,fEndBinNumberOfKryptonClusters);

  fHMaxQofKryptonClusterLast1000 = new TH1F(namefMaxQLast1000,namefMaxQLast1000,fEndBinMaxQ-fStartBinMaxQ,fStartBinMaxQ,fEndBinMaxQ);

  fHMaxQofKryptonClusterSelection = new TH1F(namefMaxQSelection,namefMaxQSelection,fEndBinMaxQ-fStartBinMaxQ,fStartBinMaxQ,fEndBinMaxQ);

  fHistogramsInitialized=kTRUE;
}

void AliHLTTPCKryptonClusterFinder::ResetHistograms(TString histoName){

  if (histoName.CompareTo("KryptonSpectrumFullPatch")==0){
    fHKryptonSpectrumFullPatch->Reset();
  }
  if (histoName.CompareTo("KryptonSpectrumSelection")==0){
    fHKryptonSpectrumSelection->Reset();
  }
  if (histoName.CompareTo("NumberOfKryptonClusters")==0){
    fHNumberOfKryptonClusters->Reset();
  }
  if (histoName.CompareTo("NumberOfKryptonClustersSelection")==0){
    fHNumberOfKryptonClustersSelection->Reset();
  }
  if (histoName.CompareTo("MaxQ")==0){
    fHMaxQofKryptonClusterLast1000->Reset();
  }
  if (histoName.CompareTo("MaxQSelection")==0){
    fHMaxQofKryptonClusterSelection->Reset();
  }
  if (histoName.CompareTo("All")==0){
    fHKryptonSpectrumFullPatch->Reset();
    fHKryptonSpectrumSelection->Reset();
    fHNumberOfKryptonClusters->Reset();
    fHNumberOfKryptonClustersSelection->Reset();
    fHMaxQofKryptonClusterLast1000->Reset();
    fHMaxQofKryptonClusterSelection->Reset();
  }
}

void AliHLTTPCKryptonClusterFinder::SetSelection(Int_t minRow, Int_t maxRow){
  fSelectionMinRowNumber=minRow;
  fSelectionMaxRowNumber=maxRow;
}

void AliHLTTPCKryptonClusterFinder::GetHistogramObjectArray(TObjArray & histos){
  //  TObjArray histos;
  histos.Add(fHKryptonSpectrumFullPatch);
  histos.Add(fHKryptonSpectrumSelection);
  histos.Add(fHNumberOfKryptonClusters);
  histos.Add(fHNumberOfKryptonClustersSelection);
  histos.Add(fHMaxQofKryptonClusterLast1000);
  histos.Add(fHMaxQofKryptonClusterSelection);
  //  return histos;
}

void AliHLTTPCKryptonClusterFinder::WriteHistograms(TString filename){

  TFile file(filename,"recreate");
  fHKryptonSpectrumFullPatch->Write();
  fHKryptonSpectrumSelection->Write();
  fHNumberOfKryptonClusters->Write();
  fHNumberOfKryptonClustersSelection->Write();
  fHMaxQofKryptonClusterLast1000->Write();
  fHMaxQofKryptonClusterSelection->Write();
  file.Close();
}
