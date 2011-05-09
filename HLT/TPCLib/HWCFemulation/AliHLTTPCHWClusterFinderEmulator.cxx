// @(#) $Id: AliHLTTPCHWClusterFinderEmulator.cxx 45447 2010-11-14 21:05:13Z sgorbuno $

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

//  @file   AliHLTTPCHWClusterFinderEmulator.cxx
//  @author Kenneth Aamodt, Kalliopi Kanaki
//  @date   
//  @brief  Cluster Finder for the TPC
//  @note 

#include "AliHLTTPCHWClusterFinderEmulator.h"
#include "AliHLTTPCClusterMCLabel.h"


//#include <sys/time.h>
#include <algorithm>
//#include <cmath>

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTTPCMapping.h"


AliHLTTPCHWClusterFinderEmulator::AliHLTTPCHWClusterFinderEmulator()
  :
  fDeconvTime(0),
  fDeconvPad(0),  
  fSlice(0),
  fPatch(0),
  fMapping(0)
{
  //constructor 
}

AliHLTTPCHWClusterFinderEmulator::~AliHLTTPCHWClusterFinderEmulator()
{   
  //destructor
  delete fMapping;
}

AliHLTTPCHWClusterFinderEmulator::AliHLTTPCHWClusterFinderEmulator(const AliHLTTPCHWClusterFinderEmulator&)
  :
  fDeconvTime(0),
  fDeconvPad(0),  
  fSlice(0),
  fPatch(0),
  fMapping(0)
{
  // dummy
}

AliHLTTPCHWClusterFinderEmulator& AliHLTTPCHWClusterFinderEmulator::operator=(const AliHLTTPCHWClusterFinderEmulator&)
{
  // dummy
  return *this;
}

void AliHLTTPCHWClusterFinderEmulator::Init( int slice, int patch )
{
  fSlice = slice;
  fPatch = patch;
  delete fMapping;
  fMapping = new AliHLTTPCMapping(fPatch);
}

 
struct AliHLTTPCHWCFClusterFragment
{
  AliHLTTPCHWCFClusterFragment():fQ(0),fT(0),fT2(0),fP(0),fP2(0),fTMean(0),fMC(){}
  AliHLTUInt64_t fQ, fT, fT2, fP, fP2, fTMean;
  vector<AliHLTInt32_t> fMC;
};

int AliHLTTPCHWClusterFinderEmulator::FindClusters( const AliHLTUInt32_t *rawEvent,
						    AliHLTUInt64_t rawEventSize32,
						    AliHLTUInt32_t *output,
						    AliHLTUInt64_t &outputSize32,
						    const AliHLTInt32_t *mcLabels,
						    AliHLTTPCClusterMCData *outputMC
						    )
{
  enum{ kWaitForNextChannel, kReadChannel, kForwardRCU, kStop };

  AliHLTTPCHWCFClusterFragment emptyFragment;
  emptyFragment.fQ = 0;
  emptyFragment.fT = 0;
  emptyFragment.fT2 = 0;
  emptyFragment.fP = 0;
  emptyFragment.fP2 = 0;
  emptyFragment.fTMean = 0;

  AliHLTTPCClusterMCWeight emptyWeight = {-1,0};

  AliHLTTPCClusterMCLabel emptyMC = {{emptyWeight,emptyWeight,emptyWeight}};


  AliHLTUInt64_t maxOutputSize32 = outputSize32;
  outputSize32 = 0;
  if( outputMC ) outputMC->fCount = 0;   

  if( !fMapping ) return -1;    
  if( !rawEvent ) return -2;    

  // Initialise 

  AliHLTTPCHWCFClusterFragment *arrFragments = new AliHLTTPCHWCFClusterFragment[kMaxNTimeBins*2];
  if( !arrFragments ) return -3;

  AliHLTTPCHWCFClusterFragment *searchRange = arrFragments;
  AliHLTTPCHWCFClusterFragment *insertRange = arrFragments + kMaxNTimeBins;

  Int_t searchNFragments = 0;
  Int_t insertNFragments = 0;
  Int_t insertRow = 0;
  Int_t insertPad = 0;
  Int_t insertBranch = -1;
  Int_t channelNWords = 0;
  Int_t channelIWord = 0;
  Int_t bunchNSignalsLeft = 0;
  Int_t bunchTime = 0;
  Int_t iMCLabel = 0;
  Int_t status = kWaitForNextChannel;  

  int ret = 0;

  // Read the data, word by word 

  for( AliHLTUInt64_t  iWord=0; iWord<=rawEventSize32; iWord++ ){
    
    Int_t newRow = insertRow;
    Int_t newPad = insertPad;
    Bool_t newBranch = insertBranch;

    Int_t nRangesToFlush = 0;

    AliHLTUInt32_t word = 0; 

    if( iWord==rawEventSize32 ){ // end of data
      nRangesToFlush = 2;
      status = kStop;
    } else {
      word = rawEvent[iWord];
    }
    
    UShort_t flag = (word >> 30);
   
    if( flag > 1 && (status==kWaitForNextChannel || status==kReadChannel) ){
      
      //  RCU trailer received, stop the clusterfinder and forward the trailer
      //cout<<"RCU trailer received, stop the clusterfinder and forward the trailer"<<endl;
      nRangesToFlush = 2;
      status = kForwardRCU;

    } else if( flag==1 && (status==kWaitForNextChannel || status==kReadChannel) ){ 

      // header of a new channel
      Short_t  hwAddress = word & 0xFFF;
      newBranch = (hwAddress >> 11) & 0x1;
      newRow = fMapping->GetRow(hwAddress);
      newPad = fMapping->GetPad(hwAddress);      
    
      //cout<<"header of a new channel: "<<newRow<<" "<<newPad<<endl;

     // merge and flush search range, flush insert range in case there is a gap in pads

      nRangesToFlush = 1;
      if( newRow!=insertRow || newPad!=insertPad+1 || newBranch!=insertBranch ) nRangesToFlush = 2;

      channelNWords = (word >> 16) & 0x3FF; // payload size in 10-bit words
      channelIWord = 0;
      status = kReadChannel;     
      if( (word >> 29) & 0x1 ) status = kWaitForNextChannel; // there were readout errors      
    
    } else if( flag==0 && status==kReadChannel ){ 

      // bunch data, read three 10-bit words

      for( int ishift=20; ishift>=0; ishift-=10 ){ 
	
	if( ++channelIWord > channelNWords ){ // look for next channel	
	  //cout<<"wrong N channel words"<<endl;
	  status = kWaitForNextChannel;
	  break;
	}
	
	AliHLTTPCHWCFClusterFragment &f = insertRange[insertNFragments];
	
	AliHLTUInt64_t word10 = (word >> ishift) & 0x3FF;
	
	if( bunchNSignalsLeft <=0 ){ // start new bunch	
	  f = emptyFragment;
	  bunchNSignalsLeft = ((Int_t)word10) - 1;
	  bunchTime = -1;
	} else { // add data to the current bunch	
	  if( bunchTime <0 ){    // time has not been read so far	
	    bunchTime = word10;
	    bunchNSignalsLeft--;
	    f.fTMean = bunchTime + ( 1 - bunchNSignalsLeft )/2;
	  } else { // add the signal
	    AliHLTUInt64_t q = word10<<kFixedPoint;
	    f.fQ += q;
	    f.fT += q*bunchTime;
	    f.fT2+= q*bunchTime*bunchTime;
	    f.fP += q*insertPad;
	    f.fP2+= q*insertPad*insertPad;
	    if( mcLabels ) f.fMC.push_back(mcLabels[iMCLabel++]);
	    bunchNSignalsLeft--;
	    bunchTime--;
	    if( bunchNSignalsLeft<=0 ){
	      //cout<<"finish bunch "<<f.fP/f.fQ<<" "<<f.fT/f.fQ<<" "<<f.fQ<<endl;
	      insertNFragments++; // current bunch is finished
	    }
	  } 
	}
      }
    } 
    
    
    
    if( nRangesToFlush>0 ){
      
      if( bunchNSignalsLeft>0 ){ // Should not happen, the data is corrupted. Finish the last bunch
	insertRange[insertNFragments].fTMean += bunchNSignalsLeft /2; 
	insertNFragments++;
	bunchNSignalsLeft = 0;
      }
      
      // do 2D merging and write the output
      
      for( int iPad=0; iPad<nRangesToFlush; iPad++ ){        
	
	int iSearch = 0;
	
	//cout<<"Flush row "<<insertRow<<" pad "<<insertPad-1+iPad<<":"<<endl;
	for( int iInsert=0; iInsert<=insertNFragments && iSearch<searchNFragments; iInsert++){
	  
	  // last iteration (iInsert==insertNFragments) to store fragments from the search range
	  
	  AliHLTTPCHWCFClusterFragment &f = insertRange[iInsert];
	  
	  while(iSearch<searchNFragments){
	    
	    AliHLTTPCHWCFClusterFragment &s = searchRange[iSearch];
	    
	    if( (iInsert==insertNFragments) || (s.fTMean > f.fTMean) ){ // store the fragment
	      //cout<<"store "<< s.fP/s.fQ<<" "<<s.fT/s.fQ<<" "<<s.fQ<<endl;
	      if( outputSize32+5 <= maxOutputSize32 ){
		AliHLTUInt32_t *c = &output[outputSize32];
		AliHLTUInt64_t q = s.fQ>>kFixedPoint;
		if( q>0 ){
		  c[0] = (((AliHLTUInt32_t) 0x3)<<30) + ((insertRow &0x3f)<<24) + ((s.fQ>>(kFixedPoint-6))&0xFFFFFF);		  
		  *((AliHLTFloat32_t*)&c[1]) = ((AliHLTFloat32_t)(s.fP/q))/(1<<kFixedPoint);
		  *((AliHLTFloat32_t*)&c[2]) = ((AliHLTFloat32_t)(s.fT/q))/(1<<kFixedPoint);
		  *((AliHLTFloat32_t*)&c[3]) = ((AliHLTFloat32_t)(s.fT2/q))/(1<<kFixedPoint);
		  *((AliHLTFloat32_t*)&c[4]) = ((AliHLTFloat32_t)(s.fP2/q))/(1<<kFixedPoint);
		  outputSize32+=5;
		  if( mcLabels && outputMC ){		    
		    AliHLTTPCClusterMCLabel &mc = outputMC->fLabels[outputMC->fCount];
		    mc = emptyMC;
		    sort(s.fMC.begin(), s.fMC.end() );
		    int ilab = -1;
		    int jlab = -1;
		    for( UInt_t i=0; i<s.fMC.size(); i++ ){
		      if( s.fMC[i]==-1 ) continue; 
		      if( ilab>=0 && s.fMC[i]==jlab ){
			mc.fClusterID[ilab].fWeight++;
		      } else {
			if( ilab>=2 ) break;
			ilab++;
			jlab = s.fMC[i];
			mc.fClusterID[ilab].fMCID = jlab;
			mc.fClusterID[ilab].fWeight = 1;
		      }
		    }
		    outputMC->fCount++;
		  }
		}
	      } else ret = -4;// No space in the output buffer

	    } else if( s.fTMean == f.fTMean ){ // merge two fragments
	      //cout<<"merge "<< s.fP/s.fQ<<" "<<s.fT/s.fQ<<" "<<s.fQ<<endl;

	      f.fQ += s.fQ;
	      f.fT += s.fT;
	      f.fT2 += s.fT2;
	      f.fP += s.fP;
	      f.fP2 += s.fP2;
	      f.fMC.insert(f.fMC.end(), s.fMC.begin(), s.fMC.end());
	    } else break; // go to the next fragment in the insert range

	    iSearch++;
	  }
	}
	
	// swap insert and search ranges

	searchNFragments = insertNFragments;
	insertNFragments = 0;
	AliHLTTPCHWCFClusterFragment *tmp = searchRange;
	searchRange = insertRange;
	insertRange = tmp;
      }
    }

    if( status==kForwardRCU && flag>1 && outputSize32<maxOutputSize32 ){
      output[outputSize32++] = word;      
    }

    insertRow = newRow;
    insertPad = newPad;
    insertBranch = newBranch;
  }
  
  delete[] arrFragments;

  return ret;
}

