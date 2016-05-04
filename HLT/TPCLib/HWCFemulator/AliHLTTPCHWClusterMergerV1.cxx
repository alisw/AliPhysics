// $Id: AliHLTTPCHWClusterMergerV1.cxx 53494 2011-12-09 06:24:43Z sgorbuno $

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Sergey Gorbunov                                       *
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

//  @file   AliHLTTPCHWClusterMergerV1.cxx
//  @author Matthias Richter, Sergey Gorbunov
//  @date   2011-11-25
//  @brief  Merger class for HLT TPC Hardware clusters
//          Handles merging of branch border clusters


#include "AliHLTTPCHWClusterMergerV1.h"
#include "AliHLTTPCHWCFSupport.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCGeometry.h"
#include "AliHLTTPCRawCluster.h"
#include <algorithm>

/** ROOT macro for the implementation of ROOT specific class methods */

AliHLTTPCHWClusterMergerV1::AliHLTTPCHWClusterMergerV1()
  : AliHLTLogging()
  , fNRows(0)
  , fNRowPads(0)
  , fNBorders(0)
  , fMapping(0)
  , fBorders()
  , fpData(NULL)
  , fRawClusterBlocks(NULL)
  , fMCBlocks(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  // constructor
  
  const int kNPatchesTotal = fkNSlices*fkNPatches;

  fRawClusterBlocks = new AliHLTComponentBlockData* [kNPatchesTotal];
  fMCBlocks = new AliHLTComponentBlockData* [kNPatchesTotal];
  for( int i=0; i<kNPatchesTotal; i++){
    fRawClusterBlocks[i] = NULL;
    fMCBlocks[i] = NULL;
  }
}

AliHLTTPCHWClusterMergerV1::~AliHLTTPCHWClusterMergerV1()
{
  // destructor
  delete[] fMapping;
  delete[] fRawClusterBlocks;
  delete[] fMCBlocks;
}

Int_t AliHLTTPCHWClusterMergerV1::Init( Bool_t processingRCU2Data )
{
  // initialisation
  
  Int_t iResult=0;

  delete[] fMapping;
  fMapping = 0;
  fNRows = AliHLTTPCGeometry::GetNRows();
  fNRowPads = 0;  
  fNBorders = 0;
  for( int i=0; i<fNRows; i++ ){
    int nPads = AliHLTTPCGeometry::GetNPads(i);	
    if( fNRowPads<nPads ) fNRowPads = nPads;
  }
  int nPadsTotal = fNRows*fNRowPads;
  fMapping = new AliHLTInt16_t [nPadsTotal];

  if( !fMapping ){
    HLTError("Can not allocate memory: %d bytes", nPadsTotal*sizeof(AliHLTInt16_t));
    return -ENOMEM;
  }
  for( int i=0; i<nPadsTotal; i++ ){
    fMapping[i] = -1;
  }

  AliHLTTPCHWCFSupport support; 
  support.SetProcessingRCU2Data( processingRCU2Data );

  for( int iPart=0; iPart<6; iPart++ ){
    const AliHLTUInt32_t *m = support.GetMapping(0,iPart);
    if( !m ){
      HLTError("Can not get mapping for partition %d", iPart);
      iResult = -1;
      continue;
    }
    int nHWAddress=m[0];
    for( int iHW=0; iHW<nHWAddress; iHW++ ){
      AliHLTUInt32_t configWord = m[iHW+1];
      int row = (configWord>>8) & 0x3F;
      int pad =  configWord & 0xFF;
      bool border = (configWord>>14) & 0x1;
      if( !border ) continue;
      row+=AliHLTTPCGeometry::GetFirstRow(iPart);
      if( row>=fNRows ) continue;
      if( pad>=AliHLTTPCGeometry::GetNPads(row) ) continue;
      fMapping[row*fNRowPads + pad] = -2-iPart;
    }
  }

  fBorders.clear();
  for( int row=0; row<fNRows; row++ ){
    for( int pad=0; pad<fNRowPads-1; pad++ ){
      AliHLTInt16_t *m = fMapping + row*fNRowPads;
      if( m[pad]<-1 && m[pad+1]<-1 ){		
	fBorders.push_back( AliBorderParam( pad+1., -m[pad]-2) ); // pad, patch
	fBorders.push_back( AliBorderParam( pad+1., -m[pad+1]-2) ); // pad, patch

	m[pad] = fNBorders;	
	for( int i=1; i<fkMergeWidth; i++ ){	  
	  if( pad-i<0 || m[pad-i]!=-1 ) break;
	  m[pad-i] = fNBorders;
	}
	m[pad+1] = fNBorders+1;
	for( int i=1; i<fkMergeWidth; i++ ){	  
	  if( pad+1+i>=fNRowPads || m[pad+1+i]!=-1 ) break;
	  m[pad+1+i] = fNBorders+1;	
	}	
	fNBorders+=2;
      }
    }
  } 

  Clear();

  return iResult;
}



Int_t AliHLTTPCHWClusterMergerV1::SetDataBlock(  AliHLTComponentBlockData *block)
{
  //
  // set block of clusters or MC labels
  //

  if( !block ){    
    HLTError("Input NULL pointer to data block");
    return -1;
  }
  
  Int_t blockType=0;
  
  if( block->fDataType == (AliHLTTPCDefinitions::fgkRawClustersDataType | kAliHLTDataOriginTPC) ) blockType=1;
  else if( block->fDataType == (AliHLTTPCDefinitions::fgkAliHLTDataTypeClusterMCInfo | kAliHLTDataOriginTPC) ) blockType=2;
  
  if( blockType==0 ){
    HLTError("Unexpected type of input block:  %d", block->fDataType);
    return -1;
  }
  
  UInt_t minSlice     = AliHLTTPCDefinitions::GetMinSliceNr(*block); 
  UInt_t minPartition = AliHLTTPCDefinitions::GetMinPatchNr(*block);
  UInt_t maxSlice     = AliHLTTPCDefinitions::GetMaxSliceNr(*block); 
  UInt_t maxPartition = AliHLTTPCDefinitions::GetMaxPatchNr(*block);
  if( maxSlice!=minSlice || maxPartition!=minPartition ){
    HLTError("Unexpected input block (slices %d-%d, patches %d-%d): only one TPC partition per block is expected", minSlice, maxSlice, minPartition, maxPartition );
    return -1;
  }

  if( minSlice>=(UInt_t)fkNSlices || minPartition >=(UInt_t)fkNPatches ){
    HLTError("Wrong Slice/Patch number of input block: slice %d, patch %d", minSlice, minPartition );
    return -1;
  }

  Int_t iSlicePatch = minSlice*fkNPatches + minPartition;
  
  if( blockType == 1 ) fRawClusterBlocks[iSlicePatch] = block;
  else if( blockType == 2 )  fMCBlocks[iSlicePatch] = block;

  return 0;
}


void AliHLTTPCHWClusterMergerV1::Clear()
{
  /// cleanup
  
  //  fClusters.~vector<AliClusterRecord>();
  // new(&fClusters) vector<AliClusterRecord>;

  fpData = NULL;
  for( int i=0; i<fkNSlices*fkNPatches; i++){
    fRawClusterBlocks[i] = NULL;
    fMCBlocks[i] = NULL;
  }
}


int AliHLTTPCHWClusterMergerV1::Merge()
{
  /// merge clusters
   
  if( !fMapping ){
    HLTError("Mapping is not initialised");
    return 0;
  }

  if( !fpData ){    
    HLTError("Pointer to input data is not set");
    return -1;
  }

  
  int iResult = 0;
  int count = 0;

  vector<AliBorderRecord> *fBorderClusters = new vector<AliBorderRecord> [fNBorders];

  for( int iSlice=0; iSlice<fkNSlices; iSlice++){
    
    for( int iBorder=0; iBorder<fNBorders; iBorder++) fBorderClusters[iBorder].clear();
    
    for( int iPatch=0; iPatch<fkNPatches; iPatch++){
      int iSlicePatch =  iSlice*fkNPatches+iPatch;
      AliHLTComponentBlockData *block = fRawClusterBlocks[iSlicePatch];
      if( !block ) continue;
      AliHLTTPCRawClusterData* clusters= (AliHLTTPCRawClusterData*)( fpData + block->fOffset);
      if(!clusters) continue;
  
      AliHLTComponentBlockData *mcblock = fMCBlocks[iSlicePatch];
      AliHLTTPCClusterMCData* mclabels = 0;
      if( mcblock ) mclabels = (AliHLTTPCClusterMCData*)(fpData + mcblock->fOffset );
      if( mclabels && mclabels->fCount != clusters->fCount ){
	HLTError("Slice %d, patch %d: Number of MC labels (%d) not equal to N clusters (%d)", iSlice, iPatch,  mclabels->fCount, clusters->fCount );
	mclabels->fCount = 0;
	mcblock->fSize = sizeof(AliHLTTPCClusterMCData);
	mclabels = 0;
	fMCBlocks[iSlicePatch] = 0;
      }
      
      for( UInt_t iCluster=0; iCluster<clusters->fCount; iCluster++){
	AliHLTTPCRawCluster &cluster = clusters->fClusters[iCluster];

	// check if the cluster is at the branch border    
	Int_t patchRow = cluster.GetPadRow();	
	int sliceRow=patchRow+AliHLTTPCGeometry::GetFirstRow(iPatch);
 	float pad = cluster.GetPad();
	int iPad = (int) pad;

	int iBorder = -1;

	if( sliceRow>=0 && sliceRow <fNRows && iPad>=0 && iPad < fNRowPads ){
	  iBorder = fMapping[sliceRow*fNRowPads+iPad];    
	  if( iBorder>=0 ){
	    float dPad = pad - fBorders[iBorder].fPadPosition;
	    if( cluster.GetSigmaPad2()>1.e-4 ){
	      if( dPad*dPad > 12.*cluster.GetSigmaPad2() ) iBorder = -1;
	    } else {
	      if( fabs(dPad)>1. ) iBorder = -1;
	    }    
	  }
	}
	if( iBorder>=0 ){
	  AliHLTTPCClusterMCLabel *mc=0;
	  if( mclabels ) mc = &( mclabels->fLabels[iCluster] );
	  fBorderClusters[iBorder].push_back( AliBorderRecord( &cluster, mc, (int)cluster.GetTime()) );
	}
      }
    } // patches
  
    
    for( int iBorder=0; iBorder<fNBorders; iBorder+=2){      
      vector<AliBorderRecord> &border1 = fBorderClusters[iBorder];
      vector<AliBorderRecord> &border2 = fBorderClusters[iBorder+1];
      UInt_t n1 = border1.size();
      UInt_t n2 = border2.size();      
      if( n1==0 || n2==0 ) continue;

      // sort 
      
      std::sort(border1.begin(),border1.end(), CompareTime);
      std::sort(border2.begin(),border2.end(), CompareTime);
      
      // merge

      UInt_t ib1 = 0;
      for( UInt_t ib2 = 0; (ib1<n1) && (ib2<n2); ib2++ ){
	
	// search first cluster at border1 to merge
	
	while( ib1<n1 && border1[ib1].fTimeBin>border2[ib2].fTimeBin+fkMergeTimeWindow ){
	  ib1++;
	}
	
	if( ib1>= n1 ) break;
	if( border1[ib1].fTimeBin < border2[ib2].fTimeBin-fkMergeTimeWindow) continue;
      
	// merge c2->c1 
	
	AliHLTTPCRawCluster *c1 = border1[ib1].fCluster;
	AliHLTTPCRawCluster *c2 = border2[ib2].fCluster;
	AliHLTTPCClusterMCLabel *mc1 = border1[ib1].fMC;
	AliHLTTPCClusterMCLabel *mc2 = border2[ib2].fMC;


	if( c1->GetPad()<-1 || c2->GetPad()<-1 ) continue;

	float w1 = c1->GetCharge();
	float w2 = c2->GetCharge();
	if( w1<=0 ) w1 = 1;
	if( w2<=0 ) w2 = 1;
	float w = 1./(w1+w2);
	w1*=w;
	w2*=w;
	
	c1->SetCharge( c1->GetCharge() + c2->GetCharge() );
	if( c1->GetQMax() < c2->GetQMax() ) c1->SetQMax( c2->GetQMax() );
	
	c1->SetSigmaPad2(  w1*c1->GetSigmaPad2() + w2*c2->GetSigmaPad2()
			 + (c1->GetPad() - c2->GetPad())*(c1->GetPad() - c2->GetPad())*w1*w2 );
	
	c1->SetSigmaTime2( w1*c1->GetSigmaTime2() + w2*c2->GetSigmaTime2()
			+ (c1->GetTime() - c2->GetTime())*(c1->GetTime() - c2->GetTime())*w1*w2 );
      
	c1->SetPad( w1*c1->GetPad() + w2*c2->GetPad() );
	c1->SetTime( w1*c1->GetTime() + w2*c2->GetTime() );
      
	// merge MC labels
	if( mc1 && mc2 ){

	  AliHLTTPCClusterMCWeight labels[6] = {
	    mc1->fClusterID[0],
	    mc1->fClusterID[1],
	    mc1->fClusterID[2],
	    mc2->fClusterID[0],
	    mc2->fClusterID[1],
	    mc2->fClusterID[2]
	  };

	  sort(labels, labels+6, CompareMCLabels);
	  for( unsigned int i=1; i<6; i++ ){
	    if(labels[i-1].fMCID==labels[i].fMCID ){
	      labels[i].fWeight+=labels[i-1].fWeight;
	      labels[i-1].fWeight = 0;
	    }
	  }
	
	  sort(labels, labels+6, CompareMCWeights );
	  
	  for( unsigned int i=0; i<3; i++ ){
	    mc1->fClusterID[i] = labels[i];
	  }

	} // MC labels

	// wipe c2
	c2->SetPad(-10);

	count++;
	// do not merge to c1 anymore
	ib1++;
      }    
    } // iBorder
    

    // remove merged clusters from data

    for( int iPatch=0; iPatch<fkNPatches; iPatch++){
      int iSlicePatch =  iSlice*fkNPatches+iPatch;
      AliHLTComponentBlockData *block = fRawClusterBlocks[iSlicePatch];
      if( !block ) continue;
      AliHLTTPCRawClusterData* clusters= (AliHLTTPCRawClusterData*)(fpData + block->fOffset);
      if(!clusters) continue;
      AliHLTUInt32_t nClustersOrig = clusters->fCount;
      
      AliHLTComponentBlockData *mcblock = fMCBlocks[iSlicePatch];
      AliHLTTPCClusterMCData* mclabels = 0;
      if( mcblock ) mclabels = (AliHLTTPCClusterMCData*)(fpData + mcblock->fOffset);
      if( mclabels && mclabels->fCount != nClustersOrig ) mclabels = 0;      

      clusters->fCount=0;
      for( UInt_t  iCluster=0; iCluster<nClustersOrig; iCluster++){
	AliHLTTPCRawCluster &cluster = clusters->fClusters[iCluster];
	if( cluster.GetPad()<-1 ) continue; // the cluster has been merged
	if( clusters->fCount != iCluster ){
	  clusters->fClusters[clusters->fCount] = cluster;
	  if( mclabels ) mclabels->fLabels[clusters->fCount] = mclabels->fLabels[iCluster];
	}
	clusters->fCount++;
      }
      
      block->fSize = sizeof(AliHLTTPCRawClusterData) + clusters->fCount*sizeof(AliHLTTPCRawCluster);
      if( mclabels ){
	mclabels->fCount = clusters->fCount;
	mcblock->fSize = sizeof(AliHLTTPCClusterMCData) + mclabels->fCount*sizeof(AliHLTTPCClusterMCLabel);	
      }
    }

  } // iSlice

  delete[] fBorderClusters;

  if (iResult<0) return iResult;
  return count;
}

