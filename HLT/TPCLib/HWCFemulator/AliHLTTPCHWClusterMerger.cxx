// $Id$

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

//  @file   AliHLTTPCHWClusterMerger.cxx
//  @author Matthias Richter, Sergey Gorbunov
//  @date   2011-11-25
//  @brief  Merger class for HLT TPC Hardware clusters
//          Handles merging of branch border clusters

#include <algorithm>

#include "AliHLTTPCHWClusterMerger.h"
#include "AliHLTTPCGeometry.h"
#include "AliHLTTPCHWCFSupport.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCHWClusterMerger)

AliHLTTPCHWClusterMerger::AliHLTTPCHWClusterMerger()
  : AliHLTLogging()
  , fMapping(0)
  , fNRows(0)
  , fNRowPads(0)
  , fNBorders(0)
  , fBorders(0)
  , fBorderNClusters(0)
  , fBorderFirstCluster(0)
  , fBorderClusters(0)
  , fBorderNClustersTotal(0)
  , fClusters()
  , fRemovedClusterIds()
  , fIter()
  , fEnd()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  // constructor
}

AliHLTTPCHWClusterMerger::~AliHLTTPCHWClusterMerger()
{
  // destructor
  delete[] fMapping;
  delete [] fBorders;
  delete[] fBorderNClusters;
  delete[] fBorderFirstCluster;
  delete[] fBorderClusters;
}

void AliHLTTPCHWClusterMerger::Init()
{
  // initialisation
  
  delete[] fMapping;
  delete[] fBorders;
  delete[] fBorderNClusters;
  delete[] fBorderFirstCluster;
  delete[] fBorderClusters;
  fMapping = 0;
  fBorders = 0;
  fBorderNClusters = 0;
  fBorderFirstCluster = 0;
  fBorderClusters = 0;
  fBorderNClustersTotal = 0;
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
    return;
  }
  for( int i=0; i<nPadsTotal; i++ ){
    fMapping[i] = -1;
  }

  AliHLTTPCHWCFSupport support; 

  for( int iPart=0; iPart<6; iPart++ ){
    const AliHLTUInt32_t *m = support.GetMapping(0,iPart);
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
      fMapping[row*fNRowPads + pad] = -2;
    }
  }

  std::vector<AliHLTFloat32_t> vBorders;
  for( int row=0; row<fNRows; row++ ){
    for( int pad=0; pad<fNRowPads-1; pad++ ){
      AliHLTInt16_t *m = fMapping + row*fNRowPads;
      if( m[pad]==-2 && m[pad+1]==-2 ){	
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
	vBorders.push_back(pad+1.);
	vBorders.push_back(pad+1.);
      }
    }
  } 

  //cout<<" NBorders = "<<fNBorders/2<<endl;

  /*  
  cout<<"Borders:"<<endl;
  for( int row=0; row<fNRows; row++ ){
    for( int pad=0; pad<fNRowPads; pad++ ){
      AliHLTInt16_t *m = fMapping + row*fNRowPads;
      if( m[pad]>=0 ) {
	cout<<row<<" "<<pad<<" "<<m[pad]%2<<endl;
      }
    }
  }
  */

  fBorders = new AliHLTFloat32_t [fNBorders];
  if( !fBorders ){
    HLTError("Can not allocate memory: %d bytes", fNBorders*sizeof(AliHLTFloat32_t));
    fNBorders = 0;
    return;
  }
  for( int i=0; i<fNBorders; i++ ){
    fBorders[i] = vBorders[i];
  }
  fBorderNClusters = new int[fkNSlices*fNBorders];
  if( !fBorderNClusters ){
    HLTError("Can not allocate memory: %d bytes", fkNSlices*fNBorders*sizeof(int));    
    return;
  }

  fBorderFirstCluster = new int[fkNSlices*fNBorders];
  if( !fBorderFirstCluster ){
    HLTError("Can not allocate memory: %d bytes", fkNSlices*fNBorders*sizeof(int));    
    return;
  }

  for( int i=0; i<fkNSlices*fNBorders; i++){
    fBorderNClusters[i] = 0;  
    fBorderFirstCluster[i] = 0;  
  }
}

bool AliHLTTPCHWClusterMerger::CheckCandidate(int /*slice*/, int partition, int partitionrow, float pad, float /*time*/, float sigmaPad2
) 
{
  /// check cluster if it is a candidate for merging
  int slicerow=partitionrow+AliHLTTPCGeometry::GetFirstRow(partition);
  
  if( !fMapping ) Init();
  if( !fMapping ) return 0;
  if( slicerow <0 || slicerow>= fNRows ) return 0;
  int iPad = (int) pad;
  if( iPad<0 || iPad>=fNRowPads ) return 0;
  int atBorder = fMapping[slicerow*fNRowPads+iPad];    
  if( atBorder<0 ) return 0;
  float dPad = pad - fBorders[atBorder];
  if( sigmaPad2>1.e-4 ){
    if( dPad*dPad > 12.*sigmaPad2 ) return 0;
  } else {
    if( fabs(dPad)>1. ) return 0;
  }
  return 1;
}

int AliHLTTPCHWClusterMerger::AddCandidate(int slice,
					   int partition,
					   short partitionrow,
					   float pad,
					   float time,
					   float sigmaPad2,
					   float sigmaTime2,
					   unsigned short charge,
					   unsigned short qmax,
					   unsigned short flags,
					   AliHLTUInt32_t id,
					   const AliHLTTPCClusterMCLabel *mc 
					   )
{
  /// add a candidate for merging and register in the index grid
    
  int iBorder = -1;
  int iPad = (int) pad;
  if( !fMapping ) Init();
 
  int slicerow=partitionrow+AliHLTTPCGeometry::GetFirstRow(partition);
 
  if (id!=~AliHLTUInt32_t(0)) {
    if (slice<0) {
      slice=AliHLTTPCGeometry::CluID2Slice(id);
    } else if ((unsigned)slice!=AliHLTTPCGeometry::CluID2Slice(id)) {
      HLTError("cluster id 0x%08x is not consistent with specified slice %d", id, slice);
    }
    if (partition<0) {
      partition=AliHLTTPCGeometry::CluID2Partition(id);
    } else if ((unsigned)partition!=AliHLTTPCGeometry::CluID2Partition(id)) {
      HLTError("cluster id 0x%08x is not consistent with specified partition %d", id, partition);
    }
  }

  if( slice<0 || slice>=fkNSlices ){
    HLTError("cluster id 0x%08x has wrong slice number %d", id, slice);    
  }
  
  if( partition<0 || partition>=6 ){
    HLTError("cluster id 0x%08x has wrong partition number %d", id, partition);
  }
  
  if( slice>=0 && slice<fkNSlices && slicerow>=0 && slicerow <fNRows && fMapping && iPad>=0 && iPad < fNRowPads ){
    iBorder = fMapping[slicerow*fNRowPads+iPad];
    
    if( iBorder>=0 ){
      float dPad = pad - fBorders[iBorder];
      if( sigmaPad2>1.e-4 ){
	if( dPad*dPad > 12.*sigmaPad2 ) iBorder = -1;
      } else {
	if( fabs(dPad)>1. ) iBorder = -1;
      }    
    }
    if( iBorder>=0 ) iBorder = slice*fNBorders + iBorder;
  }

  fClusters.push_back(AliClusterRecord(slice, partition, iBorder, -1, id,
				       AliHLTTPCRawCluster(partitionrow, pad, time, sigmaPad2, sigmaTime2, charge, qmax, flags),
				       mc!=NULL?*mc:AliHLTTPCClusterMCLabel() ));

  if( iBorder>=0 ){
    fBorderNClusters[iBorder]++;
    fBorderNClustersTotal++;
  }
  return fClusters.size()-1;
}

int AliHLTTPCHWClusterMerger::FillIndex()
{
  /// loop over cached raw clusters and fill index

  delete[] fBorderClusters;
  fBorderClusters = 0;
  fBorderClusters = new AliBorderRecord[fBorderNClustersTotal];
  
  fBorderFirstCluster[0] = 0;
  for(int i=1; i<fkNSlices*fNBorders; i++ ){
    fBorderFirstCluster[i] = fBorderFirstCluster[i-1] + fBorderNClusters[i-1];
    fBorderNClusters[i-1] = 0;
  }
  fBorderNClusters[fkNSlices*fNBorders-1]=0;
  for (AliHLTUInt32_t pos=0; pos<fClusters.size(); pos++) {
    int iBorder = fClusters[pos].GetBorder();
    if( iBorder<0 ) continue;
    AliBorderRecord &b = fBorderClusters[fBorderFirstCluster[iBorder]+fBorderNClusters[iBorder]];
    b.fClusterRecordID = pos;
    b.fTimeBin = (int)fClusters[pos].GetCluster().GetTime();
    fBorderNClusters[iBorder]++;
  }
  /*
  cout<<"Border clusters: "<<endl;
  for( int iSlice=0; iSlice<fkNSlices; iSlice++){
    for( int ib = 0; ib<fNBorders; ib++ ){
      int ind = iSlice*fNBorders+ib;
      cout<<iSlice<<" "<<ib<<" :"<<fBorderFirstCluster[ind]<<", "<<fBorderNClusters[ind]<<endl;    
    }
  }
  */
  /*
  for( int ib=1; ib<fkNSlices*fNBorders; ib++ ){
    if( fBorderFirstCluster[ib] != fBorderFirstCluster[ib-1]+fBorderNClusters[ib-1] ){
      cout<<"Something wrong with cluster borders !!! "<<endl;
    }
  }
  */
  return 0;
}

void AliHLTTPCHWClusterMerger::Clear()
{
  /// cleanup
  fClusters.clear();
  fRemovedClusterIds.clear();
  
  delete[] fBorderClusters;
  fBorderClusters = 0;
  for( int i=0; i<fkNSlices*fNBorders; i++){
    fBorderNClusters[i] = 0;  
    fBorderFirstCluster[i] = 0; 
  }
  fBorderNClustersTotal = 0;
}


int AliHLTTPCHWClusterMerger::Merge()
{
  /// merge clusters
  /// first all candidates are filled into the index grid, looping over all clusters
  /// in the grid automatically results in time ordered clusters per row (ordering with
  /// respect to cells, within a cell two clusters can be reversed)
  
  if( !fMapping ) Init();
  if( !fMapping ) return 0;
  //static int sLeft = 0, sRight = 0, sMerged = 0;
  int iResult=0;
  int count=0;
  if ((iResult=FillIndex())<0) {
    return iResult;
  }
  // sort 
  for( int i=0; i<fkNSlices*fNBorders; i++){
    std::sort(fBorderClusters+fBorderFirstCluster[i],fBorderClusters+fBorderFirstCluster[i]+fBorderNClusters[i], CompareTime);
  }

  for( int iBorder=0; iBorder<fkNSlices*fNBorders; iBorder+=2){
    AliBorderRecord *border1 = fBorderClusters+fBorderFirstCluster[iBorder];
    AliBorderRecord *border2 = fBorderClusters+fBorderFirstCluster[iBorder+1];
    int n1 = fBorderNClusters[iBorder];
    int n2 = fBorderNClusters[iBorder+1];
    /*
sLeft+=n1;
    sRight+=n2;
    bool prn = (n1>0) && (n2>0);
    if( prn ){
      cout<<"Border "<<iBorder/2<<", left: "<<n1<<" clusters :"<<endl;
      for( int ib=0; ib<n1; ib++ ){
	cout<<ib<<": "<<border1[ib].fTimeBin<<endl;
      }
      cout<<"Border "<<iBorder/2<<", right: "<<n2<<" clusters :"<<endl;
      for( int ib=0; ib<n2; ib++ ){
	cout<<ib<<": "<<border2[ib].fTimeBin<<endl;
      }
    }
    */
    int ib1 = 0;
    for( int ib2 = 0; (ib1<n1) && (ib2<n2); ib2++ ){
      // search first cluster at border1 to merge
      
      while( ib1<n1 && border1[ib1].fTimeBin>border2[ib2].fTimeBin+fkMergeTimeWindow ){
	ib1++;
      }
      
      if( ib1>= n1 ) break;
      if( border1[ib1].fTimeBin < border2[ib2].fTimeBin-fkMergeTimeWindow) continue;

      // merge c2->c1 
      //if( prn ) cout<<"merge "<<ib2<<"->"<<ib1<<endl;
      AliClusterRecord &c1 = fClusters[border1[ib1].fClusterRecordID];
      AliClusterRecord &c2 = fClusters[border2[ib2].fClusterRecordID];
      
      float w1 = c1.Cluster().fCharge;
      float w2 = c2.Cluster().fCharge;
      if( w1<=0 ) w1 = 1;
      if( w2<=0 ) w2 = 1;
      float w = 1./(w1+w2);
      w1*=w;
      w2*=w;
      
      c1.Cluster().fCharge += c2.Cluster().fCharge;
      if( c1.Cluster().fQMax < c2.Cluster().fQMax ) c1.Cluster().fQMax = c2.Cluster().fQMax;
      
      c1.Cluster().fSigmaPad2 = 
	w1*c1.Cluster().fSigmaPad2 + w2*c2.Cluster().fSigmaPad2
	+ (c1.Cluster().fPad - c2.Cluster().fPad)*(c1.Cluster().fPad - c2.Cluster().fPad)*w1*w2;
      
      c1.Cluster().fSigmaTime2 = 
	w1*c1.Cluster().fSigmaTime2 + w2*c2.Cluster().fSigmaTime2
	+ (c1.Cluster().fTime - c2.Cluster().fTime)*(c1.Cluster().fTime - c2.Cluster().fTime)*w1*w2;
      
      c1.Cluster().fPad  = w1*c1.Cluster().fPad + w2*c2.Cluster().fPad;
      c1.Cluster().fTime = w1*c1.Cluster().fTime + w2*c2.Cluster().fTime;
      
      c1.Cluster().fFlags |= c2.Cluster().fFlags;
      
      // merge MC labels
      
      AliHLTTPCClusterMCWeight labels[6] = {
	c1.GetMCLabel().fClusterID[0],
	c1.GetMCLabel().fClusterID[1],
	c1.GetMCLabel().fClusterID[2],
	c2.GetMCLabel().fClusterID[0],
	c2.GetMCLabel().fClusterID[1],
	c2.GetMCLabel().fClusterID[2]
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
	c1.MCLabel().fClusterID[i] = labels[i];
      }

      // wipe c2
      fRemovedClusterIds.push_back(c2.GetId());
      HLTDebug("merging %d into %d", border2[ib2].fClusterRecordID, border1[ib1].fClusterRecordID);
      //cout<<"Set merged flag at position "<<border2[ib2].fClusterRecordID<<endl;
      c2.SetMergedTo(border1[ib1].fClusterRecordID);      
      count++;
      // do not merge c1 anymore
      ib1++;    
    }    
  } // iBorder
  
  //cout<<"Merged "<<count<<" clusters "<<" out of "<<fClusters.size()<<endl;
  //sMerged+= count;
  // cout<<"L "<<sLeft<<" r "<<sRight<<" m "<<sMerged<<endl;
  if (iResult<0) return iResult;
  return count;
}

