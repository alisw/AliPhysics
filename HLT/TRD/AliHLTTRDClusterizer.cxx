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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// HLT TRD cluster finder                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliHLTTRDClusterizer.h"
#include "AliHLTTRDCluster.h"
#include "AliTRDgeometry.h"
#include "AliTRDcluster.h"
#include "AliTRDReconstructor.h"
#include <TClonesArray.h>

ClassImp(AliHLTTRDClusterizer)

//_____________________________________________________________________________
AliHLTTRDClusterizer::AliHLTTRDClusterizer(const AliTRDReconstructor *const rec)
  :AliTRDclusterizer(rec)
  ,fClMemBlock(NULL)
  ,fTrMemBlock(NULL)
  ,fTrMemCurrPtr(NULL)
  ,fLastDet(-1)
  ,fClusters(NULL)
  ,fAddedSize(0)
{
  //
  // AliHLTTRDClusterizer default constructor
  //
}

//_____________________________________________________________________________
AliHLTTRDClusterizer::AliHLTTRDClusterizer(const Text_t *const name, const Text_t *const title, const AliTRDReconstructor *const rec)
  : AliTRDclusterizer(name,title,rec)
  ,fClMemBlock(NULL)
  ,fTrMemBlock(NULL)
  ,fTrMemCurrPtr(NULL)
  ,fLastDet(-1)
  ,fClusters(NULL)
  ,fAddedSize(0)
{
  //
  // AliHLTTRDClusterizer constructor
  //
}

//_____________________________________________________________________________
AliHLTTRDClusterizer::AliHLTTRDClusterizer(const AliHLTTRDClusterizer& c)
  : AliTRDclusterizer(c)
  ,fClMemBlock(NULL)
  ,fTrMemBlock(NULL)
  ,fTrMemCurrPtr(NULL)
  ,fLastDet(-1)
  ,fClusters(NULL)
  ,fAddedSize(0)
{
  //
  // AliHLTTRDClusterizer copy constructor
  //
}

//_____________________________________________________________________________
AliHLTTRDClusterizer& AliHLTTRDClusterizer::operator=(const AliHLTTRDClusterizer& c)
{
  //
  // Assignment operator
  //
  
  if(this!=&c) 
    c.Copy(*this);
  return *this;
}

//_____________________________________________________________________________
void AliHLTTRDClusterizer::Copy(TObject& c) const
{
  //
  // Copy function
  //

  ((AliHLTTRDClusterizer&)c).fClMemBlock = NULL;
  ((AliHLTTRDClusterizer&)c).fTrMemBlock = NULL;
  ((AliHLTTRDClusterizer&)c).fTrMemCurrPtr = NULL;
}

//_____________________________________________________________________________
void AliHLTTRDClusterizer::AddClusterToArray(AliTRDcluster* cluster)
{
  //
  // Add a cluster to the array
  //

  if(fLastDet!=cluster->GetDetector()){
    fLastDet = cluster->GetDetector();
    fClusters = new(GetClMemBlock()+fAddedSize) AliHLTTRDClustersArray(fLastDet);
    fAddedSize += sizeof(AliHLTTRDClustersArray);
  }
  new(&fClusters->fCluster[fClusters->fCount]) AliHLTTRDClustersArray::cluster_type(cluster);
  fClusters->fCount++;
  fAddedSize += sizeof(AliHLTTRDClustersArray::cluster_type);
}

//_____________________________________________________________________________
void AliHLTTRDClusterizer::AddTrackletsToArray()
{
  //
  // Add the online tracklets of this chamber to the array
  //

  // memcpy(&(((UInt_t*)GetTrMemBlock())[fNoOfTracklets]),fTrackletContainer[0],256*sizeof(UInt_t));
  // memcpy(&(((UInt_t*)GetTrMemBlock())[fNoOfTracklets+256]),fTrackletContainer[1],256*sizeof(UInt_t));

  // fNoOfTracklets += 512;

  /*
  // this way of accessing tracklets is obsolete
  // and did not give any results since long time
  // removing it now to allow for removal of useless
  // legacy code in the clusterizer
  //
  // if tracklets are needed here, they should be taken
  // from TrackletsArray()
  //
  // jkl

  UInt_t* trackletword;
  AliHLTTRDTrackletWordArray* trklArr = new(fTrMemCurrPtr) AliHLTTRDTrackletWordArray(fDet);
  fTrMemCurrPtr += sizeof(AliHLTTRDTrackletWordArray);
  for(Int_t side=0; side<2; side++)
    {
      Int_t trkl=0;
      trackletword=fTrackletContainer[side];
      while(trackletword[trkl]>0){
	trkl++;
      }
      memcpy(fTrMemCurrPtr,fTrackletContainer[side],trkl*sizeof(UInt_t));
      fTrMemCurrPtr += trkl*sizeof(UInt_t);
      trklArr->fCount += trkl;
    }
  */

  // fTrackletContainer[0]+=256;
  // fTrackletContainer[1]+=256;
  // fNoOfTracklets += 512;

}
