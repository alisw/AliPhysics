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
$Log$
Revision 1.1.4.2  2000/04/10 11:34:02  kowal2

Clusters handling in a new data structure

*/


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Projection Chamber clusters objects                                //
//
//  Origin: Marian Ivanov , GSI Darmstadt
//                                                                           //
//                                                                           //
//                                                                          //
///////////////////////////////////////////////////////////////////////////////
#include "AliTPC.h"
#include "AliTPCCluster.h"
#include "TClonesArray.h"
#include "TDirectory.h"


const Int_t kDefSize = 1;  //defalut size


ClassImp(AliTPCClustersRow) 


//*****************************************************************************
//
//_____________________________________________________________________________
AliTPCClustersRow::AliTPCClustersRow() 
{  
  //
  //default constructor
  fNclusters=0;
  fClusters =  new TClonesArray("AliTPCcluster",kDefSize); 
}
//_____________________________________________________________________________
AliTPCClustersRow::AliTPCClustersRow(Int_t size) 
{  
  fNclusters=0;
  fClusters = new TClonesArray("AliTPCcluster",size);
}

//_____________________________________________________________________________
const  AliTPCcluster* AliTPCClustersRow::operator[](Int_t i)
{
  //
  // return cluster at internal position i
  //
  if (fClusters==0) return 0;
  void * cl = fClusters->UncheckedAt(i);
  if (cl==0) return 0;
  return  (AliTPCcluster*)cl;
}
//_____________________________________________________________________________
void  AliTPCClustersRow::Sort()
{
  // sort cluster 
  if (fClusters) fClusters->Sort();
}

//_____________________________________________________________________________
void AliTPCClustersRow::InsertCluster(const AliTPCcluster * c) 
{ 
  //
  // Add a simulated cluster copy to the list
  //
  if(!fClusters) fClusters=new TClonesArray("AliTPCcluster",1000);
  TClonesArray &lclusters = *fClusters;
  new(lclusters[fNclusters++]) AliTPCcluster(*c);
}

//_____________________________________________________________________________
Int_t AliTPCClustersRow::Find(Double_t y) const 
{
  //
  // return index of cluster nearest to given y position
  //
  AliTPCcluster* cl;
  cl=(AliTPCcluster*)fClusters->UncheckedAt(0);
  if (y <= cl->fY) return 0;  
  cl=(AliTPCcluster*)fClusters->UncheckedAt(fNclusters-1);
  if (y > cl->fY) return fNclusters; 
  //if (y <= clusters[0]->fY) return 0;
  //if (y > clusters[num_of_clusters-1]->fY) return num_of_clusters;  
  Int_t b=0, e=fNclusters-1, m=(b+e)/2;
  //  Int_t b=0, e=num_of_clusters-1, m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    cl = (AliTPCcluster*)fClusters->UncheckedAt(m);
    if (y > cl->fY) b=m+1;
    //    if (y > clusters[m]->fY) b=m+1;
    else e=m; 
  }
  return m;
}


//_____________________________________________________________________________

ClassImp(AliTPCClustersArray) 

AliTPCClustersArray::AliTPCClustersArray()
{
  fParam = 0;
}

AliTPCClustersArray::~AliTPCClustersArray()
{
  //
  //object is only owner of fParam
  //
  if (fParam) delete fParam;
}

const AliTPCClustersRow * AliTPCClustersArray::GetRow(Int_t sector,Int_t row)
{
  //
  //return clusters ((AliTPCClustersRow *) per given sector and padrow
  //
  if (fParam==0) return 0;
  Int_t index = fParam->GetIndex(sector,row);  
  return (AliTPCClustersRow *)(*this)[index];
}

Bool_t  AliTPCClustersArray::LoadRow(Int_t sector,Int_t row)
{
  //
  //return clusters ((AliTPCClustersRow *) per given sector and padrow
  //
  if (fParam==0) return 0;
  Int_t index = fParam->GetIndex(sector,row);  
  return LoadSegment(index);
}

Bool_t  AliTPCClustersArray::StoreRow(Int_t sector,Int_t row)
{
  //
  //return clusters ((AliTPCClustersRow *) per given sector and padrow
  //
  if (fParam==0) return 0;
  Int_t index = fParam->GetIndex(sector,row);  
  StoreSegment(index);
  return kTRUE;
}



Bool_t AliTPCClustersArray::Setup(AliTPCParam *param)
{
  //
  //setup  function to adjust array parameters
  //
  if (param==0) return kFALSE;
  fParam = new AliTPCParam(*param);
  return MakeArray(fParam->GetNRowsTotal());
}
