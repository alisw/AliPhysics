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
Revision 1.2.4.1  2000/06/09 07:09:29  kowal2

Clustering and tracking classes are splitted from the simulation ones

Revision 1.2  2000/04/17 09:37:33  kowal2
removed obsolete AliTPCDigitsDisplay.C

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
#include "AliTPCParam.h" 
#include "AliSegmentArray.h" 
#include "AliCluster.h"
#include "AliClusters.h"
#include "AliClustersArray.h" 
#include "AliTPCClustersRow.h" 

#include "AliTPCClustersArray.h"
#include "TClonesArray.h"
#include "TDirectory.h"



//_____________________________________________________________________________

ClassImp(AliTPCClustersArray) 

AliTPCClustersArray::AliTPCClustersArray()
{
  fParam = 0;
  SetClass("AliTPCClustersRow");
}

AliTPCClustersArray::~AliTPCClustersArray()
{
  //
}



AliTPCClustersRow * AliTPCClustersArray::GetRow(Int_t sector,Int_t row) 
{
  //
  //return clusters ((AliTPCClustersRow *) per given sector and padrow
  //
  if (fParam==0) return 0;
  Int_t index = ((AliTPCParam*)fParam)->GetIndex(sector,row);  
  return (AliTPCClustersRow *)(*this)[index];
}

AliTPCClustersRow *  AliTPCClustersArray::CreateRow(Int_t sector, Int_t row)
{
  //
  //create digits row  
  //
  //if row just exist - delete it
  AliTPCParam * param = (AliTPCParam*)fParam;
  Int_t index = param->GetIndex(sector,row);  
  AliTPCClustersRow * clusters = (AliTPCClustersRow *)(*this)[index];
  if (clusters !=0) delete clusters;

  clusters = (AliTPCClustersRow *) AddSegment(index);
  if (clusters == 0) return 0;
  return clusters;
}

AliTPCClustersRow * AliTPCClustersArray::LoadRow(Int_t sector,Int_t row)
{
  //
  //return clusters ((AliTPCClustersRow *) per given sector and padrow
  //
  if (fParam==0) return 0;
  Int_t index = ((AliTPCParam*)fParam)->GetIndex(sector,row);  
  return (  AliTPCClustersRow *) LoadSegment(index);
}

Bool_t  AliTPCClustersArray::StoreRow(Int_t sector,Int_t row)
{
  //
  //return clusters ((AliTPCClustersRow *) per given sector and padrow
  //
  if (fParam==0) return 0;
  Int_t index = ((AliTPCParam*)fParam)->GetIndex(sector,row);  
  StoreSegment(index);
  return kTRUE;
}

Bool_t  AliTPCClustersArray::ClearRow(Int_t sector,Int_t row)
{
  //
  //return clusters ((AliTPCDigitsRow *) per given sector and padrow
  //
  if (fParam==0) return 0;
  Int_t index = ((AliTPCParam*)fParam)->GetIndex(sector,row);  
  ClearSegment(index);
  return kTRUE;
}



Bool_t AliTPCClustersArray::Setup(const AliDetectorParam *param)
{
  //
  //setup  function to adjust array parameters
  //
  if (param==0) return kFALSE;
  fParam = (AliDetectorParam *)param;
  return MakeArray(((AliTPCParam*)fParam)->GetNRowsTotal());

}
Bool_t AliTPCClustersArray::Update()
{
  //
  //setup  function to adjust array parameters
  //
  if (fParam ==0 ) return kFALSE;
  if (fTree!=0) return MakeDictionary( ((AliTPCParam*)fParam)->GetNRowsTotal()) ;
  ((AliTPCParam*)fParam)->Update();
  return MakeArray(((AliTPCParam*)fParam)->GetNRowsTotal());
}


/*
void AliTPCClustersArray::MakeTree()
{
  //  AliSegmentID  segment;
  if (fClusterType==0) {
    Error("AliTPCCLustersArray", "cluster type isn't adjusted");
    return;
  }
  AliClusters * psegment = (AliClusters *)NewSegment();  
  psegment->SetClass(fClusterType->GetName());  
  psegment->SetArray(100);
  if (fTree) delete fTree;
  fTree = new TTree("Segment Tree","Tree with segments");
  fBranch = fTree->Branch("Segment",psegment->IsA()->GetName(),&psegment,64000,1);
  delete psegment;
}              
*/
AliSegmentID * AliTPCClustersArray::NewSegment()
{
  //
  //create object according class information 
  if (fClusterType==0) {
    Error("AliTPCCLustersArray", "cluster type isn't adjusted");
    return 0;
  }
  AliSegmentID *segment=AliSegmentArray::NewSegment();
  ((AliTPCClustersRow*)segment)->SetClass(fClusterType->GetName()); 
  ((AliTPCClustersRow*)segment)->SetArray(100);
  return segment;
}
