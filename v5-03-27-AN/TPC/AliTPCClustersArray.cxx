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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Projection Chamber clusters objects                                //
//
//  Origin: Marian Ivanov , GSI Darmstadt
//                                                                           //
//                                                                           //
//                                                                          //
///////////////////////////////////////////////////////////////////////////////
#include "AliTPCParam.h" 
#include "AliSegmentArray.h" 
#include "AliComplexCluster.h"
#include "AliClusters.h"
#include "AliClustersArray.h" 
#include "AliTPCClustersRow.h" 

#include "AliTPCClustersArray.h"
#include "TClonesArray.h"
#include "TDirectory.h"
#include <TClass.h>



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
