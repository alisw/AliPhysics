// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Kalliopi Kanaki<Kalliopi.Kanaki@ift.uib.no>           *
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

/** @file   AliHLTTPCOfflineCluster.cxx
    @author Kalliopi Kanaki
    @date   
    @brief  Cluster converter from offline to HLT format and back
*/

#include "AliHLTTPCOfflineCluster.h"
#include "AliHLTTPCTransform.h"
#include "TObjArray.h"
#include <sys/time.h>

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTPCOfflineCluster)

AliHLTTPCOfflineCluster::AliHLTTPCOfflineCluster()
{
//constructor  
}

AliHLTTPCOfflineCluster::AliHLTTPCOfflineCluster(const AliHLTTPCSpacePointData& /*hltCluster*/){
  
}

AliHLTTPCOfflineCluster::AliHLTTPCOfflineCluster(const AliTPCclusterMI& /*offCluster*/){
 
}

AliHLTTPCOfflineCluster& AliHLTTPCOfflineCluster::operator=(const AliTPCclusterMI& /*offCluster*/){
  return *this;
}

AliHLTTPCOfflineCluster::~AliHLTTPCOfflineCluster(){


  //destructor  
}

AliTPCclusterMI* AliHLTTPCOfflineCluster::ConvertHLTToOffline(AliHLTTPCSpacePointData spacePoint){
   
   AliTPCclusterMI *offCluster = new AliTPCclusterMI();
   
   offCluster->SetPad(spacePoint.fX);	          // pad   
   offCluster->SetRow((Int_t)spacePoint.fPadRow); // row 
   offCluster->SetTimeBin(spacePoint.fZ);         // time bin
   offCluster->SetQ(spacePoint.fCharge);     	  // charge
   offCluster->SetMax(spacePoint.fQMax);     	  // max Q (amplitude)
   //offCluster->SetDetector(0);	     	  // detector/slice
   //offCluster->SetType(0);		     	  // default from constructor
   //offCluster->IsUsed(0);		     	  // default from constructor
   //offCluster->SetInfo(NULL);		     	  // default from constructor
     
   return offCluster;
}


AliHLTTPCSpacePointData AliHLTTPCOfflineCluster::ConvertOfflineToHLT(AliTPCclusterMI *offCluster){

     
   AliHLTTPCSpacePointData spacePoint = { 0.,0.,0.,0,0,0.,0.,0,0,kFALSE,0 };
       
   spacePoint.fX      = offCluster->GetPad();
   spacePoint.fPadRow = offCluster->GetRow();
   spacePoint.fZ      = offCluster->GetTimeBin();
   spacePoint.fCharge = (UInt_t)offCluster->GetQ();
   spacePoint.fQMax   = (UInt_t)offCluster->GetMax();

   return spacePoint;
   
}
