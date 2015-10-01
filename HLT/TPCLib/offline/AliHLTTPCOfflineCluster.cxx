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
   
   offCluster->SetRow((Int_t)spacePoint.fPadRow); // padrow
   offCluster->SetX(spacePoint.fX);	          // X coordinate in slice system
   offCluster->SetY(spacePoint.fY);               // Y coordinate in slice system
   offCluster->SetZ(spacePoint.fZ);               // Z coordinate in slice system
   offCluster->SetQ(spacePoint.fCharge);     	  // charge
   offCluster->SetMax(spacePoint.fQMax);     	  // max Q (amplitude)
   offCluster->SetSigmaY2(spacePoint.fSigmaY2);   // error in Y direction
   offCluster->SetSigmaZ2(spacePoint.fSigmaZ2);   // error in Z direction
   //offCluster->SetDetector(0);	     	  // detector/slice
   //offCluster->SetType(0);		     	  // default from constructor
   //offCluster->IsUsed(0);		     	  // default from constructor
   //offCluster->SetInfo(NULL);		     	  // default from constructor
        
   return offCluster;
}


AliHLTTPCSpacePointData AliHLTTPCOfflineCluster::ConvertOfflineToHLT(AliTPCclusterMI *offCluster){

     
   AliHLTTPCSpacePointData spacePoint;
       
   spacePoint.fPadRow  = offCluster->GetRow();
   spacePoint.fX       = offCluster->GetX(); // these are in the detector system
   spacePoint.fY       = offCluster->GetY(); // the HLT clusters have to be transformed to the slice system
   spacePoint.fZ       = offCluster->GetZ(); // for the X, Y, Z to be consistent with our definitions
   spacePoint.fCharge  = (UInt_t)offCluster->GetQ();
   spacePoint.fQMax    = (UInt_t)offCluster->GetMax();
   spacePoint.fSigmaY2 = offCluster->GetSigmaY2();
   spacePoint.fSigmaZ2 = offCluster->GetSigmaZ2();
   
   return spacePoint;
   
}
