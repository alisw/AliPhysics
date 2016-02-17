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

//-------------------------------------------------------------------------
//  Implementation of the ITS Upgrade tracker mother class.              
//-------------------------------------------------------------------------
#include <TTree.h>
#include <Riostream.h> 

#include "AliITStrackerU.h"
#include "AliESDEvent.h"


ClassImp(AliITStrackerU)
//_________________________________________________________________________
AliITStrackerU::AliITStrackerU():
AliTracker()
{
  //
  // Default constructor
  //

 Info("AliITStrackerU","Default Costructor");
  
}

AliITStrackerU::~AliITStrackerU()
{
 //
 // Default destructor
 //  
}
//_________________________________________________________________________
Int_t AliITStrackerU::Clusters2Tracks(AliESDEvent */*event*/)
{
  //
  // To be implemented 
  //
  
  Info("Clusters2Tracks","To be implemented");
  return 0;
}
//_________________________________________________________________________
Int_t AliITStrackerU::PropagateBack(AliESDEvent * /*event*/)
{
  //
  // To be implemented 
  //
  
 Info("PropagateBack","To be implemented");
  return 0;
}
//_________________________________________________________________________
Int_t AliITStrackerU::RefitInward(AliESDEvent * /*event*/)
{
  //
  // To be implemented 
  //
  
  Info("RefitInward","To be implemented");
  return 0;
}
//_________________________________________________________________________
Int_t AliITStrackerU::LoadClusters(TTree * /*treeR*/)
{
  //
  // To be implemented 
  //
  
  Info("LoadClusters","To be implemented");
  return 0;
} 
//_________________________________________________________________________
void AliITStrackerU::UnloadClusters()
{
  //
  // To be implemented 
  //
  
  Info("UnloadClusters","To be implemented");
} 
//_________________________________________________________________________
AliCluster * AliITStrackerU::GetCluster(Int_t /*index*/) const
{
  //
  // To be implemented 
  //
  
  Info("GetCluster","To be implemented");
  return 0x0;
} 
