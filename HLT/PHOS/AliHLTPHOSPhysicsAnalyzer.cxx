
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Ãystein Djuvsland <oysteind@ift.uib.no>                        *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#include "AliHLTPHOSPhysicsAnalyzer.h"
#include "TVector3.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "AliPHOSGeometry.h"
#include "Rtypes.h"
#include "AliHLTPHOSCommonDefs.h"
#include "AliHLTPHOSClusterDataStruct.h"
 

ClassImp(AliHLTPHOSPhysicsAnalyzer);

AliHLTPHOSPhysicsAnalyzer::AliHLTPHOSPhysicsAnalyzer():fClustersPtr(0), fRootHistPtr(0), fPHOSRadius(0)
						    
						       
{
  //Constructor

  AliPHOSGeometry *geom=AliPHOSGeometry::GetInstance("noCPV");

  fPHOSRadius = geom->GetIPtoCrystalSurface();
  
  for(UInt_t i = 0; i < N_MODULES; i++)
    {
      fRotParametersCos[i] = cos((geom->GetPHOSAngle(i+1))*2*TMath::Pi()/360);

      fRotParametersSin[i] = sin((geom->GetPHOSAngle(i+1))*2*TMath::Pi()/360);
    }
}

AliHLTPHOSPhysicsAnalyzer::AliHLTPHOSPhysicsAnalyzer(const AliHLTPHOSPhysicsAnalyzer &):fClustersPtr(0), fRootHistPtr(0), fPHOSRadius(0)

{
  //Cooy constructor

  AliPHOSGeometry *geom=AliPHOSGeometry::GetInstance("noCPV");

  fPHOSRadius = geom->GetIPtoCrystalSurface();
  
  for(UInt_t i = 0; i < N_MODULES; i++)
    {
      fRotParametersCos[i] = cos((geom->GetPHOSAngle(i+1))*2*TMath::Pi()/360);
      
      fRotParametersSin[i] = sin((geom->GetPHOSAngle(i+1))*2*TMath::Pi()/360);
    }

}


AliHLTPHOSPhysicsAnalyzer::~AliHLTPHOSPhysicsAnalyzer()
{
  //Destructor

  fClustersPtr = 0;
  fRootHistPtr = 0;
}

void
AliHLTPHOSPhysicsAnalyzer::LocalPosition(AliHLTPHOSClusterDataStruct* clusterPtr, Float_t* locPositionPtr)
{
  //Get local position for a cluster

  locPositionPtr[0] = clusterPtr->fLocalPositionPtr[0];
  locPositionPtr[1] = clusterPtr->fLocalPositionPtr[1];

}

void
AliHLTPHOSPhysicsAnalyzer::GlobalPosition(AliHLTPHOSClusterDataStruct* clusterPtr, Float_t* positionPtr)
{
  //Get global position for a cluster
  
  Float_t tempPosX = 0;

  Int_t module = clusterPtr->fPHOSModule;

  tempPosX = kCRYSTAL_SIZE*(clusterPtr->fLocalPositionPtr[0]-N_XCOLUMNS_MOD/2) + kCRYSTAL_SIZE/2;

  positionPtr[0] = tempPosX*fRotParametersSin[module] + fPHOSRadius*fRotParametersCos[module];

  positionPtr[1] = tempPosX*fRotParametersCos[module] - fPHOSRadius*fRotParametersSin[module];

  positionPtr[2] = kCRYSTAL_SIZE*(clusterPtr->fLocalPositionPtr[1]-N_ZROWS_MOD/2) + kCRYSTAL_SIZE/2;

}

void
AliHLTPHOSPhysicsAnalyzer::GlobalPosition(Float_t* locPositionPtr, Float_t* positionPtr, Int_t module)
{ 
  //Get global position from local postion and module number

  positionPtr[0] = kCRYSTAL_SIZE*(locPositionPtr[0]-N_XCOLUMNS_MOD/2)*fRotParametersCos[module-1] + fPHOSRadius*fRotParametersSin[module-1];

  positionPtr[1] = kCRYSTAL_SIZE*(locPositionPtr[0]-N_XCOLUMNS_MOD/2)*fRotParametersSin[module-1] - fPHOSRadius*fRotParametersCos[module-1];
  
  positionPtr[2] = kCRYSTAL_SIZE*(locPositionPtr[1]-N_ZROWS_MOD);

}

void
AliHLTPHOSPhysicsAnalyzer::WriteHistogram(Char_t* fileName)
{
  //Write the histogram

  TFile *outfile = new TFile(fileName,"recreate");  
  
  fRootHistPtr->Write();
  
  outfile->Close();
  
  delete outfile;
  outfile = 0;

}


