
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Øystein Djuvsland <oysteind@ift.uib.no>                        *
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
#include <cmath>
#include "../../PHOS/AliPHOSGeometry.h"
#include <iostream>

ClassImp(AliHLTPHOSPhysicsAnalyzer);

AliHLTPHOSPhysicsAnalyzer::AliHLTPHOSPhysicsAnalyzer():fClustersPtr(NULL)
						    
						       
{

  AliPHOSGeometry *geom=AliPHOSGeometry::GetInstance("noCPV");

  fPHOSRadius = geom->GetIPtoCrystalSurface();
  
  for(Int_t i = 0; i < N_MODULES; i++)
    {
      fRotParametersCos[i] = cos((geom->GetPHOSAngle(i+1))*2*TMath::Pi()/360);

      fRotParametersSin[i] = sin((geom->GetPHOSAngle(i+1))*2*TMath::Pi()/360);
    }
}

AliHLTPHOSPhysicsAnalyzer::AliHLTPHOSPhysicsAnalyzer(const AliHLTPHOSPhysicsAnalyzer &):fClustersPtr(NULL)

{
  AliPHOSGeometry *geom=AliPHOSGeometry::GetInstance("noCPV");

  fPHOSRadius = geom->GetIPtoCrystalSurface();
  
  for(Int_t i = 0; i < N_MODULES; i++)
    {
      fRotParametersCos[i] = cos((geom->GetPHOSAngle(i+1))*2*TMath::Pi()/360);
      
      fRotParametersSin[i] = sin((geom->GetPHOSAngle(i+1))*2*TMath::Pi()/360);
    }
  cout << "Copy constructor not tested!\n";
}


AliHLTPHOSPhysicsAnalyzer::~AliHLTPHOSPhysicsAnalyzer()
{

}

void
AliHLTPHOSPhysicsAnalyzer::LocalPosition(AliHLTPHOSClusterDataStruct* clusterPtr, Float_t* locPositionPtr)
{

  locPositionPtr[0] = clusterPtr->fLocalPositionPtr[0];
  locPositionPtr[1] = clusterPtr->fLocalPositionPtr[1];

}

void
AliHLTPHOSPhysicsAnalyzer::GlobalPosition(AliHLTPHOSClusterDataStruct* clusterPtr, Float_t* positionPtr)
{
  
  Float_t tempPosX = 0;

  Int_t module = clusterPtr->fPHOSModule;

  tempPosX = CRYSTAL_SIZE*(clusterPtr->fLocalPositionPtr[0]-N_COLUMNS_MOD/2) + CRYSTAL_SIZE/2;

  positionPtr[0] = tempPosX*fRotParametersSin[module] + fPHOSRadius*fRotParametersCos[module];

  positionPtr[1] = tempPosX*fRotParametersCos[module] - fPHOSRadius*fRotParametersSin[module];

  positionPtr[2] = CRYSTAL_SIZE*(clusterPtr->fLocalPositionPtr[1]-N_ROWS_MOD/2) + CRYSTAL_SIZE/2;

}

void
AliHLTPHOSPhysicsAnalyzer::GlobalPosition(Float_t* locPositionPtr, Float_t* positionPtr, Int_t module)
{ 
  
  positionPtr[0] = CRYSTAL_SIZE*(locPositionPtr[0]-N_COLUMNS_MOD/2)*fRotParametersCos[module-1] + fPHOSRadius*fRotParametersSin[module-1];

  positionPtr[1] = CRYSTAL_SIZE*(locPositionPtr[0]-N_COLUMNS_MOD/2)*fRotParametersSin[module-1] - fPHOSRadius*fRotParametersCos[module-1];
  
  positionPtr[2] = CRYSTAL_SIZE*(locPositionPtr[1]-N_ROWS_MOD);

}

void
AliHLTPHOSPhysicsAnalyzer::WriteHistogram(Char_t* fileName)
{
  TFile *outfile = new TFile(fileName,"recreate");  
  
  fRootHistPtr->Write();
  
  outfile->Close();

}


