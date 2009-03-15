// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** 
 * @file   AliHLTPHOSPhysicsAnalyzer.cxx
 * @author Oystein Djuvsland
 * @date 
 * @brief  Physics analysis base class  */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSPhysicsAnalyzer.h"
#include "TVector3.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "AliPHOSGeometry.h"
#include "Rtypes.h"
#include "AliHLTPHOSCommonDefs.h"
#include "AliHLTPHOSRecPointDataStruct.h"
 

ClassImp(AliHLTPHOSPhysicsAnalyzer);

AliHLTPHOSPhysicsAnalyzer::AliHLTPHOSPhysicsAnalyzer():
  fRecPointsPtr(0), 
  fRootHistPtr(0), 
  fPHOSRadius(0)
{
  //Constructor
  //See header file for documentation
  AliPHOSGeometry *geom=AliPHOSGeometry::GetInstance("noCPV");

  //  fPHOSRadius = geom->GetIPtoCrystalSurface();
  fPHOSRadius = geom->GetIPtoCrystalSurface();

  for(int i = 0; i < NMODULES; i++)
    {
//       fRotParametersCos[i] = cos((geom->GetPHOSAngle(i+1))*2*TMath::Pi()/360);
//       fRotParametersSin[i] = sin((geom->GetPHOSAngle(i+1))*2*TMath::Pi()/360);
      fRotParametersCos[i] = cos((geom->GetPHOSAngle(i))*2*TMath::Pi()/360);
      fRotParametersSin[i] = sin((geom->GetPHOSAngle(i))*2*TMath::Pi()/360);

    }
}
/*
AliHLTPHOSPhysicsAnalyzer::AliHLTPHOSPhysicsAnalyzer(const AliHLTPHOSPhysicsAnalyzer &):fClustersPtr(0), fRootHistPtr(0), fPHOSRadius(0)

{
  //Cooy constructor
  //See header file for documentation
  AliPHOSGeometry *geom=AliPHOSGeometry::GetInstance("noCPV");

  fPHOSRadius = geom->GetIPtoCrystalSurface();
  
  for(UInt_t i = 0; i < N_MODULES; i++)
    {
      fRotParametersCos[i] = cos((geom->GetPHOSAngle(i+1))*2*TMath::Pi()/360);
      
      fRotParametersSin[i] = sin((geom->GetPHOSAngle(i+1))*2*TMath::Pi()/360);
    }

}
*/

AliHLTPHOSPhysicsAnalyzer::~AliHLTPHOSPhysicsAnalyzer()
{
  //Destructor
  //See header file for documentation
  fRecPointsPtr = 0;
  fRootHistPtr = 0;
}

void
AliHLTPHOSPhysicsAnalyzer::LocalPosition(AliHLTPHOSRecPointDataStruct* /*recPointPtr*/, Float_t* /*locPositionPtr*/)
{
  //Get local position for a recPoint

  //  locPositionPtr[0] = recPointPtr->fLocalPositionPtr[0];
  //locPositionPtr[1] = recPointPtr->fLocalPositionPtr[1];

}

void
AliHLTPHOSPhysicsAnalyzer::GlobalPosition(AliHLTPHOSRecPointDataStruct* recPointPtr, Float_t* positionPtr)
{
  //Get global position for a recPoint
  //See header file for documentation
  Float_t tempPosX = 0;
 
  
  Int_t module = recPointPtr->fModule;

  tempPosX = kCRYSTALSIZE*(recPointPtr->fX-NXCOLUMNSMOD/2) + kCRYSTALSIZE/2;

  positionPtr[0] = tempPosX*fRotParametersSin[module] + fPHOSRadius*fRotParametersCos[module];

  positionPtr[1] = tempPosX*fRotParametersCos[module] - fPHOSRadius*fRotParametersSin[module];
 
  positionPtr[2] = kCRYSTALSIZE*(recPointPtr->fZ-NZROWSMOD/2) + kCRYSTALSIZE/2;

}

void
AliHLTPHOSPhysicsAnalyzer::GlobalPosition(Float_t* locPositionPtr, Float_t* positionPtr, Int_t module)
{ 
  //Get global position from local postion and module number
  //See header file for documentation
  positionPtr[0] = kCRYSTALSIZE*(locPositionPtr[0]-NXCOLUMNSMOD/2)*fRotParametersCos[module] + fPHOSRadius*fRotParametersSin[module];

  positionPtr[1] = kCRYSTALSIZE*(locPositionPtr[0]-NXCOLUMNSMOD/2)*fRotParametersSin[module] - fPHOSRadius*fRotParametersCos[module];
  
  positionPtr[2] = kCRYSTALSIZE*(locPositionPtr[1]-NZROWSMOD);

}

void
AliHLTPHOSPhysicsAnalyzer::WriteHistogram(const Char_t* fileName)
{
  //Write the histogram
  //See header file for documentation
  TFile *outfile = new TFile(fileName,"recreate");  
  
  fRootHistPtr->Write();
  
  outfile->Close();
  
  delete outfile;
  outfile = 0;

}

void
AliHLTPHOSPhysicsAnalyzer::Analyze(AliHLTPHOSRecPointContainerStruct* /*recPointsArrayPtr*/, Int_t /*nRecPoints*/)
{
  //comment
  return;
}

