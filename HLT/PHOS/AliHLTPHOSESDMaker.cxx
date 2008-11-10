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
 * @file   AliHLTPHOSESDMaker.cxx
 * @author Oystein Djuvsland
 * @date 
 * @brief  ESD maker for PHOS HLT 
 */

// see header file for class documentation
// or
// refer to README to build package
// or 
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSESDMaker.h"
#include "AliHLTPHOSCaloClusterContainerStruct.h"
#include "AliHLTPHOSBase.h"
#include "AliHLTPHOSCaloClusterDataStruct.h"
#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"
#include "TClonesArray.h"

ClassImp(AliHLTPHOSESDMaker);

AliHLTPHOSESDMaker::AliHLTPHOSESDMaker() : 
  AliHLTPHOSBase(),
  fNCaloClusters(0),
  fCaloClustersPtr(0),
  fESDEventPtr(0),
  fCaloClusterContainerPtr(0)
{
  //See header file for documentation
  fCaloClustersPtr = new TClonesArray("AliESDCaloCluster", 0);
}

AliHLTPHOSESDMaker::~AliHLTPHOSESDMaker()
{
  //See header file for documentation
  if(fCaloClustersPtr)
    {
      //fCaloClustersPtr->Delete();
      //fCaloClustersPtr = 0;
    }
}

Int_t
AliHLTPHOSESDMaker::FillESDCaloClusters()
{
  //See header file for documentation

//   AliESDCaloCluster *caloCluster = 0;
//   AliHLTPHOSCaloClusterDataStruct* caloClusterStruct = 0;

  for(UInt_t i = 0; i < fCaloClusterContainerPtr->fNCaloClusters; i++)
    {

//       caloCluster = (AliESDCaloCluster*)fCaloClustersPtr->At(i + fNCaloClusters);
//       caloClusterStruct = &(fCaloClusterContainerPtr->fCaloClusterArray[i]);
//       caloCluster->SetID(caloClusterStruct->fID);
//       caloCluster->SetClusterType(caloClusterStruct->fClusterType);
//       caloCluster->SetPosition((Float_t*)&caloClusterStruct->fGlobalPos[0]);
//       caloCluster->SetE(caloClusterStruct->fEnergy);
//       caloCluster->SetClusterDisp(caloClusterStruct->fDispersion);
//       caloCluster->SetClusterChi2(caloClusterStruct->fFitQuality);
//       caloCluster->SetPid((Float_t*)&caloClusterStruct->fPID[0]);
//       caloCluster->SetM20(caloClusterStruct->fM20);
//       caloCluster->SetM02(caloClusterStruct->fM02);
//       // PT   caloCluster->SetM11(caloClusterStruct->fM11);
//       caloCluster->SetNExMax(caloClusterStruct->fNExMax);
//       caloCluster->SetEmcCpvDistance(caloClusterStruct->fEmcCpvDistance);
//       caloCluster->SetDistanceToBadChannel(caloClusterStruct->fDistToBadChannel);
//       caloCluster->SetNCells(caloClusterStruct->fNCells);
//       caloCluster->SetCellsAbsId(caloClusterStruct->fCellsAbsId);
//       caloCluster->SetCellsAmplitudeFraction(caloClusterStruct->fCellsAmpFraction);
//       fNCaloClusters++;
    }

  return 0;
}

Int_t
AliHLTPHOSESDMaker::FillESDEvent()
{
  //See header file for documentation
  
//   AliESDCaloCluster *caloCluster = 0;
//   AliHLTPHOSCaloClusterDataStruct* caloClusterStruct = 0;
  for(UInt_t i = 0; i < fCaloClusterContainerPtr->fNCaloClusters; i++)
    {
//       //      caloCluster = (AliESDCaloCluster*)fCaloClustersPtr->New(i + fNCaloClusters);
//       caloCluster = (AliESDCaloCluster*)fCaloClustersPtr->New(i);
//       caloClusterStruct = &(fCaloClusterContainerPtr->fCaloClusterArray[i]);

//       caloCluster->SetID(caloClusterStruct->fID);
//       caloCluster->SetClusterType(caloClusterStruct->fClusterType);
//       caloCluster->SetPosition((Float_t*)&caloClusterStruct->fGlobalPos[0]);
//       caloCluster->SetE(caloClusterStruct->fEnergy);
//       caloCluster->SetClusterDisp(caloClusterStruct->fDispersion);
//       caloCluster->SetClusterChi2(caloClusterStruct->fFitQuality);
//       caloCluster->SetPid((Float_t*)&caloClusterStruct->fPID[0]);
//       caloCluster->SetM20(caloClusterStruct->fM20);
//       caloCluster->SetM02(caloClusterStruct->fM02);
//       //     caloCluster->SetM11(caloClusterStruct->fM11);
//       caloCluster->SetNExMax(caloClusterStruct->fNExMax);
//       caloCluster->SetEmcCpvDistance(caloClusterStruct->fEmcCpvDistance);
//       caloCluster->SetDistanceToBadChannel(caloClusterStruct->fDistToBadChannel);
//       caloCluster->SetNCells(caloClusterStruct->fNCells);
//       caloCluster->SetCellsAbsId(caloClusterStruct->fCellsAbsId);
//       caloCluster->SetCellsAmplitudeFraction(caloClusterStruct->fCellsAmpFraction);
//       fESDEventPtr->AddCaloCluster(caloCluster);
//       fNCaloClusters++;  
    }
  
  return 0;
}
Int_t
AliHLTPHOSESDMaker::FillESDEvent(AliHLTPHOSCaloClusterContainerStruct* caloClusterContainerPtr)
{
  //See header file for documentation
  
//   AliESDCaloCluster *caloCluster = 0;
//   AliHLTPHOSCaloClusterDataStruct* caloClusterStruct = 0;
  
  //cout << "ESD: # of clusters: " << caloClusterContainerPtr->fNCaloClusters << endl; 
  for(UInt_t i = 0; i < caloClusterContainerPtr->fNCaloClusters; i++)
    {
      // caloCluster = (AliESDCaloCluster*)fCaloClustersPtr->New(i + fNCaloClusters);
//       caloCluster = (AliESDCaloCluster*)fCaloClustersPtr->New(i + fNCaloClusters);
//       caloClusterStruct = &(caloClusterContainerPtr->fCaloClusterArray[i]);
//       caloCluster->SetID(caloClusterStruct->fID);
//       caloCluster->SetClusterType(caloClusterStruct->fClusterType);
//       //      caloCluster->SetPosition((Float_t*)&caloClusterStruct->fGlobalPos[0]);
//       caloCluster->SetPosition((Float_t*)caloClusterStruct->fGlobalPos);
//       caloCluster->SetE(caloClusterStruct->fEnergy);
//       // cout << "\t\t ESD: Cluster energy: " << caloClusterStruct->fEnergy << endl;
// //       cout << "\t\t ESD: Position: x = " << caloClusterStruct->fGlobalPos[0] << " - y = " << caloClusterStruct->fGlobalPos[1] << " - z = " << caloClusterStruct->fGlobalPos[2] << endl;
//       caloCluster->SetClusterDisp(caloClusterStruct->fDispersion);
//       caloCluster->SetClusterChi2(caloClusterStruct->fFitQuality);
//       caloCluster->SetPid((Float_t*)&caloClusterStruct->fPID[0]);
//       caloCluster->SetM20(caloClusterStruct->fM20);
//       caloCluster->SetM02(caloClusterStruct->fM02);
//       caloCluster->SetNExMax(caloClusterStruct->fNExMax);
//       caloCluster->SetEmcCpvDistance(caloClusterStruct->fEmcCpvDistance);
//       caloCluster->SetDistanceToBadChannel(caloClusterStruct->fDistToBadChannel);
//       caloCluster->SetNCells(caloClusterStruct->fNCells);
//       caloCluster->SetCellsAbsId(caloClusterStruct->fCellsAbsId);
//       caloCluster->SetCellsAmplitudeFraction(caloClusterStruct->fCellsAmpFraction);
//       fESDEventPtr->AddCaloCluster(caloCluster);
      
//       fNCaloClusters++;  
    }
  
  return 0;
}

void 
AliHLTPHOSESDMaker::ResetESD()
{
  fNCaloClusters = 0;
  // fCaloClustersPtr->Delete();
}

