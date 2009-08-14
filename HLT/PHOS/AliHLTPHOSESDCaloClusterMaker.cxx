
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
 * @file   AliHLTPHOSESDCaloClusterMaker.cxx
 * @author Oystein Djuvsland
 * @date 
 * @brief  ESD Calo Cluster maker for PHOS HLT 
 */

// see header file for class documentation
// or
// refer to README to build package
// or 
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSESDCaloClusterMaker.h"
#include "AliHLTPHOSCaloClusterHeaderStruct.h"
#include "AliHLTPHOSBase.h"
#include "AliHLTPHOSCaloClusterDataStruct.h"
#include "AliESDCaloCluster.h"
#include "TClonesArray.h"
#include "AliHLTPHOSCaloClusterReader.h"
#include "TH1F.h"
#include "TFile.h"
#include "TNtuple.h"

ClassImp(AliHLTPHOSESDCaloClusterMaker);



AliHLTPHOSESDCaloClusterMaker::AliHLTPHOSESDCaloClusterMaker() : 
  AliHLTPHOSBase(),
  fNCaloClusters(0),
  fClusterReaderPtr(0),
  fESDCaloClusterPtr(0)
{
  //See header file for documentation
  fClusterReaderPtr = new AliHLTPHOSCaloClusterReader();
  fESDCaloClusterPtr = new AliESDCaloCluster();

}

AliHLTPHOSESDCaloClusterMaker::~AliHLTPHOSESDCaloClusterMaker()
{
  //See header file for documentation

}

Int_t
AliHLTPHOSESDCaloClusterMaker::FillESDCaloClusters(TClonesArray* esdClustersPtr, AliHLTPHOSCaloClusterHeaderStruct* caloClusterHeaderPtr)
{
  //See header file for documentation

  //  AliESDCaloCluster* caloCluster = 0;
  AliHLTPHOSCaloClusterDataStruct* caloClusterStructPtr = 0;
  fClusterReaderPtr->SetMemory(caloClusterHeaderPtr);

  //  fESDCaloClusterPtr = new AliESDCaloCluster();
  //  fNCaloClusters = 0;
  fNCaloClusters = esdClustersPtr->GetEntries();
  int count = 0;
  while((caloClusterStructPtr = fClusterReaderPtr->NextCluster()) != 0)
    {
      new((*esdClustersPtr)[fNCaloClusters]) AliESDCaloCluster();
      fESDCaloClusterPtr = static_cast<AliESDCaloCluster*>((*esdClustersPtr)[fNCaloClusters]);

      fESDCaloClusterPtr->SetID(caloClusterStructPtr->fID);
      fESDCaloClusterPtr->SetClusterType(caloClusterStructPtr->fClusterType);
      fESDCaloClusterPtr->SetPosition((Float_t*)&caloClusterStructPtr->fGlobalPos[0]);
      fESDCaloClusterPtr->SetE(caloClusterStructPtr->fEnergy);
      fESDCaloClusterPtr->SetTOF(caloClusterStructPtr->fTOF);
      fESDCaloClusterPtr->SetClusterDisp(caloClusterStructPtr->fDispersion);
      fESDCaloClusterPtr->SetClusterChi2(caloClusterStructPtr->fFitQuality);
      fESDCaloClusterPtr->SetPid((Float_t*)&caloClusterStructPtr->fPID[0]);
      fESDCaloClusterPtr->SetM20(caloClusterStructPtr->fM20);
      fESDCaloClusterPtr->SetM02(caloClusterStructPtr->fM02);
      fESDCaloClusterPtr->SetNExMax(caloClusterStructPtr->fNExMax);
      fESDCaloClusterPtr->SetEmcCpvDistance(caloClusterStructPtr->fEmcCpvDistance);
      fESDCaloClusterPtr->SetDistanceToBadChannel(caloClusterStructPtr->fDistToBadChannel);
      fESDCaloClusterPtr->SetNCells(caloClusterStructPtr->fNCells);
      //       fESDCaloClusterPtr->SetCellsAbsId(caloClusterStructPtr->fCellsAbsId);
      //       fESDCaloClusterPtr->SetCellsAmplitudeFraction(caloClusterStructPtr->fCellsAmpFraction);
      fNCaloClusters++;
      count++;
    }
  //  delete fESDCaloClusterPtr;
  return count;

}
