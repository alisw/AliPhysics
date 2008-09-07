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
 * @file   AliHLTPHOSHistogramProducer.cxx
 * @author Oystein Djuvsland
 * @date 
 * @brief  Histogram producer for PHOS HLT 
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSHistogramProducer.h"
#include "AliHLTPHOSBase.h"
#include "AliHLTPHOSCaloClusterContainerStruct.h"
#include "AliHLTPHOSCaloClusterDataStruct.h"

#include "TH1D.h"
#include "TNtuple.h"

AliHLTPHOSHistogramProducer::AliHLTPHOSHistogramProducer():
  AliHLTPHOSBase(),
  fClusterEnergiesHistPtr(0),
  fMultiplicitiesHistPtr(0),
  fClusterNtuplePtr(0),
  fFillClusterEnergies(false),
  fFillMultiplicities(false),
  fFillNtuple(false),
  fMaxNtupleEntries(1000000000)
{
  //comment
}

AliHLTPHOSHistogramProducer::~AliHLTPHOSHistogramProducer()
{
  //comment
  if(fClusterEnergiesHistPtr)
    {
      delete fClusterEnergiesHistPtr;
      fClusterEnergiesHistPtr = 0;
    }
  if(fMultiplicitiesHistPtr)
    {
      delete fMultiplicitiesHistPtr;
      fMultiplicitiesHistPtr = 0;
    }
  if(fClusterNtuplePtr)
    {
      delete fClusterNtuplePtr;
      fClusterNtuplePtr = 0;
    }
}

Int_t 
AliHLTPHOSHistogramProducer::Fill(AliHLTPHOSCaloClusterContainerStruct* clusterContainerPtr)
{
  //comment
  AliHLTPHOSCaloClusterDataStruct* tmpClusterPtr = 0;

  for(Int_t i = 0; i < clusterContainerPtr->fNCaloClusters; i++)
    {
      tmpClusterPtr = &(clusterContainerPtr->fCaloClusterArray[i]);
      if(fFillClusterEnergies)
	{
	  fClusterEnergiesHistPtr->Fill(tmpClusterPtr->fEnergy);
	}
      if(fFillMultiplicities)
	{
	  fMultiplicitiesHistPtr->Fill(tmpClusterPtr->fNCells);
	}
      if(fFillNtuple)
	{
	  if(fClusterNtuplePtr->GetEntries() > fMaxNtupleEntries)
	    {
	      delete fClusterNtuplePtr;
	      fClusterNtuplePtr = new TNtuple("cluster_ntuple", "Cluster Ntuple", "energy:multiplicity:x:y:z");
	    }
	  fClusterNtuplePtr->Fill(tmpClusterPtr->fEnergy, tmpClusterPtr->fNCells, tmpClusterPtr->fGlobalPos[0], tmpClusterPtr->fGlobalPos[1], tmpClusterPtr->fGlobalPos[2]); 
	}
    }
  return 0;
}

Int_t
AliHLTPHOSHistogramProducer::InitializeObjects()
{
  //comment
  if(fFillClusterEnergies)
    {
      fClusterEnergiesHistPtr = new TH1D("energy_dist", "Energy spectrum - all clusters", 2000, 0, 99);
    }
  if(fFillMultiplicities)
    {
      fMultiplicitiesHistPtr = new TH1D("multiplicity_dist", "Multiplicities - all clusters", 100, 0, 99);
    }
  if(fFillNtuple)
    {
      fClusterNtuplePtr = new TNtuple("cluster_ntuple", "Cluster Ntuple", "energy:multiplicity:x:y:z");
    }
  return 0;
}
