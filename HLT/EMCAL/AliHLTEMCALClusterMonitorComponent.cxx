/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Freancesco Blanco                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTEMCALClusterMonitorComponent.h"
#include "AliHLTEMCALClusterMonitor.h"

#include "TFile.h"
#include "TString.h"
#include "AliESDVZERO.h"


/** 
 * @file   AliHLTEMCALClusterMonitorComponent.cxx
 * @author Francesco Blanco
 * @date   
 * @brief  A component to pusk back histograms for EMCAL HLT
 */

//FIXME
AliHLTEMCALClusterMonitorComponent gAliHLTEMCALClusterMonitorComponent;

AliHLTEMCALClusterMonitorComponent::AliHLTEMCALClusterMonitorComponent() :
				  AliHLTCaloProcessor(),
				  fRootFileName("histofile_local.root"),
				  fPushFraction(10),
				  fLocalEventCount(0),
				  fBeVerbose(0),
				  fHistoMakerPtr(0)
{
	//see header file for documentation
}


AliHLTEMCALClusterMonitorComponent::~AliHLTEMCALClusterMonitorComponent()
{
	//see header file for documentation
}


int
AliHLTEMCALClusterMonitorComponent::DoInit(int argc, const char** argv )
{
	//see header file for documentation

	fHistoMakerPtr = new AliHLTEMCALClusterMonitor();
	SetupCTPData();
	for(int i = 0; i < argc; i++)
	{
		if(!strcmp("-roothistofilename", argv[i]))
			fRootFileName = argv[i+1];

		if(!strcmp("-pushfraction", argv[i]))
			fPushFraction = atoi(argv[i+1]);

		if(!strcmp("-beverbose", argv[i]))
			fBeVerbose = atoi(argv[i+1]);

	}

	if (fBeVerbose) cout << "\nI-CLUSTERMONITORCOMPONENT: local root file name is: " << fRootFileName << endl;

	return 0;
}


int 
AliHLTEMCALClusterMonitorComponent::Deinit()
{ 
	//see header file for documentation
	if(fHistoMakerPtr)
	{
		delete fHistoMakerPtr;
		fHistoMakerPtr = 0;
	}

	return 0;
}

const char*
AliHLTEMCALClusterMonitorComponent::GetComponentID()
{
	//see header file for documentation
	return "EmcalClusterMonitor";
}


void
AliHLTEMCALClusterMonitorComponent::GetInputDataTypes(vector<AliHLTComponentDataType>& list)
{ 
	//see header file for documentation
	list.clear();
	list.push_back(kAliHLTDataTypeCaloCluster);
	list.push_back(kAliHLTDataTypeESDContent | kAliHLTDataOriginVZERO);
}

AliHLTComponentDataType 
AliHLTEMCALClusterMonitorComponent::GetOutputDataType()
{
	//see header file for documentation
	return kAliHLTDataTypeHistogram | kAliHLTDataOriginEMCAL;
}

void 
AliHLTEMCALClusterMonitorComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
	//see header file for documentation
	constBase = 0;
	// to be reviewed later
	inputMultiplier = 100;
}

int 
AliHLTEMCALClusterMonitorComponent::DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
		AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* /*outputPtr*/, AliHLTUInt32_t& /*size*/,
		std::vector<AliHLTComponentBlockData>& /*outputBlocks*/)
{

	const AliHLTComponentBlockData* iter = 0;
	unsigned long ndx;

	UInt_t specification = 0;
	for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
	{
	
		AliHLTCaloClusterDataStruct* caloClusterPtr = 0;

		iter = blocks+ndx;
		
		if (fBeVerbose) PrintComponentDataTypeInfo(iter->fDataType);

		if (iter->fDataType == kAliHLTDataTypeCaloCluster) caloClusterPtr = reinterpret_cast<AliHLTCaloClusterDataStruct*>(iter->fPtr);
		else {
			if(fBeVerbose) HLTWarning("\nI-CLUSTERMONITORCOMPONENT: Data block does not contain cluster type \n");
		}

		Int_t nClusters = iter->fSize/sizeof(AliHLTCaloClusterDataStruct);
		const AliESDVZERO* esdVZERO = dynamic_cast<const AliESDVZERO*>(GetFirstInputObject(kAliHLTDataTypeESDContent | kAliHLTDataOriginVZERO, "AliESDVZERO"));

		specification |= iter->fSpecification;

		fHistoMakerPtr->MakeHisto(caloClusterPtr, nClusters, esdVZERO);

	}

	fLocalEventCount++;

	TFile rootHistFile(fRootFileName,"recreate");
	
	fHistoMakerPtr->GetHistograms()->Write();
	
	if (fLocalEventCount%fPushFraction == 0) {
	  
	  if (fBeVerbose) cout << "\nI-CLUSTERMONITORCOMPONENT: pushback done at " << fLocalEventCount << " events " << endl;
	  
	  PushHistograms(fHistoMakerPtr->GetHistograms());
	}
	
	return 0;
}

void AliHLTEMCALClusterMonitorComponent::PushHistograms(TCollection* list)
{
  HLTDebug("Collection has %d histograms", list->GetEntries());
  TIter next(list);

  TH1* histo = 0;
  while ((histo = static_cast<TH1*>(next()))) {
    if (histo->GetEntries() > 0) {
      HLTDebug("Pushing histogram %s", histo->GetName());
      Int_t nbytes = PushBack(histo, kAliHLTDataTypeHistogram | kAliHLTDataOriginEMCAL, 0);
      if (nbytes > 0)
        histo->Reset();
    } else {
      HLTDebug("Not pushing histogram %s, because it has 0 entries", histo->GetName());
    }
  }
}

AliHLTComponent*
AliHLTEMCALClusterMonitorComponent::Spawn()
{
	//see header file for documentation
	return new AliHLTEMCALClusterMonitorComponent();
}
