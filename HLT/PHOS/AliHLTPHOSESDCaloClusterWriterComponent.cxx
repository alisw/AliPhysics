
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Oystein Djuvsland <oysteind@ift.uib.no>                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#include <iostream>

#include "AliHLTPHOSESDCaloClusterWriterComponent.h"
#include "AliESDCaloCluster.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TFile.h"

/** @file   AliHLTPHOSESDCaloClusterWriterComponent.cxx
    @author Oystein Djuvsland
    @date   
    @brief  An ESD calo cluster writer component for PHOS HLT
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
using namespace std;
#endif

AliHLTPHOSESDCaloClusterWriterComponent gAliHLTPHOSESDCaloClusterWriterComponent;

AliHLTPHOSESDCaloClusterWriterComponent::AliHLTPHOSESDCaloClusterWriterComponent(): 
  AliHLTPHOSProcessor(), 
  fOutfile(0),
  fOutfileName(0),
  fWriteModulo(1000),
  fESDCaloClusterTreePtr(0),
  fESDCaloClustersPtr(0)
{
  //See headerfile for documentation
}

AliHLTPHOSESDCaloClusterWriterComponent::~AliHLTPHOSESDCaloClusterWriterComponent()
{
  //See headerfile for documentation

  if(fESDCaloClustersPtr)
    {
      fESDCaloClustersPtr->Write();
      delete fESDCaloClustersPtr;
      fESDCaloClustersPtr = 0;
    }
  if (fOutfile)
    {
      fOutfile->Close();
      delete fOutfile;
      fOutfile = 0;
    }
}


int
AliHLTPHOSESDCaloClusterWriterComponent::Deinit()
{
  //See headerfile for documentation

  if(fESDCaloClustersPtr)
    {
      fESDCaloClustersPtr->Write();
      delete fESDCaloClustersPtr;
      fESDCaloClustersPtr = 0;
    }
  if (fOutfile)
    {
      fOutfile->Close();
      delete fOutfile;
      fOutfile = 0;
    }
    return 0;
}

const Char_t*
AliHLTPHOSESDCaloClusterWriterComponent::GetComponentID()
{
  //See headerfile for documentation

  return "PhosEsdCaloClusterWriter";
}

void
AliHLTPHOSESDCaloClusterWriterComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  //See headerfile for documentation
  list.clear();
  list.push_back(AliHLTPHOSDefinitions::fgkESDCaloClusterDataType);

}

AliHLTComponentDataType
AliHLTPHOSESDCaloClusterWriterComponent::GetOutputDataType()
{
  //See headerfile for documentation
  return kAliHLTAnyDataType;
}

void
AliHLTPHOSESDCaloClusterWriterComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  //See headerfile for documentation

  constBase = 30;
  inputMultiplier = 2;
}
 
Int_t 
AliHLTPHOSESDCaloClusterWriterComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/ ) 
{
  // see header file for class documentation
 
  const TObject* iter = 0;
  //  const AliHLTComponentBlockData* iter = 0;

  //UInt_t specification = 0;
  //  iter = GetFirstInputObject(AliHLTPHOSDefinitions::fgkESDCaloClusterDataType);
  //cout << iter << endl;
  for(iter = GetFirstInputObject(AliHLTPHOSDefinitions::fgkESDCaloClusterDataType); iter != NULL; iter = GetNextInputObject())
      {


	//fESDCaloClustersPtr->AddAll(reinterpret_cast<const TClonesArray*>(iter));
	fESDCaloClustersPtr = reinterpret_cast<TClonesArray*>(const_cast<TObject*>(iter)->Clone());
	
      }



    //  iter = GetFirstInputBlock(AliHLTPHOSDefinitions::fgkESDCaloClusterDataType);

//   if(iter)
//     {
//       specification = specification|iter->fSpecification;
//      cout << (reinterpret_cast<TClonesArray*>(iter->fPtr))->GetEntries() << endl;
//       //      fESDCaloClustersPtr->AddAll(reinterpret_cast<TClonesArray*>(iter->fPtr));
//       while((iter = GetNextInputBlock()))
// 	{
// 	  //	  fESDCaloClustersPtr->AddAll(reinterpret_cast<TClonesArray*>(iter->fPtr));
// 	}
      
//     }
//  HLTError("Number of clusters: %d", fESDCaloClustersPtr->GetEntries());
  fESDCaloClusterTreePtr->Fill();

  fPhosEventCount++;

  if(fPhosEventCount % fWriteModulo == 0) 
    {
      WriteTree();
    }

  //fESDCaloClustersPtr->Delete();

  return 0;
}

int 
AliHLTPHOSESDCaloClusterWriterComponent::WriteTree()
{

  // See headerfile for documentation
  char tmpFilename[128];
  sprintf(tmpFilename, "%s_%d.root", fOutfileName, fPhosEventCount/fWriteModulo);

  fOutfile = TFile::Open(tmpFilename, "RECREATE");
  if(fOutfile == 0)
    {
      HLTError("Could not open file %s for writing", tmpFilename);
      return -1;
    }
  fESDCaloClusterTreePtr->Write();
  fOutfile->Close();
  fESDCaloClusterTreePtr->Reset();
  return 0;
}

int
AliHLTPHOSESDCaloClusterWriterComponent::DoInit(int argc, const char** argv )
{
  //See headerfile for documentation
  fOutfileName = new char[128];
  char tmpFilename[128];
  ScanArgumentsModule(argc, argv);
  for (int i = 0; i < argc; i++)
    {
      if(!strcmp("-filename", argv[i]))
	{
	  sprintf(fOutfileName, "/tmp/%s", argv[i+1]);
	  sprintf(tmpFilename, "%s_0.root", fOutfileName);
	}
      if(!strcmp("-writemodulo", argv[i]))
	{
	  fWriteModulo = atoi(argv[i+1]);
	}
    }
  //  fESDCaloClustersPtr = new TClonesArray(AliESDCaloCluster::Class(), 0);
  fESDCaloClusterTreePtr = new TTree("caloClusterTree", "Tree containing AliESDCaloClusters reconstructed in HLT");
  fESDCaloClusterTreePtr->Branch("CaloClusters", &fESDCaloClustersPtr);

    //fESDCaloClusterTreePtr->SetBranchAddress("CaloClusters", &fESDCaloClustersPtr);

  return 0;
}

AliHLTComponent*
AliHLTPHOSESDCaloClusterWriterComponent::Spawn()
{
  //See headerfile for documentation

  return new AliHLTPHOSESDCaloClusterWriterComponent();
}
 
