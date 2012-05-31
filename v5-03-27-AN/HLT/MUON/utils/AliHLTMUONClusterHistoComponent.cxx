//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Author: Arshad Ahmad Masoodi <Arshad.Ahmad@cern.ch>            *
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

/** @file   AliHLTMUONClusterHistoComponent.cxx
    @author Arshad Ahmad <Arshad.Ahmad@cern.ch>
    @date  27 Nov 2009
    @brief  Component for onlinehistograms
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTMUONClusterHistoComponent.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliHLTDataTypes.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONClustersBlockStruct.h"
#include "AliHLTMUONDataBlockReader.h"
#include <TFile.h>
#include <TString.h>
#include "TObjString.h"
#include "TObjArray.h"


/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTMUONClusterHistoComponent);

AliHLTMUONClusterHistoComponent::AliHLTMUONClusterHistoComponent() :
  AliHLTMUONProcessor(),
  fChargePerClusterBending(NULL),
  fChargePerClusterNonBending(NULL),
  fNumberOfClusters(NULL),
  fPlotChargePerClusterBending(kTRUE),
  fPlotChargePerClusterNonBending(kTRUE),
  fPlotNClusters(kTRUE)
{
  // see header file for class documentation
}

AliHLTMUONClusterHistoComponent::~AliHLTMUONClusterHistoComponent()
{
  // see header file for class documentation
  
 if (fChargePerClusterBending != NULL)		delete fChargePerClusterBending;
 if (fChargePerClusterNonBending != NULL)	delete fChargePerClusterNonBending;
 if (fNumberOfClusters != NULL) 	 	delete fNumberOfClusters;
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTMUONClusterHistoComponent::GetComponentID()
{
  // see header file for class documentation
  
  return AliHLTMUONConstants::ClusterHistogrammerId();
//   return "MUONClusterHistogrammer";

}

void AliHLTMUONClusterHistoComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back(AliHLTMUONConstants::ClusterBlockDataType() );
  //list.push_back( AliHLTMUONConstants::RecHitsBlockDataType() );
  //list.push_back( AliHLTMUONConstants::TriggerRecordsBlockDataType() );
}

AliHLTComponentDataType AliHLTMUONClusterHistoComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeHistogram;
}


int AliHLTMUONClusterHistoComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)

{
  // see header file for class documentation
  tgtList.clear();
  tgtList.push_back( AliHLTMUONConstants::HistogramDataType() );
  return tgtList.size();
}

void AliHLTMUONClusterHistoComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  
  // The total constant size is the size of the TH1F object, the space needed for
  // the arrays of data points and finally we add 4k for any extra streamer meta data
  // etc. added by ROOT.
  constBase = sizeof(TH1F)*3 + sizeof(Double_t)*(400*2+101) + 4*1024*3;
  
  // The multiplier can be zero since the histograms generated are a constant size
  // independent of the input data size.
  inputMultiplier = 0;
}

AliHLTComponent* AliHLTMUONClusterHistoComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTMUONClusterHistoComponent;
}

int AliHLTMUONClusterHistoComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  
  fPlotChargePerClusterBending=kTRUE;
  fPlotChargePerClusterNonBending=kTRUE;
  fPlotNClusters=kTRUE;
   
  if(fPlotChargePerClusterBending)
    fChargePerClusterBending = new TH1F("fChargePerClusterBending","Total Charge of clusters ",400,0,4000);
  
  if(fPlotChargePerClusterNonBending)
    fChargePerClusterNonBending = new TH1F("fChargePerClusterNonBending","Total Charge of clusters ",400,0,4000);
  
  if(fPlotNClusters)
    fNumberOfClusters = new TH1F("fNumberOfClusters","Total Number of Clusters",101,0,100);
  
  int iResult=0;
  TString configuration="";
  TString argument="";
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (!configuration.IsNull()) configuration+=" ";
    configuration+=argument;
  }
  
  if (!configuration.IsNull()) {
    iResult=Configure(configuration.Data());
  }

  return iResult;
}

int AliHLTMUONClusterHistoComponent::DoDeinit()
{
  // see header file for class documentation
  if (fChargePerClusterBending != NULL)
  {
    delete fChargePerClusterBending;
    fChargePerClusterBending = NULL;
  }
  if (fChargePerClusterNonBending != NULL)
  {
    delete fChargePerClusterNonBending;
    fChargePerClusterNonBending = NULL;
  }
  if (fNumberOfClusters != NULL)
  {
    delete fNumberOfClusters;
    fNumberOfClusters = NULL;
  }
  return 0;
}

int AliHLTMUONClusterHistoComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/,
                                                const AliHLTComponentBlockData* /*blocks*/,
                                                AliHLTComponentTriggerData& /*trigData*/,
                                                AliHLTUInt8_t* /*outputPtr*/,
                                                AliHLTUInt32_t& /*size*/,
                                                AliHLTComponentBlockDataList& /*outputBlocks*/)
{
  AliHLTUInt32_t specification = 0x0;
  
  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) ) return 0;
     
  for (const AliHLTComponentBlockData* block = GetFirstInputBlock(AliHLTMUONConstants::ClusterBlockDataType());
       block != NULL;block = GetNextInputBlock()){
    
    if (block->fDataType != AliHLTMUONConstants::ClusterBlockDataType()) continue;
    specification |= block->fSpecification;  // The specification bit pattern should indicate all the DDLs that contributed.
    
    HLTDebug("Handling block with fDataType = '%s', fPtr = %p and fSize = %u bytes.",
	     DataType2Text(block->fDataType).c_str(), block->fPtr, block->fSize
	     );
    
    AliHLTMUONClustersBlockReader clusterBlock(block->fPtr, block->fSize);
    if (not BlockStructureOk(clusterBlock))
    {
      //FIXME: Need to inherit AliHLTMUONProcessor DoInit functionality properly for
      // the following Dump feature to work properly. Like in AliHLTMUONRawDataHistoComponent.
      //if (DumpDataOnError()) DumpEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks);
      continue;
    }
    
    const AliHLTMUONClusterStruct *cluster=clusterBlock.GetArray();
    for (AliHLTUInt32_t i = 0; i < clusterBlock.Nentries(); ++i){
      
      //commented for future use
      //AliHLTInt32_t detelement=cluster->fDetElemId;  // Detector ID number from AliRoot geometry

      //AliHLTUInt16_t nchannelsB=cluster->fNchannelsB; // Number of channels/pads in the cluster in bending plane.
      //AliHLTUInt16_t nchannelsNB=cluster->fNchannelsNB; // Number of channels/pads in the cluster in non-bending plane.
      

      AliHLTFloat32_t BCharge= cluster->fChargeB; // Cluster charge in bending plane. Can be -1 if invalid or uncomputed.
      AliHLTFloat32_t NBCharge= cluster->fChargeNB; // Cluster charge in bending plane. Can be -1 if invalid or uncomputed.
      
      if(fPlotChargePerClusterBending)
	fChargePerClusterBending->Fill(BCharge);
      if(fPlotChargePerClusterNonBending)	
	fChargePerClusterNonBending->Fill(NBCharge);
      cluster++;
    }  // forloop clusterBlock.Nentries
    if(fPlotNClusters and clusterBlock.Nentries()>0)
	    fNumberOfClusters->Fill(clusterBlock.Nentries());
    
  }//forloop clusterblock
  
  if( fPlotChargePerClusterBending ) PushBack( fChargePerClusterBending, AliHLTMUONConstants::HistogramDataType(), specification);
  if( fPlotChargePerClusterNonBending ) PushBack( fChargePerClusterBending, AliHLTMUONConstants::HistogramDataType(), specification);
  if( fPlotNClusters ) PushBack( fNumberOfClusters, AliHLTMUONConstants::HistogramDataType(), specification);

  return 0;
}


int AliHLTMUONClusterHistoComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // see header file for class documentation
  int iResult=0;
  const char* path="HLT/ConfigMUON/MUONHistoComponent";
  const char* defaultNotify="";
  if (cdbEntry) {
    path=cdbEntry;
    defaultNotify=" (default)";
  }
  if (path) {
    HLTInfo("reconfigure from entry %s%s, chain id %s", path, defaultNotify,(chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
    AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
        if (pEntry) {
        TObjString* pString=dynamic_cast<TObjString*>(pEntry->GetObject());
        if (pString) {
        HLTInfo("received configuration object string: \'%s\'", pString->GetString().Data());
        iResult=Configure(pString->GetString().Data());
        } else {
        HLTError("configuration object \"%s\" has wrong type, required TObjString", path);
      }
        } else {
      HLTError("can not fetch object \"%s\" from CDB", path);
    }
  }

  return iResult;
}


int AliHLTMUONClusterHistoComponent::Configure(const char* arguments)
{
  
  int iResult=0;
  if (!arguments) return iResult;
  
  TString allArgs=arguments;
  TString argument;
  
  TObjArray* pTokens=allArgs.Tokenize(" ");
  
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;
      
    }
    delete pTokens;
  }
  
  return iResult;
}
