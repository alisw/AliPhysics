//-*- Mode: C++ -*-
// $Id$
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Chiara Oppedisano <Chiara.Oppedisano@to.infn.it>      *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file    AliHLTZDCESDRecoComponent.cxx
    @author  Chiara Oppedisano <Chiara.Oppedisano@to.infn.it>
    @brief   ZDC reconstruction component
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTErrorGuard.h"
#include "AliHLTSystem.h"
#include "AliHLTZDCESDRecoComponent.h"
#include "AliHLTDefinitions.h"
#include "AliRawReaderMemory.h"
#include "AliGeomManager.h"
#include "AliGRPObject.h"
#include "AliZDCReconstructor.h"
#include "AliZDCRecoParam.h"
#include "AliZDCRecoParampp.h"
#include "AliZDCRecoParamPbPb.h"
#include "AliESDEvent.h"
#include "AliESDZDC.h"
#include "AliCDBManager.h"
#include "AliHLTDataTypes.h"
#include <cstdlib>
#include <cerrno>

ClassImp(AliHLTZDCESDRecoComponent)
    
//_____________________________________________________________________
AliHLTZDCESDRecoComponent::AliHLTZDCESDRecoComponent()
  : AliHLTProcessor(),    
    fRawReader(NULL),
    fReconstructor(NULL)
{
}

//_____________________________________________________________________
AliHLTZDCESDRecoComponent::~AliHLTZDCESDRecoComponent()
{
    // see header file for class documentation
}

//_____________________________________________________________________
const char* AliHLTZDCESDRecoComponent::GetComponentID()
{
    // see header file for class documentation
    return "ZDCESDReco"; // The ID of this component
}

//_____________________________________________________________________
void AliHLTZDCESDRecoComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
    // see header file for class documentation
      list.push_back(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginZDC);
}

//_____________________________________________________________________
AliHLTComponentDataType AliHLTZDCESDRecoComponent::GetOutputDataType()
{
      return kAliHLTDataTypeESDContent|kAliHLTDataOriginZDC;
}

//_____________________________________________________________________
void AliHLTZDCESDRecoComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
    // see header file for class documentation
    constBase = 1000;
    inputMultiplier = 0;
}

//_____________________________________________________________________
AliHLTComponent* AliHLTZDCESDRecoComponent::Spawn()
{
    // see header file for class documentation
    return new AliHLTZDCESDRecoComponent;
}

//_____________________________________________________________________
int AliHLTZDCESDRecoComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;

  // -- Load GeomManager
  if(AliGeomManager::GetGeometry()==NULL){
    AliGeomManager::LoadGeometry();
  }
  
  // read configuration object : HLT/ConfigZDC/
  TString cdbPath="HLT/ConfigZDC/";
  cdbPath+=GetComponentID();
  iResult=ConfigureFromCDBTObjString(cdbPath);

  // init stage 3: read the component arguments
  if (iResult>=0) {
    iResult=ConfigureFromArgumentString(argc, argv);
  }

  // init stage 1: default values for all data members
  
  TObject* pOCDBEntry = LoadAndExtractOCDBObject("GRP/GRP/Data");
  AliGRPObject* pGRP = pOCDBEntry?dynamic_cast<AliGRPObject*>(pOCDBEntry):NULL;
  Float_t  beamEnergy=0.;
  if(pGRP) beamEnergy = pGRP->GetBeamEnergy(); 

  TString beamType="";
  if(pGRP) beamType = pGRP->GetBeamType();

  //!!! set the beam type and the energy explicitly

  beamType="A-A"; 
  beamEnergy = 2760.;

  // implement the component initialization
  do {
    if (iResult<0) break;

    fReconstructor = new AliZDCReconstructor;
    if (!fReconstructor) {
      iResult=-ENOMEM;
      break;
    }

    fRawReader = new AliRawReaderMemory;
    if (!fRawReader) {
      iResult=-ENOMEM;
      break;
    }

    // implement further initialization
  } while (0);

  if (iResult<0) {
    // implement cleanup

    if (fRawReader) delete fRawReader;
    fRawReader = NULL;
    
    if (fReconstructor) delete fReconstructor;
    fReconstructor = NULL;
  }

  if (iResult>=0) {

   if((beamType.CompareTo("p-p"))==0) fReconstructor->SetRecoParam(AliZDCRecoParampp::GetLowFluxParam());
    else if((beamType.CompareTo("A-A"))==0) fReconstructor->SetRecoParam(AliZDCRecoParamPbPb::GetHighFluxParam(beamEnergy));
    else HLTWarning(" Beam type not known by ZDC!");
  

    fReconstructor->Init(beamType, beamEnergy);
  }

  return iResult;
}

//_____________________________________________________________________
int AliHLTZDCESDRecoComponent::ScanConfigurationArgument(int /*argc*/, const char** argv)
{
  // Scan configuration arguments
  // Return the number of processed arguments
  //        -EPROTO if argument format error (e.g. number expected but not found)
  //
  // The AliHLTComponent base class implements a parsing loop for argument strings and
  // arrays of strings which is invoked by ConfigureFromArgumentString/ConfigureFromCDBTObjString
  // The component needs to implement ScanConfigurationArgument in order to decode the arguments.

  int i=0;
  TString argument=argv[i];

  if (argument.IsNull()) return 0;

  return 0;
}

//_____________________________________________________________________
int AliHLTZDCESDRecoComponent::DoDeinit()
{
    if(fRawReader) delete fRawReader;
    fRawReader=NULL;
    if(fReconstructor) delete fReconstructor;
    fReconstructor=NULL;
    return 0;
}

//_____________________________________________________________________
int AliHLTZDCESDRecoComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/,
				       AliHLTComponentTriggerData& /*trigData*/)
{
  // event processing function
  int iResult=0;
  
  // check if this is a data event, there are a couple of special events
  // which should be ignored for normal processing
  if (!IsDataEvent()) return 0;
  
  // get ZDC raw input data block and set up the rawreader
  const AliHLTComponentBlockData* pBlock = GetFirstInputBlock(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginZDC);
  if (!pBlock) {
    ALIHLTERRORGUARD(1, "No ZDC input block at event %d", GetEventCount());
    return 0;
  }

  // add input block to raw reader
  if (!fRawReader->SetMemory((UChar_t*) pBlock->fPtr, pBlock->fSize )){
    HLTError("Could not add buffer of data block  %s, 0x%08x to rawreader",
	     DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification);
    iResult = -1;
  }
  
  TTree *clusterTree = new TTree();
  if (!clusterTree) {
    iResult=-ENOMEM;
  }

  if (iResult >= 0) {

    // set ZDC EquipmentID
    fRawReader->SetEquipmentID(3840);

    fReconstructor->Reconstruct(fRawReader, clusterTree);

    fReconstructor->FillZDCintoESD(clusterTree, (AliESDEvent *) NULL);

    // send AliESDZDC
    PushBack(static_cast<TObject*>(fReconstructor->GetZDCESDData()), 
	     kAliHLTDataTypeESDContent|kAliHLTDataOriginZDC,0);
   
  }

  delete clusterTree;

  // clear the rawreader
  fRawReader->ClearBuffers();    
  
  return iResult;
}


//_____________________________________________________________________
int AliHLTZDCESDRecoComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // reconfigure the component from the specified CDB entry, or default CDB entry
  // function is invoked by the framework if a reconfigure command was received.
  // 
  int iResult=0;
  TString cdbPath;
  if (cdbEntry) {
    cdbPath=cdbEntry;
  } else {
    cdbPath="HLT/ConfigZDC/";
    cdbPath+=GetComponentID();
  }
  AliInfoClass(Form("reconfigure '%s' from entry %s%s", chainId, cdbPath.Data(), cdbEntry?"":" (default)"));
  iResult=ConfigureFromCDBTObjString(cdbPath);

  return iResult;
}


//_____________________________________________________________________
int AliHLTZDCESDRecoComponent::ReadPreprocessorValues(const char* /*modules*/)
{
  // see header file for class documentation
  return 0;
}
