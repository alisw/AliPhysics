// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTTPCOfflineClustererComponent.cxx
    @author Matthias Richter & Jacek Otwinowski 
    @date   
    @brief  Wrapper component to the TPC offline cluster finder
*/

#include "AliHLTTPCOfflineClustererComponent.h"
#include "AliHLTTPCDefinitions.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliTPCRecoParam.h"
#include "AliTPCParam.h"
#include "AliTPCParamSR.h"
#include "AliRawReaderMemory.h"
#include "AliTPCclustererMI.h"
#include "AliTPCClustersRow.h"
#include "AliTracker.h"
#include "AliDAQ.h"
#include "TString.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TObjString.h"
#include "TTree.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCOfflineClustererComponent)

AliHLTTPCOfflineClustererComponent::AliHLTTPCOfflineClustererComponent() 
: AliHLTProcessor(),
fOutputPercentage(100),
fGeometryFileName(""),
fTPCRecoParam(0),
fTPCGeomParam(0),
fRawReader(0),
fClusterer(0)
{
  // Default constructor
  fGeometryFileName = getenv("ALICE_ROOT");
  fGeometryFileName += "/HLT/TPCLib/offline/geometry.root";
}

AliHLTTPCOfflineClustererComponent::~AliHLTTPCOfflineClustererComponent()
{
  // Destructor
}

const char* AliHLTTPCOfflineClustererComponent::GetComponentID()
{
  // Return component ID
  return "TPCOfflineClusterer";
}

void AliHLTTPCOfflineClustererComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // Get the list of input data types
  list.push_back(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTPC);
}

AliHLTComponentDataType AliHLTTPCOfflineClustererComponent::GetOutputDataType()
{
  // Return output data type
  // TObjArray/TClonesArray of clusters
  return  kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC/*AliHLTTPCDefinitions::fgkOfflineClustersDataType*/;

}

void AliHLTTPCOfflineClustererComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // Get output data size
  constBase = 20000;
  inputMultiplier = ((double)fOutputPercentage)/100.0;
}

AliHLTComponent* AliHLTTPCOfflineClustererComponent::Spawn()
{
  // Return a new instance of the class 
  return new AliHLTTPCOfflineClustererComponent;
}

int AliHLTTPCOfflineClustererComponent::DoInit( int argc, const char** argv )
{
  // Perfom initialisation. Checks arguments (argc,argv)  
  // and make initialisation. Returns 0 if success.  
  //
#ifdef HAVE_NOT_TPCOFFLINE_REC
  HLTError("AliRoot version > v4-13-Release required");
  return -ENOSYS;
#endif //HAVE_NOT_TPCOFFLINE_REC

  int iResult=0;
  HLTInfo("parsing %d arguments", argc);

  TString argument="";
  TString configuration=""; 
  int bMissingParam=0;

  // loop over input parameters
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

    if (argument.CompareTo("-geometry")==0) {
      if ((bMissingParam=(++i>=argc))) break;

      HLTInfo("got \'-geometry\' argument: %s", argv[i]);
      fGeometryFileName = argv[i];
      HLTInfo("Geometry file is: %s", fGeometryFileName.c_str());

      // the remaining arguments are treated as configuration
    } else {
      if (!configuration.IsNull()) configuration+=" ";
      configuration+=argument;
    }
  } // end loop

  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }

  if (iResult>=0 && !configuration.IsNull()) {
    iResult=Configure(configuration.Data());
  } else {
    iResult=Reconfigure(NULL, NULL);
  }

  //
  // initialisation
  //
   
  // Load geometry
  AliGeomManager::LoadGeometry(fGeometryFileName.c_str());
  if((AliGeomManager::GetGeometry()) == 0) {
    HLTError("Cannot load geometry from file %s",fGeometryFileName.c_str());
    iResult=-EINVAL;
  }
 
  // Raw Reader
  fRawReader = new AliRawReaderMemory;

  // TPC reconstruction parameters
  fTPCRecoParam = AliTPCRecoParam::GetHLTParam();
  if (fTPCRecoParam) {
    fTPCRecoParam->SetClusterSharing(kTRUE);
  }

  // TPC geometry parameters
  fTPCGeomParam = new AliTPCParamSR;
  if (fTPCGeomParam) {
    fTPCGeomParam->ReadGeoMatrices();
  }

  // Init clusterer
  fClusterer = new AliTPCclustererMI(fTPCGeomParam,fTPCRecoParam);
#ifndef HAVE_NOT_TPCOFFLINE_REC
  fClusterer->StoreInClonesArray(kTRUE); // output clusters stored in one TClonesArray
#endif

  if (!fRawReader || !fClusterer || !fTPCRecoParam || !fTPCGeomParam) {
    HLTError("failed creating internal objects");
    iResult=-ENOMEM;
  }
  return iResult;
}

int AliHLTTPCOfflineClustererComponent::DoDeinit()
{
  // Deinitialisation of the component
  if (fTPCRecoParam) delete fTPCRecoParam; fTPCRecoParam=0;
  if (fTPCGeomParam) delete fTPCGeomParam; fTPCGeomParam=0;
  if (fRawReader) delete fRawReader; fRawReader=0;
  if (fClusterer) delete fClusterer; fClusterer=0;

  return 0;
}

int AliHLTTPCOfflineClustererComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  // see header file for class documentation
  HLTInfo("processing data");

  int iResult=0;
  int commonMinSlice=-1;
  int commonMaxSlice=-1;
  int commonMinPatch=-1;
  int commonMaxPatch=-1;

  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTPC);
       pBlock!=NULL && iResult>=0;
       pBlock=GetNextInputBlock()) {
    int slice=AliHLTTPCDefinitions::GetMinSliceNr(pBlock->fSpecification);
    int patch=AliHLTTPCDefinitions::GetMinPatchNr(pBlock->fSpecification);

    if (slice!=AliHLTTPCDefinitions::GetMaxSliceNr(pBlock->fSpecification) ||
	    patch!=AliHLTTPCDefinitions::GetMaxPatchNr(pBlock->fSpecification)) {
          HLTError("ambiguous readout partition (specification 0x%08x), skipping input block", pBlock->fSpecification);
          break;
    }
    if (slice<0 || slice>35 || patch<0 || patch>5) {
      HLTError("invalid readout partition %d/%d (specification 0x%08x, skipping input block", slice, patch,  pBlock->fSpecification);
      break;
    }

    if (commonMinPatch<0) {
      commonMinSlice=slice;
      commonMaxSlice=slice;
      commonMinPatch=patch;
      commonMaxPatch=patch;
    } else {
      if (commonMinSlice<slice) commonMinSlice=slice;
      if (commonMaxSlice>slice) commonMinSlice=slice;
      if (commonMinPatch<patch) commonMinPatch=patch;
      if (commonMaxPatch>patch) commonMinPatch=patch;
    }

    if (fRawReader && fClusterer) {
      // setup raw reader and cluster finder
      int ddlId=AliDAQ::DdlIDOffset("TPC");
      if (patch<2) {
	ddlId+=2*slice+patch;
      } else {
	ddlId+=72;
	ddlId+=4*slice+patch-2;	  
      }

#ifdef HAVE_NOT_ALIRAWREADERMEMORY_ADDBUFFER
      // AliRawReaderMemory does not support mulriple blocks, 
      // process on block by block basis
      fRawReader->SetMemory( reinterpret_cast<UChar_t*>( pBlock->fPtr ), pBlock->fSize );
      fRawReader->SetEquipmentID(ddlId);
      iResult=RunClusterer(pBlock->fSpecification);
#else //!HAVE_NOT_ALIRAWREADERMEMORY_ADDBUFFER
      // add all raw data blocks to one clusterer
      fRawReader->AddBuffer( reinterpret_cast<UChar_t*>( pBlock->fPtr ), pBlock->fSize, ddlId );
#endif //HAVE_NOT_ALIRAWREADERMEMORY_ADDBUFFER
    }
  }

#ifndef HAVE_NOT_ALIRAWREADERMEMORY_ADDBUFFER
  iResult=RunClusterer(AliHLTTPCDefinitions::EncodeDataSpecification(commonMinSlice, 
								     commonMaxSlice,
								     commonMinPatch,
								     commonMaxPatch));
  fRawReader->ClearBuffers();
#endif //HAVE_NOT_ALIRAWREADERMEMORY_ADDBUFFER

  return iResult;
}

int AliHLTTPCOfflineClustererComponent::RunClusterer(AliHLTUInt32_t outspec)
{
  // see header file for class documentation
  int iResult=0;
  {
    if (fRawReader && fClusterer) {
      // run the cluster finder
      fClusterer->Digits2Clusters(fRawReader);

      Int_t nbClusters = 0;
      TClonesArray* outClrow=NULL;
#ifndef HAVE_NOT_TPCOFFLINE_REC
      outClrow=fClusterer->GetOutputClonesArray();

#endif //HAVE_NOT_TPCOFFLINE_REC
     
      if (outClrow) {
      nbClusters = outClrow->GetEntriesFast() ;

      // insert TClonesArray of clusters into output stream
      PushBack(outClrow, kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC/*AliHLTTPCDefinitions::fgkOfflineClustersDataType*/, outspec);

      // clear array 
      outClrow->Delete();
      }
      HLTInfo("processing done: Number of clusters %d (specification 0x%08x)", nbClusters, outspec);

    } else {
      HLTError("component not initialized");
      iResult=-EFAULT;
    }
  }
  return iResult;
}

int AliHLTTPCOfflineClustererComponent::Configure(const char* arguments)
{
  // see header file for class documentation
  int iResult=0;
  if (!arguments) return iResult;

  TString allArgs=arguments;
  TString argument;
  int bMissingParam=0;

  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;

      if (argument.CompareTo("-something")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;

      } else {
	HLTError("unknown argument %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    }
    delete pTokens;
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTTPCOfflineClustererComponent::Reconfigure(const char* /*cdbEntry*/, const char* /*chainId*/)
{
  // see header file for class documentation
  int iResult=0;
  // CDB stuff needs to be implemented
  return iResult;
}
