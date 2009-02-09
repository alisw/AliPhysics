/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
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

/** @file   AliHLTEMCALClusterizerComponent.cxx
    @author Mateusz Ploskon
    @date   
    @brief  EMCAL clusterizer component for HLT. */

#if __GNUC__== 3
using namespace std;
#endif

#include "AliHLTEMCALClusterizerComponent.h"
#include "AliHLTEMCALDefinitions.h"
#include "AliHLTEMCALUtils.h"

#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TTree.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliRawReaderMemory.h"

// this is a global object used for automatic component registration, do not use this
AliHLTEMCALClusterizerComponent gAliHLTEMCALClusterizerComponent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTEMCALClusterizerComponent)

AliHLTEMCALClusterizerComponent::AliHLTEMCALClusterizerComponent()
  : AliHLTProcessor()
  , fOutputPercentage(100)
  , fStorageDBpath("local://$ALICE_ROOT/OCDB")
  , fCDB(NULL)
  , fGeometryFileName("")
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTEMCALClusterizerComponent::AliHLTEMCALClusterizerComponent(const AliHLTEMCALClusterizerComponent &/*c*/)
  : AliHLTProcessor()
  , fOutputPercentage(100)
  , fStorageDBpath("local://$ALICE_ROOT/OCDB")
  , fCDB(NULL)
  , fGeometryFileName("")
{
  // may not use the copy contructor
  HLTError("May not use.");
}

AliHLTEMCALClusterizerComponent& AliHLTEMCALClusterizerComponent::operator=(const AliHLTEMCALClusterizerComponent&)
{
  // may not use the copy contructor
  HLTError("May not use.");
  return *this;
}

AliHLTEMCALClusterizerComponent::~AliHLTEMCALClusterizerComponent()
{
  // see header file for class documentation
}

AliHLTComponentDataType AliHLTEMCALClusterizerComponent::GetOutputDataType()
{
  //return AliHLTEMCALDefinitions::fgkClusterDataType | AliHLTEMCALDefinitions::fgkDigitDataType;
  return AliHLTEMCALDefinitions::fgkClusterDataType;
}

const char* AliHLTEMCALClusterizerComponent::GetComponentID()
{ 
  return "EMCALClusterizer";
}

void AliHLTEMCALClusterizerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) 
{
  list.push_back(kAliHLTAnyDataType);
}

void AliHLTEMCALClusterizerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) 
{
  constBase = 0;
  inputMultiplier = ((double)fOutputPercentage)/100.0;
}

int AliHLTEMCALClusterizerComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;
  HLTInfo("parsing %d arguments", argc);

  TString argument="";
  TString configuration=""; 
  int bMissingParam=0;
  bool bHaveMandatory1=false;
  bool bHaveMandatory2=false;

  char *cpErr = 0;

  for (int i=0; i<argc && iResult>=0; i++) 
    {
      argument=argv[i];
      if (argument.IsNull()) continue;

    // -mandatory1
    if (argument.CompareTo("-cdb")==0) 
      {
	bHaveMandatory1|=1;
	if ((bMissingParam=(++i>=argc))) break;
	HLTInfo("got \'-cdb\' argument: %s", argv[i]);
	fStorageDBpath = argv[i];
	HLTInfo("CDB path is: %s", fStorageDBpath.c_str());
	// -mandatory2
      } 
    else if (argument.CompareTo("-geometry")==0) 
      {
	bHaveMandatory2|=1;
	if ((bMissingParam=(++i>=argc))) break;
	HLTInfo("got \'-geometry\' argument");
	fGeometryFileName = argv[i];
	HLTInfo("Geometry file is: %s", fGeometryFileName.c_str());
	// -optional1
      } 
    else if (argument.CompareTo("-output_percentage")==0) 
      {
	if ((bMissingParam=(++i>=argc))) break;
	HLTInfo("got \'-output_percentage\' argument: %s", argv[i]);
	fOutputPercentage = strtoul(argv[i], &cpErr, 0);
	if ( *cpErr )
	  {
	    HLTError("Unable to convert ouput_percentage to a number %s", argv[i]);
	    return -EINVAL;
	  }
      } 
    else 
      {
	// the remaining arguments are treated as configuration
	if (!configuration.IsNull()) configuration+=" ";
	configuration+=argument;
      }
    }

  if (bMissingParam) 
    {
      HLTError("missing parameter for argument %s", argument.Data());
      iResult=-EINVAL;
    }

  if (iResult>=0 && !bHaveMandatory1) 
    {
      HLTError("mandatory argument \'-cdb\' missing");
      iResult=-EPROTO;
    }

  if (iResult>=0 && !bHaveMandatory2) 
    {
      HLTError("mandatory argument \'-geometry\' missing");
      iResult=-EPROTO;
    }
  
  if (iResult>=0 && !configuration.IsNull()) 
    {
      iResult=Configure(configuration.Data());
    } 
  else 
    {
      iResult=Reconfigure(NULL, NULL);
    }
  
  // Initialize here
  // raw reader

  // geometry
  if (AliHLTEMCALUtils::GetGeometry() == NULL)
    {
      HLTError("unable to init geometry");
      iResult=-EPROTO;      
    }

  // OCDB

  // clusterizer and raw utils
  if (AliHLTEMCALUtils::GetRawUtils() == NULL)
    {
      HLTError("unable to init rawutils");
      iResult=-EPROTO;      
    }

  if (AliHLTEMCALUtils::GetRawUtils() == NULL)
    {
      HLTError("unable to init rawutils");
      iResult=-EPROTO;      
    }

  if (AliHLTEMCALUtils::GetRecParam() == NULL)
    {
      HLTError("unable to init reco params");
      iResult=-EPROTO;      
    }      

  if (AliHLTEMCALUtils::GetClusterizer() == NULL)
    {
      HLTError("unable to init clusterizer");
      iResult=-EPROTO;      
    }      

  return iResult;
}

int AliHLTEMCALClusterizerComponent::DoDeinit()
{
  // see header file for class documentation
  AliHLTEMCALUtils::Cleanup();
  HLTInfo("processing cleanup");
  return 0;
}

int AliHLTEMCALClusterizerComponent::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
				      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
				      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks ) {
  // see header file for class documentation
  // check if the input data are there at all - empty events possible
  HLTDebug("processing data");
  if (evtData.fStructSize==0 && blocks==NULL && trigData.fStructSize==0 &&
      outputPtr==0 && size==0)
  {
    outputBlocks.clear();
    // this is just to get rid of the warning "unused parameter"
    return 0;
  }

  UChar_t *buffer = 0;
  ULong_t  buffSize = 0;
  Int_t    eqID = 4608; //fixed for the moment - get this one from the input data
  TString  recoOption = "";
  AliHLTUInt32_t dBlockSpecification = 0;

  // Loop over all input blocks in the event
  for ( unsigned long i = 0; i < evtData.fBlockCnt; i++ )
    {
      //if (blocks[i].fDataType != (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginEMCAL))
      if (blocks[i].fDataType != (kAliHLTDataTypeDDLRaw)) // for the moment until we have the EMCAL in the publisher setup
	{
 	  HLTWarning("COMPARE FAILED received type=%d expected-or-type=%d",
		     blocks[i].fDataType, kAliHLTDataTypeDDLRaw | kAliHLTDataOriginEMCAL);
	  continue;
	}

      dBlockSpecification = blocks[i].fSpecification;
      unsigned long blockSize = blocks[i].fSize;
      buffSize += blockSize;
    }

  buffer = (UChar_t*)calloc(buffSize, 1);
  AliHLTUInt8_t *pBuf = (AliHLTUInt8_t *)buffer;
  if (buffer == NULL)
    {
      HLTError("Unable to allocate %lu bytes", buffSize);
      return -1;
    }

  // Make the memory continous
  unsigned long copied = 0;
  for ( unsigned long i = 0; i < evtData.fBlockCnt; i++ )
    {
      // we process only the raw data from TRD
      if ( blocks[i].fDataType != (kAliHLTDataTypeDDLRaw | kAliHLTDataOriginEMCAL) )
	continue;

      void *pos = (void*)(pBuf + copied);
      void *copyret = memcpy(pos, blocks[i].fPtr, blocks[i].fSize);
      if (copyret < 0)
	{
	  HLTError("MEMORY Unable to copy %lu bytes", blocks[i].fSize);
	  return -1;
	}
      copied += blocks[i].fSize;
    }

  // do it all function RAW->ClusterTree
  //   TTree* AliHLTEMCALUtils::RawBuffer2Clusters(UChar_t *buffer, ULong_t buffSize, 
  // 					      Int_t eqID, 
  // 					      Option_t* sDigitsOption)

  AliHLTEMCALUtils::ResetReconstructionTrees();
  TTree* outputTree = AliHLTEMCALUtils::RawBuffer2Clusters(buffer, buffSize, eqID, recoOption.Data());
  PushBack(outputTree, AliHLTEMCALDefinitions::fgkClusterDataType);

  TTree *digitsTree = AliHLTEMCALUtils::GetDigitsTree();
  PushBack(digitsTree, AliHLTEMCALDefinitions::fgkDigitDataType);  
  
  free(buffer);
  
  return 0;
}

int AliHLTEMCALClusterizerComponent::Configure(const char* arguments)
{
  // see header file for class documentation
  int iResult=0;
  if (!arguments) return iResult;
  HLTInfo("parsing configuration string \'%s\'", arguments);

  TString allArgs = arguments;
  TString argument;
  int bMissingParam=0;

  TObjArray* pTokens = allArgs.Tokenize(" ");
  if (pTokens) 
    {
      for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) 
	{
	  argument=((TObjString*)pTokens->At(i))->GetString();
	  if (argument.IsNull()) continue;
	  HLTInfo("processing argument %\n", argument.Data());
	  // -config1
	  if (argument.CompareTo("-cdb")==0) 
	    {
	      if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	      HLTInfo("got \'-cdb\': %s", ((TObjString*)pTokens->At(i))->GetString().Data());	      
	      // -config2
	    } 
	  else if (argument.CompareTo("-geometry")==0) 
	    {
	      if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	      HLTInfo("got \'-geometry\'");
	    } 
	  else if (argument.CompareTo("-output_percentage")==0) 
	    {
	      if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	      HLTInfo("got \'-output_percentage\'");
	    } 
	  else 
	    {
	      HLTError("unknown argument %s", argument.Data());
	      iResult=-EINVAL;
	      break;
	    }
	}
      delete pTokens;
    }

  if (bMissingParam) 
    {
      HLTError("missing parameter for argument %s", argument.Data());
      iResult=-EINVAL;
    }

  return iResult;
}

int AliHLTEMCALClusterizerComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // see header file for class documentation
  int iResult=0;
  const char* path="HLT/ConfigEMCAL/EMCALClusterizerComponent";
  const char* defaultNotify="";
  if (cdbEntry) 
    {
      path=cdbEntry;
      defaultNotify=" (default)";
    }
  if (path) 
    {
      HLTInfo("reconfigure from entry %s%s, chain id %s", path, defaultNotify,(chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
      AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
      if (pEntry) 
	{
	  TObjString* pString=dynamic_cast<TObjString*>(pEntry->GetObject());
	  if (pString) 
	    {
	      HLTInfo("received configuration object string: \'%s\'", pString->GetString().Data());
	      iResult=Configure(pString->GetString().Data());
	    } else 
	    {
	      HLTError("configuration object \"%s\" has wrong type, required TObjString", path);
	    }
	} 
      else 
	{
	  HLTError("can not fetch object \"%s\" from CDB", path);
	}
    }
  return iResult;
}

int AliHLTEMCALClusterizerComponent::ReadPreprocessorValues(const char* modules)
{
  // see header file for class documentation
  int iResult=0;
  TString detectors(modules!=NULL?modules:"");
  HLTInfo("read preprocessor values for detector(s): %s", detectors.IsNull()?"none":detectors.Data());
  return iResult;
}
