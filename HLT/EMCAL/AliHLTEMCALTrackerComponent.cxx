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

/** @file   AliHLTEMCALTrackerComponent.cxx
    @author Mateusz Ploskon
    @date   
    @brief  EMCAL tracker component for HLT. */

#if __GNUC__== 3
using namespace std;
#endif

#include "AliHLTEMCALTrackerComponent.h"
#include "AliHLTEMCALDefinitions.h"
#include "AliHLTEMCALUtils.h"

#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TTree.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliRawReaderMemory.h"
#include "AliESDEvent.h"

// this is a global object used for automatic component registration, do not use this
AliHLTEMCALTrackerComponent gAliHLTEMCALTrackerComponent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTEMCALTrackerComponent)

AliHLTEMCALTrackerComponent::AliHLTEMCALTrackerComponent()
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

AliHLTEMCALTrackerComponent::AliHLTEMCALTrackerComponent(const AliHLTEMCALTrackerComponent &/*c*/)
  : AliHLTProcessor()
  , fOutputPercentage(100)
  , fStorageDBpath("local://$ALICE_ROOT/OCDB")
  , fCDB(NULL)
  , fGeometryFileName("")
{
  // may not use the copy contructor
  HLTError("May not use.");
}

AliHLTEMCALTrackerComponent& AliHLTEMCALTrackerComponent::operator=(const AliHLTEMCALTrackerComponent&)
{
  // may not use the copy contructor
  HLTError("May not use.");
  return *this;
}

AliHLTEMCALTrackerComponent::~AliHLTEMCALTrackerComponent()
{
  // see header file for class documentation
}

AliHLTComponentDataType AliHLTEMCALTrackerComponent::GetOutputDataType()
{
  //return AliHLTEMCALDefinitions::fgkClusterDataType | AliHLTEMCALDefinitions::fgkDigitDataType;
  return AliHLTEMCALDefinitions::fgkClusterDataType;
}

const char* AliHLTEMCALTrackerComponent::GetComponentID()
{ 
  return "EMCALTracker";
}

void AliHLTEMCALTrackerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) 
{
  list.push_back(kAliHLTAnyDataType);
}

void AliHLTEMCALTrackerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) 
{
  constBase = 0;
  inputMultiplier = ((double)fOutputPercentage)/100.0;
}

int AliHLTEMCALTrackerComponent::DoInit( int argc, const char** argv )
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

  // tracker and raw utils
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

  return iResult;
}

int AliHLTEMCALTrackerComponent::DoDeinit()
{
  // see header file for class documentation
  AliHLTEMCALUtils::Cleanup();
  HLTInfo("processing cleanup");
  return 0;
}

int AliHLTEMCALTrackerComponent::DoEvent( const AliHLTComponentEventData & /*evtData*/,
					  AliHLTComponentTriggerData & /*trigData*/ )
{
  //
  // see header file for class documentation
  //

  // check if the input data are there at all - empty events possible
  
  HLTDebug("HLT::TRDTracker::DoEvent", "BLOCKS", "NofBlocks %lu", GetNumberOfInputBlocks() );

  //implement a usage of the following
  //   AliHLTUInt32_t triggerDataStructSize = trigData.fStructSize;
  //   AliHLTUInt32_t triggerDataSize = trigData.fDataSize;
  //   void *triggerData = trigData.fData;
  //HLTDebug("Struct size %d Data size %d Data location 0x%x", trigData.fStructSize, trigData.fDataSize, (UInt_t*)trigData.fData);

  //   another way to check the blocks
  //   AliHLTComponentBlockData *dblock = (AliHLTComponentBlockData *)GetFirstInputBlock( AliHLTEMCALDefinitions::fgkClusterDataType );
  //   if (dblock == 0)
  //     {
  //       HLTError(Form("First Input Block not found! 0x%x", dblock));
  //       return -1;
  //     }

  // all those should be received by the component
  AliESDEvent *esd = 0; // we assume we receive this one from a global merger component
  TTree *clustersTree = 0;
  TTree *digitsTree = 0;

  Int_t ibForce = 0;
  TObject  *tobjin = 0;
  tobjin = (TObject *)GetFirstInputObject( AliHLTEMCALDefinitions::fgkClusterDataType, "TTree", ibForce);
  if (tobjin)
    {
      clustersTree = (TTree*)tobjin;
    }

  tobjin = (TObject *)GetFirstInputObject( AliHLTEMCALDefinitions::fgkDigitDataType, "TTree", ibForce);
  if (tobjin)
    {
      digitsTree = (TTree*)tobjin;
    }

  // the data type used here is a prototype
  // esd eventually should come from TPC alone or TPC+TRD matched tracks
  // check the AliHLTEMCALDefinitions for extra comments
  tobjin = (TObject *)GetFirstInputObject( AliHLTEMCALDefinitions::fgkESDDataType, "TTree", ibForce);
  if (tobjin)
    {
      esd = (AliESDEvent*)tobjin;
    }

  if (digitsTree != 0 && clustersTree != 0 && esd != 0)
    {
      Bool_t retValue = AliHLTEMCALUtils::FillESD(digitsTree, clustersTree, esd);
      
      if (retValue == kTRUE)
	{
	  PushBack(esd, AliHLTEMCALDefinitions::fgkEMCALESDDataType);
	}
      else
	{
	  HLTWarning("Fill ESD failed.");
	}
    }
  
  if (clustersTree)
    delete clustersTree;

  if (digitsTree)
    delete digitsTree;
  
  if (esd)
    delete esd;

  return 0;
}

int AliHLTEMCALTrackerComponent::Configure(const char* arguments)
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

int AliHLTEMCALTrackerComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // see header file for class documentation
  int iResult=0;
  const char* path="HLT/ConfigEMCAL/EMCALTrackerComponent";
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

int AliHLTEMCALTrackerComponent::ReadPreprocessorValues(const char* modules)
{
  // see header file for class documentation
  int iResult=0;
  TString detectors(modules!=NULL?modules:"");
  HLTInfo("read preprocessor values for detector(s): %s", detectors.IsNull()?"none":detectors.Data());
  return iResult;
}
