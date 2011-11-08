// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
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

//  @file   AliHLTComponent.cxx
//  @author Matthias Richter, Timm Steinbeck
//  @date   
//  @brief  Base class implementation for HLT components. */
//  @note   The class is both used in Online (PubSub) and Offline (AliRoot)
//          context


#if __GNUC__>= 3
using namespace std;
#endif

//#include "AliHLTStdIncludes.h"
#include "AliHLTComponent.h"
#include "AliHLTComponentHandler.h"
#include "AliHLTMessage.h"
#include "AliHLTCTPData.h"
#include "AliHLTErrorGuard.h"
#include "AliRawDataHeader.h"
#include "TString.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjectTable.h"
#include "TClass.h"
#include "TStopwatch.h"
#include "TFormula.h"
#include "TUUID.h"
#include "TMD5.h"
#include "TRandom3.h"
#include "AliHLTMemoryFile.h"
#include "AliHLTMisc.h"
#include <cassert>
#include <ctime>
#include <stdint.h>

/**
 * default compression level for ROOT objects
 */
#define ALIHLTCOMPONENT_DEFAULT_OBJECT_COMPRESSION 5
#define ALIHLTCOMPONENT_STATTIME_SCALER 1000000

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTComponent);

/** stopwatch macro using the stopwatch guard */
#define ALIHLTCOMPONENT_STOPWATCH(type) AliHLTStopwatchGuard swguard(fpStopwatches!=NULL?reinterpret_cast<TStopwatch*>(fpStopwatches->At((int)type)):NULL)
//#define ALIHLTCOMPONENT_STOPWATCH(type) 

/** stopwatch macro for operations of the base class */
#define ALIHLTCOMPONENT_BASE_STOPWATCH() ALIHLTCOMPONENT_STOPWATCH(kSWBase)
/** stopwatch macro for operations of the detector algorithm (DA) */
#define ALIHLTCOMPONENT_DA_STOPWATCH() ALIHLTCOMPONENT_STOPWATCH(kSWDA)

AliHLTComponent::AliHLTComponent()
  :
  fEnvironment(),
  fCurrentEvent(0),
  fEventCount(-1),
  fFailedEvents(0),
  fCurrentEventData(),
  fpInputBlocks(NULL),
  fCurrentInputBlock(-1),
  fSearchDataType(kAliHLTVoidDataType),
  fClassName(),
  fpInputObjects(NULL),
  fpOutputBuffer(NULL),
  fOutputBufferSize(0),
  fOutputBufferFilled(0),
  fOutputBlocks(),
  fpStopwatches(new TObjArray(kSWTypeCount)),
  fMemFiles(),
  fpRunDesc(NULL),
  fCDBSetRunNoFunc(false),
  fChainId(),
  fChainIdCrc(0),
  fpBenchmark(NULL),
  fFlags(0),
  fEventType(gkAliEventTypeUnknown),
  fComponentArgs(),
  fEventDoneData(NULL),
  fEventDoneDataSize(0),
  fCompressionLevel(ALIHLTCOMPONENT_DEFAULT_OBJECT_COMPRESSION)
  , fLastObjectSize(0)
  , fpCTPData(NULL)
  , fPushbackPeriod(0)
  , fLastPushBackTime(-1),
  fEventModulo(-1)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  memset(&fEnvironment, 0, sizeof(AliHLTAnalysisEnvironment));
  if (fgpComponentHandler)
    fgpComponentHandler->ScheduleRegister(this);
  //SetLocalLoggingLevel(kHLTLogDefault);
}

AliHLTComponent::~AliHLTComponent()
{
  // see header file for function documentation
  if (fpBenchmark) delete fpBenchmark;
  fpBenchmark=NULL;

  CleanupInputObjects();
  if (fpStopwatches!=NULL) delete fpStopwatches;
  fpStopwatches=NULL;
  AliHLTMemoryFilePList::iterator element=fMemFiles.begin();
  while (element!=fMemFiles.end()) {
    if (*element) {
      if ((*element)->IsClosed()==0) {
	HLTWarning("memory file has not been closed, possible data loss or incomplete buffer");
	// close but do not flush as we dont know whether the buffer is still valid
	(*element)->CloseMemoryFile(0);
      }
      delete *element;
      *element=NULL;
    }
    element++;
  }
  if (fpRunDesc) {
    delete fpRunDesc;
    fpRunDesc=NULL;
  }
  if (fEventDoneData)
    delete [] reinterpret_cast<AliHLTUInt8_t*>( fEventDoneData );
  fEventDoneData=NULL;

  if (fpCTPData) {
    delete fpCTPData;
  }
  fpCTPData=NULL;
}

AliHLTComponentHandler* AliHLTComponent::fgpComponentHandler=NULL;

int AliHLTComponent::SetGlobalComponentHandler(AliHLTComponentHandler* pCH, int bOverwrite) 
{
  // see header file for function documentation
  int iResult=0;
  if (fgpComponentHandler==NULL || bOverwrite!=0)
    fgpComponentHandler=pCH;
  else
    iResult=-EPERM;
  return iResult;
}

int AliHLTComponent::UnsetGlobalComponentHandler() 
{
  // see header file for function documentation
  return SetGlobalComponentHandler(NULL,1);
}

int AliHLTComponent::SetComponentEnvironment(const AliHLTAnalysisEnvironment* comenv, void* environParam)
{
  // see header file for function documentation
  HLTLogKeyword(fChainId.c_str());
  int iResult=0;
  if (comenv) {
    memset(&fEnvironment, 0, sizeof(AliHLTAnalysisEnvironment));
    memcpy(&fEnvironment, comenv, comenv->fStructSize<sizeof(AliHLTAnalysisEnvironment)?comenv->fStructSize:sizeof(AliHLTAnalysisEnvironment));
    fEnvironment.fStructSize=sizeof(AliHLTAnalysisEnvironment);
    fEnvironment.fParam=environParam;
  }
  return iResult;
}

int AliHLTComponent::Init(const AliHLTAnalysisEnvironment* comenv, void* environParam, int argc, const char** argv )
{
  // see header file for function documentation
  HLTLogKeyword(fChainId.c_str());
  int iResult=0;
  if (comenv) {
    SetComponentEnvironment(comenv, environParam);
  }
  fPushbackPeriod=0;
  fLastPushBackTime=-1;

  fComponentArgs="";
  const char** pArguments=NULL;
  int iNofChildArgs=0;
  TString argument="";
  int bMissingParam=0;
  if (argc>0) {
    pArguments=new const char*[argc];
    if (pArguments) {
      for (int i=0; i<argc && iResult>=0; i++) {
	if (fComponentArgs.size()>0) fComponentArgs+=" ";
	fComponentArgs+=argv[i];
	argument=argv[i];
	if (argument.IsNull()) continue;

	// benchmark
	if (argument.CompareTo("-benchmark")==0) {

	  // -loglevel=
	} else if (argument.BeginsWith("-loglevel=")) {
	  TString parameter=argument.ReplaceAll("-loglevel=", "");
	  parameter.Remove(TString::kLeading, ' '); // remove all blanks
	  if (parameter.BeginsWith("0x") &&
	      parameter.Replace(0,2,"",0).IsHex()) {
	    unsigned int loglevel=kHLTLogNone;
	    sscanf(parameter.Data(),"%x", &loglevel);
	    SetLocalLoggingLevel((AliHLTComponentLogSeverity)loglevel);
	  } else {
	    HLTError("wrong parameter for argument %s, hex number expected", argument.Data());
	    iResult=-EINVAL;
	  }
	  // -object-compression=
	} else if (argument.BeginsWith("-object-compression=")) {
	  argument.ReplaceAll("-object-compression=", "");
	  if (argument.IsDigit()) {
	    fCompressionLevel=argument.Atoi();
	    if (fCompressionLevel<0 || fCompressionLevel>9) {
	      HLTWarning("invalid compression level %d, setting to default %d", fCompressionLevel, ALIHLTCOMPONENT_DEFAULT_OBJECT_COMPRESSION);
	      fCompressionLevel=ALIHLTCOMPONENT_DEFAULT_OBJECT_COMPRESSION;
	    }
	  } else {
	    HLTError("wrong parameter for argument -object-compression, number expected");
	  }
	  // -pushback-period=
	} else if (argument.BeginsWith("-pushback-period=")) {
	  argument.ReplaceAll("-pushback-period=", "");
	  if (argument.IsDigit()) {
	    fPushbackPeriod=argument.Atoi();
	  } else {
	    HLTError("wrong parameter for argument -pushback-period, number expected");
	  }
	  // -event-modulo
	} else if (argument.BeginsWith("-event-modulo=")) {
	  argument.ReplaceAll("-event-modulo=", "");
	  if (argument.IsDigit()) {
	    fEventModulo=argument.Atoi();
	    if (fEventModulo < 1)
	    {
	      fEventModulo = -1;
	      HLTError("number passed in -event-modulo must be a positive integer greater or equal to 1.");
	    }
	  } else {
	    HLTError("wrong parameter for argument -event-modulo, integer number expected");
	  }
	  // -disable-component-stat
	} else if (argument.CompareTo("-disable-component-stat")==0) {
	  fFlags|=kDisableComponentStat;
	} else {
	  pArguments[iNofChildArgs++]=argv[i];
	}
      }
    } else {
      iResult=-ENOMEM;
    }
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  if (iResult>=0) {
    iResult=CheckOCDBEntries();
  }
  if (iResult>=0) {
    iResult=DoInit(iNofChildArgs, pArguments);
  }
  if (iResult>=0) {
    fEventCount=0;

    // find out if the component wants to get the steering events
    // explicitly
    AliHLTComponentDataTypeList inputDt;
    GetInputDataTypes(inputDt);
    bool bRequireSteeringBlocks=false;
    for (AliHLTComponentDataTypeList::iterator dt=inputDt.begin();
	 dt!=inputDt.end() && !bRequireSteeringBlocks;
	 dt++) {
      bRequireSteeringBlocks|=MatchExactly(*dt,kAliHLTDataTypeSOR);
      bRequireSteeringBlocks|=MatchExactly(*dt,kAliHLTDataTypeRunType);
      bRequireSteeringBlocks|=MatchExactly(*dt,kAliHLTDataTypeEOR);
      bRequireSteeringBlocks|=MatchExactly(*dt,kAliHLTDataTypeDDL);
      bRequireSteeringBlocks|=MatchExactly(*dt,kAliHLTDataTypeComponentStatistics);
    }
    if (bRequireSteeringBlocks) fFlags|=kRequireSteeringBlocks;
  }
  if (pArguments) delete [] pArguments;

#if defined(HLT_COMPONENT_STATISTICS)
  // benchmarking stopwatch for the component statistics
  fpBenchmark=new TStopwatch;
#endif // HLT_COMPONENT_STATISTICS

  return iResult;
}

int AliHLTComponent::Deinit()
{
  // see header file for function documentation
  HLTLogKeyword(fChainId.c_str());
  int iResult=0;
  iResult=DoDeinit();
  if (fpRunDesc) {
    // TODO: the warning should be kept, but the condition is wrong since the
    // AliHLTRunDesc is set before the SOR event in the SetRunDescription
    // method. A couple of state flags should be defined but that is a bit more
    // work to do. For the moment disable the warning (2009-07-01)
    // 2009-09-08: now, the info is not cleared in the ProcessEvent, because it
    // might be needed by components during the event processing.
    //HLTWarning("did not receive EOR for run %d", fpRunDesc->fRunNo);
    AliHLTRunDesc* pRunDesc=fpRunDesc;
    fpRunDesc=NULL;
    delete pRunDesc;
  }
  if (fpCTPData) {
    delete fpCTPData;
  }
  fpCTPData=NULL;

  fEventCount=0;
  fFlags=0;
  return iResult;
}

int AliHLTComponent::InitCDB(const char* cdbPath, AliHLTComponentHandler* pHandler)
{
  // see header file for function documentation
  int iResult=0;
  HLTInfo("Using CDB: %s", cdbPath);
  if (pHandler) {
  // I have to think about separating the library handling from the
  // component handler. Requiring the component handler here is not
  // the cleanest solution.
  // We presume the library already to be loaded, which is the case
  // because it is loaded in the initialization of the logging functionality
  //
  // find the symbol
  AliHLTMiscInitCDB_t pFunc=(AliHLTMiscInitCDB_t)pHandler->FindSymbol(ALIHLTMISC_LIBRARY, ALIHLTMISC_INIT_CDB);
  if (pFunc) {
    TString path;
    if (cdbPath && cdbPath[0]!=0) {
      path=cdbPath;
      // very temporary fix, have to check for other formats
      if (!path.BeginsWith("local://")) {
	path="local://";
	path+=cdbPath;
      }
    }
    if ((iResult=(*pFunc)(path.Data()))>=0) {
      if (!(fCDBSetRunNoFunc=pHandler->FindSymbol(ALIHLTMISC_LIBRARY, ALIHLTMISC_SET_CDB_RUNNO))) {
	Message(NULL, kHLTLogWarning, "AliHLTComponent::InitCDB", "init CDB",
		"can not find function to set CDB run no");
      }
    }
  } else {
    Message(NULL, kHLTLogError, "AliHLTComponent::InitCDB", "init CDB",
	    "can not find initialization function");
    iResult=-ENOSYS;
  }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTComponent::SetCDBRunNo(int runNo)
{
  // see header file for function documentation
  if (!fCDBSetRunNoFunc) return 0;
  return (*((AliHLTMiscSetCDBRunNo_t)fCDBSetRunNoFunc))(runNo);
}

int AliHLTComponent::SetRunDescription(const AliHLTRunDesc* desc, const char* /*runType*/)
{
  // see header file for function documentation
  if (!desc) return -EINVAL;
  if (desc->fStructSize!=sizeof(AliHLTRunDesc)) {
    HLTError("invalid size of RunDesc struct (%ul)", desc->fStructSize);
    return -EINVAL;
  }

  if (!fpRunDesc) {
    fpRunDesc=new AliHLTRunDesc;
    if (!fpRunDesc) return -ENOMEM;
    *fpRunDesc=kAliHLTVoidRunDesc;
  }

  if (fpRunDesc->fRunNo!=kAliHLTVoidRunNo && fpRunDesc->fRunNo!=desc->fRunNo) {
    HLTWarning("Run description has already been set");
  }
  *fpRunDesc=*desc;
  SetCDBRunNo(fpRunDesc->fRunNo);
  // TODO: we have to decide about the runType
  return 0;
}

int AliHLTComponent::SetComponentDescription(const char* desc)
{
  // see header file for function documentation
  int iResult=0;
  if (!desc) return 0;

  TString descriptor=desc;
  TObjArray* pTokens=descriptor.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntriesFast() && iResult>=0; i++) {
      TString argument=pTokens->At(i++)->GetName();
      if (!argument || argument.IsNull()) continue;

      // chainid
      if (argument.BeginsWith("chainid")) {
	argument.ReplaceAll("chainid", "");
	if (argument.BeginsWith("=")) {
	  fChainId=argument.Replace(0,1,"");
	  fChainIdCrc=CalculateChecksum((const AliHLTUInt8_t*)fChainId.c_str(), fChainId.length());
	  HLTDebug("setting component description: chain id %s crc 0x%8x", fChainId.c_str(), fChainIdCrc);
	} else {
	  fChainId="";
	}
      } else {
	HLTWarning("unknown component description %s", argument.Data());
      }
    }
    delete pTokens;
  }
  
  return iResult;
}

int AliHLTComponent::ConfigureFromArgumentString(int argc, const char** argv)
{
  // Configure from an array of argument strings
  // Function supports individual arguments in the argv array and arguments
  // separated by blanks.
  //
  // Each argv entry can contain blanks, quotes and newlines. Newlines are interpreted
  // as blanks. Enclosing quotes deactivate blank as delimiter.
  // The separated arguments are stored in an array of strings, and the pointers to
  // those strings in an array of pointers. The latter is used in the custom argument
  // scan of the component.

  int iResult=0;
  vector<string> stringarray;   // array of argument copies
  // array of pointers to the argument copies
  // note: not necessarily in sync with the entries in stringarray
  // can contain any pointer to valid arguments in arbitrary sequence
  vector<const char*> ptrarray;

  TString argument="";
  int i=0;
  for (i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsWhitespace()) continue;

    // special handling for single component arguments ending with
    // a sequence of blanks. All characters until the first occurence
    // of a blank are removed. If the remainder contains only whitespaces
    // the argument is a single argument, just having whitespaces at the end.
    argument.Remove(0, argument.First(' '));
    if (argument.IsWhitespace()) {
      ptrarray.push_back(argv[i]);
      continue;
    }

    // outer loop checks for enclosing quotes
    // extra blank to insert blank token before leading quotes, then
    // quoted arguments are always the even ones
    argument=" ";
    argument+=argv[i];
    // insert blank in consecutive quotes to correctly tokenize
    argument.ReplaceAll("''", "' '");
    // replace newlines by blanks
    argument.ReplaceAll("\n", " ");
    if (argument.IsNull()) continue;
    TObjArray* pTokensQuote=argument.Tokenize("'");
    if (pTokensQuote) {
      if (pTokensQuote->GetEntriesFast()>0) {
	for (int k=0; k<pTokensQuote->GetEntriesFast(); k++) {
	  argument=pTokensQuote->At(k)->GetName();
	  if (argument.IsWhitespace()) continue;
	  if (k%2) {
	    // every second entry is enclosed by quotes and thus
	    // one single argument
	    stringarray.push_back(argument.Data());
	    ptrarray.push_back(stringarray.back().c_str());
	  } else {
    TObjArray* pTokens=argument.Tokenize(" ");
    if (pTokens) {
      if (pTokens->GetEntriesFast()>0) {
	for (int n=0; n<pTokens->GetEntriesFast(); n++) {
	  TString data=pTokens->At(n)->GetName();
	  if (!data.IsNull() && !data.IsWhitespace()) {
	    stringarray.push_back(data.Data());
	    ptrarray.push_back(stringarray.back().c_str());
	  }
	}
      }
      delete pTokens;
    }
	  }
	}
      }
      delete pTokensQuote;
    }
  }

  for (i=0; (unsigned)i<ptrarray.size() && iResult>=0;) {
    int result=ScanConfigurationArgument(ptrarray.size()-i, &ptrarray[i]);
    if (result==0) {
      HLTWarning("unknown component argument %s", ptrarray[i]);
      i++;
    } else if (result>0) {
      i+=result;
    } else {
      iResult=result;
      if (iResult==-EINVAL) {
	HLTError("unknown argument %s", ptrarray[i]);
      } else if (iResult==-EPROTO) {
	HLTError("missing/wrong parameter for argument %s (%s)", ptrarray[i], (ptrarray.size()>(unsigned)i+1)?ptrarray[i+1]:"missing");
      } else {
	HLTError("scan of argument %s failed (%d)", ptrarray[i], iResult);
      }
    }
  }

  return iResult;
}

int AliHLTComponent::ConfigureFromCDBTObjString(const char* entries, const char* key)
{
  // load a list of OCDB objects and configure from the objects
  // can either be a TObjString or a TMap with a TObjString:TObjString key-value pair
  int iResult=0;
  TString arguments;
  TString confEntries=entries;
  TObjArray* pTokens=confEntries.Tokenize(" ");
  if (pTokens) {
    for (int n=0; n<pTokens->GetEntriesFast(); n++) {
      const char* path=pTokens->At(n)->GetName();
      const char* chainId=GetChainId();
      HLTInfo("configure from entry \"%s\"%s%s, chain id %s", path, key?" key ":"",key?key:"", (chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
      TObject* pOCDBObject = LoadAndExtractOCDBObject(path, key);
      if (pOCDBObject) {
	TObjString* pString=dynamic_cast<TObjString*>(pOCDBObject);
	if (!pString) {
	  TMap* pMap=dynamic_cast<TMap*>(pOCDBObject);
	  if (pMap) {
	    // this is the case where no key has been specified and the OCDB
	    // object is a TMap, search for the default key
	    TObject* pObject=pMap->GetValue("default");
	    if (pObject && (pString=dynamic_cast<TObjString*>(pObject))!=NULL) {
	      HLTInfo("using default key of TMap of configuration object \"%s\"", path);
	    } else {
	      HLTError("no default key available in TMap of configuration object \"%s\"", path);
	      iResult=-ENOENT;
	      break;
	    }
	  }
	}

	if (pString) {
	  HLTInfo("received configuration object string: \'%s\'", pString->GetName());
	  arguments+=pString->GetName();
	  arguments+=" ";
	} else {
	  HLTError("configuration object \"%s\"%s%s has wrong type, required TObjString", path, key?" key ":"",key?key:"");
	  iResult=-EINVAL;
	}
      } else {
	HLTError("can not fetch object \"%s\" from OCDB", path);
	iResult=-ENOENT;
      }
    }
    delete pTokens;
  }
  if (iResult>=0 && !arguments.IsNull())  {
    const char* array=arguments.Data();
    iResult=ConfigureFromArgumentString(1, &array);
  }
  return iResult;
}

TObject* AliHLTComponent::LoadAndExtractOCDBObject(const char* path, const char* key) const
{
  // see header file for function documentation
  AliCDBEntry* pEntry=AliHLTMisc::Instance().LoadOCDBEntry(path, GetRunNo());
  if (!pEntry) return NULL;
  TObject* pObject=AliHLTMisc::Instance().ExtractObject(pEntry);
  TMap* pMap=dynamic_cast<TMap*>(pObject);
  if (pMap && key) {
    pObject=pMap->GetValue(key);
    if (!pObject) {
      pObject=pMap->GetValue("default");
      if (pObject) {
	HLTWarning("can not find object for key \"%s\" in TMap of configuration object \"%s\", using key \"default\"", key, path);
      }
    }
    if (!pObject) {
      HLTError("can not find object for key \"%s\" in TMap of configuration object \"%s\"", key, path);
      return NULL;
    }
  }
  return pObject;
}

int AliHLTComponent::DoInit( int /*argc*/, const char** /*argv*/)
{
  // default implementation, childs can overload
  HLTLogKeyword("dummy");
  return 0;
}

int AliHLTComponent::DoDeinit()
{
  // default implementation, childs can overload
  HLTLogKeyword("dummy");
  return 0;
}

int AliHLTComponent::Reconfigure(const char* /*cdbEntry*/, const char* /*chainId*/)
{
  // default implementation, childs can overload
  HLTLogKeyword("dummy");
  return 0;
}

int AliHLTComponent::ReadPreprocessorValues(const char* /*modules*/)
{
  // default implementation, childs can overload
  HLTLogKeyword("dummy");
  return 0;
}

int AliHLTComponent::ScanConfigurationArgument(int /*argc*/, const char** /*argv*/)
{
  // default implementation, childs can overload
  HLTLogKeyword("dummy");
  HLTWarning("The function needs to be implemented by the component");
  return 0;
}

int AliHLTComponent::StartOfRun()
{
  // default implementation, childs can overload
  HLTLogKeyword("dummy");
  return 0;
}

int AliHLTComponent::EndOfRun()
{
  // default implementation, childs can overload
  HLTLogKeyword("dummy");
  return 0;
}


int AliHLTComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& /*tgtList*/)
{
  // default implementation, childs can overload
  HLTLogKeyword("dummy");
  return 0;
}

void AliHLTComponent::GetOCDBObjectDescription( TMap* const /*targetArray*/)
{
  // default implementation, childs can overload
  HLTLogKeyword("dummy");
}

int AliHLTComponent::CheckOCDBEntries(const TMap* const externList)
{
  // check the availability of the OCDB entry descriptions in the TMap
  //  key : complete OCDB path of the entry
  //  value : auxiliary object - short description
  // if the external map was not provided the function invokes
  // interface function GetOCDBObjectDescription() to retrieve the list.
  int iResult=0;
  if (externList) {
    iResult=AliHLTMisc::Instance().CheckOCDBEntries(externList);
  } else {
    TMap* pMap=new TMap;
    if (pMap) {
      pMap->SetOwnerKeyValue(kTRUE);
      GetOCDBObjectDescription(pMap);
      iResult=AliHLTMisc::Instance().CheckOCDBEntries(pMap);
      delete pMap;
      pMap=NULL;
    }
  }

  return iResult;
}

void AliHLTComponent::DataType2Text( const AliHLTComponentDataType& type, char output[kAliHLTComponentDataTypefIDsize+kAliHLTComponentDataTypefOriginSize+2] ) const
{
  // see header file for function documentation
  memset( output, 0, kAliHLTComponentDataTypefIDsize+kAliHLTComponentDataTypefOriginSize+2 );
  strncat( output, type.fOrigin, kAliHLTComponentDataTypefOriginSize );
  strncat( output, ":", 1 );
  strncat( output, type.fID, kAliHLTComponentDataTypefIDsize );
}

string AliHLTComponent::DataType2Text( const AliHLTComponentDataType& type, int mode)
{
  // see header file for function documentation
  string out("");

  // 'typeid' 'origin'
  // aligned to 8 and 4 chars respectively, blocks enclosed in quotes and
  // separated by blank e.g.
  // 'DDL_RAW ' 'TPC '
  if (mode==3) {
    int i=0;
    char tmp[8];
    out+="'";
    for (i=0; i<kAliHLTComponentDataTypefIDsize; i++) {
      unsigned char* puc=(unsigned char*)type.fID;
      if (puc[i]<32)
	sprintf(tmp, "\\%x", type.fID[i]);
      else
	sprintf(tmp, "%c", type.fID[i]);
      out+=tmp;
    }
    out+="' '";
    for (i=0; i<kAliHLTComponentDataTypefOriginSize; i++) {
      unsigned char* puc=(unsigned char*)type.fOrigin;
      if ((puc[i])<32)
	sprintf(tmp, "\\%x", type.fOrigin[i]);
      else
	sprintf(tmp, "%c", type.fOrigin[i]);
      out+=tmp;
    }
    out+="'";
    return out;
  }

  // origin typeid as numbers separated by colon e.g.
  // aligned to 8 and 4 chars respectively, all characters separated by
  // quotes, e.g.
  // '84'80'67'32':'68'68'76'95'82'65'87'32'
  if (mode==2) {
    int i=0;
    char tmp[8];
    for (i=0; i<kAliHLTComponentDataTypefOriginSize; i++) {
      sprintf(tmp, "'%d", type.fOrigin[i]);
      out+=tmp;
    }
    out+="':'";
    for (i=0; i<kAliHLTComponentDataTypefIDsize; i++) {
      sprintf(tmp, "%d'", type.fID[i]);
      out+=tmp;
    }
    return out;
  }

  // origin typeid separated by colon e.g.
  // aligned to 8 and 4 chars respectively, all characters separated by
  // quotes, e.g.
  // 'T'P'C' ':'D'D'L'_'R'A'W' '
  if (mode==1) {
    int i=0;
    char tmp[8];
    for (i=0; i<kAliHLTComponentDataTypefOriginSize; i++) {
      unsigned char* puc=(unsigned char*)type.fOrigin;
      if ((puc[i])<32)
	sprintf(tmp, "'\\%x", type.fOrigin[i]);
      else
	sprintf(tmp, "'%c", type.fOrigin[i]);
      out+=tmp;
    }
    out+="':'";
    for (i=0; i<kAliHLTComponentDataTypefIDsize; i++) {
      unsigned char* puc=(unsigned char*)type.fID;
      if (puc[i]<32)
	sprintf(tmp, "\\%x'", type.fID[i]);
      else
	sprintf(tmp, "%c'", type.fID[i]);
      out+=tmp;
    }
    return out;
  }

  // origin typeid
  // aligned to 8 and 4 chars respectively, separated by colon e.g.
  // TPC :DDL_RAW 
  if (type==kAliHLTVoidDataType) {
    out="VOID:VOID";
  } else {
    // some gymnastics in order to avoid a '0' which is part of either or both
    // ID and origin terminating the whole string. Unfortunately, string doesn't
    // stop appending at the '0' if the number of elements to append was 
    // explicitely specified
    string tmp("");
    tmp.append(type.fOrigin, kAliHLTComponentDataTypefOriginSize);
    out.append(tmp.c_str());
    out.append(":");
    tmp="";
    tmp.append(type.fID, kAliHLTComponentDataTypefIDsize);
    out.append(tmp.c_str());
  }
  return out;
}


void* AliHLTComponent::AllocMemory( unsigned long size ) 
{
  // see header file for function documentation
  if (fEnvironment.fAllocMemoryFunc)
    return (*fEnvironment.fAllocMemoryFunc)(fEnvironment.fParam, size );
  HLTFatal("no memory allocation handler registered");
  return NULL;
}

int AliHLTComponent::MakeOutputDataBlockList( const AliHLTComponentBlockDataList& blocks, AliHLTUInt32_t* blockCount,
					      AliHLTComponentBlockData** outputBlocks ) 
{
  // see header file for function documentation
    if ( blockCount==NULL || outputBlocks==NULL )
	return -EFAULT;
    AliHLTUInt32_t count = blocks.size();
    if ( !count )
	{
	*blockCount = 0;
	*outputBlocks = NULL;
	return 0;
	}
    *outputBlocks = reinterpret_cast<AliHLTComponentBlockData*>( AllocMemory( sizeof(AliHLTComponentBlockData)*count ) );
    if ( !*outputBlocks )
	return -ENOMEM;
    for ( unsigned long i = 0; i < count; i++ ) {
	(*outputBlocks)[i] = blocks[i];
	if (MatchExactly(blocks[i].fDataType, kAliHLTAnyDataType)) {
	  (*outputBlocks)[i].fDataType=GetOutputDataType();
	  /* data type was set to the output data type by the PubSub AliRoot
	     Wrapper component, if data type of the block was ********:****.
	     Now handled by the component base class in order to have same
	     behavior when running embedded in AliRoot
	  memset((*outputBlocks)[i].fDataType.fID, '*', kAliHLTComponentDataTypefIDsize);
	  memset((*outputBlocks)[i].fDataType.fOrigin, '*', kAliHLTComponentDataTypefOriginSize);
	  */
	}
    }
    *blockCount = count;
    return 0;

}

int AliHLTComponent::GetEventDoneData( unsigned long size, AliHLTComponentEventDoneData** edd ) const
{
  // see header file for function documentation
  if (fEnvironment.fGetEventDoneDataFunc)
    return (*fEnvironment.fGetEventDoneDataFunc)(fEnvironment.fParam, fCurrentEvent, size, edd );
  return -ENOSYS;
}

int AliHLTComponent::ReserveEventDoneData( unsigned long size )
{
  // see header file for function documentation
  int iResult=0;

  unsigned long capacity=fEventDoneDataSize;
  if (fEventDoneData) capacity-=sizeof(AliHLTComponentEventDoneData)+fEventDoneData->fDataSize;
  if (size>capacity) {
    unsigned long newSize=sizeof(AliHLTComponentEventDoneData)+size+(fEventDoneDataSize-capacity);
    AliHLTComponentEventDoneData* newEDD = reinterpret_cast<AliHLTComponentEventDoneData*>( new AliHLTUInt8_t[newSize] );
    if (!newEDD)
      return -ENOMEM;
    newEDD->fStructSize = sizeof(AliHLTComponentEventDoneData);
    newEDD->fDataSize = 0;
    newEDD->fData = reinterpret_cast<AliHLTUInt8_t*>(newEDD)+newEDD->fStructSize;
    if (fEventDoneData) {
      memcpy( newEDD->fData, fEventDoneData->fData, fEventDoneData->fDataSize );
      newEDD->fDataSize = fEventDoneData->fDataSize;
      delete [] reinterpret_cast<AliHLTUInt8_t*>( fEventDoneData );
    }
    fEventDoneData = newEDD;
    fEventDoneDataSize = newSize;
  }
  return iResult;

}

int AliHLTComponent::PushEventDoneData( AliHLTUInt32_t eddDataWord )
{
  // see header file for function documentation
  if (!fEventDoneData)
    return -ENOMEM;
  if (fEventDoneData->fDataSize+sizeof(AliHLTUInt32_t)>fEventDoneDataSize)
    return -ENOSPC;
  *reinterpret_cast<AliHLTUInt32_t*>((reinterpret_cast<AliHLTUInt8_t*>(fEventDoneData->fData)+fEventDoneData->fDataSize)) = eddDataWord;
  fEventDoneData->fDataSize += sizeof(AliHLTUInt32_t);
  return 0;
}

void AliHLTComponent::ReleaseEventDoneData()
{
   // see header file for function documentation
 if (fEventDoneData)
    delete [] reinterpret_cast<AliHLTUInt8_t*>( fEventDoneData );
  fEventDoneData = NULL;
  fEventDoneDataSize = 0;
}


int AliHLTComponent::FindMatchingDataTypes(AliHLTComponent* pConsumer, AliHLTComponentDataTypeList* tgtList) 
{
  // see header file for function documentation
  int iResult=0;
  if (pConsumer) {
    AliHLTComponentDataTypeList itypes;
    AliHLTComponentDataTypeList otypes;
    otypes.push_back(GetOutputDataType());
    if (MatchExactly(otypes[0],kAliHLTMultipleDataType)) {
      otypes.clear();
      int count=0;
      if ((count=GetOutputDataTypes(otypes))>0) {
      } else if (GetComponentType()!=kSink) {
	HLTWarning("component %s indicates multiple output data types but GetOutputDataTypes returns %d", GetComponentID(), count);
      }
    }
    pConsumer->GetInputDataTypes(itypes);
    AliHLTComponentDataTypeList::iterator otype=otypes.begin();
    for (;otype!=otypes.end();otype++) {
      //PrintDataTypeContent((*otype), "publisher \'%s\'");
      if ((*otype)==(kAliHLTAnyDataType|kAliHLTDataOriginPrivate)) {
	if (tgtList) tgtList->push_back(*otype);
	iResult++;
	continue;
      }
      
      AliHLTComponentDataTypeList::iterator itype=itypes.begin();
      for ( ; itype!=itypes.end() && (*itype)!=(*otype) ; itype++) {/* empty body */};
      //if (itype!=itypes.end()) PrintDataTypeContent(*itype, "consumer \'%s\'");
      if (itype!=itypes.end()) {
	if (tgtList) tgtList->push_back(*otype);
	iResult++;
      }
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

void AliHLTComponent::PrintDataTypeContent(AliHLTComponentDataType& dt, const char* format)
{
  // see header file for function documentation
  const char* fmt="\'%s\'";
  if (format) fmt=format;
  AliHLTLogging::Message(NULL, kHLTLogNone, NULL , NULL, Form(fmt, (DataType2Text(dt)).c_str()));
  AliHLTLogging::Message(NULL, kHLTLogNone, NULL , NULL, 
			 Form("%x %x %x %x %x %x %x %x : %x %x %x %x", 
			      dt.fID[0],
			      dt.fID[1],
			      dt.fID[2],
			      dt.fID[3],
			      dt.fID[4],
			      dt.fID[5],
			      dt.fID[6],
			      dt.fID[7],
			      dt.fOrigin[0],
			      dt.fOrigin[1],
			      dt.fOrigin[2],
			      dt.fOrigin[3]));
}

void AliHLTComponent::FillBlockData( AliHLTComponentBlockData& blockData )
{
  // see header file for function documentation
  blockData.fStructSize = sizeof(blockData);
  FillShmData( blockData.fShmKey );
  blockData.fOffset = ~(AliHLTUInt32_t)0;
  blockData.fPtr = NULL;
  blockData.fSize = 0;
  FillDataType( blockData.fDataType );
  blockData.fSpecification = kAliHLTVoidDataSpec;
}

void AliHLTComponent::FillShmData( AliHLTComponentShmData& shmData )
{
  // see header file for function documentation
  shmData.fStructSize = sizeof(shmData);
  shmData.fShmType = gkAliHLTComponentInvalidShmType;
  shmData.fShmID = gkAliHLTComponentInvalidShmID;
}

void AliHLTComponent::FillDataType( AliHLTComponentDataType& dataType )
{
  // see header file for function documentation
  dataType=kAliHLTAnyDataType;
}

void AliHLTComponent::CopyDataType(AliHLTComponentDataType& tgtdt, const AliHLTComponentDataType& srcdt) 
{
  // see header file for function documentation
  memcpy(&tgtdt.fID[0], &srcdt.fID[0], kAliHLTComponentDataTypefIDsize);
  memcpy(&tgtdt.fOrigin[0], &srcdt.fOrigin[0], kAliHLTComponentDataTypefOriginSize);
}

void AliHLTComponent::SetDataType(AliHLTComponentDataType& tgtdt, const char* id, const char* origin) 
{
  // see header file for function documentation
  tgtdt.fStructSize=sizeof(AliHLTComponentDataType);
  if (id) {
    memset(&tgtdt.fID[0], 0, kAliHLTComponentDataTypefIDsize);
    strncpy(&tgtdt.fID[0], id, strlen(id)<(size_t)kAliHLTComponentDataTypefIDsize?strlen(id):kAliHLTComponentDataTypefIDsize);
  }
  if (origin) {
    memset(&tgtdt.fOrigin[0], 0, kAliHLTComponentDataTypefOriginSize);
    strncpy(&tgtdt.fOrigin[0], origin, strlen(origin)<(size_t)kAliHLTComponentDataTypefOriginSize?strlen(origin):kAliHLTComponentDataTypefOriginSize);
  }
}

void AliHLTComponent::SetDataType(AliHLTComponentDataType& dt, AliHLTUInt64_t id, AliHLTUInt32_t origin)
{
  // see header file for function documentation
  dt.fStructSize=sizeof(AliHLTComponentDataType);
  assert(kAliHLTComponentDataTypefIDsize==sizeof(id));
  assert(kAliHLTComponentDataTypefOriginSize==sizeof(origin));
  memcpy(&dt.fID, &id, kAliHLTComponentDataTypefIDsize);
  memcpy(&dt.fOrigin, &origin, kAliHLTComponentDataTypefOriginSize);
}

void AliHLTComponent::FillEventData(AliHLTComponentEventData& evtData)
{
  // see header file for function documentation
  memset(&evtData, 0, sizeof(AliHLTComponentEventData));
  evtData.fStructSize=sizeof(AliHLTComponentEventData);
  evtData.fEventID=kAliHLTVoidEventID;
}

void AliHLTComponent::PrintComponentDataTypeInfo(const AliHLTComponentDataType& dt) 
{
  // see header file for function documentation
  TString msg;
  msg.Form("AliHLTComponentDataType(%d): ID=\"", dt.fStructSize);
  for ( int i = 0; i < kAliHLTComponentDataTypefIDsize; i++ ) {
   if (dt.fID[i]!=0) msg+=dt.fID[i];
   else msg+="\\0";
  }
  msg+="\" Origin=\"";
  for ( int i = 0; i < kAliHLTComponentDataTypefOriginSize; i++ ) {
   if (dt.fOrigin[i]!=0) msg+=dt.fOrigin[i];
   else msg+="\\0";
  }
  msg+="\"";
  AliHLTLogging::Message(NULL, kHLTLogNone, NULL , NULL, msg.Data());
}

int AliHLTComponent::GetEventCount() const
{
  // see header file for function documentation
  return fEventCount;
}

int AliHLTComponent::IncrementEventCounter()
{
  // see header file for function documentation
  if (fEventCount>=0) fEventCount++;
  return fEventCount;
}

int AliHLTComponent::GetNumberOfInputBlocks() const
{
  // see header file for function documentation
  if (fpInputBlocks!=NULL) {
    return fCurrentEventData.fBlockCnt;
  }
  return 0;
}

AliHLTEventID_t AliHLTComponent::GetEventId() const
{
  // see header file for function documentation
  if (fpInputBlocks!=NULL) {
    return fCurrentEventData.fEventID;
  }
  return 0;
}

const TObject* AliHLTComponent::GetFirstInputObject(const AliHLTComponentDataType& dt,
						    const char* classname,
						    int bForce)
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  fSearchDataType=dt;
  if (classname) fClassName=classname;
  else fClassName.clear();
  int idx=FindInputBlock(fSearchDataType, 0, 1);
  TObject* pObj=NULL;
  if (idx>=0) {
    HLTDebug("found block %d when searching for data type %s", idx, DataType2Text(dt).c_str());
    if ((pObj=GetInputObject(idx, fClassName.c_str(), bForce))!=NULL) {
      fCurrentInputBlock=idx;
    } else {
    }
  }
  return pObj;
}

const TObject* AliHLTComponent::GetFirstInputObject(const char* dtID, 
						    const char* dtOrigin,
						    const char* classname,
						    int         bForce)
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  AliHLTComponentDataType dt;
  SetDataType(dt, dtID, dtOrigin);
  return GetFirstInputObject(dt, classname, bForce);
}

const TObject* AliHLTComponent::GetNextInputObject(int bForce)
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  int idx=FindInputBlock(fSearchDataType, fCurrentInputBlock+1, 1);
  //HLTDebug("found block %d when searching for data type %s", idx, DataType2Text(fSearchDataType).c_str());
  TObject* pObj=NULL;
  if (idx>=0) {
    if ((pObj=GetInputObject(idx, fClassName.c_str(), bForce))!=NULL) {
      fCurrentInputBlock=idx;
    }
  }
  return pObj;
}

int AliHLTComponent::FindInputBlock(const AliHLTComponentDataType& dt, int startIdx, int bObject) const
{
  // see header file for function documentation
  int iResult=-ENOENT;
  if (fpInputBlocks!=NULL) {
    int idx=startIdx<0?0:startIdx;
    for ( ; (UInt_t)idx<fCurrentEventData.fBlockCnt && iResult==-ENOENT; idx++) {
      if (dt!=fpInputBlocks[idx].fDataType) continue;

      if (bObject!=0) {
	if (fpInputBlocks[idx].fPtr==NULL) continue;
	AliHLTUInt32_t firstWord=*((AliHLTUInt32_t*)fpInputBlocks[idx].fPtr);
	if (firstWord!=fpInputBlocks[idx].fSize-sizeof(AliHLTUInt32_t)) continue;
      }
      iResult=idx;
    }
  }
  return iResult;
}

TObject* AliHLTComponent::CreateInputObject(int idx, int bForce)
{
  // see header file for function documentation
  TObject* pObj=NULL;
  if (fpInputBlocks!=NULL) {
    if ((UInt_t)idx<fCurrentEventData.fBlockCnt) {
      if (fpInputBlocks[idx].fPtr) {
	AliHLTUInt32_t firstWord=*((AliHLTUInt32_t*)fpInputBlocks[idx].fPtr);
	if (firstWord==fpInputBlocks[idx].fSize-sizeof(AliHLTUInt32_t)) {
	  HLTDebug("create object from block %d size %d", idx, fpInputBlocks[idx].fSize);
	  AliHLTMessage msg(fpInputBlocks[idx].fPtr, fpInputBlocks[idx].fSize);
	  TClass* objclass=msg.GetClass();
	  pObj=msg.ReadObject(objclass);
	  if (pObj && objclass) {
	    HLTDebug("object %p type %s created", pObj, objclass->GetName());
	  } else {
	  }
	  //} else {
       	} else if (bForce!=0) {
	  HLTError("size mismatch: block size %d, indicated %d", fpInputBlocks[idx].fSize, firstWord+sizeof(AliHLTUInt32_t));
	}
      } else {
	HLTFatal("block descriptor empty");
      }
    } else {
      HLTError("index %d out of range %d", idx, fCurrentEventData.fBlockCnt);
    }
  } else {
    HLTError("no input blocks available");
  }
  
  return pObj;
}

TObject* AliHLTComponent::GetInputObject(int idx, const char* /*classname*/, int bForce)
{
  // see header file for function documentation
  if (fpInputObjects==NULL) {
    fpInputObjects=new TObjArray(fCurrentEventData.fBlockCnt);
  }
  TObject* pObj=NULL;
  if (fpInputObjects) {
    pObj=fpInputObjects->At(idx);
    if (pObj==NULL) {
      pObj=CreateInputObject(idx, bForce);
      if (pObj) {
	fpInputObjects->AddAt(pObj, idx);
      }
    }
  } else {
    HLTFatal("memory allocation failed: TObjArray of size %d", fCurrentEventData.fBlockCnt);
  }
  return pObj;
}

int AliHLTComponent::CleanupInputObjects()
{
  // see header file for function documentation
  if (!fpInputObjects) return 0;
  TObjArray* array=fpInputObjects;
  fpInputObjects=NULL;
  for (int i=0; i<array->GetEntriesFast(); i++) {
    TObject* pObj=array->At(i);
    // grrr, garbage collection strikes back: When read via AliHLTMessage
    // (CreateInputObject), and written to a TFile afterwards, the
    // TFile::Close calls ROOOT's garbage collection. No clue why the
    // object ended up in the key list and needs to be deleted
    //
    // Matthias 09.11.2008 follow up
    // This approach doesn't actually work in all cases: the object table
    // can be switched off globally, the flag needs to be checked here as
    // well in order to avoid memory leaks.
    // This means we have to find another solution for the problem if it
    // pops up again.
    if (pObj &&
	(!TObject::GetObjectStat() || gObjectTable->PtrIsValid(pObj))) {
      delete pObj;
    }
  }
  delete array;
  return 0;
}

AliHLTComponentDataType AliHLTComponent::GetDataType(const TObject* pObject)
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  AliHLTComponentDataType dt=kAliHLTVoidDataType;
  int idx=fCurrentInputBlock;
  if (pObject) {
    if (fpInputObjects==NULL || (idx=fpInputObjects->IndexOf(pObject))>=0) {
    } else {
      HLTError("unknown object %p", pObject);
    }
  }
  if (idx>=0) {
    if ((UInt_t)idx<fCurrentEventData.fBlockCnt) {
      dt=fpInputBlocks[idx].fDataType;
    } else {
      HLTFatal("severe internal error, index out of range");
    }
  }
  return dt;
}

AliHLTUInt32_t AliHLTComponent::GetSpecification(const TObject* pObject)
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  AliHLTUInt32_t iSpec=kAliHLTVoidDataSpec;
  int idx=fCurrentInputBlock;
  if (pObject) {
    if (fpInputObjects==NULL || (idx=fpInputObjects->IndexOf(pObject))>=0) {
    } else {
      HLTError("unknown object %p", pObject);
    }
  }
  if (idx>=0) {
    if ((UInt_t)idx<fCurrentEventData.fBlockCnt) {
      iSpec=fpInputBlocks[idx].fSpecification;
    } else {
      HLTFatal("severe internal error, index out of range");
    }
  }
  return iSpec;
}

int AliHLTComponent::Forward(const TObject* pObject)
{
  // see header file for function documentation
  int iResult=0;
  int idx=fCurrentInputBlock;
  if (pObject) {
    if (fpInputObjects==NULL || (idx=fpInputObjects->IndexOf(pObject))>=0) {
    } else {
      HLTError("unknown object %p", pObject);
      iResult=-ENOENT;
    }
  }
  if (idx>=0) {
    fOutputBlocks.push_back(fpInputBlocks[idx]);
  }
  return iResult;
}

int AliHLTComponent::Forward(const AliHLTComponentBlockData* pBlock)
{
  // see header file for function documentation
  int iResult=0;
  int idx=fCurrentInputBlock;
  if (pBlock) {
    if ((idx=FindInputBlock(pBlock))>=0) {
    } else {
      HLTError("unknown Block %p", pBlock);
      iResult=-ENOENT;      
    }
  }
  if (idx>=0) {
    // check for fpInputBlocks pointer done in FindInputBlock
    fOutputBlocks.push_back(fpInputBlocks[idx]);
  }
  return iResult;
}

const AliHLTComponentBlockData* AliHLTComponent::GetFirstInputBlock(const AliHLTComponentDataType& dt)
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  fSearchDataType=dt;
  fClassName.clear();
  int idx=FindInputBlock(fSearchDataType, 0);
  const AliHLTComponentBlockData* pBlock=NULL;
  if (idx>=0) {
    // check for fpInputBlocks pointer done in FindInputBlock
    pBlock=&fpInputBlocks[idx];
    fCurrentInputBlock=idx;
  }
  return pBlock;
}

const AliHLTComponentBlockData* AliHLTComponent::GetFirstInputBlock(const char* dtID, 
								    const char* dtOrigin)
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  AliHLTComponentDataType dt;
  SetDataType(dt, dtID, dtOrigin);
  return GetFirstInputBlock(dt);
}

const AliHLTComponentBlockData* AliHLTComponent::GetInputBlock(int index) const
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  assert( 0 <= index and index < (int)fCurrentEventData.fBlockCnt );
  return &fpInputBlocks[index];
}

const AliHLTComponentBlockData* AliHLTComponent::GetNextInputBlock()
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  int idx=FindInputBlock(fSearchDataType, fCurrentInputBlock+1);
  const AliHLTComponentBlockData* pBlock=NULL;
  if (idx>=0) {
    // check for fpInputBlocks pointer done in FindInputBlock
    pBlock=&fpInputBlocks[idx];
    fCurrentInputBlock=idx;
  }
  return pBlock;
}

int AliHLTComponent::FindInputBlock(const AliHLTComponentBlockData* pBlock) const
{
  // see header file for function documentation
  int iResult=-ENOENT;
  if (fpInputBlocks!=NULL) {
    if (pBlock) {
      if (pBlock>=fpInputBlocks && pBlock<fpInputBlocks+fCurrentEventData.fBlockCnt) {
	iResult=(int)(pBlock-fpInputBlocks);
      }
    } else {
      iResult=-EINVAL;
    }
  }
  return iResult;
}

AliHLTUInt32_t AliHLTComponent::GetSpecification(const AliHLTComponentBlockData* pBlock)
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  AliHLTUInt32_t iSpec=kAliHLTVoidDataSpec;
  int idx=fCurrentInputBlock;
  if (pBlock) {
    if (fpInputObjects==NULL || (idx=FindInputBlock(pBlock))>=0) {
    } else {
      HLTError("unknown Block %p", pBlock);
    }
  }
  if (idx>=0) {
    // check for fpInputBlocks pointer done in FindInputBlock
    iSpec=fpInputBlocks[idx].fSpecification;
  }
  return iSpec;
}

int AliHLTComponent::PushBack(const TObject* pObject, const AliHLTComponentDataType& dt, AliHLTUInt32_t spec, 
			      void* pHeader, int headerSize)
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  int iResult=0;
  fLastObjectSize=0;
  if (fPushbackPeriod>0) {
    // suppress the output
    TDatime time;
    if (fLastPushBackTime<0 || (int)time.Get()-fLastPushBackTime<fPushbackPeriod) return 0;
  }
  if (pObject) {
    AliHLTMessage msg(kMESS_OBJECT);
    msg.SetCompressionLevel(fCompressionLevel);
    msg.WriteObject(pObject);
    Int_t iMsgLength=msg.Length();
    if (iMsgLength>0) {
      // Matthias Sep 2008
      // NOTE: AliHLTMessage does implement it's own SetLength method
      // which is not architecture independent. The original SetLength
      // stores the size always in network byte order.
      // I'm trying to remember the rational for that, might be that
      // it was just some lack of knowledge. Want to change this, but
      // has to be done carefullt to be backward compatible.
      msg.SetLength(); // sets the length to the first (reserved) word

      // does nothing if the level is 0
      msg.Compress();

      char *mbuf = msg.Buffer();
      if (msg.CompBuffer()) {
	msg.SetLength(); // set once more to have to byte order
	mbuf = msg.CompBuffer();
	iMsgLength = msg.CompLength();
      }
      assert(mbuf!=NULL);
      iResult=InsertOutputBlock(mbuf, iMsgLength, dt, spec, pHeader, headerSize);
      if (iResult>=0) {
	HLTDebug("object %s (%p) size %d compression %d inserted to output", pObject->ClassName(), pObject, iMsgLength, msg.GetCompressionLevel());
      }
      fLastObjectSize=iMsgLength;
    } else {
      HLTError("object serialization failed for object %p", pObject);
      iResult=-ENOMSG;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTComponent::PushBack(const TObject* pObject, const char* dtID, const char* dtOrigin, AliHLTUInt32_t spec,
			      void* pHeader, int headerSize)
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  AliHLTComponentDataType dt;
  SetDataType(dt, dtID, dtOrigin);
  return PushBack(pObject, dt, spec, pHeader, headerSize);
}

int AliHLTComponent::PushBack(const void* pBuffer, int iSize, const AliHLTComponentDataType& dt, AliHLTUInt32_t spec,
			      const void* pHeader, int headerSize)
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  if (fPushbackPeriod>0) {
    // suppress the output
    TDatime time;
    if (fLastPushBackTime<0 || (int)time.Get()-fLastPushBackTime<fPushbackPeriod) return 0;
  }

  return InsertOutputBlock(pBuffer, iSize, dt, spec, pHeader, headerSize);
}

int AliHLTComponent::PushBack(const void* pBuffer, int iSize, const char* dtID, const char* dtOrigin, AliHLTUInt32_t spec,
			      const void* pHeader, int headerSize)
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  AliHLTComponentDataType dt;
  SetDataType(dt, dtID, dtOrigin);
  return PushBack(pBuffer, iSize, dt, spec, pHeader, headerSize);
}

int AliHLTComponent::InsertOutputBlock(const void* pBuffer, int iBufferSize, const AliHLTComponentDataType& dt, AliHLTUInt32_t spec,
				       const void* pHeader, int iHeaderSize)
{
  // see header file for function documentation
  int iResult=0;
  int iBlkSize = iBufferSize + iHeaderSize;

  if ((pBuffer!=NULL && iBufferSize>0) || (pHeader!=NULL && iHeaderSize>0)) {
    if (fpOutputBuffer && iBlkSize<=(int)(fOutputBufferSize-fOutputBufferFilled)) {
      AliHLTUInt8_t* pTgt=fpOutputBuffer+fOutputBufferFilled;

      // copy header if provided but skip if the header is the target location
      // in that case it has already been copied
      if (pHeader!=NULL && pHeader!=pTgt) {
	memcpy(pTgt, pHeader, iHeaderSize);
      }

      pTgt += (AliHLTUInt8_t) iHeaderSize;

      // copy buffer if provided but skip if buffer is the target location
      // in that case it has already been copied
      if (pBuffer!=NULL && pBuffer!=pTgt) {
	memcpy(pTgt, pBuffer, iBufferSize);
	
	//AliHLTUInt32_t firstWord=*((AliHLTUInt32_t*)pBuffer);	
	//HLTDebug("copy %d bytes from %p to output buffer %p, first word %#x", iBufferSize, pBuffer, pTgt, firstWord);
      }
      //HLTDebug("buffer inserted to output: size %d data type %s spec %#x", iBlkSize, DataType2Text(dt).c_str(), spec);
    } else {
      if (fpOutputBuffer) {
	HLTError("too little space in output buffer: %d of %d, required %d", fOutputBufferSize-fOutputBufferFilled, fOutputBufferSize, iBlkSize);
      } else {
	HLTError("output buffer not available");
      }
      iResult=-ENOSPC;
    }
  }
  if (iResult>=0) {
    AliHLTComponentBlockData bd;
    FillBlockData( bd );
    bd.fOffset        = fOutputBufferFilled;
    bd.fSize          = iBlkSize;
    bd.fDataType      = dt;
    bd.fSpecification = spec;
    fOutputBlocks.push_back( bd );
    fOutputBufferFilled+=bd.fSize;
  }

  return iResult;
}

int AliHLTComponent::EstimateObjectSize(const TObject* pObject) const
{
  // see header file for function documentation
  if (!pObject) return 0;

  AliHLTMessage msg(kMESS_OBJECT);
  msg.WriteObject(pObject);
  return msg.Length();  
}

AliHLTMemoryFile* AliHLTComponent::CreateMemoryFile(int capacity, const char* dtID,
						    const char* dtOrigin,
						    AliHLTUInt32_t spec)
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  AliHLTComponentDataType dt;
  SetDataType(dt, dtID, dtOrigin);
  return CreateMemoryFile(capacity, dt, spec);
}

AliHLTMemoryFile* AliHLTComponent::CreateMemoryFile(int capacity,
						    const AliHLTComponentDataType& dt,
						    AliHLTUInt32_t spec)
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  AliHLTMemoryFile* pFile=NULL;
  if (capacity>=0 && static_cast<unsigned int>(capacity)<=fOutputBufferSize-fOutputBufferFilled){
    AliHLTUInt8_t* pTgt=fpOutputBuffer+fOutputBufferFilled;
    pFile=new AliHLTMemoryFile((char*)pTgt, capacity);
    if (pFile) {
      unsigned int nofBlocks=fOutputBlocks.size();
      if (nofBlocks+1>fMemFiles.size()) {
	fMemFiles.resize(nofBlocks+1, NULL);
      }
      if (nofBlocks<fMemFiles.size()) {
	fMemFiles[nofBlocks]=pFile;
	AliHLTComponentBlockData bd;
	FillBlockData( bd );
	bd.fOffset        = fOutputBufferFilled;
	bd.fSize          = capacity;
	bd.fDataType      = dt;
	bd.fSpecification = spec;
	fOutputBufferFilled+=bd.fSize;
	fOutputBlocks.push_back( bd );
      } else {
	HLTError("can not allocate/grow object array");
	pFile->CloseMemoryFile(0);
	delete pFile;
	pFile=NULL;
      }
    }
  } else {
    HLTError("can not create memory file of size %d (%d available)", capacity, fOutputBufferSize-fOutputBufferFilled);
  }
  return pFile;
}

AliHLTMemoryFile* AliHLTComponent::CreateMemoryFile(const char* dtID,
						    const char* dtOrigin,
						    AliHLTUInt32_t spec,
						    float capacity)
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  AliHLTComponentDataType dt;
  SetDataType(dt, dtID, dtOrigin);
  int size=fOutputBufferSize-fOutputBufferFilled;
  if (capacity<0 || capacity>1.0) {
    HLTError("invalid parameter: capacity %f", capacity);
    return NULL;
  }
  size=(int)(size*capacity);
  return CreateMemoryFile(size, dt, spec);
}

AliHLTMemoryFile* AliHLTComponent::CreateMemoryFile(const AliHLTComponentDataType& dt,
						    AliHLTUInt32_t spec,
						    float capacity)
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  int size=fOutputBufferSize-fOutputBufferFilled;
  if (capacity<0 || capacity>1.0) {
    HLTError("invalid parameter: capacity %f", capacity);
    return NULL;
  }
  size=(int)(size*capacity);
  return CreateMemoryFile(size, dt, spec);
}

int AliHLTComponent::Write(AliHLTMemoryFile* pFile, const TObject* pObject,
			   const char* key, int option)
{
  // see header file for function documentation
  int iResult=0;
  if (pFile && pObject) {
    pFile->cd();
    iResult=pObject->Write(key, option);
    if (iResult>0) {
      // success
    } else {
      iResult=-pFile->GetErrno();
      if (iResult==-ENOSPC) {
	HLTError("error writing memory file, buffer too small");
      }
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTComponent::CloseMemoryFile(AliHLTMemoryFile* pFile)
{
  // see header file for function documentation
  int iResult=0;
  if (pFile) {
    AliHLTMemoryFilePList::iterator element=fMemFiles.begin();
    int i=0;
    while (element!=fMemFiles.end() && iResult>=0) {
      if (*element && *element==pFile) {
	iResult=pFile->CloseMemoryFile();
	
	// sync memory files and descriptors
	if (iResult>=0) {
	  fOutputBlocks[i].fSize=(*element)->GetSize()+(*element)->GetHeaderSize();
	}
	delete *element;
	*element=NULL;
	return iResult;
      }
      element++; i++;
    }
    HLTError("can not find memory file %p", pFile);
    iResult=-ENOENT;
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTComponent::CreateEventDoneData(AliHLTComponentEventDoneData edd)
{
  // see header file for function documentation
  int iResult=0;

  AliHLTComponentEventDoneData* newEDD = NULL;
  
  unsigned long newSize=edd.fDataSize;
  if (fEventDoneData)
    newSize += fEventDoneData->fDataSize;

  if (newSize>fEventDoneDataSize) {
    newEDD = reinterpret_cast<AliHLTComponentEventDoneData*>( new AliHLTUInt8_t[ sizeof(AliHLTComponentEventDoneData)+newSize ] );
    if (!newEDD)
      return -ENOMEM;
    newEDD->fStructSize = sizeof(AliHLTComponentEventDoneData);
    newEDD->fDataSize = newSize;
    newEDD->fData = reinterpret_cast<AliHLTUInt8_t*>(newEDD)+newEDD->fStructSize;
    unsigned long long offset = 0;
    if (fEventDoneData) {
      memcpy( newEDD->fData, fEventDoneData->fData, fEventDoneData->fDataSize );
      offset += fEventDoneData->fDataSize;
    }
    memcpy( reinterpret_cast<AliHLTUInt8_t*>(newEDD->fData)+offset, edd.fData, edd.fDataSize );
    if (fEventDoneData)
      delete [] reinterpret_cast<AliHLTUInt8_t*>( fEventDoneData );
    fEventDoneData = newEDD;
    fEventDoneDataSize = newSize;
  }
  else if (fEventDoneData) {
    memcpy( reinterpret_cast<AliHLTUInt8_t*>(fEventDoneData->fData)+fEventDoneData->fDataSize, edd.fData, edd.fDataSize );
    fEventDoneData->fDataSize += edd.fDataSize;
  }
  else {
    HLTError("internal mismatch, fEventDoneData=%d but buffer is NULL", fEventDoneDataSize);
    iResult=-EFAULT;
  }
  return iResult;
}

namespace
{
  // helper function for std:sort, implements an operator<
  bool SortComponentStatisticsById(const AliHLTComponentStatistics& a, const AliHLTComponentStatistics& b)
  {
    return a.fId<b.fId;
  }

  // helper function for std:sort
  bool SortComponentStatisticsDescendingByLevel(const AliHLTComponentStatistics& a, const AliHLTComponentStatistics& b)
  {
    return a.fId>b.fId;
  }

  // helper class to define operator== between AliHLTComponentStatistics and AliHLTComponentStatistics.fId
  class AliHLTComponentStatisticsId {
  public:
    AliHLTComponentStatisticsId(AliHLTUInt32_t id) : fId(id) {}
    AliHLTComponentStatisticsId(const AliHLTComponentStatisticsId& src) : fId(src.fId) {}
    bool operator==(const AliHLTComponentStatistics& a) const {return a.fId==fId;}
  private:
    AliHLTComponentStatisticsId();
    AliHLTUInt32_t fId;
  };

  // operator for std::find of AliHLTComponentStatistics by id
  bool operator==(const AliHLTComponentStatistics& a, const AliHLTComponentStatisticsId& b)
  {
    return b==a;
  }

  bool AliHLTComponentStatisticsCompareIds(const AliHLTComponentStatistics& a, const AliHLTComponentStatistics& b)
  {
    return a.fId==b.fId;
  }

  // helper class to define operator== between AliHLTComponentBlockData and AliHLTComponentBlockData.fSpecification
  class AliHLTComponentBlockDataSpecification {
  public:
    AliHLTComponentBlockDataSpecification(AliHLTUInt32_t specification) : fSpecification(specification) {}
    AliHLTComponentBlockDataSpecification(const AliHLTComponentBlockDataSpecification& src) : fSpecification(src.fSpecification) {}
    bool operator==(const AliHLTComponentBlockData& bd) const {return bd.fSpecification==fSpecification;}
  private:
    AliHLTComponentBlockDataSpecification();
    AliHLTUInt32_t fSpecification;
  };

  // operator for std::find of AliHLTComponentBlockData by specification
  bool operator==(const AliHLTComponentBlockData& bd, const AliHLTComponentBlockDataSpecification& spec)
  {
    return spec==bd;
  }

  // operator for std::find
  bool operator==(const AliHLTComponentBlockData& a, const AliHLTComponentBlockData& b)
  {
    if (!MatchExactly(a.fDataType,b.fDataType)) return false;
    return a.fSpecification==b.fSpecification;
  }

} // end of namespace

int AliHLTComponent::ProcessEvent( const AliHLTComponentEventData& evtData,
				   const AliHLTComponentBlockData* blocks, 
				   AliHLTComponentTriggerData& trigData,
				   AliHLTUInt8_t* outputPtr, 
				   AliHLTUInt32_t& size,
				   AliHLTUInt32_t& outputBlockCnt, 
				   AliHLTComponentBlockData*& outputBlocks,
				   AliHLTComponentEventDoneData*& edd )
{
  // see header file for function documentation
  HLTLogKeyword(fChainId.c_str());
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  int iResult=0;
  fCurrentEvent=evtData.fEventID;
  fCurrentEventData=evtData;
  fpInputBlocks=blocks;
  fCurrentInputBlock=-1;
  fSearchDataType=kAliHLTAnyDataType;
  fpOutputBuffer=outputPtr;
  fOutputBufferSize=size;
  fOutputBufferFilled=0;
  fOutputBlocks.clear();
  outputBlockCnt=0;
  outputBlocks=NULL;

  AliHLTComponentBlockDataList forwardedBlocks;

  // optional component statistics
  AliHLTComponentStatisticsList compStats;
  bool bAddComponentTableEntry=false;
  vector<AliHLTUInt32_t> parentComponentTables;
  int processingLevel=-1;
#if defined(HLT_COMPONENT_STATISTICS)
  if ((fFlags&kDisableComponentStat)==0) {
    AliHLTComponentStatistics outputStat;
    memset(&outputStat, 0, sizeof(AliHLTComponentStatistics));
    outputStat.fStructSize=sizeof(AliHLTComponentStatistics);
    outputStat.fId=fChainIdCrc;
    if (fpBenchmark) {
      fpBenchmark->Stop();
      outputStat.fComponentCycleTime=(AliHLTUInt32_t)(fpBenchmark->RealTime()*ALIHLTCOMPONENT_STATTIME_SCALER);
      fpBenchmark->Reset();
      fpBenchmark->Start();
    }
    compStats.push_back(outputStat);
  }
#endif // HLT_COMPONENT_STATISTICS

  // data processing is skipped
  // -  if there are only steering events in the block list.
  //    For the sake of data source components data processing
  //    is not skipped if there is no block list at all or if it
  //    just contains the eventType block
  // - always skipped if the event is of type
  //   - gkAliEventTypeConfiguration
  //   - gkAliEventTypeReadPreprocessor
  const unsigned int skipModeDefault=0x1;
  const unsigned int skipModeForce=0x2;
  unsigned int bSkipDataProcessing=skipModeDefault;

  // find special events
  if (fpInputBlocks && evtData.fBlockCnt>0) {
    // first look for all special events and execute in the appropriate
    // sequence afterwords
    int indexComConfEvent=-1;
    int indexUpdtDCSEvent=-1;
    int indexSOREvent=-1;
    int indexEOREvent=-1;
    int indexECSParamBlock=-1;
    for (unsigned int i=0; i<evtData.fBlockCnt && iResult>=0; i++) {
      if (fpInputBlocks[i].fDataType==kAliHLTDataTypeSOR) {
	indexSOREvent=i;
	// the AliHLTCalibrationProcessor relies on the SOR and EOR events
	bSkipDataProcessing&=~skipModeDefault;
      } else if (fpInputBlocks[i].fDataType==kAliHLTDataTypeRunType) {
	// run type string
	// handling is not clear yet
	if (fpInputBlocks[i].fPtr) {
	  HLTDebug("got run type \"%s\"\n", fpInputBlocks[i].fPtr);
	}
      } else if (fpInputBlocks[i].fDataType==kAliHLTDataTypeEOR) {
	indexEOREvent=i;
	// the calibration processor relies on the SOR and EOR events
	bSkipDataProcessing&=~skipModeDefault;
      } else if (fpInputBlocks[i].fDataType==kAliHLTDataTypeDDL) {
	// DDL list
	// this event is most likely deprecated
      } else if (fpInputBlocks[i].fDataType==kAliHLTDataTypeComConf) {
	indexComConfEvent=i;
      } else if (fpInputBlocks[i].fDataType==kAliHLTDataTypeUpdtDCS) {
	indexUpdtDCSEvent=i;
      } else if (fpInputBlocks[i].fDataType==kAliHLTDataTypeEvent) {
	fEventType=fpInputBlocks[i].fSpecification;
	if (fEventType != gkAliEventTypeConfiguration and
	    fEventType != gkAliEventTypeReadPreprocessor
	   )
	{
	  // We can actually get the event type from the CDH if it is valid.
	  // Otherwise just use the specification of the input block.
	  const AliRawDataHeader* cdh;
	  if (ExtractTriggerData(trigData, NULL, NULL, &cdh, NULL) == 0)
	  {
	    fEventType = ExtractEventTypeFromCDH(cdh);
	  }
	}

	// skip always in case of gkAliEventTypeConfiguration
	if (fpInputBlocks[i].fSpecification==gkAliEventTypeConfiguration) bSkipDataProcessing|=skipModeForce;

	// skip always in case of gkAliEventTypeReadPreprocessor
	if (fpInputBlocks[i].fSpecification==gkAliEventTypeReadPreprocessor) bSkipDataProcessing|=skipModeForce;

	// never skip if the event type block is the only block
	if (evtData.fBlockCnt==1) bSkipDataProcessing&=~skipModeDefault;

      } else if (fpInputBlocks[i].fDataType==kAliHLTDataTypeComponentStatistics) {
	if (compStats.size()>0) {
	  AliHLTUInt8_t* pData=reinterpret_cast<AliHLTUInt8_t*>(fpInputBlocks[i].fPtr);
	  for (AliHLTUInt32_t offset=0;
	       offset+sizeof(AliHLTComponentStatistics)<=fpInputBlocks[i].fSize;
	       offset+=sizeof(AliHLTComponentStatistics)) {
	    AliHLTComponentStatistics* pStat=reinterpret_cast<AliHLTComponentStatistics*>(pData+offset);
	    if (pStat && compStats[0].fLevel<=pStat->fLevel) {
	      compStats[0].fLevel=pStat->fLevel+1;
	    }
	    if (find(compStats.begin(), compStats.end(), AliHLTComponentStatisticsId(pStat->fId))==compStats.end()) {
	      compStats.push_back(*pStat);
	    }
	  }
	}
      } else if (fpInputBlocks[i].fDataType==kAliHLTDataTypeComponentTable) {
	AliHLTComponentBlockDataList::iterator element=forwardedBlocks.begin();
	while ((element=find(element, forwardedBlocks.end(), AliHLTComponentBlockDataSpecification(fpInputBlocks[i].fSpecification)))!=forwardedBlocks.end()) {
	  if (element->fDataType==fpInputBlocks[i].fDataType) break;
	  // TODO: think about some more checks inclusing also the actual data buffers
	  // this has to handle multiplicity>1 in the online system, where all components
	  // send the information on SOR, because this event is broadcasted
	}
	if (element==forwardedBlocks.end()) {
	  forwardedBlocks.push_back(fpInputBlocks[i]);
	}
	parentComponentTables.push_back(fpInputBlocks[i].fSpecification);
	if (fpInputBlocks[i].fSize>=sizeof(AliHLTComponentTableEntry)) {
	  const AliHLTComponentTableEntry* entry=reinterpret_cast<AliHLTComponentTableEntry*>(fpInputBlocks[i].fPtr);
	  if (entry->fStructSize==sizeof(AliHLTComponentTableEntry)) {
	    if (processingLevel<0 || processingLevel<=(int)entry->fLevel) 
	      processingLevel=entry->fLevel+1;
	  }
	}
      } else if (fpInputBlocks[i].fDataType==kAliHLTDataTypeECSParam) {
	indexECSParamBlock=i;
      } else {
	// the processing function is called if there is at least one
	// non-steering data block. Steering blocks are not filtered out
	// for sake of performance 
	bSkipDataProcessing&=~skipModeDefault;
	if (compStats.size()>0) {
	  compStats[0].fInputBlockCount++;
	  compStats[0].fTotalInputSize+=fpInputBlocks[i].fSize;
	}
      }
    }

    if (indexSOREvent>=0) {
      // start of run
      bAddComponentTableEntry=true;
      compStats.clear();   // no component statistics for SOR
      if (fpRunDesc==NULL) {
	fpRunDesc=new AliHLTRunDesc;
	if (fpRunDesc) *fpRunDesc=kAliHLTVoidRunDesc;
      }
      if (fpRunDesc) {
	AliHLTRunDesc rundesc;
	if ((iResult=CopyStruct(&rundesc, sizeof(AliHLTRunDesc), indexSOREvent, "AliHLTRunDesc", "SOR"))>0) {
	  if (fpRunDesc->fRunNo==kAliHLTVoidRunNo) {
	    *fpRunDesc=rundesc;
	    HLTDebug("set run decriptor, run no %d", fpRunDesc->fRunNo);
	    SetCDBRunNo(fpRunDesc->fRunNo);
	  } else if (fpRunDesc->fRunNo!=rundesc.fRunNo) {
	    HLTWarning("already set run properties run no %d, ignoring SOR with run no %d", fpRunDesc->fRunNo, rundesc.fRunNo);
	  }
	}
      } else {
	iResult=-ENOMEM;
      }

      if (indexECSParamBlock>=0) {
	if (fpInputBlocks[indexECSParamBlock].fSize>0) {
	  const char* param=reinterpret_cast<const char*>(fpInputBlocks[indexECSParamBlock].fPtr);
	  TString paramString;
	  if (param[fpInputBlocks[indexECSParamBlock].fSize-1]!=0) {
	    HLTWarning("ECS parameter string not terminated");
	    paramString.Insert(0, param, fpInputBlocks[indexECSParamBlock].fSize);
	    paramString+="";
	  } else {
	    paramString=param;
	  }
	  ScanECSParam(paramString.Data());
	} else {
	  HLTWarning("empty ECS parameter received");
	}
      } else {
	// TODO: later on we might throw a warning here since the CTP trigger classes
	// should be mandatory
      }
    }
    if (indexEOREvent>=0) {
      fLastPushBackTime=0; // always send at EOR
      bAddComponentTableEntry=true;
      compStats.clear();   // no component statistics for EOR
      if (fpRunDesc!=NULL) {
	if (fpRunDesc) {
	  AliHLTRunDesc rundesc;
	  if ((iResult=CopyStruct(&rundesc, sizeof(AliHLTRunDesc), indexEOREvent, "AliHLTRunDesc", "SOR"))>0) {
	    if (fpRunDesc->fRunNo!=rundesc.fRunNo) {
	      HLTWarning("run no mismatch: SOR %d, EOR %d", fpRunDesc->fRunNo, rundesc.fRunNo);
	    } else {
	      HLTDebug("EOR run no %d", fpRunDesc->fRunNo);
	    }
	  }
	  // we do not unload the fpRunDesc struct here in order to have the run information
	  // available during the event processing
	  // https://savannah.cern.ch/bugs/?39711
	  // the info will be cleared in DeInit
	}
      } else {
	HLTWarning("did not receive SOR, ignoring EOR");
      }
    }
    if (fEventType==gkAliEventTypeConfiguration) {
      if (indexComConfEvent>=0) {
      TString cdbEntry;
      if (indexComConfEvent>=0 && fpInputBlocks[indexComConfEvent].fPtr!=NULL && fpInputBlocks[indexComConfEvent].fSize>0) {
	cdbEntry.Append(reinterpret_cast<const char*>(fpInputBlocks[indexComConfEvent].fPtr), fpInputBlocks[indexComConfEvent].fSize);
      }
      HLTDebug("received component configuration command: entry %s", cdbEntry.IsNull()?"none":cdbEntry.Data());
      int tmpResult=Reconfigure(cdbEntry[0]==0?NULL:cdbEntry.Data(), fChainId.c_str());
      if (tmpResult<0) {
	HLTWarning("reconfiguration of component %p (%s) failed with error code %d", this, GetComponentID(), tmpResult);
      }
      } else {
	ALIHLTERRORGUARD(1, "incomplete Configure event, missing parameter data block");
      }
    }
    if (fEventType==gkAliEventTypeReadPreprocessor) {
      if (indexUpdtDCSEvent>=0) {
      TString modules;
      if (fpInputBlocks[indexUpdtDCSEvent].fPtr!=NULL && fpInputBlocks[indexUpdtDCSEvent].fSize>0) {
	modules.Append(reinterpret_cast<const char*>(fpInputBlocks[indexUpdtDCSEvent].fPtr), fpInputBlocks[indexUpdtDCSEvent].fSize);
      }
      HLTDebug("received preprocessor update command: detectors %s", modules.IsNull()?"ALL":modules.Data());
      int tmpResult=ReadPreprocessorValues(modules[0]==0?"ALL":modules.Data());
      if (tmpResult<0) {
	HLTWarning("preprocessor update of component %p (%s) failed with error code %d", this, GetComponentID(), tmpResult);
      }
      } else {
	ALIHLTERRORGUARD(1, "incomplete ReadPreprocessor event, missing parameter data block");
      }
    }
  } else {
    // processing function needs to be called if there are no input data
    // blocks in order to make data source components working.
    bSkipDataProcessing&=~skipModeDefault;
  }

  // data processing is not skipped if the component explicitly asks
  // for the private blocks
  if ((fFlags&kRequireSteeringBlocks)!=0) bSkipDataProcessing=0;

  // data processing is not skipped for data sources
  if (GetComponentType()==AliHLTComponent::kSource) bSkipDataProcessing=0;

  if (fpCTPData) {
    // set the active triggers for this event
    fpCTPData->SetTriggers(trigData);
    // increment CTP trigger counters if available
    if (IsDataEvent()) fpCTPData->Increment(trigData);
  }
  
  // Check if the event processing should be skipped because of the
  // down scaling from the event modulo argument. Using a prime number
  // as pre divisor to pseudo-randomise the event number to get a more
  // uniform distribution.
  if (fEventModulo > 1)
  {
    bSkipDataProcessing |= ( ((AliHLTUInt64_t(fCurrentEvent) / AliHLTUInt64_t(4789)) % AliHLTUInt64_t(fEventModulo)) == 0 );
  }

  AliHLTComponentBlockDataList blockData;
  if (iResult>=0 && !bSkipDataProcessing)
  { // dont delete, sets the scope for the stopwatch guard
    // do not use ALIHLTCOMPONENT_DA_STOPWATCH(); macro
    // in order to avoid 'shadowed variable' warning
    AliHLTStopwatchGuard swguard2(fpStopwatches!=NULL?reinterpret_cast<TStopwatch*>(fpStopwatches->At((int)kSWDA)):NULL);
    AliHLTMisc::AliOnlineGuard onlineGuard;
    iResult=DoProcessing(evtData, blocks, trigData, outputPtr, size, blockData, edd);
  } // end of the scope of the stopwatch guard
  if (iResult>=0 && !bSkipDataProcessing) {
    if (fOutputBlocks.size()>0) {
      // High Level interface

      //HLTDebug("got %d block(s) via high level interface", fOutputBlocks.size());      
      // sync memory files and descriptors
      AliHLTMemoryFilePList::iterator element=fMemFiles.begin();
      int i=0;
      while (element!=fMemFiles.end() && iResult>=0) {
	if (*element) {
	  if ((*element)->IsClosed()==0) {
	    HLTWarning("memory file has not been closed, force flush");
	    iResult=CloseMemoryFile(*element);
	  }
	}
	element++; i++;
      }

      if (iResult>=0) {
	// create the descriptor list
	if (blockData.size()>0) {
	  HLTError("low level and high interface must not be mixed; use PushBack methods to insert data blocks");
	  iResult=-EFAULT;
	} else {
	  if (compStats.size()>0 && IsDataEvent()) {
	    int offset=AddComponentStatistics(fOutputBlocks, fpOutputBuffer, fOutputBufferSize, fOutputBufferFilled, compStats);
	    if (offset>0) fOutputBufferFilled+=offset;
	  }
	  if (bAddComponentTableEntry) {
	    int offset=AddComponentTableEntry(fOutputBlocks, fpOutputBuffer, fOutputBufferSize, fOutputBufferFilled, parentComponentTables, processingLevel);
	    if (offset>0) size+=offset;
	  }
	  if (forwardedBlocks.size()>0) {
	    fOutputBlocks.insert(fOutputBlocks.end(), forwardedBlocks.begin(), forwardedBlocks.end());
	  }
	  iResult=MakeOutputDataBlockList(fOutputBlocks, &outputBlockCnt, &outputBlocks);
	  size=fOutputBufferFilled;
	}
      }
    } else {
      // Low Level interface
      if (blockData.empty() && size!=0) {
	// set the size to zero because there is no announced data
	//HLTWarning("no output blocks added by component but output buffer filled %d of %d", size, fOutputBufferSize);
	size=0;
      }

      if (compStats.size()>0) {
	int offset=AddComponentStatistics(blockData, fpOutputBuffer, fOutputBufferSize, size, compStats);
	if (offset>0) size+=offset;
      }
      if (bAddComponentTableEntry) {
	int offset=AddComponentTableEntry(blockData, fpOutputBuffer, fOutputBufferSize, size, parentComponentTables, processingLevel);
	if (offset>0) size+=offset;
      }
      if (forwardedBlocks.size()>0) {
	blockData.insert(blockData.end(), forwardedBlocks.begin(), forwardedBlocks.end());
      }
      iResult=MakeOutputDataBlockList(blockData, &outputBlockCnt, &outputBlocks);
    }
    if (iResult<0) {
      HLTFatal("component %s (%p): can not convert output block descriptor list", GetComponentID(), this);
    }
  }
  if (iResult<0 || bSkipDataProcessing) {
    outputBlockCnt=0;
    outputBlocks=NULL;
  }
  CleanupInputObjects();
  if (iResult>=0 && IsDataEvent()) {
    IncrementEventCounter();
  }
  if (outputBlockCnt==0) {
    // no output blocks, set size to 0
    size=0;
  }

  // reset the internal EventData struct
  FillEventData(fCurrentEventData);

  // reset the active triggers
  if (fpCTPData) fpCTPData->SetTriggers(0);

  // set the time for the pushback period
  if (fPushbackPeriod>0) {
    // suppress the output
    TDatime time;
    if (fLastPushBackTime<0) {
      // choose a random offset at beginning to equalize traffic for multiple instances
      // of the component
      gRandom->SetSeed(fChainIdCrc);
      fLastPushBackTime=time.Get();
      fLastPushBackTime-=gRandom->Integer(fPushbackPeriod);
    } else if ((int)time.Get()-fLastPushBackTime>=fPushbackPeriod) {
      fLastPushBackTime=time.Get();
    }
  }

  return iResult;
}

int  AliHLTComponent::AddComponentStatistics(AliHLTComponentBlockDataList& blocks, 
					     AliHLTUInt8_t* buffer,
					     AliHLTUInt32_t bufferSize,
					     AliHLTUInt32_t offset,
					     AliHLTComponentStatisticsList& stats) const
{
  // see header file for function documentation
  int iResult=0;
  if ((fFlags&kDisableComponentStat)!=0) return 0;
#if defined(HLT_COMPONENT_STATISTICS)
  if (stats.size()==0) return -ENOENT;
  // check if there is space for at least one entry
  if (offset+sizeof(AliHLTComponentStatistics)>bufferSize) return 0;
  stats[0].fTotalOutputSize=offset;
  stats[0].fOutputBlockCount=blocks.size();
  if (fpBenchmark) {
    fpBenchmark->Stop();
    stats[0].fTime=(AliHLTUInt32_t)(fpBenchmark->RealTime()*ALIHLTCOMPONENT_STATTIME_SCALER);
    stats[0].fCTime=(AliHLTUInt32_t)(fpBenchmark->CpuTime()*ALIHLTCOMPONENT_STATTIME_SCALER);
    fpBenchmark->Continue();
  }

  sort(stats.begin(), stats.end(), SortComponentStatisticsDescendingByLevel);

  // shrink the number of entries if the buffer is too small
  if (offset+stats.size()*sizeof(AliHLTComponentStatistics)>bufferSize) {
    unsigned originalSize=stats.size();
    AliHLTUInt32_t removedLevel=0;
    do {
      // remove all entries of the level of the last entry
      removedLevel=stats.back().fLevel;
      AliHLTComponentStatisticsList::iterator element=stats.begin();
      element++;
      while (element!=stats.end()) {
	if (element->fLevel<=removedLevel) {
	  element=stats.erase(element);
	} else {
	  element++;
	}
      }
    } while (stats.size()>1 && 
	     (offset+stats.size()*sizeof(AliHLTComponentStatistics)>bufferSize));
    HLTWarning("too little space in output buffer to add block of %d statistics entries (size %d), available %d, removed %d entries",
	       originalSize, sizeof(AliHLTComponentStatistics), bufferSize-offset, originalSize-stats.size());
  } else {
    HLTDebug("adding block of %d statistics entries", stats.size());
  }
  assert(stats.size()>0);
  if (stats.size()==0) return 0;

  if (offset+stats.size()*sizeof(AliHLTComponentStatistics)<=bufferSize) {
    AliHLTComponentBlockData bd;
    FillBlockData( bd );
    bd.fOffset        = offset;
    bd.fSize          = stats.size()*sizeof(AliHLTComponentStatistics);
    bd.fDataType      = kAliHLTDataTypeComponentStatistics;
    bd.fSpecification = fChainIdCrc;
    memcpy(buffer+offset, &(stats[0]), bd.fSize);
    blocks.push_back(bd);
    iResult=bd.fSize;
  }
#else
  if (blocks.size() && buffer && bufferSize && offset && stats.size()) {
    // get rid of warning
  }
#endif // HLT_COMPONENT_STATISTICS
  return iResult;
}

int  AliHLTComponent::AddComponentTableEntry(AliHLTComponentBlockDataList& blocks, 
					     AliHLTUInt8_t* buffer,
					     AliHLTUInt32_t bufferSize,
					     AliHLTUInt32_t offset,
					     const vector<AliHLTUInt32_t>& parents,
					     int processingLevel) const
{
  // see header file for function documentation
  int iResult=0;
  if ((fFlags&kDisableComponentStat)!=0) return 0;
#if defined(HLT_COMPONENT_STATISTICS)
  // the payload consists of the AliHLTComponentTableEntry struct,
  // followed by a an array of 32bit crc chain ids and the component
  // description string
  unsigned int payloadSize=sizeof(AliHLTComponentTableEntry);
  payloadSize+=parents.size()*sizeof(AliHLTUInt32_t);

  // the component description has the following format:
  // chain-id{component-id:arguments}
  const char* componentId=const_cast<AliHLTComponent*>(this)->GetComponentID();
  unsigned int descriptionSize=fChainId.size()+1;
  descriptionSize+=2; // the '{}' around the component id
  descriptionSize+=strlen(componentId);
  descriptionSize+=1; // the ':' between component id and arguments
  descriptionSize+=fComponentArgs.size();

  payloadSize+=descriptionSize;
  if (buffer && (offset+payloadSize<=bufferSize)) {
    AliHLTUInt8_t* pTgt=buffer+offset;
    memset(pTgt, 0, payloadSize);

    // write entry
    AliHLTComponentTableEntry* pEntry=reinterpret_cast<AliHLTComponentTableEntry*>(pTgt);
    pEntry->fStructSize=sizeof(AliHLTComponentTableEntry);
    pEntry->fLevel=processingLevel>=0?processingLevel:0;
    pEntry->fNofParents=parents.size();
    pEntry->fSizeDescription=descriptionSize;
    pTgt=pEntry->fBuffer;

    // write array of parents
    if (parents.size()>0) {
      unsigned int copy=parents.size()*sizeof(vector<AliHLTUInt32_t>::value_type);
      memcpy(pTgt, &parents[0], parents.size()*sizeof(vector<AliHLTUInt32_t>::value_type));
      pTgt+=copy;
    }

    // write component description
    memcpy(pTgt, fChainId.c_str(), fChainId.size());
    pTgt+=fChainId.size();
    *pTgt++='{';
    memcpy(pTgt, componentId, strlen(componentId));
    pTgt+=strlen(componentId);
    *pTgt++=':';
    memcpy(pTgt, fComponentArgs.c_str(), fComponentArgs.size());
    pTgt+=fComponentArgs.size();
    *pTgt++='}';
    *pTgt++=0;

    AliHLTComponentBlockData bd;
    FillBlockData( bd );
    bd.fOffset        = offset;
    bd.fSize          = payloadSize;
    bd.fDataType      = kAliHLTDataTypeComponentTable;
    bd.fSpecification = fChainIdCrc;
    blocks.push_back(bd);
    iResult=bd.fSize;
  }
#else
  if (blocks.size() && buffer && bufferSize && offset && parents.size() && processingLevel) {
    // get rid of warning
  }
#endif // HLT_COMPONENT_STATISTICS
  return iResult;
}

AliHLTComponent::AliHLTStopwatchGuard::AliHLTStopwatchGuard()
  :
  fpStopwatch(NULL),
  fpPrec(NULL)
{
  // standard constructor (not for use)
}

AliHLTComponent::AliHLTStopwatchGuard::AliHLTStopwatchGuard(TStopwatch* pStopwatch)
  :
  fpStopwatch(pStopwatch),
  fpPrec(NULL)
{
  // constructor

  // check for already existing guard
  if (fgpCurrent) fpPrec=fgpCurrent;
  fgpCurrent=this;

  // stop the preceeding guard if it controls a different stopwatch
  int bStart=1;
  if (fpPrec && fpPrec!=this) bStart=fpPrec->Hold(fpStopwatch);

  // start the stopwatch if the current guard controls a different one
  if (fpStopwatch && bStart==1) fpStopwatch->Start(kFALSE);
}

AliHLTComponent::AliHLTStopwatchGuard::AliHLTStopwatchGuard(const AliHLTStopwatchGuard&)
  :
  fpStopwatch(NULL),
  fpPrec(NULL)
{
  //
  // copy constructor not for use
  //
}

AliHLTComponent::AliHLTStopwatchGuard& AliHLTComponent::AliHLTStopwatchGuard::operator=(const AliHLTStopwatchGuard&)
{
  //
  // assignment operator not for use
  //
  fpStopwatch=NULL;
  fpPrec=NULL;
  return *this;
}

AliHLTComponent::AliHLTStopwatchGuard* AliHLTComponent::AliHLTStopwatchGuard::fgpCurrent=NULL;

AliHLTComponent::AliHLTStopwatchGuard::~AliHLTStopwatchGuard()
{
  // destructor

  // resume the preceeding guard if it controls a different stopwatch
  int bStop=1;
  if (fpPrec && fpPrec!=this) bStop=fpPrec->Resume(fpStopwatch);

  // stop the stopwatch if the current guard controls a different one
  if (fpStopwatch && bStop==1) fpStopwatch->Stop();

  // resume to the preceeding guard
  fgpCurrent=fpPrec;
}

int AliHLTComponent::AliHLTStopwatchGuard::Hold(const TStopwatch* pSucc)
{
  // see header file for function documentation
  if (fpStopwatch!=NULL && fpStopwatch!=pSucc) fpStopwatch->Stop();
  return fpStopwatch!=pSucc?1:0;
}

int AliHLTComponent::AliHLTStopwatchGuard::Resume(const TStopwatch* pSucc)
{
  // see header file for function documentation
  if (fpStopwatch!=NULL && fpStopwatch!=pSucc) fpStopwatch->Start(kFALSE);
  return fpStopwatch!=pSucc?1:0;
}

int AliHLTComponent::SetStopwatch(TObject* pSW, AliHLTStopwatchType type) 
{
  // see header file for function documentation
  int iResult=0;
  if (pSW!=NULL && type<kSWTypeCount) {
    if (fpStopwatches) {
      TObject* pObj=fpStopwatches->At((int)type);
      if (pSW==NULL        // explicit reset
	  || pObj==NULL) { // explicit set
	fpStopwatches->AddAt(pSW, (int)type);
      } else if (pObj!=pSW) {
	HLTWarning("stopwatch %d already set, reset first", (int)type);
	iResult=-EBUSY;
      }
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTComponent::SetStopwatches(TObjArray* pStopwatches)
{
  // see header file for function documentation
  if (pStopwatches==NULL) return -EINVAL;

  int iResult=0;
  for (int i=0 ; i<(int)kSWTypeCount && pStopwatches->GetEntriesFast(); i++)
    SetStopwatch(pStopwatches->At(i), (AliHLTStopwatchType)i);
  return iResult;
}

AliHLTUInt32_t AliHLTComponent::GetRunNo() const
{
  // see header file for function documentation

  // 2010-02-11 OCDB is now the reliable source for the run number, it is
  // initialized either by aliroot or the external interface depending
  // on the environment. It turned out that the rundescriptor is not set
  // in the aliroot mode, resulting in an invalid run number. However this
  // did not cause problems until now. OCDB initialization has been revised
  // already in 10/2009. This was a remnant.
  // Have to check whether we get rid of the rundescriptor at some point.
  if (fpRunDesc==NULL) return AliHLTMisc::Instance().GetCDBRunNo();
  if (fpRunDesc->fRunNo!=(unsigned)AliHLTMisc::Instance().GetCDBRunNo()) {
    HLTWarning("run number mismatch: ocdb %d   run descriptor %d", AliHLTMisc::Instance().GetCDBRunNo(), fpRunDesc->fRunNo);
  }
  return fpRunDesc->fRunNo;
}

AliHLTUInt32_t AliHLTComponent::GetRunType() const
{
  // see header file for function documentation
  if (fpRunDesc==NULL) return kAliHLTVoidRunType;
  return fpRunDesc->fRunType;
}

AliHLTUInt32_t    AliHLTComponent::GetTimeStamp() const
{
  // see header file for function documentation
  if (fCurrentEventData.fEventCreation_s) {
    return  fCurrentEventData.fEventCreation_s;
  }
  // using the actual UTC if the time stamp was not set by the framework
  return static_cast<AliHLTUInt32_t>(time(NULL));
}

AliHLTUInt32_t    AliHLTComponent::GetPeriodNumber() const
{
  // see header file for function documentation
  return (GetEventId()>>36)&0xfffffff;
}

AliHLTUInt32_t    AliHLTComponent::GetOrbitNumber() const
{
  // see header file for function documentation
  return (GetEventId()>>12)&0xffffff;
}

AliHLTUInt16_t    AliHLTComponent::GetBunchCrossNumber() const
{
  // see header file for function documentation
  return GetEventId()&0xfff;
}

bool AliHLTComponent::IsDataEvent(AliHLTUInt32_t* pTgt) const
{
  // see header file for function documentation
  if (pTgt) *pTgt=fEventType;
  return (fEventType==gkAliEventTypeData ||
	  fEventType==gkAliEventTypeDataReplay);
}

int AliHLTComponent::CopyStruct(void* pStruct, unsigned int iStructSize, unsigned int iBlockNo,
				const char* structname, const char* eventname)
{
  // see header file for function documentation
  int iResult=0;
  if (pStruct!=NULL && iStructSize>sizeof(AliHLTUInt32_t)) {
    if (fpInputBlocks!=NULL && iBlockNo<fCurrentEventData.fBlockCnt) {
      AliHLTUInt32_t* pTgt=(AliHLTUInt32_t*)pStruct;
      if (fpInputBlocks[iBlockNo].fPtr && fpInputBlocks[iBlockNo].fSize) {
	AliHLTUInt32_t copy=*((AliHLTUInt32_t*)fpInputBlocks[iBlockNo].fPtr);
	if (fpInputBlocks[iBlockNo].fSize!=copy) {
	  HLTWarning("%s event: mismatch of block size (%d) and structure size (%d)", eventname, fpInputBlocks[iBlockNo].fSize, copy);
	  if (copy>fpInputBlocks[iBlockNo].fSize) copy=fpInputBlocks[iBlockNo].fSize;
	}
	if (copy!=iStructSize) {
	  HLTWarning("%s event: mismatch in %s version (data type version %d)", eventname, structname, ALIHLT_DATA_TYPES_VERSION);
	  if (copy>iStructSize) {
	    copy=iStructSize;
	  } else {
	    memset(pTgt, 0, iStructSize);
	  }
	}
	memcpy(pTgt, fpInputBlocks[iBlockNo].fPtr, copy);
	*pTgt=iStructSize;
	iResult=copy;
      } else {
	HLTWarning("%s event: missing data block", eventname);
      }
    } else {
      iResult=-ENODATA;
    }
  } else {
    HLTError("invalid struct");
    iResult=-EINVAL;
  }
  return iResult;
}

AliHLTUInt32_t AliHLTComponent::CalculateChecksum(const AliHLTUInt8_t* buffer, int size)
{
  // see header file for function documentation
  AliHLTUInt32_t  remainder = 0; 
  const AliHLTUInt8_t crcwidth=(8*sizeof(AliHLTUInt32_t));
  const AliHLTUInt32_t topbit=1 << (crcwidth-1);
  const AliHLTUInt32_t polynomial=0xD8;  /* 11011 followed by 0's */

  // code from
  // http://www.netrino.com/Embedded-Systems/How-To/CRC-Calculation-C-Code

  /*
   * Perform modulo-2 division, a byte at a time.
   */
  for (int byte = 0; byte < size; ++byte)
    {
      /*
       * Bring the next byte into the remainder.
       */
      remainder ^= (buffer[byte] << (crcwidth - 8));

      /*
       * Perform modulo-2 division, a bit at a time.
       */
      for (uint8_t bit = 8; bit > 0; --bit)
        {
	  /*
	   * Try to divide the current data bit.
	   */
	  if (remainder & topbit)
            {
	      remainder = (remainder << 1) ^ polynomial;
            }
	  else
            {
	      remainder = (remainder << 1);
            }
        }
    }

  /*
   * The final remainder is the CRC result.
   */
  return (remainder);
}

int AliHLTComponent::ExtractComponentTableEntry(const AliHLTUInt8_t* pBuffer, AliHLTUInt32_t size,
						string& retChainId, string& retCompId, string& retCompArgs,
						vector<AliHLTUInt32_t>& parents, int& level)
{
  // see header file for function documentation
  retChainId.clear();
  retCompId.clear();
  retCompArgs.clear();
  parents.clear();
  level=-1;
  if (!pBuffer || size==0) return 0;

  const AliHLTComponentTableEntry* pEntry=reinterpret_cast<const AliHLTComponentTableEntry*>(pBuffer);
  if (size<8/* the initial size of the structure*/ ||
      pEntry==NULL || pEntry->fStructSize<8) return -ENOMSG;

  if (pEntry->fStructSize!=sizeof(AliHLTComponentTableEntry)) return -EBADF;
  level=pEntry->fLevel;
  const AliHLTUInt32_t* pParents=reinterpret_cast<const AliHLTUInt32_t*>(pEntry->fBuffer);
  const AliHLTUInt8_t* pEnd=pBuffer+size;

  if (pParents+pEntry->fNofParents>=reinterpret_cast<const AliHLTUInt32_t*>(pEnd)) return -ENODEV;
  for (unsigned int i=0; i<pEntry->fNofParents; i++, pParents++) {
    parents.push_back(*pParents);
  }

  const char* pDescription=reinterpret_cast<const char*>(pParents);
  if (pDescription+pEntry->fSizeDescription>=reinterpret_cast<const char*>(pEnd) ||
      *(pDescription+pEntry->fSizeDescription)!=0) {
    return -EBADF;
  }

  TString descriptor=reinterpret_cast<const char*>(pDescription);
  TString chainId;
  TString compId;
  TString compArgs;
  TObjArray* pTokens=descriptor.Tokenize("{");
  if (pTokens) {
    int n=0;
    if (pTokens->GetEntriesFast()>n) {
      retChainId=pTokens->At(n++)->GetName();
    }
    if (pTokens->GetEntriesFast()>n) {
      compId=pTokens->At(n++)->GetName();
    }
    delete pTokens;
  }
  if (!compId.IsNull() && (pTokens=compId.Tokenize(":"))!=NULL) {
    int n=0;
    if (pTokens->GetEntriesFast()>n) {
      compId=pTokens->At(n++)->GetName();
    }
    if (pTokens->GetEntriesFast()>n) {
      compArgs=pTokens->At(n++)->GetName();
    }
    delete pTokens;
  }
  compId.ReplaceAll("}", "");
  compArgs.ReplaceAll("}", "");

  retCompId=compId;
  retCompArgs=compArgs;

  if (retChainId.size()==0) return -ENODATA;

  return 1;
}

int AliHLTComponent::ExtractTriggerData(
    const AliHLTComponentTriggerData& trigData,
    const AliHLTUInt8_t (**attributes)[gkAliHLTBlockDAttributeCount],
    AliHLTUInt64_t* status,
    const AliRawDataHeader** cdh,
    AliHLTReadoutList* readoutlist,
    bool printErrors
  )
{
  // see header file for function documentation
  
  // Check that the trigger data structure is the correct size.
  if (trigData.fStructSize != sizeof(AliHLTComponentTriggerData))
  {
    if (printErrors)
    {
      AliHLTLogging log;
      log.LoggingVarargs(kHLTLogError, Class_Name(), FUNCTIONNAME(), __FILE__, __LINE__,
          "Invalid trigger structure size: %d but expected %d.", trigData.fStructSize, sizeof(AliHLTComponentTriggerData)
        );
    }
    return -ENOENT;
  }
  
  // Get the size of the AliHLTEventTriggerData structure without the readout list part.
  // The way we do this here should also handle memory alignment correctly.
  AliHLTEventTriggerData* dummy = NULL;
  size_t sizeWithoutReadout = (char*)(&dummy->fReadoutList) - (char*)(dummy);
  
  // Check that the trigger data pointer points to data of a size we can handle.
  // Either it is the size of AliHLTEventTriggerData or the size of the old
  // version of AliHLTEventTriggerData using AliHLTEventDDLV0.
  if (trigData.fDataSize != sizeof(AliHLTEventTriggerData) and
      trigData.fDataSize != sizeWithoutReadout + sizeof(AliHLTEventDDLV0)
     )
  {
    if (printErrors)
    {
      AliHLTLogging log;
      log.LoggingVarargs(kHLTLogError, Class_Name(), FUNCTIONNAME(), __FILE__, __LINE__,
          "Invalid trigger data size: %d but expected %d.", trigData.fDataSize, sizeof(AliHLTEventTriggerData)
        );
    }
    return -EBADF;
  }
  
  AliHLTEventTriggerData* evtData = reinterpret_cast<AliHLTEventTriggerData*>(trigData.fData);
  assert(evtData != NULL);
  
  // Check that the CDH has 8 words.
  if (cdh != NULL and evtData->fCommonHeaderWordCnt != 8)
  {
    if (printErrors)
    {
      AliHLTLogging log;
      log.LoggingVarargs(kHLTLogError, Class_Name(), FUNCTIONNAME(), __FILE__, __LINE__,
          "Common Data Header (CDH) has wrong number of data words: %d but expected %d",
          evtData->fCommonHeaderWordCnt, sizeof(AliRawDataHeader)/sizeof(AliHLTUInt32_t)
        );
    }
    return -EBADMSG;
  }
  
  // Check that the readout list has the correct count of words. i.e. something we can handle,
  if (readoutlist != NULL and
      evtData->fReadoutList.fCount != (unsigned)gkAliHLTDDLListSizeV0 and
      evtData->fReadoutList.fCount != (unsigned)gkAliHLTDDLListSizeV1
     )
  {
    if (printErrors)
    {
      AliHLTLogging log;
      log.LoggingVarargs(kHLTLogError, Class_Name(), FUNCTIONNAME(), __FILE__, __LINE__,
          "Readout list structure has wrong number of data words: %d but expected %d",
          evtData->fReadoutList.fCount, gkAliHLTDDLListSize
        );
    }
    return -EPROTO;
  }
  
  if (attributes != NULL)
  {
    *attributes = &evtData->fAttributes;
  }
  if (status != NULL)
  {
    *status = evtData->fHLTStatus;
  }
  if (cdh != NULL)
  {
    const AliRawDataHeader* cdhptr = reinterpret_cast<const AliRawDataHeader*>(&evtData->fCommonHeader[0]);
    *cdh = cdhptr;
  }
  if (readoutlist != NULL)
  {
    *readoutlist = AliHLTReadoutList(evtData->fReadoutList);
  }
  return 0;
}

AliHLTUInt32_t AliHLTComponent::ExtractEventTypeFromCDH(const AliRawDataHeader* cdh)
{
  // see header file for function documentation
  
  UChar_t l1msg = cdh->GetL1TriggerMessage();
  if ((l1msg & 0x1) == 0x0) return gkAliEventTypeData;
  // The L2SwC bit must be one if we got here, i.e. l1msg & 0x1 == 0x1.
  if (((l1msg >> 2) & 0xF) == 0xE) return gkAliEventTypeStartOfRun;
  if (((l1msg >> 2) & 0xF) == 0xF) return gkAliEventTypeEndOfRun;
  // Check the C1T bit to see if this is a calibration event,
  // if not then it must be some other software trigger event.
  if (((l1msg >> 6) & 0x1) == 0x1) return gkAliEventTypeCalibration;
  return gkAliEventTypeSoftware;
}

int AliHLTComponent::LoggingVarargs(AliHLTComponentLogSeverity severity, 
				    const char* originClass, const char* originFunc,
				    const char* file, int line, ... ) const
{
  // see header file for function documentation
  int iResult=0;

  va_list args;
  va_start(args, line);

  // logging function needs to be const in order to be called from const member functions
  // without problems. But at this point we face the problem with virtual members which
  // are not necessarily const.
  AliHLTComponent* nonconst=const_cast<AliHLTComponent*>(this);
  AliHLTLogging::SetLogString(this, ", %p", "%s (%s_pfmt_): ", 
			      fChainId[0]!=0?fChainId.c_str():nonconst->GetComponentID(),
			      nonconst->GetComponentID());
  iResult=SendMessage(severity, originClass, originFunc, file, line, AliHLTLogging::BuildLogString(NULL, args, true /*append*/));
  va_end(args);

  return iResult;
}

TUUID AliHLTComponent::GenerateGUID()
{
  // Generates a globally unique identifier.

  // Start by creating a new UUID. We cannot use the one automatically generated
  // by ROOT because the algorithm used will not guarantee unique IDs when generating
  // these UUIDs at a high rate in parallel.
  TUUID uuid;
  // We then use the generated UUID to form part of the random number seeds which
  // will be used to generate a proper random UUID. For good measure we use a MD5
  // hash also. Note that we want to use the TUUID class because it will combine the
  // host address information into the UUID. Using gSystem->GetHostByName() apparently
  // can cause problems on Windows machines with a firewall, because it always tries
  // to contact a DNS. The TUUID class handles this case appropriately.
  union
  {
    UChar_t buf[16];
    UShort_t word[8];
    UInt_t dword[4];
  };
  uuid.GetUUID(buf);
  TMD5 md5;
  md5.Update(buf, sizeof(buf));
  TMD5 md52 = md5;
  md5.Final(buf);
  dword[0] += gSystem->GetUid();
  dword[1] += gSystem->GetGid();
  dword[2] += gSystem->GetPid();
  for (int i = 0; i < 4; ++i)
  {
    gRandom->SetSeed(dword[i]);
    dword[i] = gRandom->Integer(0xFFFFFFFF);
  }
  md52.Update(buf, sizeof(buf));
  md52.Final(buf);
  // To keep to the standard we need to set the version and reserved bits.
  word[3] = (word[3] & 0x0FFF) | 0x4000;
  buf[8] = (buf[8] & 0x3F) | 0x80;

  // Create the name of the new class and file.
  char uuidstr[64];
  sprintf(uuidstr, "%08x-%04x-%04x-%02x%02x-%02x%02x%02x%02x%02x%02x",
    dword[0], word[2], word[3], buf[8], buf[9], buf[10], buf[11], buf[12], buf[13], buf[14], buf[15]
  );

  uuid.SetUUID(uuidstr);
  return uuid;
}

int AliHLTComponent::ScanECSParam(const char* ecsParam)
{
  // see header file for function documentation

  // format of the parameter string from ECS
  // <command>;<parameterkey>=<parametervalue>;<parameterkey>=<parametervalue>;...
  // search for a subset of the parameterkeys
  //   RUN_TYPE=
  //   RUN_NUMBER=
  //   HLT_IN_DDL_LIST=
  //   CTP_TRIGGER_CLASS=
  //   DATA_FORMAT_VERSION=
  //   BEAM_TYPE=
  //   HLT_OUT_DDL_LIST=
  //   HLT_TRIGGER_CODE=
  //   DETECTOR_LIST=
  //   HLT_MODE=
  // The command apears not to be sent by the online framework
  int iResult=0;
  TString string=ecsParam;
  TObjArray* parameter=string.Tokenize(";");
  if (parameter) {
    for (int i=0; i<parameter->GetEntriesFast(); i++) {
      TString entry=parameter->At(i)->GetName();
      HLTDebug("scanning ECS entry: %s", entry.Data());
      TObjArray* entryParams=entry.Tokenize("=");
      if (entryParams) {
	if (entryParams->GetEntriesFast()>1) {
	  if ((((TObjString*)entryParams->At(0))->GetString()).CompareTo("CTP_TRIGGER_CLASS")==0) {
	    int result=InitCTPTriggerClasses(entryParams->At(1)->GetName());
	    if (iResult>=0 && result<0) iResult=result;
	  } else {
	    // TODO: scan the other parameters
	    // e.g. consistency check of run number
	  }
	}
	delete entryParams;
      }
    }
    delete parameter;
  }

  return iResult;
}

int AliHLTComponent::SetupCTPData()
{
  // see header file for function documentation
  if (fpCTPData) delete fpCTPData;
  fpCTPData=new AliHLTCTPData;
  if (!fpCTPData) return -ENOMEM;
  return 0;
}

int AliHLTComponent::InitCTPTriggerClasses(const char* ctpString)
{
  // see header file for function documentation
  if (!fpCTPData) return 0; // silently accept as the component has to announce that it want's the CTP info
  return fpCTPData->InitCTPTriggerClasses(ctpString);
}

bool AliHLTComponent::EvaluateCTPTriggerClass(const char* expression, AliHLTComponentTriggerData& trigData) const
{
  // see header file for function documentation
  if (!fpCTPData) {
    static bool bWarningThrown=false;
    if (!bWarningThrown) HLTError("Trigger classes not initialized, use SetupCTPData from DoInit()");
    bWarningThrown=true;
    return false;
  }

  return fpCTPData->EvaluateCTPTriggerClass(expression, trigData);
}

int AliHLTComponent::CheckCTPTrigger(const char* name) const
{
  // see header file for function documentation
  if (!fpCTPData) {
    static bool bWarningThrown=false;
    if (!bWarningThrown) HLTError("Trigger classes not initialized, use SetupCTPData from DoInit()");
    bWarningThrown=true;
    return false;
  }

  return fpCTPData->CheckTrigger(name);
}

Double_t AliHLTComponent::GetBz()
{
  // Returns Bz.
  return AliHLTMisc::Instance().GetBz();
}

Double_t AliHLTComponent::GetBz(const Double_t *r)
{
  // Returns Bz (kG) at the point "r" .
  return AliHLTMisc::Instance().GetBz(r);
}

void AliHLTComponent::GetBxByBz(const Double_t r[3], Double_t b[3])
{
  // Returns Bx, By and Bz (kG) at the point "r" .
  AliHLTMisc::Instance().GetBxByBz(r, b);
}
