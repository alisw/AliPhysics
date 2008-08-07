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

/** @file   AliHLTComponent.cxx
    @author Matthias Richter, Timm Steinbeck
    @date   
    @brief  Base class implementation for HLT components. */

#if __GNUC__>= 3
using namespace std;
#endif

//#include "AliHLTStdIncludes.h"
#include "AliHLTComponent.h"
#include "AliHLTComponentHandler.h"
#include "AliHLTMessage.h"
#include "TString.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjectTable.h"
#include "TClass.h"
#include "TStopwatch.h"
#include "AliHLTMemoryFile.h"
#include "AliHLTMisc.h"
#include <cassert>

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
  fpDDLList(NULL),
  fCDBSetRunNoFunc(false),
  fChainId(),
  fChainIdCrc(0),
  fpBenchmark(NULL),
  fRequireSteeringBlocks(false),
  fEventType(gkAliEventTypeUnknown)
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

int AliHLTComponent::Init(const AliHLTAnalysisEnvironment* comenv, void* environParam, int argc, const char** argv )
{
  // see header file for function documentation
  HLTLogKeyword(GetComponentID());
  int iResult=0;
  if (comenv) {
    memset(&fEnvironment, 0, sizeof(AliHLTAnalysisEnvironment));
    memcpy(&fEnvironment, comenv, comenv->fStructSize<sizeof(AliHLTAnalysisEnvironment)?comenv->fStructSize:sizeof(AliHLTAnalysisEnvironment));
    fEnvironment.fStructSize=sizeof(AliHLTAnalysisEnvironment);
    fEnvironment.fParam=environParam;
  }
  const char** pArguments=NULL;
  int iNofChildArgs=0;
  TString argument="";
  int bMissingParam=0;
  if (argc>0) {
    pArguments=new const char*[argc];
    if (pArguments) {
      for (int i=0; i<argc && iResult>=0; i++) {
	argument=argv[i];
	if (argument.IsNull()) continue;

	// benchmark
	if (argument.CompareTo("benchmark")==0) {

	  // loglevel
	} else if (argument.CompareTo("loglevel")==0) {
	  if ((bMissingParam=(++i>=argc))) break;
	  TString parameter(argv[i]);
	  parameter.Remove(TString::kLeading, ' '); // remove all blanks
	  if (parameter.BeginsWith("0x") &&
	      parameter.Replace(0,2,"",0).IsHex()) {
	    AliHLTComponentLogSeverity loglevel=kHLTLogNone;
	    sscanf(parameter.Data(),"%x", (unsigned int*)&loglevel);
	    SetLocalLoggingLevel(loglevel);
	  } else {
	    HLTError("wrong parameter for argument %s, hex number expected", argument.Data());
	    iResult=-EINVAL;
	  }
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
    iResult=DoInit(iNofChildArgs, pArguments);
  }
  if (iResult>=0) {
    fEventCount=0;

    // find out if the component wants to get the steering events
    // explicitly
    AliHLTComponentDataTypeList inputDt;
    GetInputDataTypes(inputDt);
    for (AliHLTComponentDataTypeList::iterator dt=inputDt.begin();
	 dt!=inputDt.end() && !fRequireSteeringBlocks;
	 dt++) {
      fRequireSteeringBlocks|=MatchExactly(*dt,kAliHLTDataTypeSOR);
      fRequireSteeringBlocks|=MatchExactly(*dt,kAliHLTDataTypeRunType);
      fRequireSteeringBlocks|=MatchExactly(*dt,kAliHLTDataTypeEOR);
      fRequireSteeringBlocks|=MatchExactly(*dt,kAliHLTDataTypeDDL);
      fRequireSteeringBlocks|=MatchExactly(*dt,kAliHLTDataTypeComponentStatistics);
    }
  }
  if (pArguments) delete [] pArguments;

#if defined(__DEBUG) || defined(HLT_COMPONENT_STATISTICS)
  // benchmarking stopwatch for the component statistics
  fpBenchmark=new TStopwatch;
  if (fpBenchmark) {
    fpBenchmark->Start();
  }
#endif

  return iResult;
}

int AliHLTComponent::Deinit()
{
  // see header file for function documentation
  HLTLogKeyword(GetComponentID());
  int iResult=0;
  iResult=DoDeinit();
  if (fpRunDesc) {
    HLTWarning("did not receive EOR for run %d", fpRunDesc->fRunNo);
    AliHLTRunDesc* pRunDesc=fpRunDesc;
    fpRunDesc=NULL;
    delete pRunDesc;
  }
  fEventCount=0;
  return iResult;
}

int AliHLTComponent::InitCDB(const char* cdbPath, AliHLTComponentHandler* pHandler)
{
  // see header file for function documentation
  int iResult=0;
  if (pHandler) {
  // I have to think about separating the library handling from the
  // component handler. Requiring the component hanlder here is not
  // the cleanest solution.
  // We presume the library already to be loaded
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
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      TString argument=((TObjString*)pTokens->At(i++))->GetString();
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
  }
  
  return iResult;
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

void AliHLTComponent::DataType2Text( const AliHLTComponentDataType& type, char output[kAliHLTComponentDataTypefIDsize+kAliHLTComponentDataTypefOriginSize+2] ) const
{
  // see header file for function documentation
  memset( output, 0, kAliHLTComponentDataTypefIDsize+kAliHLTComponentDataTypefOriginSize+2 );
  strncat( output, type.fOrigin, kAliHLTComponentDataTypefOriginSize );
  strcat( output, ":" );
  strncat( output, type.fID, kAliHLTComponentDataTypefIDsize );
}

string AliHLTComponent::DataType2Text( const AliHLTComponentDataType& type, int mode)
{
  // see header file for function documentation
  string out("");

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

int AliHLTComponent::GetEventDoneData( unsigned long size, AliHLTComponentEventDoneData** edd ) 
{
  // see header file for function documentation
  if (fEnvironment.fGetEventDoneDataFunc)
    return (*fEnvironment.fGetEventDoneDataFunc)(fEnvironment.fParam, fCurrentEvent, size, edd );
  return -ENOSYS;
}

int AliHLTComponent::FindMatchingDataTypes(AliHLTComponent* pConsumer, AliHLTComponentDataTypeList* tgtList) 
{
  // see header file for function documentation
  int iResult=0;
  if (pConsumer) {
    AliHLTComponentDataTypeList itypes;
    AliHLTComponentDataTypeList otypes;
    otypes.push_back(GetOutputDataType());
    if (otypes[0]==kAliHLTMultipleDataType) {
      otypes.clear();
      int count=0;
      if ((count=GetOutputDataTypes(otypes))>0) {
      } else if (GetComponentType()!=kSink) {
	HLTWarning("component %s indicates multiple output data types but GetOutputDataTypes returns %d", GetComponentID(), count);
      }
    }
    ((AliHLTComponent*)pConsumer)->GetInputDataTypes(itypes);
    AliHLTComponentDataTypeList::iterator otype=otypes.begin();
    for (;otype!=otypes.end();otype++) {
      //PrintDataTypeContent((*otype), "publisher \'%s\'");
      if ((*otype)==(kAliHLTAnyDataType|kAliHLTDataOriginPrivate)) {
	if (tgtList) tgtList->push_back(*otype);
	iResult++;
	continue;
      }
      
      AliHLTComponentDataTypeList::iterator itype=itypes.begin();
      for (;itype!=itypes.end() && (*itype)!=(*otype); itype++);
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
  for (int i=0; i<array->GetEntries(); i++) {
    TObject* pObj=array->At(i);
    // grrr, garbage collection strikes back: When read via AliHLTMessage
    // (CreateInputObject), and written to a TFile afterwards, the
    // TFile::Close calls ROOOT's garbage collection. No clue why the
    // object ended up in the key list and needs to be deleted
    if (pObj && gObjectTable->PtrIsValid(pObj)) delete pObj;
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
    if (fpInputObjects==NULL || (idx=FindInputBlock(pBlock))>=0) {
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

int AliHLTComponent::PushBack(TObject* pObject, const AliHLTComponentDataType& dt, AliHLTUInt32_t spec, 
			      void* pHeader, int headerSize)
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  int iResult=0;
  if (pObject) {
    AliHLTMessage msg(kMESS_OBJECT);
    msg.WriteObject(pObject);
    Int_t iMsgLength=msg.Length();
    if (iMsgLength>0) {
      msg.SetLength(); // sets the length to the first (reserved) word
      iResult=InsertOutputBlock(msg.Buffer(), iMsgLength, dt, spec, pHeader, headerSize);
      if (iResult>=0) {
	HLTDebug("object %s (%p) size %d inserted to output", pObject->ClassName(), pObject, iMsgLength);
      }
    } else {
      HLTError("object serialization failed for object %p", pObject);
      iResult=-ENOMSG;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTComponent::PushBack(TObject* pObject, const char* dtID, const char* dtOrigin, AliHLTUInt32_t spec,
			      void* pHeader, int headerSize)
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  AliHLTComponentDataType dt;
  SetDataType(dt, dtID, dtOrigin);
  return PushBack(pObject, dt, spec, pHeader, headerSize);
}

int AliHLTComponent::PushBack(void* pBuffer, int iSize, const AliHLTComponentDataType& dt, AliHLTUInt32_t spec,
			      void* pHeader, int headerSize)
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  return InsertOutputBlock(pBuffer, iSize, dt, spec, pHeader, headerSize);
}

int AliHLTComponent::PushBack(void* pBuffer, int iSize, const char* dtID, const char* dtOrigin, AliHLTUInt32_t spec,
			      void* pHeader, int headerSize)
{
  // see header file for function documentation
  ALIHLTCOMPONENT_BASE_STOPWATCH();
  AliHLTComponentDataType dt;
  SetDataType(dt, dtID, dtOrigin);
  return PushBack(pBuffer, iSize, dt, spec, pHeader, headerSize);
}

int AliHLTComponent::InsertOutputBlock(void* pBuffer, int iBufferSize, const AliHLTComponentDataType& dt, AliHLTUInt32_t spec,
			      void* pHeader, int iHeaderSize)
{
  // see header file for function documentation
  int iResult=0;
  int iBlkSize = iBufferSize + iHeaderSize;
  if (pBuffer) {
    if (fpOutputBuffer && iBlkSize<=(int)(fOutputBufferSize-fOutputBufferFilled)) {
      AliHLTUInt8_t* pTgt=fpOutputBuffer+fOutputBufferFilled;
      AliHLTComponentBlockData bd;
      FillBlockData( bd );
      bd.fOffset        = fOutputBufferFilled;
      bd.fSize          = iBlkSize;
      bd.fDataType      = dt;
      bd.fSpecification = spec;
      if (pHeader!=NULL && pHeader!=pTgt) {
	memcpy(pTgt, pHeader, iHeaderSize);
      }

      pTgt += (AliHLTUInt8_t) iHeaderSize;

      if (pBuffer!=NULL && pBuffer!=pTgt) {
	memcpy(pTgt, pBuffer, iBufferSize);
	
	//AliHLTUInt32_t firstWord=*((AliHLTUInt32_t*)pBuffer);	
	//HLTDebug("copy %d bytes from %p to output buffer %p, first word %#x", iBufferSize, pBuffer, pTgt, firstWord);
      }
      fOutputBufferFilled+=bd.fSize;
      fOutputBlocks.push_back( bd );
      //HLTDebug("buffer inserted to output: size %d data type %s spec %#x", iBlkSize, DataType2Text(dt).c_str(), spec);
    } else {
      if (fpOutputBuffer) {
	HLTError("too little space in output buffer: %d, required %d", fOutputBufferSize-fOutputBufferFilled, iBlkSize);
      } else {
	HLTError("output buffer not available");
      }
      iResult=-ENOSPC;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTComponent::EstimateObjectSize(TObject* pObject) const
{
  // see header file for function documentation
  if (!pObject) return -EINVAL;
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

int AliHLTComponent::CreateEventDoneData(AliHLTComponentEventDoneData /*edd*/)
{
  // see header file for function documentation
  int iResult=-ENOSYS;
  //#warning  function not yet implemented
  HLTWarning("function not yet implemented");
  return iResult;
}

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
  HLTLogKeyword(GetComponentID());
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
  AliHLTComponentStatisticsList compStats;
#if defined(__DEBUG) || defined(HLT_COMPONENT_STATISTICS)
  AliHLTComponentStatistics outputStat;
  memset(&outputStat, 0, sizeof(AliHLTComponentStatistics));
  outputStat.fId=fChainIdCrc;
  compStats.push_back(outputStat);
  if (fpBenchmark) {
    fpBenchmark->Reset();
    fpBenchmark->Start();
  }
#endif

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
	    compStats.push_back(*pStat);
	  }
	}
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
    }
    if (indexEOREvent>=0) {
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
	  AliHLTRunDesc* pRunDesc=fpRunDesc;
	  fpRunDesc=NULL;
	  delete pRunDesc;
	}
      } else {
	HLTWarning("did not receive SOR, ignoring EOR");
      }
    }
    if (indexComConfEvent>=0 || fEventType==gkAliEventTypeConfiguration) {
      TString cdbEntry;
      if (indexComConfEvent>=0 && fpInputBlocks[indexComConfEvent].fPtr!=NULL && fpInputBlocks[indexComConfEvent].fSize>0) {
	cdbEntry.Append(reinterpret_cast<const char*>(fpInputBlocks[indexComConfEvent].fPtr), fpInputBlocks[indexComConfEvent].fSize);
      }
      HLTDebug("received component configuration command: entry %s", cdbEntry.IsNull()?"none":cdbEntry.Data());
      int tmpResult=Reconfigure(cdbEntry[0]==0?NULL:cdbEntry.Data(), fChainId.c_str());
      if (tmpResult<0) {
	HLTWarning("reconfiguration of component %p (%s) failed with error code %d", this, GetComponentID(), tmpResult);
      }
    }
    if (indexUpdtDCSEvent>=0 || fEventType==gkAliEventTypeReadPreprocessor) {
      TString modules;
      if (fpInputBlocks[indexUpdtDCSEvent].fPtr!=NULL && fpInputBlocks[indexUpdtDCSEvent].fSize>0) {
	modules.Append(reinterpret_cast<const char*>(fpInputBlocks[indexUpdtDCSEvent].fPtr), fpInputBlocks[indexUpdtDCSEvent].fSize);
      }
      HLTDebug("received preprocessor update command: detectors %s", modules.IsNull()?"ALL":modules.Data());
      int tmpResult=ReadPreprocessorValues(modules[0]==0?"ALL":modules.Data());
      if (tmpResult<0) {
	HLTWarning("preprocessor update of component %p (%s) failed with error code %d", this, GetComponentID(), tmpResult);
      }
    }
  } else {
    // processing function needs to be called if there are no input data
    // blocks in order to make data source components working.
    bSkipDataProcessing&=~skipModeDefault;
  }

  // data processing is not skipped if the component explicitly asks
  // for the private blocks
  if (fRequireSteeringBlocks) bSkipDataProcessing=0;

  AliHLTComponentBlockDataList blockData;
  if (iResult>=0 && !bSkipDataProcessing)
  { // dont delete, sets the scope for the stopwatch guard
    // do not use ALIHLTCOMPONENT_DA_STOPWATCH(); macro
    // in order to avoid 'shadowed variable' warning
    AliHLTStopwatchGuard swguard2(fpStopwatches!=NULL?reinterpret_cast<TStopwatch*>(fpStopwatches->At((int)kSWDA)):NULL);
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
	  if (compStats.size()>0) {
	    int offset=AddComponentStatistics(fOutputBlocks, fpOutputBuffer, fOutputBufferSize, fOutputBufferFilled, compStats);
	    if (offset>0) fOutputBufferFilled+=offset;
	  }
	  iResult=MakeOutputDataBlockList(fOutputBlocks, &outputBlockCnt, &outputBlocks);
	  size=fOutputBufferFilled;
	}
      }
    } else {
      // Low Level interface
      if (compStats.size()>0) {
	int offset=AddComponentStatistics(blockData, fpOutputBuffer, fOutputBufferSize, size, compStats);
	if (offset>0) size+=offset;
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
  if (stats.size()==0) return -ENOENT;
  stats[0].fTotalOutputSize=offset;
  stats[0].fOutputBlockCount=blocks.size();
  if (fpBenchmark) {
    stats[0].fTime=(AliHLTUInt32_t)(fpBenchmark->RealTime()*1000000);
    stats[0].fCTime=(AliHLTUInt32_t)(fpBenchmark->CpuTime()*1000000);
  }
  if (offset+stats.size()*sizeof(AliHLTComponentStatistics)<=bufferSize) {
    AliHLTComponentBlockData bd;
    FillBlockData( bd );
    bd.fOffset        = offset;
    bd.fSize          = stats.size()*sizeof(AliHLTComponentStatistics);
    bd.fDataType      = kAliHLTDataTypeComponentStatistics;
    bd.fSpecification = kAliHLTVoidDataSpec;
    unsigned int master=0;
    for (unsigned int i=1; i<blocks.size(); i++) {
      if (blocks[i].fSize>blocks[master].fSize && 
	  !MatchExactly(blocks[i].fDataType, kAliHLTVoidDataType|kAliHLTDataOriginPrivate))
	master=i;
    }
    if (blocks.size()>0 && !MatchExactly(blocks[master].fDataType, kAliHLTVoidDataType|kAliHLTDataOriginPrivate)) {
      // take the data origin of the biggest block as specification
      // this is similar to the treatment in the HOMER interface. For traditional
      // reasons, the bytes are swapped there on a little endian architecture, so
      // we do it as well.
      memcpy(&bd.fSpecification, &blocks[master].fDataType.fOrigin, sizeof(bd.fSpecification));
#ifdef R__BYTESWAP // set on little endian architectures
      bd.fSpecification=((bd.fSpecification & 0xFFULL) << 24) | 
	((bd.fSpecification & 0xFF00ULL) << 8) | 
	((bd.fSpecification & 0xFF0000ULL) >> 8) | 
	((bd.fSpecification & 0xFF000000ULL) >> 24);
#endif
    }
    memcpy(buffer+offset, &(stats[0]), bd.fSize);
    blocks.push_back(bd);
    iResult=bd.fSize;
  }
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

int AliHLTComponent::AliHLTStopwatchGuard::Hold(TStopwatch* pSucc)
{
  // see header file for function documentation
  if (fpStopwatch!=NULL && fpStopwatch!=pSucc) fpStopwatch->Stop();
  return fpStopwatch!=pSucc?1:0;
}

int AliHLTComponent::AliHLTStopwatchGuard::Resume(TStopwatch* pSucc)
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
  for (int i=0 ; i<(int)kSWTypeCount && pStopwatches->GetEntries(); i++)
    SetStopwatch(pStopwatches->At(i), (AliHLTStopwatchType)i);
  return iResult;
}

AliHLTUInt32_t AliHLTComponent::GetRunNo() const
{
  // see header file for function documentation
  if (fpRunDesc==NULL) return kAliHLTVoidRunNo;
  return fpRunDesc->fRunNo;
}

AliHLTUInt32_t AliHLTComponent::GetRunType() const
{
  // see header file for function documentation
  if (fpRunDesc==NULL) return kAliHLTVoidRunType;
  return fpRunDesc->fRunType;
}

bool AliHLTComponent::IsDataEvent(AliHLTUInt32_t* pTgt)
{
  // see header file for function documentation
  if (pTgt) *pTgt=fEventType;
  return (fEventType==gkAliEventTypeData ||
	  fEventType==gkAliEventTypeDataReplay ||
	  fEventType==gkAliEventTypeCalibration);
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

void AliHLTComponent::SetDDLBit(AliHLTEventDDL &list, Int_t ddlId, Bool_t state ) const
{
  // see header file for function documentation
  
  // -- Detector offset
  Int_t ddlIdBase =  TMath::FloorNint( (Double_t) ddlId / 256.0 );
  
  // -- Word Base = 1. word of detector ( TPC has 8 words, TOF 3 ) 
  Int_t wordBase = 0;

  if ( ddlIdBase <= 3 )
    wordBase = ddlIdBase;
  else if ( ddlIdBase > 3 && ddlIdBase < 5 )
    wordBase = ddlIdBase + 7;
  else 
    wordBase = ddlIdBase + 9;

  // -- Bit index in Word
  Int_t bitIdx = ddlId % 32;

  // -- Index of word
  Int_t wordIdx = wordBase;

  // -- if TPC (3) or TOD (5) add word idx
  if ( ( ddlIdBase == 3 ) || ( ddlIdBase == 5 ) ) {
    wordIdx += TMath::FloorNint( (Double_t) ( ddlId - ( ddlIdBase * 256 ) ) / 32.0 );
  }

  // -- Set -- 'OR' word with bit mask;
  if ( state )
    list.fList[wordIdx] |= ( 0x00000001 << bitIdx );
  // -- Unset -- 'AND' word with bit mask;
  else
    list.fList[wordIdx] &= ( 0xFFFFFFFF ^ ( 0x00000001 << bitIdx ) );
}

Int_t AliHLTComponent::GetFirstUsedDDLWord(AliHLTEventDDL &list) const
{
  // see header file for function documentation

  Int_t iResult = -1;

  for ( Int_t wordNdx = 0 ; wordNdx < gkAliHLTDDLListSize ; wordNdx++ ) {

    if ( list.fList[wordNdx] != 0 && iResult == -1 ) {
      // check for special cases TPC and TOF
      if ( wordNdx > 3 && wordNdx <= 10 ) {
	wordNdx = 10;
	iResult = 3;
      }
      else if ( wordNdx > 12 && wordNdx <= 14 ) {
	wordNdx = 14;
	iResult = 12;
      }
      else
	iResult = wordNdx;
    }
    else if ( list.fList[wordNdx] != 0 && iResult >= 0 ) {
      HLTError( "DDLIDs for minimum of TWO detectors ( %d, %d ) set, this function works only for ONE detector.", iResult, wordNdx );
      iResult = -1;
      break;
    }
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
