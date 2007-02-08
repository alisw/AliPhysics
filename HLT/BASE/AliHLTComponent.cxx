// $Id$

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Matthias Richter <Matthias.Richter@ift.uib.no>                *
 *          Timm Steinbeck <timm@kip.uni-heidelberg.de>                   *
 *          for The ALICE Off-line Project.                               *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTComponent.cxx
    @author Matthias Richter, Timm Steinbeck
    @date   
    @brief  Base class implementation for HLT components. */

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTStdIncludes.h"
#include "AliHLTComponent.h"
#include "AliHLTComponentHandler.h"
#include "TString.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTComponent)

AliHLTComponent::AliHLTComponent()
  :
  fEnvironment(),
  fCurrentEvent(0),
  fEventCount(-1)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  memset(&fEnvironment, 0, sizeof(AliHLTComponentEnvironment));
  if (fgpComponentHandler)
    fgpComponentHandler->ScheduleRegister(this);
}

AliHLTComponent::AliHLTComponent(const AliHLTComponent&)
  :
  fEnvironment(),
  fCurrentEvent(0),
  fEventCount(-1)
{
  // see header file for class documentation
  HLTFatal("copy constructor untested");
}

AliHLTComponent& AliHLTComponent::operator=(const AliHLTComponent&)
{ 
  // see header file for class documentation
  HLTFatal("assignment operator untested");
  return *this;
}

AliHLTComponent::~AliHLTComponent()
{
  // see header file for function documentation
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

int AliHLTComponent::Init( AliHLTComponentEnvironment* environ, void* environParam, int argc, const char** argv )
{
  // see header file for function documentation
  int iResult=0;
  if (environ) {
    memcpy(&fEnvironment, environ, sizeof(AliHLTComponentEnvironment));
    fEnvironment.fParam=environParam;
  }
  iResult=DoInit(argc, argv);
  if (iResult>=0) fEventCount=0;
  return iResult;
}

int AliHLTComponent::Deinit()
{
  // see header file for function documentation
  int iResult=0;
  iResult=DoDeinit();
  return iResult;
}

int AliHLTComponent::DoInit( int argc, const char** argv )
{
  // see header file for function documentation
  if (argc==0 && argv==NULL) {
    // this is currently just to get rid of the warning "unused parameter"
  }
  return 0;
}

int AliHLTComponent::DoDeinit()
{
  // see header file for function documentation
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

string AliHLTComponent::DataType2Text( const AliHLTComponentDataType& type )
{
  // see header file for function documentation
  string out("");
  
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

int AliHLTComponent::MakeOutputDataBlockList( const vector<AliHLTComponentBlockData>& blocks, AliHLTUInt32_t* blockCount,
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
	if (blocks[i].fDataType==kAliHLTAnyDataType) {
	  memset((*outputBlocks)[i].fDataType.fID, '*', kAliHLTComponentDataTypefIDsize);
	  memset((*outputBlocks)[i].fDataType.fOrigin, '*', kAliHLTComponentDataTypefOriginSize);
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

int AliHLTComponent::FindMatchingDataTypes(AliHLTComponent* pConsumer, vector<AliHLTComponentDataType>* tgtList) 
{
  // see header file for function documentation
  int iResult=0;
  if (pConsumer) {
    vector<AliHLTComponentDataType> ctlist;
    ((AliHLTComponent*)pConsumer)->GetInputDataTypes(ctlist);
    vector<AliHLTComponentDataType>::iterator type=ctlist.begin();
    while (type!=ctlist.end() && iResult==0) {
      if ((*type)==GetOutputDataType() ||
	  (*type)==kAliHLTAnyDataType) {
	if (tgtList) tgtList->push_back(*type);
	iResult++;
	// this loop has to be changed in case of multiple output types
	break;
      }
      type++;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

void AliHLTComponent::FillBlockData( AliHLTComponentBlockData& blockData ) const
{
  // see header file for function documentation
  blockData.fStructSize = sizeof(blockData);
  FillShmData( blockData.fShmKey );
  blockData.fOffset = ~(AliHLTUInt32_t)0;
  blockData.fPtr = NULL;
  blockData.fSize = 0;
  FillDataType( blockData.fDataType );
  blockData.fSpecification = ~(AliHLTUInt32_t)0;
}

void AliHLTComponent::FillShmData( AliHLTComponentShmData& shmData ) const
{
  // see header file for function documentation
  shmData.fStructSize = sizeof(shmData);
  shmData.fShmType = gkAliHLTComponentInvalidShmType;
  shmData.fShmID = gkAliHLTComponentInvalidShmID;
}

void AliHLTComponent::FillDataType( AliHLTComponentDataType& dataType ) const
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
  tgtdt.fStructSize = sizeof(AliHLTComponentDataType);
  memset(&tgtdt.fID[0], 0, kAliHLTComponentDataTypefIDsize);
  memset(&tgtdt.fOrigin[0], 0, kAliHLTComponentDataTypefOriginSize);

  if ((int)strlen(id)>kAliHLTComponentDataTypefIDsize) {
    HLTWarning("data type id %s is too long, truncated to %d", id, kAliHLTComponentDataTypefIDsize);
  }
  strncpy(&tgtdt.fID[0], id, kAliHLTComponentDataTypefIDsize);

  if ((int)strlen(origin)>kAliHLTComponentDataTypefOriginSize) {
    HLTWarning("data type origin %s is too long, truncated to %d", origin, kAliHLTComponentDataTypefOriginSize);
  }
  strncpy(&tgtdt.fOrigin[0], origin, kAliHLTComponentDataTypefOriginSize);
}

void AliHLTComponent::FillEventData(AliHLTComponentEventData& evtData)
{
  // see header file for function documentation
  memset(&evtData, 0, sizeof(AliHLTComponentEventData));
  evtData.fStructSize=sizeof(AliHLTComponentEventData);
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
  int iResult=0;
  fCurrentEvent=evtData.fEventID;
  iResult=DoProcessing(evtData, blocks, trigData, outputPtr, size, outputBlockCnt, outputBlocks, edd);
  IncrementEventCounter();
  return iResult;
}
