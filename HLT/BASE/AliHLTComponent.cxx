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
#include "AliHLTSystem.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTComponent)

AliHLTComponent::AliHLTComponent()
  :
  fEnvironment(),
  fCurrentEvent(0)
{ 
  memset(&fEnvironment, 0, sizeof(AliHLTComponentEnvironment));
  if (fpComponentHandler)
    fpComponentHandler->ScheduleRegister(this);
}

AliHLTComponent::~AliHLTComponent()
{
}

AliHLTComponentHandler* AliHLTComponent::fpComponentHandler=NULL;

int AliHLTComponent::SetGlobalComponentHandler(AliHLTComponentHandler* pCH, int bOverwrite) 
{
  int iResult=0;
  if (fpComponentHandler==NULL || bOverwrite!=0)
    fpComponentHandler=pCH;
  else
    iResult=-EPERM;
  return iResult;
}

int AliHLTComponent::UnsetGlobalComponentHandler() {
  return SetGlobalComponentHandler(NULL,1);
}

int AliHLTComponent::Init( AliHLTComponentEnvironment* environ, void* environ_param, int argc, const char** argv )
{
  int iResult=0;
  if (environ) {
    memcpy(&fEnvironment, environ, sizeof(AliHLTComponentEnvironment));
    fEnvironment.fParam=environ_param;
  }
  iResult=DoInit(argc, argv);
  return iResult;
}

int AliHLTComponent::Deinit()
{
  int iResult=0;
  iResult=DoDeinit();
  return iResult;
}

int AliHLTComponent::DoInit( int argc, const char** argv )
{
  if (argc==0 && argv==NULL) {
    // this is currently just to get rid of the warning "unused parameter"
  }
  return 0;
}

int AliHLTComponent::DoDeinit()
{
  return 0;
}

void AliHLTComponent::DataType2Text( const AliHLTComponentDataType& type, char output[14] ) {
memset( output, 0, 14 );
strncat( output, type.fOrigin, 4 );
strcat( output, ":" );
strncat( output, type.fID, 8 );
}

void* AliHLTComponent::AllocMemory( unsigned long size ) {
  if (fEnvironment.fAllocMemoryFunc)
    return (*fEnvironment.fAllocMemoryFunc)(fEnvironment.fParam, size );
  return NULL;
}

int AliHLTComponent::MakeOutputDataBlockList( const vector<AliHLTComponentBlockData>& blocks, AliHLTUInt32_t* blockCount,
					      AliHLTComponentBlockData** outputBlocks ) {
    if ( !blockCount || !outputBlocks )
	return EFAULT;
    AliHLTUInt32_t count = blocks.size();
    if ( !count )
	{
	*blockCount = 0;
	*outputBlocks = NULL;
	return 0;
	}
    *outputBlocks = reinterpret_cast<AliHLTComponentBlockData*>( AllocMemory( sizeof(AliHLTComponentBlockData)*count ) );
    if ( !*outputBlocks )
	return ENOMEM;
    for ( unsigned long i = 0; i < count; i++ )
	(*outputBlocks)[i] = blocks[i];
    *blockCount = count;
    return 0;

}

int AliHLTComponent::GetEventDoneData( unsigned long size, AliHLTComponentEventDoneData** edd ) {
  if (fEnvironment.fGetEventDoneDataFunc)
    return (*fEnvironment.fGetEventDoneDataFunc)(fEnvironment.fParam, fCurrentEvent, size, edd );
  return -ENOSYS;
}

int AliHLTComponent::FindMatchingDataTypes(AliHLTComponent* pConsumer, vector<AliHLTComponentDataType>* tgtList) 
{
  int iResult=0;
  if (pConsumer) {
    vector<AliHLTComponentDataType> ctlist;
    ((AliHLTComponent*)pConsumer)->GetInputDataTypes(ctlist);
    vector<AliHLTComponentDataType>::iterator type=ctlist.begin();
    while (type!=ctlist.end() && iResult==0) {
      if ((*type)==GetOutputDataType()) {
	if (tgtList) tgtList->push_back(*type);
	iResult++;
	break;
      }
      type++;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}
