// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

/** @file   AliHLTControlTask.cxx
    @author Matthias Richter
    @date   
    @brief  Special task to produce the control events.
*/

#include "AliHLTControlTask.h"
#include "AliHLTComponentHandler.h"
#include <cassert>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTControlTask)

AliHLTControlTask::AliHLTControlTask()
  : AliHLTTask()
  , fBlocks()
  , fpData(NULL)
  , fSize(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTControlTask::~AliHLTControlTask()
{
  // see header file for class documentation
  ResetBlocks();
}

int AliHLTControlTask::CreateComponent(AliHLTConfiguration* /*pConf*/, AliHLTComponentHandler* pCH, AliHLTComponent*& pComponent) const
{
  // see header file for class documentation
  int iResult=0;
  if ((pComponent=new AliHLTControlEventComponent(this))) {
    const AliHLTAnalysisEnvironment* pEnv=pCH->GetEnvironment();
    const char* argv[]={
      "-disable-component-stat"
    };
    int argc=sizeof(argv)/sizeof(const char*);
    if ((iResult=pComponent->Init(pEnv, NULL, argc, argv))>=0) {
      //HLTDebug("component %s (%p) created", pComponent->GetComponentID(), pComponent); 
    } else {
      HLTError("Initialization of component \"%s\" failed with error %d", pComponent->GetComponentID(), iResult);
    }
    return iResult;
  }
  return -ENOMEM;
}

void AliHLTControlTask::SetBlocks(const AliHLTComponentBlockDataList& list)
{
  // see header file for class documentation
  fBlocks.assign(list.begin(), list.end());
  AliHLTComponentBlockDataList::iterator element=fBlocks.begin();
  for (;element!=fBlocks.end(); element++) fSize+=element->fSize;

  // allocate buffer for the payload of all blocks
  fpData=new AliHLTUInt8_t[fSize];
  AliHLTUInt8_t offset=0;

  // copy and redirect
  for (element=fBlocks.begin();element!=fBlocks.end(); element++) {
    memcpy(fpData+offset, element->fPtr, element->fSize);
    element->fPtr=fpData+offset;
    offset+=element->fSize;
  }
}

void AliHLTControlTask::ResetBlocks()
{
  // see header file for class documentation
  fBlocks.clear();
  if (fpData) delete [] fpData;
  fpData=NULL;
  fSize=0;
}

AliHLTControlTask::AliHLTControlEventComponent::AliHLTControlEventComponent(const AliHLTControlTask* pParent)
  :
  fpParent(pParent)
{
  // see header file for class documentation
  assert(pParent);
}

AliHLTControlTask::AliHLTControlEventComponent::~AliHLTControlEventComponent()
{
  // see header file for class documentation
}

AliHLTComponentDataType AliHLTControlTask::AliHLTControlEventComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTMultipleDataType;
}

int AliHLTControlTask::AliHLTControlEventComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  // see header file for class documentation
  tgtList.clear();
  tgtList.push_back(kAliHLTDataTypeSOR);
  tgtList.push_back(kAliHLTDataTypeEOR);
  return tgtList.size();
}

void AliHLTControlTask::AliHLTControlEventComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  if (fpParent && fpParent->fSize>0) constBase=fpParent->fSize;
  else constBase=0;
  inputMultiplier=0;
}

int AliHLTControlTask::AliHLTControlEventComponent::GetEvent(const AliHLTComponentEventData& /*evtData*/,
							     AliHLTComponentTriggerData& /*trigData*/,
							     AliHLTUInt8_t* outputPtr, 
							     AliHLTUInt32_t& size,
							     vector<AliHLTComponentBlockData>& outputBlocks )
{
  // see header file for class documentation
  if (!fpParent) return -ENODEV;
  const AliHLTControlTask* pParent=fpParent;

  AliHLTUInt32_t capacity=size;
  size=0;
  if (capacity<pParent->fSize) {
    return -ENOSPC;
  }

  // return if no event has been set
  if (pParent->fpData==NULL ||
      pParent->fBlocks.size()==0) {
    //HLTInfo("no control event to send");
    return 0;
  }

  for (unsigned int i=0; i<pParent->fBlocks.size(); i++) {
    HLTDebug("publishing control block %s", DataType2Text(pParent->fBlocks[i].fDataType).c_str());
    memcpy(outputPtr+size, pParent->fBlocks[i].fPtr, pParent->fBlocks[i].fSize);
    AliHLTComponentBlockData bd;
    FillBlockData(bd);
    bd.fOffset=size;
    bd.fSize=pParent->fBlocks[i].fSize;
    bd.fDataType=pParent->fBlocks[i].fDataType;
    bd.fSpecification=pParent->fBlocks[i].fSpecification;
    outputBlocks.push_back( bd );
    size+=bd.fSize;
  }

  return size;
}
