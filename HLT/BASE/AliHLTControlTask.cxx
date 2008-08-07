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
  :
  fEvent(kAliHLTVoidDataType),
  fSpecification(kAliHLTVoidDataSpec),
  fpData(NULL),
  fSize(0)
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
}

int AliHLTControlTask::CreateComponent(AliHLTConfiguration* /*pConf*/, AliHLTComponentHandler* pCH, AliHLTComponent*& pComponent) const
{
  // see header file for class documentation
  int iResult=0;
  if ((pComponent=new AliHLTControlEventComponent(this))) {
    const AliHLTAnalysisEnvironment* pEnv=pCH->GetEnvironment();
    if ((iResult=pComponent->Init(pEnv, NULL, 0, NULL))>=0) {
      //HLTDebug("component %s (%p) created", pComponent->GetComponentID(), pComponent); 
    } else {
      HLTError("Initialization of component \"%s\" failed with error %d", pComponent->GetComponentID(), iResult);
    }
    return iResult;
  }
  return -ENOMEM;
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
  else constBase=sizeof(AliHLTRunDesc);
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
  if (size<pParent->fSize) {
    return -ENOSPC;
  }

  if (pParent->fpData && pParent->fSize) {
    memcpy(outputPtr, pParent->fpData, pParent->fSize);
  }

  // return if no event has been set
  if (pParent->fEvent==kAliHLTVoidDataType) {
    HLTInfo("no control event to send");
    return 0;
  }

  HLTInfo("publishing control event %s", DataType2Text(pParent->fEvent).c_str());
  AliHLTComponentBlockData bd;
  FillBlockData(bd);
  bd.fOffset=0;
  bd.fSize=pParent->fSize;
  bd.fDataType=pParent->fEvent;
  bd.fSpecification=pParent->fSpecification;
  outputBlocks.push_back( bd );

  return bd.fSize;
}
