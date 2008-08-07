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

/** @file   AliHLTOUTTask.cxx
    @author Matthias Richter
    @date   
    @brief  A special HLTOUT sibling working as a data sink in chains
*/

#include "AliHLTOUTTask.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOUTTask)

AliHLTOUTTask::AliHLTOUTTask(const char* chains)
  :
  AliHLTOUT(),
  AliHLTTask(),
  fpDummyTask(NULL),
  fpDummyConfiguration(NULL),
  fBlockDescList()
{ 
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  TString taskname=chains;
  taskname.ReplaceAll(" ", "_");
  taskname+="_hltouttask";

  // This is just a trick to use the existing Task handling, especially the
  // ProcessTask function.
  // 1. The 'BlockFilter' just forwards the data blocks of all specified chains
  // 2. A dummy task is added to the target list in order to set the data segments
  //    after forwarding. In case of an empty target list the data is just discarded.
  // That's maybe not the most elegant solution but far less awkward as one would
  // expect.
  fpConfiguration=new AliHLTConfiguration(taskname.Data(), "BlockFilter", chains, NULL);
  TString dummyname=chains;
  dummyname.ReplaceAll(" ", "_");
  dummyname+="_never_used_dummy_";
  fpDummyConfiguration=new AliHLTConfiguration(dummyname.Data(), "BlockFilter", taskname.Data(), NULL);
  if (fpDummyConfiguration) {
    fpDummyTask=new AliHLTTask(fpDummyConfiguration);
    SetTarget(fpDummyTask);
  }
}

AliHLTOUTTask::~AliHLTOUTTask()
{
  // see header file for class documentation

  // UnsetTarget called automatically
  if (fpDummyTask) delete fpDummyTask;
  fpDummyTask=NULL;

  if (fpDummyConfiguration) delete fpDummyConfiguration;
  fpDummyConfiguration=NULL;

  if (fpConfiguration) delete fpConfiguration;
  fpConfiguration=NULL;
}

int AliHLTOUTTask::CustomInit(AliHLTComponentHandler* pCH)
{
  // see header file for class documentation
  if (!fpDummyTask) return -ENOENT;
  return fpDummyTask->Init(NULL, pCH);
}

int AliHLTOUTTask::CustomCleanup()
{
  // see header file for class documentation
  if (!fpDummyTask) return -ENOENT;
  return fpDummyTask->Deinit();
}

int AliHLTOUTTask::GenerateIndex()
{
  // see header file for class documentation
  int iResult=0;
  if (fpDataBuffer) {
    if (!fpDataBuffer->FindConsumer(fpDummyTask->GetComponent())) {
      // in order to subscribe to the buffers the dummy consumer
      // needs to be set. The dummy task is not in the AliHLTSystem chain
      // and therefor not automatically set.
      fpDataBuffer->SetConsumer(fpDummyTask->GetComponent());
    }
    fBlockDescList.clear();
    if (fpDataBuffer->GetNofSegments()>0) {
      if ((iResult=fpDataBuffer->Subscribe(fpDummyTask->GetComponent(), fBlockDescList))>=0) {
	for (unsigned int i=0; i<fBlockDescList.size(); i++) {
	  AliHLTOUTBlockDescriptor desc(fBlockDescList[i].fDataType, fBlockDescList[i].fSpecification, i, this);
	  //HLTDebug("adding block %d: %s %#x", i, AliHLTComponent::DataType2Text(fBlockDescList[i].fDataType).c_str(), fBlockDescList[i].fSpecification);
	  iResult=AddBlockDescriptor(desc);
	}
      } else {
	HLTError("failed to subscribe to data buffer");
      }
    }
  } else {
    // 2008-08-07 this is not a failure condition
    // If the chain has not been processed because LocalReconstruction
    // is not enabled, the task will be empty
    //HLTWarning("no data buffer available");
  }
  return iResult;
}

int AliHLTOUTTask::GetDataBuffer(AliHLTUInt32_t index, const AliHLTUInt8_t* &pBuffer, 
				 AliHLTUInt32_t& size)
{
  // see header file for class documentation
  int iResult=0;
  if (index>=fBlockDescList.size()) return -ENOENT;
  pBuffer=reinterpret_cast<AliHLTUInt8_t*>(fBlockDescList[index].fPtr);
  size=fBlockDescList[index].fSize;
  return iResult;
}

AliHLTOUT::AliHLTOUTByteOrder AliHLTOUTTask::CheckBlockByteOrder(AliHLTUInt32_t /*index*/)
{
  // see header file for class documentation
  return kInvalidByteOrder;
}

int AliHLTOUTTask::CheckBlockAlignment(AliHLTUInt32_t /*index*/, AliHLTOUT::AliHLTOUTDataType /*type*/)
{
  // see header file for class documentation
  int iResult=0;
  return iResult;
}

int AliHLTOUTTask::ResetInput()
{
  // see header file for class documentation
  int iResult=0;
  if (!fpDataBuffer) return 0;

  if (fBlockDescList.size()==0 && fpDataBuffer->GetNofPendingConsumers()>0) {
    // not yet subscribed
    fpDataBuffer->Subscribe(fpDummyTask->GetComponent(), fBlockDescList);
  }

  for (unsigned int i=0; i<fBlockDescList.size(); i++) {
    fpDataBuffer->Release(&(fBlockDescList[i]), fpDummyTask->GetComponent(), this);
  }
  fBlockDescList.clear();
  return iResult;
}

const char* AliHLTOUTTask::GetSourceChains() const
{
  // see header file for class documentation
  if (!fpConfiguration) return "";
  return fpConfiguration->GetSourceSettings();
}
