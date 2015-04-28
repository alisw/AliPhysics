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

/** @file   AliHLTDumpTask.cxx
    @author Matthias Richter
    @date   
    @brief  Base class for data sinks with direct buffer access.
*/

#include "AliHLTDumpTask.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTDumpTask)

AliHLTDumpTask::AliHLTDumpTask(const char* chains)
  :
  AliHLTTask(),
  fpDummyTask(NULL),
  fpDummyConfiguration(NULL),
  fBlocks()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  if (chains && chains[0]!=0) SetChains(chains);
}

AliHLTDumpTask::~AliHLTDumpTask()
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

int AliHLTDumpTask::SetChains(const char* chains)
{
  // see header file for class documentation
  if (!chains || chains[0]==0) {
    HLTError("invalid chain id string");
    return -EINVAL;
  }
  TString taskname=chains;
  taskname.ReplaceAll(" ", "_");
  taskname+="_hltDumpTask";

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
    fpDummyTask->SetDependency(this);
  }
  return 0;
}

int AliHLTDumpTask::CustomInit(AliHLTComponentHandler* pCH)
{
  // see header file for class documentation
  if (!fpDummyTask) return -ENOENT;
  return fpDummyTask->Init(NULL, pCH);
}

int AliHLTDumpTask::CustomCleanup()
{
  // see header file for class documentation
  if (!fpDummyTask) return -ENOENT;
  return fpDummyTask->Deinit();
}

const char* AliHLTDumpTask::GetSourceChains() const
{
  // see header file for class documentation
  if (!fpConfiguration) return "";
  return fpConfiguration->GetSourceSettings();
}

const AliHLTComponentBlockDataList& AliHLTDumpTask::GetDataBlocks()
{
  // see header file for class documentation
  if (fBlocks.size()>0) return fBlocks;
  if (fpDataBuffer) {
    if (!fpDataBuffer->FindConsumer(fpDummyTask->GetComponent())) {
      // in order to subscribe to the buffers the dummy consumer
      // needs to be set. The dummy task is not in the AliHLTSystem chain
      // and therefor not automatically set.
      fpDataBuffer->SetConsumer(fpDummyTask->GetComponent());
    }
    fBlocks.clear();
    if (fpDataBuffer->GetNofSegments()>0 
	&& fpDataBuffer->FindConsumer(fpDummyTask->GetComponent(), 0 /*search only among pending consumers*/)>0) {
      if (fpDataBuffer->Subscribe(fpDummyTask->GetComponent(), fBlocks)>=0) {
	return fBlocks;
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
  fBlocks.clear();
  return fBlocks;
}

int AliHLTDumpTask::ReleaseDataBlocks()
{
  // see header file for class documentation
  int iResult=0;
  if (!fpDataBuffer) return 0;

  if (fBlocks.size()==0 && fpDataBuffer->GetNofPendingConsumers()>0) {
    // There are data blocks in the parents which this task has not yet
    // subscribed to. The subscription takes place in GetDataBlocks.
    // However, this method is not necessarily called.
    //
    // In order to switch buffer states correctly, first let the dummy
    // task as the only consumer subscribe to all those buffers and
    // release them further down. Basically, the buffers are arranged
    // in a different internal list, which is the only state they can be
    // released from. This approach has been chosen to implement the
    // DumpTask having no real consumers but at the same time has to
    // behave like a normal task in AliHLTTask::ProcessTask
    fpDataBuffer->Subscribe(fpDummyTask->GetComponent(), fBlocks);
  }

  for (unsigned int i=0; i<fBlocks.size(); i++) {
    fpDataBuffer->Release(&(fBlocks[i]), fpDummyTask->GetComponent(), this);
  }
  fBlocks.clear();
  return iResult;
}
