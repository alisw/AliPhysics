// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

/** @file   testAliHLT_C_Component_WrapperInterface.C
    @author Matthias Richter
    @date   
    @brief  Test program for the old wrapper interface, frozen since Jul 08
 */

#include "AliHLTProcessor.h"
#include "AliHLTModuleAgent.h"
#include "AliHLT_C_Component_WrapperInterface.h"


const char* gDummy="dummy";
AliHLTUInt32_t gRunNo=kAliHLTVoidRunNo;
const char* gChainId="<void>";
AliHLTComponentDataType gInputDt=kAliHLTVoidDataType;

class TestProcessor : public AliHLTProcessor
{
public:
  TestProcessor();
  ~TestProcessor();

  const char* GetComponentID() {return "TestProcessor";}

  void GetInputDataTypes( vector<AliHLTComponentDataType>& list)
  {list.push_back(kAliHLTAnyDataType);}

  AliHLTComponentDataType GetOutputDataType() {return kAliHLTAnyDataType;}

  void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
  {constBase=100; inputMultiplier=2.0;}
  
  AliHLTComponent* Spawn() {return new TestProcessor;}
private:
  int DoInit( int argc, const char** argv ) {
    if (fState!=kCreated) {
      HLTError("wrong state (%d): component already initialized", fState);
      return -EBUSY;
    }

    if (argc>0 && argv) {
    }
    gRunNo=GetRunNo();
    gChainId=GetChainId();

    fState=kInitialized;
    return 0;
  }
  int DoDeinit() {
    if (fState!=kInitialized) {
      HLTError("wrong state (%d): required %d kInitialized", fState, kInitialized);
      return -ENODEV;
    }

    return 0;
  }

  int DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/) {
    int iResult=0;
    if (fState!=kInitialized) {
      HLTError("wrong state (%d): component not initialized", fState);
      return -ENODEV;
    }

    HLTInfo("processing event");
    for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock();
	 pBlock!=NULL; 
	 pBlock=GetNextInputBlock()) {
      gInputDt=pBlock->fDataType;
      AliHLTComponentDataType dt;
      SetDataType(dt, "-OUTPUT-", "MYCO");
      iResult=PushBack(pBlock->fPtr, pBlock->fSize/2, dt, ~pBlock->fSpecification);
      Forward();
    }

    return iResult;
  }

  enum {
    kCreated =0,
    kInitialized,
    kProcessing,
  };

  int fState; //!transient
};

TestProcessor::TestProcessor()
  : 
  AliHLTProcessor(),
  fState(kCreated)
 
{
}

TestProcessor::~TestProcessor()
{
}

class TestAgent : public AliHLTModuleAgent
{
public:
  TestAgent() : AliHLTModuleAgent("TEST") {}
  ~TestAgent() {}

  int RegisterComponents(AliHLTComponentHandler* pHandler) const
  {
    pHandler->AddComponent(new TestProcessor);
    return 0;
  }
  
};

TestAgent gAgent;

int Logging( void* /*param*/, 
	     AliHLTComponentLogSeverity severity,
	     const char* origin,
	     const char* keyword,
	     const char* message)
{
  cout << "Logging: "<< severity << " " << origin << " " << keyword << " " << message << endl;
  return 0;
}

void* AllocMemory( void* param, unsigned long size )
{
  if (param!=&gDummy) {
    cerr << "AllocMemoryFunc callback with wrong parameter " << endl;
    abort();
  }
  if (size==0) return NULL;
  return new AliHLTUInt8_t[size];
}

const char* GetComponentDescription(void* param)
{
  if (param!=&gDummy) {
    cerr << "GetComponentDescription callback with wrong parameter " << endl;
    abort();
  }
  return "-chainid=test";
}

int main(int /*argc*/, const char** /*argv*/)
{
  int iResult=0;
  AliHLTComponentEnvironment environment;
  memset(&environment, 0, sizeof(environment));
  environment.fStructSize=sizeof(environment);
  environment.fAllocMemoryFunc=AllocMemory;
  environment.fLoggingFunc=Logging;
  //environment.fGetComponentDescription=GetComponentDescription;
  if ((iResult=AliHLT_C_Component_InitSystem( &environment ))<0) {
    cerr << "error: AliHLT_C_Component_InitSystem failed with " << iResult << endl;
    return iResult;
  }

  if ((iResult=AliHLT_C_Component_LoadLibrary("../util/.libs/libAliHLTUtil.so"))<0) {
    cerr << "error: AliHLT_C_Component_LoadLibrary failed with " << iResult << endl;
    return iResult;
  }

  AliHLTRunDesc desc;
  memset(&desc, 0, sizeof(desc));
  desc.fStructSize=sizeof(desc);
  desc.fRunNo=0xbeef;
  if ((iResult=AliHLT_C_SetRunDescription(&desc, "dummy run"))<0) {
    cerr << "error: AliHLT_C_Component_SetRunDescription failed with " << iResult << endl;
    return iResult;
  }

  AliHLTComponentHandle handle;
  if ((iResult=AliHLT_C_CreateComponent("TestProcessor", &gDummy, 0, NULL, &handle ))<0) {
    cerr << "error: AliHLT_C_Component_CreateComponent failed with " << iResult << endl;
    return iResult;
  }

  if (gRunNo!=0xbeef) {
    cerr << "error: propagation of run number failed " << hex << gRunNo << " vs. 0xbeef" << endl;
    return -1;
  }

  // can be used in the new interface again
//   if (strcmp(gChainId, "test")) {
//     cerr << "propagation of chain id failed: '" << gChainId << "' vs. test" << endl;
//     return -1;
//   }


  const char* inputData="some data to be copied";
  AliHLTComponentBlockData inputBlock;
  memset(&inputBlock, 0, sizeof(inputBlock));
  inputBlock.fStructSize=sizeof(inputBlock);
  inputBlock.fPtr=(void*)inputData;
  inputBlock.fSize=strlen(inputData);
  inputBlock.fDataType=kAliHLTDataTypeDDLRaw;
  inputBlock.fSpecification=0xdead;

  AliHLTComponentEventData evtData;
  memset(&evtData, 0, sizeof(evtData));
  evtData.fStructSize=sizeof(evtData);
  evtData.fBlockCnt=1;

  AliHLTComponentTriggerData trigData;
  memset(&trigData, 0, sizeof(trigData));
  trigData.fStructSize=sizeof(trigData);

  AliHLTUInt8_t outputPtr[100];
  AliHLTUInt32_t size=sizeof(outputPtr);

  AliHLTUInt32_t outputBlockCnt=0;
  AliHLTComponentBlockData* outputBlocks=NULL;
  AliHLTComponentEventDoneData* edd=NULL;

  if ((iResult=AliHLT_C_ProcessEvent( handle, &evtData, &inputBlock, &trigData, outputPtr,
				      &size, &outputBlockCnt, &outputBlocks, &edd ))<0) {
    cerr << "error: AliHLT_C_Component_ProcessEvent failed with " << iResult << endl;
    return iResult;
  }

  if (outputBlockCnt<2) {
    cerr << "error: mismatch in output block count, expecting >2, got " << outputBlockCnt << endl;
    return -1;
  }

  if (outputBlocks==NULL) {
    cerr << "error: did not get output block array " << endl;
    return -1;
  }

  AliHLTComponentDataType outdt;
  AliHLTComponent::SetDataType(outdt, "-OUTPUT-", "MYCO");
  bool bHaveForwarded=false;
  bool bHaveCopied=false;
  for (unsigned int i=0; i<outputBlockCnt; i++) {
    if (outputBlocks[i].fDataType==kAliHLTDataTypeDDLRaw) {
      if (outputBlocks[i].fPtr!=inputData ||
	  outputBlocks[i].fSize!=strlen(inputData)) {
	cerr << "error: failed comparing forwarded input block" << endl;
	return -1;
      }
      bHaveForwarded=true;
    } else if (outputBlocks[i].fDataType==outdt) {
      if (outputBlocks[i].fSize!=strlen(inputData)/2) {
	cerr << "error: wrong size of copied block" << endl;
	return -1;
      }
      if (memcmp(inputData, outputPtr+outputBlocks[i].fOffset, outputBlocks[i].fSize)) {
	cerr << "error: failed comparing copied block" << endl;
	return -1;
      }
      bHaveCopied=true;
    }
  }

  if (!bHaveForwarded) {
    cerr << "error: did not get forwarded data block" << endl;
    return -1;
  }

  if (!bHaveCopied) {
    cerr << "error: did not get copied data block" << endl;
    return -1;
  }

  AliHLT_C_DestroyComponent(handle);

  if ((iResult=AliHLT_C_Component_DeinitSystem( ))<0) {
    cerr << "AliHLT_C_Component_DeinitSystem failed with " << iResult << endl;
    return iResult;
  }

  return 0;
}
