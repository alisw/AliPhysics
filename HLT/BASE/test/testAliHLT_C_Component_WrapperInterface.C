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
    @brief  Test program for the external wrapper interface
 */

#include "AliHLTProcessor.h"
#include "AliHLTModuleAgent.h"
#include "AliHLT_C_Component_WrapperInterface.h"


const char* gDummy="dummy";
AliHLTUInt32_t gRunNo=kAliHLTVoidRunNo;
const char* gChainId="<void>";

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
    if (argc>0 && argv) {
    }
    gRunNo=GetRunNo();
    gChainId=GetChainId();
    return 0;
  }
  int DoDeinit() {
    return 0;
  }

  int DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData) {
    return 0;
  }
};

TestProcessor::TestProcessor()
  : AliHLTProcessor() 
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

int Logging( void* param, 
	     AliHLTComponentLogSeverity severity,
	     const char* origin,
	     const char* keyword,
	     const char* message)
{
  cout << "Logging: "<< severity << " " << origin << " " << keyword << " " << message << endl;
  return 0;
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
  environment.fLoggingFunc=Logging;
  environment.fGetComponentDescription=GetComponentDescription;
  if ((iResult=AliHLT_C_Component_InitSystem( &environment ))<0) {
    cerr << "AliHLT_C_Component_InitSystem failed with " << iResult << endl;
    return iResult;
  }

  if ((iResult=AliHLT_C_Component_LoadLibrary("libAliHLTUtil.so"))<0) {
    cerr << "AliHLT_C_Component_LoadLibrary failed with " << iResult << endl;
    return iResult;
  }

  AliHLTRunDesc desc;
  memset(&desc, 0, sizeof(desc));
  desc.fStructSize=sizeof(desc);
  desc.fRunNo=0xbeef;
  if ((iResult=AliHLT_C_SetRunDescription(&desc, "dummy run"))<0) {
    cerr << "AliHLT_C_Component_SetRunDescription failed with " << iResult << endl;
    return iResult;
  }

  AliHLTComponentHandle handle;
  if ((iResult=AliHLT_C_CreateComponent("TestProcessor", &gDummy, 0, NULL, &handle ))<0) {
    cerr << "AliHLT_C_Component_CreateComponent failed with " << iResult << endl;
    return iResult;
  }

  if (gRunNo!=0xbeef) {
    cerr << "propagation of run number failed " << hex << gRunNo << " vs. 0xbeef" << endl;
    return -1;
  }

  if (strcmp(gChainId, "test")) {
    cerr << "propagation of chain id failed: '" << gChainId << "' vs. test" << endl;
    return -1;
  }

  if ((iResult=AliHLT_C_Component_DeinitSystem( ))<0) {
    cerr << "AliHLT_C_Component_DeinitSystem failed with " << iResult << endl;
    return iResult;
  }

  return 0;
}
