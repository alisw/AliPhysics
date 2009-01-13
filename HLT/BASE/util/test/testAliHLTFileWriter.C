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

/** @file   testAliHLTFileWriter.C
    @author Matthias Richter
    @date   
    @brief  Test macro/program for the AliHLTFileWriter
 */

#ifndef __CINT__
#include <ostream>
#include <vector>
#include "AliHLTComponent.h"
#include "AliHLTComponentHandler.h"
#include "AliHLTConfiguration.h"
#include "TString.h"
#include "TFile.h"
#include "TSystem.h"
#else
#define EINVAL 22
#define EPROTO 71
#endif //__CINT__

typedef AliHLTUInt32_t TestElement_t;
const char* testDirectory="/tmp/testAliHLTFileWriter";

int SetEvent(unsigned eventNo, vector<TestElement_t>& testData)
{
  unsigned shift=(sizeof(TestElement_t)-1)*8;
  AliHLTUInt32_t mask=0xf<<shift;
  for (unsigned i=0; i<testData.size(); i++) {
    testData[i]&=~mask;
    testData[i]|=eventNo<<shift;
  }
  return 0;
}

int TestCycle(AliHLTComponent* pComponent, const char* arguments, AliHLTComponentBlockDataList blocks,
	      vector<TestElement_t>& testData, unsigned int nofEvents=1)
{
  int iResult=0;
  vector<char*> argList;
  if ((iResult=AliHLTConfiguration::InterpreteString(arguments, argList))<0) {
    cerr << "failed: chopping arguments \"" << arguments << "\"" << endl;
  }

  // notepad for generalized unit test class
  // - test component for outputsize==0 if no blocks are added

  // init FileWriter component
  if ((iResult=pComponent->Init(NULL, NULL, argList.size(), const_cast<const char**>(&(argList[0]))))!=0) {
    cerr << "failed: initializing component" << endl;
    return iResult;
  }

  // process event
  AliHLTComponentBlockData eventTypeBlock;
  AliHLTComponent::FillBlockData(eventTypeBlock);
  eventTypeBlock.fDataType=kAliHLTDataTypeEvent;
  eventTypeBlock.fSpecification=gkAliEventTypeData;
  blocks.push_back(eventTypeBlock);

  AliHLTComponentEventData evtData;
  memset(&evtData, 0, sizeof(AliHLTComponentEventData));
  evtData.fStructSize=sizeof(AliHLTComponentEventData);
  evtData.fEventID=0;
  evtData.fBlockCnt=blocks.size();
  
  AliHLTComponentTriggerData trigData;
  memset(&trigData, 0, sizeof(AliHLTComponentTriggerData));
  trigData.fStructSize=sizeof(AliHLTComponentTriggerData);

  AliHLTUInt32_t size=0;
  AliHLTUInt32_t outputBlockCnt=0;
  AliHLTComponentBlockData* outputBlocks=NULL;
  AliHLTComponentEventDoneData* edd=NULL;

  for (unsigned eventNo=0; eventNo<nofEvents; eventNo++) {
    SetEvent(eventNo, testData);
    evtData.fEventID=eventNo;
    if ((iResult=pComponent->ProcessEvent(evtData, &(blocks[0]), trigData, NULL, size, outputBlockCnt, outputBlocks, edd))!=0) {
      cerr << "failed: event processing" << endl;
      break;
    }
  }
  blocks.pop_back();

  // cleanup FileWriter component
  if ((iResult=pComponent->Deinit())!=0) {
    cerr << "failed: deinit component" << endl;
  }

  // cleanup argument vector
  for (vector<char*>::iterator element=argList.begin();
       element!=argList.end(); element++) {
    delete *element;
  }
  argList.clear();

  return iResult;
}

int BuildFileName(unsigned eventNo, unsigned blockID,
		  const AliHLTComponentDataType& dataType,
		  const AliHLTUInt32_t /*specification*/,
		  bool bConcatenateEvents, bool bConcatenateBlocks,
		  TString& filename)
{
  int iResult=0;
  filename="";

  filename+=testDirectory;
  if (!filename.EndsWith("/")) {
    filename+="/";
  }

  filename+="event";
  if (!bConcatenateEvents) {
    filename+=Form("_0x%08x", eventNo);
  }
  if (!bConcatenateBlocks && !bConcatenateEvents) {
    filename+=Form("_0x%02x", blockID);
    if (dataType!=kAliHLTVoidDataType) {
      filename+="_";
      filename+=AliHLTComponent::DataType2Text(dataType).data();
    }
    //filename+=Form("", specification);
  }
  return iResult;
}

int CheckOutputFiles(AliHLTComponentBlockDataList blocks, vector<TestElement_t>& testData, 
		     bool bConcatenateEvents, bool bConcatenateBlocks,
		     unsigned int nofEvents=1)
{
  int iResult=0;
  AliHLTUInt32_t fileOffset=0;
  TString lastFileName;
  TFile* pFile=NULL;
  TArrayC buffer;
  for (unsigned eventNo=0; eventNo<nofEvents && iResult>=0; eventNo++) {
    SetEvent(eventNo, testData);
    for (unsigned blockNo=0; blockNo<blocks.size(); blockNo++) {
      TString fileName;
      BuildFileName(eventNo, blockNo, blocks[blockNo].fDataType, blocks[blockNo].fSpecification, bConcatenateEvents, bConcatenateBlocks, fileName);
      if (pFile && !fileName.IsNull() && fileName!=lastFileName) {
	pFile->Close();
	delete pFile;
	pFile=NULL;
	fileOffset=0;
      }
      lastFileName=fileName;
      fileName+="?filetype=raw";
      if (!pFile) pFile=new TFile(fileName);
      if (!pFile || pFile->IsZombie()) {
	cerr << "failed: can not open file " << fileName.ReplaceAll("?filetype=raw", "") << endl;
	iResult=-ENOENT;
	break;
      }
      if (buffer.GetSize()<(int)blocks[blockNo].fSize) {
	buffer.Set(blocks[blockNo].fSize);
      }
      if ((pFile->ReadBuffer((char*)buffer.GetArray(), blocks[blockNo].fSize))!=0) {
	cerr << "failed: reading " << blocks[blockNo].fSize << " at offset " << fileOffset << " from file " << fileName.ReplaceAll("?filetype=raw", "") << endl;
	iResult=-EIO;
	break;
      }
      AliHLTUInt8_t* origin=(AliHLTUInt8_t*)blocks[blockNo].fPtr;
      AliHLTUInt8_t* processed=(AliHLTUInt8_t*)buffer.GetArray();
      if (memcmp(origin, processed, blocks[blockNo].fSize)!=0) {
	cerr << "failed: data missmatch in event " << eventNo << " block " << blockNo << endl;
	iResult=-EFAULT;
	break;
      }
      cout << "checked: file " << fileName.ReplaceAll("?filetype=raw", "") << " offset " << fileOffset << endl;
      fileOffset+=blocks[blockNo].fSize;
    }
  }
  if (pFile) {
    pFile->Close();
    delete pFile;
    pFile=NULL;
  }
  return iResult;
}

int testAliHLTFileWriter()
{
  int iResult=0;
  AliHLTComponentHandler chandler;
  if ((iResult=chandler.LoadLibrary("../.libs/libAliHLTUtil.so"))<0) {
    cerr << "failed: loading libAliHLTUtil" << endl;
    return iResult;
  }

  // create FileWriter component
  AliHLTComponent* pComponent=NULL;
  if ((iResult=chandler.CreateComponent("FileWriter", pComponent))<0) {
    cerr << "failed: can not create component \"FileWriter\"" << endl;
    return iResult;
  }

  AliHLTComponentBlockDataList blocks;
  const unsigned testDataSize=500;
  vector<TestElement_t> testData;
  testData.assign(testDataSize, 0);
  unsigned int i=0;
  for (; i<testDataSize; i++) testData[i]=i;

  const unsigned elementsPerBlock=100;
  AliHLTComponentBlockData bd;
  AliHLTComponent::FillBlockData(bd);
  for (i=0; i*elementsPerBlock<testData.size(); i++) {
    bd.fPtr=&(testData[i*elementsPerBlock]);
    bd.fSize=elementsPerBlock*sizeof(TestElement_t);
    bd.fSpecification=i;
    AliHLTComponent::SetDataType(bd.fDataType, "DUMMYDAT", "TEST");
    blocks.push_back(bd);
  }
  
  const char* testArguments[]={
    "",
    "-concatenate-blocks ",
    "-concatenate-blocks -concatenate-events ",
    // not currently working   "-burst-buffer 10000 ",
    // not currently working   "-burst-buffer 5000 ",
    // not currently working   "-burst-buffer 10000 -concatenate-blocks ",
    "-burst-buffer 10000 -concatenate-blocks -concatenate-events ",
    "-burst-buffer 5000 -concatenate-blocks -concatenate-events ",
    NULL
  };

  for (const char** currentTest=testArguments; *currentTest!=NULL; currentTest++) {
    TString arguments=*currentTest;
      arguments+=" -directory "; arguments+=testDirectory;
    if (!gSystem->AccessPathName(testDirectory)) {
      TString command="rm -r "; command+=testDirectory;
      gSystem->Exec(command);
    }
    cout << arguments.Data() << endl;
    if ((iResult=TestCycle(pComponent, arguments.Data(), blocks, testData, 4))<0) {
      break;
    }
    if ((iResult=CheckOutputFiles(blocks, testData, arguments.Contains("-concatenate-events "), arguments.Contains("-concatenate-blocks"), 4))<0) {
      break;
    }
  };

  delete pComponent;

  return iResult;
}

int main(int /*argc*/, const char** /*argv*/)
{
  int iResult=0;
  if ((iResult=testAliHLTFileWriter())<0) {
  }
  return iResult;
}
