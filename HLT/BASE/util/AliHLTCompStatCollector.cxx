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

/** @file   AliHLTCompStatCollector.cxx
    @author Matthias Richter
    @date   
    @brief  Collector component for the component statistics information.
*/

#include "AliHLTCompStatCollector.h"
#include "TStopwatch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2C.h"
#include "TTree.h"
#include "TString.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTCompStatCollector)

AliHLTCompStatCollector::AliHLTCompStatCollector()
  :
  AliHLTProcessor(),
  fpTimer(NULL),
  fpStatTree(NULL),
  fCycleTime(0),
  fNofSets(0),
  fArraySize(1000),
  fPosition(0),
  fpLevelArray(NULL),
  fpSpecArray(NULL),
  fpBlockNoArray(NULL),
  fpIdArray(NULL),
  fpTimeArray(NULL),
  fpCTimeArray(NULL),
  fpInputBlockCountArray(NULL),
  fpTotalInputSizeArray(NULL),
  fpOutputBlockCountArray(NULL),
  fpTotalOutputSizeArray(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTCompStatCollector::~AliHLTCompStatCollector()
{
  // see header file for class documentation
  ClearAll();
}

void AliHLTCompStatCollector::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.push_back(kAliHLTDataTypeComponentStatistics);
}

AliHLTComponentDataType AliHLTCompStatCollector::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTMultipleDataType;
}

int AliHLTCompStatCollector::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  // see header file for class documentation
  tgtList.clear();
  tgtList.push_back(kAliHLTDataTypeHistogram);
  tgtList.push_back(kAliHLTDataTypeTTree);
  return tgtList.size();
}

void AliHLTCompStatCollector::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  constBase=1000;
  inputMultiplier=100.0;
}

int AliHLTCompStatCollector::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;
  TString argument="";
  int bMissingParam=0;
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

    // -array-size
    if (argument.CompareTo("-array-size")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      TString param=argv[argc];
      if (param.IsDigit()) {
	fArraySize=param.Atoi();
      } else {
	HLTError("expecting number as parameter for option '-array-size'");
	iResult=-EINVAL;
      }

    } else {
      iResult=-EINVAL;
    }
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }

  if (iResult>=0) {
    fpLevelArray=new UInt_t[fArraySize];
    fpSpecArray=new UInt_t[fArraySize];
    fpBlockNoArray=new UInt_t[fArraySize];
    fpIdArray=new UInt_t[fArraySize];
    fpTimeArray=new UInt_t[fArraySize];
    fpCTimeArray=new UInt_t[fArraySize];
    fpInputBlockCountArray=new UInt_t[fArraySize];
    fpTotalInputSizeArray=new UInt_t[fArraySize];
    fpOutputBlockCountArray=new UInt_t[fArraySize];
    fpTotalOutputSizeArray=new UInt_t[fArraySize];

    fpStatTree=new TTree("CompStat", "HLT component statistics");
    if (fpStatTree) {
      fpStatTree->SetDirectory(0);
      fpStatTree->Branch("cycleTime",        &fCycleTime, "cycleTime/F");
      fpStatTree->Branch("nofSets",          &fNofSets, "nofSets/I");
      fpStatTree->Branch("Level",            fpLevelArray, "Level[nofSets]/i");
      fpStatTree->Branch("Specification",    fpSpecArray, "Specification[nofSets]/i");
      fpStatTree->Branch("BlockNo",          fpBlockNoArray, "BlockNo[nofSets]/i");
      fpStatTree->Branch("Id",               fpIdArray, "Id[nofSets]/i");
      fpStatTree->Branch("Time",             fpTimeArray, "Time[nofSets]/i");
      fpStatTree->Branch("CTime",            fpCTimeArray, "CTime[nofSets]/i");
      fpStatTree->Branch("InputBlockCount",  fpInputBlockCountArray, "InputBlockCount[nofSets]/i");
      fpStatTree->Branch("TotalInputSize",   fpTotalInputSizeArray, "TotalInputSize[nofSets]/i");
      fpStatTree->Branch("OutputBlockCount", fpOutputBlockCountArray, "OutputBlockCount[nofSets]/i");
      fpStatTree->Branch("TotalOutputSize",  fpTotalOutputSizeArray, "TotalOutputSize[nofSets]/i");
    }
  }
  return iResult;
}

int AliHLTCompStatCollector::DoDeinit( )
{
  // see header file for class documentation
  ClearAll();

  return 0;
}

int AliHLTCompStatCollector::DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  // see header file for class documentation
  int iResult=0;

  ResetFillingVariables();
  if (fpTimer) {
    fCycleTime=fpTimer->RealTime()*1000000;
  }

  int blockNo=0;
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeComponentStatistics);
       pBlock && iResult>=0;
       pBlock=GetNextInputBlock(), blockNo++) {
    unsigned int current=fPosition;
    iResult=FillVariablesSorted(pBlock->fPtr, pBlock->fSize);
    for (; current<fPosition; current++) {
      fpSpecArray[current]=pBlock->fSpecification;
      fpBlockNoArray[current]=blockNo;
    }
    // indicate availability of component statistic block
    iResult=1;
  }

  if (iResult>0) {
    fNofSets=fPosition;
    fpStatTree->Fill();
    iResult=PushBack(fpStatTree, kAliHLTDataTypeTTree);

    // init the timer for the next cycle
    if (!fpTimer)  fpTimer=new TStopwatch;
    if (fpTimer) {
      fpTimer->Reset();
      fpTimer->Start();
    }
  }

  if (iResult>0) iResult=-1;
  return iResult;
}

void AliHLTCompStatCollector::ResetFillingVariables()
{
  // see header file for class documentation
  fCycleTime=0;
  fNofSets=0;
  fPosition=0;
  memset(fpLevelArray, 0, sizeof(UInt_t)*fArraySize);
  memset(fpSpecArray, 0, sizeof(UInt_t)*fArraySize);
  memset(fpBlockNoArray, 0, sizeof(UInt_t)*fArraySize);
  memset(fpIdArray, 0, sizeof(UInt_t)*fArraySize);
  memset(fpTimeArray, 0, sizeof(UInt_t)*fArraySize);
  memset(fpCTimeArray, 0, sizeof(UInt_t)*fArraySize);
  memset(fpInputBlockCountArray, 0, sizeof(UInt_t)*fArraySize);
  memset(fpTotalInputSizeArray, 0, sizeof(UInt_t)*fArraySize);
  memset(fpOutputBlockCountArray, 0, sizeof(UInt_t)*fArraySize);
  memset(fpTotalOutputSizeArray, 0, sizeof(UInt_t)*fArraySize);
}

int AliHLTCompStatCollector::FillVariablesSorted(void* ptr, int size)
{
  // see header file for class documentation
  int iResult=0;
  if (size%sizeof(AliHLTComponentStatistics)) {
    HLTError("data block is not aligned to the size of the AliHLTComponentStatistics struct");
    return -EINVAL;
  }
  AliHLTComponentStatistics* pStat=reinterpret_cast<AliHLTComponentStatistics*>(ptr);
  UInt_t nofStats=size/sizeof(AliHLTComponentStatistics);
  vector<int> indexList;
  UInt_t i=0;
  for (i=0; i<nofStats; i++) {
    vector<int>::iterator element=indexList.begin();
    for (; element!=indexList.end(); element++) {
      if (pStat[i].fLevel>pStat[*element].fLevel) {
	break;
      }
    }
    indexList.insert(element, i);
  }

  i=fPosition;
  for (vector<int>::iterator element=indexList.begin();
       element!=indexList.end();
       element++, i++) {
    if (i<fArraySize) {
      fpLevelArray[i]=pStat[*element].fLevel;
      fpIdArray[i]=pStat[*element].fId;
      fpTimeArray[i]=pStat[*element].fTime;
      fpCTimeArray[i]=pStat[*element].fCTime;
      fpInputBlockCountArray[i]=pStat[*element].fInputBlockCount;
      fpTotalInputSizeArray[i]=pStat[*element].fTotalInputSize;
      fpOutputBlockCountArray[i]=pStat[*element].fOutputBlockCount;
      fpTotalOutputSizeArray[i]=pStat[*element].fTotalOutputSize;
    }
  }

  if (i>=fArraySize) {
    HLTWarning("too little space in branch variables to fill %d statistics blocks, available %d at position %d", i, fArraySize, fPosition);
    fPosition=fArraySize;
  } else {
    fPosition=i;
  }
  
  return iResult;
}

void AliHLTCompStatCollector::ClearAll()
{
  if (fpTimer) delete fpTimer; fpTimer=NULL;
  if (fpStatTree) delete fpStatTree; fpStatTree=NULL;
  if (fpLevelArray) delete fpLevelArray; fpLevelArray=NULL;
  if (fpSpecArray) delete fpSpecArray; fpSpecArray=NULL;
  if (fpBlockNoArray) delete fpBlockNoArray; fpBlockNoArray=NULL;
  if (fpIdArray) delete fpIdArray; fpIdArray=NULL;
  if (fpTimeArray) delete fpTimeArray; fpTimeArray=NULL;
  if (fpCTimeArray) delete fpCTimeArray; fpCTimeArray=NULL;
  if (fpInputBlockCountArray) delete fpInputBlockCountArray; fpInputBlockCountArray=NULL;
  if (fpTotalInputSizeArray) delete fpTotalInputSizeArray; fpTotalInputSizeArray=NULL;
  if (fpOutputBlockCountArray) delete fpOutputBlockCountArray; fpOutputBlockCountArray=NULL;
  if (fpTotalOutputSizeArray) delete fpTotalOutputSizeArray; fpTotalOutputSizeArray=NULL;
}
