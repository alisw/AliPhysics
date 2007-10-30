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

/** @file   AliHLTOUTComponent.cxx
    @author Matthias Richter
    @date   
    @brief  The HLTOUT data sink component similar to HLTOUT nodes */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
using namespace std;
#endif

#include <cassert>
#include <iostream>
#include "AliHLTOUTComponent.h"
#include "AliHLTOUT.h"
#include "AliHLTHOMERWriter.h"
#include "AliDAQ.h" // equipment Ids
#include "AliRawDataHeader.h" // Common Data Header 
#include <TDatime.h> // seed for TRandom
#include <TRandom.h> // random int generation for DDL no

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOUTComponent)

AliHLTOUTComponent::AliHLTOUTComponent()
  :
  AliHLTOfflineDataSink(),
  fWriters(),
  fNofDDLs(10),
  fIdFirstDDL(4864), // 0x13<<8
  fWriteDigits(kTRUE),
  fWriteRaw(kTRUE),
  fBuffer(),
  fpLibManager(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  // I guess DDL definitions should never change any more
  assert(fNofDDLs==AliDAQ::NumberOfDdls("HLT"));
  fNofDDLs=AliDAQ::NumberOfDdls("HLT");
  assert(fIdFirstDDL==AliDAQ::DdlIDOffset("HLT"));
  fIdFirstDDL=AliDAQ::DdlIDOffset("HLT");
}

AliHLTOUTComponent::~AliHLTOUTComponent()
{
  // see header file for class documentation
  if (fpLibManager) delete fpLibManager;
  fpLibManager=NULL;
}

const char* AliHLTOUTComponent::GetComponentID()
{
  // see header file for class documentation
  return "HLTOUT";
}

void AliHLTOUTComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back(kAliHLTAnyDataType);
}

AliHLTComponent* AliHLTOUTComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTOUTComponent;
}

int AliHLTOUTComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;
  TString argument="";
  int bMissingParam=0;
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

    {
      HLTError("unknown argument %s", argument.Data());
      break;
    }
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  if (iResult>=0) {
  }

  fpLibManager=new AliHLTHOMERLibManager;
  if (fpLibManager) {
    int writerNo=0;
    for (writerNo=0; writerNo<fNofDDLs; writerNo++) {
      AliHLTMonitoringWriter* pWriter=fpLibManager->OpenWriter();
      if (pWriter) {
	fWriters.push_back(pWriter);
      } else {
	iResult=-ENOMEM;
	break;
      }
    }
  } else {
    iResult=-ENOMEM;
  }

  return iResult;
}

int AliHLTOUTComponent::DoDeinit()
{
  // see header file for class documentation
  int iResult=0;

  if (fpLibManager) {
    AliHLTMonitoringWriterPVector::iterator element=fWriters.begin();
    while (element!= fWriters.end()) {
      assert(*element);
      if (*element!=NULL) fpLibManager->DeleteWriter(dynamic_cast<AliHLTHOMERWriter*>(*element));
      element=fWriters.erase(element);
    }
  }
  
  return iResult;
}

int AliHLTOUTComponent::DumpEvent( const AliHLTComponentEventData& evtData,
			 const AliHLTComponentBlockData* blocks, 
			 AliHLTComponentTriggerData& /*trigData*/ )
{
  // see header file for class documentation
  int iResult=0;
  HLTInfo("write %d output blocks", evtData.fBlockCnt);
  fWriters.clear();
  if (iResult>=0) {
    homer_uint64 homerHeader[kCount_64b_Words];
    HOMERBlockDescriptor homerDescriptor(homerHeader);
    for (int n=0; n<(int)evtData.fBlockCnt; n++ ) {
      memset( homerHeader, 0, sizeof(homer_uint64)*kCount_64b_Words );
      homerDescriptor.Initialize();
      homerDescriptor.SetType(reinterpret_cast<homer_uint64>(blocks[n].fDataType.fID));
      homerDescriptor.SetSubType1(reinterpret_cast<homer_uint64>(blocks[n].fDataType.fOrigin));
      homerDescriptor.SetSubType2(static_cast<homer_uint64>(blocks[n].fSpecification));
      int writerNo=ShuffleWriters(fWriters, blocks[n].fSize);
      assert(writerNo>=0 && writerNo<fWriters.size());
      fWriters[writerNo]->AddBlock(&homerDescriptor, blocks[n].fPtr);
    }
  }

  return iResult;
}

int AliHLTOUTComponent::FillESD(int eventNo, AliRunLoader* runLoader, AliESDEvent* /*esd*/)
{
  // see header file for class documentation
  int iResult=0;
  if (fWriters.size()==0) return 0;
  
  // search for the writer with the biggest data volume in order to allocate the
  // output buffer of sufficient size
  AliHLTMonitoringWriterPVector::iterator writer=fWriters.begin();
  vector<int> sorted;
  for (int i=0; i<fWriters.size(); i++) {
    assert(fWriters[i]);
    if (fWriters[i]) {
      if (sorted.size()>=0 && fWriters[i]->GetTotalMemorySize()>fWriters[sorted[0]]->GetTotalMemorySize()) {
	sorted.insert(sorted.begin(), i);
      } else {
	sorted.push_back(i);
      }
    }
    writer++;
  }

  vector<int>::iterator ddlno=sorted.begin();
  while (ddlno!=sorted.end()) {
    const AliHLTUInt8_t* pBuffer=NULL;
    int bufferSize=0;
    
    if ((bufferSize=FillOutputBuffer(eventNo, *writer, pBuffer))>0) {
      if (fWriteDigits) WriteDigits(eventNo, runLoader, *ddlno, pBuffer, bufferSize);
      if (fWriteRaw) WriteRawFile(eventNo, runLoader, *ddlno, pBuffer, bufferSize);
    }
    ddlno++;
  }
  return iResult;
}

int AliHLTOUTComponent::ShuffleWriters(AliHLTMonitoringWriterPVector &list, AliHLTUInt32_t size)
{
  // see header file for class documentation
  int iResult=-ENOENT;
  assert(list.size()>0);
  if (list.size()==0) return iResult;
  vector<int> writers;
  int i=0;
  for (i=0; i<list.size(); i++) {
    if (list[i]->GetTotalMemorySize()==0)
      writers.push_back(i);
    else if (iResult<0 ||
	     list[i]->GetTotalMemorySize()<list[iResult]->GetTotalMemorySize())
      iResult=i;
      
  }
  if (writers.size()>0) {
    iResult=writers[0];
    if (writers.size()>0) {
      // shuffle among the empty writers
      TDatime dt;
      TRandom rand;
      rand.SetSeed(dt.Get()*(iResult+1));
      i=rand.Integer(writers.size()-1);
      assert(i>0 && i<writers.size()-1);
      iResult=writers[i];
    }
  } else {
    // take the writer with the least data volume
    assert(iResult>=0);
  }
  return iResult;
}

int AliHLTOUTComponent::FillOutputBuffer(int eventNo, AliHLTMonitoringWriter* pWriter, const AliHLTUInt8_t* &pBuffer)
{
  // see header file for class documentation
  int iResult=0;
  int bufferSize=0;

  // space for common data header
  bufferSize+=sizeof(AliRawDataHeader);
  assert(sizeof(AliRawDataHeader)==24);

  // space for HLT event header
  bufferSize+=sizeof(AliHLTOUT::AliHLTOUTEventHeader);

  // space for payload from the writer
  if (pWriter) bufferSize+=pWriter->GetTotalMemorySize();

  if (bufferSize>fBuffer.size())
    fBuffer.resize(bufferSize);

  if (bufferSize<=fBuffer.size()) {
    AliRawDataHeader* pCDH=reinterpret_cast<AliRawDataHeader*>(&fBuffer[0]);
    AliHLTOUT::AliHLTOUTEventHeader* pHLTH=reinterpret_cast<AliHLTOUT::AliHLTOUTEventHeader*>(&fBuffer[sizeof(AliRawDataHeader)]);
    memset(pCDH, 0, sizeof(AliRawDataHeader));
    memset(pHLTH, 0, sizeof(AliHLTOUT::AliHLTOUTEventHeader));
    pHLTH->fVersion=1;
    if (pWriter) {
      // copy payload
      pWriter->Copy(&fBuffer[sizeof(AliRawDataHeader)+sizeof(AliHLTOUT::AliHLTOUTEventHeader)], 0, 0, 0, 0);
      pHLTH->fLength=pWriter->GetTotalMemorySize();
      // set status bit to indicate HLT payload
      pCDH->fStatusMiniEventID|=0x1<<(AliHLTOUT::fgkCDHStatusFlagsOffset+AliHLTOUT::fgkCDHFlagsHLTPayload);
    }
    pHLTH->fLength+=sizeof(AliHLTOUT::AliHLTOUTEventHeader);
    pHLTH->fEventID=eventNo;

    pCDH->fSize=sizeof(AliRawDataHeader)+pHLTH->fLength;
    
    pBuffer=&fBuffer[0];
  } else {
    pBuffer=NULL;
    iResult=-ENOMEM;
  }

  return iResult;
}

int AliHLTOUTComponent::WriteDigits(int eventNo, AliRunLoader* runLoader, int hltddl, const AliHLTUInt8_t* pBuffer, int bufferSize)
{
  // see header file for class documentation
  int iResult=0;
  return iResult;
}

int AliHLTOUTComponent::WriteRawFile(int eventNo, AliRunLoader* runLoader, int hltddl, const AliHLTUInt8_t* pBuffer, int bufferSize)
{
  // see header file for class documentation
  int iResult=0;
  const char* fileName=AliDAQ::DdlFileName("HLT", hltddl);
  assert(fileName!=NULL);
  TString filePath(fileName);
  if (fileName) {
    ios::openmode filemode=(ios::openmode)0;
    ofstream rawfile(filePath.Data(), filemode);
    if (rawfile.good()) {
      if (pBuffer && bufferSize>0) {
	rawfile.write(reinterpret_cast<const char*>(pBuffer), bufferSize);
      } else {
	HLTWarning("writing zero length raw data file %s");
      }
      HLTDebug("wrote %d byte(s) to file %s", bufferSize, filePath.Data());
    } else {
      HLTError("can not open file %s for writing", filePath.Data());
      iResult=-EBADF;
    }
    rawfile.close();
  }
  return iResult;
}
