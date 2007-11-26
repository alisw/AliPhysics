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

/** @file   AliHLTTPCDigitDumpComponent.cxx
    @author Matthias Richter
    @date   
    @brief  Special file writer converting TPC digit input to ASCII. */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <cassert>
#include "AliHLTTPCDigitDumpComponent.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCDigitReaderRaw.h"
#include "AliHLTTPCDefinitions.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCDigitDumpComponent)

AliHLTTPCDigitDumpComponent::AliHLTTPCDigitDumpComponent()
  :
  AliHLTFileWriter(),
  fRawreaderMode(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCDigitDumpComponent::~AliHLTTPCDigitDumpComponent()
{
  // see header file for class documentation
}

const char* AliHLTTPCDigitDumpComponent::GetComponentID()
{
  // see header file for class documentation
  return "TPCDigitDump";
}

void AliHLTTPCDigitDumpComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back(kAliHLTAnyDataType);
}

AliHLTComponent* AliHLTTPCDigitDumpComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCDigitDumpComponent;
}

int AliHLTTPCDigitDumpComponent::InitWriter()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTPCDigitDumpComponent::ScanArgument(int argc, const char** argv)
{
  // see header file for class documentation
  int iResult=0;
  TString argument="";
  bool bMissingParam=0;
  int i=0;
  for (; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

    // -rawreadermode
    if (argument.CompareTo("-rawreadermode")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      int mode=AliHLTTPCDigitReaderRaw::DecodeMode(argv[i]);
      if (mode<0) {
	HLTError("invalid rawreadermode specifier '%s'", argv[i]);
	iResult=-EINVAL;
      } else {
	fRawreaderMode=static_cast<unsigned>(mode);
      }
      break;
    }
  }

  if (bMissingParam) {
    iResult=-EPROTO;
  }
  if (iResult>=0) iResult=i+1;

  return iResult;
}

int AliHLTTPCDigitDumpComponent::CloseWriter()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTPCDigitDumpComponent::DumpEvent( const AliHLTComponentEventData& evtData,
					    const AliHLTComponentBlockData* blocks, 
					    AliHLTComponentTriggerData& /*trigData*/ )
{
  // see header file for class documentation
  int iResult=0;
  int iPrintedSlice=-1;
  int iPrintedPart=-1;
  int blockno=0;
  HLTDebug("%d blocks", evtData.fBlockCnt);
  const AliHLTComponentBlockData* pDesc=NULL;

  for (pDesc=GetFirstInputBlock(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTPC); pDesc!=NULL; pDesc=GetNextInputBlock(), blockno++) {
    // TODO: Matthias 26.11.2007
    // I have a very strange behavior with the processing of that log message leading to a crash in vsnprintf
    // Maybe it's just trivial, but then it should not crash. Needs more investigation
    //HLTDebug("event %d block %d: %s 0x%08x size %d", evtData.fEventID, blockno, DataType2Text(pDesc->fDataType).c_str(), pDesc->fSpecification, pDesc->fSize);
    HLTDebug("dump digits for block %d specification 0x%08x", blockno, pDesc->fSpecification);
    TString filename;
    iResult=BuildFileName(evtData.fEventID, blockno, pDesc->fDataType, pDesc->fSpecification, filename);
    ios::openmode filemode=(ios::openmode)0;
    if (fCurrentFileName.CompareTo(filename)==0) {
      // append to the file
      filemode=ios::app;
    } else {
      // store the file for the next block
      fCurrentFileName=filename;
    }
    if (iResult>=0) {
      ofstream dump(filename.Data(), filemode);
      if (dump.good()) {
	int part=AliHLTTPCDefinitions::GetMinPatchNr(*pDesc);
	assert(part==AliHLTTPCDefinitions::GetMaxPatchNr(*pDesc));
	int slice=AliHLTTPCDefinitions::GetMinSliceNr(*pDesc);
	assert(slice==AliHLTTPCDefinitions::GetMaxSliceNr(*pDesc));
	int firstRow=AliHLTTPCTransform::GetFirstRow(part);
	int lastRow=AliHLTTPCTransform::GetLastRow(part);
	AliHLTTPCDigitReaderRaw reader(fRawreaderMode);
	reader.InitBlock(pDesc->fPtr,pDesc->fSize,firstRow,lastRow,part,slice);

	int iPrintedRow=-1;
	int iPrintedPad=-1;
	int iLastTime=-1;
	while (reader.Next()) {
	  if (iPrintedSlice!=slice || iPrintedPart!=part) {
	    iPrintedSlice=slice;
	    iPrintedPart=part;
	    dump << endl;
	    dump << "====================================================================" << endl;
	    dump << "    Slice: " << iPrintedSlice << "   Pad: " << iPrintedPad;
	  }
	  if (iPrintedRow!=reader.GetRow()) {
	    iPrintedRow=reader.GetRow();
	    dump << endl;
	    dump << "--------------------------------------------------------------------" << endl;
	    dump << "Row: " << iPrintedRow << endl;
	  }
	  if (iPrintedPad!=reader.GetPad()) {
	    dump << endl;
	    iPrintedPad=reader.GetPad();
	    dump << "    Pad: " << iPrintedPad;
	  }
	  if (iLastTime!=reader.GetTime()+1 && iLastTime!=reader.GetTime()-1 ) {
	    dump << endl;
	    dump << "        Time: " << reader.GetTime();
	  }
	  iLastTime=reader.GetTime();
	  dump << "  " << reader.GetSignal();
	}
	dump << endl << endl;
      } else {
	HLTError("can not open file %s for writing", filename.Data());
	iResult=-EBADF;
      }
      dump.close();
    }
  }
  return iResult;
}
