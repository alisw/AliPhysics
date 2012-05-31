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
#include "AliHLTTPCDigitReader.h"
#include "AliHLTTPCDigitReaderUnpacked.h"
#include "AliHLTTPCDigitReader32Bit.h"
#include "AliHLTTPCDefinitions.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCDigitDumpComponent)

AliHLTTPCDigitDumpComponent::AliHLTTPCDigitDumpComponent()
  :
  AliHLTFileWriter(),
  fDigitReaderType(kDigitReader32Bit),
  fRcuTrailerSize(2),
  fUnsorted(true),
  fbBulkMode(true),
  fpReader(NULL),
  f32BitFormat(kFALSE)
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
  int iResult=0;
  switch (fDigitReaderType) {
  case kDigitReaderUnpacked:
    HLTInfo("create DigitReaderUnpacked");
    fpReader=new AliHLTTPCDigitReaderUnpacked; 
    break;
  case kDigitReader32Bit:
    HLTInfo("create DigitReader32Bit");
    fpReader=new AliHLTTPCDigitReader32Bit();
    f32BitFormat = kTRUE;
    break;
  }
  if (!fpReader) {
    HLTError("can not create digit reader of type %d", fDigitReaderType);
    iResult=-EFAULT;
  } else {
    fpReader->SetUnsorted(fUnsorted);
  }
  return iResult;
}

int AliHLTTPCDigitDumpComponent::ScanArgument(int argc, const char** argv)
{
  // see header file for class documentation
  int iResult=0;
  TString argument="";
  bool bMissingParam=0;
  int i=0;
  do {
    if (i>=argc || (argument=argv[i]).IsNull()) continue;

    // -rawreadermode
    if (argument.CompareTo("-rawreadermode")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      HLTWarning("argument '-rawreadermode' deprecated");
      break;
    }

    // -digitreader
    if (argument.CompareTo("-digitreader")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      TString param=argv[i];
      if (param.CompareTo("unpacked", TString::kIgnoreCase)==0) {
	fDigitReaderType=kDigitReaderUnpacked;
      } else if (param.CompareTo("packed", TString::kIgnoreCase)==0) {
	HLTWarning("argument 'packed' is deprecated, falling back to DigitReader32Bit");
	fDigitReaderType=kDigitReader32Bit;
      } else if (param.CompareTo("raw", TString::kIgnoreCase)==0) {
	HLTWarning("argument 'raw' is deprecated, falling back to DigitReader32Bit");
	fDigitReaderType=kDigitReader32Bit;
      } else if (param.CompareTo("decoder", TString::kIgnoreCase)==0) {
	HLTWarning("argument 'decoder' is deprecated, falling back to DigitReader32Bit");
	fDigitReaderType=kDigitReader32Bit;
      } else if (param.CompareTo("32bit", TString::kIgnoreCase)==0) {
	fDigitReaderType=kDigitReader32Bit;
      } else {
	HLTError("unknown digit reader type %s", param.Data());
	iResult=-EINVAL;
      }

      break;
    }

    // -rcutrailersize
    if (argument.CompareTo("-rcutrailersize")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      char *endptr=NULL;
      fRcuTrailerSize=strtoul(argv[i], &endptr, 0);
      if (/*endptr ||*/ fRcuTrailerSize<1) {
	HLTError("invalid parameter '%s', %s", argv[i], endptr==NULL?"number >= 1 expected":"can not convert string to number");
	iResult=-EINVAL;
      }
      break;
    }

    // -unsorted
    if (argument.CompareTo("-unsorted")==0) {
      fUnsorted=true;
      break;
    }

    // -sorted
    if (argument.CompareTo("-sorted")==0) {
      fUnsorted=false;
      break;
    }

    // -bulk
    if (argument.CompareTo("-bulk")==0) {
      fbBulkMode=true;
      break;
    }

    // -stream
    if (argument.CompareTo("-stream")==0) {
      fbBulkMode=false;
      break;
    }
  } while (0); // just use the do/while here to have the option of breaking

  if (bMissingParam) iResult=-EPROTO;
  else if (iResult>=0) iResult=i;

  return iResult;
}

int AliHLTTPCDigitDumpComponent::CloseWriter()
{
  // see header file for class documentation
  if (fpReader) delete fpReader;
  fpReader=NULL;
  return 0;
}

int AliHLTTPCDigitDumpComponent::DumpEvent( const AliHLTComponentEventData& evtData,
					    const AliHLTComponentBlockData* /*blocks*/, 
					    AliHLTComponentTriggerData& /*trigData*/ )
{
  // see header file for class documentation
  int iResult=0;
  int iPrintedSlice=-1;
  int iPrintedPart=-1;
  int blockno=0;
  const AliHLTComponentBlockData* pDesc=NULL;

  AliHLTTPCDigitReader* pReader=fpReader;
  if (!pReader) return -ENODEV;

  for (pDesc=GetFirstInputBlock(kAliHLTAnyDataType); pDesc!=NULL; pDesc=GetNextInputBlock(), blockno++) {
    HLTDebug("event %Lu block %d: %s 0x%08x size %d", evtData.fEventID, blockno, DataType2Text(pDesc->fDataType).c_str(), pDesc->fSpecification, pDesc->fSize);

    if (fDigitReaderType==kDigitReaderUnpacked && pDesc->fDataType!=AliHLTTPCDefinitions::fgkUnpackedRawDataType) continue;
    else if (fDigitReaderType!=kDigitReaderUnpacked && pDesc->fDataType!=(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTPC)) continue;

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

	iResult=pReader->InitBlock(pDesc->fPtr,pDesc->fSize,firstRow,lastRow,part,slice);

	int iPrintedRow=-1;
	int iPrintedPad=-1;
	int iLastTime=-1;
	if (fbBulkMode) {
	  while (pReader->NextChannel()) {
	    if (PrintHeaders(slice, iPrintedSlice, part, iPrintedPart, pReader, iPrintedRow, iPrintedPad, dump)) {
	      iLastTime=-1;
	    }
	    while (pReader->NextBunch()) {
	      int bunchLength=pReader->GetBunchSize();
	      
	      // Kenneth: 20-04-09. The following if have been added because of inconsistency in the 40 bit decoder and the 32 bit decoder.
	      // GetSignals() in the 40 bit decoder returns an array of UInt_t while the 32 bit one returns UShort_t
	      if(f32BitFormat == kTRUE){
		const  UShort_t* bunchData=pReader->GetSignalsShort();
		
		// bunch data is printed in 'reverse' order in order to produce
		// the same output as in stream reading mode
		dump << "                     Time " << pReader->GetTime()+bunchLength-1 << ":  ";
		for (int bin=bunchLength-1; bin>=0; bin--) {
		  dump << "  " << bunchData[bin];
		}
		dump << "    -> Time: " << pReader->GetTime() << endl;
	      }
	      else{
		const  UInt_t* bunchData=pReader->GetSignals();
		dump << "                     Time " << pReader->GetTime()+bunchLength-1 << ":  ";
		for (int bin=0; bin<bunchLength; bin++) {
		  dump << "  " << bunchData[bin];
		}
		dump << "    -> Time: " << pReader->GetTime() << endl;
	      }
	    }
	  }
	  dump << endl;
	} else {
	while (pReader->Next()) {
	  if ((iPrintedSlice!=-1 && iLastTime!=-1 && iLastTime!=pReader->GetTime()+1 && iLastTime!=pReader->GetTime()-1)) {
	    dump << "    -> Time: " << iLastTime << endl;
	  } else if ((iPrintedPad!=-1 && iPrintedPad!=pReader->GetPad()) ||
		     (iPrintedRow!=-1 && iPrintedRow!=pReader->GetRow())) {
	    dump << "    -> Time: " << iLastTime << endl;
	    //dump << endl;
	  }

	  if (PrintHeaders(slice, iPrintedSlice, part, iPrintedPart, pReader, iPrintedRow, iPrintedPad, dump)) {
	    iLastTime=-1;
	  }
	  if (iLastTime==-1 || (iLastTime!=pReader->GetTime()+1 && iLastTime!=pReader->GetTime()-1)) {
	    dump << "                     Time " << pReader->GetTime() << ":  ";
	  }
	  iLastTime=pReader->GetTime();
	  dump << "  " << pReader->GetSignal();
	}
	if (iLastTime>=0) dump << "    -> Time: " << iLastTime << endl << endl;
	}
      } else {
	HLTError("can not open file %s for writing", filename.Data());
	iResult=-EBADF;
      }
      dump.close();
    }
    pReader->Reset();
  }
  return iResult;
}

int AliHLTTPCDigitDumpComponent::PrintHeaders(int slice, int &iPrintedSlice,
					      int part, int &iPrintedPart,
					      AliHLTTPCDigitReader* pReader,
					      int &iPrintedRow, int &iPrintedPad,
					      ofstream &dump) const
{
  // see header file for class documentation
  int iResult=0;
  assert(pReader);
  if (iPrintedSlice!=slice || iPrintedPart!=part) {
    iPrintedSlice=slice;
    iPrintedPart=part;
    dump << "====================================================================" << endl;
    dump << "    Slice: " << iPrintedSlice << "   Partition: " << iPrintedPart << endl;
    iPrintedRow=-1;
  }
  if (iPrintedRow!=pReader->GetRow()) {
    iPrintedRow=pReader->GetRow();
    dump << "--------------------------------------------------------------------" << endl;
    dump << "Row: " << iPrintedRow << endl;
    iPrintedPad=-1;
  }
  if (iPrintedPad!=pReader->GetPad()) {
    iPrintedPad=pReader->GetPad();
    dump << "Row: " << iPrintedRow << "  Pad: " << iPrintedPad << "  HW address: " << pReader->GetAltroBlockHWaddr() << endl;
    iResult=1;
  }

  return iResult;
}
