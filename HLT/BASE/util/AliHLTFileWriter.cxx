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

/** @file   AliHLTFileWriter.cxx
    @author Matthias Richter
    @date   
    @brief  HLT file writer component implementation. */

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTFileWriter.h"
#include "AliHLTBlockDataCollection.h"
#include <TObjArray.h>
#include <TObjString.h>
#include <TSystem.h>
//#include <TMath.h>
//#include <TFile.h>
#include <cassert>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTFileWriter)

AliHLTFileWriter::AliHLTFileWriter()
  :
  AliHLTDataSink(),
  fpBlockDataCollection(NULL),
  fBaseName(""),
  fExtension(""),
  fDirectory(""),
  fSubDirFormat(""),
  fIdFormat("_0x%08x"),
  fSpecFormat(""),
  fBlcknoFormat("_0x%02x"),
  fCurrentFileName(""),
  fMode(0)
  , fpBurstBuffer(NULL)
  , fBurstBufferSize(0)
  , fBurstBlocks()
  , fBurstBlockEvents()
  , fPublisherConfName()
  , fPublisherConfEvent(-1)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTFileWriter::~AliHLTFileWriter()
{
  // see header file for class documentation

  // file list and file name list are owner of their objects and
  // delete all the objects
}

int AliHLTFileWriter::SetDefaults()
{
  // see header file for class documentation
  fBaseName="";
  fExtension="";
  fDirectory="";
  fSubDirFormat="";
  fIdFormat="_0x%08x";
  fSpecFormat="";
  fBlcknoFormat="_0x%02x";
  fCurrentFileName="";
  fMode=0;
  fpBurstBuffer=NULL;
  fBurstBufferSize=0;
  fBurstBlocks.clear();
  fBurstBlockEvents.clear();
  fPublisherConfEvent=-1;
  return 0;
}

const char* AliHLTFileWriter::GetComponentID()
{
  // see header file for class documentation
  return "FileWriter";
}

void AliHLTFileWriter::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back(kAliHLTAllDataTypes);
}

AliHLTComponent* AliHLTFileWriter::Spawn()
{
  // see header file for class documentation
  return new AliHLTFileWriter;
}

int AliHLTFileWriter::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;
  fpBlockDataCollection=new AliHLTBlockDataCollection;
  TString argument="";
  int bMissingParam=0;
  char* cpErr=NULL;
  int i=0;
  for (; i<argc && iResult>=0; i++) {
    cpErr=NULL;
    argument=argv[i];
    if (argument.IsNull()) continue;

    // -basename
    if (argument.CompareTo("-datafile")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fBaseName=argv[i];
      TObjArray* pTokens=fBaseName.Tokenize(".");
      if (pTokens) {
	int iEntries=pTokens->GetEntriesFast();
	if (iEntries>1) {
	  int n=0;
	  fBaseName=((TObjString*)pTokens->At(n++))->GetString();
	  while (n<iEntries-1) {
	    fBaseName+="." + ((TObjString*)pTokens->At(n++))->GetString();
	  }
	  fExtension=((TObjString*)pTokens->At(n))->GetString();
	}
	delete pTokens;
      }

      // -directory
    } else if (argument.CompareTo("-directory")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fDirectory=argv[i];

      // -subdir
    } else if (argument.BeginsWith("-subdir")) {
      argument.ReplaceAll("-subdir", "");
      if (argument.BeginsWith("=")) {
	fSubDirFormat=argument.Replace(0,1,"");
	if (strchr(fSubDirFormat.Data(), '%')==NULL) {
	  fSubDirFormat+="%lu";
	}
      } else {
	fSubDirFormat="event%03lu";
      }
      // no additional eventno in the filename unless set again
      // the sub dir contains the id
      fIdFormat="";

      // -idfmt
    } else if (argument.BeginsWith("-idfmt")) {
      argument.ReplaceAll("-idfmt", "");
      if (argument.BeginsWith("=")) {
	fIdFormat=argument.Replace(0,1,"");
      }

      // -specfmt
    } else if (argument.BeginsWith("-specfmt")) {
      argument.ReplaceAll("-specfmt", "");
      if (argument.BeginsWith("=")) {
	fSpecFormat=argument.Replace(0,1,"");
      } else {
	fSpecFormat="_0x%08x";
      }

      // -blocknofmt
    } else if (argument.BeginsWith("-blcknofmt") || 
	       argument.BeginsWith("-blocknofmt")) {
      // for the sake of backward compatibility we consider also the
      // old argument with typo for a while
      argument.ReplaceAll("-blcknofmt", "");
      argument.ReplaceAll("-blocknofmt", "");
      if (argument.BeginsWith("=")) {
	fBlcknoFormat=argument.Replace(0,1,"");
      } else {
	fBlcknoFormat="_0x%02x";
      }

      // -enumeration
    } else if (argument.CompareTo("-enumerate")==0) {
      SetMode(kEnumerate);

      // -concatenate-blocks
    } else if (argument.CompareTo("-concatenate-blocks")==0) {
      SetMode(kConcatenateBlocks);

      // -concatenate-events
    } else if (argument.CompareTo("-concatenate-events")==0) {
      SetMode(kConcatenateEvents);

      // -publisher-conf
    } else if (argument.CompareTo("-publisher-conf")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fPublisherConfName=argv[i];

      // -write-all-events
    } else if (argument.CompareTo("-write-all-events")==0) {
      SetMode(kWriteAllEvents);

      // -write-all-blocks
    } else if (argument.CompareTo("-write-all-blocks")==0) {
      SetMode(kWriteAllBlocks);

      // -write-all
    } else if (argument.CompareTo("-write-all")==0) {
      SetMode(kWriteAllEvents);
      SetMode(kWriteAllBlocks);

      // -burst-buffer
    } else if (argument.CompareTo("-burst-buffer")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fBurstBufferSize = strtoul( argv[i], &cpErr ,0);
      if ( *cpErr ) break;

      // -skip-datatype
    } else if(argument.CompareTo("-skip-datatype")==0){
      SetMode(kSkipDataType);

      // check for selection arguments (AliHLTBlockDataCollection)
    } else if (fpBlockDataCollection && 
	       (iResult=fpBlockDataCollection->ScanArgument(argc-i, &argv[i]))>0) {
	i+=iResult-1;
	iResult=0;
    } else {
      if ((iResult=ScanArgument(argc-i, &argv[i]))==-EINVAL) {
	HLTError("unknown argument %s", argument.Data());
	break;
      } else if (iResult==-EPROTO) {
	bMissingParam=1;
	break;
      } else if (iResult>0) {
	i+=iResult-1;
	iResult=0;
      }
    }
  }

  if (cpErr && *cpErr) {
    HLTError("Cannot convert specifier '%s' for argument '%s'", argv[i], argument.Data());
    iResult=-EINVAL;
  } else if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  if (fpBlockDataCollection &&
      (iResult<0 || fpBlockDataCollection->IsEmpty())) {
    delete fpBlockDataCollection;
    fpBlockDataCollection=NULL;
  }
  if (iResult>=0 && fBurstBufferSize>0) {
    if (!CheckMode(kConcatenateBlocks) || !CheckMode(kConcatenateEvents)) {
      HLTError("burst write currently only supported for mode kConcatenateBlocks AND kConcatenateEvents");
      iResult=-EINVAL;
    } else {
    fpBurstBuffer=new AliHLTUInt8_t[fBurstBufferSize];
    if (!fpBurstBuffer) {
      iResult=-ENOMEM;
      fBurstBufferSize=0;
    }
    }
  }

  if (!fPublisherConfName.IsNull()) {
    if (CheckMode(kConcatenateBlocks) || CheckMode(kConcatenateEvents)) {
      fPublisherConfName="";
      HLTWarning("option 'concatenate blocks/events' collides with writing of FilePublisher configuration, ignoring option '-publisher-conf'");
    }
  }

  if (iResult>=0) {
    iResult=InitWriter();
    if (!fDirectory.IsNull()) {
      gSystem->mkdir(fDirectory);
    }
  }

  return iResult;
}

int AliHLTFileWriter::InitWriter()
{
  // see header file for class documentation
  
  // fCurrentFileName is used in dump event, just touched her to avoid
  // coding convention violation RC11. The function can not be declared
  // const since it is just the default implementation, overloaded
  // virtual function might not be const
  fCurrentFileName="";
  return 0; // note: this doesn't mean 'error'
}

int AliHLTFileWriter::ScanArgument(int /*argc*/, const char** /*argv*/)
{
  // see header file for class documentation

  // there are no other arguments than the standard ones
  // fCurrentFileName is used in dump event, just touched her to avoid
  // coding convention violation RC11. The function can not be declared
  // const since it is just the default implementation, overloaded
  // virtual function might not be const
  fCurrentFileName="";
  return -EINVAL;
}

int AliHLTFileWriter::DoDeinit()
{
  // see header file for class documentation
  int iResult=0;
  if (fpBurstBuffer) {
    if ((iResult=BurstWrite())<0) {
      HLTError("failed BurstWrite");
    }
    delete [] fpBurstBuffer;
    fpBurstBuffer=NULL;
    fBurstBufferSize=0;
    fBurstBlocks.clear();
    fBurstBlockEvents.clear();
  }

  iResult=CloseWriter();
  ClearMode(kEnumerate);

  if (fpBlockDataCollection) delete fpBlockDataCollection;
  fpBlockDataCollection=NULL;

  SetDefaults();
  return iResult;
}

int AliHLTFileWriter::CloseWriter()
{
  // see header file for class documentation

  // fCurrentFileName is used in dump event, just touched her to avoid
  // coding convention violation RC11. The function can not be declared
  // const since it is just the default implementation, overloaded
  // virtual function might not be const
  fCurrentFileName="";
  return 0;
}

int AliHLTFileWriter::DumpEvent( const AliHLTComponentEventData& evtData,
			 AliHLTComponentTriggerData& /*trigData*/ )
{
  // see header file for class documentation
  int iResult=0;
  if (!IsDataEvent() && !CheckMode(kWriteAllEvents)) return 0;

  if (CheckMode(kConcatenateEvents)==0) {
    // reset the current file name in order to open a new file
    // for the first block. If events are concatenated, the current
    // file name stays in order to be opended in append mode.
    fCurrentFileName="";
  }
  const AliHLTComponentBlockData* pDesc=NULL;

  int blockno=0;
  for (pDesc=GetFirstInputBlock(); pDesc!=NULL; pDesc=GetNextInputBlock(), blockno++) {
    if (fpBlockDataCollection) {
      if (!fpBlockDataCollection->IsSelected(*pDesc)) continue;
    } else if (pDesc->fDataType==(kAliHLTAnyDataType|kAliHLTDataOriginPrivate) && !CheckMode(kWriteAllBlocks))
      continue;
    HLTDebug("block %d out of %d", blockno, evtData.fBlockCnt);
    iResult=ScheduleBlock(blockno, evtData.fEventID, pDesc);
  }
  return iResult;
}

int AliHLTFileWriter::BuildFileName(const AliHLTEventID_t eventID, const int blockID,
				    const AliHLTComponentDataType& dataType,
				    const AliHLTUInt32_t specification,
				    TString& filename)
{
  // see header file for class documentation
  int iResult=0;
  //HLTDebug("build file name for event %d block %d", eventID, blockID);
  filename="";

  AliHLTUInt32_t eventType=gkAliEventTypeUnknown;

  if (!fDirectory.IsNull()) {
    filename+=fDirectory;
    if (!filename.EndsWith("/"))
      filename+="/";
  }
  if (!fSubDirFormat.IsNull()) {
    filename+=Form(fSubDirFormat, eventID);
    if (!filename.EndsWith("/"))
      filename+="/";
  }
  if (filename.EndsWith("/")) {
    gSystem->mkdir(filename);
  }

  IsDataEvent(&eventType);
  if (CheckMode(kWriteAllEvents) && !CheckMode(kConcatenateEvents)) {
    if (eventType==gkAliEventTypeStartOfRun) filename+="SOR_";
    else if (eventType==gkAliEventTypeEndOfRun) filename+="EOR_";
  }

  if (!fBaseName.IsNull())
    filename+=fBaseName;
  else
    filename+="event";
  if (!CheckMode(kConcatenateEvents)) {
    if (!CheckMode(kEnumerate)) {
      if (eventID!=kAliHLTVoidEventID && !fIdFormat.IsNull()) {
	filename+=Form(fIdFormat, eventID);
      }
    } else {
      filename+=Form("_%d", GetEventCount());
    }
  }
  if (blockID>=0 && !CheckMode(kConcatenateBlocks)) {
    if (!fBlcknoFormat.IsNull())
      filename+=Form(fBlcknoFormat, blockID);
    if (dataType!=kAliHLTVoidDataType &&
	!CheckMode(kSkipDataType)) {
      filename+="_";
      filename+=AliHLTComponent::DataType2Text(dataType).data();
    }
    if (!fSpecFormat.IsNull())
      filename+=Form(fSpecFormat, specification);
  }
  if (!fExtension.IsNull())
    filename+="." + fExtension;
  filename.ReplaceAll(" ", "");
  return iResult;
}

int AliHLTFileWriter::SetMode(Short_t mode) 
{
  // see header file for class documentation
  fMode|=mode;
  //HLTDebug("mode set to 0x%x", fMode);
  return fMode;
}

int AliHLTFileWriter::ClearMode(Short_t mode)
{
  // see header file for class documentation
  fMode&=~mode;
  //HLTDebug("mode set to 0x%x", fMode);
  return fMode;
}

int AliHLTFileWriter::CheckMode(Short_t mode) const
{
  // see header file for class documentation

  //HLTDebug("check mode 0x%x for flag 0x%x: %d", fMode, mode, (fMode&mode)!=0);
  return (fMode&mode)!=0;
}

int AliHLTFileWriter::ScheduleBlock(int blockno, const AliHLTEventID_t& eventID,
				    const AliHLTComponentBlockData* pDesc)
{
  // see header file for class documentation
  int iResult=0;
  if (fpBurstBuffer==NULL ||
      (fBurstBlocks.size()==0 && pDesc->fSize>fBurstBufferSize) ) {
    return WriteBlock(blockno, eventID, pDesc);
  }
  AliHLTComponentBlockData bd=*pDesc;
  bd.fPtr=NULL;
  if (fBurstBlocks.size()>0) {
    bd.fOffset=fBurstBlocks.back().fOffset+fBurstBlocks.back().fSize;
  } else {
    bd.fOffset=0;
  }
  if (bd.fOffset+bd.fSize>fBurstBufferSize) {
    if ((iResult=BurstWrite())>=0) {
      iResult=WriteBlock(blockno, eventID, pDesc);
    }
  } else {
    memcpy(fpBurstBuffer+bd.fOffset, pDesc->fPtr, bd.fSize);
    fBurstBlocks.push_back(bd);
    fBurstBlockEvents.push_back(eventID);
  }

  return iResult;
}

int AliHLTFileWriter::BurstWrite()
{
  // see header file for class documentation
  int iResult=0;
  if (fBurstBlocks.size()==0) return 0;
  assert(fBurstBlocks.size()==fBurstBlockEvents.size());
  HLTDebug("writing %d postponed blocks", fBurstBlocks.size());
  int blockno=0;
  AliHLTComponentBlockDataList::iterator block=fBurstBlocks.begin();
  AliHLTComponentBlockDataList::iterator firstBlock=block;
  vector<AliHLTEventID_t>::iterator event=fBurstBlockEvents.begin();
  if (CheckMode(kConcatenateEvents)) {
    block=fBurstBlocks.end()-1;
    event=fBurstBlockEvents.end()-1;
  }
  for (; block!=fBurstBlocks.end() && iResult>=0; block++, event++, blockno++) {
    if (event!=fBurstBlockEvents.begin() && *event!=*(event-1)) {
      blockno=0;
    }
    if (CheckMode(kConcatenateEvents)) {
      // all blocks in the burst buffer are written in one go
      // just the block descriptor is updated appropriately
      (*block).fSize+=(*block).fOffset;
      (*block).fPtr=fpBurstBuffer;
    } else if (CheckMode(kConcatenateBlocks)) {
      // all blocks of the same event are written in one go
      // just the block descriptor is updated appropriately
      if (event+1==fBurstBlockEvents.end() ||
	  *event!=*(event+1)) {
	(*block).fSize+=(*block).fOffset-(*firstBlock).fOffset;
	(*block).fPtr=fpBurstBuffer+(*firstBlock).fOffset;
	firstBlock=block+1;
      } else {
	// coninue if it wasn't the last block of this event
	continue;
      }
    } else {
      (*block).fPtr=fpBurstBuffer+(*block).fOffset;
    }
    (*block).fOffset=0;
    iResult=WriteBlock(blockno, *event, &(*block));
  }
  fBurstBlocks.clear();
  fBurstBlockEvents.clear();

  return iResult;
}

int AliHLTFileWriter::WriteBlock(int blockno, const AliHLTEventID_t& eventID,
				 const AliHLTComponentBlockData* pDesc)
{
  // see header file for class documentation
  int iResult=0;
  TString filename;
  HLTDebug("dataspec 0x%x", pDesc->fSpecification);
  iResult=BuildFileName(eventID, blockno, pDesc->fDataType, pDesc->fSpecification, filename);
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
      dump.write((static_cast<const char*>(pDesc->fPtr)), pDesc->fSize);
      HLTDebug("wrote %d byte(s) to file %s", pDesc->fSize, filename.Data());
    } else {
      HLTError("can not open file %s for writing", filename.Data());
      iResult=-EBADF;
    }
    dump.close();
  }
  if (iResult>=0 && !fPublisherConfName.IsNull()) {
    if (!CheckMode(kConcatenateBlocks) &&
	!CheckMode(kConcatenateEvents)) {
      // append if not the first entry
      if (fPublisherConfEvent>=0) filemode=ios::app;
      else filemode=(ios::openmode)0;
      ofstream conf(fPublisherConfName.Data(), filemode);
      if (conf.good()) {
	if (fPublisherConfEvent>=0 &&
	    fPublisherConfEvent!=GetEventCount()) {
	  conf << "-nextevent " << endl;
	}
	fPublisherConfEvent=GetEventCount();
	conf << "-datatype ";
	conf << DataType2Text(pDesc->fDataType, 3);
	conf << " -dataspec ";
	TString specstr; specstr.Form("0x%08x", pDesc->fSpecification);
	conf << specstr;
	conf << " -datafile ";
	conf << filename;
	conf << endl;
      } else {
	fPublisherConfName="";
	HLTError("can not open file %s for writing of configuration commands", fPublisherConfName.Data());
      }
      conf.close();
    } else {
	fPublisherConfName="";
	HLTWarning("option 'concatenate blocks/events' collides with writing of FilePublisher configuration, disable ...");
    }
  }
  return iResult;
}
