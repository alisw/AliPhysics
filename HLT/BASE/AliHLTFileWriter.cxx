// $Id$

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Matthias Richter <Matthias.Richter@ift.uib.no>                *
 *          for The ALICE Off-line Project.                               *
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
#include <TObjArray.h>
#include <TObjString.h>
#include <TMath.h>
#include <TFile.h>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTFileWriter)

AliHLTFileWriter::AliHLTFileWriter()
  :
  fBaseName(""),
  fExtension(""),
  fDirectory(""),
  fCurrentFileName(""),
  fMode(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTFileWriter::AliHLTFileWriter(const AliHLTFileWriter&)
  :
  fBaseName(""),
  fExtension(""),
  fDirectory(""),
  fCurrentFileName(""),
  fMode(0)
{
  // see header file for class documentation
  HLTFatal("copy constructor untested");
}

AliHLTFileWriter& AliHLTFileWriter::operator=(const AliHLTFileWriter&)
{ 
  // see header file for class documentation
  HLTFatal("assignment operator untested");
  return *this;
}

AliHLTFileWriter::~AliHLTFileWriter()
{
  // see header file for class documentation

  // file list and file name list are owner of their objects and
  // delete all the objects
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
  list.push_back(kAliHLTAnyDataType);
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
  TString argument="";
  int bMissingParam=0;
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

    // -basename
    if (argument.CompareTo("-datafile")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fBaseName=argv[i];
      TObjArray* pTokens=fBaseName.Tokenize(".");
      if (pTokens) {
	int iEntries=pTokens->GetEntries();
	if (iEntries>1) {
	  int i=0;
	  fBaseName=((TObjString*)pTokens->At(i++))->GetString();
	  while (i<iEntries-1) {
	    fBaseName+="." + ((TObjString*)pTokens->At(i++))->GetString();
	  }
	  fExtension=((TObjString*)pTokens->At(i))->GetString();
	}
	delete pTokens;
      }

      // -directory
    } else if (argument.CompareTo("-directory")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fDirectory=argv[i];

      // -enumeration
    } else if (argument.CompareTo("-enumerate")==0) {
      SetMode(kEnumerate);

      // -concatenate-blocks
    } else if (argument.CompareTo("-concatenate-blocks")==0) {
      SetMode(kConcatenateBlocks);

      // -concatenate-events
    } else if (argument.CompareTo("-concatenate-events")==0) {
      SetMode(kConcatenateEvents);

    } else {
      if ((iResult=ScanArgument(argc-i, &argv[i]))==-EINVAL) {
	HLTError("unknown argument %s", argument.Data());
	break;
      } else if (iResult==-EPROTO) {
	bMissingParam=1;
	break;
      } else if (iResult>=0) {
	i+=iResult;
	iResult=0;
      }
    }
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  if (iResult>=0) {
    iResult=InitWriter();
  }

  return iResult;
}

int AliHLTFileWriter::InitWriter()
{
  // see header file for class documentation
  return 0; // note: this doesn't mean 'error'
}

int AliHLTFileWriter::ScanArgument(int argc, const char** argv)
{
  // see header file for class documentation

  // there are no other arguments than the standard ones
  return -EINVAL;
}

int AliHLTFileWriter::DoDeinit()
{
  // see header file for class documentation
  int iResult=CloseWriter();
  ClearMode(kEnumerate);
  return iResult;
}

int AliHLTFileWriter::CloseWriter()
{
  // see header file for class documentation
  return 0; // note: this doesn't mean 'error'
}

int AliHLTFileWriter::DumpEvent( const AliHLTComponentEventData& evtData,
			 const AliHLTComponentBlockData* blocks, 
			 AliHLTComponentTriggerData& trigData )
{
  // see header file for class documentation
  int iResult=0;
  if (CheckMode(kConcatenateEvents)==0) {
    // reset the current file name in order to open a new file
    // for the first block. If events are concatenated, the current
    // file name stays in order to be opended in append mode.
    fCurrentFileName="";
  }
  for (int n=0; n<(int)evtData.fBlockCnt; n++ ) {
    //HLTDebug("block %d out of %d", n, evtData.fBlockCnt);
    TString filename;
    HLTDebug("dataspec 0x%x", blocks[n].fSpecification);
    iResult=BuildFileName(evtData.fEventID, n, blocks[n].fDataType, filename);
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
	dump.write((const char*)blocks[n].fPtr, blocks[n].fSize);
	HLTDebug("wrote %d byte(s) to file %s", blocks[n].fSize, filename.Data());
      } else {
	HLTError("can not open file %s for writing", filename.Data());
	iResult=-EBADF;
      }
      dump.close();
    }
  }
  return iResult;
}

int AliHLTFileWriter::BuildFileName(const AliHLTEventID_t eventID, const int blockID, const AliHLTComponentDataType& dataType, TString& filename)
{
  // see header file for class documentation
  int iResult=0;
  //HLTDebug("build file name for event %d block %d", eventID, blockID);
  filename="";
  if (!fDirectory.IsNull()) {
    filename+=fDirectory;
    if (!fDirectory.EndsWith("/"))
      filename+="/";
  }
  if (!fBaseName.IsNull())
    filename+=fBaseName;
  else
    filename+="event";
  if (!CheckMode(kConcatenateEvents)) {
    if (!CheckMode(kEnumerate)) {
      if (eventID!=kAliHLTVoidEventID) {
	filename+=Form("_0x%08x", eventID);
      }
    } else {
      filename+=Form("_%d", GetEventCount());
    }
  }
  if (blockID>=0 && !CheckMode(kConcatenateBlocks)) {
    filename+=Form("_0x%x", blockID);
    if (dataType!=kAliHLTVoidDataType) {
      filename+="_";
      filename+=AliHLTComponent::DataType2Text(dataType).data();
    }
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

int AliHLTFileWriter::CheckMode(Short_t mode)
{
  // see header file for class documentation

  //HLTDebug("check mode 0x%x for flag 0x%x: %d", fMode, mode, (fMode&mode)!=0);
  return (fMode&mode)!=0;
}
