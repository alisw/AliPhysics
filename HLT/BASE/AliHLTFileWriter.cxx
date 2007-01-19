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
#include <TObjString.h>
#include <TMath.h>
#include <TFile.h>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTFileWriter)

AliHLTFileWriter::AliHLTFileWriter()
  :
  fBaseName(""),
  fDirectory(""),
  fEnumeration(""),
  fbSeparate(0)
{
}

AliHLTFileWriter::AliHLTFileWriter(const AliHLTFileWriter&)
  :
  fBaseName(""),
  fDirectory(""),
  fEnumeration(""),
  fbSeparate(0)
{
  HLTFatal("copy constructor untested");
}

AliHLTFileWriter& AliHLTFileWriter::operator=(const AliHLTFileWriter&)
{ 
  HLTFatal("assignment operator untested");
  return *this;
}

AliHLTFileWriter::~AliHLTFileWriter()
{
  // file list and file name list are owner of their objects and
  // delete all the objects
}

const char* AliHLTFileWriter::GetComponentID()
{
  return "FileWriter";
}

void AliHLTFileWriter::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  list.clear();
  list.push_back(kAliHLTAnyDataType);
}

AliHLTComponent* AliHLTFileWriter::Spawn()
{
  return new AliHLTFileWriter;
}

int AliHLTFileWriter::DoInit( int argc, const char** argv )
{
  int iResult=0;
  TString argument="";
  int bMissingParam=0;
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];

    // -basename
    if (argument.CompareTo("-datafile")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fBaseName=argv[i];

      // -directory
    } else if (argument.CompareTo("-directory")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fDirectory=argv[i];

      // -enumeration
    } else if (argument.CompareTo("-enumeration")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fEnumeration=argv[i];

      // -separate
    } else if (argument.CompareTo("-enumeration")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      TString parameter(argv[i]);
      fbSeparate=parameter.Atoi();
      if (fbSeparate < 0 || fbSeparate > 1) {
	HLTError("invalid parameter for argument %s", argument.Data());
      }

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
  return iResult;
}

int AliHLTFileWriter::ScanArgument(int argc, const char** argv)
{
  // there are no other arguments than the standard ones
  return -EINVAL;
}

int AliHLTFileWriter::DoDeinit()
{
  int iResult=0;
  return iResult;
}

int AliHLTFileWriter::DumpEvent( const AliHLTComponentEventData& evtData,
			 const AliHLTComponentBlockData* blocks, 
			 AliHLTComponentTriggerData& trigData )
{
  int iResult=0;
  for (int n=0; n<(int)evtData.fBlockCnt; n++ ) {
    //HLTDebug("block %d out of %d", n, evtData.fBlockCnt);
    TString filename;
    iResult=BuildFileName(evtData.fEventID, n, blocks[n].fDataType, filename);
    if (iResult>=0) {
      ofstream dump(filename.Data());
      if (dump.good()) {
	dump.write((const char*)blocks[n].fPtr, blocks[n].fSize);
	HLTDebug("wrote %d byte(s) to file %s", blocks[n].fSize, filename.Data());
      }
      dump.close();
    }
  }
  return iResult;
}

int AliHLTFileWriter::BuildFileName(const AliHLTEventID_t eventID, const int blockID, const AliHLTComponentDataType& dataType, TString& filename)
{
  int iResult=0;
  //HLTDebug("build file name for event %d block %d", eventID, blockID);
  filename.Form("event_%#08x_0x%x_", eventID, blockID);
  filename+=AliHLTComponent::DataType2Text(dataType).data();
  return iResult;
}
