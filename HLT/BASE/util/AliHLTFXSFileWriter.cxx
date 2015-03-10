// $Id$

///**************************************************************************
///* This file is property of and copyright by the                          * 
///* ALICE Experiment at CERN, All rights reserved.                         *
///*                                                                        *
///* Permission to use, copy, modify and distribute this software and its   *
///* documentation strictly for non-commercial purposes is hereby granted   *
///* without fee, provided that the above copyright notice appears in all   *
///* copies and that both the copyright notice and this permission notice   *
///* appear in the supporting documentation. The authors make no claims     *
///* about the suitability of this software for any purpose. It is          *
///* provided "as is" without express or implied warranty.                  *
///**************************************************************************/

/// @file   AliHLTFXSFileWriter.cxx
/// @author Timo Breitner
/// @date   
/// @brief  HLT FXS file writer component implementation.
#include <fstream>
#include <TSystem.h>
#include "AliHLTFXSFileWriter.h"

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTFXSFileWriter)

const Int_t AliHLTFXSFileWriter::gkCalibObjectMaxVersion=65000;

AliHLTFXSFileWriter::AliHLTFXSFileWriter()
  : AliHLTDataSink()
  , fDirectory("FXS")
{
  // An HLT data sink component which writes FXS data to file(s).
  //
  // Component ID: \b FXSFileWriter      <br>
  // Library: \b libAliHLTUtil.so     <br>
  // Input Data Types: ::kAliHLTDataTypeFXSCalib <br>
  // Output Data Types: none <br>
}

AliHLTFXSFileWriter::~AliHLTFXSFileWriter()
{
  // destructor
}

const char* AliHLTFXSFileWriter::GetComponentID()
{
  // overloaded from AliHLTComponent
  return "FXSFileWriter";
}

void AliHLTFXSFileWriter::GetInputDataTypes( AliHLTComponentDataTypeList& list)
{
  // overloaded from AliHLTComponent
  list.clear();
  list.push_back(kAliHLTDataTypeFXSCalib);
}

AliHLTComponent* AliHLTFXSFileWriter::Spawn()
{
  // overloaded from AliHLTComponent
  return new AliHLTFXSFileWriter;
}

int AliHLTFXSFileWriter::DoInit( int argc, const char** argv )
{
  // overloaded from AliHLTComponent: initialization
  TString argument="";
  int bMissingParam=0;
  int i=0;
  for (; i<argc; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

    if (argument.CompareTo("-directory")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fDirectory=argv[i];
    }
  }

  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    return -EINVAL;
  }

  if (!fDirectory.IsNull()) {
    if ( gSystem->mkdir(fDirectory) != 0 )
      return -EPERM;
  }

  return 0;
}

int AliHLTFXSFileWriter::DumpEvent( const AliHLTComponentEventData& evtData,
			 AliHLTComponentTriggerData& /*trigData*/ )
{
  // overloaded from AliHLTDataSink: event processing
  int iResult=0;

  const AliHLTComponentBlockData* pDesc=NULL;

  int blockno=0;
  for (pDesc=GetFirstInputBlock(); pDesc!=NULL; pDesc=GetNextInputBlock(), blockno++) {
    if (pDesc->fDataType!=(kAliHLTDataTypeFXSCalib|kAliHLTDataOriginAny))
      continue;
    HLTDebug("block %d out of %d", blockno, evtData.fBlockCnt);
    if(pDesc->fSize <= sizeof(AliHLTFXSHeader)){
      HLTError("Block size (%d) not larger than FXS header size (%d).",pDesc->fSize, sizeof(AliHLTFXSHeader));
      return -EINVAL;
    }
    AliHLTFXSHeader* header=reinterpret_cast<AliHLTFXSHeader*>(pDesc->fPtr);
    const char* data=static_cast<const char*>(pDesc->fPtr) + sizeof(AliHLTFXSHeader);
    int datasize=pDesc->fSize - sizeof(AliHLTFXSHeader);

    TString runNumber;
    runNumber += header->fRunNumber;
    TString detector = header->fOrigin;
    TString fileID = header->fFileID;
    TString ddlNumber;
    for (int ddl=0; ddl<gkAliHLTFXSHeaderfDDLNumberSize; ++ddl) {
      ddlNumber += header->fDDLNumber[ddl];
    }
    TString sep="/";

    TString basedir = fDirectory + sep + runNumber + sep + detector + sep + ddlNumber;
    TString basefile = basedir + sep + fileID;

    gSystem->mkdir(basedir, kTRUE);

    TString targetfile;

    for(int version=0; version<gkCalibObjectMaxVersion; ++version){
      TString tmpfile=basefile+ "_";
      tmpfile += version;
      ifstream test(tmpfile.Data());
      if(! test.good()){
	targetfile=tmpfile;
	break;
      }
      test.close();
    }

    if(targetfile.IsNull()){
      HLTError("Maximum version (%d) reached, cannot store any more objects.", gkCalibObjectMaxVersion);
      return -1;
    }

    ofstream dump(targetfile.Data(), ofstream::out );
    if (dump.good()) {
      dump.write(data, datasize);
      HLTDebug("wrote %d byte(s) to file %s", datasize, targetfile.Data());
    } else {
      HLTError("can not open file %s for writing", targetfile.Data());
      iResult=-EBADF;
    }
    dump.close();
  }
  return iResult;
}
