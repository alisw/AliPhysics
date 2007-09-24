// @(#) $Id$

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

/** @file   AliHLTRootFileWriterComponent.cxx
    @author Matthias Richter
    @date   
    @brief  Base class for writer components to store data in a ROOT file

                                                                          */

#include "AliHLTRootFileWriterComponent.h"
#include "TFile.h"
#include "TString.h"

/** the global object for component registration */
AliHLTRootFileWriterComponent gAliHLTRootFileWriterComponent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTRootFileWriterComponent)

AliHLTRootFileWriterComponent::AliHLTRootFileWriterComponent()
  :
  AliHLTFileWriter(),
  fEventID(kAliHLTVoidEventID),
  fCurrentFile(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  // all blocks of one event go into the same file
  SetMode(kConcatenateBlocks);
}

AliHLTRootFileWriterComponent::AliHLTRootFileWriterComponent(const AliHLTRootFileWriterComponent&)
  :
  AliHLTFileWriter(),
  fEventID(kAliHLTVoidEventID),
  fCurrentFile(NULL)
{
  // see header file for class documentation
  HLTFatal("copy constructor untested");
}

AliHLTRootFileWriterComponent& AliHLTRootFileWriterComponent::operator=(const AliHLTRootFileWriterComponent&)
{
  // see header file for class documentation
  HLTFatal("assignment operator untested");
  return *this;
}

AliHLTRootFileWriterComponent::~AliHLTRootFileWriterComponent()
{
  // see header file for class documentation
}

int AliHLTRootFileWriterComponent::CloseWriter()
{
  // see header file for class documentation
  if (fCurrentFile!=NULL) {
    HLTDebug("close root file");
    TFile* pFile=fCurrentFile; fCurrentFile=NULL;
    pFile->Close(); delete pFile;
  }
  return 0;
}

int AliHLTRootFileWriterComponent::DumpEvent( const AliHLTComponentEventData& evtData,
					    const AliHLTComponentBlockData* blocks, 
					    AliHLTComponentTriggerData& trigData )
{
  // see header file for class documentation
  int iResult=0;
  if (evtData.fStructSize==0 && blocks==NULL && trigData.fStructSize==0) {
    // this is just to get rid of the warning "unused parameter"
  }
  const TObject* pObj=GetFirstInputObject(kAliHLTAnyDataType);
  HLTDebug("got first object %p", pObj);
  int count=0;
  while (pObj && iResult>=0) {
    iResult=WriteObject(evtData.fEventID, pObj);
    if (iResult) {
      count++;
      HLTDebug("wrote object of class %s, data type %s", pObj->ClassName(), (DataType2Text(GetDataType(pObj)).c_str())); 
    }
    pObj=GetNextInputObject();
  }
  HLTDebug("wrote %d object(s) from %d input blocks to file", count, GetNumberOfInputBlocks());
  return iResult;
}

int AliHLTRootFileWriterComponent::ScanArgument(int argc, const char** argv)
{
  // see header file for class documentation
  // no other arguments known
  if (argc==0 && argv==NULL) {
    // this is just to get rid of the warning "unused parameter"
  }
  int iResult=-EINVAL;
  return iResult;
}

int AliHLTRootFileWriterComponent::WriteObject(const AliHLTEventID_t eventID, const TObject *pOb)
{
  // see header file for class documentation
  int iResult=0;
  if (pOb) {
    HLTDebug("write object %p (%s)", pOb, pOb->GetName());
    if (CheckMode(kConcatenateEvents) && eventID != fEventID &&
	fCurrentFile!=NULL && eventID!=kAliHLTVoidEventID) {
      TFile* pFile=fCurrentFile; fCurrentFile=NULL;
      pFile->Close(); delete pFile;
    }
    if (fCurrentFile==NULL) {
      fCurrentFile=OpenFile(eventID, 0);
      if (fCurrentFile) fEventID=eventID;
    }
    if (fCurrentFile) {
      fCurrentFile->cd();
      pOb->Write();
    } else {
      iResult=-EBADF;
    }
  }
  return iResult;
}

TFile* AliHLTRootFileWriterComponent::OpenFile(const AliHLTEventID_t eventID, const int blockID, const char* option)
{
  // see header file for class documentation
  TFile* pFile=NULL;
  TString filename("");
  if ((BuildFileName(eventID, blockID, kAliHLTVoidDataType, 0, filename))>=0 && filename.IsNull()==0) {
    pFile=new TFile(filename, option);
    if (pFile) {
      if (pFile->IsZombie()) {
	delete pFile;
	pFile=NULL;
	HLTError("can not open ROOT file %s", filename.Data());
      }
    } else {
      HLTFatal("memory allocation failed");
    }
  } else {
    HLTError("failed to build a new file name for event %#8x", eventID);
  }
  return pFile;
}
