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

/** @file   AliHLTTRDEsdWriterComponent.cxx
    @author Mateusz Ploskon
    @date   
    @brief  Writer component to store tracks of the HLT TRD

                                                                          */
#include "AliHLTTRDEsdWriterComponent.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "TTree.h"
#include "AliHLTTRDDefinitions.h"

/** global instance for component registration */
AliHLTTRDEsdWriterComponent gTRDEsdWriter;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTRDEsdWriterComponent)

AliHLTTRDEsdWriterComponent::AliHLTTRDEsdWriterComponent()
  :
  fTree(NULL),
  fESD(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTRDEsdWriterComponent::AliHLTTRDEsdWriterComponent(const AliHLTTRDEsdWriterComponent&)
  :
  fTree(NULL),
  fESD(NULL)
{
}

AliHLTTRDEsdWriterComponent& AliHLTTRDEsdWriterComponent::operator=(const AliHLTTRDEsdWriterComponent&)
{
  return *this;
}

void AliHLTTRDEsdWriterComponent::GetInputDataTypes( vector<AliHLTComponent_DataType>& list)
{
  // Get the list of input data  
  list.clear(); // We do not have any requirements for our input data type(s).
  list.push_back( AliHLTTRDDefinitions::fgkTRDSATracksDataType );
}

AliHLTTRDEsdWriterComponent::~AliHLTTRDEsdWriterComponent()
{
  // see header file for class documentation
}

int AliHLTTRDEsdWriterComponent::InitWriter()
{
  // see header file for class documentation
  int iResult=0;
  fESD = new AliESDEvent;
  if (fESD) {
    fESD->CreateStdContent();
    fTree = new TTree("esdTree", "Tree with HLT ESD objects");
    if (fTree) {
      fESD->WriteToTree(fTree);
    }
  }
  if (fTree==NULL) {
    iResult=-ENOMEM;
  }
  return iResult;
}

int AliHLTTRDEsdWriterComponent::CloseWriter()
{
  // see header file for class documentation
  int iResult=0;
  if (fTree) {
    WriteObject(kAliHLTVoidEventID, fTree);
    TTree* pTree=fTree;
    fTree=NULL;
    delete pTree;
  } else {
    HLTWarning("not initialized");
  }
  iResult=AliHLTRootFileWriterComponent::CloseWriter();
  return iResult;
}

int AliHLTTRDEsdWriterComponent::DumpEvent( const AliHLTComponentEventData& evtData,
					    const AliHLTComponentBlockData* blocks, 
					    AliHLTComponentTriggerData& trigData )
{
  // see header file for class documentation
  int iResult=0;
  TTree* pTree=fTree;

  AliHLTUInt32_t fDblock_Specification = 0;

  //implement a usage of the following
  //   AliHLTUInt32_t triggerDataStructSize = trigData.fStructSize;
  //   AliHLTUInt32_t triggerDataSize = trigData.fDataSize;
  //   void *triggerData = trigData.fData;
  Logging( kHLTLogDebug, "HLT::TRDEsdWriter::DumpEvent", "Trigger data received", 
	   "Struct size %d Data size %d Data location 0x%x", trigData.fStructSize, trigData.fDataSize, (UInt_t*)trigData.fData);

  //AliHLTComponentBlockData *dblock = (AliHLTComponentBlockData *)GetFirstInputBlock( AliHLTTRDDefinitions::fgkTRDSAEsdDataType );
  AliHLTComponentBlockData *dblock = (AliHLTComponentBlockData *)GetFirstInputBlock( AliHLTTRDDefinitions::fgkTRDSATracksDataType );
  if (dblock != 0)
    {
      fDblock_Specification = dblock->fSpecification;
    }
  else
    {
      Logging( kHLTLogWarning, "HLT::TRDEsdWriter::DumpEvent", "DATAIN", "First Input Block not found! 0x%x", dblock);
      return -1;
    }

  int ibForce = 1;
  TObject *tobjin = (TObject *)GetFirstInputObject( AliHLTTRDDefinitions::fgkTRDSATracksDataType, "AliESDtrack", ibForce);
  Logging( kHLTLogInfo, "HLT::TRDEsdWriter::DumpEvent", "1stBLOCK", "Pointer = 0x%x", tobjin);

  AliESDtrack* track = (AliESDtrack *)tobjin;
  if (!track)
    {
      Logging( kHLTLogWarning, "HLT::TRDEsdWriter::DumpEvent", "DATAIN", "First Input Block not a ESDtrack! 0x%x", tobjin);
      return -1;
    }

  Int_t nTracks = 0;
  while (tobjin != 0)
    {
      if (track != 0)
	{
	  //Logging( kHLTLogInfo, "HLT::TRDEsdWriter::DumpEvent", "Track found", "0x%x", track);	  
	  Logging( kHLTLogInfo, "HLT::TRDEsdWriter::DumpEvent", "DONE", "Track %d 0x%x Pt %1.2f", nTracks, track, track->Pt());
	  fESD->AddTrack(track);
	  nTracks++;
	}

      track = 0;
      tobjin = 0;
      tobjin = (TObject *)GetNextInputObject( ibForce );
      //Logging( kHLTLogInfo, "HLT::TRDEsdWriter::DumpEvent", "nextBLOCK", "Pointer = 0x%x", tobjin);
      track = (AliESDtrack *)tobjin;
    }

  Logging( kHLTLogInfo, "HLT::TRDEsdWriter::DumpEvent", "Fill", "Ntracks: %d", nTracks);	  
  pTree->Fill();
  fESD->Reset();

  return iResult;
}

int AliHLTTRDEsdWriterComponent::ScanArgument(int argc, const char** argv)
{
  // see header file for class documentation
  int iResult=AliHLTRootFileWriterComponent::ScanArgument(argc, argv);
  return iResult;
}
