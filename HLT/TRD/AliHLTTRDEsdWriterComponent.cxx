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
#include "AliHLTTRDUtils.h"			\

/** global instance for component registration */
AliHLTTRDEsdWriterComponent gTRDEsdWriter;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTRDEsdWriterComponent)

AliHLTTRDEsdWriterComponent::AliHLTTRDEsdWriterComponent()
  :
  AliHLTRootFileWriterComponent(),
  fTree(NULL),
  fOutputPercentage(100),
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
  AliHLTRootFileWriterComponent(),
  fTree(NULL),
  fOutputPercentage(100),
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
					    const AliHLTComponentBlockData*  blocks, 
					    AliHLTComponentTriggerData& /*trigData*/ )
{
int result=0;

for ( unsigned long iBlock = 0; iBlock < evtData.fBlockCnt; iBlock++ )
    {
    // HLTDebug("i am a debug message"); // y does it not print out debug messages???
    /*HLTInfo("Block # %i/%i; Event 0x%08LX (%Lu)",
		    iBlock, evtData.fBlockCnt,
		    evtData.fEventID, evtData.fEventID);*/
    

    TClonesArray* tracksArray = NULL;
    const AliHLTComponentBlockData &block = blocks[iBlock];
    tracksArray = new TClonesArray("AliTRDtrackV1");
    
    //HLTInfo("BLOCK fPtr 0x%x, fOffset %i, fSize %i, fSpec 0x%x, fDataType %s", block.fPtr, block.fOffset, block.fSize, block.fSpecification, DataType2Text(block.fDataType).c_str()); //HLTInfo instead of HLTDebug, because debug gives no output... -> strange

    AliHLTTRDUtils::ReadTracks(tracksArray, block.fPtr, block.fSize); 

    // give out number of tracklets in tracksArray
    Int_t nbEntries = tracksArray->GetEntries();
    HLTInfo(" %i TRDtracks in tracksArray", nbEntries); 
    
    }

 /*  AliESDtrack* track = (AliESDtrack *)tobjin;
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
  fESD->Reset(); */

  return result;
}

// int AliHLTTRDEsdWriterComponent::ScanArgument(int argc, const char** argv)
// {
//   // see header file for class documentation
//   int iResult=AliHLTRootFileWriterComponent::ScanArgument(argc, argv);
//   return iResult;
// }

int AliHLTTRDEsdWriterComponent::DoEvent(	const AliHLTComponent_EventData& /*evtData*/,
						const AliHLTComponent_BlockData* /*blocks*/,
						AliHLTComponent_TriggerData& /*trigData*/,
						AliHLTUInt8_t* /*outputPtr*/,
						AliHLTUInt32_t& /*size*/,
						vector<AliHLTComponent_BlockData>& /*outputBlocks*/)
{
HLTDebug("ignor me");
return 0;

}

Int_t AliHLTTRDEsdWriterComponent::ScanArgument( int argc, const char** argv )
{
  // perform initialization. We check whether our relative output size is specified in the arguments.
  int i = 0;
  char* cpErr;
  HLTDebug("argv[%d] == %s", i, argv[i] );
      if ( !strcmp( argv[i], "output_percentage" ) )
	{
	  if ( i+1>=argc )
	    {
	      HLTError("Missing output_percentage parameter");
	      return ENOTSUP;
	    }
	  HLTDebug("argv[%d+1] == %s", i, argv[i+1] );
	  fOutputPercentage = strtoul( argv[i+1], &cpErr, 0 );
	  if ( *cpErr )
	    {
	      HLTError("Cannot convert output_percentage parameter '%s'", argv[i+1] );
	      return EINVAL;
	    }
	  HLTInfo("Output percentage set to %lu %%", fOutputPercentage );
	
    }
    //AliHLTRootFileWriterComponent::ScanArgument(argc, argv);

  return 0;
}


void AliHLTTRDEsdWriterComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // Get the output data size
  constBase = 0;
  inputMultiplier = ((double)fOutputPercentage)/100.0;
}
