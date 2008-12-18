// $Id$

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Matthias Richter <Matthias.Richter@ift.uib.no>                *
 *          Timm Steinbeck <timm@kip.uni-heidelberg.de>                   *
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

/** @file   AliHLTTRDCalibrationComponent.cxx
    @author Timm Steinbeck, Matthias Richter
    @date   
    @brief  A TRDCalibration processing component for the HLT. */

#if __GNUC__ >= 3
using namespace std;
#endif

#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"

#include "AliHLTTRDCalibrationComponent.h"
#include "AliHLTTRDDefinitions.h"
#include "AliHLTTRDTrack.h"

#include "AliCDBManager.h"
#include "AliTRDclusterizerHLT.h"
#include "AliRawReaderMemory.h"
#include "AliTRDCalibraFillHisto.h"

#include <cstdlib>
#include <cerrno>
#include <string>

// this is a global object used for automatic component registration, do not use this
AliHLTTRDCalibrationComponent gAliHLTTRDCalibrationComponent;

ClassImp(AliHLTTRDCalibrationComponent);
   
AliHLTTRDCalibrationComponent::AliHLTTRDCalibrationComponent():
  AliHLTCalibrationProcessor(),
  fTRDCalibraFillHisto(NULL),
  fUseHLTTracks(kFALSE),
  fOutputPercentage(100), // By default we copy to the output exactly what we got as input  
  fStrorageDBpath("local://$ALICE_ROOT"),
  fCDB(NULL)
{
  // Default constructor
}

AliHLTTRDCalibrationComponent::~AliHLTTRDCalibrationComponent()
{
  // Destructor
  ;
}

const char* AliHLTTRDCalibrationComponent::GetComponentID()
{
  // Return the component ID const char *
  return "TRDCalibration"; // The ID of this component
}

void AliHLTTRDCalibrationComponent::GetInputDataTypes( vector<AliHLTComponent_DataType>& list)
{
  // Get the list of input data
  list.clear(); // We do not have any requirements for our input data type(s).
  list.push_back( AliHLTTRDDefinitions::fgkTRDSATracksDataType );
}

AliHLTComponent_DataType AliHLTTRDCalibrationComponent::GetOutputDataType()
{
  // Get the output data type
  return AliHLTTRDDefinitions::fgkCalibrationDataType;
}

void AliHLTTRDCalibrationComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // Get the output data size
  constBase = 0;
  inputMultiplier = ((double)fOutputPercentage)/100.0;
}

AliHLTComponent* AliHLTTRDCalibrationComponent::Spawn()
{
  // Spawn function, return new instance of this class
  return new AliHLTTRDCalibrationComponent;
};

Int_t AliHLTTRDCalibrationComponent::ScanArgument( int argc, const char** argv )
{
  // perform initialization. We check whether our relative output size is specified in the arguments.
  fOutputPercentage = 100;
  int i = 0;
  char* cpErr;
  while ( i < argc )
    {
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
	  i += 2;
	  continue;
	}
      else if ( strcmp( argv[i], "-cdb" ) == 0)
	{
	  if ( i+1 >= argc )
	    {
	      HLTError("Missing -cdb argument");
	      return ENOTSUP;	      
	    }
	  fStrorageDBpath = argv[i+1];
	  HLTInfo("DB storage is %s", fStrorageDBpath.c_str() );	  
	  i += 2;
	  continue;
	}  
      else if ( strcmp( argv[i], "-useHLTTracks" ) == 0){
	fUseHLTTracks = kTRUE;
	i++;
	HLTInfo("Expecting block of AliHLTTracks as input. Using low-level interface");
      }
      else{
	HLTError("Unknown option '%s'", argv[i] );
	return EINVAL;
      }
    }
  return 0;
}

Int_t AliHLTTRDCalibrationComponent::InitCalibration()
{
  //init the calibration
  fCDB = AliCDBManager::Instance();
  if (!fCDB)
    {
      HLTError("Could not get CDB instance, fCDB 0x%x", fCDB);
    }
  else
    {
      fCDB->SetRun(0); // THIS HAS TO BE RETRIEVED !!!
      fCDB->SetDefaultStorage(fStrorageDBpath.c_str());
      HLTDebug("fCDB 0x%x", fCDB);
    }
  fTRDCalibraFillHisto = AliTRDCalibraFillHisto::Instance();
  fTRDCalibraFillHisto->SetHisto2d(); // choose to use histograms
  fTRDCalibraFillHisto->SetCH2dOn();  // choose to calibrate the gain
  fTRDCalibraFillHisto->SetPH2dOn();  // choose to calibrate the drift velocity
  fTRDCalibraFillHisto->SetPRF2dOn(); // choose to look at the PRF
  fTRDCalibraFillHisto->Init2Dhistos(); // initialise the histos
  return 0;
}

Int_t AliHLTTRDCalibrationComponent::DeinitCalibration()
{
  HLTDebug("DeinitCalibration");
  
  // Deinitialization of the component

  if (fCDB)
    {
      HLTDebug("destroy fCDB");
      fCDB->Destroy();
      fCDB = 0;
    }
  return 0;
}

Int_t AliHLTTRDCalibrationComponent::ProcessCalibration(const AliHLTComponent_EventData& evtData,
							const AliHLTComponent_BlockData* blocks,
							AliHLTComponent_TriggerData& /*trigData*/,
							AliHLTUInt8_t* /*outputPtr*/,
							AliHLTUInt32_t& /*size*/,
							vector<AliHLTComponent_BlockData>& /*outputBlocks*/)
{
  // Process an event
  //   Logging( kHLTLogInfo, "HLT::TRDCalibration::ProcessCalibration", "Output percentage set", "Output percentage set to %lu %%", fOutputPercentage );
  HLTDebug("NofBlocks %lu", evtData.fBlockCnt );
  // Process an event

  // Loop over all input blocks in the event
  vector<AliHLTComponent_DataType> expectedDataTypes;
  GetInputDataTypes(expectedDataTypes);
  for ( unsigned long iBlock = 0; iBlock < evtData.fBlockCnt; iBlock++ )
    {
      const AliHLTComponentBlockData &block = blocks[iBlock];
      AliHLTComponentDataType inputDataType = block.fDataType;
      Bool_t correctDataType = kFALSE;
      
      for(UInt_t i = 0; i < expectedDataTypes.size(); i++)
	if( expectedDataTypes.at(i) == inputDataType)
	  correctDataType = kTRUE;
      if (!correctDataType)
	{
	  HLTDebug( "Block # %i/%i; Event 0x%08LX (%Lu) Wrong received datatype: %s - Skipping",
		    iBlock, evtData.fBlockCnt,
		    evtData.fEventID, evtData.fEventID, 
		    DataType2Text(inputDataType).c_str());
	  continue;
	}
      else 
	HLTDebug("We get the right data type: Block # %i/%i; Event 0x%08LX (%Lu) Received datatype: %s",
		    iBlock, evtData.fBlockCnt-1,
		    evtData.fEventID, evtData.fEventID, 
		    DataType2Text(inputDataType).c_str());

      TClonesArray* tracksArray = NULL;
      int ibForce = 0;
      if (fUseHLTTracks)
	{
	  tracksArray = new TClonesArray("AliTRDtrackV1");
	  HLTDebug("BLOCK fPtr 0x%x, fOffset %i, fSize %i, fSpec 0x%x, fDataType %s", block.fPtr, block.fOffset, block.fSize, block.fSpecification, DataType2Text(block.fDataType).c_str());
	  ReadTracks(tracksArray, block.fPtr, block.fSize);
	}
      else
	{
	  TObject *objIn = NULL;
	  if (iBlock == 0)
	    objIn = (TObject *)GetFirstInputObject( AliHLTTRDDefinitions::fgkTRDSATracksDataType, "TClonesArray", ibForce);
	  else
	    objIn = (TObject *)GetNextInputObject( ibForce );
	  HLTDebug("1stBLOCK, Pointer = 0x%x", objIn);
	  if (objIn){
	    tracksArray = (TClonesArray* )objIn;
	  }
	}

      if (tracksArray)
	{
	  Int_t nbEntries = tracksArray->GetEntries();
	  HLTDebug(" %i TRDtracks in tracksArray", nbEntries);
	  AliTRDtrackV1* trdTrack = 0x0;
	  for (Int_t i = 0; i < nbEntries; i++){
	    HLTDebug("%i/%i: ", i+1, nbEntries);
	    trdTrack = (AliTRDtrackV1*)tracksArray->At(i);
	    trdTrack->Print();
	    fTRDCalibraFillHisto->UpdateHistogramsV1(trdTrack);
	  }
	}
      

      TObjArray *outArray = FormOutput();
      if (outArray){
	PushBack(outArray, AliHLTTRDDefinitions::fgkCalibrationDataType);
	delete outArray;
      }
      
    }
  return 0;
  
}

/**
 * Read track to the TClonesArray from the memory 
 */
//============================================================================
Int_t AliHLTTRDCalibrationComponent::ReadTracks(TClonesArray *outArray, void* inputPtr, AliHLTUInt32_t size)
{
  AliHLTUInt8_t* iterPtr = (AliHLTUInt8_t* )inputPtr;
  
  cout << "\nReading tracks from the Memory\n ============= \n";
  AliHLTTRDTrack * hltTrack;
  AliHLTUInt32_t trackSize = 0, curSize = 0;
  Int_t counter=0;
  
  while (curSize < size)
    {
      hltTrack = (AliHLTTRDTrack*) iterPtr;
      HLTDebug("curSize %i, size %i",curSize, size);
      
      trackSize = hltTrack->GetSize();
      HLTDebug("GetSize() %i", trackSize);

      hltTrack->ReadTrackletsFromMemory(iterPtr + sizeof(AliHLTTRDTrack));

      AliTRDtrackV1* curTRDTrack = new((*outArray)[counter]) AliTRDtrackV1();
      hltTrack->ExportTRDTrack(curTRDTrack);
      
      curSize += trackSize; 
      iterPtr += trackSize;
      counter++;
    }

  //CheckTrackArray(outArray);
  return counter;
}


Int_t AliHLTTRDCalibrationComponent::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
//   //Int_t PushToFXS(TObject* pObject, const char* pDetector, const char* pFileID, const char* pDDLNumber = "");
//   //ireturn = PushToFXS(object, "TRD ", "TRDCalib", "1024 ");

  TObjArray *outArray = FormOutput();
  HLTDebug("Shipping data. Dummy. outputArray: pointer = 0x%x; NEntries = %i;", outArray, outArray->GetEntries());
  Int_t ireturn = PushToFXS(outArray, "TRD ", "TRDCalib");
  
  if (outArray){
    outArray->Delete();
    delete outArray;
  }
  
  return ireturn;
}


/**
 * Form output array of histrograms 
 */
//============================================================================
TObjArray* AliHLTTRDCalibrationComponent::FormOutput()
{
  TObjArray *outArray=new TObjArray(3);

  // gain histo
  TH2I *hCH2d = fTRDCalibraFillHisto->GetCH2d();
  outArray->Add(hCH2d);

  // drift velocity histo
  TProfile2D *hPH2d = fTRDCalibraFillHisto->GetPH2d();
  outArray->Add(hPH2d);
       
  // PRF histo
  TProfile2D *hPRF2d = fTRDCalibraFillHisto->GetPRF2d(); 
  outArray->Add(hPRF2d);

  HLTDebug("GetCH2d = 0x%x; NEntries = %i; size = %i", hCH2d, hCH2d->GetEntries(), sizeof(hCH2d));
  hCH2d->Print();
  HLTDebug("GetPH2d = 0x%x; NEntries = %i; size = %i", hPH2d, hPH2d->GetEntries(), sizeof(hPH2d));
  hPH2d->Print();
  HLTDebug("GetPRF2d = 0x%x; NEntries = %i; size = %i", hPRF2d, hPRF2d->GetEntries(), sizeof(hPRF2d));
  hPRF2d->Print();
  HLTDebug("output Array: pointer = 0x%x; NEntries = %i; size = %i", outArray, outArray->GetEntries(), sizeof(outArray));
    
  
  
  return outArray;
  
}

