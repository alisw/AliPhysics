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
#include "AliHLTTRDUtils.h"

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliRawReaderMemory.h"
#include "AliTRDCalibraFillHisto.h"
#include "AliTRDtrackV1.h"

#include <cstdlib>
#include <cerrno>
#include <string>

ClassImp(AliHLTTRDCalibrationComponent);

AliHLTTRDCalibrationComponent::AliHLTTRDCalibrationComponent()
: AliHLTCalibrationProcessor(),
  fTRDCalibraFillHisto(NULL),
  fOutputSize(50000),
  fTracksArray(NULL),
  fOutArray(NULL),
  fNevent(0),
  feveryNevent(20),
  fRecievedTimeBins(kFALSE)
{
  // Default constructor
}

AliHLTTRDCalibrationComponent::~AliHLTTRDCalibrationComponent()
{
  // Destructor
}

const char* AliHLTTRDCalibrationComponent::GetComponentID()
{
  // Return the component ID const char *
  return "TRDCalibration"; // The ID of this component
}

void AliHLTTRDCalibrationComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // Get the list of input data
  list.clear(); // We do not have any requirements for our input data type(s).
  list.push_back(AliHLTTRDDefinitions::fgkTRDSATracksDataType);
}

AliHLTComponentDataType AliHLTTRDCalibrationComponent::GetOutputDataType()
{
  // Get the output data type
  return AliHLTTRDDefinitions::fgkCalibrationDataType;
}

void AliHLTTRDCalibrationComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // Get the output data size
  constBase = fOutputSize;
  inputMultiplier = 0;
}

AliHLTComponent* AliHLTTRDCalibrationComponent::Spawn()
{
  // Spawn function, return new instance of this class
  return new AliHLTTRDCalibrationComponent;
};

Int_t AliHLTTRDCalibrationComponent::ScanArgument( int argc, const char** argv )
{
  // perform initialization. We check whether our relative output size is specified in the arguments.
  int i = 0;
  char* cpErr;
  while ( i < argc )
    {
      HLTDebug("argv[%d] == %s", i, argv[i] );
      if ( !strcmp( argv[i], "output_size" ) )
        {
          if ( i+1>=argc )
            {
              HLTError("Missing output_size parameter");
              return ENOTSUP;
            }
          HLTDebug("argv[%d+1] == %s", i, argv[i+1] );
          fOutputSize = strtoul( argv[i+1], &cpErr, 0 );
          if ( *cpErr )
            {
              HLTError("Cannot convert output_size parameter '%s'", argv[i+1] );
              return EINVAL;
            }
          HLTInfo("Output size set to %lu %%", fOutputSize );
          i += 2;
          continue;
        }
      if ( !strcmp( argv[i], "-everyNevent" ) )
        {
          if ( i+1>=argc )
            {
              HLTError("Missing everyNevent parameter");
              return ENOTSUP;
            }
          HLTDebug("argv[%d+1] == %s", i, argv[i+1] );
          fOutputSize = strtoul( argv[i+1], &cpErr, 0 );
          if ( *cpErr )
            {
              HLTError("Cannot convert everyNevent parameter '%s'", argv[i+1] );
              return EINVAL;
            }
          HLTInfo("Pushing back every %d event", feveryNevent);
          i += 2;
          continue;
        }

      else {
        HLTError("Unknown option '%s'", argv[i] );
        return EINVAL;
      }
    }
  return 0;
}

Int_t AliHLTTRDCalibrationComponent::InitCalibration()
{
  if(!AliCDBManager::Instance()->IsDefaultStorageSet()){
    HLTError("DefaultStorage is not set in CDBManager");
    return -EINVAL;
  }
  if(AliCDBManager::Instance()->GetRun()<0){
    HLTError("Run Number is not set in CDBManager");
    return -EINVAL;
  }
  HLTInfo("CDB default storage: %s; RunNo: %i", (AliCDBManager::Instance()->GetDefaultStorage()->GetBaseFolder()).Data(), AliCDBManager::Instance()->GetRun());

  fTRDCalibraFillHisto = AliTRDCalibraFillHisto::Instance();
  fTRDCalibraFillHisto->SetHisto2d(); // choose to use histograms
  fTRDCalibraFillHisto->SetCH2dOn();  // choose to calibrate the gain
  fTRDCalibraFillHisto->SetPH2dOn();  // choose to calibrate the drift velocity
  fTRDCalibraFillHisto->SetPRF2dOn(); // choose to look at the PRF
  fTRDCalibraFillHisto->SetIsHLT(); // per detector
  //fTRDCalibraFillHisto->SetDebugLevel(1);// debug
  fTRDCalibraFillHisto->SetFillWithZero(kTRUE);

  fTracksArray = new TClonesArray("AliTRDtrackV1");
  fOutArray = new TObjArray(3);

  return 0;
}

Int_t AliHLTTRDCalibrationComponent::DeinitCalibration()
{
  
  // Deinitialization of the component
  
  HLTDebug("DeinitCalibration");
  fTracksArray->Delete();
  delete fTracksArray;
  fTRDCalibraFillHisto->Destroy();
  //fOutArray->Delete();
  delete fOutArray;

  return 0;
}

Int_t AliHLTTRDCalibrationComponent::ProcessCalibration(const AliHLTComponent_EventData& evtData,
                                                        const AliHLTComponent_BlockData* blocks,
                                                        AliHLTComponent_TriggerData& /*trigData*/,
                                                        AliHLTUInt8_t* /*outputPtr*/,
                                                        AliHLTUInt32_t& /*size*/,
                                                        vector<AliHLTComponent_BlockData>& /*outputBlocks*/)
{
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
      if (!correctDataType) {
        HLTDebug( "Block # %i/%i; Event 0x%08LX (%Lu) Wrong received datatype: %s - Skipping",
                  iBlock, evtData.fBlockCnt,
                  evtData.fEventID, evtData.fEventID,
                  DataType2Text(inputDataType).c_str());
        continue;
      }
      else {
        HLTDebug("We get the right data type: Block # %i/%i; Event 0x%08LX (%Lu) Received datatype: %s; Block Size: %i",
                 iBlock, evtData.fBlockCnt-1,
                 evtData.fEventID, evtData.fEventID,
                 DataType2Text(inputDataType).c_str(),
		 block.fSize);
      }

      Int_t nTimeBins;
      AliHLTTRDUtils::ReadTracks(fTracksArray, block.fPtr, block.fSize, &nTimeBins);
      
      if(!fRecievedTimeBins){
	HLTDebug("Reading number of time bins from input block. Value is: %d", nTimeBins);
	fTRDCalibraFillHisto->Init2Dhistos(); // initialise the histos
	fTRDCalibraFillHisto->SetNumberClusters(0); // At least 1 clusters
	fTRDCalibraFillHisto->SetNumberClustersf(nTimeBins); // Not more than %d  clusters
	fRecievedTimeBins=kTRUE;
      }

      Int_t nbEntries = fTracksArray->GetEntries();
      HLTDebug(" %i TRDtracks in tracksArray", nbEntries);
      AliTRDtrackV1* trdTrack = 0x0;
      for (Int_t i = 0; i < nbEntries; i++){
	HLTDebug("%i/%i: ", i+1, nbEntries);
	trdTrack = (AliTRDtrackV1*)fTracksArray->At(i);
	trdTrack->Print();
	fTRDCalibraFillHisto->UpdateHistogramsV1(trdTrack);
      }
      
      if(!fOutArray->At(0))FormOutput();
      if (fNevent%feveryNevent==0 && fOutArray) {
        PushBack(fOutArray, AliHLTTRDDefinitions::fgkCalibrationDataType);
      }

      fTracksArray->Delete();
      fNevent++;

    }
  return 0;

}

/**
 * Form output array of histrograms
 */
//============================================================================
void AliHLTTRDCalibrationComponent::FormOutput()
{
  // gain histo
  TH2I *hCH2d = fTRDCalibraFillHisto->GetCH2d();
  fOutArray->Add(hCH2d);

  // drift velocity histo
  TProfile2D *hPH2d = fTRDCalibraFillHisto->GetPH2d();
  fOutArray->Add(hPH2d);

  // PRF histo
  TProfile2D *hPRF2d = fTRDCalibraFillHisto->GetPRF2d();
  fOutArray->Add(hPRF2d);

  HLTDebug("GetCH2d = 0x%x; NEntries = %i; size = %i", hCH2d, hCH2d->GetEntries(), sizeof(hCH2d));
  hCH2d->Print();
  HLTDebug("GetPH2d = 0x%x; NEntries = %i; size = %i", hPH2d, hPH2d->GetEntries(), sizeof(hPH2d));
  hPH2d->Print();
  HLTDebug("GetPRF2d = 0x%x; NEntries = %i; size = %i", hPRF2d, hPRF2d->GetEntries(), sizeof(hPRF2d));
  hPRF2d->Print();
  HLTDebug("output Array: pointer = 0x%x; NEntries = %i; size = %i", fOutArray, fOutArray->GetEntries(), sizeof(fOutArray));
}

Int_t AliHLTTRDCalibrationComponent::ShipDataToFXS(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  //fTRDCalibraFillHisto->DestroyDebugStreamer();
  PushToFXS((TObject*)fOutArray, "TRD", "GAINDRIFTPRF");
}
