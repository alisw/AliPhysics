#if !defined(__CINT__) || defined(__MAKECINT__)

// Root include files
#include "TClonesArray.h"

// AliRoot include files
#include "AliDAQ.h"
#include "AliHitMap.h"

#include "AliRawReader.h"
#include "AliRawReaderFile.h"

#include "AliTOFrawData.h"
#include "AliTOFRawMap.h"
#include "AliTOFRawStream.h"

#endif

extern gBenchmark;

void AliTOFRawDataReadBuffer(Int_t iEvent=0);

void AliTOFRawDataReadBuffer(Int_t iEvent)
{

  //
  // To read TOF raw data
  //

  Int_t Volume[5];
  AliTOFHitData *HitData;
  AliTOFHitDataBuffer *DataBuffer;
  AliTOFHitDataBuffer *PackedDataBuffer;

  /* create a tree for decoded data */
  TTree DataTree("DataTree", "Decoded Data");
  DataTree.Branch("HitData", "AliTOFHitData", &HitData);

  /* create a tree for decoded packed data */
  TTree PackedDataTree("PackedDataTree", "Decoded Packed Data");
  PackedDataTree.Branch("HitData", "AliTOFHitData", &HitData);
  
  AliRawReaderFile reader(iEvent);
  AliTOFRawStream stream(&reader);

  reader.RewindEvents();

  gBenchmark->Reset();
  /* loop over events */
  for (Int_t iEvent = 0; reader.NextEvent(); iEvent++) {
    printf("processing event %d\n", iEvent);

    /* reset buffers (actually not needed)*/
    stream.ResetBuffers();
    
    /* decode all DDLs */
    gBenchmark->Start("time");
    stream.DecodeDDL(0, AliDAQ::NumberOfDdls("TOF") - 1,0);
    gBenchmark->Stop("time");
    
    /* loop over DDLs */
    for (Int_t iDDL = 0; iDDL < AliDAQ::NumberOfDdls("TOF"); iDDL++){
      
      /* read decoded data */
      DataBuffer = stream.GetDataBuffer(iDDL);
      PackedDataBuffer = stream.GetPackedDataBuffer(iDDL);
      
      /* get buffer entries */
      Int_t nDBEntries = DataBuffer->GetEntries();
      Int_t nPDBEntries = PackedDataBuffer->GetEntries();

      /* read data buffer hits */
      for (Int_t iHit = 0; iHit < nDBEntries; iHit++){
      	HitData = DataBuffer->GetHit(iHit);
	HitData->SetDDLID(iDDL);
	/* add volume information to hit */
	stream.EquipmentId2VolumeId(HitData, HitData->GetVolume());
	DataTree.Fill();
      }
      /* reset buffer */
      DataBuffer->Reset();
      
      /* read data buffer hits */
      for (Int_t iHit = 0; iHit < nPDBEntries; iHit++){
      	HitData = PackedDataBuffer->GetHit(iHit);
	HitData->SetDDLID(iDDL);
	/* add volume information to hit */
	stream.EquipmentId2VolumeId(HitData, HitData->GetVolume());
	PackedDataTree.Fill();
      }
      /* reset buffer */
      PackedDataBuffer->Reset();
    }
    
  }
  gBenchmark->Print("time");

  TFile fileOut("TOF_rawQA.root", "RECREATE");
  DataTree.Write();
  PackedDataTree.Write();
  fileOut.Close();
  

}
