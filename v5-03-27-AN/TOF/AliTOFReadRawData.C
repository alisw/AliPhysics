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

void AliTOFReadRawData(const char *input, const Int_t evIn, const Char_t *outFile, Int_t runNumber=137366,Int_t decoderVersion=2);

void AliTOFReadRawData(const char *input, const Int_t evIn, const Char_t *outFile, Int_t runNumber, Int_t decoderVersion)
{

  //
  // To read TOF raw data
  //

  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("raw://");
  cdb->SetRun(runNumber);

  AliGeomManager::LoadGeometry();
  AliGeomManager::ApplyAlignObjsFromCDB("ITS TPC TRD TOF HMPID PHOS EMCAL T0 VZERO FMD PMD MUON ZDC");

  AliTOFcalib *fTOFcalib = new AliTOFcalib();
  fTOFcalib->Init();

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

  AliRawReader *reader;  
  TString fileName(input);
  if (fileName.EndsWith("/")) {
    reader = new AliRawReaderFile(fileName);
  } else if (fileName.EndsWith(".root")) {
    reader = new AliRawReaderRoot(fileName);
  } else if (!fileName.IsNull()) {
    reader = new AliRawReaderDate(fileName);
    reader->SelectEvents(7);
  }

  TClonesArray *clonesRawData = new TClonesArray("AliTOFrawData",1000);
  AliTOFRawStream stream;

  gBenchmark->Reset();

  Int_t countNoise = 0;

  /* loop over events */
  for (Int_t iEvent = 0; reader->NextEvent(); iEvent++) {
    if (iEvent!=evIn) continue;
      printf("processing event %d\n", iEvent);

    stream.Clear();
    stream.SetRawReader(reader);

    //reader->RewindEvents();

    UInt_t timeStamp = reader->GetTimestamp();

    gBenchmark->Start("time");

    /* loop over DDLs */
    for (Int_t iDDL = 0; iDDL < AliDAQ::NumberOfDdls("TOF"); iDDL++){

      reader->Reset();

      /* reset buffers (actually not needed)*/
      stream.ResetBuffers();
    
      switch(decoderVersion) {
      case 1:
	stream.LoadRawDataBuffers(iDDL,0);
	break;
      case 2:
	stream.LoadRawDataBuffersV2(iDDL,0);
	break;
      case 0:
	stream.LoadRawData(iDDL);
	break;
      }

      clonesRawData = (TClonesArray*)stream.GetRawData();
      for (Int_t iRawData = 0; iRawData<clonesRawData->GetEntriesFast(); iRawData++) {

	AliTOFrawData *tofRawDatum = (AliTOFrawData*)clonesRawData->UncheckedAt(iRawData);

	if (tofRawDatum->GetTOF()==-1) continue;
	Int_t detectorIndex[5]={-1,-1,-1,-1,-1};
	stream.EquipmentId2VolumeId(iDDL, tofRawDatum->GetTRM(), tofRawDatum->GetTRMchain(),
				    tofRawDatum->GetTDC(), tofRawDatum->GetTDCchannel(), detectorIndex);
	Int_t dummy = detectorIndex[3];
	detectorIndex[3] = detectorIndex[4];
	detectorIndex[4] = dummy;

	/* check valid index */
	if (detectorIndex[0]==-1||detectorIndex[1]==-1||detectorIndex[2]==-1||detectorIndex[3]==-1||detectorIndex[4]==-1) continue;

	// Do not reconstruct anything in the holes
	if (detectorIndex[0]==13 || detectorIndex[0]==14 || detectorIndex[0]==15 ) { // sectors with holes
	  if (detectorIndex[1]==2) { // plate with holes
	    inholes++;
	    continue;
	  }
	}


	Int_t index = AliTOFGeometry::GetIndex(detectorIndex);

	/* check channel enabled */
	//if (!fTOFcalib->IsChannelEnabled(index)) continue; // bad channels

	/* get cluster info */
	Float_t tot = tofRawDatum->GetTOT() * AliTOFGeometry::ToTBinWidth() * 1.e-3; /* ns */
	Float_t timeRaw = tofRawDatum->GetTOF() * AliTOFGeometry::TdcBinWidth(); /* ps */
	Float_t time = tofRawDatum->GetTOF() * AliTOFGeometry::TdcBinWidth(); /* ps */
	Float_t deltaBC = tofRawDatum->GetDeltaBC(); // deltaBC
	Float_t l0l1 = tofRawDatum->GetL0L1Latency(); // L0-L1 latency
	Float_t xHit, yHit, zHit;
	Int_t sector, strip, padx, padz;
	ofstream ftxt;

	/* get correction */
	Float_t corr = fTOFcalib->GetTimeCorrection(index, tot, deltaBC, l0l1, timeStamp); /* ps */

	/* apply time correction */
	time -= corr;
	
	//if (iEvent==12)
#if 0	
	cout << " " << index
	     << " " << fTOFcalib->IsChannelEnabled(index)
	     << " " << tot
	     << " " << timeRaw
	     << " " << time
	     << endl;
#endif
	
	/* geometry */
	AliTOFGeometry geo;
	for (Int_t j=0; j<5; j++) detectorIndex[j] = -1;
	geo.GetVolumeIndices(index, detectorIndex);
	xHit  = geo.GetX(detectorIndex);
	yHit  = geo.GetY(detectorIndex);
	zHit  = geo.GetZ(detectorIndex);
	
	//if (fTOFcalib->IsChannelEnabled(index) == 1) {
	  sector = detectorIndex[0];
	  strip  = detectorIndex[2];
	  padx   = detectorIndex[4];
	  padz   = detectorIndex[3];
	  
	  //printf("x %f, y %f, z %f\n",xHit,yHit,zHit);
	  //printf("sector %d, strip %d, padx %d, padz %d\n",sector,strip,padx,padz);
#if 1  
	  ftxt.open(outFile,ios::app); 
	  ftxt << time << " " << timeRaw << " " << tot << " " << sector 
	       << " " << strip << " " << padx << " " << padz 
	       << " " << xHit << " " << yHit << " " << zHit << endl;
	  ftxt.close();  
#endif
	  //}
	  //else countNoise++;

      } // closed loop on TOF raw data per current DDL file
     
      clonesRawData->Clear();




      /*
      // read decoded data
      DataBuffer = stream.GetDataBuffer(iDDL);
      PackedDataBuffer = stream.GetPackedDataBuffer(iDDL);
      
      // get buffer entries
      Int_t nDBEntries = DataBuffer->GetEntries();
      Int_t nPDBEntries = PackedDataBuffer->GetEntries();

      // read data buffer hits
      for (Int_t iHit = 0; iHit < nDBEntries; iHit++){
      	HitData = DataBuffer->GetHit(iHit);
	HitData->SetDDLID(iDDL);
	// add volume information to hit
	stream.EquipmentId2VolumeId(HitData, HitData->GetVolume());
	DataTree.Fill();
      }
      // reset buffer
      DataBuffer->Reset();
      
      // read data buffer hits
      for (Int_t iHit = 0; iHit < nPDBEntries; iHit++){
      	HitData = PackedDataBuffer->GetHit(iHit);
	HitData->SetDDLID(iDDL);
	// add volume information to hit
	stream.EquipmentId2VolumeId(HitData, HitData->GetVolume());
	PackedDataTree.Fill();
      }

      // reset buffer
      PackedDataBuffer->Reset();
      */

    } // loop on DDLs
    
  }
  printf("noise %d\n",countNoise);
  gBenchmark->Print("time");

  /*
  TFile fileOut("TOF_rawQA.root", "RECREATE");
  DataTree.Write();
  PackedDataTree.Write();
  fileOut.Close();
  */

}
