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

void AliTOFRawDataRead(Int_t iEvent=0);

void AliTOFRawDataRead(Int_t iEvent)
{
  //
  // To read TOF raw data
  //

  AliTOFrawData *tofRawDatum=new AliTOFrawData();
  TTree *PackedDataTree= new TTree("PackedDataTree", "Decoded Packed Data");
  PackedDataTree->Branch("HitData", "AliTOFrawData", &tofRawDatum);

  TClonesArray *clonesRawData = new TClonesArray("AliTOFrawData",1000);
  Int_t fPackedDigits=0;

  Int_t detectorIndex[5] = {-1, -1, -1, -1, -1};
  Int_t dummy = -1;

  AliRawReader *reader = new AliRawReaderFile(iEvent);
  reader->RewindEvents();

  ofstream ftxt;
  ftxt.open("TOFrawDataReading.txt",ios::app);

  while (reader->NextEvent()) {

    AliTOFRawMap *rawMap = new AliTOFRawMap(clonesRawData);

    Int_t slot[4] = {-1, -1, -1, -1};

   for (Int_t indexDDL = 0; indexDDL < AliDAQ::NumberOfDdls("TOF"); indexDDL++) {

     reader->Reset();
     AliTOFRawStream stream(reader);
     stream.LoadRawData(indexDDL);

     clonesRawData = (TClonesArray*)stream.GetRawData();

     for (Int_t iRawData = 0; iRawData<clonesRawData->GetEntriesFast(); iRawData++) {

       tofRawDatum = (AliTOFrawData*)clonesRawData->UncheckedAt(iRawData);

       if (tofRawDatum->GetTOT()==-1 || tofRawDatum->GetTOF()==-1) continue;

       if (indexDDL<10) ftxt << "  " << indexDDL;
       else             ftxt << " " << indexDDL;
       if (tofRawDatum->GetTRM()<10) ftxt << "  " << tofRawDatum->GetTRM();
       else                          ftxt << " " << tofRawDatum->GetTRM();
       ftxt << "  " << tofRawDatum->GetTRMchain();
       if (tofRawDatum->GetTDC()<10) ftxt << "  " << tofRawDatum->GetTDC();
       else                          ftxt << " " << tofRawDatum->GetTDC();
       ftxt << "  " << tofRawDatum->GetTDCchannel();

       stream.EquipmentId2VolumeId(indexDDL, tofRawDatum->GetTRM(), tofRawDatum->GetTRMchain(),
				   tofRawDatum->GetTDC(), tofRawDatum->GetTDCchannel(), detectorIndex);
       dummy = detectorIndex[3];
       detectorIndex[3] = detectorIndex[4];
       detectorIndex[4] = dummy;

       if (detectorIndex[0]<10) ftxt  << "  ->  " << detectorIndex[0];
       else                     ftxt  << "  -> " << detectorIndex[0];
       ftxt << "  " << detectorIndex[1];
       if (detectorIndex[2]<10) ftxt << "  " << detectorIndex[2];
       else                     ftxt << " " << detectorIndex[2];
       ftxt << "  " << detectorIndex[3];
       if (detectorIndex[4]<10) ftxt << "  " << detectorIndex[4];
       else                     ftxt << " " << detectorIndex[4];

       if (tofRawDatum->GetTOT()<10)                                            ftxt << "        " << tofRawDatum->GetTOT();
       else if (tofRawDatum->GetTOT()>=10 && tofRawDatum->GetTOT()<100)         ftxt << "       " << tofRawDatum->GetTOT();
       else if (tofRawDatum->GetTOT()>=100 && tofRawDatum->GetTOT()<1000)       ftxt << "      " << tofRawDatum->GetTOT();
       else if (tofRawDatum->GetTOT()>=1000 && tofRawDatum->GetTOT()<10000)     ftxt << "     " << tofRawDatum->GetTOT();
       else if (tofRawDatum->GetTOT()>=10000 && tofRawDatum->GetTOT()<100000)   ftxt << "    " << tofRawDatum->GetTOT();
       else if (tofRawDatum->GetTOT()>=100000 && tofRawDatum->GetTOT()<1000000) ftxt << "   " << tofRawDatum->GetTOT();
       else                                                                     ftxt << "  " << tofRawDatum->GetTOT();
       if (tofRawDatum->GetTOF()<10)                                            ftxt << "        " << tofRawDatum->GetTOF() << endl;
       else if (tofRawDatum->GetTOF()>=10 && tofRawDatum->GetTOF()<100)         ftxt << "       " << tofRawDatum->GetTOF() << endl;
       else if (tofRawDatum->GetTOF()>=100 && tofRawDatum->GetTOF()<1000)       ftxt << "      " << tofRawDatum->GetTOF() << endl;
       else if (tofRawDatum->GetTOF()>=1000 && tofRawDatum->GetTOF()<10000)     ftxt << "     " << tofRawDatum->GetTOF() << endl;
       else if (tofRawDatum->GetTOF()>=10000 && tofRawDatum->GetTOF()<100000)   ftxt << "    " << tofRawDatum->GetTOF() << endl;
       else if (tofRawDatum->GetTOF()>=100000 && tofRawDatum->GetTOF()<1000000) ftxt << "   " << tofRawDatum->GetTOF() << endl;
       else                                                                     ftxt << "  " << tofRawDatum->GetTOF() << endl;

       PackedDataTree->Fill();
     } // end loop

   } // endl loop on DDL files

   iEvent++;

  } // end while loop on event

  ftxt.close();
  TFile fileOut("TOF_rawQA_OldDecoder.root", "RECREATE");
  PackedDataTree->Write();
  fileOut.Close();

}

