// gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/RAW -I$ALICE_ROOT/TOF")
// .L AliTOFRawDataRead.C++
// AliTOFRawDataRead()

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

  TClonesArray *fTOFrawData = new TClonesArray("AliTOFrawData",1000);
  Int_t fPackedDigits=0;

  Int_t indexDDL = 0;
  Int_t detectorIndex[5] = {-1, -1, -1, -1, -1};
  Int_t dummy = -1;

  AliRawReaderFile reader(iEvent);
  reader.RewindEvents();

  ofstream ftxt;
  ftxt.open("TOFrawDataReading.txt",ios::app);

  while (reader.NextEvent()) {

    AliTOFRawMap *rawMap = new AliTOFRawMap(fTOFrawData);

    Int_t slot[4] = {-1, -1, -1, -1};

   for (indexDDL = 0; indexDDL < AliDAQ::NumberOfDdls("TOF"); indexDDL++) {

     rawMap->Clear();
     fTOFrawData->Clear();
     fPackedDigits = 0;

     printf(" \n \n \n DRM number %2i \n \n \n ", indexDDL);

     reader.Reset();
     AliTOFRawStream stream(&reader);

     reader.Select("TOF", indexDDL, indexDDL);
     Bool_t signal = kFALSE;

     while(stream.Next()) {

       signal = (stream.GetSector()!=-1 &&
		 stream.GetPlate()!=-1 &&
		 stream.GetStrip()!=-1 &&
		 stream.GetPadZ()!=-1 &&
		 stream.GetPadX()!=-1);

       if (signal) {
	 printf("  %2i  %1i  %2i  %1i  %2i  \n", stream.GetSector(), stream.GetPlate(), stream.GetStrip(), stream.GetPadZ(), stream.GetPadX());

	 slot[0] = stream.GetTRM();
	 slot[1] = stream.GetTRMchain();
	 slot[2] = stream.GetTDC();
	 slot[3] = stream.GetTDCchannel();

	 if (rawMap->TestHit(slot) != kEmpty) {

	   AliTOFrawData *rawDigit = static_cast<AliTOFrawData*>(rawMap->GetHit(slot));
	   rawDigit->Update(stream.GetTofBin(), stream.GetToTbin(), stream.GetLeadingEdge(), stream.GetTrailingEdge(), stream.GetPSbit(), stream.GetACQ(), stream.GetErrorFlag());

	 }
	 else {

	   TClonesArray &arrayTofRawData =  *fTOFrawData;
	   new (arrayTofRawData[fPackedDigits++]) AliTOFrawData(stream.GetTRM(), stream.GetTRMchain(), stream.GetTDC(), stream.GetTDCchannel(),
								stream.GetTofBin(), stream.GetToTbin(), stream.GetLeadingEdge(), stream.GetTrailingEdge(),
								stream.GetPSbit(), stream.GetACQ(), stream.GetErrorFlag());

	   rawMap->SetHit(slot);

	   printf("  %6i \n", rawMap->GetHitIndex(slot));

	 }

       } // end if (signal)

     } // closed loop on TOF raw data per current DDL file

     printf("\n \n \n end of reading DRM number %2i\n", indexDDL);
     printf("                                     packed data %5i\n \n \n ", fTOFrawData->GetEntriesFast());

     for (Int_t iRawData = 0; iRawData<fTOFrawData->GetEntriesFast(); iRawData++) {

       AliTOFrawData *tofRawDatum = (AliTOFrawData*)fTOFrawData->UncheckedAt(iRawData);

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

       if (tofRawDatum->GetTOT()<10)                                    ftxt << "        " << tofRawDatum->GetTOT();
       else if (tofRawDatum->GetTOT()>=10 && tofRawDatum->GetTOT()<100) ftxt << "       " << tofRawDatum->GetTOT();
       else                                                             ftxt << "      " << tofRawDatum->GetTOT();
       if (tofRawDatum->GetTOF()<10)                                      ftxt << "      " << tofRawDatum->GetTOF() << endl;
       else if (tofRawDatum->GetTOF()>=10 && tofRawDatum->GetTOF()<100)   ftxt << "     " << tofRawDatum->GetTOF() << endl;
       else if (tofRawDatum->GetTOF()>=100 && tofRawDatum->GetTOF()<1000) ftxt << "    " << tofRawDatum->GetTOF() << endl;
       else                                                               ftxt << "   " << tofRawDatum->GetTOF() << endl;

     } // end loop

   } // endl loop on DDL files

   iEvent++;

  } // end while loop on event

  ftxt.close();

}
