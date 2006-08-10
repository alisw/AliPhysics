void AliTOFRawDataRead(Int_t iEvent=0)
{
  //
  // To read TOF raw data
  //

  Int_t ii = 0;
  Int_t indexDDL = 0;
  Int_t detectorIndex[5] = {-1, -1, -1, -1, -1};

  AliRawReaderFile reader(iEvent);
  reader.RewindEvents();

  while (reader.NextEvent()) {

    for (indexDDL = 0; indexDDL < AliDAQ::NumberOfDdls("TOF"); indexDDL++) {

      reader.Reset();
      AliTOFRawStream stream(&reader);
      reader.Select("TOF", indexDDL, indexDDL);

      //FILE *fpw = fopen("TOFrawDataRead.txt","w");

      while(stream.Next()) {

	for (ii=0; ii<5; ii++) detectorIndex[ii] = -1;

	detectorIndex[0] = (Int_t)stream.GetSector();
	detectorIndex[1] = (Int_t)stream.GetPlate();
	detectorIndex[2] = (Int_t)stream.GetStrip();
	detectorIndex[3] = (Int_t)stream.GetPadZ();
	detectorIndex[4] = (Int_t)stream.GetPadX();

	if (detectorIndex[0]==-1 ||
	    detectorIndex[1]==-1 ||
	    detectorIndex[2]==-1 ||
	    detectorIndex[3]==-1 ||
	    detectorIndex[4]==-1) continue;
	else {

	  printf("%2i  %2i  %2i  %2i  %2i      %2i  %1i  %2i  %1i  %2i    %7i  %8i\n",
		 stream.GetDDL(),stream.GetTRM(),stream.GetTRMchain(),
		 stream.GetTDC(),stream.GetTDCchannel(),
		 stream.GetSector(),stream.GetPlate(),
		 stream.GetStrip(),stream.GetPadZ(),stream.GetPadX(),
		 stream.GetToTbin(),stream.GetTofBin());

	} // end else

      } // end while loop on next stream

    } // endl loop on DDL files

    iEvent++;

  } // end while loop on event

}
