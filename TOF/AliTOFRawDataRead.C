void AliTOFRawDataRead(Int_t iEvent=0)
{
  //
  // To read TOF raw data
  //

  Int_t indexDDL = 0;

  AliRawReaderFile reader(ievt);

  for (indexDDL = 0; indexDDL < 72; indexDDL++) {

    reader.Reset();
    AliTOFRawStream stream(&reader);
    reader.Select("TOF", indexDDL, indexDDL);

    //FILE *fpw = fopen("TOFrawDataRead.txt","w");

    while(stream.Next()) {

      if (stream.GetSector()==0) {

	printf("%2i %2i %2i %2i         %2i %1i %2i %1i %2i    %7i %8i\n",
	       stream.GetDDL(),stream.GetTRM(),
	       stream.GetTDC(),stream.GetChannel(),
	       stream.GetSector(),stream.GetPlate(),
	       stream.GetStrip(),stream.GetPadZ(),stream.GetPadX(),
	       stream.GetADCbin(),stream.GetTofBin());
	/*
	  Int_t iSector = stream.GetSector();
	  Int_t iPlate  = stream.GetPlate();
	  Int_t iStrip  = stream.GetStrip();
	  Int_t iPadZ   = stream.GetPadZ();
	  Int_t iPadX   = stream.GetPadX();
	  Int_t iTof    = stream.GetTofBin();
	  Int_t iAdc    = stream.GetADCbin();
	  
	  fprintf(fpw,"%2i %1i %2i %1i %2i %8i %7i\n",
	  iSector, iPlate, iStrip, iPadZ, iPadX, iTof, iAdc);
	*/

      }

    }

  }

}
