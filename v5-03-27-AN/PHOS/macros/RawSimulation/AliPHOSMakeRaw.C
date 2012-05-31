/* $Id$ */
MakeRaw()
{
  // Create DDL files manually by firing each crystal.
  // The DDL files are stored into directory raw0
  // Yuri Kharlov. 15 February 2006
  
  const Int_t kAdcThreshold = 1; // Lower ADC threshold to write to raw data

  Int_t prevDDL = -1;

  // Create a shaper pulse object
  AliPHOSPulseGenerator pulse;

  Int_t *adcValuesLow = new Int_t[pulse.GetRawFormatTimeBins()];
  Int_t *adcValuesHigh= new Int_t[pulse.GetRawFormatTimeBins()];

  const Int_t maxDDL = 20;
  AliAltroBuffer  *buffer[maxDDL];
  AliAltroMapping *mapping[maxDDL];

  // loop over digits
  Float_t energy = 2.;
  Float_t time   = 15.6e-09;
  Int_t module   = 1;
  for (Int_t iX=0; iX<64; iX++) {
    for (Int_t iZ=0; iZ<56; iZ++) {
      
      Int_t iRCU = -111;
      if     ( 0<=iX&&iX<32 &&  0<=iZ&&iZ<28) iRCU=0;
      else if( 0<=iX&&iX<32 && 28<=iZ&&iZ<56) iRCU=1;
      else if(32<=iX&&iX<64 &&  0<=iZ&&iZ<28) iRCU=2;
      else if(32<=iX&&iX<64 && 28<=iZ&&iZ<56) iRCU=3;
      Int_t iDDL = 4 * (module - 1) + iRCU;
      
      // new DDL
      if (iDDL != prevDDL) {
	if (buffer[iDDL] == 0) {
	  // open new file and write dummy header
	  TString fileName = AliDAQ::DdlFileName("PHOS",iDDL);
	  
	  TString path = gSystem->Getenv("ALICE_ROOT");
	  path += "/PHOS/mapping/RCU";
	  path += iRCU;
	  path += ".data";
	  
	  // 	printf("DDL %d, mapping %s\n",iDDL,path.Data());
	  mapping[iDDL] = new AliCaloAltroMapping(path.Data());
	  buffer[iDDL]  = new AliAltroBuffer(fileName.Data(),mapping[iDDL]);
	  buffer[iDDL]->WriteDataHeader(kTRUE, kFALSE);  //Dummy;
	}
	prevDDL = iDDL;
      }
      
      // if a signal is out of time range, write only trailer
      if (time > pulse.GetRawFormatTimeMax()*0.5 ) {
	AliInfo("Signal is out of time range.\n");
	buffer[iDDL]->FillBuffer(0);
	buffer[iDDL]->FillBuffer(pulse.GetRawFormatTimeBins() );  // time bin
	buffer[iDDL]->FillBuffer(3);                               // bunch length
	buffer[iDDL]->WriteTrailer(3, iZ, iX, 0);  // trailer
      }
      pulse.SetAmplitude(energy);
      pulse.SetTZero(time);
      pulse.MakeSamples();
      pulse.GetSamples(adcValuesHigh, adcValuesLow) ; 
      //       printf("digit E=%.4f GeV, t=%g s, (mod,iZ,iX)=(%d,%d,%d)\n",
      //  	     energy,time,module-1,iZ,iX);
      buffer[iDDL]->WriteChannel(iZ, iX, 0, 
				 pulse.GetRawFormatTimeBins(), 
				 adcValuesLow , kAdcThreshold);
      buffer[iDDL]->WriteChannel(iZ, iX, 1, 
				 pulse.GetRawFormatTimeBins(), 
				 adcValuesHigh, kAdcThreshold);
    }
  }
  // write real header and close last file
  for (Int_t iDDL=0; iDDL<maxDDL; iDDL++) {
    if (buffer[iDDL]) {
      buffer[iDDL]->Flush();
      buffer[iDDL]->WriteDataHeader(kFALSE, kFALSE);
      delete buffer[iDDL];
      if (mapping[iDDL]) delete mapping[iDDL];
    }
  }

  // move DDL files to directory raw0
  gSystem->Exec("mkdir -p raw0; mv *.ddl raw0");
}
//-----------------------------------------------------------------------------
PlotRawDDL()
{
  // Macro to read raw data from the DDL files 
  // and fill histogram with fired channels.
  // Yuri Kharlov. 19 October 2006

  AliRawReaderFile* rf = new AliRawReaderFile(".");
  AliCaloRawStream in(rf,"PHOS");
  in.SetOldRCUFormat(kFALSE); // Set to kTRUE for beam test 2006 only
  Int_t iBin = 0;
  Int_t iEvent=0;

  TFile *file = new TFile("PlotRawDDL.root","recreate");
  TH2F *hXY = new TH2F("hXY","X,Y of raw data",64,0.,64.,56,0.,56.);
  TH2F *hXYmodule[5];
  for (Int_t iMod=0; iMod<5; iMod++) {
    TString name="hXY";
    TString title="X,Y of raw data in module ";
    hXYmodule[iMod] = new TH2F(name+iMod,title+iMod,64,0.,64.,56,0.,56.);
  }

  // Loop over events
  while (rf->NextEvent()) {
    Int_t runNum = rf->GetRunNumber();
    cout << "Run "<<runNum<<", event "<<iEvent<<" of type " 
	 << rf->GetType() << endl;

    // Loop over channels within this event
    
    while ( in.Next() ) { 
      if (in.IsNewHWAddress()) {
	  hXY->Fill((Double_t)in.GetRow(),(Double_t)in.GetColumn());
	  hXYmodule[in.GetModule()]->Fill((Double_t)in.GetRow(),(Double_t)in.GetColumn());
      }
    }
    iEvent++;
  }
  file->Write();
  file->Close();
}
