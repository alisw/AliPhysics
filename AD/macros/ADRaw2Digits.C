//_____________________________________________________//
//                                                     //
//    This macro reads AD DDL Raw Data and          //
//    converts it into Digits                          //
//                                                     //
//____________________________________________________ //


void ADRaw2Digits(Int_t nEvent = 5) 
{
     for(int iEvent=0; iEvent<nEvent; iEvent++) {
		printf("=========== EVENT  %d ===========\n",iEvent);
		TString FileName = "raw.root";
		RawStreamEvent(FileName.Data(),iEvent); }
}


Bool_t RawStreamEvent(TString fileName = "./",Int_t iEvent = 0)
{
  // Reads DDL data from fileName

  TStopwatch timer;
  timer.Start();

// Creates a TreeD to dump Digits

  AliRunLoader* rl = AliRunLoader::Open("galice.root");    

  AliADLoader* loader = (AliADLoader*) rl->GetLoader("ADLoader");
  
  if(!loader) {
    AliError("no AD loader found");
    return kFALSE; }

  TTree* treeD  = loader->TreeD();
  if(!treeD) {
      loader->MakeTree("D");
      treeD = loader->TreeD(); }
        
  AliADdigit  digit;
  AliADdigit* pdigit = &digit;
  const Int_t kBufferSize = 4000;
   
  treeD->Branch("AD", "AliADdigit",  &pdigit, kBufferSize);

  AliRawReader* rawReader = 0x0;
  rawReader = new AliRawReaderRoot(fileName); // DDL files
  
  AliADRawStream* rawStream  = new AliADRawStream(rawReader);    
  
  rawReader->GotoEvent(iEvent);  
     
  rawStream->Next();
  
  printf("Data size is %d\n",rawReader->GetDataSize());
  
  for(Int_t i=0; i<16; i++) {
  
      Short_t chargeADC[21];
      for(Int_t iClock=0; iClock < 21; ++iClock) {
	chargeADC[iClock] = rawStream->GetPedestal(i,iClock);
      }

      new(pdigit) AliADdigit(i, (Int_t)rawStream->GetTime(i), (Int_t)rawStream->GetWidth(i), rawStream->GetIntegratorFlag(i, 10),chargeADC); 
      treeD->Fill();
  }
 
// Checks if everything is OK by printing results 

   for(int i=0;i<16;i++) {
 	printf("Channel %d : %d %d \n",i,rawStream->GetADC(i),rawStream->GetTime(i)); }
   //treeD->Print(); printf(" \n"); 
   	
  loader->WriteDigits("OVERWRITE");
  loader->UnloadDigits();
  
  	
	
  delete rawReader;
  delete rawStream;

  timer.Stop();
  timer.Print();
}


