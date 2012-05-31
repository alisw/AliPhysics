//_____________________________________________________//
//                                                     //
//    This macro reads VZERO DDL Raw Data and          //
//    converts it into Digits                          //
//                                                     //
//____________________________________________________ //


void VZERORaw2Digits(Int_t nEvent = 1) 
{
     for(int i=0;i<nEvent;i++) {
		printf("=========== EVENT  %d ===========\n",i);
		TString DirName = "raw";
		DirName += i;
		RawStreamEvent(DirName.Data()); }
}


Bool_t RawStreamEvent(TString fileName = "./")
{
  // Reads DDL data from fileName

  TStopwatch timer;
  timer.Start();

// Creates a TreeD to dump Digits

  AliRunLoader* rl = AliRunLoader::Open("galice.root");    

  AliVZEROLoader* loader = (AliVZEROLoader*) rl->GetLoader("VZEROLoader");
  
  if(!loader) {
    AliError("no VZERO loader found");
    return kFALSE; }

  TTree* treeD  = loader->TreeD();
  if(!treeD) {
      loader->MakeTree("D");
      treeD = loader->TreeD(); }
        
  AliVZEROdigit  digit;
  AliVZEROdigit* pdigit = &digit;
  const Int_t kBufferSize = 4000;
   
  treeD->Branch("VZERO", "AliVZERODigit",  &pdigit, kBufferSize);

  AliRawReader* rawReader = 0x0;
  rawReader = new AliRawReaderFile(fileName); // DDL files
  
  AliVZERORawStream* rawStream  = new AliVZERORawStream(rawReader);    
     
  rawReader->NextEvent();
    
  rawStream->Next();
  
  for(Int_t i=0; i<64; i++) {
      new(pdigit) AliVZEROdigit(i, (Int_t)rawStream->GetADC(i), (Int_t)rawStream->GetTime(i)); 
      treeD->Fill();
  }
 
// Checks if everything is OK by printing results 

//   for(int i=0;i<64;i++) {
// 	printf("Channel %d : %d %d \n",i,rawStream->GetADC(i),rawStream->GetTime(i)); }
//   treeD->Print(); printf(" \n"); 
   	
  loader->WriteDigits("OVERWRITE");
  loader->UnloadDigits();	
	
  delete rawReader;
  delete rawStream;

  timer.Stop();
  timer.Print();
}


