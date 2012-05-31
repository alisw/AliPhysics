//_____________________________________________________//
//                                                     //
//    This macro reads ACORDE DDL Raw Data and          //
//    converts it into Digits                          //
//                                                     //
//____________________________________________________ //


void ACORDERaw2Digits(Int_t nEvents = 1, char* fileName = "rawdata.root")
{
  // Reads DDL data from fileName

  TStopwatch timer;
  timer.Start();

// Creates a TreeD to dump Digits

  AliRunLoader* rl = AliRunLoader::Open("galice.root");    

  AliACORDELoader* loader = (AliACORDELoader*) rl->GetLoader("ACORDELoader");
  
  if(!loader) {
    AliError("no ACORDE loader found");
    return kFALSE; }

  TTree* treeD  = loader->TreeD();
  if(!treeD) {
      loader->MakeTree("D");
      treeD = loader->TreeD(); }
        
  AliACORDEdigit  digit;
  AliACORDEdigit* pdigit = &digit;
  const Int_t kBufferSize = 4000;
   
  treeD->Branch("ACORDE", "AliACORDEdigit",  &pdigit, kBufferSize);

  AliRawReader* rawReader = 0x0;
//  rawReader = new AliRawReaderFile(fileName); // DDL files
  rawReader = new AliRawReaderRoot(fileName); // DDL files

  AliACORDERawStream* rawStream  = new AliACORDERawStream(rawReader);    

  for (Int_t i=0; i<nEvents; i++) {
    printf("=========== EVENT  %d ===========\n",i);
    if (!rawReader->NextEvent())
      break;

    rawStream->Reset();
    if (!rawStream->Next())
      break;
    printf("Data size is %d\n",rawStream->DataSize());
  /*  
  for(Int_t i=0; i<64; i++) {
      new(pdigit) AliACORDEdigit(i, (Int_t)rawStream->GetADC(i), (Int_t)rawStream->GetTime(i)); 
      treeD->Fill();
  }
  */ 
// Checks if everything is OK by printing results 

//   for(int i=0;i<64;i++) {
// 	printf("Channel %d : %d %d \n",i,rawStream->GetADC(i),rawStream->GetTime(i)); }
//   treeD->Print(); printf(" \n"); 
  }

  loader->WriteDigits("OVERWRITE");
  loader->UnloadDigits();	
	
  delete rawReader;
  delete rawStream;

  timer.Stop();
  timer.Print();
}
