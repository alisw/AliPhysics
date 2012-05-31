//_____________________________________________________//
//                                                     //
//    This macro reads ACORDE DDL Raw Data and          //
//    converts it into Digits                          //
//                                                     //
//____________________________________________________ //


void ACORDERaw2Digits(Int_t nEvents = 1, char* fileName = "alien:///alice/data/2008/LHC08a_ACORDE/000016788/raw/08000016788014.20.root")
{
  // Reads DDL data from fileName

  TStopwatch timer;
  timer.Start();

  TGrid::Connect("alien://");

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
    for (Int_t j=0; j<4; j++)
      printf(" %x",rawStream->GetWord(j));
    printf("\n");
  }

  delete rawReader;
  delete rawStream;

  timer.Stop();
  timer.Print();
}
