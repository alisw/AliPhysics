void AliTRDDigits2Raw() {

  const char * inFile = "galice.root";
  AliRunLoader *rl = AliRunLoader::Open(inFile,"Event","read");

  rl->GetEvent(0);
  AliLoader *trdloader=rl->GetLoader("TRDLoader");
  trdloader->LoadDigits();
  TTree *digitsTree=trdloader->TreeD();

  digitsTree->Print();

  AliTRDrawData *raw = new AliTRDrawData();
  raw->SetDebug(2);

  if (raw->Digits2Raw(digitsTree)) {
    printf("Digits2Raw done! \n");
  } else {
    printf("Digits2Raw returned FALSE! \n");
  }

  delete rl;

}

