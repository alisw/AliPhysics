// ----------------------------------------------------//
//                                                     //
//    This macro does Raw to SDigits                   //
//                                                     //
// ----------------------------------------------------//


void AliPMDRaw2SDigits(Int_t nevt = 1) 
{
  TStopwatch timer;
  timer.Start();

  // Open the AliRoot file

  AliRunLoader *runLoader = AliRunLoader::Open("galice.root");
  if (!runLoader)
    {
      cerr<<"Can't load RunLoader"<<endl;
      return 1;
    }

  AliPMDRawToSDigits pmdr2sd;

  Int_t ievt = 0;

  for (ievt = 0; ievt < nevt; ievt++)
    {
      runLoader->GetEvent(ievt);
      AliRawReaderFile reader(ievt);
      pmdr2sd.Raw2SDigits(runLoader, &reader);
    }

  timer.Stop();
  timer.Print();
}

