
void ReadACORDEDigits()

{

  // get run loader loader
   AliRunLoader* rl = 
     AliRunLoader::Open("galice.root",
			AliConfig::GetDefaultEventFolderName(),"read");
   if (rl == 0x0) {
     cerr<<"Can not open session for file galice.root\n";
     return;
   }
   
   rl->LoadgAlice();
   gAlice = rl->GetAliRun();

   // get acorde loader
   AliACORDE* acorde = (AliACORDE*)gAlice->GetDetector("ACORDE");
   AliLoader* aloader =rl->GetLoader("ACORDELoader");

   // loop over events
   Int_t nevt = rl->GetNumberOfEvents();
   cout << " There are " << nevt << " evts " << endl;
   for (Int_t ievt=0; ievt<nevt;ievt++) {
     rl->GetEvent(ievt);
     aloader->LoadDigits("READ");
     TTree* treeD = aloader->TreeD();
     TClonesArray *adigits = new TClonesArray ("AliACORDEdigit", 1000);
     treeD->GetBranch("ACORDEdigit")->SetAddress(&adigits);

     // loop over entries
     Int_t nent = treeD->GetEntries();
     cout << " There are " << nent << " entries in event " << ievt << endl;

     for (Int_t ient=0;ient<nent;ient++) {
       acorde->ResetDigits();
       treeD->GetEvent(ient);
       Int_t ndig = adigits->GetEntriesFast();
       cout << " There are " << ndig << " digits in entry " << ient << endl;
       for (Int_t idig=0;idig<ndig;idig++) {
	 AliACORDEdigit* digit = (AliACORDEdigit*) adigits->At(idig);
	 Int_t mod = digit->GetModule();
	 Float_t time = digit->GetTime();
	 cout << " Digit " << idig 
	      <<" : module="<< mod << " time=" << time << endl;
	 for (Int_t i=0;i<3;i++)
	   cout << " track " << i << " is " << digit->GetTrack(i) << endl;
       } // end loop over digits
     } // end loop over entries
     aloader->UnloadDigits();
   } // end loop over events
}
