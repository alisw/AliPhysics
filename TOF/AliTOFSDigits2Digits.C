Int_t AliTOFSDigits2Digits(Int_t numberOfEvents = 0) {

  /////////////////////////////////////////////////////////////////////////
  //
  // Creates TOF digits from the summable digits for all event in the header file
  //
  // Use case:
  // start root
  // // load the macro
  // root[0] .L AliTOFSDigits2Digits.C
  // root[1] AliTOFSDigits2Digits()
  //
  // By default, it creates digits for all the events in the header file.
  //
  // If you want create digits only the firts event
  // you can use the following lines:
  //
  // root[0] .L AliTOFSDigits2Digits.C
  // root[1] AliTOFSDigits2Digits(1)
  //
  // Created by: F. Pierella
  // Updated to the new I/O: I.Belikov (Jouri.Belikov@cern.ch)
  //
  // Report problems to decaro@sa.infn.it
  //
  /////////////////////////////////////////////////////////////////////////

  Int_t rc = 0;

  if (gAlice)
    {
      delete AliRunLoader::GetRunLoader();
      delete gAlice;
      gAlice = 0x0;
    }

  AliRunLoader* rl = AliRunLoader::Open("galice.root");
  if (rl == 0x0)
    {
      cerr<<"Can not open session"<<endl;
      rc = 1;   
      return rc;
    }

  if (rl->LoadgAlice())
    {
      cerr<<"Error occured while loading gAlice \n";
      rc = 2;
      return rc;
    }

  AliLoader *tofl = rl->GetLoader("TOFLoader");
   if (tofl == 0x0)
     {
       cerr<<"Can not get the TOF Loader \n";
       rc = 3;
       return rc;
     }

   gAlice=rl->GetAliRun();
   if (!gAlice)
     {
       cerr<<"Can't get gAlice !\n";
       rc = 4;
       return rc;
     }

   tofl->LoadSDigits("read");
   tofl->LoadDigits("recreate");

   Int_t totndig=0;   // total number of digits
   Int_t tottracks=0; // total number of tracks contributing to totndig

   if (numberOfEvents==0) numberOfEvents=rl->GetNumberOfEvents();

   TClonesArray *fSDigits=new TClonesArray("AliTOFSDigit",  1000); 
   TClonesArray *fDigits =new TClonesArray("AliTOFdigit",  1000);
   TClonesArray &da=*fDigits; 

   for (Int_t ievent = 0; ievent < numberOfEvents; ievent++) {
     rl->GetEvent(ievent);
     
     TTree *sTree=tofl->TreeS();
     if (sTree == 0)
       {
	 cerr<<"Can't get the sdigit tree !\n";
	 rc = 5;
	 return rc;
       }
     TBranch *branch=sTree->GetBranch("TOF");
     if (!branch)
       {
	 cerr<<"Cant' get the branch !\n";
	 rc = 6;
	 return rc;
       }
     branch->SetAddress(&fSDigits);

    TTree *dTree=tofl->TreeD();
    if (dTree == 0) {
       tofl->MakeTree("D");
       dTree=tofl->TreeD();
    }
    branch=dTree->GetBranch("TOF");
    if (!branch) dTree->Branch("TOF",&fDigits);    
    else branch->SetAddress(&fDigits);

    Int_t nEntries = sTree->GetEntries();                                
    for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {
      sTree->GetEvent(iEntry);

      Int_t ndig = fSDigits->GetEntriesFast();
      cout << "----------------<AliTOFSDigits2Digits>---------------- \n";
      cout << "Found " << ndig << " TOF SDigits for event " << ievent << endl;
      cout << "------------------------------------------------------ \n";

      for (Int_t k = 0; k < ndig; k++) {
	Int_t    vol[5];       // location for a digit

	// Get the information for this digit
	AliTOFSDigit *tofsdigit = (AliTOFSDigit *)fSDigits->UncheckedAt(k);

	Int_t nslot=tofsdigit->GetNDigits(); // get the number of slots
	                                     // for current sdigit

	// TOF sdigit volumes (always the same for all slots)
	Int_t sector    = tofsdigit->GetSector(); // range [0-17]
	Int_t plate     = tofsdigit->GetPlate();  // range [0- 4]
	Int_t strip     = tofsdigit->GetStrip();  // range [0-19]
	Int_t padz      = tofsdigit->GetPadz();   // range [0- 1]
	Int_t padx      = tofsdigit->GetPadx();   // range [0-47]

	vol[0] = sector;
	vol[1] = plate;
	vol[2] = strip;
	vol[3] = padx;
	vol[4] = padz;

	//--------------------- QA section ----------------------
	// in the while, I perform QA
	Bool_t isSDigitBad = (sector<0 || sector>17 || 
                               plate<0 || plate >4  || 
                                padz<0 || padz>1    || 
                                padx<0 || padx>47);
	
	if (isSDigitBad)
	  {
	    cout << "<AliTOFSDigits2Digits>  strange sdigit found \n";
	    rc = 7;
	    return rc;
	  }
	//-------------------------------------------------------

	// start loop on number of slots for current sdigit
	for (Int_t islot = 0; islot < nslot; islot++) {
	  Float_t  digit[2];           // TOF digit variables
          const Int_t kMAXDIGITS = 3;  // number 3 is a legacy from AliDigit object

	  Int_t tracknum[kMAXDIGITS];  //contributing tracks for the current slot

	  Float_t tdc=tofsdigit->GetTdc(islot); digit[0]=tdc;
	  Float_t adc=tofsdigit->GetAdc(islot); digit[1]=adc;

	  tracknum[0]=tofsdigit->GetTrack(islot,0);
	  tracknum[1]=tofsdigit->GetTrack(islot,1);
	  tracknum[2]=tofsdigit->GetTrack(islot,2);

	  for (Int_t i = 0; i < kMAXDIGITS; i++) {
	    tottracks++;
	    // search for the first empty location
	    if(tracknum[i]==-1){
	      tottracks--;
	      break;
	    }
	  }

	  // adding a TOF digit for each slot
          {
          Int_t ndigits=da.GetEntriesFast();
	  new (da[ndigits]) AliTOFdigit(tracknum, vol, digit);
          }          
	  totndig++;
	}

      } // end loop on sdigits
      fSDigits->Clear();

    } // end loop on entries

    dTree->Fill();
    tofl->WriteDigits("OVERWRITE");

    // free used memory
    fDigits->Clear();

  } // end loop on events 

  delete fSDigits;
  delete fDigits;

  tofl->UnloadDigits();
  tofl->UnloadSDigits();
  rl->UnloadHeader();
  rl->UnloadgAlice();

  cout << "---------------------------------------------------------- \n";
  cout << "<AliTOFSDigits2Digits> Summary \n";
  cout << "contributing tracks to " << totndig << " digits: " << tottracks << endl;
  cout << "---------------------------------------------------------- \n";

  if (gAlice)
    {
      delete AliRunLoader::GetRunLoader();
      delete gAlice;
      gAlice = 0x0;
    }

  return rc;

}
