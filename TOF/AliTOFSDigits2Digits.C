Int_t AliTOFSDigits2Digits(TString infileNameSDigits, Int_t firstEvent=0,Int_t nEvents=1, Int_t ndump=15)
{
  // Create TOF digits out of TOF sdigits (no merging implemented here).
  // Input: infileNameSDigits --> TOF sdigits filename 
  //        firstEvent        --> first event to digitize
  //        nEvents           --> number of events to digitize 
  //                              including the first and the last ones
  // if nEvents==0 we sdigitize all events in infileNameSDigits
  //        ndump             --> number of dumped digits

  // Author: F. Pierella (Bologna University)
  // report problem to pierella@bo.infn.it

  // Use case: (start root)
  //root [0] .L AliTOFSDigits2Digits.C                                    
  //root [1] AliTOFSDigits2Digits("fileWithTOFSdigits.root")


  // number 3 is a legacy from AliDigit object
  const Int_t kMAXDIGITS = 3;


  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  } else {
    delete gAlice;
    gAlice = 0;
  }
  
  Int_t rc=0;
  
  
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(infileNameSDigits.Data());
  if(file){
    cout<<"headerFile already open \n";
  }
  else {
    // file is open in "update" mode
    // in order to have a writable digits tree
    if(!file)file=TFile::Open(infileNameSDigits.Data(),"update");
  }
  
  // Get AliRun object from file
  if (!gAlice) {
    gAlice = (AliRun*)file->Get("gAlice");
    if (gAlice) printf("AliRun object found on file\n");
  }

  if (nEvents == 0) nEvents = (Int_t) gAlice->TreeE()->GetEntries();
  
  AliTOF * tof = (AliTOF *) gAlice->GetDetector("TOF") ;

  if (!tof) {
    cout << "<AliTOFSDigits2Digits> No TOF detector found" << endl;
    rc = 2;
    return rc;
  }

  Int_t totndig=0;   // total number of digits
  Int_t tottracks=0; // total number of tracks contributing to totndig
  Int_t upperBoundEvNum=firstEvent+nEvents;
  for (Int_t ievent = firstEvent; ievent < upperBoundEvNum; ievent++) {

    gAlice->GetEvent(ievent);
    if (gAlice->TreeD () == 0)
      gAlice->MakeTree ("D");
    
    //Make branches
    char branchname[20];
    sprintf (branchname, "%s", tof->GetName ());
    //Make branch for digits
    tof->MakeBranch ("D");


    // get the TOF branch in TreeS for current event
    char tname[100]; sprintf(tname,"TreeS%d",ievent);
    
    TTree *sdigitstree=(TTree*)file->Get(tname);                   
    
    TClonesArray * fSDigits= new TClonesArray("AliTOFSDigit",  1000); 
    sdigitstree->GetBranch("TOF")->SetAddress(&fSDigits);           
    
    Int_t nEntries = sdigitstree->GetEntries();                                
    //cout << nEntries << endl;
    
    // Loop through all entries in the tree
    Int_t nbytes;

    for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {
      
      // Import the tree
      nbytes += sdigitstree->GetEvent(iEntry);
      
      // Get the number of sdigits
      Int_t ndig = fSDigits->GetEntriesFast();
      cout << "------------------<AliTOFSDigits2Digits>------------------" << endl;
      cout << "Found " << ndig << " TOF SDigits for event " << ievent << endl;
      cout << "------------------------------------------------------------" << endl;

      // start loop on sdigits
      for (Int_t k = 0; k < ndig; k++) {

	Int_t    vol[5];       // location for a digit
	
	// Get the information for this digit
	AliTOFSDigit *tofsdigit = (AliTOFSDigit *) fSDigits->UncheckedAt(k);

	Int_t nslot=tofsdigit->GetNDigits(); // get the number of slots
	                                     // for current sdigit

	// TOF sdigit volumes (always the same for all slots)
	Int_t sector    = tofsdigit->GetSector(); // range [1-18]
	Int_t plate     = tofsdigit->GetPlate();  // range [1- 5]
	Int_t strip     = tofsdigit->GetStrip();  // range [1-20]
	Int_t padz      = tofsdigit->GetPadz();   // range [1- 2]
	Int_t padx      = tofsdigit->GetPadx();   // range [1-48]

	vol[0] = sector;
	vol[1] = plate;
	vol[2] = strip;
	vol[3] = padx;
	vol[4] = padz;

	//--------------------- QA section ----------------------
	// in the while, I perform QA
	Bool_t isSDigitBad = (sector<1 || sector>18 || plate<1 || plate >5 || padz<1 || padz>2 || padx<1 || padx>48);
	
	if (isSDigitBad) {
	  cout << "<AliTOFSDigits2Digits>  strange sdigit found" << endl;
	  rc = 3;
	  return rc;
	}
	//-------------------------------------------------------

	//------------------- Dump section ----------------------
	if(k<ndump){
	  cout << k << "-th | " << "Sector " << sector << " | Plate " << plate << " | Strip " << strip << " | PadZ " << padz << " | PadX " << padx << endl;
	  cout << k << "-th sdigit" << endl;
	  cout << "----------------------------------------------------"<< endl;
	}
	// ------------------------------------------------------

	// start loop on number of slots for current sdigit
	for (Int_t islot = 0; islot < nslot; islot++) {
	  Float_t  digit[2];     // TOF digit variables
	  Int_t tracknum[kMAXDIGITS];     // contributing tracks for the current slot
	
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
	  tof->AddDigit(tracknum, vol, digit);
	  totndig++;
	}


      } // end loop on sdigits
      
    } // end loop on entries

    // free used memory
    fSDigits->Clear();
    fSDigits=0;

    gAlice->TreeD()->Reset();
    gAlice->TreeD()->Fill();
    //gAlice->TreeS()->Write(0,TObject::kOverwrite) ;
    gAlice->TreeD()->AutoSave();

  } // end loop on events 

  cout << "----------------------------------------------------------" << endl;
  cout << "<AliTOFSDigits2Digits> Summary" << endl;
  cout << "contributing tracks to " << totndig << " digits: " << tottracks << endl; 
  cout << "----------------------------------------------------------" << endl;

  return rc;
  

}
