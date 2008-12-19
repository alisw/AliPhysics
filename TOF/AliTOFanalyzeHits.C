Int_t AliTOFanalyzeHits(Int_t numberOfEvents=0, Bool_t drawing=kFALSE)
{

  /////////////////////////////////////////////////////////////////////////
  //
  // Analyzes TOF hits and fills QA-histograms
  // and writes the histograms in the TOF_hitsQA.root file
  //
  // Use case:
  // start aliroot
  // root [0] .L AliTOFanalyzeHits.C
  // root [1] AliTOFanalyzeHits()
  //
  // By default, it analyzes hits for all the events in the header file
  // and does not draw the histograms filled
  //
  // If you want analyze hits only the 1th event
  // you can use the following line:
  //
  // root[0] .L AliTOFanalyzeHits.C
  // root[1] AliTOFanalyzeHits(1)
  //
  // Updated to the new I/O: C. Zampolli
  //
  // Report problems to decaro@sa.infn.it
  //
  /////////////////////////////////////////////////////////////////////////

  Int_t rc = 0;
  
  // Define the histograms
  // x,y,z, rho, tof, padx, padz, sector, plate, strip, (x vs y)

  // hit-map in a plane
  TH2F *h2hitMap = new TH2F("h2hitMap","Hit Map (projection on the plane)",2500,-12500.,12500.,800,-400.,400.);
  // time of flight distribution for primaries and secondaries
  TH1F *htofp    = new TH1F("htofp","Time of Flight (primaries)",800,0.,80.);
  TH1F *htofs    = new TH1F("htofs","Time of Flight (secondaries)",800,0.,80.);
  // momentum when striking the TOF for primaries and secondaries
  TH1F *htofmomp = new TH1F("htofmomp","Momentum at TOF (primaries)",100,0.,10.);
  TH1F *htofmoms = new TH1F("htofmoms","Momentum at TOF (secondaries)",100,0.,10.);
  // TOF hit volumes
  TH1F *hsector  = new TH1F("hsector","Sector",18,0.,18.);
  TH1F *hplate   = new TH1F("hplate","Plate ", 5,0., 5.);
  TH1F *hstrip   = new TH1F("hstrip","Strip ",20,0.,20.);
  TH1F *hpadz    = new TH1F("hpadz","Pad along z ",2,0.,2.);
  TH1F *hpadx    = new TH1F("hpadx","Pad along x",48,0.,48.);
  // track length when striking the TOF (used by AliTOFT0)
  TH1F *htrackLenp= new TH1F("htrackLenp","Track Length on TOF for Primaries",800,0.,800.);

  // Histograms added to control the right TOF element numbering:
  // it should be increasing with the azimuthal and polar angles

  TH2F *hmoduleVStheta     = new TH2F("hmoduleVStheta", "hmoduleVStheta", 180,0.,180.,5,0,5);
  TH2F *hsectorVSphi       = new TH2F("hsectorVSphi", "hsectorVSphi", 360,0.,360.,18,0,18);
  TH2F *hstripVStheta      = new TH2F("hstripVStheta", "hstripVStheta", 180,0.,180.,20,0,20);
  //TH2F *hpadzVStheta       = new TH2F("hpadzVStheta", "hpadzVStheta", 180,0.,180.,2,0,2);
  TH2F *hpadxVSphi         = new TH2F("hpadxVSphi", "hpadxVSphi", 360,0.,360.,48,0,48);
  TH2F *hpadz2stripVStheta = new TH2F("hpadz2stripVStheta", "hpadz2stripVStheta", 180,0.,180.,40,0,40);
  //TH2F *hdzVSpadz2strip    = new TH2F("hdzVSpadz2strip", "hdzVSpadz2strip",40,0,40,70,-3.5,3.5);

  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }

  if (gAlice)
    {
      delete AliRunLoader::GetRunLoader();
      delete gAlice;
      gAlice = 0x0;
    }

  AliRunLoader *rl = AliRunLoader::Open("galice.root",AliConfig::GetDefaultEventFolderName(),"read");
  if (!rl)
    {
      cerr<<"Can't load RunLoader from file! \n";
      rc = 1;
      return rc;
    }

  rl->LoadgAlice();
  gAlice = rl->GetAliRun();

  if (!gAlice)
    {
      cerr << "<AliTOFanalyzeHits> AliRun object not found on file \n";
      rc = 2;
      return rc;
    }

  rl->LoadHeader();

  // Get the pointer to the TOF detector
  AliLoader *tofl = rl->GetLoader("TOFLoader");
  AliTOF *tof = (AliTOF *) gAlice->GetDetector("TOF");
  if (tof == 0x0 || tofl == 0x0) {
    cerr << "<AliTOFanalyzeHits> Can not find TOF or TOFLoader \n";
    rc = 3;
    return rc;
  }

  Int_t countHits = 0;

  if (numberOfEvents==0) numberOfEvents=(Int_t)(rl->GetNumberOfEvents());

  for (Int_t ievent=0; ievent<numberOfEvents; ievent++) {
    printf ("Processing event %d \n", ievent);
    rl->GetEvent(ievent);

    // Get the pointer Hit tree
    tofl->LoadHits();
    TTree *hitTree = tofl->TreeH();
    tof->SetTreeAddress();
    if (!hitTree) {
      cout << "<AliTOFanalyzeHits> No TreeH found" << endl;
      rc = 4;
      return rc;
    }

    rl->LoadKinematics();
    //AliStack* stack = rl->Stack(); // it is not necessary to use the stack!

    // Get the number of entries in the hit tree
    // (Number of primary particles creating a hit somewhere)
    Int_t nTrack = (Int_t) hitTree->GetEntries();
    cout << "<AliTOFanalyzeHits> Found " << nTrack 
         << " primary particles with hits \n";

    Int_t nPrimaryOnTof = 0;
    Int_t nSecondaryOnTof = 0;
    Int_t nelectron  = 0;
    Int_t npion      = 0;
    Int_t nkaon      = 0;
    Int_t nproton    = 0;    
    Int_t nmuon      = 0;
  
    // Loop through all entries in the tree
    for (Int_t iTrack = 0; iTrack < nTrack; iTrack++) {

      tof->ResetHits();
      hitTree->GetEvent(iTrack);

      // Loop through the TOF hits  
      AliTOFhitT0 *hit = (AliTOFhitT0 *) tof->FirstHit(-1);
      while (hit) {

	countHits++;
	
        Float_t x     = hit->X();
        Float_t y     = hit->Y();
        Float_t z     = hit->Z();
        Float_t phiAngle=TMath::Pi() + TMath::ATan2(-y,-x);
        Float_t rhoRadius=TMath::Sqrt(x*x+y*y);
        Float_t thetaAngle=TMath::Pi() + TMath::ATan2(-rhoRadius,-z);
        Float_t dummy=rhoRadius*phiAngle;
        h2hitMap->Fill(dummy,z);

        phiAngle*=180./TMath::Pi();
        thetaAngle*=180./TMath::Pi();

        Float_t flightTime = hit->GetTof(); // [s]
        flightTime *= 1.e+09; // convert in [ns]
        Float_t angle = hit->GetIncA();
        Float_t tofmom = hit->GetMom(); // [GeV/c]
        Float_t trackLen = hit->GetLen(); // [cm]

        // TOF hit volumes
        Int_t sector = hit->GetSector(); // range [1-18]
        Int_t plate  = hit->GetPlate();  // range [1- 5]
        Int_t strip  = hit->GetStrip();  // range [1-20]
        Int_t padz   = hit->GetPadz();   // range [1- 2]
        Int_t padx   = hit->GetPadx();   // range [1-48]
        // it is QA, then I perform QA!

        Bool_t isHitBad = (sector<0 || sector>17 || 
			    plate<0 || plate>4   || 
			     padz<0 || padz>1    || 
			     padx<0 || padx>47   ||
			   ((strip<0 || strip>14) && plate == 2) ||
			   ((strip<0 || strip>18) && (plate == 1 || plate == 3)) ||
			   ((strip<0 || strip>19) && (plate == 0 || plate == 4)));

	if (isHitBad) {
          cout << "<AliTOFanalyzeHits>  strange hit found \n";
	  cout << "sector = " << sector <<
	          " plate = " << plate <<
	          " strip = " << strip <<
	          " padx = " << padx <<
	          " padz = " << padz << endl;
	  rc = 5;
	  return rc;
        }
	
        hmoduleVStheta->Fill(thetaAngle,plate);
	hstripVStheta->Fill(thetaAngle,strip);
        hsectorVSphi->Fill(phiAngle,sector);
        //hpadzVStheta->Fill(thetaAngle,padx);
        hpadxVSphi->Fill(phiAngle,padx);

	Float_t dummy2 = 2*strip + padz;
	hpadz2stripVStheta->Fill(thetaAngle,dummy2);

	/*
	  Float_t dummy3;
	  if (hit->GetDz()<=0) dummy3 = hit->GetDz()-1.75;
	  else dummy3 = hit->GetDz()+1.75;	
	  hdzVSpadz2strip->Fill(dummy2,dummy3);
	*/

        // filling hit volume histos
        hsector->Fill(sector);
        hplate->Fill(plate);
        hstrip->Fill(strip);
        hpadx->Fill(padx);
        hpadz->Fill(padz);
/*
        Int_t track = hit->Track();
        TParticle *part = gAlice->GetMCApp()->Particle(track);

        // getting MC info for the current track
        if (part->GetFirstMother()<0){
          Int_t icod = TMath::Abs(part->GetPdgCode());
 	  switch (icod) {
	  case 211:
	    npion++;
	    break ;
	  case 321:
	    nkaon++;
	    break ;
	  case 2212:
	    nproton++;
	    break ;
	  case 11:
	    nelectron++;
	    break ;
	  case 13:
	    nmuon++;
	    break ;
	  }
	  htofp->Fill(flightTime);
	  htofmomp->Fill(tofmom);
          htrackLenp->Fill(trackLen);
        } else {
	  htofs->Fill(flightTime);
	  htofmoms->Fill(tofmom);
        }
*/
        // go to next hit
        hit = (AliTOFhitT0 *) tof->NextHit();

      }

    }

    tofl->UnloadHits();
    rl->UnloadKinematics();

    cout << "<AliTOFanalyzeHits> Found " << countHits << " hits in total \n";
    cout << npion     << " primary pions reached the TOF detector \n";
    cout << nkaon     << " primary kaons reached the TOF detector \n";
    cout << nproton   << " primary protons reached the TOF detector \n";
    cout << nelectron << " primary electrons reached the TOF detector \n";
    cout << nmuon     << " primary muons reached the TOF detector \n";

  }

  cout << "hpadx->GetEntries() = " << hpadx->GetEntries() << endl;

  rl->UnloadHeader();
  rl->UnloadgAlice();

  if (drawing) {  
    TCanvas *cHits = new TCanvas("cHits","AliTOFanalyzeHits hit volumes",50,50,900,900);
    cHits->Divide(3,2);
    cHits->cd(1);
    hsector->Draw();
    cHits->cd(2);
    hplate->Draw();
    cHits->cd(3);
    hstrip->Draw();
    cHits->cd(4);
    hpadz->Draw();
    cHits->cd(5);
    hpadx->Draw();
    
    TCanvas *chitmap = new TCanvas("chitmap","AliTOFanalyzeHits Hit Map",50,50,600,600);
    chitmap->cd();
    h2hitMap->Draw();
    
    TCanvas *ctrackLen = new TCanvas("ctrackLen","AliTOFanalyzeHits Track Length for primaries on TOF",50,50,400,400);
    ctrackLen->cd();
    htrackLenp->Draw();
    
    TCanvas *ctofmom = new TCanvas("ctofmom","AliTOFanalyzeHits flight times",50,50,700,700);
    ctofmom->Divide(2,2);
    ctofmom->cd(1);
    gPad->SetLogy();
    htofp->Draw();
    ctofmom->cd(2);
    gPad->SetLogy();
    htofs->Draw();
    ctofmom->cd(3);
    gPad->SetLogy();
    htofmomp->Draw();
    ctofmom->cd(4);
    gPad->SetLogy();
    htofmoms->Draw();
  }
  
  // save histos into file TOF_hitsQA.root
  TFile *fout = new TFile("TOF_hitsQA.root","RECREATE");
  h2hitMap->Write();
  htofp->Write();
  htofs->Write();
  htofmomp->Write();
  htofmoms->Write();
  hsector->Write();
  hplate->Write();
  hstrip->Write();
  hpadz->Write();
  hpadx->Write();
  htrackLenp->Write();

  hmoduleVStheta->Write();
  hsectorVSphi->Write();
  hstripVStheta->Write();
  //hpadzVStheta->Write();
  hpadxVSphi->Write();
  hpadz2stripVStheta->Write();
  //hdzVSpadz2strip->Write();

  fout->Close(); 

  cout << " Finished AliTOFanalizeHits \n";

  if (gAlice)
    {
      delete AliRunLoader::GetRunLoader();
      delete gAlice;
      gAlice = 0x0;
    }

  return rc;

}
