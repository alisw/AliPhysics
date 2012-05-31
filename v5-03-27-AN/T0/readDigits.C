void readDigits() 
{
 
  // Dynamically link some shared libs
  /*  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }
  */
    Float_t c = 0.0299792; // cm/ps
    Float_t channelWidth = 25;
    Float_t timeDelay = 150;
    char filename[100];
    sprintf(filename,"galice.root");
    AliRunLoader* rl = AliRunLoader::Open("galice.root",AliConfig::GetDefaultEventFolderName(),"read");
    if (rl == 0x0)
   {
     cerr<<"Can not open session for file galice.root\n";
     return;
   }
    
    rl->LoadgAlice();
    gAlice = rl->GetAliRun();
    
    AliT0* T0  = (AliT0 *)gAlice->GetDetector("T0");
    
    rl->LoadHeader();
    Int_t retval;
    AliLoader* lstart = rl->GetLoader("T0Loader");

  Int_t iNevents=rl->GetNumberOfEvents();
  cout<<"  nevents   "<<iNevents<<endl;


  TH1F *hTimediff = new TH1F("hTimediff","Time difference",100,5950,6050);
  TH1F *hBestTimeright = new TH1F("hBestTimeright","First time right",
				  100,8450.,8650);
  TH1F *hBestTimeleft = new TH1F("hBestTimeleft","First time left",
				 100,8450.,8650.);
  TH1F *hRealVertex = new TH1F("hRealVertex","Real Vertex",100,-15,15);
  
  TH1F *hVertex = new TH1F("hVertex","Z position of vertex",   100,-15,15);
  //  digits = new AliT0digit();
  AliT0digit *fDigits   ; // digits
  fDigits = new AliT0digit();

 // Event ------------------------- LOOP  
  for (Int_t j=0; j<iNevents; j++){
    rl->GetEvent(j);

    AliHeader *header = gAlice->GetHeader();
    AliGenEventHeader* genHeader = header->GenEventHeader();
    TArrayF *o = new TArrayF(3); 
    genHeader->PrimaryVertex(*o);
    Float_t zRealVertex=o->At(2);
    hRealVertex->Fill( zRealVertex);

    lstart->LoadDigits("READ");
    TTree *digitsTree =  lstart->TreeD();
    TBranch *brDigits=digitsTree->GetBranch("T0");
    AliT0digit *fDigits = new AliT0digit();
    if (brDigits) {
      brDigits->SetAddress(&fDigits);
    }else{
      cerr<<"EXEC Branch T0 digits not found"<<endl;
      return;
    }
    brDigits->GetEntry(0);
    Int_t   besttimeright = fDigits->BestTimeRight();
    Int_t   besttimeleft = fDigits->BestTimeLeft();
    Int_t   timeDiff = fDigits->TimeDiff();
    Int_t    sumMult=   fDigits->SumMult();
    hTimediff->Fill(timeDiff);
    hBestTimeright->Fill(besttimeright);
    hBestTimeleft->Fill(besttimeleft );
    Float_t vertex= (timeDiff* channelWidth - timeDelay*1000.)*c;
    cout<<j<<" "<<besttimeright<<" "<< besttimeleft<<" "<<timeDiff<<" "<<vertex<<endl;
    hVertex->Fill(vertex);
  }
  Hfile = new TFile("Figdigits.root","RECREATE","Histograms for T0 digits");
  printf("Writting histograms to root file \n");
  Hfile->cd();
  //Create a canvas, set the view range, show histograms
  gStyle->SetOptStat(111111);
  //  TCanvas *c1 = new TCanvas("c1","Alice T0 Time ",400,10,600,600);
   hTimediff->Write();
  hBestTimeright->Write(); 
  hBestTimeleft ->Write(); 
  hVertex->Write();
  hRealVertex->Write();
    
 
} // end of macro




