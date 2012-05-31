void readRec() 
{
  //read T0 RecPoints and plots histos 

  // Dynamically link some shared libs
  /*  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }
  */
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


  TH1F *hBestTimeC = new TH1F("hBestTimeC","First time Cside",
				  100,12.,13.);
  TH1F *hBestTimeA = new TH1F("hBestTimeA","First time Aside",100,12.,13.);
  TH1F *hMean= new TH1F("hMean"," T0 ",100,12.,13.);
  TH1F *hAcc = new TH1F("hAcc","rec - real vertex",100,-10,10);
				 
  TH1F *hRealVertex = new TH1F("hRealVertex","Real Vertex",100,-15,15);
  
  TH1F *hVertex = new TH1F("hVertex","Z position of vertex",   100,-15,15);
 
  TH1F *hAmp = new TH1F("hAmp"," amplitude",   100, 10,1000);
  TH1F *hTime = new TH1F("hTime"," time",   100,12000.,20000);

  TArrayI *amp = new TArrayI(24);  
  TArrayI *time = new TArrayI(24);  

  AliT0RecPoint *fRec   ; // digits
  fRec = new AliT0RecPoint();

 // Event ------------------------- LOOP  
  for (Int_t ievent=0; ievent<iNevents; ievent++){
    // for (Int_t ievent=0; ievent<990; ievent++){
    rl->GetEvent(ievent);

    AliHeader *header = gAlice->GetHeader();
    AliGenEventHeader* genHeader = header->GenEventHeader();
    TArrayF *o = new TArrayF(3); 
    genHeader->PrimaryVertex(*o);
    Float_t zRealVertex=o->At(2);
    hRealVertex->Fill( zRealVertex);

    lstart->LoadRecPoints("READ");
    TTree *recTree =  lstart->TreeR();
    TBranch *brRec=recTree->GetBranch("T0");
    AliT0RecPoint *fRec = new AliT0RecPoint();
    if (brRec) {
      brRec->SetAddress(&fRec);
    }else{
      cerr<<"EXEC Branch T0 Recpoints not found"<<endl;
      return;
    }
    brRec->GetEntry(0);
    Int_t   mean = fRec->GetMeanTime();
    hMean->Fill(mean/1000.);
    Int_t   besttimeleft = fRec->GetBestTimeLeft();
    Int_t   besttimeright = fRec->GetBestTimeRight();
    hBestTimeC->Fill(0.001 * besttimeright);
    hBestTimeA->Fill(0.001 * besttimeleft );
    Float_t vertex= fRec->GetVertex();
    if(vertex<99){
    cout<<ievent<<" "<<mean<<" real vertex "<< zRealVertex<<" vertex "<<vertex<<
      " a "<< besttimeleft<<" c "<< besttimeright<<endl;
    hAcc->Fill(zRealVertex-vertex);
    hVertex->Fill(vertex);
     for (Int_t i=0; i<24; i++){ 
      hAmp->Fill(fRec->GetAmp(i));
      hTime->Fill(fRec->GetTime(i));
      //  cout<<"time "<<fRec->GetTime(i)<<" amp "<<fRec->GetAmp(i)<<endl;
    } 
    }
  }
  Hfile = new TFile("FigRec.root","RECREATE","Histograms for T0 
digits");
  printf("Writting histograms to root file \n");
  Hfile->cd();
  //Create a canvas, set the view range, show histograms
  gStyle->SetOptStat(111111);
  //  TCanvas *c1 = new TCanvas("c1","Alice T0 Time ",400,10,600,600);
 //  hTimediff->Write();
  hBestTimeC->Write(); 
  hBestTimeA ->Write(); 
  hVertex->Write();
  hRealVertex->Write();
  hAcc->Write();  
  hMean->Write();
  hAmp->Write();
  hTime->Write();
} // end of macro




