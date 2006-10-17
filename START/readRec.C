void readRec() 
{
  //read START RecPoints and plots histos 

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
    
    AliSTART* START  = (AliSTART *)gAlice->GetDetector("START");
    
    rl->LoadHeader();
    Int_t retval;
    AliLoader* lstart = rl->GetLoader("STARTLoader");

  Int_t iNevents=rl->GetNumberOfEvents();
  cout<<"  nevents   "<<iNevents<<endl;


  TH1F *hBestTimeC = new TH1F("hBestTimeC","First time Cside",
				  100,12.,13.);
  TH1F *hBestTimeA = new TH1F("hBestTimeA","First time Aside",100,12.,13.);
  TH1F *hMean= new TH1F("hMean"," T0 ",100,12.,13.);
  TH1F *hAcc = new TH1F("hAcc","rec - real vertex",100,-10,10);
				 
  TH1F *hRealVertex = new TH1F("hRealVertex","Real Vertex",100,-15,15);
  
  TH1F *hVertex = new TH1F("hVertex","Z position of vertex",   100,-15,15);
  AliSTARTRecPoint *fRec   ; // digits
  fRec = new AliSTARTRecPoint();

 // Event ------------------------- LOOP  
  for (Int_t ievent=0; ievent<iNevents; ievent++){
    rl->GetEvent(ievent);

    AliHeader *header = gAlice->GetHeader();
    AliGenEventHeader* genHeader = header->GenEventHeader();
    TArrayF *o = new TArrayF(3); 
    genHeader->PrimaryVertex(*o);
    Float_t zRealVertex=o->At(2);
    hRealVertex->Fill( zRealVertex);

    lstart->LoadRecPoints("READ");
    TTree *recTree =  lstart->TreeR();
    TBranch *brRec=recTree->GetBranch("START");
    AliSTARTRecPoint *fRec = new AliSTARTRecPoint();
    if (brRec) {
      brRec->SetAddress(&fRec);
    }else{
      cerr<<"EXEC Branch START Recpoints not found"<<endl;
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
    cout<<ievent<<" "<<mean<<" real vertex "<< zRealVertex<<" vertex "<<vertex<<
      " a "<< besttimeleft<<" c "<< besttimeright<<endl;
    hAcc->Fill(zRealVertex-vertex/2.);
    hVertex->Fill(vertex/2.);
  }
  Hfile = new TFile("FigRec.root","RECREATE","Histograms for START 
digits");
  printf("Writting histograms to root file \n");
  Hfile->cd();
  //Create a canvas, set the view range, show histograms
  gStyle->SetOptStat(111111);
  //  TCanvas *c1 = new TCanvas("c1","Alice START Time ",400,10,600,600);
 //  hTimediff->Write();
  hBestTimeC->Write(); 
  hBestTimeA ->Write(); 
  hVertex->Write();
  hRealVertex->Write();
  hAcc->Write();  
 hMean->Write();
} // end of macro




