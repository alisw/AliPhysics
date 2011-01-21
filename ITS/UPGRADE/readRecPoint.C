void readRecPoint(){
  gSystem->Load("libITSUpgradeSim");
  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeRec");
  gROOT->SetStyle("Plain");
  Int_t nbins=100;
  Int_t xmin=0;
  Int_t xmax=50000;//00*1e-09;

  const Int_t nLayers = 6;

  AliITSsegmentationUpgrade *seg = new AliITSsegmentationUpgrade();
  if(!seg){
    printf("no segmentation info available... Exiting");
    return;
  }

  TH1D *hNel[nLayers];
  for(Int_t i=0; i< nLayers; i++ ) {
  hNel[i] = new TH1D(Form("hNel%i",i),Form("cluster charge distribution [ Layer %i] ",i),nbins,xmin,xmax);
  hNel[i]->SetXTitle("N electrons");
  }
  TH1D * type = new TH1D("hCluType"," cluster type" , 50,0,15 );

  TH2F *xyGlob = new TH2F("xyGlob"," X - Y Global coordinates ",100,-50,50,100,-50,50);
  xyGlob->SetXTitle("cm"); 
 xyGlob->SetMarkerStyle(7); 
  TH1F *zGlob  = new TH1F("zGlob", " Z Global coordinates ",200, -100,100 );
  zGlob->SetXTitle("cm"); 


  gAlice=NULL;
  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
  runLoader->LoadgAlice();

  gAlice = runLoader->GetAliRun();

  runLoader->LoadHeader();
  runLoader->LoadKinematics();
  runLoader->LoadRecPoints();

  AliITSLoader *dl = (AliITSLoader*)runLoader->GetDetectorLoader("ITS");


  TTree *clusTree = 0x0;

  TClonesArray statITSCluster("AliITSRecPoint");

  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    runLoader->GetEvent(iEvent);
    clusTree=dl->TreeR();
    TClonesArray *ITSCluster = &statITSCluster;
    TBranch* itsClusterBranch=clusTree->GetBranch("ITSRecPoints");
    if (!itsClusterBranch) {
      printf("can't get the branch with the ITS clusters ! \n");
      return;
    }
    itsClusterBranch->SetAddress(&ITSCluster);
    clusTree->GetEntry(0);   
    Double_t charge=0.;
    Int_t nCluster = ITSCluster->GetEntriesFast();
    for(Int_t i=0; i<nCluster; i++){
      AliITSRecPoint *recp = (AliITSRecPoint*)ITSCluster->UncheckedAt(i);
      Double_t xyz[3]={-1,-1,-1};
      seg->DetToGlobal(recp->GetLayer(), recp->GetDetLocalX(), recp->GetDetLocalZ(), xyz[0],xyz[1],xyz[2]) ;
      xyGlob->Fill(xyz[0],xyz[1]);
      zGlob->Fill(xyz[2]);
      charge=recp->GetQ();
      // cout<< "layer "<< recp->GetLayer() << "   local system    X "<< recp->GetDetLocalX() << " Z "<< recp->GetDetLocalZ() <<endl;  
      type->Fill(recp->GetType());
      hNel[recp->GetLayer()]->Fill(charge);
    }

  }

  TCanvas *xyCanv =  new TCanvas("xvCanvClus","RecPoint X-Y Positions",500,500);
  xyCanv->cd();
  xyGlob->Draw();
  TCanvas *zCanv =  new TCanvas("zCanvClus","RecPoint Z Positions",500,500);
  zCanv->cd();
  zGlob->Draw();
  new TCanvas();
  type->Draw();

  TCanvas *c = new TCanvas("c","Cluster charge distribution",1000,800);
  c->Divide(3,2);
  for(Int_t ip =1; ip<=6; ip++){
    c->cd(ip);
    hNel[ip-1]->Draw();
  } 
}

