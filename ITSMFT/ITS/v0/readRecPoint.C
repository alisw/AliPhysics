void readRecPoint(){
  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");
  gSystem->Load("libITSUpgradeRec");
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1111111);
  Int_t nbins=100;
  Int_t xmin=0;
  Int_t xmax=50000;//00*1e-09;


  AliITSsegmentationUpgrade *seg = new AliITSsegmentationUpgrade();
  if(!seg){
    printf("no segmentation info available... Exiting");
    return;
  }
  Int_t nLayers = seg->GetNLayers();

  TH1D **hNel, **hNsect;
  hNel = new TH1D[nLayers]; 
  hNsect = new TH1D[nLayers]; 
 // [nLayers];
  for(Int_t i=0; i< nLayers; i++ ) {
  hNel[i] = new TH1D(Form("hNel%i",i),Form("cluster charge distribution [ Layer %i] ",i),nbins,xmin,xmax);
  hNel[i]->SetXTitle("N electrons");
  hNsect[i] = new TH1D(Form("hNsect%i",i),Form("cluster entries per sector [ Layer %i] ",i),seg->GetNSectors(),-0.5,seg->GetNSectors()-0.5);
  hNsect[i]->SetXTitle("Sector Number");
  hNsect[i]->SetYTitle("# clusters");
  hNsect[i]->SetMinimum(0);
  }

  TH1D * type = new TH1D("hCluType"," cluster type" , 50,0,15 );

  TH2F *xyGlob = new TH2F("xyGlob"," X - Y Global coordinates ",100,-50,50,100,-50,50);
  xyGlob->SetXTitle("cm"); 
  xyGlob->SetMarkerStyle(7); 
  TH1F *zGlob  = new TH1F("zGlob", " Z Global coordinates ",200, -50,50 );
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

  TClonesArray statITSCluster("AliITSRecPointU");

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
      AliITSRecPointU *recp = (AliITSRecPointU*)ITSCluster->UncheckedAt(i);
      //  cout<<"layer "<<recp->GetLayer()<<endl;
      Double_t xyz[3]={-1,-1,-1};
      seg->DetToGlobal(recp->GetLayer(), recp->GetModule(), recp->GetDetLocalX(), recp->GetDetLocalZ(), xyz[0],xyz[1],xyz[2]) ;
      xyGlob->Fill(xyz[0],xyz[1]);
      zGlob->Fill(xyz[2]);
      charge=recp->GetQ();
      type->Fill(recp->GetType());
      hNel[recp->GetLayer()]->Fill(charge);
      hNsect[recp->GetLayer()]->Fill(recp->GetModule());
    }
  }
  
  Int_t size = 400;

  TCanvas *xyCanv =  new TCanvas("xvCanvClus","RecPoint X-Y Positions",10,10,size,size);
  xyCanv->cd();
  xyGlob->Draw();

  TCanvas *zCanv =  new TCanvas("zCanvClus","RecPoint Z Positions",size+20,10,size,size);
  zCanv->cd();
  zGlob->Draw();
  TCanvas *typeCanv = new TCanvas("ClusType","Cluster type distribution",10,size+40,2*size+10,size);
  type->Draw();


  TCanvas *c = new TCanvas("c","Cluster charge distribution",900,0, 1000,550);
   c->Divide(3,(nLayers/3)+1);
  for(Int_t ip =1; ip<=nLayers; ip++){
    c->cd(ip);
    hNel[ip-1]->Draw();
  }
 
  TCanvas *cS = new TCanvas("cS","clusters in sectors",900,640,1000,550);
   cS->Divide(3,(nLayers/3)+1);
  for(Int_t ip =1; ip<=nLayers; ip++){
    cS->cd(ip);
    hNsect[ip-1]->Draw();
  } 
}

