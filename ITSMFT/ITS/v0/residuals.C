


Int_t nLayers =0;
Int_t GetLayer(AliTrackReference *ref, TObject *obj);

void residuals(){

  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeRec");
  gSystem->Load("libITSUpgradeSim");

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

  nLayers = seg->GetNLayers();

  Int_t nbins=100;
  Int_t xmin = -10;
  Int_t xmax= 10;
  Int_t zmin = -3000;
  Int_t zmax= 3000;

  TH1D **hDiffx, **hDiffy, **hDiffz, **hDiffRphi, **hDiffR;
  
  hDiffx    = new TH1D[nLayers];
  hDiffy    = new TH1D[nLayers];
  hDiffz    = new TH1D[nLayers];
  hDiffRphi = new TH1D[nLayers];
  hDiffR    = new TH1D[nLayers];
    
  for(Int_t i=0; i<nLayers; i++){
    hDiffx[i] = new TH1D(Form("hDiffX%i",i),Form("#Delta X Glob (Clusters - track ref  [Layer %i])",i),100,-40,40);
    hDiffy[i] = new TH1D(Form("hDiffY%i",i),Form("#Delta Y Glob (Clusters - track ref  [Layer %i])",i),100,-40,40);
    if(i<4) hDiffz[i] = new TH1D(Form("hDiffZ%i",i),Form("#Delta Z Glob (Clusters - track ref  [Layer %i])",i),50,-300,300);
    else hDiffz[i] = new TH1D(Form("hDiffZ%i",i),Form("#Delta Z Glob (Clusters - track ref  [Layer %i])",i),100,-2000,2000);
    hDiffRphi[i] = new TH1D(Form("hDiffRphi%i",i),Form("#Delta R #phi [Layer %i]",i),200,-100,100);
    hDiffRphi[i]->SetXTitle("#mum");
    hDiffR[i]= new TH1D(Form("hR%i",i),Form("#Delta R (cluster-track ref) [Layer %i]",i),nbins,xmin/20.,xmax/20.);
    hDiffR[i]->SetXTitle("#mum");
  }

  gAlice=NULL;
  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
  runLoader->LoadgAlice();

  gAlice = runLoader->GetAliRun();

  runLoader->LoadHeader();
  runLoader->LoadKinematics();
  runLoader->LoadTrackRefs();
  runLoader->LoadRecPoints();

  AliLoader *dl = runLoader->GetDetectorLoader("ITS");

  //Trackf
  TTree *trackRefTree = 0x0; 
  TClonesArray *trackRef = new TClonesArray("AliTrackReference",1000);
  //RECPOINT
  TTree *clusTree = 0x0;
  TClonesArray statITSCluster("AliITSRecPointU");


  AliITSsegmentationUpgrade *segmentation2=segmentation2 = new AliITSsegmentationUpgrade();
  // Loop over events
  for(Int_t i=0; i<runLoader->GetNumberOfEvents(); i++){
    runLoader->GetEvent(i); 
    AliStack *stack = runLoader->Stack(); 
    trackRefTree=runLoader->TreeTR();
    TBranch *br = trackRefTree->GetBranch("TrackReferences");
    if(!br) {
      printf("no TR branch available , exiting \n");
      return;
    }
    br->SetAddress(&trackRef);

    printf("event : %i \n",i);
    // init the clusters tree 
    clusTree=dl->TreeR();
    TClonesArray *ITSCluster = &statITSCluster;
    TBranch* itsClusterBranch=clusTree->GetBranch("ITSRecPoints");
    if (!itsClusterBranch) {
      printf("can't get the branch with the ITS clusters ! \n");
      return;
    }

    itsClusterBranch->SetAddress(&ITSCluster);

    // init the trackRef tree 
    trackRefTree=runLoader->TreeTR();
    trackRefTree->SetBranchAddress("TrackReferences",&trackRef);

    //for(Int_t il=0; il<nLayers; il++){
    if (!clusTree->GetEvent(0))    continue;
    Int_t nCluster = ITSCluster->GetEntriesFast();
    for(Int_t in=0; in<nCluster; in++){
      AliITSRecPointU *recp = (AliITSRecPointU*)ITSCluster->UncheckedAt(in);
      Int_t il = recp->GetLayer();  
      Double_t xr,yr,zr = 0.;
      Double_t xz[2];
      xz[0]= recp->GetDetLocalX(); //gets fXloc
      xz[1]= recp->GetDetLocalZ();   //gets fZloc
      segmentation2->DetToGlobal(il,recp->GetModule(),xz[0], xz[1],xr,yr,zr);
      for(Int_t iLabel =0; iLabel<12; iLabel++){
	if(recp->GetTrackID(iLabel)<0) continue; 
	trackRefTree->GetEntry(stack->TreeKEntry(recp->GetTrackID(iLabel)));  	
	Int_t nref=trackRef->GetEntriesFast();
	for(Int_t iref =0; iref<nref; iref++){
	  Double_t x,y,z=0.;
	  Int_t labelTR=-999;
	  AliTrackReference *trR = (AliTrackReference*)trackRef->At(iref);
	  if(!trR) continue;
	  if(trR->DetectorId()!=AliTrackReference::kITS) continue;
	  if(!stack->IsPhysicalPrimary(recp->GetTrackID(iLabel))) continue;   
	  Int_t lay = GetLayer(trR, seg);    
	  if(TMath::Abs(lay-il)>0.5) continue;          
	  x=trR->X();
	  y=trR->Y();
	  z=trR->Z();
	  Double_t xx, zz;
	  Int_t mod;
	  segmentation2->GlobalToDet(il,x,y,x,xx,zz,mod);
	  //printf("cluster module %i  - module %i \n",recp->GetModule(),mod);
	  Double_t rTr = TMath::Sqrt(x*x+y*y); 

	  Double_t difx, dify, difz=0.;
	  Double_t dify=0.;
	  Double_t difz=0.;
	  Double_t deltaR=0.;
	  difx = xr - x;
	  dify = yr - y;
	  difz = zr - z;
	  Double_t phiGlob = TMath::ATan2(yr,xr);
	  if(yr<0) phiGlob+=TMath::TwoPi();
	  Double_t phiTr = TMath::ATan2(y,x);
	  if(y<0) phiTr+=TMath::TwoPi();
	  Double_t deltaPhi =  phiGlob-phiTr;
	      
	  Double_t deltaRphi = (phiGlob*TMath::Sqrt(xr*xr+yr*yr)-phiTr*rTr)*1e+04;
	  if(TMath::Abs(deltaRphi)/1e+04 > 3*(10000.*seg->GetCellSizeX(lay))) {
	    // printf("   Layer %i   type %i  dRphi %f  segm %f  \n",lay,recp->GetType(),deltaRphi/1e+04,seg->GetCellSizeX(lay)*10000.);
	    // printf("   Layer %i   type %i   Track (x,y,z) = (%f,%f,%f) Phi %f   Rphi %f         Cluster (x,y,z) = (%f,%f,%f)  Phi %f Rphi %f \n",lay,recp->GetType(),x,y,z,TMath::RadToDeg()*phiTr,phiTr*TMath::Sqrt(x*x+y*y),xr,yr,zr,TMath::RadToDeg()*phiGlob,phiGlob*TMath::Sqrt(xr*xr+yr*yr));
	  }              
	  deltaR = (TMath::Sqrt(xr*xr+yr*yr)-TMath::Sqrt(x*x+y*y))*1e+04;

	  hDiffx[il]->Fill(difx*1e04);
	  hDiffy[il]->Fill(dify*1e04);
	  hDiffz[il]->Fill(difz*1e04);
	  hDiffRphi[il]->Fill(deltaRphi);
	  hDiffR[il]->Fill(deltaR);
           
	}//loop entries fast 
      } 
    }//loop clusters
    //}//loop layer
  }//loop entries runloader
  
  TCanvas *c3x = new TCanvas("cGx","Delta in X global",200,10,900,900);
  c3x->Divide(3,nLayers/3);
  TCanvas *c3y = new TCanvas("cGy","Delta in Y global",200,10,900,900);
  c3y->Divide(3,nLayers/3);
  TCanvas *c3z = new TCanvas("cGz","Delta in Z global",200,10,900,900);
  c3z->Divide(3,nLayers/3);
  TCanvas *c3rphi = new TCanvas("cRphi","Delta R-Phi",200,10,900,900);
  c3rphi->Divide(3,nLayers/3);
  TCanvas *c3r = new TCanvas("cR","Delta R-Phi",200,10,900,900);
  c3r->Divide(3,nLayers/3);
  

  for(Int_t iL=0; iL< nLayers; iL++){
    c3x->cd(iL+1);
    hDiffx[iL]->Draw();
    c3y->cd(iL+1);
    hDiffy[iL]->Draw();
    c3z->cd(iL+1);
    hDiffz[iL]->Draw();
    c3rphi->cd(iL+1);
    hDiffRphi[iL]->Draw();
    c3r->cd(iL+1);
    hDiffR[iL]->Draw();
  }
 
}// main macro read


Int_t GetLayer(AliTrackReference *ref, TObject *obj){
 
  AliITSsegmentationUpgrade *seg = (AliITSsegmentationUpgrade*)obj;
  Int_t ilayer =-1;
  for(Int_t ila=0; ila<6; ila++){
    if(ila==(nLayers-1)){
      if(ref->R()>seg->GetRadius(ila)+seg->GetThickness(ila)) ilayer=(nLayers-1);
    }
    else {
      if(ref->R()>seg->GetRadius(ila)+seg->GetThickness(ila) && ref->R()<seg->GetRadius(ila+1)) ilayer=ila;
    }
  }
  return ilayer;
}

