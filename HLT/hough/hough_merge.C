void hough_merge(char *rootfile,Bool_t MC=false)
{
  gStyle->SetOptStat(0);
  
  Int_t xbin = 60;
  Int_t ybin = 60;
  Float_t xrange[2] = {-0.006 , 0.006}; //Pt 0.1->4GeV 0.006-0.006
  //Float_t yrange[2] = {-0.17 , 0.17}; //slice 2 0.55->0.88
  Float_t yrange[2] = {-0.26 , 0.26}; //slice 2 0.55->0.88
  //Float_t yrange[2] = {0.55 , 0.88}; //slice 2 0.55->0.88

  Int_t xr[2] = {0,250};
  Int_t yr[2] = {-125,125};
  TH1F *ntracks = new TH1F("ntracks","",100,0,200);
  TH2F *raw = new TH2F("raw","",250,xr[0],xr[1],250,yr[0],yr[1]);
  TH2F *road = new TH2F("road","",250,0,250,250,yr[0],yr[1]);
  TH2F *fake = new TH2F("fake","",300,0,300,250,-125,125);
  TH2F *peaks = new TH2F("peaks","",xbin,xrange[0],xrange[1],ybin,yrange[0],yrange[1]);
  peaks->SetMarkerStyle(3);
  peaks->SetMarkerColor(2);
  
  int slice = 2,patch=0,n_phi_segments=40;
  float eta[2] = {0.3,0.4};
  // float eta[2] = {0.04,0.05};
  
  AliL3HoughTransformer *a;
  AliL3HoughMaxFinder *b;
  AliL3HoughEval *c = new AliL3HoughEval(slice);
  AliL3HoughMerge *d = new AliL3HoughMerge(slice);
  TH2F *hist;
  for(int pat=0; pat<5; pat++)
    {
      a = new AliL3HoughTransformer(slice,pat,eta,n_phi_segments);
      hist = new TH2F("hist","Parameter space",xbin,xrange[0],xrange[1],ybin,yrange[0],yrange[1]);
      a->GetPixels(rootfile,raw);
      a->Transform2Circle(hist,87);
      
      b = new AliL3HoughMaxFinder("KappaPhi");
      AliL3TrackArray *tracks = b->FindMaxima(hist);
      c->SetTransformer(a);
      c->LookInsideRoad(tracks);
      d->FillTracks(tracks,pat);
      delete a;
      delete b;
    }
  return;
  

  TH2F *merge_hist = new TH2F("merge_hist","Parameter space",xbin,xrange[0],xrange[1],ybin,yrange[0],yrange[1]);
  
  d->FillHisto(merge_hist);
  
  
  b = new AliL3HoughMaxFinder(slice,pat);
  b->FindMaxima(merge_hist);
  
  AliL3TrackArray *mtrs = b->GetTracks();
  
  c->CompareMC(rootfile,mtrs,eta);
  
  AliL3Transform *transform = new AliL3Transform();
  for(int tr=0; tr<mtrs->GetNTracks(); tr++)
    {
      AliL3HoughTrack *t = (AliL3HoughTrack*)mtrs->GetCheckedTrack(tr);
      if(!t) {printf("NO TRACK11\n"); continue;}
      if(t->GetMCid() < 0) {printf("Track %d was fake, weigth %d\n",tr,t->GetNHits()); continue;}
      printf("Found pt %f weigth %d\n",t->GetPt(),t->GetNHits());
      //printf("charge %f\n",t->GetCharge());
      for(Int_t j=0; j<174; j++)
	{
	  Float_t xyz[3];
	  if(!t->GetCrossingPoint(j,xyz)) 
	    {
	      printf("Track does not cross line\n");
	      continue;
	    }
	  road->Fill(xyz[0],xyz[1],1);
	  
	}
	
      
    }
  

  TCanvas *c1 = new TCanvas("c1","",800,800);
  SetCanvasOptions(c1);
  c1->Divide(2,2);
  c1->cd(1);
  merge_hist->Draw("box");
  
  //TCanvas *c2 = new TCanvas("c2","",2);
  c1->cd(3);
  raw->Draw();

  //TCanvas *c3 = new TCanvas("c3","",2);
  c1->cd(4);
  road->SetMarkerStyle(7);
  road->Draw();

  delete b;
  delete d;
  delete c;
  //----------------------------------------------------------------------------
  if(MC)
    {
      TFile *file = new TFile(rootfile);
      file->cd();
      
      gAlice = (AliRun*)file->Get("gAlice");
      gAlice->GetEvent(0);
      
      TClonesArray *particles=gAlice->Particles();
      Int_t n=particles->GetEntriesFast();
      Float_t torad=TMath::Pi()/180;
      Float_t phi_min = slice*20 - 10;
      Float_t phi_max = slice*20 + 10;
      
      for (Int_t j=0; j<n; j++) {
	TParticle *p=(TParticle*)particles->UncheckedAt(j);
	if (p->GetFirstMother()>=0) continue;  //secondary particle
	if(p->Eta() < eta[0] || p->Eta() > eta[1]) continue;
	Double_t ptg=p->Pt(),pxg=p->Px(),pyg=p->Py(),pzg=p->Pz();
	Double_t phi_part = TMath::ATan2(pyg,pxg);
	if (phi_part < 0) phi_part += 2*TMath::Pi();
	
	if(phi_part < phi_min*torad || phi_part > phi_max*torad) {continue;}
	if(ptg<0.100) continue;
	int found = 0;
	for(int m=0; m<mtrs->GetNTracks(); m++)
	  {
	    AliL3HoughTrack *tra = (AliL3HoughTrack*)mtrs->GetCheckedTrack(m);
	    if(tra->GetMCid() != j) continue;
	    found = 1;
	    double dPt = fabs(ptg-tra->GetPt());
	    float  phi0[2] = {tra->GetPhi0(),0};
	    transform->Local2GlobalAngle(phi0,slice);
	    double dPhi0 = fabs(phi_part-phi0[0]);
	    printf("Difference Pt %f Phi0 %f\n",dPt,dPhi0);

	  }
	if(!found) printf("Track %d was not found\n",j);
	printf("Generated particle %d, with pt %f phi %f\n",p->GetPdgCode(),ptg,phi_part);
      }
      file->Close();
      delete file;
    }
  
}
