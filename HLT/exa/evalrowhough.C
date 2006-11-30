// $Id$

struct GoodTrack 
{
  Int_t label;
  Double_t eta;
  Int_t code;
  Double_t px,py,pz;
  Double_t x,y,z;
  Int_t nhits;
  Int_t sector;
};
typedef struct GoodTrack GoodTrack;

//Histograms
TNtuple *fNtuppel;
TH1F *fPtRes;
TH1F *fNGoodTracksPt;
TH1F *fNFoundTracksPt;
TH1F *fNFakeTracksPt;
TH1F *fTrackEffPt;
TH1F *fFakeTrackEffPt;
TH1F *fNGoodTracksEta;
TH1F *fNFoundTracksEta;
TH1F *fNFakeTracksEta;
TH1F *fTrackEffEta;
TH1F *fFakeTrackEffEta;
TH1F *fNGoodTracksPhi;
TH1F *fNFoundTracksPhi;
TH1F *fNFakeTracksPhi;
TH1F *fTrackEffPhi;
TH1F *fFakeTrackEffPhi;
TProfile *fFakesPt;
TProfile *fFakesPhi;
TProfile *fSlicesPt;
TProfile *fSlicesPhi;
TProfile2D *fFakesPtvsPhi;
TNtuple *fNtupleRes;
TH1F *fPtRes2;
TH1F *fEtaRes;
TH1F *fPsiRes;

TH1F *fNGoodTracksRad;
TH1F *fNFoundTracksRad;
TH1F *fNGoodTracksZ;
TH1F *fNFoundTracksZ;
TH1F *fRadEff;
TH1F *fZEff;

void evalrowhough(Char_t *path="./fitter",Char_t *offlinepath="./",Int_t ievent=0, Float_t ptmin=0.4)  
{
  CreateHistos(25,0.1,4.1);

  Char_t fname[1024];
  sprintf(fname,"%s/tracks_%d.raw",path,ievent);
  AliHLTFileHandler *tfile = new AliHLTFileHandler();
  if(!tfile->SetBinaryInput(fname)){
    LOG(AliHLTLog::kError,"AliHLTEvaluation::Setup","File Open")
      <<"Inputfile "<<fname<<" does not exist"<<ENDLOG; 
    return;
  }

  AliHLTTrackArray *fTracks = new AliHLTTrackArray("AliHLTHoughTrack");
  tfile->Binary2TrackArray(fTracks);
  //fTracks->QSort();
  tfile->CloseBinaryInput();
  delete tfile;

  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      AliHLTHoughTrack *track = fTracks->GetCheckedTrack(i);
      if(!track) continue; 
      
      track->SetEta(fTracks->GetCheckedTrack(i)->GetPseudoRapidity());
      UInt_t *ids = track->GetHitNumbers();
      Int_t slice = (ids[0]>>25)&0x7f;
      track->SetSlice(slice);
      cout<<"Track with label "<<ids<<" "<<track->GetNHits()<<" "<<track->GetMCid()<<" "<<ids[0]<<" "<<track->GetWeight()<<" "<<track->GetSlice()<<endl;
    }

  Char_t filename[1024];
  sprintf(filename,"%s/good_tracks_tpc_%d",offlinepath,ievent);
  ifstream in(filename);
  if(!in)
    {
      cerr<<"AliHLTEvaluate::GetGoodParticles : Problems opening file :"<<filename<<endl;
      return;
    }

  Int_t MaxTracks=20000;
  Int_t fGoodGen=0;
  fGoodTracks = new GoodTrack[MaxTracks];
  
  while (in>>fGoodTracks[fGoodGen].label>>fGoodTracks[fGoodGen].code>>
	 fGoodTracks[fGoodGen].px>>fGoodTracks[fGoodGen].py>>fGoodTracks[fGoodGen].pz>>
	 fGoodTracks[fGoodGen].x>>fGoodTracks[fGoodGen].y>>fGoodTracks[fGoodGen].z)
    {
      fGoodTracks[fGoodGen].nhits=-1;
      fGoodTracks[fGoodGen].sector=-1; 
      fGoodGen++;
      if (fGoodGen==MaxTracks) 
	{
	  cerr<<"AliHLTEvaluate::GetGoodParticles : Too many good tracks !\n";
	  break;
	}
    }

  Int_t fGoodGen2 = 0;
  Int_t fGoodFound = 0;
  Int_t fGoodFound2 = 0;
  Float_t fMinGoodPt = ptmin;
  Float_t fMaxGoodPt = 99999.;

  if(!fGoodTracks)
    {
      cerr<<"AliHLTEvaluate::FillEffHistos : No good tracks"<<endl;
      return;
    }
  cout<<"Comparing "<<fGoodGen<<" good tracks ..."<<endl;
  for(Int_t i=0; i<fGoodGen; i++)
    {
      //cout<<"Checking particle "<<i<<endl;
      Float_t ptg = TMath::Sqrt(fGoodTracks[i].px*fGoodTracks[i].px + fGoodTracks[i].py*fGoodTracks[i].py);
      Float_t pg = TMath::Sqrt(fGoodTracks[i].px*fGoodTracks[i].px + fGoodTracks[i].py*fGoodTracks[i].py + fGoodTracks[i].pz*fGoodTracks[i].pz);
      if(ptg < fMinGoodPt || ptg > fMaxGoodPt) continue;
      Float_t pzg=fGoodTracks[i].pz;
      Float_t dipangle=TMath::ATan2(pzg,ptg)*180./TMath::Pi();
      Float_t etag=0.5 * TMath::Log((pg+pzg)/(pg-pzg));
      Float_t phig=TMath::ATan2(fGoodTracks[i].py,fGoodTracks[i].px);
      Float_t rg = TMath::Sqrt(fGoodTracks[i].x*fGoodTracks[i].x + fGoodTracks[i].y*fGoodTracks[i].y);
      Float_t zg = TMath::Abs(fGoodTracks[i].z);
      if (phig<0) phig+=2*TMath::Pi();

      fNGoodTracksPt->Fill(ptg);
      fNGoodTracksEta->Fill(etag);
      fNGoodTracksPhi->Fill(phig);
      if((zg != 0) || (rg != 0)) { 
	fNGoodTracksRad->Fill(rg);
	fNGoodTracksZ->Fill(zg);
      }
      fGoodGen2++;
      Int_t found = 0;
      Int_t nslices = 1;
      Int_t org_slice;
      for(Int_t k=0; k<fTracks->GetNTracks(); k++)
	{
	  AliHLTHoughTrack *track = (AliHLTHoughTrack *)fTracks->GetCheckedTrack(k);
	  if(!track) continue;
	  //	  Int_t nHits = track->GetNumberOfPoints();
	  //	  if(nHits < fMinPointsOnTrack) break;
	  Int_t nHits = 1;
	  Int_t tracklabel;
	  tracklabel = track->GetMCid();
	  
	  if(TMath::Abs(tracklabel) != fGoodTracks[i].label) continue;
	  if(found == 0) {
	    found=1;
	    Float_t pt=track->GetPt();
	    org_slice = track->GetSlice();
	    Float_t eta = track->GetEta();
	    Float_t phi = track->GetPsi();
	    //cout<<eta<<" "<<phi<<" "<<etag<<" "<<phig<<" "<<track->GetEtaIndex()<<endl;
	    if(tracklabel == fGoodTracks[i].label) 
	      {
		fGoodFound++;
		
		fNFoundTracksPt->Fill(ptg); 
		fNFoundTracksEta->Fill(etag);
		fNFoundTracksPhi->Fill(phig);
		if((zg != 0) || (rg != 0)) { 
		  fNFoundTracksRad->Fill(rg);
		  fNFoundTracksZ->Fill(zg);
		}
		fNtuppel->Fill(ptg,pt,nHits);
		fPtRes->Fill((pt-ptg)/ptg*100.);
		fEtaRes->Fill(eta-etag);
		fPsiRes->Fill(phi-phig);
	      }
	    else 
	      {
		fNFakeTracksPt->Fill(ptg); 
		fNFakeTracksEta->Fill(etag);
		fNFakeTracksPhi->Fill(phig);
	      }
	    //fPtRes->Fill((pt-ptg)/ptg*100.);
	    //fNtuppel->Fill(ptg,pt,nHits);
	  }
	  else {
	    found++;
	    if(track->GetSlice() != org_slice)
	      nslices = 2;
	  }
	}
      if(!found)
	//cout<<"Track "<<fGoodTracks[i].label<<" was not found"<<endl;
	continue;
      else {
	fFakesPt->Fill(ptg,(Double_t)found-1.0);
	fFakesPhi->Fill(phig,(Double_t)found-1.0);
	fSlicesPt->Fill(ptg,(Double_t)nslices);
	fSlicesPhi->Fill(phig,(Double_t)nslices);
	//Float_t local_phi = phi-(phi/)
	//fFakesPtvsPhi->Fill(1.0/ptg,
      }
      fGoodFound2 += found;
    }

  cout<<fGoodFound<<"("<<fGoodFound2<<") tracks found out of "<<fGoodGen2<<" generated tracks"<<endl;
  cout<<" The overall efficiency is "<<(Float_t)fGoodFound/(Float_t)fGoodGen2<<endl;
  
  CalcEffHistos();
  plotptres(ptmin);
  Write2File();
}

void CreateHistos(Int_t nbin,Float_t xlow,Float_t xup)
{
  //Create the histograms 
  
  fNtuppel = new TNtuple("fNtuppel","Pt resolution","pt_gen:pt_found:nHits");
  fNtuppel->SetDirectory(0);
  fPtRes = new TH1F("fPtRes","Relative Pt resolution",30,-10.,10.); 
  fNGoodTracksPt = new TH1F("fNGoodTracksPt","Good tracks vs pt",nbin,xlow,xup);    
  fNFoundTracksPt = new TH1F("fNFoundTracksPt","Found tracks vs pt",nbin,xlow,xup);
  fNFakeTracksPt = new TH1F("fNFakeTracksPt","Fake tracks vs pt",nbin,xlow,xup);
  fTrackEffPt = new TH1F("fTrackEffPt","Tracking efficiency vs pt",nbin,xlow,xup);
  fFakeTrackEffPt = new TH1F("fFakeTrackEffPt","Efficiency for fake tracks vs pt",nbin,xlow,xup);
  
  fNGoodTracksEta = new TH1F("fNGoodTracksEta","Good tracks vs eta",nbin,-1,1);
  fNFoundTracksEta = new TH1F("fNFoundTracksEta","Found tracks vs eta",nbin,-1,1);
  fNFakeTracksEta = new TH1F("fNFakeTracksEta","Fake tracks vs eta",nbin,-1,1);
  fTrackEffEta = new TH1F("fTrackEffEta","Tracking efficienct vs eta",nbin,-1,1);
  fFakeTrackEffEta = new TH1F("fFakeTrackEffEta","Efficiency for fake tracks vs eta",nbin,-1,1);

  
  fNGoodTracksPhi = new TH1F("fNGoodTracksPhi","Good tracks vs phi",nbin,0,6.28);
  fNFoundTracksPhi = new TH1F("fNFoundTracksPhi","Found tracks vs phi",nbin,0,6.28);
  fNFakeTracksPhi = new TH1F("fNFakeTracksPhi","Fake tracks vs phi",nbin,0,6.28);
  fTrackEffPhi = new TH1F("fTrackEffPhi","Tracking efficienct vs phi",nbin,0,6.28);
  fFakeTrackEffPhi = new TH1F("fFakeTrackEffPhi","Efficiency for fake tracks vs phi",nbin,0,6.28);
  fFakesPt = new TProfile("fFakesPt","Number of ghosts vs pt",nbin,xlow,xup,0,100);
  fFakesPhi = new TProfile("fFakesPhi","Number of ghosts vs phi",nbin,0,6.28,0,100);
  fSlicesPt = new TProfile("fSlicesPt","Number of slices vs pt",nbin,xlow,xup,0,100);
  fSlicesPhi = new TProfile("fSlicesPhi","Number of slices vs phi",nbin,0,6.28,0,100);
  fPtRes2 = new TH1F("fPtRes2","Relative Pt resolution (Pt>2 GeV)",20,-100.,100.);
  fEtaRes = new TH1F("fEtaRes","Eta resolution",30,-0.1,0.1);
  fPsiRes = new TH1F("fPsiRes","Psi resolution",30,-0.1,0.1);

  fNGoodTracksRad = new TH1F("fNGoodTracksRad","Good tracks vs distance(x,y)",1,0,0.2);    
  fNFoundTracksRad = new TH1F("fNFoundTracksRad","Found tracks vs distance(x,y)",1,0,0.2);
  fRadEff = new TH1F("fRadEff","Tracking efficiency vs distance(x,y)",1,0,0.2);
  fNGoodTracksZ = new TH1F("fNGoodTracksZ","Good tracks vs distance(z)",1,0,0.2);    
  fNFoundTracksZ = new TH1F("fNFoundTracksZ","Found tracks vs distance(z)",1,0,0.2);
  fZEff = new TH1F("fZEff","Tracking efficiency vs distance(z)",1,0,0.2);
}

void Write2File()
{
  //Write histograms to file:
  
  TFile *of = TFile::Open("hough_efficiencies.root","RECREATE");
  if(!of->IsOpen())
    {
      cout<<"Problems opening rootfile"<<endl;
      return;
    }
  
  of->cd();
  fNtuppel->Write();
  fPtRes->Write();
  fNGoodTracksPt->Write();
  fNFoundTracksPt->Write();
  fNFakeTracksPt->Write();
  fTrackEffPt->Write();
  fFakeTrackEffPt->Write();
  fNGoodTracksEta->Write();
  fNFoundTracksEta->Write();
  fNFakeTracksEta->Write();
  fTrackEffEta->Write();
  fFakeTrackEffEta->Write();

  fNGoodTracksPhi->Write();
  fNFoundTracksPhi->Write();
  fNFakeTracksPhi->Write();
  fTrackEffPhi->Write();
  fFakeTrackEffPhi->Write();
  
  fFakesPt->Write();
  fFakesPhi->Write();
  fSlicesPt->Write();
  fSlicesPhi->Write();

  fPtRes2->Write();
  fEtaRes->Write();
  fPsiRes->Write();

  fNGoodTracksRad->Write();
  fNFoundTracksRad->Write(); 
  fRadEff->Write();
  fNGoodTracksZ->Write();
  fNFoundTracksZ->Write();
  fZEff->Write();

  of->Close();
}

void CalcEffHistos()
{  
  fNFoundTracksPt->Sumw2(); fNGoodTracksPt->Sumw2();
  fTrackEffPt->Divide(fNFoundTracksPt,fNGoodTracksPt,1,1.,"b");
  fFakeTrackEffPt->Divide(fNFakeTracksPt,fNGoodTracksPt,1,1.,"b");
  fTrackEffPt->SetMaximum(1.4);
  fTrackEffPt->SetXTitle("P_{T} [GeV]");
  fTrackEffPt->SetLineWidth(2);
  fFakeTrackEffPt->SetFillStyle(3013);
  fTrackEffPt->SetLineColor(4);
  fFakeTrackEffPt->SetFillColor(2);

  fNFoundTracksEta->Sumw2(); fNGoodTracksEta->Sumw2();
  fTrackEffEta->Divide(fNFoundTracksEta,fNGoodTracksEta,1,1.,"b");
  fFakeTrackEffEta->Divide(fNFakeTracksEta,fNGoodTracksEta,1,1.,"b");
  fTrackEffEta->SetMaximum(1.4);
  fTrackEffEta->SetXTitle("#eta");
  fTrackEffEta->SetLineWidth(2);
  fFakeTrackEffEta->SetFillStyle(3013);
  fTrackEffEta->SetLineColor(4);
  fFakeTrackEffEta->SetFillColor(2);
       
  fNFoundTracksPhi->Sumw2(); fNGoodTracksPhi->Sumw2();
  fTrackEffPhi->Divide(fNFoundTracksPhi,fNGoodTracksPhi,1,1.,"b");
  fFakeTrackEffPhi->Divide(fNFakeTracksPhi,fNGoodTracksPhi,1,1.,"b");
  fTrackEffPhi->SetMaximum(1.4);
  fTrackEffPhi->SetXTitle("#psi [rad]");
  fTrackEffPhi->SetLineWidth(2);
  fFakeTrackEffPhi->SetFillStyle(3013);
  fTrackEffPhi->SetLineColor(4);
  fFakeTrackEffPhi->SetFillColor(2);

  fEtaRes->SetXTitle("#eta");
  fPsiRes->SetXTitle("#psi [rad]");

  fRadEff->Divide(fNFoundTracksRad,fNGoodTracksRad,1,1.,"b");
  fRadEff->SetXTitle("Distance [cm]");
  fZEff->Divide(fNFoundTracksZ,fNGoodTracksZ,1,1.,"b");
  fZEff->SetXTitle("Distance [cm]");
}

void plotptres(Float_t ptmin)
{
  //Plot the relative pt resolution vs pt.
  
  const Int_t n=8;
  
  Double_t hltx[9] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,2.0};
  Double_t hltxerr[n];
  Double_t hltavx[n];
  Double_t hlty[n];
  Double_t hltyerr[n];
  Char_t string[1024];

  for(Int_t i=0; i<n; i++)
    {
      TH1F *hist = new TH1F("hist","",100,-50.*(float)(i+1),50.*(float)(i+1));
      sprintf(string,"pt_gen > %f && pt_gen <= %f",hltx[i],hltx[i+1]);
      fNtuppel->Draw("(pt_found-pt_gen)/pt_gen*100>>hist",string,"goff");
      TF1 *f1 = new TF1("f1","gaus");
      hist->Fit("f1","E");
      hlty[i] = f1->GetParameter(2);
      hltyerr[i] = f1->GetParError(2);
      //hlty[i] = hist->GetRMS();
      //hltyerr[i] = 0;
      hltavx[i] = (hltx[i]+hltx[i+1])/2.0;
      hltxerr[i] = (hltx[i+1]-hltx[i])/2.0;
      cout<<string<<" "<<i<<" "<<hlty[i]<<" "<<hltyerr[i]<<endl;
      delete f1;
      delete hist;
    }

  sprintf(string,"pt_gen > %f",2.0);
  fNtuppel->Draw("(pt_found-pt_gen)/pt_gen*100>>fPtRes2",string,"goff");
  TF1 *f1 = new TF1("f1","gaus");
  fPtRes2->Fit("f1","E");
  fPtRes2->SetXTitle("#Delta P_{t} / P_{t} [%]");
  fPtRes2->SetYTitle("#Delta N");
  delete f1;

  TGraphErrors *gr1 = new TGraphErrors(n,hltavx,hlty,hltxerr,hltyerr);
  TCanvas *c1 = new TCanvas("c1","",1);
  fPtRes = c1->DrawFrame(0.1,0,2,15);
  gr1->SetTitle("");
  gr1->Draw("P");
  gr1->SetMarkerStyle(20);
  gr1->SetLineWidth(2);
  gr1->SetMarkerSize(1.3);
  f1 = new TF1("f1","pol1",ptmin,2.0);
  gr1->Fit("f1","ER");
  fPtRes->SetXTitle("p_{t} [GeV]");
  fPtRes->SetYTitle("#Delta P_{t} / P_{t} [%]");
  delete f1;
}
