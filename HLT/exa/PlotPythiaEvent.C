void PlotPythiaEvent(char* path="/data1/AliRoot/pp/pythia/14000")
{
  //gSystem->Load("$ROOTSYS/lib/libEG.so");                   // Root Event Generator interface

  TH1F *hVr = new TH1F("hVr","Secondary vertex distribution (r direction)",100,0,100);
  TH1F *hVr2 = new TH1F("hVr2","Secondary vertex distribution near primary vertex",100,0,50);
  
  TH1F *hVz = new TH1F("hVz","Secondary vertex distribution",100,-100,100);
  TH1F *hVz10 = new TH1F("hVz10","Secondary vertex distribution",10,-10,10);
  TH1F *hVz20 = new TH1F("hVz20","Secondary vertex distribution",200,-20,20);
  TH1F *hVz30 = new TH1F("hVz30","Secondary vertex distribution",300,-30,30);
  TH1F *hVz40 = new TH1F("hVz40","Secondary vertex distribution",400,-40,40);
  TH1F *hVz50 = new TH1F("hVz50","Secondary vertex distribution",500,-50,50);
  TH1F *hVz60 = new TH1F("hVz60","Secondary vertex distribution",600,-60,60);
  TH1F *hVz70 = new TH1F("hVz70","Secondary vertex distribution",700,-70,70);
  TH1F *hVz80 = new TH1F("hVz80","Secondary vertex distribution",800,-80,80);
  TH1F *hVz90 = new TH1F("hVz90","Secondary vertex distribution",900,-90,90);
  TH1F *hVz100 = new TH1F("hVz100","Secondary vertex distribution",100,-100,100);

  TH1F *hPt = new TH1F("hPt","Transversal momentum",50,0,10);
  TH1F *hPt10 = new TH1F("hPt10","Transversal momentum for particles within the vertex cut",50,0,10);
  TH1F *hPt20 = new TH1F("hPt20","Transversal momentum for particles within the vertex cut",50,0,10);
  TH1F *hPt30 = new TH1F("hPt30","Transversal momentum for particles within the vertex cut",50,0,10);
  TH1F *hPt40 = new TH1F("hPt40","Transversal momentum for particles within the vertex cut",50,0,10);
  TH1F *hPt50 = new TH1F("hPt50","Transversal momentum for particles within the vertex cut",50,0,10);
  TH1F *hPt60 = new TH1F("hPt60","Transversal momentum for particles within the vertex cut",50,0,10);
  TH1F *hPt70 = new TH1F("hPt70","Transversal momentum for particles within the vertex cut",50,0,10);
  TH1F *hPt80 = new TH1F("hPt80","Transversal momentum for particles within the vertex cut",50,0,10);
  TH1F *hPt90 = new TH1F("hPt90","Transversal momentum for particles within the vertex cut",50,0,10);
  TH1F *hPt100 = new TH1F("hPt100","Transversal momentum for particles within the vertex cut",50,0,10);

/*TH1F *hPt3 = new TH1F("hPt3","Transversal momentum for particles outside the vertex cut",100,0,10);
  hPt3->SetXTitle("GeV/c");
  hPt3->SetYTitle("N");

  TH1F *hPt4 = new TH1F("hPt4","hPt2/hPt",100,0,10);
  hPt4->SetXTitle("GeV/c");
  hPt4->SetYTitle("N");

  TH1F *hPt5 = new TH1F("hPt5","hPt3/hPt",100,0,10);
  hPt5->SetXTitle("GeV/c");
  hPt5->SetYTitle("N");*/

  Char_t fname[256];
  Int_t nFiles = 4;
  for (Int_t iFile = 0; iFile < nFiles; iFile++) {

    sprintf(fname,"%s/pythia6_%d.root",path,iFile);
    TFile *f = new TFile(fname);
    cout << "Opened file " << fname << endl;

    TClonesArray* arr = new TClonesArray("TParticle");

    TTree *tree = (TTree *) f->Get("pythia");
    tree->GetBranch("particles")->SetAddress(&arr);
    Int_t nEntries = tree->GetEntries();
    
    TParticle *p;
    TParticle d1;
    //TParticle d2;

    for (Int_t i = 0; i < nEntries; i++) {
      tree->GetEvent(i);
      TClonesArray &part = *arr;
      Int_t nPart = part.GetEntriesFast();
      for (Int_t iPart = 0; iPart < nPart; iPart++) {
	p = (TParticle*)part[iPart];
	if (p->GetPdgCode()==3122){
	  
          Int_t mother = p->GetFirstMother();

	  Int_t fdaughter = p->GetFirstDaughter();
	  Int_t ldaughter = p->GetLastDaughter();
          d1 = (TParticle*)part[fdaughter-1];
          //if (!d1) continue;
          d2 = (TParticle*)part[ldaughter-1];
	  //cout << fdaughter << " " << ldaughter << endl;
	  //cout << d1->GetPdgCode() << " " << d2->GetPdgCode() << endl;

	  Double_t vx = d2->Vx();
	  Double_t vy = d2->Vy();
	  Double_t vz = d2->Vz();
	  Double_t vt = TMath::Sqrt((vx*vx)+(vy*vy));
	  Double_t vr = TMath::Sqrt((vx*vx)+(vy*vy)+(vz*vz));
	  Double_t pt = p->Pt();

	  hVr->Fill(vr);
	  if (vr < 10) hVr2->Fill(vr);
	  hVz->Fill(vz);
	  hPt->Fill(pt);
	  if (TMath::Abs(vz) < 10) 
	    {
	      hVz10->Fill(vz);
	      hPt10->Fill(pt);
	    }
	  if (TMath::Abs(vz) < 20)
	    {
	      hVz20->Fill(vz);
	      hPt20->Fill(pt);
	    }
	  if (TMath::Abs(vz) < 30)
	    {
	      hVz30->Fill(vz);
	      hPt30->Fill(pt);
	    }
	  if (TMath::Abs(vz) < 40)
	    {
	      hVz40->Fill(vz);
	      hPt40->Fill(pt);
	    }
	  if (TMath::Abs(vz) < 50)
	    {
	      hVz50->Fill(vz);
	      hPt50->Fill(pt);
	    }
	  if (TMath::Abs(vz) < 60)
	    {
	      hVz60->Fill(vz);
	      hPt60->Fill(pt);
	    }
	  if (TMath::Abs(vz) < 70)
	    {
	      hVz70->Fill(vz);
	      hPt70->Fill(pt);
	    }
	  if (TMath::Abs(vz) < 80)
	    {
	      hVz80->Fill(vz);
	      hPt80->Fill(pt);
	    }
	  if (TMath::Abs(vz) < 90)
	    {
	      hVz90->Fill(vz);
	      hPt90->Fill(pt);
	    }
	  if (TMath::Abs(vz) < 100)
	    {
	      hVz100->Fill(vz);
	      hPt100->Fill(pt);
	    }

	}
	
      }
      
    }
    
  }
  
  //TCanvas *c = new TCanvas("c","c",50,50,600,800);
  //c->Divide(1,2);
  //c->cd(1);
  //gPad->SetLogy();
  //hVz->Draw();
  //hPt->Draw();
  //c->cd(2);
  //gPad->SetLogy();
  //hVz2->Draw();
  //hPt3->Draw();

  //hPt4->Divide(hPt2,hPt);
  //hPt5->Divide(hPt3,hPt);
  //TCanvas *c2 = new TCanvas("c2","c2",50,50,600,800);
  //c2->Divide(1,2);
  //c2->cd(1);
  //hPt4->Draw();
  //c2->cd(2);
  //hPt5->Draw();

  TFile *outfile = new TFile("/tmp/results.root","recreate");
  hVr->Write();
  hVr2->Write();

  hVz->Write();
  hVz10->Write();
  hVz20->Write();
  hVz30->Write();
  hVz40->Write();
  hVz50->Write();
  hVz60->Write();
  hVz70->Write();
  hVz80->Write();
  hVz90->Write();
  hVz100->Write();

  hPt->Write(); 
  hPt10->Write();
  hPt20->Write();
  hPt30->Write();
  hPt40->Write();
  hPt50->Write();
  hPt60->Write();
  hPt70->Write();
  hPt80->Write();
  hPt90->Write();
  hPt100->Write();

  outfile->Close();
  f->Close();

}

