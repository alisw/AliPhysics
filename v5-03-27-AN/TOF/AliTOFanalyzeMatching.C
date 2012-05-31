void AliTOFanalyzeMatching(const char* datafile)
{

  //
  // Matching efficiency and contamination
  // for different particle species
  // (pions, kaons and protons).
  // All histos are saved into a separate file
  // datafile is assumed to be the file name containing
  // the results of the matching in TNtuple format.
  //
  // Author: F. Pierella | pierella@bo.infn.it
  //
  // Use case:
  // start root
  // root[0] .L AliTOFanalyzeMatching.C
  // root[1] AliTOFanalyzeMatching("matchingNtuple.root")

  // output (histos!) filename
  char outFileName[100];
  strcpy(outFileName,"histo");
  strcat(outFileName,datafile);
  
  // dummy histos (for normalization)
  TH1F* hpitot= new TH1F("hpitot","",12,0.,3.);
  TH1F* hkatot= new TH1F("hkatot","",12,0.,3.);
  TH1F* hprtot= new TH1F("hprtot","",12,0.,3.);
  
  TH1F* hpimatched= new TH1F("hpimatched","",12,0.,3.);
  TH1F* hkamatched= new TH1F("hkamatched","",12,0.,3.);
  TH1F* hprmatched= new TH1F("hprmatched","",12,0.,3.);
  
  
  // matching efficiency histos
  TH1F* hpimatcheff= new TH1F("hpimatcheff","Matching efficiency for pions",12,0.,3.);
  TH1F* hkamatcheff= new TH1F("hkamatcheff","Matching efficiency for kaons",12,0.,3.);
  TH1F* hprmatcheff= new TH1F("hprmatcheff","Matching efficiency for protons",12,0.,3.);

  // matching contamination histos
  TH1F* hpimatchcon= new TH1F("hpimatchcon","Matching contamination for pions",12,0.,3.);
  TH1F* hkamatchcon= new TH1F("hkamatchcon","Matching contamination for kaons",12,0.,3.);
  TH1F* hprmatchcon= new TH1F("hprmatchcon","Matching contamination for protons",12,0.,3.);
  
  
  TFile *file = TFile::Open(datafile,"old");
  TNtuple* fNtuple= (TNtuple*)file->Get("Ntuple"); // get ntuple from file
  Int_t nvar = fNtuple->GetNvar(); cout <<"N of var.="<< nvar << endl;
  fNtuple->GetEvent(0);
  
  file->cd();
  Int_t nparticles = (Int_t)fNtuple->GetEntries();
  
  for (Int_t i=0; i < nparticles; i++) {
    fNtuple->GetEvent(i);
    Int_t event=fNtuple->GetLeaf("event")->GetValue();   
    Int_t pdgcode=fNtuple->GetLeaf("ipart")->GetValue(); 
    Int_t matc=fNtuple->GetLeaf("matc")->GetValue(0);
    Float_t px=fNtuple->GetLeaf("pxvtx")->GetValue(0);
    Float_t py=fNtuple->GetLeaf("pyvtx")->GetValue(0);
    Float_t pz=fNtuple->GetLeaf("pzvtx")->GetValue(0);
    
    Float_t pvtx=TMath::Sqrt(px*px+py*py+pz*pz);  
    Float_t ptvtx=TMath::Sqrt(px*px+py*py);
    Int_t abspdgcode=TMath::Abs(pdgcode);
    

    // N (1+2+3+4+(-4)) cases
    if(matc>=1 || matc==-4){
      switch(abspdgcode){
      case 211:
	hpitot->Fill(pvtx);
	break;
      case 321:
	hkatot->Fill(pvtx);
	break;
      case 2212:
	hprtot->Fill(pvtx);
	break;
      }
    }


    // N_matched (3+4) cases
    if(matc==3 || matc==4){
      switch(abspdgcode){
      case 211:
	hpimatched->Fill(pvtx);
	break;
      case 321:
	hkamatched->Fill(pvtx);
	break;
      case 2212:
	hprmatched->Fill(pvtx);
	break;
      }
    }


    // N_t (3) case
    if(matc==3){
      switch(abspdgcode){
      case 211:
	hpimatcheff->Fill(pvtx);
	break;
      case 321:
	hkamatcheff->Fill(pvtx);
	break;
      case 2212:
	hprmatcheff->Fill(pvtx);
	break;
      }
    }

    // N_w (4) case
    if(matc==4){
      switch(abspdgcode){
      case 211:
	hpimatchcon->Fill(pvtx);
	break;
      case 321:
	hkamatchcon->Fill(pvtx);
	break;
      case 2212:
	hprmatchcon->Fill(pvtx);
	break;
      }
    }

  }

  // histo normalization
  // efficiency
  hpimatcheff->Divide(hpitot);
  hkamatcheff->Divide(hkatot);
  hprmatcheff->Divide(hprtot);

  // contamination
  hpimatchcon->Divide(hpimatched);
  hkamatchcon->Divide(hkamatched);
  hprmatchcon->Divide(hprmatched);


  TFile *houtfile = new TFile(outFileName,"recreate");
  houtfile->cd();

  hpitot->Write();
  hkatot->Write();
  hprtot->Write();
    
  hpimatched->Write();
  hkamatched->Write();
  hprmatched->Write();
    
  hpimatcheff->Write();
  hkamatcheff->Write();
  hprmatcheff->Write();

  hpimatchcon->Write();
  hkamatchcon->Write();
  hprmatchcon->Write();
  houtfile->Close();   

  cout << "File " << outFileName << " with histos has been created" << endl;
}
