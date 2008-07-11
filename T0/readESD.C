void readESD() 
{
  //read T0 ESD from reconstructed cosmic data 

  char filename[100];
  sprintf(filename,"");
  TFile *file = new TFile("AliESDs.root");
  
  //  TH1F *hMean= new TH1F("hMean"," T0 ", 10000,1.249e6,1.259e6);
  TH1F *hMean= new TH1F("hMean"," T0 ",6000,0,6000);
  TH1F *hVertex = new TH1F("hVertex","Z position of vertex",   100,-120,-110);
  
  //  TH1F *hAmp = new TH1F("hAmp"," amplitude",   2500, 10,2510);
  //  TH1F *hTime = new TH1F("hTime"," time",   500,1250000,1250500);
  TH1F *hAmp[24];  TH1F * hTime[24]; 
  Char_t  buf1[20], buf2[20];
  
  for(Int_t ic=0; ic<24; ic++) 
    {
      sprintf(buf1,"Amp%i",ic+1);
      hAmp[ic] = new TH1F(buf1,"Amp",50, 0,50 );
      sprintf(buf2,"Time%i",ic+1);
      //      hTime[ic] = new TH1F("hTime"," time",   10000,1.249e6,1.259e6);
          hTime[ic] = new TH1F("hTime"," time",   6000,0,6000);
    }
  
  Double_t *amp, *time;
  TTree *esdTree =  (TTree*)file->Get("esdTree");
  esdTree->ls();

  TBranch *brRec=esdTree->GetBranch("AliESDTZERO.");
  AliESDTZERO *fesdT0 = new AliESDTZERO();
  if (brRec) {
    brRec->SetAddress(&fesdT0);
  }else{
    cerr<<"EXEC Branch T0  not found"<<endl;
    return;
  }
    //   brRec->GetEntry(0);
  
  // Event ------------------------- LOOP  
  for (Int_t ievent=0; ievent<brRec->GetEntries(); ievent++){
    brRec->GetEntry(ievent);
    Int_t   mean = fesdT0->GetT0();
    hMean->Fill(mean);
    Float_t vertex= fesdT0->GetT0zVertex();
      if (ievent <100)  cout<<ievent<<" "<<mean<<" vertex "<<vertex<<endl;
    if(vertex<999){
      hVertex->Fill(vertex);
      amp=fesdT0->GetT0amplitude();
      time=fesdT0->GetT0time();
      for (Int_t i=0; i<24; i++){ 
	hAmp[i]->Fill(amp[i]);
	hTime[i]->Fill(time[i]);
	if (ievent <100)   cout<<ievent<<" "<<i<<" time "<<time[i]<<" amp "<<amp[i]<<endl;
      } 
    }
  }
  
  Hfile = new TFile("FigESD_25081.root","RECREATE","Histograms for T0 
digits");
  printf("Writting histograms to root file \n");
  Hfile->cd();
  //Create a canvas, set the view range, show histograms
  //  g
  //  gStyle->SetOptStat(111111);
  hVertex->Write();
  hMean->Write();
  for (Int_t i=0; i<24; i++){ 
    hAmp[i]->Write();
    hTime[i]->Write();
  }
 
  TCanvas *c1 = new TCanvas("c1", "Time C side",0,48,1280,951);
  c1->Divide(4,3);
  // gStyle->SetOptFit(1111);
  for (Int_t i=0; i<12; i++)
    {
      c1->cd(i+1);
      hTime[i]->Draw();
    }
  TCanvas *c2 = new TCanvas("c2", "Time A side",0,48,1280,951);
  c2->Divide(4,3);
  gStyle->SetOptFit(1111);
  for (Int_t i=12; i<24; i++)
    {
      c2->cd(i+1-12);
      hTime[i]->Draw();
    }
  
  TCanvas *c3 = new TCanvas("c3", "Time mean & vertex",0,48,1280,951);
  c3->Divide(1,2);
  c3->cd(1);
  hVertex->Draw();
  c3->cd(2);
  hMean->Draw();


  

} // end of macro
