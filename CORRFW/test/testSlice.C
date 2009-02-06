
void testSlice() {
  
  gStyle->SetPalette(1);

  gSystem->Load("libANALYSIS");
  gSystem->Load("libCORRFW");

  TFile * f = TFile::Open("output.root");

  AliCFContainer * c_in = (AliCFContainer*)f->Get("container");
  
  printf("====\n");
  printf("old container properties\n");
  printf("nvar=%d\t nstep=%d\t nbins[0]=%d\t nbins[1]=%d\n",c_in->GetNVar(),c_in->GetNStep(),c_in->GetNBins(0),c_in->GetNBins(1));
  printf("====\n");


  TCanvas * can = new TCanvas("can","",1300,850);
  can->Divide(3,2);
  Short_t iCan=1;

  can->cd(iCan++);
  TH2D* h_in_step0 = c_in->ShowProjection(0,1,0);
  h_in_step0->SetTitle("input container - step 0");
  h_in_step0->Draw("text colz");

  can->cd(iCan++);
  TH2D* h_in_step1 = c_in->ShowProjection(0,1,1) ;
  h_in_step1->SetTitle("input container - step 1");
  h_in_step1->Draw("text colz");
  

  Int_t* vars = new Int_t[1 /*2*/];
  vars[0]=0;
  //vars[1]=1;

  Double_t epsilon=1.e-07;

  Double_t varMin[2]={0.0 ,            0.0  };
  Double_t varMax[2]={2.0 - epsilon ,  1.0  };


  can->cd(iCan++);
  TH1D* h_slice_step0 = c_in->ShowSlice(0,varMin,varMax,0);
  h_slice_step0->SetTitle("test slice - step 0");
  h_slice_step0->Draw();

  printf("Creating Slice...\n");
  AliCFContainer* c_out = c_in->MakeSlice(1,vars,varMin,varMax);

  printf("====\n");
  printf("new container properties\n");
  printf("nvar=%d\t nstep=%d\t nbins[0]=%d\n",c_out->GetNVar(),c_out->GetNStep(),c_out->GetNBins(0));
  printf("range = %d  -> %d\n",((AliCFGridSparse*)c_out->GetGrid(0))->GetGrid()->GetAxis(0)->GetFirst(),((AliCFGridSparse*)c_out->GetGrid(0))->GetGrid()->GetAxis(0)->GetLast());
  printf("====\n");
  
  can->cd(iCan++);
  TH1D* h_out_step0 = c_out->ShowProjection(0,0);
  h_out_step0->SetTitle("output container - step 0");
  h_out_step0->Draw();

  can->cd(iCan++);
  TH1D* h_out_step1 = c_out->ShowProjection(0,1);
  h_out_step1->SetTitle("output container - step 1");
  h_out_step1->Draw();



  AliCFEffGrid * eff = new AliCFEffGrid("eff","",*c_out);
  eff->CalculateEfficiency(1,0);

  printf("====\n");
  printf("efficiency map properties\n");
  printf("nvar=%d\t nbins[0]=%d\n",eff->GetNVar(),eff->GetNBins(0));
  printf("====\n");
  
  can->cd(iCan++);
  eff->Project(0)->Draw();

}
