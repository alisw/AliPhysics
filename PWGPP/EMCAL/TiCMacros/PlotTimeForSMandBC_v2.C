void PlotTimeForSMandBC_v2(TString period="LHC15n"){
  //  TFile *file=new TFile("AnalysisResults_step2Merged.root");
  TFile *file=new TFile("AnalysisResults.root");
  TList *list=file->Get("chistolist");
  //  list->ls();
  TH2F *hSMBC[80];
  TH2F *hSM[20];
  TH2F *hSMBCLG[80];
  TH2F *hSMLG[20];

  char text[256];
  for(Int_t j=0;j<20;j++){//SM loop
    for(Int_t i=0;i<4;i++){//BC loop
      sprintf(text,"SupMod%dBC%d",j,i);
      hSMBC[j+i*20]=(TH2F *)list->FindObject(text);
      sprintf(text,"SupMod%dBC%dLG",j,i);
      hSMBCLG[j+i*20]=(TH2F *)list->FindObject(text);
    }
    sprintf(text,"SupMod%d",j);
    hSM[j]=(TH2F *)list->FindObject(text);
    sprintf(text,"SupMod%dLG",j);
    hSMLG[j]=(TH2F *)list->FindObject(text);
  }

  TH2F *hAll=hSM[0]->Clone("hAll");
  hAll->Sumw2();
  for(Int_t j=1;j<20;j++){//SM loop
    hSM[j]->Sumw2();
    hAll->Add(hSM[j],1.);
  }
  TH1F *hAll_pro=(TH1F*)hAll->ProjectionY("hAll_pro",1,100);
  TH1F *hAll_pro1=(TH1F*)hAll->ProjectionY("hAll_pro1",1,5);
  TH1F *hAll_pro2=(TH1F*)hAll->ProjectionY("hAll_pro2",6,10);
  TH1F *hAll_pro3=(TH1F*)hAll->ProjectionY("hAll_pro3",11,25);
  TH1F *hAll_pro4=(TH1F*)hAll->ProjectionY("hAll_pro4",26,100);


    
  TCanvas *c1=new TCanvas("c1","c1",800,600); 
  c1->Divide(5,4);
  for(Int_t i=0;i<4;i++){//BC loop
    for(Int_t j=0;j<20;j++){//SM loop
      c1->cd(j+1)->SetLogz();
      if(hSMBC[j+i*20])
	hSMBC[j+i*20]->Draw("colz");
    }
    sprintf(text,"plots/AllSM_pass2_%s_BC%d.pdf",period.Data(),i);
    c1->Print(text);
  }
  for(Int_t j=0;j<20;j++){//SM loop
    c1->cd(j+1)->SetLogz();
    if(hSM[j])
      hSM[j]->Draw("colz");
  }
  sprintf(text,"plots/AllSMmerged_pass2_%s.pdf",period.Data());
  c1->Print(text);

  TCanvas *c2=new TCanvas("c2","c2",800,600); 
  c2->Divide(5,4);
  for(Int_t i=0;i<4;i++){//BC loop
    for(Int_t j=0;j<20;j++){//SM loop
      c2->cd(j+1)->SetLogz();
      if(hSMBCLG[j+i*20])
	hSMBCLG[j+i*20]->Draw("colz");
    }
    sprintf(text,"plots/AllSM_pass2_%s_BC%d_LG.pdf",period.Data(),i);
    c2->Print(text);
  }
  for(Int_t j=0;j<20;j++){//SM loop
    c2->cd(j+1)->SetLogz();
    if(hSMLG[j])
      hSMLG[j]->Draw("colz");
  }
  sprintf(text,"plots/AllSMmerged_pass2_%s_LG.pdf",period.Data());
  c2->Print(text);

  //cout<<  hAll->GetNbinsX()<<endl;//100 bins in a range 0-20 GeV
  //1 bin is 0.2 GeV/c
  TCanvas *c3=new TCanvas("c3","c3",800,600);
  gPad->SetLogy();
  gStyle->SetOptFit();
  hAll_pro->SetMarkerStyle(20);
  hAll_pro->Fit("gaus","","",-20.,20.);
  hAll_pro->Draw();
  sprintf(text,"plots/AllmergedFit_%s.pdf",period.Data());
  c3->Print(text);

  TCanvas *c3a=new TCanvas("c3a","c3a",800,600);
  gPad->SetLogz();
  //  gStyle->SetOptFit();
  //hAll_pro->SetMarkerStyle(20);
  //hAll_pro->Fit("gaus","","",-20.,20.);
  hAll->Draw("colz");
  sprintf(text,"plots/Allmerged_2D_%s.pdf",period.Data());
  c3a->Print(text);

  
  hAll_pro1->SetMarkerStyle(20);
  hAll_pro2->SetMarkerStyle(20);
  hAll_pro3->SetMarkerStyle(20);
  hAll_pro4->SetMarkerStyle(20);

  TCanvas *c4=new TCanvas("c4","c4",800,600);
  c4->Divide(2,2);
  gStyle->SetOptFit();
  c4->cd(1);
  gPad->SetLogy();
  hAll_pro1->Fit("gaus","","",-20.,20.);
  hAll_pro1->Draw();
  c4->cd(2);
  gPad->SetLogy();
  hAll_pro2->Fit("gaus","","",-20.,20.);
  hAll_pro2->Draw();
  c4->cd(3);
  gPad->SetLogy();
  hAll_pro3->Fit("gaus","","",-20.,20.);
  hAll_pro3->Draw();
  c4->cd(4);
  gPad->SetLogy();
  hAll_pro4->Fit("gaus","","",-20.,20.);
  hAll_pro4->Draw();
  sprintf(text,"plots/AllmergedFit_pTbins_%s.pdf",period.Data());
  c4->Print(text);

  
}
