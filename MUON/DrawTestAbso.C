void DrawTestAbso(Text_t *FileName = "MUONTestAbso")
{
  //  FileName : name of the input file (.root)
  //  Draw histograms to test the absorber correction (mean energy loss + Branson correction)
  //  MUONTestAbso.root can be produced by the mean of the macro MUONTestAbso.C


  strcat(FileName,".root");
  TFile *hfile = (TFile*)gROOT->FindObject(FileName); 
  if (hfile) hfile->Close(); 
  TFile *fill = new TFile(FileName);  
  
  gROOT->Reset();
  gStyle->SetOptStat(11111111);
  gStyle->SetOptFit(1);
  
  
  TCanvas *c1 = new TCanvas("c1",FileName, 250, 20, 800, 820);
  c1->GetFrame()->SetBorderSize(6);
  c1->GetFrame()->SetBorderMode(-1);

  c1->Divide(1,2);

  TPad *pad1 = new TPad("pad1", "The pad with the Histo", 0.01, 0.01, 0.99, 0.99, 21);
  c1->cd(1);
  pad1->Draw();
  pad1->cd();
  pad1->SetGridx();
  pad1->SetGridy();
  pad1->SetFillColor(22);
  pad1->GetFrame()->SetBorderMode(-1);
  hDeltaP1->Draw();

  TPad *pad2 = new TPad("pad2", "The pad with the Histo", 0.01, 0.01, 0.99, 0.99, 21);
  c1->cd(2);
  pad2->Draw();
  pad2->cd();
  pad2->SetGridx();
  pad2->SetGridy();
  pad2->SetFillColor(22);
  pad2->GetFrame()->SetBorderMode(-1);
  hDeltaAngle1->Draw();


  TCanvas *c2 = new TCanvas("c2",FileName, 200, 10, 800, 820);
  c2->GetFrame()->SetBorderSize(6);
  c2->GetFrame()->SetBorderMode(-1);

  c2->Divide(1,2);

  TPad *pad3 = new TPad("pad3", "The pad with the Histo", 0.01, 0.01, 0.99, 0.99, 21);
  c2->cd(1);
  pad3->Draw();
  pad3->cd();
  pad3->SetGridx();
  pad3->SetGridy();
  pad3->SetFillColor(22);
  pad3->GetFrame()->SetBorderMode(-1);
  hDeltaP2->Draw();

  TPad *pad4 = new TPad("pad4", "The pad with the Histo", 0.01, 0.01, 0.99, 0.99, 21);
  c2->cd(2);
  pad4->Draw();
  pad4->cd();
  pad4->SetGridx();
  pad4->SetGridy();
  pad4->SetFillColor(22);
  pad4->GetFrame()->SetBorderMode(-1);
  hDeltaAngle2->Draw();

  TCanvas *c3 = new TCanvas("c3",FileName, 150, 5, 800, 700);

  c3->GetFrame()->SetBorderSize(6);
  c3->GetFrame()->SetBorderMode(-1);

  TPad *pad5 = new TPad("pad5", "The pad with the Histo", 0.01, 0.01, 0.99, 0.99, 21);
  pad5->Draw();
  pad5->cd();
  pad5->SetGridx();
  pad5->SetGridy();
  pad5->SetFillColor(22);
  pad5->GetFrame()->SetBorderMode(-1);

  g1= new TF1("g1","gaus",9.3,9.8) ; // 9.25
  hInvMassRes->Fit("g1","RQ");
  hInvMassRes->GetXaxis()->SetTitleFont(20);
  hInvMassRes->SetXTitle("Mass (GeV/c^2!)");
  hInvMassRes->SetFillColor(63);
  hInvMassRes->Draw();

}
