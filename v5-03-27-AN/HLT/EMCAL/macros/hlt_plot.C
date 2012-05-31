void hlt_plot() {
  
  gStyle->SetPalette(1);
  gStyle->SetOptStat("nemruo");

  TFile *f = new TFile("histofile_local.root");  
  TProfile2D *h[4];
  TH2F *g[4];

  TCanvas *c[3];

  Char_t outname[256];

  for (int i=0; i<4; i++) {
    sprintf(outname, "can%d",i);
    c[i]= new TCanvas(outname,"",0,0,600,400);
    c[i]->Divide(2,2);    
 }
 
  // TCanvas *c1 = new TCanvas("c1","",0,0,600,400);
  // c1->Divide(2,2);

  // TCanvas *c2 = new TCanvas("c2","",0,0,600,400);
  // c2->Divide(2,2);

  // TCanvas *c3 = new TCanvas("c3","",0,0,600,400);
  // c3->Divide(2,2);
  
  
  for (int i=0; i<4; i++) {

    c[0]->cd(i+1);
    h[i] = (TProfile2D*)f->Get(Form("fAmp%d", i));
    h[i]->SetMaximum(1023);
    cout << "fAmp" << i << "entries: " << h[i]->GetEntries()<<endl; 
    h[i]->Draw("colz");

    c[1]->cd(i+1);
    h[i] = (TProfile2D*)f->Get(Form("fTime%d", i));
    h[i]->SetMaximum(15);
    cout << "fTime" << i << "entries: " << h[i]->GetEntries()<<endl;
    h[i]->Draw("colz");
    
    c[2]->cd(i+1);
    g[i] = (TH2F*)f->Get(Form("fAT%d", i));
    cout << "fAT" << i << "entries: " << g[i]->GetEntries()<<endl;
    g[i]->Draw("colz");


    c[3]->cd(i+1);
    g[i] = (TH2F*)f->Get(Form("fCellVsEne%d", i));
    cout << "fCellVsEne" << i << "entries: " << g[i]->GetEntries()<<endl;
    g[i]->Draw("colz");


  }
  
  sprintf(outname, "fAmp.gif");
  c[0]->SaveAs(outname);
  
  sprintf(outname, "fTime.gif");
  c[1]->SaveAs(outname);

  sprintf(outname, "fAT.gif");
  c[2]->SaveAs(outname);

  sprintf(outname, "fCellVsEne.gif");
  c[3]->SaveAs(outname);



 TCanvas *c4 = new TCanvas("c4","",0,0,600,400);
 fCellEne->Draw();
 fClusCellEne->SetLineColor(kRed);
 fClusCellEne->Draw("same");
 sprintf(outname, "ClusCellEne.gif");
 c4->SaveAs(outname);

}
