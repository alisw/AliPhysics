plotDiscreteWeight(Int_t len=5)
{
  AliQuenchingWeights afq;
  afq.InitMult();
  afq.PlotDiscreteWeights();
}

plotDiscreteWeightMacro(Int_t len=5)
{
  AliQuenchingWeights afq;
  afq.InitMult();

  TCanvas *c = new TCanvas("cdisc","Discrete Weight for Multiple Scattering",0,0,500,400);
  TH2F *hframe = new TH2F("hdisc","",2,0,1.1,2,0,1);
  hframe->SetStats(0);
  hframe->SetXTitle("#hat{q} [GeV^{2}/fm]");
  hframe->SetYTitle("Probability #Delta E = 0 , p_{0}");
  hframe->Draw();

  TGraph *gq=new TGraph(20);
  Int_t i=0;
  for(Double_t q=0.1;q<=1.05;q+=0.05){
    Double_t disc,cont;
    afq.CalcMult(1,1.0, q, len, cont, disc);
    //cout << " " << q << " " << disc << endl;
    gq->SetPoint(i,q,disc);i++;
  }
  gq->SetMarkerStyle(20);
  gq->Draw("pl");

  TGraph *gg=new TGraph(20);
  Int_t i=0;
  for(Double_t q=0.05;q<=1.05;q+=0.05){
    Double_t disc,cont;
    afq.CalcMult(2,1.0, q, 5., cont, disc);
    //cout << " " << q << " " << disc << endl;
    gg->SetPoint(i,q,disc);i++;
  }
  gg->SetMarkerStyle(24);
  gg->Draw("pl");

  TLegend *l1a = new TLegend(0.5,0.6,.95,0.8);
  l1a->SetFillStyle(0);
  l1a->SetBorderSize(0);
  Char_t label[100];
  sprintf(label,"L = %d fm",len);
  l1a->AddEntry(gq,label,"");
  l1a->AddEntry(gq,"quark","pl");
  l1a->AddEntry(gg,"gluon","pl");
  l1a->Draw();

  c->Update();

}
