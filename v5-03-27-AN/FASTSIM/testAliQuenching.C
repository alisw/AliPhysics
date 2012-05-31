testAliQuenching()
{
  AliQuenchingWeights afq;
  afq.InitMult();
  afq.SetQTransport(1.);
  //afq.SetECMethod(0);
  afq.SampleEnergyLoss();

  Int_t lengths[6]={2,4,6,8,10,12};
  TCanvas *c = new TCanvas("cELW","Energy Loss Weights",0,0,800,500);
  c->Divide(3,2);
  for(Int_t len=0;len<6;len++){
    c->cd(len+1);
    TH1F *h = new TH1F(*afq.GetHisto(2,lengths[len]));
    if(h) h->DrawCopy();
  }
  c->Update();

  afq.PlotDiscreteWeights(4);
  afq.PlotAvgELoss(4);
  afq.PlotAvgELossVsPt(1.,4);
}


