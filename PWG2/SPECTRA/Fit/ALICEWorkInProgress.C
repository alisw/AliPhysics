void ALICEWorkInProgress(TCanvas *c,TString today="11/05/2010", TString label = "ALICE performance"){

  TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo",0.72,0.75,0.82,0.89);
  myPadLogo->SetFillColor(0); 
  myPadLogo->SetBorderMode(0);
  myPadLogo->SetBorderSize(2);
  myPadLogo->SetFrameBorderMode(0);
  myPadLogo->SetLeftMargin(0.0);
  myPadLogo->SetTopMargin(0.0);
  myPadLogo->SetBottomMargin(0.0);
  myPadLogo->SetRightMargin(0.0);
  myPadLogo->SetFillStyle(0);
  myPadLogo->Draw();
  myPadLogo->cd();
  TASImage *myAliceLogo = new TASImage("~/WORK/ALICE/ANALYSIS/macros/alice_logo.png");
  myAliceLogo->Draw();
  c->cd();
  TPaveText* t1=new TPaveText(0.65,0.7,0.89,0.75,"NDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->AddText(0.,0.,label);
  t1->SetTextColor(kRed);
  t1->SetTextFont(42);
  t1->Draw();
  TPaveText* t2=new TPaveText(0.65,0.65,0.89,0.7,"NDC");
  t2->SetFillStyle(0);
  t2->SetBorderSize(0);
  t2->SetTextColor(kRed);
  t2->SetTextFont(52);
  t2->AddText(0.,0.,today.Data());
  t2->Draw();
}
