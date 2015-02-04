/// \file AliTPCExBdraw.C
//////
/// ~~~{.cpp}
/// char *storage = "local://OCDBres"
/// Int_t RunNumber=0;
/// AliCDBManager::Instance()->SetDefaultStorage(storage);
/// AliCDBManager::Instance()->SetRun(RunNumber) 
/// AliTPCExBFirst * exb = AliTPCcalibDB::Instance()->GetExB();
/// // See example macro $ALICE_ROOT/TPC/macros/AliTPCExBdraw.C
/// .L $ALICE_ROOT/TPC/macros/AliTPCExBdraw.C 
/// DrawPrim(0)
/// DrawLaser(0)
/// ~~~

void DrawPrim(Double_t theta,Float_t magf=5){
  /// Draw  "primary track distortion"
  ///
  /// Double_t theta=0.5;

  TF1 *fdrfi_rpi0=new TF1("fdrfi_rpi0",Form("AliTPCExB::GetDrphi(x,0*pi/2,%f*x)",theta),20,250);
  TF1 *fdrfi_rpi1=new TF1("fdrfi_rpi1",Form("AliTPCExB::GetDrphi(x,1*pi/2,%f*x)",theta),20,250);
  TF1 *fdrfi_rpi2=new TF1("fdrfi_rpi2",Form("AliTPCExB::GetDrphi(x,2*pi/2,%f*x)",theta),20,250);
  TF1 *fdrfi_rpi3=new TF1("fdrfi_rpi3",Form("AliTPCExB::GetDrphi(x,3*pi/2,%f*x)",theta),20,250);

  fdrfi_rpi0->GetXaxis()->SetTitle("r (cm)");
  fdrfi_rpi0->GetYaxis()->SetTitle("drd#phi (cm)");
  
  fdrfi_rpi0->SetLineColor(2);
  fdrfi_rpi1->SetLineColor(3);
  fdrfi_rpi2->SetLineColor(4);
  fdrfi_rpi3->SetLineColor(5);
  char leg0[1000];
  char leg1[1000];
  char leg2[1000];
  char leg3[1000];
  Float_t conv= 1/(magf*0.299792458e-3);
  fdrfi_rpi0->GetHistogram()->Fit("pol2","same","",90,250);
  sprintf(leg0,"#phi=0*#pi/2      Track res: dy=%.2f mm d#phi=%.2f mrad dC = %.6f 1/GeV",
	  10*pol2->GetParameter(0),1000*pol2->GetParameter(1), conv*pol2->GetParameter(2));
  fdrfi_rpi1->GetHistogram()->Fit("pol2","same","",90,250);
  sprintf(leg1,"#phi=1*#pi/2      Track res: dy=%.2f mm d#phi=%.2f mrad dC = %.6f 1/GeV",
	  10*pol2->GetParameter(0),1000*pol2->GetParameter(1),conv*pol2->GetParameter(2));
  fdrfi_rpi2->GetHistogram()->Fit("pol2","same","",90,250);
  sprintf(leg2,"#phi=2*#pi/2      Track res: dy=%.2f mm d#phi=%.2f mrad dC = %.6f 1/GeV",
	  10*pol2->GetParameter(0),1000*pol2->GetParameter(1),conv*pol2->GetParameter(2));
  fdrfi_rpi3->GetHistogram()->Fit("pol2","same","",90,250);
  sprintf(leg3,"#phi=3*#pi/2      Track res: dy=%.2f mm d#phi=%.2f mrad dC = %.6f 1/GeV",
	  10*pol2->GetParameter(0),1000*pol2->GetParameter(1),conv*pol2->GetParameter(2));
  fdrfi_rpi0->SetMaximum(fdrfi_rpi0->GetMaximum()*2.0);
  fdrfi_rpi0->Draw();
  fdrfi_rpi1->Draw("same");
  fdrfi_rpi2->Draw("same");
  fdrfi_rpi3->Draw("same");

  TLegend *legend = new TLegend(0.10,0.50,0.8,0.90, Form("ExB distortion along track (tan(#theta)=%f)",theta));
  legend->AddEntry(fdrfi_rpi0, leg0, "l");
  legend->AddEntry(fdrfi_rpi1, leg1, "l");
  legend->AddEntry(fdrfi_rpi2, leg2, "l");
  legend->AddEntry(fdrfi_rpi3, leg3, "l");
  legend->Draw();
}



void DrawLaser(Double_t fipi2,Float_t magf=5){
  /// Draw  "primary track distortion"
  ///
  /// Double_t fipi2=0.5;

  TF2 *fdistout = new TF1("fdistout","[0]+[1]*(x-250)+[2]*(x-250)*(x-250)",90,250);

  TF1 *fdrfi_rpi0=new TF1("fdrfi_rpi0",Form("AliTPCExB::GetDrphi(x,%f*pi/2,250-20)",fipi2),20,250);
  TF1 *fdrfi_rpi1=new TF1("fdrfi_rpi1",Form("AliTPCExB::GetDrphi(x,%f*pi/2,250-90)",fipi2),20,250);
  TF1 *fdrfi_rpi2=new TF1("fdrfi_rpi2",Form("AliTPCExB::GetDrphi(x,%f*pi/2,250-160)",fipi2),20,250);
  TF1 *fdrfi_rpi3=new TF1("fdrfi_rpi3",Form("AliTPCExB::GetDrphi(x,%f*pi/2,250-230)",fipi2),20,250);

  fdrfi_rpi0->GetXaxis()->SetTitle("r (cm)");
  fdrfi_rpi0->GetYaxis()->SetTitle("drd#phi (cm)");
  
  fdrfi_rpi0->SetLineColor(2);
  fdrfi_rpi1->SetLineColor(3);
  fdrfi_rpi2->SetLineColor(4);
  fdrfi_rpi3->SetLineColor(5);
  char leg0[1000];
  char leg1[1000];
  char leg2[1000];
  char leg3[1000];
  Float_t conv= 1/(magf*0.299792458e-3);
  fdrfi_rpi0->GetHistogram()->Fit("fdistout","same","",90,250);
  sprintf(leg0,"laser 0      Track res: dy=%.2f mm d#phi=%.2f mrad dC = %.6f 1/GeV",
	  10*fdistout->GetParameter(0),1000*fdistout->GetParameter(1), conv*fdistout->GetParameter(2));
  fdrfi_rpi1->GetHistogram()->Fit("fdistout","same","",90,250);
  sprintf(leg1,"laser 1      Track res: dy=%.2f mm d#phi=%.2f mrad dC = %.6f 1/GeV",
	  10*fdistout->GetParameter(0),1000*fdistout->GetParameter(1),conv*fdistout->GetParameter(2));
  fdrfi_rpi2->GetHistogram()->Fit("fdistout","same","",90,250);
  sprintf(leg2,"laser 2      Track res: dy=%.2f mm d#phi=%.2f mrad dC = %.6f 1/GeV",
	  10*fdistout->GetParameter(0),1000*fdistout->GetParameter(1),conv*fdistout->GetParameter(2));
  fdrfi_rpi3->GetHistogram()->Fit("fdistout","same","",90,250);
  sprintf(leg3,"laser 3      Track res: dy=%.2f mm d#phi=%.2f mrad dC = %.6f 1/GeV",
	  10*fdistout->GetParameter(0),1000*fdistout->GetParameter(1),conv*fdistout->GetParameter(2));
  fdrfi_rpi0->SetMaximum(fdrfi_rpi3->GetMaximum()*2.0);
  fdrfi_rpi0->Draw();
  fdrfi_rpi1->Draw("same");
  fdrfi_rpi2->Draw("same");
  fdrfi_rpi3->Draw("same");

  TLegend *legend = new TLegend(0.10,0.50,0.8,0.90, Form("ExB distortion along track (fi=%f*#pi/2))",fipi2));
  legend->AddEntry(fdrfi_rpi0, leg0, "l");
  legend->AddEntry(fdrfi_rpi1, leg1, "l");
  legend->AddEntry(fdrfi_rpi2, leg2, "l");
  legend->AddEntry(fdrfi_rpi3, leg3, "l");
  legend->Draw();
}



void DrawDpt(Float_t d1pt, Float_t ptmin, Float_t ptmax){
  ///

}
