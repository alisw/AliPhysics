
void Draw(Double_t theta){
  //
  // Draw  "primary track distortion"
  //
  //Double_t theta=0.5;
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
  
  fdrfi_rpi0->Draw();
  fdrfi_rpi1->Draw("same");
  fdrfi_rpi2->Draw("same");
  fdrfi_rpi3->Draw("same");

  TLegend *legend = new TLegend(0.15,0.70,0.6,0.85, Form("ExB distortion alog track (tan(#theta)=%f)",theta));
  legend->AddEntry(fdrfi_rpi0, "#phi=0", "l");
  legend->AddEntry(fdrfi_rpi1, "#phi=#pi/2", "l");
  legend->AddEntry(fdrfi_rpi2, "#phi=pi", "l");
  legend->AddEntry(fdrfi_rpi3, "#phi=3*#pi/2", "l");
  legend->Draw();
}
