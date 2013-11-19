void EmpiricalGEMQFit(){
  //
  // Empirical fit of the relative Q resolution for the iron source with 4 GEM layers
  // Fit:
  // Assumption : 
  // 1.) Effective electron transparency proportiaonal to the effective ion transparency
  // 2.) RMS of the effecitive gain can be expressed as linear function of U1/U3,U2/U4 and U3/U4
  //
  // Agreement with data within the expect error estimate -relative agreement 
  // RMS (sigma_{meas}/sigma_{fit}) ~ 3%
  // 
  TFile *fgem = TFile::Open("GEMScansIKF.root");
  TTree * tree= (TTree*)f.Get("NeCO2N2");
  tree->SetMarkerStyle(25);
  TStatToolkit toolkit;
  Double_t chi20=0;
  Int_t    npoints=0;
  Int_t npointsMax=10000;
  TVectorD param0,param1,param2,param3;
  TMatrixD covar0,covar1,covar2,covar3;
  
  tree->SetAlias("UGEMA","(UGEM1+UGEM2+UGEM3+UGEM4)");
  TString fstringFast="";
  fstringFast+="1/IB++";                // fraction of the 
  fstringFast+="(UGEM1/UGEMA)/IB++";    // fraction of the gain
  fstringFast+="(UGEM2/UGEMA)/IB++";    // fraction of the gain
  fstringFast+="(UGEM3/UGEMA)/IB++";    // fraction of the gain
  //
  TString *strResolFit = TStatToolkit::FitPlane(tree,"Sigma^2:Sigma^2", fstringFast.Data(),cutFit, chi20,npoints,param0,covar0,-1,0, npointsMax, 0);
  strResolFit->Tokenize("++")->Print();
  tree->SetAlias("fitSigma2",strResolFit->Data());
  //
  gStyle->SetOptTitle(0);
  //
  TCanvas *canvas = new TCanvas("canvasEmp","canvasEmp",600,500);
  canvas->SetBottomMargin(0.15);
  canvas->SetRightMargin(0.1);
  canvas->SetTicks(1,1);
  tree->SetMarkerSize(0.7)
  tree->Draw("Sigma:sqrt(fitSigma2):sqrt(1/IB)",cutFit,"colz");  
  TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
  htemp->GetXaxis()->SetTitle("#sigma_{fit}(%)");
  htemp->GetYaxis()->SetTitle("#sigma_{meas}(%)");
  htemp->GetZaxis()->SetTitle("#sqrt{IBF}");
  htemp->SetTitle("Fe resolution");
  htemp->GetXaxis()->SetLimits(8,20);
  htemp->GetYaxis()->SetLimits(8,20);
  htemp->Draw("colz");
  TLatex latex;
  latex.DrawLatex(8.5,18.5,"#sigma=#sqrt{p_{0}+p_{1}#frac{1}{IB}+p_{2}#frac{U1}{U4xIB}+p_{3}#frac{U2}{U4xIB} + p_{4}#frac{U3}{U4xIB}}");
  latex.DrawLatex(8.5,17,TString::Format("p_{0}=%1.f",param0[0]));
  latex.DrawLatex(8.5,16,TString::Format("p_{1}=%1.f",param0[1]));
  latex.DrawLatex(8.5,15,TString::Format("p_{2}=%1.f",param0[2]));
  latex.DrawLatex(8.5,14,TString::Format("p_{3}=%1.f",param0[3]));
  latex.DrawLatex(8.5,13,TString::Format("p_{4}=%1.f",param0[4]));

  canvasFit->SaveAs("canvasFEResolutionFit.pdf");
  canvasFit->SaveAs("canvasFEResolutionFit.png");
}
