void SetStyles(TH1 *histo,int marker, int color,char *xtitle, char *ytitle);
TH1* bayneseffdiv(TH1* numerator, TH1* denominator,Char_t* name) ;
void PlotBackgroundClusters(TString filename="rootFiles/LHC11a4_bis/Et.ESD.simPbPb.EMCAL.LHC11a4_bis.root", Bool_t peripheral = kFALSE, int bin = 10, int binLast = 10, TString det = "EMCal"){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TString outname = "";
  TString outnamebin = Form("%iTo%i",bin,binLast);

  TString suffix = "";
  if(peripheral) suffix ="Peripheral";

  TFile *f = TFile::Open(filename, "READ");
  
  TList *l = (TList*)f->Get("out1");
  TH1F  *fHistChargedTracksCut = l->FindObject(("fHistChargedTracksCut"+suffix).Data());
  TH1F  *fHistChargedTracksAccepted = l->FindObject(("fHistChargedTracksAccepted"+suffix).Data());
  TH1F  *fHistGammasCut = l->FindObject(("fHistGammasTracksCut"+suffix).Data());
  TH1F  *fHistGammasAccepted = l->FindObject(("fHistGammasTracksAccepted"+suffix).Data());

  int rebin = 1;
  fHistChargedTracksCut->Rebin(rebin);
  fHistChargedTracksAccepted->Rebin(rebin);
  fHistGammasCut->Rebin(rebin);
  fHistGammasAccepted->Rebin(rebin);

  int bin1 = 1;
  int bin2 = fHistChargedTracksCut->GetXaxis()->FindBin(0.5);
  int bin3 = fHistChargedTracksCut->GetNbinsX();
  float nGammasCutBelow = fHistGammasCut->Integral(bin1,bin2);
  float nGammasCutAbove = fHistGammasCut->Integral(bin2+1,bin3);
  float nGammasCut = fHistGammasCut->Integral();
  float nTracksCutBelow = fHistChargedTracksCut->Integral(bin1,bin2);
  float nTracksCutAbove = fHistChargedTracksCut->Integral(bin2+1,bin3);
  float nTracksCut = fHistChargedTracksCut->Integral();
  cout<<"Gamma Fraction above: "<<nGammasCutAbove/nGammasCut<<" below: "<<nGammasCutBelow/nGammasCut<<endl;
  cout<<"Track Fraction above: "<<nTracksCutAbove/nTracksCut<<" below: "<<nTracksCutBelow/nTracksCut<<endl;

  TH1F *hTotalCut = fHistChargedTracksCut->Clone("hTotalCut");
  hTotalCut->Add(fHistGammasCut);
  TH1F *hTotalAccepted = fHistChargedTracksAccepted->Clone("hTotalAccepted");
  hTotalAccepted->Add(fHistGammasAccepted);
  TH1F *hTotalClusters = fHistChargedTracksCut->Clone("hTotalClusters");
  hTotalClusters->Add(fHistChargedTracksAccepted);
  hTotalClusters->Add(fHistGammasCut);
  hTotalClusters->Add(fHistGammasAccepted);
  TH1F *hTotalCharged = fHistChargedTracksCut->Clone("hTotalCharged");
  hTotalCharged->Add(fHistChargedTracksAccepted);
  TH1F *hTotalGamma = fHistGammasCut->Clone("hTotalGamma");
  hTotalGamma->Add(fHistGammasAccepted);

  TH1F *hSignalCut = bayneseffdiv(fHistGammasCut,hTotalCut,"hSignalCut");
  TH1F *hBkgdCut = bayneseffdiv(fHistChargedTracksCut,hTotalCut,"hBkgdCut");
  TH1F *hSignalAcc = bayneseffdiv(fHistGammasAccepted,hTotalAccepted,"hSignalAcc");
  TH1F *hBkgdAcc = bayneseffdiv(fHistChargedTracksAccepted,hTotalAccepted,"hBkgdAcc");
  TH1F *hFracCharged = bayneseffdiv(hTotalCharged,hTotalClusters,"hFracCharged");
  TH1F *hFracGamma = bayneseffdiv(hTotalGamma,hTotalClusters,"hFracGamma");
  SetStyles(hSignalCut,21,2,"E_{T}","fraction cut");
  SetStyles(hSignalAcc,25,2,"E_{T}","fraction cut");
  SetStyles(hBkgdCut,20,4,"E_{T}","fraction cut");
  SetStyles(hBkgdAcc,24,4,"E_{T}","fraction cut");
  SetStyles(hFracCharged,24,4,"E_{T}","fraction cut");
  SetStyles(hFracGamma,25,2,"E_{T}","fraction cut");

  float ECut = 0.5;
  TLine *lineEDep = new TLine(ECut,0,ECut,1);
  lineEDep->SetLineWidth(2);
  TLine *lineEDep2D = new TLine(ECut,0,ECut,3);
  lineEDep2D->SetLineWidth(2);

  TLegend *leg = new TLegend(0.357383,0.47043,0.458054,0.75);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(hSignalCut,"Fraction of cut particles that are signals");
  leg->AddEntry(hSignalAcc,"Fraction of accepted particles that are signals");
  leg->AddEntry(hBkgdCut,"Fraction of cut particles that are bkgd");
  leg->AddEntry(hBkgdAcc,"Fraction of accepted particles that are background");
  leg->SetTextSize(0.061828);

  TCanvas *c1 = new TCanvas("c1","c1",600,400);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  c1->SetBorderSize(0);
  c1->SetFillColor(0);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);
  hSignalCut->SetMaximum(1.0);
  hSignalCut->Draw();
  hSignalAcc->Draw("same");
  hBkgdCut->Draw("same");
  hBkgdAcc->Draw("same");
  lineEDep->Draw();
  leg->Draw();

  c1->SaveAs("/tmp/SignalBkgdEmcal.png");

  TLegend *leg1 = new TLegend(0.357383,0.47043,0.458054,0.75);
  leg1->SetFillStyle(0);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->AddEntry(hFracCharged,"Fraction of clusters from charged particles");
  leg1->AddEntry(hFracGamma,"Fraction of clusters from gammas");
  leg1->SetTextSize(0.061828);

  TCanvas *c2 = new TCanvas("c2","c2",600,400);
  c2->SetTopMargin(0.02);
  c2->SetRightMargin(0.02);
  c2->SetBorderSize(0);
  c2->SetFillColor(0);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetFrameFillColor(0);
  c2->SetFrameBorderMode(0);
//   hFracCharged->SetMaximum(1.0);
//   hFracCharged->Draw();
//   hFracGamma->Draw("same");
//   leg1->Draw();

//   c2->SaveAs("/tmp/FractionsEmcal.png");


  TH2F  *fHistMatchedTracksEvspTBkgd = l->FindObject(("fHistMatchedTracksEvspTBkgd"+suffix).Data());
  TProfile * profBkgd = fHistMatchedTracksEvspTBkgd->ProfileX();
  TH2F  *fHistMatchedTracksEvspTSignal = l->FindObject(("fHistMatchedTracksEvspTSignal"+suffix).Data());
  TProfile * profSignal = fHistMatchedTracksEvspTSignal->ProfileX();
  TF1 *func = new TF1("func","x",0,3);
  fHistMatchedTracksEvspTBkgd->GetXaxis()->SetTitle("p_{T}");
  fHistMatchedTracksEvspTBkgd->GetYaxis()->SetTitle("E_{cluster}");
  fHistMatchedTracksEvspTSignal->GetXaxis()->SetTitle("p_{T}");
  fHistMatchedTracksEvspTSignal->GetYaxis()->SetTitle("E_{cluster}");


  TF1 *funcAvgSig = new TF1("funcAvgSig","[0]*x+[1]",0,3);
  funcAvgSig->SetParameter(0,1);
  funcAvgSig->SetParameter(1,0.3);

  TF1 *funcCut = new TF1("funcCut","[0]*x+[1]",0,3);
  funcCut->SetParameter(0,9.54979e-01);
  funcCut->SetParameter(1,-2.47491e-01);


  TH3F  *fHistMatchedTracksEvspTBkgdvsMult = l->FindObject("fHistMatchedTracksEvspTBkgdMult");
  TH3F  *fHistMatchedTracksEvspTSignalvsMult = l->FindObject("fHistMatchedTracksEvspTSignalMult");
  //DoProjectProfile2D(const char* name, const char* title, TAxis* projX, TAxis* projY, bool originalRange, bool useUF, bool useOF) const
  fHistMatchedTracksEvspTBkgdvsMult->GetZaxis()->SetRange(bin,binLast);
  TH2D *hBkgd2D = (TProfile2D*) fHistMatchedTracksEvspTBkgdvsMult->Project3D("yx");
  TProfile * profBkgd2D = hBkgd2D->ProfileX();
  fHistMatchedTracksEvspTSignalvsMult->GetZaxis()->SetRange(bin,binLast);
  TH2D *hSignal2D = (TProfile2D*) fHistMatchedTracksEvspTSignalvsMult->Project3D("yx");
  TProfile * profSignal2D = hSignal2D->ProfileX();

  TH2F  *fHistChargedTracksCutMult = l->FindObject("fHistChargedTracksCutMult");
  TH2F  *fHistChargedTracksAcceptedMult = l->FindObject("fHistChargedTracksAcceptedMult");
  TH2F  *fHistGammasCutMult = l->FindObject("fHistGammasTracksCutMult");
  TH2F  *fHistGammasAcceptedMult = l->FindObject("fHistGammasTracksAcceptedMult");

  TH1D  *fHistChargedTracksCutPeri = fHistChargedTracksCutMult->ProjectionX("fHistChargedTracksCutPeri",bin,binLast);
  TH1D  *fHistChargedTracksAcceptedPeri = fHistChargedTracksAcceptedMult->ProjectionX("fHistChargedTracksAcceptedPeri",bin,binLast);
  TH1D  *fHistGammasCutPeri = fHistGammasCutMult->ProjectionX("fHistGammasTracksCutPeri",bin,binLast);
  TH1D  *fHistGammasAcceptedPeri = fHistGammasAcceptedMult->ProjectionX("fHistGammasTracksAcceptedPeri",bin,binLast);

  int rebin = 1;
  fHistChargedTracksCutPeri->Rebin(rebin);
  fHistChargedTracksAcceptedPeri->Rebin(rebin);
  fHistGammasCutPeri->Rebin(rebin);
  fHistGammasAcceptedPeri->Rebin(rebin);

  TH1D *hTotalCutPeri = fHistChargedTracksCutPeri->Clone("hTotalCutPeri");
  hTotalCutPeri->Add(fHistGammasCutPeri);
  TH1D *hTotalAcceptedPeri = fHistChargedTracksAcceptedPeri->Clone("hTotalAcceptedPeri");
  hTotalAcceptedPeri->Add(fHistGammasAcceptedPeri);
  TH1D *hTotalClustersPeri = fHistChargedTracksCutPeri->Clone("hTotalClustersPeri");
  hTotalClustersPeri->Add(fHistChargedTracksAcceptedPeri);
  hTotalClustersPeri->Add(fHistGammasCutPeri);
  hTotalClustersPeri->Add(fHistGammasAcceptedPeri);
  TH1D *hTotalChargedPeri = fHistChargedTracksCutPeri->Clone("hTotalChargedPeri");
  hTotalChargedPeri->Add(fHistChargedTracksAcceptedPeri);
  TH1D *hTotalGammaPeri = fHistGammasCutPeri->Clone("hTotalGammaPeri");
  hTotalGammaPeri->Add(fHistGammasAcceptedPeri);

  TH1D *hSignalCutPeri = bayneseffdiv(fHistGammasCutPeri,hTotalCutPeri,"hSignalCutPeri");
  TH1D *hBkgdCutPeri = bayneseffdiv(fHistChargedTracksCutPeri,hTotalCutPeri,"hBkgdCutPeri");
  TH1D *hSignalAccPeri = bayneseffdiv(fHistGammasAcceptedPeri,hTotalAcceptedPeri,"hSignalAccPeri");
  TH1D *hBkgdAccPeri = bayneseffdiv(fHistChargedTracksAcceptedPeri,hTotalAcceptedPeri,"hBkgdAccPeri");
  TH1D *hFracChargedPeri = bayneseffdiv(hTotalChargedPeri,hTotalClustersPeri,"hFracChargedPeri");
  TH1D *hFracGammaPeri = bayneseffdiv(hTotalGammaPeri,hTotalClustersPeri,"hFracGammaPeri");
  SetStyles(hSignalCutPeri,21,2,"E_{T}","fraction cut");
  SetStyles(hSignalAccPeri,25,2,"E_{T}","fraction cut");
  SetStyles(hBkgdCutPeri,20,4,"E_{T}","fraction cut");
  SetStyles(hBkgdAccPeri,24,4,"E_{T}","fraction cut");
  SetStyles(hFracChargedPeri,24,4,"E_{T}","fraction cut");
  SetStyles(hFracGammaPeri,25,2,"E_{T}","fraction cut");

  TLegend *legPeri = new TLegend(0.357383,0.47043,0.458054,0.75);
  legPeri->SetFillStyle(0);
  legPeri->SetFillColor(0);
  legPeri->SetBorderSize(0);
  legPeri->SetTextSize(0.03);
  legPeri->AddEntry(hSignalCutPeri,"Fraction of cut particles that are signals");
  legPeri->AddEntry(hSignalAccPeri,"Fraction of accepted particles that are signals");
  legPeri->AddEntry(hBkgdCutPeri,"Fraction of cut particles that are bkgd");
  legPeri->AddEntry(hBkgdAccPeri,"Fraction of accepted particles that are background");
  legPeri->SetTextSize(0.061828);

  //profSignal->Fit(funcAvgSig,"","",ECut,3);
  funcAvgSig->SetLineStyle(2);
  funcAvgSig->SetLineColor(1);
  funcCut->SetLineStyle(2);
  funcCut->SetLineColor(1);

  TCanvas *c3 = new TCanvas("c3","c3",600,400);
  c3->SetTopMargin(0.02);
  c3->SetRightMargin(0.02);
  c3->SetBorderSize(0);
  c3->SetFillColor(0);
  c3->SetFillColor(0);
  c3->SetBorderMode(0);
  c3->SetFrameFillColor(0);
  c3->SetFrameBorderMode(0);
  fHistMatchedTracksEvspTSignal->Draw("colz");
  profSignal->Draw("esame");
  profBkgd->Draw("esame");
  func->Draw("same");
  //funcCut->Draw("same");
  lineEDep2D->Draw();


  TCanvas *c10 = new TCanvas("c10","c10",600,400);
  c10->SetTopMargin(0.02);
  c10->SetRightMargin(0.02);
  c10->SetBorderSize(0);
  c10->SetFillColor(0);
  c10->SetFillColor(0);
  c10->SetBorderMode(0);
  c10->SetFrameFillColor(0);
  c10->SetFrameBorderMode(0);
  int firstbin = fHistMatchedTracksEvspTBkgd->GetXaxis()->FindBin(0.75);
  TProfile * profE = fHistMatchedTracksEvspTBkgd->ProfileY("Test",firstbin,firstbin+1);
  profE->Draw();

  TCanvas *c4 = new TCanvas("c4","c4",600,400);
  c4->SetTopMargin(0.02);
  c4->SetRightMargin(0.02);
  c4->SetBorderSize(0);
  c4->SetFillColor(0);
  c4->SetFillColor(0);
  c4->SetBorderMode(0);
  c4->SetFrameFillColor(0);
  c4->SetFrameBorderMode(0);
  fHistMatchedTracksEvspTBkgd->Draw("colz");
  profSignal->Draw("esame");
  profBkgd->Draw("esame");
  func->Draw("same");
  //funcCut->Draw("same");
  //funcAvgSig->Draw("same");
  lineEDep2D->Draw();

  TCanvas *c7 = new TCanvas("c7","c7",600,400);
  c7->SetTopMargin(0.02);
  c7->SetRightMargin(0.02);
  c7->SetBorderSize(0);
  c7->SetFillColor(0);
  c7->SetFillColor(0);
  c7->SetBorderMode(0);
  c7->SetFrameFillColor(0);
  c7->SetFrameBorderMode(0);
  hBkgd2D->Draw("colz");
  profBkgd2D->Draw("same");
  func->Draw("same");
  lineEDep2D->Draw();
  //funcCut->Draw("same");
  //funcAvgSig->Draw("same");
  //c7->SaveAs(Form("/tmp/Bkgd%i.png",bin));
  outname = "/tmp/TrackMatchingBkgd2D"+det+outnamebin+".png";
  c7->SaveAs(outname.Data());

  TCanvas *c8 = new TCanvas("c8","c8",600,400);
  c8->SetTopMargin(0.02);
  c8->SetRightMargin(0.02);
  c8->SetBorderSize(0);
  c8->SetFillColor(0);
  c8->SetFillColor(0);
  c8->SetBorderMode(0);
  c8->SetFrameFillColor(0);
  c8->SetFrameBorderMode(0);
  hSignal2D->Draw("colz");
  profSignal2D->Draw("same");
  func->Draw("same");
  lineEDep2D->Draw();
  //funcCut->Draw("same");
  //funcAvgSig->Draw("same");
  //c8->SaveAs(Form("/tmp/Signal%i.png",bin));
  outname = "/tmp/TrackMatchingSignal2D"+det+outnamebin+".png";
  c8->SaveAs(outname.Data());


  TCanvas *c9 = new TCanvas("c9","c9",600,400);
  c9->SetTopMargin(0.02);
  c9->SetRightMargin(0.02);
  c9->SetBorderSize(0);
  c9->SetFillColor(0);
  c9->SetFillColor(0);
  c9->SetBorderMode(0);
  c9->SetFrameFillColor(0);
  c9->SetFrameBorderMode(0);
  hSignalCutPeri->SetMaximum(1.0);
  hSignalCutPeri->Draw();
  hSignalAccPeri->Draw("same");
  hBkgdCutPeri->Draw("same");
  hBkgdAccPeri->Draw("same");
  lineEDep->Draw();
  legPeri->Draw();
  outname = "/tmp/TrackMatchingForCuts"+det+outnamebin+".png";
  c9->SaveAs(outname.Data());
  return;

  TCanvas *c5 = new TCanvas("c5","c5",600,400);
  c5->SetTopMargin(0.02);
  c5->SetRightMargin(0.02);
  c5->SetBorderSize(0);
  c5->SetFillColor(0);
  c5->SetFillColor(0);
  c5->SetBorderMode(0);
  c5->SetFrameFillColor(0);
  c5->SetFrameBorderMode(0);
  TH2F *hSignalOverBkgd = (TH2F*) fHistMatchedTracksEvspTSignal->Clone("SignalOverBackground");
  hSignalOverBkgd->Divide(fHistMatchedTracksEvspTBkgd);
  hSignalOverBkgd->Draw("colz");
  func->Draw("same");
  profSignal->Draw("esame");
  profBkgd->Draw("esame");

  TCanvas *c6 = new TCanvas("c6","c6",600,400);
  c6->SetTopMargin(0.02);
  c6->SetRightMargin(0.02);
  c6->SetBorderSize(0);
  c6->SetFillColor(0);
  c6->SetFillColor(0);
  c6->SetBorderMode(0);
  c6->SetFrameFillColor(0);
  c6->SetFrameBorderMode(0);
  TH2F *hBkgdOverSignal = (TH2F*) fHistMatchedTracksEvspTBkgd->Clone("BackgroundOverSignal");
  hBkgdOverSignal->Divide(fHistMatchedTracksEvspTSignal);
  hBkgdOverSignal->Draw("colz");
  func->Draw("same");
  profSignal->Draw("esame");
  profBkgd->Draw("esame");

  //c3->SaveAs("/tmp/SignalBkgdEmcal.png");

}

void SetStyles(TH1 *histo,int marker, int color,char *xtitle, char *ytitle){
  histo->SetMarkerStyle(marker);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  histo->GetXaxis()->SetTitle(xtitle);
  histo->GetYaxis()->SetTitle(ytitle);
}


TH1* bayneseffdiv(TH1* numerator, TH1* denominator,Char_t* name) 
{
    if(!numerator){
      cerr<<"Error:  numerator does not exist!"<<endl;
      return NULL;
    }
    if(!denominator){
      cerr<<"Error:  denominator does not exist!"<<endl;
      return NULL;
    }
    TH1* result = (TH1*)numerator->Clone(name);
    Int_t nbins = numerator->GetNbinsX();
    for (Int_t ibin=0; ibin<= nbins+1; ++ibin) {
      Double_t numeratorVal = numerator->GetBinContent(ibin);
      Double_t denominatorVal = denominator->GetBinContent(ibin);
      // Check if the errors are right or the thing is scaled
      Double_t numeratorValErr = numerator->GetBinError(ibin);
      if (!(numeratorValErr==0. || numeratorVal ==0.) ) {
	Double_t rescale = numeratorValErr*numeratorValErr/numeratorVal;
	numeratorVal /= rescale;
      }
      Double_t denominatorValErr = denominator->GetBinError(ibin);
      if (!(denominatorValErr==0. || denominatorVal==0. )) {
	Double_t rescale = denominatorValErr*denominatorValErr/denominatorVal;
	denominatorVal /= rescale;
      }
      Double_t quotient = 0.;
      if (denominatorVal!=0.) {
	quotient = numeratorVal/denominatorVal;
      }
      Double_t quotientErr=0;
      quotientErr = TMath::Sqrt(
				(numeratorVal+1.0)/(denominatorVal+2.0)*
				((numeratorVal+2.0)/(denominatorVal+3.0)-(numeratorVal+1.0)/(denominatorVal+2.0)));
      result->SetBinContent(ibin,quotient);
      result->SetBinError(ibin,quotientErr);
      //cout<<"Setting bin "<<ibin<<" to "<<quotient<<" "<<numeratorVal<<"/"<<denominatorVal<<endl;
    }
    return result;
}
