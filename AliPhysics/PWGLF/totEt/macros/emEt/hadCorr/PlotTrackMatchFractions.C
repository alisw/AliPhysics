void SetStyles(TH1 *histo,int marker, int color,char *xtitle, char *ytitle);
TH1* bayneseffdiv(TH1* numerator, TH1* denominator,Char_t* name) ;
void PlotTrackMatchFractions(TString filename="rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.root", TString det = "EMCal"){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TString outname = "/tmp/TrackMatchFractionsGammas"+det+".png";
  TString outname2 = "/tmp/TrackMatchFractionsHadrons"+det+".png";
  //TString outnamebin = Form("%iTo%i",bin,binLast);

  TFile *f = TFile::Open(filename, "READ");
  
  TList *l = (TList*)f->Get("out1");
  TH1F *fHistGammasCut = l->FindObject("fHistGammasTracksCut");
  TH1F *fHistGammasAccepted = l->FindObject("fHistGammasTracksAccepted");
  TH1F *fHistChargedTracksCut = l->FindObject("fHistChargedTracksCut");
  TH1F *fHistChargedTracksAccepted = l->FindObject("fHistChargedTracksAccepted");
  TH1F *fHistGammasTotal = fHistGammasCut->Clone("fHistGammasTotal");
  fHistGammasTotal->Add(fHistGammasAccepted);
  TH1F *fHistGammasFracCut = bayneseffdiv(fHistGammasCut,fHistGammasTotal,"fHistGammasFracCut");
  TH1F *fHistGammasFracAcc = bayneseffdiv(fHistGammasAccepted,fHistGammasTotal,"fHistGammasFracCut");
  SetStyles(fHistGammasFracCut,20,2,"E_{T}","Fraction");
  SetStyles(fHistGammasFracAcc,25,4,"E_{T}","Fraction");
  TH1F *fHistChargedTracksTotal = fHistChargedTracksCut->Clone("fHistChargedTracksTotal");
  fHistChargedTracksTotal->Add(fHistChargedTracksAccepted);
  TH1F *fHistChargedTracksFracCut = bayneseffdiv(fHistChargedTracksCut,fHistChargedTracksTotal,"fHistChargedTracksFracCut");
  TH1F *fHistChargedTracksFracAcc = bayneseffdiv(fHistChargedTracksAccepted,fHistChargedTracksTotal,"fHistChargedTracksFracCut");
  SetStyles(fHistChargedTracksFracCut,20,2,"E_{T}","Fraction");
  SetStyles(fHistChargedTracksFracAcc,25,4,"E_{T}","Fraction");
  TLegend *leg = new TLegend(0.357383,0.47043,0.458054,0.75);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(fHistGammasFracCut,"Fraction of #gamma cut");
  leg->AddEntry(fHistGammasFracAcc,"Fraction of #gamma accepted");
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
  fHistGammasFracCut->SetMaximum(1.0);
  fHistGammasFracCut->SetMinimum(0.0);
  fHistGammasFracCut->Draw();
  fHistGammasFracAcc->Draw("same");
  leg->Draw();
  c1->SaveAs(outname.Data());

  TLegend *leg2 = new TLegend(0.357383,0.47043,0.458054,0.75);
  leg2->SetFillStyle(0);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.03);
  leg2->AddEntry(fHistGammasFracCut,"Fraction of hadrons cut");
  leg2->AddEntry(fHistGammasFracAcc,"Fraction of hadrons accepted");
  leg2->SetTextSize(0.061828);
  TCanvas *c2 = new TCanvas("c2","c2",600,400);
  c2->SetTopMargin(0.02);
  c2->SetRightMargin(0.02);
  c2->SetBorderSize(0);
  c2->SetFillColor(0);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetFrameFillColor(0);
  c2->SetFrameBorderMode(0);
  fHistChargedTracksFracCut->SetMaximum(1.0);
  fHistChargedTracksFracCut->SetMinimum(0.0);
  fHistChargedTracksFracCut->Draw();
  fHistChargedTracksFracAcc->Draw("same");
  leg2->Draw();

  c2->SaveAs(outname2.Data());

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
