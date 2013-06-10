void SetStyles(TH1 *histo,int marker, int color,char *xtitle, char *ytitle, Bool_t scale = kFALSE);

void PlotGammaEfficiency(TString filename="rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.root"){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile *f = TFile::Open(filename, "READ");
  
  TList *l = (TList*)f->Get("out1");
  TH1F  *fHistGammasGenerated = l->FindObject("fHistGammasGenerated");
  TH1F  *fHistGammasFound = l->FindObject("fHistGammasFound");
  TH1F *hEfficiency = bayneseffdiv(fHistGammasFound,fHistGammasGenerated,"Efficiency");
  TH2F  *fHistGammasGeneratedMult = l->FindObject("fHistGammasGeneratedCent");
  TH2F  *fHistGammasFoundMult = l->FindObject("fHistGammasFoundCent");
  hEfficiency->GetXaxis()->SetTitle("Energy");
  hEfficiency->GetYaxis()->SetTitle("efficiency");
  hEfficiency->GetXaxis()->SetRange(1,hEfficiency->FindBin(3));

  TCanvas *c = new TCanvas("c","c",600,400);
  c->SetTopMargin(0.02);
  c->SetRightMargin(0.02);
  c->SetBorderSize(0);
  c->SetFillColor(0);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);

  hEfficiency->Draw();
  c->SaveAs("/tmp/GammaEfficiency.png");

  TCanvas *c1 = new TCanvas("c1","c1",600,400);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  c1->SetBorderSize(0);
  c1->SetFillColor(0);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);
  TH1D *fHistGammasGeneratedCent[22] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL, NULL};
  TH1D *fHistGammasFoundCent[22] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL, NULL};
  TH1D *fEff[22] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL, NULL};
  int colors[] = {TColor::kRed,TColor::kRed, TColor::kOrange-3, TColor::kOrange-3, TColor::kGreen+3,
                  TColor::kGreen+3, TColor::kBlue, TColor::kBlue, TColor::kViolet, TColor::kViolet,
                  TColor::kMagenta+3,TColor::kRed,TColor::kRed, TColor::kOrange-3, TColor::kOrange-3, TColor::kGreen+3,
                  TColor::kGreen+3, TColor::kBlue, TColor::kBlue, TColor::kViolet, TColor::kViolet,
                  TColor::kMagenta+3};
  int markers[] = {20,24,21,25,22, 26,23,32,33,27, 29,20,24,21,25,22, 26,23,32,33,27, 29};
  TLegend *leg1 = new TLegend(0.224832,0.153226,0.39094,0.647849);
  leg1->SetFillStyle(0);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetTextSize(0.0456989);
  for(int i=0;i<20;i+=4){
    fHistGammasGeneratedCent[i] = fHistGammasGeneratedMult->ProjectionX(Form("All%i",i),i+1,i+2);
    fHistGammasFoundCent[i] = fHistGammasFoundMult->ProjectionX(Form("Reco%i",i),i+1,i+2);
    fEff[i] = (TH1D*) bayneseffdiv(fHistGammasFoundCent[i],fHistGammasGeneratedCent[i],Form("eff%i",i));
    SetStyles(fEff[i],markers[i],colors[i],"energy","efficiency",1.0);
    leg1->AddEntry(fEff[i],Form("centbin %i",i));
    if(i==0){ 
      fEff[i]->SetMaximum(0.8);
      fEff[i]->GetXaxis()->SetRange(1,fEff[i]->FindBin(3.0));
      fEff[i]->Draw();
    }
    else{ fEff[i]->Draw("same");}

  }
  leg1->Draw();
//   fHistGammasFound->Draw();

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
    TH1F* result = (TH1F*)numerator->Clone(name);
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
void SetStyles(TH1 *histo,int marker, int color,char *xtitle, char *ytitle, Bool_t scale){
  histo->Sumw2();
  histo->SetMarkerStyle(marker);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  histo->GetXaxis()->SetTitle(xtitle);
  histo->GetYaxis()->SetTitle(ytitle);
}
