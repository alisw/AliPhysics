TH1* bayneseffdiv(TH1* numerator, TH1* denominator,Char_t* name) ;

void SetStyles(TH1 *histo,int marker, int color,char *xtitle, char *ytitle, Bool_t scale = kFALSE);
void PlotTrackMatchingEfficiency(TString filename="rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run137366.root",TString run = "000000",TString prod = "LHC11a10a_bis"){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile *f = TFile::Open(filename, "READ");
  TList *l = (TList*)f->Get("out1");
  TH1F  *fHistHadronDepositsAll = l->FindObject("fHistHadronDepositsAll");
  TH1F  *fHistHadronDepositsReco = l->FindObject("fHistHadronDepositsReco");
  TH1F *eff = bayneseffdiv(fHistHadronDepositsReco,fHistHadronDepositsAll,"eff");
  TCanvas *c1 = new TCanvas("c1","c1",600,400);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  c1->SetBorderSize(0);
  c1->SetFillColor(0);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);
  eff->SetMinimum(0.0);
  eff->SetMaximum(1.0);
  eff->GetXaxis()->SetRange(1,eff->GetXaxis()->FindBin(3.0));
  eff->GetXaxis()->SetTitle("p_{T}");
  eff->GetYaxis()->SetTitle("track matching efficiency");
  eff->Draw();
  TString name = "/tmp/TrackMatchingEfficiency"+prod+run+".png";
  c1->SaveAs(name.Data());
}
void SetStyles(TH1 *histo,int marker, int color,char *xtitle, char *ytitle, Bool_t scale){
  histo->Sumw2();
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
