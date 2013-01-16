void PlotGammaEfficiency(TString filename="rootFiles/LHC11a4_bis/Et.ESD.simPbPb.EMCAL.LHC11a4_bis.root"){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile *f = TFile::Open(filename, "READ");
  
  TList *l = (TList*)f->Get("out1");
  TH1F  *fHistGammasGenerated = l->FindObject("fHistGammasGenerated");
  TH1F  *fHistGammasFound = l->FindObject("fHistGammasFound");
  TH1F *hEfficiency = bayneseffdiv(fHistGammasFound,fHistGammasGenerated,"Efficiency");
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

}


TH1F* bayneseffdiv(TH1F* numerator, TH1F* denominator,Char_t* name) 
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
