void SetStyles(TH1 *histo,int marker, int color,char *xtitle, char *ytitle, Bool_t scale = kFALSE);

void PlotGammaEfficiency(Bool_t isPhos = kFALSE, Int_t cutset = 0){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TString det = "";
    if(isPhos){
      det = "PHOS";
    }
    else{
      det = "EMCal";
    }
    if(cutset==0){
      if(isPhos){
	filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.PHOS.LHC11a10a_bis.Run139465.root";
      }
      else{
	filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.root";
      }
    }
    if(cutset==1){
      if(isPhos){
	filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.PHOSOldHadMethod.LHC11a10a_bis.root";
      }
      else{
	filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCalOldHadMethod.LHC11a10a_bis.root";
      }
    }
    if(cutset==2){
      if(isPhos){
	filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.PHOS.LHC11a10a_bis.root";
      }
      else{
	filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.root";
      }
    }



  TFile *f = TFile::Open(filename, "READ");
  
  TList *l = (TList*)f->Get("out1");
  TH1F  *fHistGammasGenerated = l->FindObject("fHistGammasGenerated");
  TH1F  *fHistGammasFound = l->FindObject("fHistGammasFound");
  TH1F *hEfficiency = bayneseffdiv(fHistGammasFound,fHistGammasGenerated,"Efficiency");
  TH2F  *fHistGammasGeneratedMult = l->FindObject("fHistGammasGeneratedCent");
  //TH2F  *fHistGammasFoundMult = l->FindObject("fHistGammasFoundCent");
  //TH2F  *fHistGammasFoundOutOfAccMult = l->FindObject("fHistGammasFoundOutOfAccCent");
  //    TH2F *fHistGammasFoundRecoEnergyCent;
  //  TH2F *fHistGammasFoundOutOfAccRecoEnergyCent;
  TH2F  *fHistGammasFoundMult = l->FindObject("fHistGammasFoundRecoEnergyCent");
  TH2F  *fHistGammasFoundOutOfAccMult = l->FindObject("fHistGammasFoundOutOfAccRecoEnergyCent");
   fHistGammasFoundMult->Add(fHistGammasFoundOutOfAccMult);
  hEfficiency->GetXaxis()->SetTitle("Energy");
  hEfficiency->GetYaxis()->SetTitle("efficiency");
  hEfficiency->GetXaxis()->SetRange(1,hEfficiency->FindBin(3));

  TCanvas *c = new TCanvas("c","c",600,400);
  c->SetTopMargin(0.0187166);
  c->SetRightMargin(0.0201342);
  c->SetLeftMargin(0.112416);
  c->SetBottomMargin(0.15508);
  c->SetBorderSize(0);
  c->SetFillColor(0);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);

  hEfficiency->Draw();
  c->SaveAs("/tmp/GammaEfficiency.png");

  TCanvas *c1 = new TCanvas("c1","c1",600,400);
  c1->SetTopMargin(0.0187166);
  c1->SetRightMargin(0.0201342);
  c1->SetLeftMargin(0.112416);
  c1->SetBottomMargin(0.15508);
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
  TLegend *leg1 = new TLegend(0.129195,0.692513,0.29698,0.970588);
  leg1->SetFillStyle(0);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetTextSize(0.0534759);
  Int_t binwidth = 4;
  for(int i=0;i<20;i+=binwidth){
    Int_t maxbin = i+1+binwidth;
    fHistGammasGeneratedCent[i] = fHistGammasGeneratedMult->ProjectionX(Form("All%i",i),i+1,maxbin);
    fHistGammasFoundCent[i] = fHistGammasFoundMult->ProjectionX(Form("Reco%i",i),i+1,maxbin);
    fEff[i] = (TH1D*) bayneseffdiv(fHistGammasFoundCent[i],fHistGammasGeneratedCent[i],Form("eff%i",i));
    SetStyles(fEff[i],markers[i],colors[i],"energy","efficiency",1.0);
    TString legentry = Form("%i-%i",i*5,(maxbin-1)*5);
    legentry +="%";
    leg1->AddEntry(fEff[i],legentry.Data());
    if(i==0){ 
      fEff[i]->SetMaximum(1.5);
      fEff[i]->SetMinimum(0.0);
      fEff[i]->GetXaxis()->SetRange(1,fEff[i]->FindBin(3.0));
      fEff[i]->GetXaxis()->SetTitleSize(0.07);
      fEff[i]->GetYaxis()->SetTitleSize(0.07);
      fEff[i]->GetYaxis()->SetTitleOffset(0.7);
      fEff[i]->GetXaxis()->SetLabelSize(0.07);
      fEff[i]->GetYaxis()->SetLabelSize(0.07);
      fEff[i]->Draw();
    }
    else{ fEff[i]->Draw("same");}

  }
  leg1->Draw();
  TString outname = "/tmp/"+det+"GammaEfficiencyCent.eps";
    c1->SaveAs(outname.Data());
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
      if(quotient>1) quotientErr = 0;
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
