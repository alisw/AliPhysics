TH1* bayneseffdiv(TH1* numerator, TH1* denominator,Char_t* name);
TH2* bayneseffdiv2D(TH2* numerator, TH2* denominator,Char_t* name) ;
TH2F* flip(TH2F* orig,Char_t* name){
  float nbinx  = orig->GetXaxis()->GetNbins();
  float xlow = orig->GetXaxis()->GetBinLowEdge(1); 
  float xhigh = orig->GetXaxis()->GetBinLowEdge(nbinx+1); 
  float nbiny  = orig->GetYaxis()->GetNbins();
  float ylow= orig->GetYaxis()->GetBinLowEdge(1); 
  float yhigh = orig->GetYaxis()->GetBinLowEdge(nbiny+1); 
  TH2F *output = new TH2F(name,orig->GetTitle(),nbiny,ylow,yhigh,nbinx,xlow,xhigh);
  for(int i = 1;i<=nbinx;i++){
    for(int j = 1;j<=nbiny;j++){
      output->SetBinContent(j,i,orig->GetBinContent(i,j));
    }
  }
  return output;
}
void SetStyles(TH1 *histo,int marker, int color){
  histo->Sumw2();
  histo->SetMarkerStyle(marker);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  //histo->GetXaxis()->SetTitle(xtitle);
  //histo->GetYaxis()->SetTitle(ytitle);
}
Int_t colors[] = {0,TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack, 
		    TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack, 
		    TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack, 
		    TColor::kRed, TColor::kOrange, TColor::kGreen+3, TColor::kBlue, TColor::kBlack};
Int_t markers[] = {20,21,22,23,33, 24,25,26,32,27, 20,21,22,23,33, 24,25,26,32,27};

void PlotNeutronDeposits(TString filename = "rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.root"){
  gROOT->LoadMacro("macros/loadLibraries.C");
  loadLibraries();
   TFile *f = TFile::Open(filename, "READ");
    TList *l = dynamic_cast<TList*>(f->Get("out1"));
    TString histoname = "fHistNeutronsEtVsCent";
    TH2F *histo = l->FindObject(histoname.Data());
    histoname = "fHistNeutronsNumVsCent";
    TH2F *histoNum = l->FindObject(histoname.Data());
    TString outfilename = "/tmp/";
    TString calocorrfilename = "calocorrections.";
    TString det = "";
    if(filename.Contains("EMC")){
      outfilename +="EMCAL";
      calocorrfilename +="EMCAL";
      det += "Emcal";
    }
    else{
      outfilename +="PHOS";
      calocorrfilename +="PHOS";
      det += "Phos";
    }
    calocorrfilename +=".root";
    TString textfilename = "Neutrons"+det+".dat";
  TH2F  *fHistGammasGeneratedMult = l->FindObject("fHistGammasGeneratedCent");
  TH2F  *fHistGammasFoundMult = l->FindObject("fHistGammasFoundCent");
  TH2F  *fHistGammasFoundOutOfAccCent = l->FindObject("fHistGammasFoundOutOfAccCent");
  fHistGammasFoundMult->Add(fHistGammasFoundOutOfAccCent);
  // TH2F *gammaEff2D = (TH2F*) bayneseffdiv2D(fHistGammasFoundMult,fHistGammasGeneratedMult,"gammaEff2D");

    TCanvas *c1 = new TCanvas("c1","<E_{T}>",600,400);
    c1->SetTopMargin(0.02);
    c1->SetRightMargin(0.02);
    c1->SetBorderSize(0);
    c1->SetFillColor(0);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetFrameFillColor(0);
    c1->SetFrameBorderMode(0);
    TH1D *prof = histo->ProfileY();
    histo->Draw("colz");
    //prof->Draw("same");

    TCanvas *c3 = new TCanvas("c3","<E_{T}>",600,400);
    c3->SetTopMargin(0.02);
    c3->SetRightMargin(0.02);
    c3->SetBorderSize(0);
    c3->SetFillColor(0);
    c3->SetFillColor(0);
    c3->SetBorderMode(0);
    c3->SetFrameFillColor(0);
    c3->SetFrameBorderMode(0);
    //for(int i=1;i<=18;i++){
      
    //}
    prof->Draw();
    prof->GetYaxis()->SetTitle("<E_{T}>");
    prof->GetXaxis()->SetTitle("centrality bin");
    TString outfilename1 = outfilename+"Et.png";
    c3->SaveAs(outfilename1.Data());

    TCanvas *c2 = new TCanvas("c2","<Num>",600,400);
    c2->SetTopMargin(0.02);
    c2->SetRightMargin(0.02);
    c2->SetBorderSize(0);
    c2->SetFillColor(0);
    c2->SetFillColor(0);
    c2->SetBorderMode(0);
    c2->SetFrameFillColor(0);
    c2->SetFrameBorderMode(0);
    TH1D *profNum = histoNum->ProfileY();
    histoNum->Draw("colz");
    //profNum->Draw("same");

    TCanvas *c4 = new TCanvas("c4","<Num>",600,400);
    c4->SetTopMargin(0.02);
    c4->SetRightMargin(0.02);
    c4->SetBorderSize(0);
    c4->SetFillColor(0);
    c4->SetFillColor(0);
    c4->SetBorderMode(0);
    c4->SetFrameFillColor(0);
    c4->SetFrameBorderMode(0);
    profNum->Draw("same");
    profNum->GetYaxis()->SetTitle("<N_{neutrons}>");
    profNum->GetXaxis()->SetTitle("centrality bin");
    TString outfilename2 = outfilename+"Num.png";
    c4->SaveAs(outfilename2.Data());

    TFile *calocorrfile = TFile::Open(calocorrfilename, "READ");
    TString recoEffName = "ReCorrections";
    recoEffName += det;

    AliAnalysisEtRecEffCorrection *recoEff = (AliAnalysisEtRecEffCorrection *) calocorrfile->Get(recoEffName.Data());

    cout<<"I am here"<<endl;
    ofstream myfile;
    myfile.open (textfilename.Data());
    for(int i=0;i<20;i++){
      Float_t et = prof->GetBinContent(i+1);
      Float_t num = profNum->GetBinContent(i+1);
      Float_t eff = recoEff->ReconstructionEfficiency(et,i);
      float neutroncorr = 0;
      if(eff>0) neutroncorr = et*num/eff;
      float error = 0;
      cout<<"et*num/eff = "<<et<<"*"<<num<<"/"<<eff<<" "<< neutroncorr<<endl;
      myfile<<neutroncorr<<" "<<error<<endl;
    }
    myfile.close();
//     return;
//     TCanvas *c5 = new TCanvas("c5","<Num>",600,400);
//     c5->SetTopMargin(0.02);
//     c5->SetRightMargin(0.02);
//     c5->SetBorderSize(0);
//     c5->SetFillColor(0);
//     c5->SetFillColor(0);
//     c5->SetBorderMode(0);
//     c5->SetFrameFillColor(0);
//     c5->SetFrameBorderMode(0);


    TCanvas *c6 = new TCanvas("c6","<Num>",600,400);
    c6->SetTopMargin(0.02);
    c6->SetRightMargin(0.02);
    c6->SetBorderSize(0);
    c6->SetFillColor(0);
    c6->SetFillColor(0);
    c6->SetBorderMode(0);
    c6->SetFrameFillColor(0);
    c6->SetFrameBorderMode(0);
    c6->SetLogz();
    TH2F *histoNumAlt = flip(histoNum,"histoNumAlt");
    TH1D *profNumAlt = histoNum->ProfileY();
    histoNumAlt->Draw("colz");
    profNum->Draw("same");

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


TH2* bayneseffdiv2D(TH2* numerator, TH2* denominator,Char_t* name) 
{
  if(!numerator){
    cerr<<"Error:  numerator does not exist!"<<endl;
    return NULL;
  }
  if(!denominator){
    cerr<<"Error:  denominator does not exist!"<<endl;
    return NULL;
  }
  TH2* result = (TH2*)numerator->Clone(name);
  Int_t nbinsX = numerator->GetNbinsX();
  Int_t nbinsY = numerator->GetNbinsY();
  for (Int_t ibin=0; ibin<= nbinsX+1; ++ibin) {
    for (Int_t jbin=0; jbin<= nbinsY+1; ++jbin) {
      Double_t numeratorVal = numerator->GetBinContent(ibin,jbin);
      Double_t denominatorVal = denominator->GetBinContent(ibin,jbin);
      // Check if the errors are right or the thing is scaled
      Double_t numeratorValErr = numerator->GetBinError(ibin,jbin);
      if (!(numeratorValErr==0. || numeratorVal ==0.) ) {
	Double_t rescale = numeratorValErr*numeratorValErr/numeratorVal;
	numeratorVal /= rescale;
      }
      Double_t denominatorValErr = denominator->GetBinError(ibin,jbin);
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
      result->SetBinContent(ibin,jbin,quotient);
      result->SetBinError(ibin,jbin,quotientErr);
      //cout<<"Setting bin "<<ibin<<" to "<<quotient<<" "<<numeratorVal<<"/"<<denominatorVal<<endl;
    }
  }
  return result;
}
