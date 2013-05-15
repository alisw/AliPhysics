void SetStyles(TH1 *histo,int marker, int color,char *xtitle = "p", char *ytitle = "E");
TFile *junkfile = NULL;

TFile *fNominal = NULL;
TFile *fLow = NULL;
TFile *fHigh = NULL;


TH1F* bayneseffdiv(TH1F* numerator, TH1F* denominator,Char_t* name);
TH1F *GetHisto(Bool_t isPHOS,int var){
//   TString myfilename = filename;
//   if(var==1){myfilename = filenameLow;}
//   if(var==2){myfilename = filenameHigh;}
  TString listname = "";
  TString histoName = "";
  TFile *f = NULL;
  switch(var){
  case 1:
    histoName = "Efficiency1";
    listname = "List1";
    f = fLow;
    break;
  case 2:
    histoName = "Efficiency2";
    listname = "List2";
    f = fHigh;
    break;
  default:
    histoName = "Efficiency0";
    listname = "List0";
    f = fNominal;
  }
  //cout<<"My file name "<<f->GetName()<<endl;
  //TFile *f = new TFile(filename.Data());  
  f->cd();
  //TList *l = (TList*)f->Get("out1");
  f->ls();
  //l->SetName(listname.Data());
  //cout<<"List name "<<l->GetName()<<endl;
  //TH1F  *fHistGammasGenerated = l->FindObject("fHistGammasGenerated")->Clone("GeneratedLHC11b1a");
  //TH1F  *fHistGammasFound = l->FindObject("fHistGammasFound")->Clone("FoundLHC11b1a");
  TH1F  *fHistGammasGeneratedLHC11b1a = (TH1F*)f->Get("fHistGammasGenerated");
  TH1F  *fHistGammasFoundLHC11b1a = (TH1F*)f->Get("fHistGammasFound");
  junkfile->cd();
  TH1F *hEfficiency = bayneseffdiv(fHistGammasFoundLHC11b1a,fHistGammasGeneratedLHC11b1a,histoName.Data());
//   delete fHistGammasGenerated;
//   delete fHistGammasFound;
//   delete l;
  hEfficiency->Rebin(2);
  hEfficiency->Scale(0.5);
  hEfficiency->SetName(histoName.Data());
  //cout<<"Histo name "<<hEfficiency->GetName()<<" entries "<<hEfficiency->GetEntries()<<" "<<fHistGammasGenerated->GetEntries()<<endl;
  hEfficiency->GetXaxis()->SetRange(1,hEfficiency->FindBin(3));
  //f->Close();
  //delete f;
  return hEfficiency;
}
void PlotGammaEfficiencyErrors(Bool_t isPHOS = kFALSE){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);


TString det = "EMCal";
if(isPHOS) det = "PHOS";
//  TString run = "";
//  switch(runnum){
//  case 1:
//    run = "118506";//bad file
//    break;
//  case 2:
//    run = "121040";
//    break;
//  case 3:
//    run = "118518";//low stats
//    break;
//  case 4:
//    run = "121039";//really low stats
//    break;
//  case 5:
//    run = "118558";//bad file
//    break;
//  }
// TString run = "118558";//bad file - adding files didn't work.
//TString run = "121040";
//TString run = "121039";
 TString outfilename = "junk.root";
 junkfile = new TFile(outfilename.Data(),"RECREATE");
 TString filename = "rootFiles/LHC11b1a/Efficiency."+det+".LHC11b1a.Combined.root";
 TString filenameHigh = "rootFiles/LHC11b1b/Efficiency."+det+".LHC11b1b.Combined.root";
 TString filenameLow = "rootFiles/LHC11b1c/Efficiency."+det+".LHC11b1c.Combined.root";
// TString filenameHigh = "rootFiles/LHC11b1b/Et.ESD.sim."+det+".LHC11b1b.Run"+run+".root";
// TString filename = "rootFiles/LHC11b1a/Et.ESD.sim."+det+".LHC11b1a.Run"+run+".root";
// TString filenameLow = "rootFiles/LHC11b1c/Et.ESD.sim."+det+".LHC11b1c.Run"+run+".root";

fNominal = new TFile(filename.Data());  
fLow = new TFile(filenameLow.Data());  
fHigh = new TFile(filenameHigh.Data());  


  hEfficiency = GetHisto(isPHOS,0);
  SetStyles(hEfficiency,20,1,"Energy","efficiency");
  hEfficiencyLow = GetHisto(isPHOS,1);
  SetStyles(hEfficiencyLow,23,2,"Energy","efficiency");
  hEfficiencyHigh = GetHisto(isPHOS,2);
  SetStyles(hEfficiencyHigh,22,4,"Energy","efficiency");

  junkfile->Write();
  //junkfile->Close();
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
  hEfficiencyLow->Draw("same");
  hEfficiencyHigh->Draw("same");
  TString pic1name = "/tmp/GammaEfficiency"+det+".eps";
  c->SaveAs(pic1name.Data());

  TH1F *hRatioLow = (TH1F*) hEfficiencyLow->Clone("hEfficiencyLowRatio");
  TH1F *hRatioHigh = (TH1F*) hEfficiencyHigh->Clone("hEfficiencyHighRatio");
  hRatioHigh->Divide(hEfficiency);
  hRatioLow->Divide(hEfficiency);

  TCanvas *c1 = new TCanvas("c1","c1",600,400);
  c1->SetTopMargin(0.02);
  c1->SetRightMargin(0.02);
  c1->SetBorderSize(0);
  c1->SetFillColor(0);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);
  hRatioLow->SetMaximum(1.2);
  hRatioLow->SetMinimum(0.8);
  hRatioLow->GetYaxis()->SetTitle("efficiency/nominal efficiency");
  TF1 *funcHigh = new TF1("funcHigh","[0]",0,3);
  TF1 *funcLow = new TF1("funcLow","[0]",0,3);
  funcHigh->SetLineColor(hRatioHigh->GetMarkerColor());
  funcLow->SetLineColor(hRatioLow->GetMarkerColor());
  TLegend *leg = new TLegend(0.144295,0.163978,0.234899,0.346774);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(hRatioLow,"Low material budget");
  leg->AddEntry(hRatioHigh,"High material budget");
  leg->SetTextSize(0.0537634);
  hRatioLow->Fit(funcLow);
  hRatioHigh->Fit(funcHigh);
  hRatioLow->Draw();
  hRatioHigh->Draw("same");
  leg->Draw();
  pic1name = "/tmp/GammaEfficiencyError"+det+".eps";
  c1->SaveAs(pic1name.Data());

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
    //result->SetName(name);
    return result;
}

void SetStyles(TH1 *histo,int marker, int color,char *xtitle, char *ytitle){
  histo->SetMarkerStyle(marker);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  histo->GetXaxis()->SetTitle(xtitle);
  histo->GetYaxis()->SetTitle(ytitle);
  histo->Sumw2();
}
