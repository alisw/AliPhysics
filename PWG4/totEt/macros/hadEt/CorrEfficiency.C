TH1D* bayneseffdiv(TH1D* numerator, TH1D* denominator,Char_t* name) 
{
    if(!numerator){
      cerr<<"Error:  numerator does not exist!"<<endl;
      return NULL;
    }
    if(!denominator){
      cerr<<"Error:  denominator does not exist!"<<endl;
      return NULL;
    }
    TH1D* result = (TH1D*)numerator->Clone(name);
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



TH1D *GetHisto(float cut = 0.12, char *name, int mycase, bool eta, int color, int marker,bool TPC, bool ITS){
  //TFile *file = new TFile("Et.ESD.new.sim.merged.root");
  TFile *file = new TFile("Et.ESD.new.sim.LHC10d4.pp.merged.root");
  TList *list = file->FindObject("out2");
  char *myname = "ITS";
  if(TPC&&!ITS) myname = "TPC";
  if(TPC&&ITS) myname = "TPCITS";
  cout<<"Using tracks from "<<myname<<" for efficiency"<<endl;
  TH2F *numeratorParent; 
  switch(mycase){
  case 0:
    numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"PiPlus")))->Clone("RecoHadron");
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"PiMinus")));
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"KMinus")));
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"KPlus")));
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"Proton")));
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"AntiProton")));
    break;
  case 1://pion
    numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"PiPlus")))->Clone("RecoPion");
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"PiMinus")));
    break;
  case 2://kaon
    numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"KPlus")))->Clone("RecoKaon");
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"KMinus")));
    break;
  case 3://proton
    numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"Proton")))->Clone("RecoProton");
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"AntiProton")));
    break;
  case 4://electron
    numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"EPlus")))->Clone("RecoElectron");
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"EMinus")));
    break;
  }
  TH2F *denominatorParent; 
  switch(mycase){
  case 0:
    denominatorParent = (TH2F*)((TH2F*) out2->FindObject("EtNSimulatedChargedHadron"))->Clone("RecoHadron");
//     denominatorParent = (TH2F*)((TH2F*) out2->FindObject("EtNSimulatedPiPlus"))->Clone("RecoHadron");
//     denominatorParent->Add((TH2F*) out2->FindObject("EtNSimulatedPiMinus"));
//     denominatorParent->Add((TH2F*) out2->FindObject("EtNSimulatedKMinus"));
//     denominatorParent->Add((TH2F*) out2->FindObject("EtNSimulatedKPlus"));
//     denominatorParent->Add((TH2F*) out2->FindObject("EtNSimulatedProton"));
//     denominatorParent->Add((TH2F*) out2->FindObject("EtNSimulatedAntiProton"));
    break;
  case 1://pion
    denominatorParent = (TH2F*)((TH2F*) out2->FindObject("EtNSimulatedPiPlus"))->Clone("RecoPion");
    denominatorParent->Add((TH2F*) out2->FindObject("EtNSimulatedPiMinus"));
    break;
  case 2://kaon
    denominatorParent = (TH2F*)((TH2F*) out2->FindObject("EtNSimulatedKPlus"))->Clone("RecoKaon");
    denominatorParent->Add((TH2F*) out2->FindObject("EtNSimulatedKMinus"));
    break;
  case 3://proton
    denominatorParent = (TH2F*)((TH2F*) out2->FindObject("EtNSimulatedProton"))->Clone("RecoProton");
    denominatorParent->Add((TH2F*) out2->FindObject("EtNSimulatedAntiProton"));
    break;
  case 4://electron
    denominatorParent = (TH2F*)((TH2F*) out2->FindObject("EtNSimulatedEPlus"))->Clone("RecoElectron");
    denominatorParent->Add((TH2F*) out2->FindObject("EtNSimulatedEMinus"));
    break;
  }
  numeratorParent->Sumw2();
  denominatorParent->Sumw2();
  TH1D *denominator;
  TH1D *numerator;
  if(eta){
    int lowbin = numeratorParent->GetYaxis()->FindBin(-cut+.001);//make sure we don't accv0entally get the wrong bin
    int highbin = numeratorParent->GetYaxis()->FindBin(cut-.001);
    cout<<"Projecting from "<<numeratorParent->GetYaxis()->GetBinLowEdge(lowbin)<<" to "<<numeratorParent->GetYaxis()->GetBinLowEdge(highbin+1)<<endl;
    denominator = denominatorParent->ProjectionX(Form("garbage%s",name),lowbin,highbin);
    numerator = numeratorParent->ProjectionX(name,lowbin,highbin);
  }
  else{
    int lowbin = denominatorParent->GetXaxis()->FindBin(cut);//make sure we don't accidentally get the wrong bin
    int highbin = denominatorParent->GetXaxis()->GetNbins();
    cout<<"Projecting from "<<denominatorParent->GetXaxis()->GetBinLowEdge(lowbin)<<" to "<<denominatorParent->GetXaxis()->GetBinLowEdge(highbin+1)<<endl;
    numerator = numeratorParent->ProjectionY(name,lowbin,highbin);
    denominator = denominatorParent->ProjectionY(Form("denominator%s",name),lowbin,highbin);
  }
  delete numeratorParent;
  delete denominatorParent;
  //numerator->Divide(denominator);
  TH1D *result = bayneseffdiv((TH1D*) numerator,(TH1D*)denominator,name);
  //result->Rebin(2);
  //result->Scale(0.5);
  result->SetYTitle("Efficiency");
  result->GetYaxis()->SetTitleOffset(0.8);
  result->GetXaxis()->SetTitleOffset(0.8);
  result->GetYaxis()->SetLabelSize(0.05);
  result->GetXaxis()->SetLabelSize(0.05);
  result->GetYaxis()->SetTitleSize(0.05);
  result->GetXaxis()->SetTitleSize(0.05);
  result->SetMarkerColor(color);
  result->SetLineColor(color);
  result->SetMarkerStyle(marker);
  //result->Draw("e");
  return result;

}

void CorrEfficiency(char *prodname= "LHC10d4", char *shortprodname = "LHC10d4 PYTHIA D6T 7 TeV p+p", bool TPC = true,bool ITS = true){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TCanvas *c = new TCanvas("c","c",600,400);
  c->SetTopMargin(0.02);
  c->SetRightMargin(0.02);
  c->SetBorderSize(0);
  c->SetFillColor(0);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);
  //c->SetLogx();

  int colortotal = 1;
  int colorpi = 2;
  int colork = 3;
  int colorp = 4;
  int phosmarker = 20;
  int emcalmarker = 24;
  float ptcut1 = 0.05;
  float ptcut2 = 0.1;
  TH1D *PHOStotal = GetHisto(0.12,"PHOStotal",0,true,colortotal,phosmarker,TPC,ITS);
  TH1D *PHOSpi = GetHisto(0.12,"PHOSpi",1,true,colorpi,phosmarker,TPC,ITS);
  TH1D *PHOSp = GetHisto(0.12,"PHOSp",2,true,colork,phosmarker,TPC,ITS);
  TH1D *PHOSk = GetHisto(0.12,"PHOSk",3,true,colorp,phosmarker,TPC,ITS);
  if(ITS&&!TPC){PHOStotal->GetXaxis()->SetRange(PHOStotal->GetXaxis()->FindBin(0.05),PHOStotal->GetXaxis()->FindBin(1.0));}
  else{PHOStotal->GetXaxis()->SetRange(PHOStotal->GetXaxis()->FindBin(0.0),PHOStotal->GetXaxis()->FindBin(3.0));}
  PHOStotal->SetMinimum(0.0);
  PHOStotal->SetMaximum(1.0);
  //parameters[centbin][0]*exp(-pow(parameters[centbin][1]/pt,parameters[centbin][2]))
  TF1 *func = new TF1("func","[0]*exp(-pow([1]/x,[2]))");
  func->SetParameter(0,.9);
  func->SetParameter(1,.05);
  func->SetParLimits(1,1e-3,1);
  func->SetParameter(2,.1);
  //PHOStotal->Fit(func);
  PHOStotal->Draw();
  PHOSpi->Draw("same");
  PHOSp->Draw("same");
  PHOSk->Draw("same");
  TH1D *EMCALtotal = GetHisto(0.7,"EMCALtotal",0,true,colortotal,emcalmarker,TPC,ITS);
  TH1D *EMCALpi = GetHisto(0.7,"EMCALpi",1,true,colorpi,emcalmarker,TPC,ITS);
  TH1D *EMCALp = GetHisto(0.7,"EMCALp",2,true,colork,emcalmarker,TPC,ITS);
  TH1D *EMCALk = GetHisto(0.7,"EMCALk",3,true,colorp,emcalmarker,TPC,ITS);
  EMCALtotal->Draw("same");
  EMCALpi->Draw("same");
  EMCALp->Draw("same");
  EMCALk->Draw("same");


  TLegend *leg = new  TLegend(0.22651,0.247312,0.370805,0.438172);
  leg->AddEntry(PHOStotal,"#pi,K,p");
  leg->AddEntry(PHOSpi,"#pi^{#pm}");
  leg->AddEntry(PHOSk,"K^{#pm}");
  leg->AddEntry(PHOSp,"p,#bar{p}");
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.06);
 leg->Draw();

  TLine *line = new TLine(0.2,0.0,0.2,1.0);
  line->Draw();
  line->SetLineWidth(3.0);
  //line->SetLineColor(TColor::kYellow);
  line->SetLineStyle(2);
  TLatex *tex = new TLatex(0.497269,0.0513196,prodname);
  tex->SetTextSize(0.0537634);
  tex->Draw();
  TLatex *tex3 = new TLatex(1.16186,0.28348,"Closed symbols |#eta|<0.12 (PHOS)");
  tex3->SetTextSize(0.0537634);
  tex3->Draw();
  TLatex *tex4 = new TLatex(1.16186,0.213221,"Open symbols |#eta|<0.70 (EMCal)");
  tex4->SetTextSize(0.0537634);
  tex4->Draw();
  TLatex *tex2 = new TLatex(0.241937,0.448436,"TPC cut-off 150 MeV/c");
  tex2->SetTextSize(0.0537634);
  tex2->Draw();
  if(TPC){
    c->SaveAs("pics/CorrEfficiency.eps");
    c->SaveAs("pics/CorrEfficiency.png");
  }
  else{
    c->SaveAs("pics/CorrEfficiencyTPCITS.eps");
    c->SaveAs("pics/CorrEfficiencyTPCITS.png");
  }
}
