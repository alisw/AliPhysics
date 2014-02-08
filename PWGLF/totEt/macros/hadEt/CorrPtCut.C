//Christine Nattrass, University of Tennessee at Knoxville
//This macro is for calculating the correction for the pT cut-off and its systematic error
//Uses the output of AliAnalysisTaskHadEt
//This is not actually what gets used in the correction class AliAnalysisHadEtCorrections - that is done in the macro GetCorrections.C - but this is useful for making plots and playing around with different options

float mean=0;
float highbound=0;
float lowbound=0;
float syserr = 0;

TH1D *GetHisto(float ptcut = 0.15, char *name, char *filename, float etacut){
  TFile *file = new TFile(filename);
  TList *list = file->FindObject("out2");
  //TH2F *allhad = ((TH2F*) out2->FindObject("EtSimulatedAllHadron"))->Clone("allhad");
  TH2F *allhad = ((TH2F*) out2->FindObject("EtSimulatedChargedHadron"))->Clone("allhad");
  TH2F *ptlow = ((TH2F*) out2->FindObject("EtSimulatedChargedHadronAssumingNoPt"))->Clone("ptlow");
  TH2F *pthigh;
  if(ptcut>0.14){//TPC cut off
    (TH2F*)pthigh =(TH2F*) ((TH2F*) out2->FindObject("EtSimulatedChargedHadronAssumingPtTPCCut"))->Clone("pthigh");
  }
  else{
    (TH2F*)pthigh =(TH2F*) ((TH2F*) out2->FindObject("EtSimulatedChargedHadronAssumingPtITSCut"))->Clone("pthigh");
  }

  int lowbin = allhad->GetXaxis()->FindBin(0.0);//make sure we don't accidentally get the wrong bin
  int highbin = allhad->GetXaxis()->FindBin(ptcut);
  int nbins = allhad->GetXaxis()->GetNbins();
  cout<<"Projecting from "<<allhad->GetXaxis()->GetBinLowEdge(lowbin)<<" to "<<allhad->GetXaxis()->GetBinLowEdge(highbin+1)<<endl;
  cout<<"Projecting from "<<allhad->GetXaxis()->GetBinLowEdge(lowbin)<<" to "<<allhad->GetXaxis()->GetBinLowEdge(nbins)<<endl;

  //allhad->Sumw2();


  TH1D *numerator = allhad->ProjectionY("name",lowbin,highbin);
  TH1D *denominator = allhad->ProjectionY("denominator",lowbin,nbins);
  TH1D *numeratorLow = ptlow->ProjectionY("nameLow",lowbin,highbin);
  TH1D *denominatorLow = allhad->ProjectionY("denominatorLow",highbin,nbins);
  denominatorLow->Add(ptlow);
  TH1D *numeratorHigh = pthigh->ProjectionY("nameHigh",lowbin,highbin);
  TH1D *denominatorHigh = allhad->ProjectionY("denominatorHigh",highbin,nbins);
  denominatorHigh->Add(pthigh);

  numerator->Divide(denominator);
  numeratorLow->Divide(denominatorLow);
  numeratorHigh->Divide(denominatorHigh);

  TF1 *funcLow = new TF1("funcLow","[0]",-.7,.7);
  funcLow->SetParameter(0,0.01);
  numeratorLow->Fit(funcLow);
  TF1 *func = new TF1("func","[0]",-.7,.7);
  func->SetParameter(0,0.02);
  numerator->Fit(func);
  TF1 *funcHigh = new TF1("funcHigh","[0]",-.7,.7);
  funcHigh->SetParameter(0,0.02);
  numeratorHigh->Fit(funcHigh);

  mean = 1.0-func->GetParameter(0);
  lowbound = 1.0-funcHigh->GetParameter(0);
  highbound = 1.0-funcLow->GetParameter(0);
  cout<<"fpTcut = "<<mean<<","<<lowbound<<","<<highbound<<endl;
  cout<<"1/fpTcut = "<<1.0/mean<<","<<1.0/lowbound<<","<<1.0/highbound<<endl;
  //cout<<"fpTcut = "<<mean<<"-"<<mean-lowbound<<"+"<<highbound-mean<<endl;
  syserr = highbound-mean;
  if(mean-lowbound>syserr) syserr = mean-lowbound;
  cout<<Form("%2.4f^{+%2.4f}_{-%2.4f}",mean,highbound-mean,mean-lowbound)<<endl;
  cout<<"latex here ";
  cout<<Form("%2.4f \\pm %2.4f",mean,syserr)<<endl;
  cout<<"1/fpTcut = "<<1.0/mean<<"+"<<1.0/lowbound-1.0/mean<<"-"<<1.0/mean-1.0/highbound<<endl;
  numerator->SetYTitle("E_{T}^{had, p_{T}<cut-off}/E_{T}^{had, all p_{T}}");
  numerator->GetYaxis()->SetTitleOffset(1.);
  numerator->GetYaxis()->SetTitleSize(0.08);
  numerator->GetYaxis()->SetLabelSize(0.05);
  numerator->GetXaxis()->SetTitleSize(0.08);
  numerator->GetXaxis()->SetLabelSize(0.05);
  numerator->GetXaxis()->SetTitleOffset(.6);
  //numerator->Rebin(2);
  //numerator->Scale(0.5);
  //numerator->Draw("e");
  return numerator;

}
void CorrPtCut(char *prodname = "LHC10d4 PYTHIA D6T 7 TeV p+p", char *shortprodname = "LHC10d4", char *filename="Et.ESD.new.sim.LHC10d4.pp.merged.root"){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TCanvas *c = new TCanvas("c","c",500,400);
  c->SetTopMargin(0.04);
  c->SetRightMargin(0.04);
  c->SetLeftMargin(0.181452);
  c->SetBottomMargin(0.134409);
  c->SetBorderSize(0);
  c->SetFillColor(0);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);
  float etacut = 0.7;
  cout<<"Pt cut = 150 MeV/c"<<endl;
  TH1D *High = GetHisto(0.15-.001,"High",filename,etacut);
  float tpcHigh=highbound;
  float tpcLow=lowbound;
  float tpcsyserr = syserr;
  float tpcmean = mean;
  float x1 = High->GetXaxis()->GetBinLowEdge(1);
  //TBox *tpcBox = new TBox(-x1*.99,1.0-tpcLow,x1*.99,1.0-tpcHigh);
  TBox *tpcBox = new TBox(-x1*.99,1.0-(mean-syserr),x1*.99,1.0-(mean+syserr));
  tpcBox->SetFillColor(5);
  tpcBox->SetLineColor(0);
  tpcBox->SetFillStyle(1001);
  cout<<"Pt cut = 100 MeV/c"<<endl;
  TH1D *Low = GetHisto(0.1-.001,"Low",filename,etacut);
  float itsHigh=highbound;
  float itsLow=lowbound;
  float itssyserr = syserr;
  float itsmean = mean;

  cout<<Form("dataset & %2.4f \\pm %2.4f &  %2.4f \\pm %2.4f \\",itsmean,itssyserr,tpcmean,tpcsyserr)<<endl;
  float x = Low->GetXaxis()->GetBinLowEdge(1);
  //TBox *itsBox = new TBox(-x*.99,1.0-itsLow,x*.99,1.0-itsHigh);
  TBox *itsBox = new TBox(-x1*.99,1.0-(mean-syserr),x1*.99,1.0-(mean+syserr));
  itsBox->SetFillColor(5);
  itsBox->SetLineColor(0);
  itsBox->SetFillStyle(1001);
  cout<<"Pt cut = 50 MeV/c"<<endl;
  TH1D *Lowest = GetHisto(0.05-.001,"Lowest",filename,etacut);
  TF1 *funcLow = new TF1("funcLow","[0]",-.7,.7);
  funcLow->SetParameter(0,0.01);
  Low->Fit(funcLow);
  TF1 *funcHigh = new TF1("funcHigh","[0]",-.7,.7);
  funcHigh->SetParameter(0,0.02);
  High->Fit(funcLow);
  High->SetMaximum(0.06);
  High->SetMinimum(0.0);
  High->SetMarkerColor(2);
  Low->SetMarkerColor(4);
  High->SetLineColor(2);
  Low->SetLineColor(4);
  High->SetMinimum(0.0);
  High->SetMarkerStyle(20);
  Low->SetMarkerStyle(21);
  Lowest->SetMarkerStyle(22);
  High->Draw();
  tpcBox->Draw("f");
  High->Draw("same");
  itsBox->Draw("f");
  //return;
  Low->Draw("same");
  //Lowest->Draw("same");
  TLatex *tex = new TLatex(-0.723444,0.0373593+0.019,prodname);
  tex->SetTextSize(0.0537634);
  tex->Draw();
  TLegend *leg = new TLegend(0.217742,0.696237,0.477823,0.873656);
  leg->AddEntry(High,"p_{T} cut-off = 0.15 GeV/c");
  leg->AddEntry(Low,"p_{T} cut-off = 0.1 GeV/c");
  //leg->AddEntry(Lowest,"p_{T} cut-off = 0.05 GeV/c");
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.0537634);
  leg->Draw();
  c->SaveAs(Form("pics/%s/fptcut.eps",shortprodname));
  c->SaveAs(Form("pics/%s/fptcut.png",shortprodname));
  c->SaveAs(Form("pics/%s/fptcut.pdf",shortprodname));
}
