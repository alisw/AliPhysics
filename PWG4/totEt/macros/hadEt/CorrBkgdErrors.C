TH1D *GetHisto(float etacut = 0.12, bool eta = true, char *name, bool TPC,bool ITS, int mycase = 0, int color=1, int marker = 20, char *filename="Et.ESD.new.sim.merged.root",bool reweight = false,float kaonFactor=1.0, float lambdaFactor = 1.0, float baryonEnhancement = 1.0){
  TFile *file = new TFile(filename);
  TList *list = file->FindObject("out2");
  char *reweightname = "";
  if(reweight) reweightname = "Reweighted";
  char *myname = "ITS";
  if(TPC){
    if(ITS) myname = "TPCITS";
    else{    myname = "TPC";}
  }
  TH2F *signal = ((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedPiPlus",myname)))->Clone("signal");
  signal->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedPiMinus",myname)));
  signal->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedKMinus",myname)),kaonFactor);
  signal->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedKPlus",myname)),kaonFactor);
  signal->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedProton",myname)),baryonEnhancement);
  signal->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedAntiProton",myname)),baryonEnhancement);
  signal->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sUnidentifiedAssumingPion",myname)));

  //Et of all unidentified hadrons (plus hadrons identified as pions) calculated assuming their true mass
  TH2F *bkgd;
  switch(mycase){
  case 0:
    bkgd = (TH2F*)((TH2F*) out2->FindObject(Form("EtReconstructed%sMisidentifiedElectrons",myname)))->Clone("bkgd");
    bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sLambdaDaughters%s",myname,reweightname)),baryonEnhancement*lambdaFactor);
    bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sAntiLambdaDaughters%s",myname,reweightname)),baryonEnhancement*lambdaFactor);
    bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sK0SDaughters%s",myname,reweightname)),kaonFactor);
    bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sXiDaughters",myname)),baryonEnhancement);
    bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sAntiXiDaughters",myname)),baryonEnhancement);
    bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sOmegaDaughters",myname)),baryonEnhancement);
    bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sAntiOmegaDaughters",myname)),baryonEnhancement);
    break;
  case 1:
    bkgd = (TH2F*)((TH2F*) out2->FindObject(Form("EtReconstructed%sMisidentifiedElectrons",myname)))->Clone("bkgd");
    break;
  case 2:
    bkgd = (TH2F*)((TH2F*) out2->FindObject(Form("EtReconstructed%sLambdaDaughters%s",myname,reweightname)))->Clone("bkgd");
    bkgd->Scale(baryonEnhancement*lambdaFactor);
    bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sAntiLambdaDaughters%s",myname,reweightname)),baryonEnhancement*lambdaFactor);
    break;
  case 3:
    bkgd = (TH2F*)((TH2F*) out2->FindObject(Form("EtReconstructed%sK0SDaughters%s",myname,reweightname)))->Clone("bkgd");
    bkgd->Scale(kaonFactor);
    break;
  case 4:
    bkgd = (TH2F*)((TH2F*) out2->FindObject(Form("EtReconstructed%sLambdaDaughters%s",myname,reweightname)))->Clone("bkgd");
    bkgd->Scale(baryonEnhancement*lambdaFactor);
    bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sAntiLambdaDaughters%s",myname,reweightname)),baryonEnhancement*lambdaFactor);
    bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sK0SDaughters%s",myname,reweightname)),kaonFactor);
    break;
  case 5:
    bkgd = (TH2F*)((TH2F*) out2->FindObject(Form("EtReconstructed%sXiDaughters",myname)))->Clone("bkgd");
    bkgd->Scale(baryonEnhancement);
    bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sAntiXiDaughters",myname)),baryonEnhancement);
    break;
  case 6:
    bkgd = (TH2F*)((TH2F*) out2->FindObject(Form("EtReconstructed%sOmegaDaughters",myname)))->Clone("bkgd");
    bkgd->Scale(baryonEnhancement);
    bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sAntiOmegaDaughters",myname)),baryonEnhancement);
    break;
  case 7:
    bkgd = (TH2F*)((TH2F*) out2->FindObject(Form("EtReconstructed%sXiDaughters",myname)))->Clone("bkgd");
    bkgd->Scale(baryonEnhancement);
    bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sAntiXiDaughters",myname)),baryonEnhancement);
    bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sOmegaDaughters",myname)),baryonEnhancement);
    bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sAntiOmegaDaughters",myname)),baryonEnhancement);
    break;
  case 8:
    bkgd = (TH2F*)((TH2F*) out2->FindObject(Form("EtReconstructed%sLambdaDaughters%s",myname,reweightname)))->Clone("bkgd");
    bkgd->Scale(baryonEnhancement*lambdaFactor);
    break;
  case 9:
    bkgd = (TH2F*)((TH2F*) out2->FindObject(Form("EtReconstructed%sAntiLambdaDaughters%s",myname,reweightname)))->Clone("bkgd");
    bkgd->Scale(baryonEnhancement*lambdaFactor);
    break;
  }
  TH1D *denominator;
  TH1D *numerator;
  if(eta){
    int lowbin = bkgd->GetXaxis()->FindBin(etacut+.001);//make sure we don't accidentally get the wrong bin
    int highbin = bkgd->GetXaxis()->GetNbins();
    cout<<"Projecting from "<<bkgd->GetXaxis()->GetBinLowEdge(lowbin)<<" to "<<bkgd->GetXaxis()->GetBinLowEdge(highbin+1)<<endl;
    denominator = signal->ProjectionY("name",lowbin,highbin);
    numerator = bkgd->ProjectionY(name,lowbin,highbin);
  }
  else{
    int lowbin = bkgd->GetYaxis()->FindBin(-etacut+.001);//make sure we don't accidentally get the wrong bin
    int highbin = bkgd->GetYaxis()->FindBin(etacut-.001);
    cout<<"Projecting from "<<bkgd->GetYaxis()->GetBinLowEdge(lowbin)<<" to "<<bkgd->GetYaxis()->GetBinLowEdge(highbin+1)<<endl;
    denominator = signal->ProjectionX("name",lowbin,highbin);
    numerator = bkgd->ProjectionX(name,lowbin,highbin);
  }
  numerator->Divide(denominator);
  numerator->SetYTitle("Ratio of E_{T}^{background}/E_{T}^{had, meas.}");
  numerator->GetYaxis()->SetTitleOffset(1.2);
  numerator->SetLineColor(color);
  numerator->SetMarkerColor(color);
  numerator->SetMarkerStyle(marker);
  return numerator;

}

void CorrBkgdErrors(bool TPC = true, bool ITS=true, bool reweight = true,float kaonFactor=1.0, float lambdaFactor = 1.0, float baryonEnhancement = 1.0){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TCanvas *c = new TCanvas("c","c",400,400);
  c->SetTopMargin(0.04);
  c->SetRightMargin(0.04);
  c->SetBorderSize(0);
  c->SetFillColor(0);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);
  c->SetLeftMargin(0.159091);
  char * filename = "Et.ESD.new.sim.LHC10d4.pp.merged.root";
  TH1D *All = GetHisto(0.1,true,"All",TPC,ITS,0,1,20,filename,reweight,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *Electrons = GetHisto(0.1,true,"Electrons",TPC,ITS,1,2,21,filename,reweight,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *V0s = GetHisto(0.1,true,"V0s",TPC,ITS,4,4,22,filename,reweight,kaonFactor,lambdaFactor,baryonEnhancement);
  //case, color, marker
  TH1D *K0S = GetHisto(0.1,true,"K0S",TPC,ITS,3,TColor::kOrange+8,33,filename,reweight,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *Lambda = GetHisto(0.1,true,"K0S",TPC,ITS,8,TColor::kMagenta+3,29,filename,reweight,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *AntiLambda = GetHisto(0.1,true,"K0S",TPC,ITS,9,TColor::kMagenta+3,30,filename,reweight,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *Multistrange = GetHisto(0.1,true,"Multistrange",TPC,ITS,7,TColor::kGreen+2,23,filename,reweight,kaonFactor,lambdaFactor,baryonEnhancement);

  TH1D *Allpt = GetHisto(0.7,false,"Allpt",TPC,ITS,0,1,20,filename,reweight,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *Electronspt = GetHisto(0.7,false,"Electronspt",TPC,ITS,1,2,21,filename,reweight,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *V0spt = GetHisto(0.7,false,"V0spt",TPC,ITS,4,4,22,filename,reweight,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *K0Spt = GetHisto(0.1,false,"K0S",TPC,ITS,3,TColor::kOrange+8,33,filename,reweight,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *Lambdapt = GetHisto(0.1,false,"K0S",TPC,ITS,8,TColor::kMagenta+3,29,filename,reweight,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *AntiLambdapt = GetHisto(0.1,false,"K0S",TPC,ITS,9,TColor::kMagenta+3,30,filename,reweight,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *Multistrangept = GetHisto(0.7,false,"Multistrangept",TPC,ITS,7,TColor::kGreen+2,23,filename,reweight,kaonFactor,lambdaFactor,baryonEnhancement);
  if(ITS){
    All->SetMaximum(0.025);
  }
  else{
    All->SetMaximum(0.2);
  }
  All->SetMinimum(0.0);
  All->SetMarkerStyle(20);
  All->GetYaxis()->SetTitleOffset(1.8);
  TF1 *func = new TF1("func","[0]",-.7,.7);
  func->SetParameter(0,0.02);
  TF1 *funcLam = new TF1("funcLam","[0]",-.7,.7);
  funcLam->SetParameter(0,0.001);
  funcLam->SetLineColor(Lambda->GetMarkerColor());
  TF1 *funcAlam = new TF1("funcAlam","[0]",-.7,.7);
  funcAlam->SetParameter(0,0.003);
  funcAlam->SetLineColor(AntiLambda->GetMarkerColor());
  TF1 *funcK0 = new TF1("funcK0","[0]",-.7,.7);
  funcK0->SetParameter(0,0.013);
  funcK0->SetLineColor(K0S->GetMarkerColor());
  All->Fit(func);
  Lambda->Fit(funcLam);
  AntiLambda->Fit(funcAlam);
  K0S->Fit(funcK0);
  All->Draw();
  Electrons->Draw("same");
  V0s->Draw("same");
  Multistrange->Draw("same");
  K0S->Draw("same");
  Lambda->Draw("same");
  AntiLambda->Draw("same");
  TLatex *tex = new TLatex(0.161478,1.0835,"LHC10d15: p+p, Pythia6 Perugia-0");
  tex->SetTextSize(0.0537634);
  tex->Draw();
  TLegend *leg = new TLegend(0.636364,0.612903,0.893939,0.962366);
  leg->AddEntry(All,"All");
  leg->AddEntry(Electrons,"Electrons");
  leg->AddEntry(V0s,"V0s");
  leg->AddEntry(K0S,"K_{S}^{0}");
  leg->AddEntry(Lambda,"#Lambda");
  leg->AddEntry(AntiLambda,"#bar{#Lambda}");
  leg->AddEntry(Multistrange,"Multistrange");
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();
  TLatex *tex = new TLatex(-0.711139,0.0157696,Form("%2.5f#pm%2.5f",func->GetParameter(0),func->GetParError(0)));
  tex->Draw();
  TLatex *texLam = new TLatex(-0.711139,0.00201613,Form("%2.5f#pm%2.5f",funcLam->GetParameter(0),funcLam->GetParError(0)));
  texLam->SetTextColor(Lambda->GetMarkerColor());
  texLam->Draw();
  TLatex *texAlam = new TLatex(-0.711139,0.00365716,Form("%2.5f#pm%2.5f",funcAlam->GetParameter(0),funcAlam->GetParError(0)));
  texAlam->SetTextColor(AntiLambda->GetMarkerColor());
  texAlam->Draw();
  TLatex *texK0 = new TLatex(-0.711139,0.008,Form("%2.5f#pm%2.5f",funcK0->GetParameter(0),funcK0->GetParError(0)));
  texK0->SetTextColor(K0S->GetMarkerColor());
  texK0->Draw();
  TCanvas *c1 = new TCanvas("c1","c1",600,400);
  c1->SetTopMargin(0.04);
  c1->SetRightMargin(0.04);
  c1->SetBorderSize(0);
  c1->SetFillColor(0);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);
  if(ITS){
    Allpt->SetMaximum(0.04);
  }
  else{
    Allpt->SetMaximum(0.2);
  }
  if(TPC)  Allpt->GetXaxis()->SetRange(Allpt->GetXaxis()->FindBin(0.0),Allpt->GetXaxis()->FindBin(4.));
  else{  Allpt->GetXaxis()->SetRange(Allpt->GetXaxis()->FindBin(0.0),Allpt->GetXaxis()->FindBin(1.));}
  Allpt->SetMinimum(0.0);
  Allpt->SetMarkerStyle(20);
  Allpt->Draw();
  Electronspt->Draw("same");
  V0spt->Draw("same");
  K0Spt->Draw("same");
  Lambdapt->Draw("same");
  AntiLambdapt->Draw("same");
  Multistrangept->Draw("same");
  TLatex *texpt = new TLatex(0.161478,1.0835,"LHC10d15: p+p, Pythia6 Perugia-0");
  texpt->SetTextSize(0.0537634);
  texpt->Draw();
  TLegend *legpt = new TLegend(0.634228,0.430108,0.892617,0.905914);
  legpt->AddEntry(Allpt,"All");
  legpt->AddEntry(Electronspt,"Electrons");
  legpt->AddEntry(V0spt,"V0s");
  legpt->AddEntry(K0Spt,"K_{S}^{0}");
  legpt->AddEntry(Lambdapt,"#Lambda");
  legpt->AddEntry(AntiLambdapt,"#bar{#Lambda}");
  legpt->AddEntry(Multistrangept,"Multistrange");
  legpt->SetFillStyle(0);
  legpt->SetFillColor(0);
  legpt->SetBorderSize(0);
  legpt->Draw();


  char TPCnameeps[200];
  char TPCnamepng[200];
  char ITSnameeps[200];
  char ITSnamepng[200];
  TString *None = new TString("");
  TString *Factors = None;
  if(kaonFactor!=1.0||lambdaFactor!=1.0||baryonEnhancement!=1.0){
    Factors = new TString(Form("Lambda%2.1fKaon%2.1fBaryon%2.1f",lambdaFactor,kaonFactor,baryonEnhancement));
  }
  if(TPC){
    sprintf(TPCnameeps,"pics/bkgdComponentsErrorsTPC%s.eps",Factors->Data());
    sprintf(TPCnamepng,"pics/bkgdComponentsErrorsTPC%s.png",Factors->Data());
    c->SaveAs(TPCnameeps);
    c->SaveAs(TPCnamepng);
    sprintf(TPCnameeps,"pics/bkgdComponentsTPC%s.eps",Factors->Data());
    sprintf(TPCnamepng,"pics/bkgdComponentsTPC%s.png",Factors->Data());
    c1->SaveAs(TPCnameeps);
    c1->SaveAs(TPCnamepng);
  }
  else{
    sprintf(ITSnameeps,"pics/bkgdComponentsErrorsITS%s.eps",Factors->Data());
    sprintf(ITSnamepng,"pics/bkgdComponentsErrorsITS%s.png",Factors->Data());
    c->SaveAs(ITSnameeps);
    c->SaveAs(ITSnamepng);
    sprintf(ITSnameeps,"pics/bkgdComponentsITS%s.eps",Factors->Data());
    sprintf(ITSnamepng,"pics/bkgdComponentsITS%s.png",Factors->Data());
    c1->SaveAs(ITSnameeps);
    c1->SaveAs(ITSnamepng);
  }

}
