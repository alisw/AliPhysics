
//Christine Nattrass, University of Tennessee at Knoxville
//This macro is for calculating the correction for the neutral energy considered part of HadEt not recorded by the tracking detectors, that from lambdas, antilambdas, K0S, K0L, neutrons, and antineutrons.
//Since PYTHIA does not get the spectra of lambdas, antilambdas, and K0S correct, this is not a very good way of determining the correction used for data and as such should be used with great caution.
//Uses the output of AliAnalysisTaskHadEt
//This is not actually what gets used in the correction class AliAnalysisHadEtCorrections - that is done in the macro GetCorrections.C - but this is useful for making plots and playing around with different options



TH1D *GetHisto(float cut = 0.12, char *name, int mycase, bool eta, int color, int marker, bool hadronic, bool reweight,float kaonFactor=1.0, float lambdaFactor = 1.0, float baryonEnhancement = 1.0){
  //TFile *file = new TFile("Et.ESD.new.sim.merged.root");
  TFile *file = new TFile("Et.ESD.new.sim.LHC10d4.pp.merged.root");
  TList *list = file->FindObject("out2");
  char *reweightname = "";
  if(reweight) reweightname = "Reweighted";
  TH2F *numeratorParent; 
  switch(mycase){
  case 0:
    numeratorParent= (TH2F*)((TH2F*) out2->FindObject(Form("EtSimulatedLambda",reweightname)))->Clone("v0");
    numeratorParent->Scale(lambdaFactor*baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedAntiLambda",reweightname)),lambdaFactor*baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedK0S%s",reweightname)),kaonFactor);
    break;
  case 1:
    numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtSimulatedK0L%s",reweightname)))->Clone("Knnbar");
    numeratorParent->Scale(kaonFactor);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedNeutron"),baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiNeutron"),baryonEnhancement);
    break;
  case 2:
    numeratorParent = (TH2F*)((TH2F*) out2->FindObject("EtSimulatedOmega"))->Clone("ch2ndary");
    numeratorParent->Scale(baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiOmega"),baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedXi"),baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi"),baryonEnhancement);
    //numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedSigma"));
    //numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiSigma"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedXi0"),baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi0"),baryonEnhancement);
    break;
  case 3:
    numeratorParent= (TH2F*)((TH2F*) out2->FindObject(Form("EtSimulatedLambda%s",reweightname)))->Clone("allneutral");
    numeratorParent->Scale(lambdaFactor*baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedAntiLambda%s",reweightname)),lambdaFactor*baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedK0S%s",reweightname)),kaonFactor);
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedK0L%s",reweightname)),kaonFactor);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedNeutron"),baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiNeutron"),baryonEnhancement);
    break;
  case 4:
    numeratorParent= (TH2F*)((TH2F*) out2->FindObject(Form("EtSimulatedLambda%s",reweightname)))->Clone("allneutral");
    numeratorParent->Scale(lambdaFactor*baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedAntiLambda%s",reweightname)),lambdaFactor*baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedK0S%s",reweightname)),kaonFactor);
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedK0L%s",reweightname)),kaonFactor);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedNeutron"),baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiNeutron"),baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedOmega"),baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiOmega"),baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedXi"),baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi"),baryonEnhancement);
    //numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedSigma"));
    //numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiSigma"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedXi0"),baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi0"),baryonEnhancement);
    break;
  case 5:
    numeratorParent= (TH2F*)((TH2F*) out2->FindObject("EtSimulatedXi"))->Clone("allxi");
    numeratorParent->Scale(baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi"),baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedXi0"),baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi0"),baryonEnhancement);
    break;
  case 6:
    numeratorParent = (TH2F*)((TH2F*) out2->FindObject("EtSimulatedOmega"))->Clone("allomega");
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiOmega"),baryonEnhancement);
    break;
  case 7:
    numeratorParent = (TH2F*)((TH2F*) out2->FindObject("EtSimulatedSigma"))->Clone("allsigma");
    numeratorParent->Scale(baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiSigma"),baryonEnhancement);
    break;
  case 8:
    numeratorParent= (TH2F*)((TH2F*) out2->FindObject(Form("EtSimulatedLambda%s",reweightname)))->Clone("allneutral");
    numeratorParent->Scale(baryonEnhancement);
    numeratorParent->Scale(lambdaFactor*baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedAntiLambda%s",reweightname)),lambdaFactor*baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedK0S%s",reweightname)),kaonFactor);
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedK0L%s",reweightname)),kaonFactor);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedNeutron"),baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiNeutron"),baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedOmega"),baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiOmega"),baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedXi"),baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi"),baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedXi0"),baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi0"),baryonEnhancement);
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedGamma"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedEta"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedPi0"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedOmega0"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedEPlus"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedEMinus"));
    break;
  case 9:
    numeratorParent= (TH2F*)((TH2F*) out2->FindObject("EtSimulatedGamma"))->Clone("allem");
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedEta"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedPi0"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedOmega0"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedEPlus"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedEMinus"));
    break;
  case 10:
    numeratorParent= (TH2F*)((TH2F*) out2->FindObject("EtSimulatedGamma"))->Clone("gamma");
    break;
  case 11:
    numeratorParent= (TH2F*)((TH2F*) out2->FindObject("EtSimulatedPi0"))->Clone("pi0");
    break;
  case 12:
    numeratorParent= (TH2F*)((TH2F*) out2->FindObject("EtSimulatedEta"))->Clone("eta");
    break;
  case 13:
    numeratorParent= (TH2F*)((TH2F*) out2->FindObject("EtSimulatedOmega0"))->Clone("Omega0");
    break;
  case 14:
    numeratorParent= (TH2F*)((TH2F*) out2->FindObject("EtSimulatedEPlus"))->Clone("electron");
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedEMinus"));
  case 15:
    numeratorParent=(TH2F*) ((TH2F*) out2->FindObject("EtSimulatedPiPlus"))->Clone("chpi");
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedPiMinus"));
    break;
  }

  TH2F *allhad;
  //allhad=(TH2F*) ((TH2F*) out2->FindObject("EtSimulatedAllHadron"))->Clone("id");
  allhad=(TH2F*) ((TH2F*) out2->FindObject("EtSimulatedPiPlus"))->Clone("id");
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedPiMinus"));
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedKMinus"),kaonFactor);
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedKPlus"),kaonFactor);
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedProton"),baryonEnhancement);
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedAntiProton"),baryonEnhancement);
  allhad->Add((TH2F*) out2->FindObject(Form("EtSimulatedLambda%s",reweightname)),lambdaFactor*baryonEnhancement);
  allhad->Add((TH2F*) out2->FindObject(Form("EtSimulatedAntiLambda%s",reweightname)),lambdaFactor*baryonEnhancement);
  allhad->Add((TH2F*) out2->FindObject(Form("EtSimulatedK0S%s",reweightname)),kaonFactor);
  allhad->Add((TH2F*) out2->FindObject(Form("EtSimulatedK0L%s",reweightname)),kaonFactor);
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedNeutron"),baryonEnhancement);
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedAntiNeutron"),baryonEnhancement);
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedEPlus"));
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedEMinus"));
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedOmega"),baryonEnhancement);
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedAntiOmega"),baryonEnhancement);
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedXi"),baryonEnhancement);
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi"),baryonEnhancement);
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedSigma"),baryonEnhancement);
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedAntiSigma"),baryonEnhancement);
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedXi0"),baryonEnhancement);
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi0"),baryonEnhancement);

  if(hadronic){//if we are getting the correction for the hadronic only case...    
    allhad->Add((TH2F*) out2->FindObject("EtSimulatedGamma"));
    allhad->Add((TH2F*) out2->FindObject("EtSimulatedEta"));
    allhad->Add((TH2F*) out2->FindObject("EtSimulatedPi0"));
    allhad->Add((TH2F*) out2->FindObject("EtSimulatedOmega0"));
    allhad->Add((TH2F*) out2->FindObject("EtSimulatedEPlus"));
    allhad->Add((TH2F*) out2->FindObject("EtSimulatedEMinus"));
  }

  numeratorParent->Sumw2();
  allhad->Sumw2();
  TH1D *denominator;
  TH1D *numerator;
  if(eta){
    int lowbin = numeratorParent->GetYaxis()->FindBin(-cut+.001);//make sure we don't accv0entally get the wrong bin
    int highbin = numeratorParent->GetYaxis()->FindBin(cut-.001);
    cout<<"Projecting from "<<numeratorParent->GetYaxis()->GetBinLowEdge(lowbin)<<" to "<<numeratorParent->GetYaxis()->GetBinLowEdge(highbin+1)<<endl;
    denominator = allhad->ProjectionX("name",lowbin,highbin);
    numerator = numeratorParent->ProjectionX("numerator",lowbin,highbin);
  }
  else{
    int lowbin = allhad->GetXaxis()->FindBin(cut);//make sure we don't accidentally get the wrong bin
    int highbin = allhad->GetXaxis()->GetNbins();
    cout<<"Projecting from "<<allhad->GetXaxis()->GetBinLowEdge(lowbin)<<" to "<<allhad->GetXaxis()->GetBinLowEdge(highbin+1)<<endl;
    numerator = numeratorParent->ProjectionY("name",lowbin,highbin);
    denominator = allhad->ProjectionY("denominator",lowbin,highbin);
  }
  numerator->Divide(denominator);
  if(hadronic){
    numerator->SetYTitle("E_{T}^{sample}/E_{T}^{total}");
  }
  else{
    numerator->SetYTitle("E_{T}^{had,sample}/E_{T}^{had,total}");
  }
  numerator->GetYaxis()->SetTitleOffset(1.2);
  numerator->SetMarkerColor(color);
  numerator->SetLineColor(color);
  numerator->SetMarkerStyle(marker);
  //numerator->Draw("e");
  return numerator;

}

void CorrNeutral(char *prodname = "LHC10d4 PYTHIA D6T 7 TeV p+p", char *shortprodname = "LHC10d4", bool hadronic = false, bool reweighted = false, float kaonFactor=1.0, float lambdaFactor = 1.0, float baryonEnhancement = 1.0){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
//   TCanvas *c = new TCanvas("c","c",800,400);
//   c->SetTopMargin(0.0);
//   c->SetRightMargin(0.0);
//   c->SetBorderSize(0);
//   c->SetFillColor(0);
//   c->SetFillColor(0);
//   c->SetBorderMode(0);
//   c->SetFrameFillColor(0);
//   c->SetFrameBorderMode(0);
//   c->Divide(2);
//  TPad *ptpad = c->cd(1);
  TCanvas *ptpad = new TCanvas("ptpad","ptpad",400,400);
  ptpad->SetTopMargin(0.04);
  ptpad->SetRightMargin(0.04);
  ptpad->SetLeftMargin(0.149288);
  ptpad->SetBorderSize(0);
  ptpad->SetFillColor(0);
  ptpad->SetFillColor(0);
  ptpad->SetBorderMode(0);
  ptpad->SetFrameFillColor(0);
  ptpad->SetFrameBorderMode(0);

  int phosmarker = 20;
  int emcalmarker = 24;
  float ptcut1 = 0.05;
  float ptcut2 = 0.1;

  int colortotal = 1;
  int casetotal = 4;
  if(hadronic) casetotal = 8;
  TH1D *PHOStotal = GetHisto(0.12,"PHOStotal",casetotal,true,colortotal,phosmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *EMCALtotal = GetHisto(0.7,"EMCALtotal",casetotal,true,colortotal,emcalmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *pt1total = GetHisto(ptcut2,"pt1total",casetotal,false,colortotal,phosmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *pt2total = GetHisto(ptcut1,"pt2total",casetotal,false,colortotal,emcalmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);

  int colorallneutral = 2;
  TH1D *PHOSallneutral = GetHisto(0.12,"PHOSallneutral",3,true,colorallneutral,phosmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *EMCALallneutral = GetHisto(0.7,"EMCALallneutral",3,true,colorallneutral,emcalmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *pt1allneutral = GetHisto(ptcut2,"pt1allneutral",3,false,colorallneutral,phosmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *pt2allneutral = GetHisto(ptcut1,"pt2allneutral",3,false,colorallneutral,emcalmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);

  int colorchargedsecondary = TColor::kViolet-3;
  TH1D *PHOSchargedsecondary = GetHisto(0.12,"PHOSchargedsecondary",2,true,colorchargedsecondary,phosmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *EMCALchargedsecondary = GetHisto(0.7,"EMCALchargedsecondary",2,true,colorchargedsecondary,emcalmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *pt1chargedsecondary = GetHisto(ptcut2,"pt1chargedsecondary",2,false,colorchargedsecondary,phosmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *pt2chargedsecondary = GetHisto(ptcut1,"pt2chargedsecondary",2,false,colorchargedsecondary,emcalmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);

  int colorneutralUndet = 4;
  TH1D *PHOSneutralUndet = GetHisto(0.12,"PHOSneutralUndet",1,true,colorneutralUndet,phosmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *EMCALneutralUndet = GetHisto(0.7,"EMCALneutralUndet",1,true,colorneutralUndet,emcalmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *pt1neutralUndet = GetHisto(ptcut2,"pt1neutralUndet",1,false,colorneutralUndet,phosmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *pt2neutralUndet = GetHisto(ptcut1,"pt2neutralUndet",1,false,colorneutralUndet,emcalmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);

  int colorv0 = TColor::kGreen+2;
  TH1D *PHOSv0 = GetHisto(0.12,"PHOSv0",0,true,colorv0,phosmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *EMCALv0 = GetHisto(0.7,"EMCALv0",0,true,colorv0,emcalmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *pt1v0 = GetHisto(ptcut2,"pt1v0",0,false,colorv0,phosmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *pt2v0 = GetHisto(ptcut1,"pt2v0",0,false,colorv0,emcalmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);

  int colorem = TColor::kCyan;
  TH1D *PHOSem = GetHisto(0.12,"PHOSem",9,true,colorem,phosmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *EMCALem = GetHisto(0.7,"EMCALem",9,true,colorem,emcalmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *pt1em = GetHisto(ptcut2,"pt1em",9,false,colorem,phosmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  TH1D *pt2em = GetHisto(ptcut1,"pt2em",9,false,colorem,emcalmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);

  PHOStotal->SetMaximum(0.5);
  PHOStotal->SetMinimum(0.0);
  if(hadronic){
    PHOStotal->SetMaximum(0.95);
  }
  PHOStotal->SetAxisRange(0.0,4);
  PHOStotal->GetXaxis()->SetLabelSize(0.05);
  PHOStotal->GetYaxis()->SetLabelSize(0.045);
  PHOStotal->GetXaxis()->SetTitleSize(0.05);
  PHOStotal->GetYaxis()->SetTitleSize(0.06);

  PHOStotal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  PHOStotal->Draw();
  EMCALtotal->Draw("same");
  PHOSallneutral->Draw("same");
  EMCALallneutral->Draw("same");
  PHOSchargedsecondary->Draw("same");
  EMCALchargedsecondary->Draw("same");
  PHOSneutralUndet->Draw("same");
  EMCALneutralUndet->Draw("same");
  PHOSv0->Draw("same");
  EMCALv0->Draw("same");
  if(hadronic){
    PHOSem->Draw("same");
    EMCALem->Draw("same");
  }
  TLatex *tex = new TLatex(0.161478,1.0835,prodname);
  tex->SetTextSize(0.0537634);
  tex->Draw();
  TLegend *leg = new TLegend(0.193781,0.796248,0.450272,0.944371);
  leg->AddEntry(PHOStotal,"|#eta|<0.12");
  leg->AddEntry(EMCALtotal,"|#eta|<0.70");
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();
  TLegend *leg2 = new TLegend(0.518321,0.612903,0.774812,0.955343);
  leg2->AddEntry(PHOStotal,"Total");
  leg2->AddEntry(PHOSallneutral,"#Lambda,#bar{#Lambda},K^{0}_{S},K^{0}_{L},n,#bar{n}");
  leg2->AddEntry(PHOSneutralUndet,"K^{0}_{L},n,#bar{n}");
  leg2->AddEntry(PHOSv0,"#Lambda,#bar{#Lambda},K^{0}_{S}");
  leg2->AddEntry(PHOSchargedsecondary,"#Xi,#Omega");
  if(hadronic) leg2->AddEntry(PHOSem,"e^{#pm},#gamma,#eta,#pi^{0},#omega");
  leg2->SetFillStyle(0);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.0548607);
  leg2->Draw();


  TCanvas *etapad = new TCanvas("etapad","etapad",400,400);//  TPad *etapad = c->cd(2);
  etapad->SetTopMargin(0.04);
  etapad->SetRightMargin(0.04);
  etapad->SetBorderSize(0);
  etapad->SetFillColor(0);
  etapad->SetFillColor(0);
  etapad->SetBorderMode(0);
  etapad->SetFrameFillColor(0);
  etapad->SetFrameBorderMode(0);
  etapad->SetLeftMargin(0.149288);



  pt1total->GetXaxis()->SetLabelSize(0.05);
  pt1total->GetYaxis()->SetLabelSize(0.045);
  pt1total->GetXaxis()->SetTitleSize(0.05);
  pt1total->GetYaxis()->SetTitleSize(0.06);
  pt1total->SetMinimum(0.0);
  pt1total->SetMaximum(0.5);
  if(hadronic){
    pt1total->SetMaximum(0.7);
  }

  TF1 *func = new TF1("func","[0]",-.7,.7);
  func->SetParameter(0,0.2);
  pt1total->Fit(func);
  int nbins = pt1total->GetNbinsX();
  float weight = 0.0;
  float totalwidth = 0.0;
  float weightPhos = 0.0;
  float totalwidthPhos = 0.0;
  for(int i=1;i<=nbins;i++){
    float width = pt1total->GetBinWidth(i);
    weight += width * pt1total->GetBinContent(i);
    totalwidth +=width;
    if(TMath::Abs(pt1total->GetBinCenter(i))<0.12){//if within the phos acceptance
    weightPhos += width * pt1total->GetBinContent(i);
    totalwidthPhos +=width;
    }
  }
  weight = weight/totalwidth;
  weightPhos = weightPhos/totalwidthPhos;
  cout<<"weight = "<<weight<<" weight phos "<<weightPhos<<endl;
  pt1total->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  pt1total->Draw();
  pt2total->Draw("same");
  pt1allneutral->Draw("same");
  pt2allneutral->Draw("same");
  pt1chargedsecondary->Draw("same");
  pt2chargedsecondary->Draw("same");
  pt1neutralUndet->Draw("same");
  pt2neutralUndet->Draw("same");
  pt1v0->Draw("same");
  pt2v0->Draw("same");
  if(hadronic){
    pt1em->Draw("same");
    pt2em->Draw("same");
  }
  TLatex *tex = new TLatex(-.65,.23,Form("%2.5f#pm%2.5f",func->GetParameter(0),func->GetParError(0)));
  tex->Draw();

  TLegend *leg3 = new TLegend(0.539259,0.801734,0.79575,0.949857);
  leg3->AddEntry(pt1total,"p_{T} cut = 0.1");
  leg3->AddEntry(pt2total,"p_{T} cut = 0.05");
  leg3->SetFillStyle(0);
  leg3->SetFillColor(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.0548607);
  leg3->Draw();
  TLegend *leg4 = new TLegend(0.199016,0.785275,0.455507,0.955343);
  leg4->AddEntry(PHOSv0,"#Lambda,#bar{#Lambda},K^{0}_{S}");
  //leg4->AddEntry(PHOSchargedsecondary,"#Sigma,#bar{#Sigma},#Xi,#bar{#Xi},#Xi^{0},#bar{#Xi^0},#Omega,#bar{#Omega}");
  leg4->AddEntry(PHOSchargedsecondary,"#Xi,#Omega");
  if(hadronic) leg4->AddEntry(PHOSem,"e^{#pm},#gamma,#eta,#pi^{0},#omega");
  leg4->SetFillStyle(0);
  leg4->SetFillColor(0);
  leg4->SetBorderSize(0);
  leg4->SetTextSize(0.0548607);
  //leg4->Draw();

  float y = 0.0237534;
  if(hadronic) y = 0.158129;
  TLatex *tex = new TLatex(-0.719565,y,prodname);
  tex->SetTextSize(0.0537634);
  tex->Draw();

  char ptnameeps[200];
  char ptnamepng[200];
  char ptnamepdf[200];
  char ptnameC[200];
  char etanameeps[200];
  char etanamepng[200];
  char etanamepdf[200];
  char etanameC[200];
  TString *Total = new TString("total");
  TString *Neutral = new TString("neutral");
  TString *Cut = Neutral;
  if(hadronic) Cut = Total;
  TString *None = new TString("");
  TString *Factors = None;
  if(kaonFactor!=1.0||lambdaFactor!=1.0||baryonEnhancement!=1.0){
    Factors = new TString(Form("Lambda%2.2fKaon%2.2fBaryon%2.2f",lambdaFactor,kaonFactor,baryonEnhancement));
  }
  sprintf(ptnameeps,"pics/f%spt%s.eps",Cut->Data(),Factors->Data());
  sprintf(ptnamepng,"pics/f%spt%s.png",Cut->Data(),Factors->Data());
  sprintf(ptnamepdf,"pics/f%spt%s.pdf",Cut->Data(),Factors->Data());
  sprintf(ptnameC,"pics/f%spt%s.C",Cut->Data(),Factors->Data());
  sprintf(etanameeps,"pics/f%seta%s.eps",Cut->Data(),Factors->Data());
  sprintf(etanamepng,"pics/f%seta%s.png",Cut->Data(),Factors->Data());
  sprintf(etanamepdf,"pics/f%seta%s.pdf",Cut->Data(),Factors->Data());
  sprintf(etanameC,"pics/f%seta%s.C",Cut->Data(),Factors->Data());
  ptpad->SaveAs(ptnameeps);
  ptpad->SaveAs(ptnamepng);
  ptpad->SaveAs(ptnamepdf);
  ptpad->SaveAs(ptnameC);
  etapad->SaveAs(etanameeps);
  etapad->SaveAs(etanamepng);
  etapad->SaveAs(etanamepdf);
  etapad->SaveAs(etanameC);

  TCanvas *c2 = new TCanvas("c2","c2",500,400);
  c2->SetTopMargin(0.03);
  c2->SetRightMargin(0.03);
  c2->SetBorderSize(0);
  c2->SetFillColor(0);
  c2->SetFillColor(0);
  c2->SetBorderMode(0);
  c2->SetFrameFillColor(0);
  c2->SetFrameBorderMode(0);


  int colorchargedall = 1;
  TH1D *chargedall = GetHisto(ptcut2,"chargedall",2,false,colorchargedall,phosmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  int colorxi =TColor::kGreen+2;
  TH1D *xi = GetHisto(ptcut2,"xi",5,false,colorxi,phosmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  int coloromega = 4;
  TH1D *omega = GetHisto(ptcut2,"omega",6,false,coloromega,phosmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  int colorsigma = 2;
  TH1D *sigma = GetHisto(ptcut2,"sigma",7,false,colorsigma,phosmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  chargedall->SetMinimum(0.0);
  chargedall->SetMaximum(0.02);
  chargedall->Draw();
  xi->Draw("same");
  omega->Draw("same");
  //sigma->Draw("same");
  TLegend *leg5 = new TLegend(0.199016,0.785275,0.455507,0.955343);
  leg5->AddEntry(chargedall,"#Sigma,#Xi,#Xi^{0},#Omega");
  //leg5->AddEntry(sigma,"#Sigma");
  leg5->AddEntry(xi,"#Xi,#Xi^{0}");
  leg5->AddEntry(omega,"#Omega");
  leg5->SetFillStyle(0);
  leg5->SetFillColor(0);
  leg5->SetBorderSize(0);
  leg5->SetTextSize(0.0548607);
  leg5->Draw();
  if(!hadronic){return;}


  TCanvas *empad = new TCanvas("empad","empad",400,400);
  empad->SetTopMargin(0.04);
  empad->SetRightMargin(0.04);
  empad->SetLeftMargin(0.149288);
  empad->SetBorderSize(0);
  empad->SetFillColor(0);
  empad->SetFillColor(0);
  empad->SetBorderMode(0);
  empad->SetFrameFillColor(0);
  empad->SetFrameBorderMode(0);
  //TH1D *EMCALem = GetHisto(0.7,"EMCALem",9,true,colorem,emcalmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);
  pt1em->SetMarkerColor(1);
  pt1em->SetLineColor(1);
  TH1D *EMCALGamma = GetHisto(ptcut2,"EMCALGamma",10,false,4,emcalmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);//0.2%
  TH1D *EMCALPi0 = GetHisto(ptcut2,"EMCALPi0",11,false,2,phosmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);//24%
  TH1D *EMCALPiCh = GetHisto(ptcut2,"EMCALPiCh",15,false,2,emcalmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);//24%
  TH1D *EMCALEta = GetHisto(ptcut2,"EMCALEta",12,false,TColor::kViolet-3,emcalmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);//1.7%
  TH1D *EMCALOmega = GetHisto(ptcut2,"EMCALOmega",13,false,TColor::kCyan,emcalmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);//0.24%
  TH1D *EMCALElectron = GetHisto(ptcut2,"EMCALElectronFrog",14,false,TColor::kGreen+2,emcalmarker,hadronic,reweighted,kaonFactor,lambdaFactor,baryonEnhancement);//0.25%
  //EMCALElectron->Draw();return;
  //EMCALPi0->Draw();return;
  pt1em->SetMinimum(0.0);
  pt1em->SetMaximum(0.3);
  pt1em->Draw();
  //EMCALGamma->Draw();
  EMCALPiCh->Scale(0.5);
  EMCALGamma->Draw("same");
  EMCALPi0->Draw("same");
  EMCALPiCh->Draw("same");
  EMCALEta->Draw("same");
  EMCALOmega->Draw("same");
  EMCALElectron->Draw("same");
  TLegend *leg20 = new TLegend(0.194444,0.215054,0.449495,0.430108);
  leg20->AddEntry(pt1em,"#gamma#eta#pi^{0}#omega e^{#pm}");
  //leg20->AddEntry(EMCALGamma,"#gamma");
  leg20->AddEntry(EMCALEta,"#eta");
  //leg20->AddEntry(EMCALOmega,"#omega");
  leg20->AddEntry(EMCALPi0,"#pi^{0}");
  leg20->AddEntry(EMCALPiCh,"#pi^{#pm}/2");
  leg20->AddEntry(EMCALElectron,"e^{#pm}");
  leg20->SetFillStyle(0);
  leg20->SetFillColor(0);
  leg20->SetBorderSize(0);
  leg20->SetTextSize(0.0548607);
  leg20->Draw();
  empad->SaveAs("pics/ftotalEmEt.eps");

  TCanvas *percentagepad = new TCanvas("percentagepad","percentagepad",400,400);
  percentagepad->SetTopMargin(0.04);
  percentagepad->SetRightMargin(0.04);
  percentagepad->SetLeftMargin(0.149288);
  percentagepad->SetBorderSize(0);
  percentagepad->SetFillColor(0);
  percentagepad->SetFillColor(0);
  percentagepad->SetBorderMode(0);
  percentagepad->SetFrameFillColor(0);
  percentagepad->SetFrameBorderMode(0);
  TH1D *percentage = EMCALPi0->Clone("Percentage");
  percentage->GetYaxis()->SetTitle("percentage of E_{T}^{em} from #pi^0");
  percentage->Divide(pt1em);
  TF1 *funcPercent = new TF1("funcPercent","[0]",-.7,.7);
  funcPercent->SetParameter(0,0.91);
  percentage->Fit(funcPercent);
  percentage->Scale(100.0);
  TF1 *funcPercent = new TF1("funcPercent","[0]",-.7,.7);
  funcPercent->SetParameter(0,0.91);
  percentage->Fit(funcPercent);

  percentage->Draw();
  percentagepad->SaveAs("pics/ftotalpercentage.eps");
}
