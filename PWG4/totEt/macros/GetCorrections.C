//Create by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
//This macro takes an input file created by AliAnalysisHadEtMonteCarlo and creates an AliAnalysisHadEtCorrections for the determination of the corrected Et
// #include "TFile.h"
// #include <iostream>
// #include "AliAnalysisHadEtCorrections.h"

// #include <iostream>
// #include <TROOT.h> 
// #include <TSystem.h>
// #include "TStopwatch.h"

Float_t CorrNeutral(float ptcut, char *prodname, char *shortprodname, bool TPC, char *infilename);
TH1D *GetHistoCorrNeutral(float cut, char *name, int mycase, bool eta, int color, int marker, char *infilename);

Float_t CorrPtCut(float ptcut, char *prodname = "Enter Production Name", char *shortprodname = "EnterProductionName", char *filename="Et.ESD.new.sim.merged.root");
TH1D *GetHistoCorrPtCut(float ptcut = 0.15, char *name, char *filename);



void GetCorrections(char *prodname = "Enter Production Name", char *shortprodname = "EnterProductionName", bool TPC = true, char *infilename="Et.ESD.new.sim.merged.root", char *outfilename = "junk.root"){
    TStopwatch timer;
    timer.Start();
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libXMLIO.so");

    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libAOD.so");

    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");

    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
   gROOT->ProcessLine(".L AliAnalysisHadEtCorrections.cxx+g");


   TFile *outfile = new TFile(outfilename,"RECREATE");
   AliAnalysisHadEtCorrections *hadCorrectionEMCAL = new AliAnalysisHadEtCorrections();
   hadCorrectionEMCAL->SetName("hadCorrectionEMCAL");
   hadCorrectionEMCAL->SetEtaCut(0.7);
   float etacut = hadCorrectionEMCAL->GetEtaCut();
   cout<<"eta cut is "<<etacut<<endl;
   cout<<"My name is "<<hadCorrectionEMCAL->GetName()<<endl;
   hadCorrectionEMCAL->SetAcceptanceCorrectionFull(1.0);
   cout<<"Warning:  Acceptance corrections will have to be updated to include real acceptance maps of the EMCAL and the PHOS"<<endl;
   hadCorrectionEMCAL->SetAcceptanceCorrectionPHOS(360.0/60.0);
   hadCorrectionEMCAL->SetAcceptanceCorrectionEMCAL(360.0/60.0);
   float ptcut = 0.1;
   float neutralCorr = CorrNeutral(ptcut,prodname,shortprodname,TPC,infilename);
   hadCorrectionEMCAL->SetNeutralCorrection(neutralCorr);
   cout<<"Warning:  Setting neutral correction error bars to STAR value of +/-2%.  Use for development purposes only!"<<endl;
   hadCorrectionEMCAL->SetNeutralCorrectionLowBound(neutralCorr*0.98);
   hadCorrectionEMCAL->SetNeutralCorrectionHighBound(neutralCorr*1.02);

   //=====NEED TO SET NotHadronicCorrection!!============
   cout<<"Warning:  Have not set fNotHadronicCorrection!!"<<endl;
   

   float ptcutITS = CorrPtCut(0.1,prodname,shortprodname,infilename);
   hadCorrectionEMCAL->SetpTCutCorrectionITS(ptcutITS);
   float ptcutTPC = CorrPtCut(0.15,prodname,shortprodname,infilename);
   hadCorrectionEMCAL->SetpTCutCorrectionTPC(ptcutTPC);
   cout<<"Warning:  Setting pt cut correction error bars to STAR value of +/-3%.  Use for development purposes only!"<<endl;
   hadCorrectionEMCAL->SetpTCutCorrectionITSLowBound(ptcutITS*0.97);
   hadCorrectionEMCAL->SetpTCutCorrectionTPCLowBound(ptcutTPC*0.97);
   hadCorrectionEMCAL->SetpTCutCorrectionITSHighBound(ptcutITS*1.03);
   hadCorrectionEMCAL->SetpTCutCorrectionTPCHighBound(ptcutTPC*1.03);


   //Write the output
   outfile->cd();
   hadCorrectionEMCAL->Write();
   outfile->Write();
   outfile->Close();

  timer.Stop();
  timer.Print();
}

Float_t CorrNeutral(float ptcut, char *prodname, char *shortprodname, bool TPC, char *infilename){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TCanvas *c = new TCanvas("c","c",400,400);
  c->SetTopMargin(0.0);
  c->SetRightMargin(0.0);
  c->SetBorderSize(0);
  c->SetFillColor(0);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);

  TPad *ptpad = c->cd(1);
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

  char prefix[100];
  sprintf(prefix,"%s%2.1f",shortprodname,ptcut);

  char histoname[100];
  sprintf(histoname,"%stotal",histoname);
  int colortotal = 1;
  TH1D *total = GetHistoCorrNeutral(ptcut,histoname,4,false,colortotal,phosmarker,infilename);

  int colorallneutral = 2;
  TH1D *allneutral = GetHistoCorrNeutral(ptcut,"allneutral",3,false,colorallneutral,phosmarker,infilename);

  int colorchargedsecondary = TColor::kViolet-3;
  TH1D *chargedsecondary = GetHistoCorrNeutral(ptcut,"chargedsecondary",2,false,colorchargedsecondary,phosmarker,infilename);

  int colorneutralUndet = 4;
  TH1D *neutralUndet = GetHistoCorrNeutral(ptcut,"neutralUndet",1,false,colorneutralUndet,phosmarker,infilename);

  int colorv0 = TColor::kGreen+2;
  TH1D *v0 = GetHistoCorrNeutral(ptcut,"v0",0,false,colorv0,phosmarker,infilename);

  TF1 *func = new TF1("func","[0]",-.7,.7);
  func->SetParameter(0,0.2);
  total->Fit(func);

  //total->SetAxisRange(0.0,4);
  total->GetXaxis()->SetLabelSize(0.05);
  total->GetYaxis()->SetLabelSize(0.045);
  total->GetXaxis()->SetTitleSize(0.05);
  total->GetYaxis()->SetTitleSize(0.06);
  total->SetMaximum(0.3);
  total->SetMinimum(0.0);
  total->Draw();
  allneutral->Draw("same");
  chargedsecondary->Draw("same");
  neutralUndet->Draw("same");
  v0->Draw("same");

  TLatex *tex = new TLatex(0.161478,1.0835,prodname);
  tex->SetTextSize(0.0537634);
  tex->Draw();
  TLegend *leg2 = new TLegend(0.518321,0.746873,0.774812,0.955343);
  leg2->AddEntry(total,"Total");
  leg2->AddEntry(allneutral,"#Lambda,#bar{#Lambda},K^{0}_{S},K^{0}_{L},n,#bar{n}");
  leg2->AddEntry(neutralUndet,"K^{0}_{L},n,#bar{n}");
  leg2->SetFillStyle(0);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.0548607);
  leg2->Draw();
  char epsname[100];
  char pngname[100];
  sprintf(epsname,"pics/%s/fneutral.eps",shortprodname);
  sprintf(pngname,"pics/%s/fneutral.png",shortprodname);

  c->SaveAs(epsname);
  c->SaveAs(pngname);

  float corr = func->GetParameter(0);
  cout<<"Neutral correction: "<<1.0/(1.0-corr)<<endl;
  return 1.0/(1.0-corr);

}
TH1D *GetHistoCorrNeutral(float cut, char *name, int mycase, bool eta, int color, int marker, char *infilename){
  TFile *file = new TFile(infilename);
  TList *list = file->FindObject("out2");
  TH2F *numeratorParent; 
  switch(mycase){
  case 0:
    numeratorParent= (TH2F*)((TH2F*) out2->FindObject("EtSimulatedLambda"))->Clone("v0");
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiLambda"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedK0S"));
    break;
  case 1:
    numeratorParent = (TH2F*)((TH2F*) out2->FindObject("EtSimulatedK0L"))->Clone("Knnbar");
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedNeutron"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiNeutron"));
    break;
  case 2:
    numeratorParent = (TH2F*)((TH2F*) out2->FindObject("EtSimulatedOmega"))->Clone("ch2ndary");
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiOmega"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedXi"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedXi0"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi0"));
    break;
  case 3:
    numeratorParent= (TH2F*)((TH2F*) out2->FindObject("EtSimulatedLambda"))->Clone("allneutral");
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiLambda"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedK0S"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedK0L"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedNeutron"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiNeutron"));
    break;
  case 4:
    numeratorParent= (TH2F*)((TH2F*) out2->FindObject("EtSimulatedLambda"))->Clone("allneutral");
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiLambda"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedK0S"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedK0L"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedNeutron"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiNeutron"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedOmega"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiOmega"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedXi"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedXi0"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi0"));
    break;
  case 5:
    numeratorParent= (TH2F*)((TH2F*) out2->FindObject("EtSimulatedXi"))->Clone("allxi");
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedXi0"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi0"));
    break;
  case 6:
    numeratorParent = (TH2F*)((TH2F*) out2->FindObject("EtSimulatedOmega"))->Clone("allomega");
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiOmega"));
    break;
  case 7:
    numeratorParent = (TH2F*)((TH2F*) out2->FindObject("EtSimulatedSigma"))->Clone("allsigma");
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiSigma"));
    break;
  }

  TH2F *allhad;
  allhad=(TH2F*) ((TH2F*) out2->FindObject("EtSimulatedAllHadron"))->Clone("id");

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
  numerator->Rebin(2);
  numerator->Scale(0.5);
  numerator->SetYTitle("E_{T}^{had,sample}/E_{T}^{had,total}");
  numerator->GetYaxis()->SetTitleOffset(1.2);
  numerator->SetMarkerColor(color);
  numerator->SetLineColor(color);
  numerator->SetMarkerStyle(marker);
  return numerator;

}

TH1D *GetHistoCorrPtCut(float ptcut, char *name, char *filename){
  TFile *file = new TFile(filename);
  TList *list = file->FindObject("out2");
  TH2F *allhad = ((TH2F*) out2->FindObject("EtSimulatedAllHadron"))->Clone("allhad");

  int lowbin = allhad->GetXaxis()->FindBin(0.0);//make sure we don't accidentally get the wrong bin
  int highbin = allhad->GetXaxis()->FindBin(ptcut);
  int nbins = allhad->GetXaxis()->GetNbins();
  cout<<"Projecting from "<<allhad->GetXaxis()->GetBinLowEdge(lowbin)<<" to "<<allhad->GetXaxis()->GetBinLowEdge(highbin+1)<<endl;
  cout<<"Projecting from "<<allhad->GetXaxis()->GetBinLowEdge(lowbin)<<" to "<<allhad->GetXaxis()->GetBinLowEdge(nbins)<<endl;

  allhad->Sumw2();

  TH1D *numerator = allhad->ProjectionY("name",lowbin,highbin);
  TH1D *denominator = allhad->ProjectionY("denominator",lowbin,nbins);
  numerator->Divide(denominator);
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

Float_t CorrPtCut(float ptcut, char *prodname, char *shortprodname, char *filename){

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



  TH1D *High = GetHistoCorrPtCut(0.15-.001,"High",filename);
  TH1D *Low = GetHistoCorrPtCut(0.1-.001,"Low",filename);
  TH1D *Lowest = GetHistoCorrPtCut(0.05-.001,"Lowest",filename);

  TF1 *func = new TF1("func","[0]",-.7,.7);
  func->SetParameter(0,0.2);
  if(ptcut<.125){//its cuts
    Low->Fit(func);
  }
  else{//tpc cuts
    High->Fit(func);
  }

  High->SetMaximum(0.04);
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
  Low->Draw("same");
  Lowest->Draw("same");



  TLatex *tex = new TLatex(-0.723444,0.0373593,prodname);
  tex->SetTextSize(0.0537634);
  tex->Draw();
  TLegend *leg = new TLegend(0.217742,0.66129,0.477823,0.873656);
  leg->AddEntry(High,"p_{T} cut-off = 0.15 GeV/c");
  leg->AddEntry(Low,"p_{T} cut-off = 0.1 GeV/c");
  leg->AddEntry(Lowest,"p_{T} cut-off = 0.05 GeV/c");
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.0537634);
  leg->Draw();

  if(ptcut<.125){//its cuts
    c->SaveAs(Form("pics/%s/fptcutITS.eps",shortprodname));
    c->SaveAs(Form("pics/%s/fptcutITS.png",shortprodname));
  }
  else{//tpc cuts
    c->SaveAs(Form("pics/%s/fptcutTPC.eps",shortprodname));
    c->SaveAs(Form("pics/%s/fptcutTPC.png",shortprodname));
  }

  float corr = func->GetParameter(0);
  cout<<"Pt cut correction: "<<1.0/(1.0-corr)<<endl;
  return 1.0/(1.0-corr);
}
