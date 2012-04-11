#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
//#include "RooUnfoldDagostini.h"
#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"
#include "RooUnfoldErrors.h"
#endif

TH1F *hSim;
TH2F *resolution;
//TClonesArray *histoarrayCopy = new TClonesArray("TH1D",200);
TClonesArray histoarray("TH1D",200);
  TString *ytitle;
TH1F *GetReconstructedEt(TF1 *inputfunc, int nevents, char *name, int lowbound, int highbound, int nbins);
TH1F *GetReconstructedEt(TH1F *input, int nevents, char *name, int lowbound, int highbound, int nbins);
TH1F *GetTrueEt(TF1 *inputfunc, int nevents, char *name, int lowbound, int highbound, int nbins);
Float_t GetChi2(TH1F *reconstructed, TH1F *simulated);

void Unfoldpp(int had = 2, bool its = true, int difftype = 0, char *infilename="rootFiles/LHC11b10a/Et.ESD.new.sim.LHC11b10a.root", char *datainfilename="rootFiles/LHC11a/Et.ESD.new.sim.LHC11a.root", bool zerolowetbins = false, int minbin = 4, int trainingcase =0){

#ifdef __CINT__
  gSystem->Load("~/alicework/RooUnfold-1.1.1/libRooUnfold");
  //gSystem->Load("libRooUnfold");
#endif
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  //================READING IN HISTOGRAMS===========================
  TFile *outfile = new TFile("junk.root","RECREATE");
  TFile *datafile = new TFile(datainfilename);
  TFile *file = new TFile(infilename);
  TList *list = file->FindObject("out2");
  char histoname[200];
  TString *hadStr;
  TString *longHadStr;
  TString *tail;
  TString *detector;
  ytitle = new TString("1/N_{eve}dN/dE_{T}");
  if(had==1){
    hadStr = new TString("Had");
    longHadStr = new TString("E_{T}^{had}");
  }
  if(had==0){
    hadStr = new TString("Tot");
    longHadStr = new TString("E_{T}^{tot}");
  }
  if(had==2){
    hadStr = new TString("PiKP");
    longHadStr = new TString("E_{T}^{#pi,K,p}");
  }
  switch(difftype){
  case 0:
    tail = new TString("ND");//Non-diffractive
    break;
  case 1:
    tail = new TString("SD");//Singly-diffractive
    break;
  case 2:
    tail = new TString("DD");//Doubly-diffractive
    break;
  default:
    tail = new TString("");//none
  }
  if(its) detector = new TString("ITS");
  else{ detector = new TString("TPC");}
  sprintf(histoname,"Sim%sEt%s",hadStr->Data(),tail->Data());
  //sprintf(histoname,"Sim%sEt",hadStr->Data());
  file->cd();
  hTemp = (TH1F*) out2->FindObject(histoname);
  outfile->cd();
  hSim = (TH1F*) hTemp->Clone(Form("%sCopy",histoname));

  file->cd();
  sprintf(histoname,"Reco%sEtFullAcceptanceITS%s",hadStr->Data(),tail->Data());
  hTemp = (TH1F*) out2->FindObject(histoname);
  outfile->cd();
  hITS = (TH1F*) hTemp->Clone(Form("%sSim",histoname));

  datafile->cd();
   sprintf(histoname,"Reco%sEtFullAcceptanceITS%s",hadStr->Data(),tail->Data());
   hTemp = (TH1F*) out2->FindObject(histoname);
   //   if(!hTemp){ cerr<<"no histogram "<<histoname<<endl; return;}
//   outfile->cd();
//   hITSData = (TH1F*) hTemp->Clone(Form("%sData",histoname));

  file->cd();
  sprintf(histoname,"Reco%sEtFullAcceptanceTPC%s",hadStr->Data(),tail->Data());
  //sprintf(histoname,"Reco%sEtFullAcceptanceTPC",hadStr->Data());
  cout<<"Numerator "<<histoname<<endl;
  hTemp = (TH1F*)  out2->FindObject(histoname);
  outfile->cd();
  hTPC = (TH1F*) hTemp->Clone(Form("%sSim",histoname));
  file->cd();

  hSim->Sumw2();
  hITS->Sumw2();
  hTPC->Sumw2();
  int rebin = 4;
//   if(had==2){
//     //hSim->Rebin(rebin);
//     hTPC->Rebin(rebin);
//     hITS->Rebin(rebin);
//   }

  outfile->cd();
  TH1F *hITSClone = (TH1F*)hITS->Clone("hITSClone");
  file->cd();

  //Ex SimTotEtMinusRecoEtFullAcceptanceITS
  //sprintf(histoname,"Sim%sEtMinusReco%sEtFullAcceptance%s",hadStr->Data(),hadStr->Data(),detector->Data());
  sprintf(histoname,"Sim%sEtVsReco%sEtFullAcceptance%s",hadStr->Data(),hadStr->Data(),detector->Data());
  hTemp2 = (TH2F*)out2->FindObject(histoname);
  outfile->cd();
  resolution = (TH2F*) hTemp2->Clone(Form("%sCopy",histoname));

  cout<<"Histo "<<histoname<<endl;
  resolution->Draw("colz");
  //if(had==2)  resolution->Rebin2D(rebin,rebin);
  int nbins =  resolution->GetXaxis()->GetNbins();
  if(zerolowetbins){
    for(int j=0;j<minbin;j++){//loop over bins in y=reconstructed et
      for(int i=0;i<nbins;i++){
	resolution->SetBinContent(i,j,0.0);
      }
    }
  }
  for(int i=0;i<nbins;i++){//loop over bins in x
    histoarray[i]=(TH1D*)resolution->ProjectionY(Form("tmp%i",i+1),i+1,i+1);
  }
  float lowrange = 0.0;
  float highrange = 100.0;
  file->cd();
   file->Close();
  outfile->cd();

  //================FILLING TRAINING HISTOGRAMS===========================
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~EXPONENTIAL~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //TF1 *func = new TF1("func","([0]+x*[2])/[1]*exp(-x/[1])",0,100);
  TF1 *func = new TF1("func","([0])/[1]*exp(-x/[1])",0,100);
  //TF1 *func = new TF1("func","([0])/[1]*exp(-x/[1])",0,100);
  func->FixParameter(0,1.0);
  func->SetParameter(0,1);
  func->SetParameter(1,1.0/2.23876e-01);
  func->SetParameter(2,1);
//   TF1 *funcLong = new TF1("funcLong","([0]+[1]*x+[2]*x*x+[3]*x*x*x+x^[4])/[5]*exp(-(x**[6])/[5])",lowrange,highrange);
//   funcLong->SetParameter(0,1.00467e-01);
//   funcLong->SetParameter(1,-2.82339e-01);
//   funcLong->SetParameter(2,-7.10366e-02);
//   funcLong->SetParameter(3,1.22634e-02);
//   funcLong->SetParameter(4,9.25757e-01);
//   funcLong->SetParameter(5,6.77688e-01);
//   funcLong->SetParameter(6,6.30298e-01);
  int nevents = 1e6;
  float lowbound = hSim->GetXaxis()->GetBinLowEdge(1);
  float highbound = hSim->GetXaxis()->GetBinLowEdge(nbins+1);
  TH1F *hSimExponential = GetTrueEt(func,nevents,Form("testtrue%i",i),lowbound,highbound,hITS->GetNbinsX());
  TH1F *hMeasuredExponential = GetReconstructedEt(func,nevents,Form("testsmeared%i",i),lowbound,highbound,hITS->GetNbinsX());

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~FLAT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TF1 *funcFlat = new TF1("funcFlat","[0]",0,100);
  funcFlat->FixParameter(0,0.01);
  TH1F *hSimFlat = GetTrueEt(funcFlat,nevents,Form("testtrueflat%i",i),lowbound,highbound,hITS->GetNbinsX());
  TH1F *hMeasuredFlat = GetReconstructedEt(funcFlat,nevents,Form("testsmearedflat%i",i),lowbound,highbound,hITS->GetNbinsX());

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~STRAIGHT LINE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TF1 *funcStraight = new TF1("funcStraight","[0]-[1]*x",0,100);
  funcStraight->FixParameter(0,1.0);
  funcStraight->FixParameter(1,0.01);
  TH1F *hSimStraight = GetTrueEt(funcStraight,nevents,Form("testtruestraight%i",i),lowbound,highbound,hITS->GetNbinsX());
  TH1F *hMeasuredStraight = GetReconstructedEt(funcStraight,nevents,Form("testsmearedstraight%i",i),lowbound,highbound,hITS->GetNbinsX());


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~GAUS LINE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TF1 *funcGaus = new TF1("funcGaus","[0]*exp(-x*x/[1]/[1]/TMath::Pi())",0,100);//[1] is the mean x
  funcGaus->FixParameter(0,1.0);
  funcGaus->FixParameter(1,15);
  TH1F *hSimGaus = GetTrueEt(funcGaus,nevents,Form("testtruegaus%i",i),lowbound,highbound,hITS->GetNbinsX());
  TH1F *hMeasuredGaus = GetReconstructedEt(funcGaus,nevents,Form("testsmearedgaus%i",i),lowbound,highbound,hITS->GetNbinsX());

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~DOUBLE EXPONENT LINE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TF1 *funcDoubleExponential = new TF1("funcDoubleExponential","[0]/[1]*exp(-x/[1])+[2]/[3]*exp(-x/[3])",0,100);//[1] is the mean x
  //TF1 *funcDoubleExponential = new TF1("funcDoubleExponential","[0]/[1]*exp(-x/[1])",0,100);//[1] is the mean x
  funcDoubleExponential->FixParameter(0,1.0);
  funcDoubleExponential->SetParameter(1,4.0);
  funcDoubleExponential->FixParameter(2,0.03);
  funcDoubleExponential->FixParameter(3,0.25);
  TH1F *hSimDoubleExponential = GetTrueEt(funcDoubleExponential,nevents,Form("testtruedoubleexponent%i",i),lowbound,highbound,hITS->GetNbinsX());
  TH1F *hMeasuredDoubleExponential = GetReconstructedEt(funcDoubleExponential,nevents,Form("testsmeareddoubleexponent%i",i),lowbound,highbound,hITS->GetNbinsX());


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~TRUE SIMULATED~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //trying to make something which is not identical to hTPC/hITS
  TH1F *hMeasuredSimMCSmeared =  GetReconstructedEt(hSim,nevents,Form("testtruemcsmeared%i",i),lowbound,highbound,hITS->GetNbinsX());
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~REWEIGHTED SIMULATED~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TH1F *hSimReweighted = (TH1F*) hSim->Clone("hSimReweighted");
  TF1 *fWeight = new TF1("fWeight","1-[0]*([1]-x)+[2]*([3]-x)**2",0,100);
  fWeight->SetParameter(0,-0.0005);
  fWeight->SetParameter(1,3);
  fWeight->SetParameter(2,-0.0005);
  fWeight->SetParameter(3,5);
  hSimReweighted->Divide(fWeight);
  TH1F *hMeasReweighted =  GetReconstructedEt(hSimReweighted,nevents,Form("testreweighted%i",i),lowbound,highbound,hITS->GetNbinsX());
  //float chi2 = GetChi2(hITS,test);

  //================TRAINING===========================
  //  RooUnfoldResponse (const TH1* measured, const TH1* truth, const TH2* response, const char* name= 0, const char* title= 0);  // create from already-filled histograms
  TH1F *hTrainingTruth;
  TH1F *hTrainingReconstructed;
  switch(trainingcase){
  case 0:
    cout<<"Training with pure MC"<<endl;
    hTrainingTruth = hSim;
    hTrainingReconstructed = hMeasuredSimMCSmeared;
//     if(its) hTrainingReconstructed = hITS;
//     else{hTrainingReconstructed = hTPC;}
    break;
  case 1:
    cout<<"Training with reweighted MC"<<endl;
    hTrainingTruth = hSimReweighted;
    hTrainingReconstructed = hMeasReweighted;
    break;
  case 2:
    cout<<"Training with exponential"<<endl;
    hTrainingTruth = hSimExponential;
    hTrainingReconstructed = hMeasuredExponential;
    break;
  case 3:
    cout<<"Training with flat distribution"<<endl;
    hTrainingTruth = hSimFlat       ;
    hTrainingReconstructed = hMeasuredFlat       ;
    break;
  case 4:
    cout<<"Training with straight line distribution"<<endl;
    hTrainingTruth = hSimStraight       ;
    hTrainingReconstructed = hMeasuredStraight    ;
    break;
  case 5:
    cout<<"Training with half Gaussian distribution"<<endl;
    hTrainingTruth = hSimGaus       ;
    hTrainingReconstructed = hMeasuredGaus    ;
    break;
  case 6:
    cout<<"Training with double exponential"<<endl;
    hTrainingTruth = hSimDoubleExponential;
    hTrainingReconstructed = hMeasuredDoubleExponential;
    break;
  }
  RooUnfoldResponse response(hTrainingReconstructed,hTrainingTruth,resolution,"UnfoldResponseFromHistograms","MyRooUnfoldResponse");
  //RooUnfoldResponse response(hTPC,hSimReweighted,resolution,"UnfoldResponseFromHistograms","MyRooUnfoldResponse");
  //RooUnfoldResponse response(hSmeared,hTrue,resolution,"UnfoldResponseFromHistograms","MyRooUnfoldResponse");

  //================NORMALIZING===========================
  cout<<"Normalizing..."<<endl;

  Float_t binwidth = hSim->GetBinWidth(1);
  Float_t neve = hSim->GetEntries();
  //true MC
  hSim->Scale(1.0/neve/binwidth);
  neve = hTPC->GetEntries();
  hTPC->Scale(1.0/neve/binwidth);
  neve = hITS->GetEntries();
  hITS->Scale(1.0/neve/binwidth);
  //Reweighted MC
  neve = hMeasReweighted->GetEntries();
  hMeasReweighted->Scale(1.0/neve/binwidth);
  neve = hSimReweighted->GetEntries();
  hSimReweighted->Scale(1.0/neve/binwidth);
  //Exponential MC
  neve = hMeasuredExponential->GetEntries();
  hMeasuredExponential->Scale(1.0/neve/binwidth);
  neve = hSimExponential->GetEntries();
  hSimExponential->Scale(1.0/neve/binwidth);
  //Flat MC
  neve = hMeasuredFlat->GetEntries();
  hMeasuredFlat->Scale(1.0/neve/binwidth);
  neve = hSimFlat->GetEntries();
  hSimFlat->Scale(1.0/neve/binwidth);
  //Straight MC
  neve = hMeasuredStraight->GetEntries();
  hMeasuredStraight->Scale(1.0/neve/binwidth);
  neve = hSimStraight->GetEntries();
  hSimStraight->Scale(1.0/neve/binwidth);
  //Gaus MC
  neve = hMeasuredGaus->GetEntries();
  hMeasuredGaus->Scale(1.0/neve/binwidth);
  neve = hSimGaus->GetEntries();
  hSimGaus->Scale(1.0/neve/binwidth);
  //DoubleExponential MC
  neve = hMeasuredDoubleExponential->GetEntries();
  hMeasuredDoubleExponential->Scale(1.0/neve/binwidth);
  neve = hSimDoubleExponential->GetEntries();
  hSimDoubleExponential->Scale(1.0/neve/binwidth);

  //TF1 *funcScale = new TF1("func","([0])/[1]*exp(-x/[1])",0,100);
  cerr<<"180"<<endl;
  cout<<"integrals plain "<<hSim->Integral("width")<<" "<<hITS->Integral("width")<<" "<<hTPC->Integral("width")<<endl;
  cout<<"integrals plain "<<hSimReweighted->Integral("width")<<" "<<hMeasReweighted->Integral("width")<<endl;
  cout<<"integrals plain "<<hSimExponential->Integral("width")<<" "<<hMeasuredExponential->Integral("width")<<endl;
  cout<<"integrals plain "<<hSimFlat->Integral("width")<<" "<<hMeasuredFlat->Integral("width")<<endl;
  cout<<"integrals plain "<<hSimStraight->Integral("width")<<" "<<hMeasuredStraight->Integral("width")<<endl;
  cout<<"integrals plain "<<hSimGaus->Integral("width")<<" "<<hMeasuredGaus->Integral("width")<<endl;
  cout<<"integrals plain "<<hSimDoubleExponential->Integral("width")<<" "<<hMeasuredDoubleExponential->Integral("width")<<endl;
  //================UNFOLDING MC========================
  TH1F *hTruth = hSim;
  TH1F *hMeasured;
  if(its) hMeasured = hITS;
  else{hMeasured = hTPC;}
  if(zerolowetbins){
    for(int i=0;i<minbin;i++){
      hMeasured->SetBinContent(i,0.0);
    }
  }
  RooUnfoldBayes   unfold (&response, hMeasured, 3);    // OR
  unfold.SetSmoothing(true);

  //================PLOTTING===========================


  TCanvas *canvas1 = new TCanvas("canvas1","canvas1",600,600);
  canvas1->SetTopMargin(0.020979);
  canvas1->SetRightMargin(0.0184564);
  canvas1->SetLeftMargin(0.0989933);
  canvas1->SetBottomMargin(0.101399);
  canvas1->SetBorderSize(0);
  canvas1->SetFillColor(0);
  canvas1->SetFillColor(0);
  canvas1->SetBorderMode(0);
  canvas1->SetFrameFillColor(0);
  canvas1->SetFrameBorderMode(0);
  //canvas1->Divide(2);
  canvas1->SetLogy();

  TH1D* hReco1= (TH1D*) unfold.Hreco();
  hTruth->SetMarkerStyle(22);
  hReco1->SetMarkerStyle(30);
  hReco1->SetLineColor(2);
  hReco1->SetMarkerColor(2);
  hITS->SetMarkerColor(3);
  hITS->SetLineColor(3);
  hTPC->SetMarkerColor(7);
  hTPC->SetLineColor(7);
  hTruth->Draw();
  hReco1->Draw("same");
  //hITS->Draw("same");
  //hTPC->Draw("same");
  cout<<" Means "<<hSim->GetMean()<<" "<<hReco1->GetMean()<<endl;
  hSim->GetXaxis()->SetRange(minbin,hSim->GetNbinsX());
  hReco1->GetXaxis()->SetRange(minbin,hReco1->GetNbinsX());
  cout<<" Truncated Means "<<hSim->GetMean()<<" "<<hReco1->GetMean()<<endl;

   TLegend *leg = new TLegend(0.645973,0.851399,0.746644,0.951049);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hSim,"Simulated E_{T}");
  leg->AddEntry(hReco1,"Unfolded result");
  hSim->SetLineColor(1);
  hSim->GetYaxis()->SetTitleOffset(1.3);
  hSim->GetYaxis()->SetTitle("1/N_{eve} dN/dE_{T}");
  leg->SetTextSize(0.0437063);
  leg->Draw();

  //canvas1->cd(2);
  TCanvas *canvas6 = new TCanvas("canvas6","canvas6",600,600);
  canvas6->SetTopMargin(0.020979);
  canvas6->SetRightMargin(0.0184564);
  canvas6->SetLeftMargin(0.0989933);
  canvas6->SetBottomMargin(0.101399);
  canvas6->SetBorderSize(0);
  canvas6->SetFillColor(0);
  canvas6->SetFillColor(0);
  canvas6->SetBorderMode(0);
  canvas6->SetFrameFillColor(0);
  canvas6->SetFrameBorderMode(0);

  TH1D *hReco1Clone = (TH1D*) hReco1->Clone("hReco1Clone");
  hReco1Clone->Divide(hTruth);
  hReco1Clone->Draw("");
  hReco1Clone->GetYaxis()->SetTitle("N_{eve}^{reco}/N_{eve}^{true}");
  hReco1Clone->GetXaxis()->SetTitle("E_{T}");
  hReco1Clone->SetMinimum(0.0);
  hReco1Clone->SetMaximum(2.0);
  hReco1Clone->GetXaxis()->SetRange(1,hReco1Clone->FindBin(15.0));
  //   canvas1->cd(3);

//   hReco1->Draw("");
  cout<<"Bin 1: true :"<< hSim->GetBinContent(1)<<"+/-"<<hSim->GetBinError(1)<<" reco:"<<hReco1->GetBinContent(1)<<"+/-"<<hSim->GetBinError(1)<<endl;

  unfold.PrintTable(cout,hSim);


//   TCanvas *canvas2 = new TCanvas("canvas2","canvas2",1200,600);
//   canvas2->SetTopMargin(0.020979);
//   canvas2->SetRightMargin(0.0184564);
//   canvas2->SetLeftMargin(0.0989933);
//   canvas2->SetBottomMargin(0.101399);
//   canvas2->SetBorderSize(0);
//   canvas2->SetFillColor(0);
//   canvas2->SetFillColor(0);
//   canvas2->SetBorderMode(0);
//   canvas2->SetFrameFillColor(0);
//   canvas2->SetFrameBorderMode(0);
//   fWeight->Draw();
//   TCanvas *canvas3 = new TCanvas("canvas3","canvas3",1200,600);
//   canvas3->SetTopMargin(0.020979);
//   canvas3->SetRightMargin(0.0184564);
//   canvas3->SetLeftMargin(0.0989933);
//   canvas3->SetBottomMargin(0.101399);
//   canvas3->SetBorderSize(0);
//   canvas3->SetFillColor(0);
//   canvas3->SetFillColor(0);
//   canvas3->SetBorderMode(0);
//   canvas3->SetFrameFillColor(0);
//   canvas3->SetFrameBorderMode(0);
//   hReco1->Draw();
  TCanvas *canvas4 = new TCanvas("canvas4","canvas4",1200,600);
  canvas4->SetTopMargin(0.020979);
  canvas4->SetRightMargin(0.0184564);
  canvas4->SetLeftMargin(0.0989933);
  canvas4->SetBottomMargin(0.101399);
  canvas4->SetBorderSize(0);
  canvas4->SetFillColor(0);
  canvas4->SetFillColor(0);
  canvas4->SetBorderMode(0);
  canvas4->SetFrameFillColor(0);
  canvas4->SetFrameBorderMode(0);
  hTrainingTruth->Draw();
  hTrainingReconstructed->Draw("same");
  hTrainingTruth->SetLineColor(1);
  hTrainingReconstructed->SetLineColor(2);


//   TCanvas *canvas5 = new TCanvas("canvas5","canvas5",1200,600);
//   canvas5->SetTopMargin(0.020979);
//   canvas5->SetRightMargin(0.0184564);
//   canvas5->SetLeftMargin(0.0989933);
//   canvas5->SetBottomMargin(0.101399);
//   canvas5->SetBorderSize(0);
//   canvas5->SetFillColor(0);
//   canvas5->SetFillColor(0);
//   canvas5->SetBorderMode(0);
//   canvas5->SetFrameFillColor(0);
//   canvas5->SetFrameBorderMode(0);


}

TH1F *GetReconstructedEt(TF1 *inputfunc, int nevents, char *name, int lowbound, int highbound, int nbins){

  TH1F *hResult = new TH1F(name,ytitle->Data(),nbins,lowbound,highbound);
  hResult->Sumw2();
  for(int i=0;i<nevents;i++){
    float et = inputfunc->GetRandom(lowbound,highbound);
    int mybin = resolution->GetXaxis()->FindBin(et);
    if(mybin>0 && mybin<=nbins){//then this is in our range...
      float myres = ((TH1D*)histoarray.At(mybin-1))->GetRandom();
      //float etreco = (1-myres)*et;
      float etreco = myres;
      //cout<<"et true "<<et<<" reco "<<etreco<<endl;
      hResult->Fill(etreco);
    }
  }
  //float binwidth = hResult->GetBinWidth(1);
  //hResult->Scale(1.0/nevents/binwidth);
  return hResult;
}
TH1F *GetReconstructedEt(TH1F *input, int nevents, char *name, int lowbound, int highbound, int nbins){

  TH1F *hResult = new TH1F(name,ytitle->Data(),nbins,lowbound,highbound);
  hResult->Sumw2();
  for(int i=0;i<nevents;i++){
    float et = input->GetRandom();
    int mybin = resolution->GetXaxis()->FindBin(et);
    if(mybin>0 && mybin<=nbins){//then this is in our range...
      float myres = ((TH1D*)histoarray.At(mybin-1))->GetRandom();
      //float etreco = (1-myres)*et;
      float etreco = myres;
      //cout<<"et true "<<et<<" reco "<<etreco<<endl;
      hResult->Fill(etreco);
    }
  }
  //float binwidth = hResult->GetBinWidth(1);
  //hResult->Scale(1.0/nevents/binwidth);
  return hResult;
}

TH1F *GetTrueEt(TF1 *inputfunc, int nevents, char *name, int lowbound, int highbound, int nbins){

  TH1F *hResult = new TH1F(name,ytitle->Data(),nbins,lowbound,highbound);
  hResult->Sumw2();
  for(int i=0;i<nevents;i++){
    float et = inputfunc->GetRandom(lowbound,highbound);
    hResult->Fill(et);
  }
  float binwidth = hResult->GetBinWidth(1);
  //hResult->Scale(1.0/nevents/binwidth);
  return hResult;
}


Float_t GetChi2(TH1F *reconstructed, TH1F *simulated){
  Float_t chi2 = 0;
  int nbins = 200;//reconstructed->GetNbinsX();
  float ndf=0;
  for(int i=4;i<=nbins;i++){
    float y = reconstructed->GetBinContent(i);
    float yerr = reconstructed->GetBinError(i);
    float f = simulated->GetBinContent(i);
    float ferr = simulated->GetBinError(i);
    if((yerr*yerr)>0){
//       cout<<"i "<< TMath::Power(y-f,2)/(yerr*yerr+ferr*ferr)<<endl;
//       chi2+= TMath::Power(y-f,2)/(yerr*yerr+ferr*ferr);
      //cout<<"i "<< TMath::Power(y-f,2)/(yerr*yerr)<<endl;
      chi2+= TMath::Power(y-f,2)/(yerr*yerr);
      //float relerry = yerr/y;
      //float relerrf = ferr/f;
      //cout<<"i "<<i<<" chi2 increment "<< TMath::Power(y-f,2)/(yerr*yerr)<<" y "<<y<<" f "<<f<<" yerr/y "<<relerry<<" ferr/f "<<relerrf<<endl;
      ndf++;
    }
  }
  ndf--;//There is one parameter
  if(ndf>0)return chi2/ndf;
  else{return 0;}
}
