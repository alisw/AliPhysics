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

TF1 *myGaussian = new TF1("Gaussian","1/[0]/TMath::Sqrt(2*TMath::Pi())*exp(-0.5*(x/[0])**2)",-5,5);
//gaus(0) is a substitute for [0]*exp(-0.5*((x-[1])/[2])**2) and (0) means start numbering parameters at 0.
TH1F *hSim;
TH2F *resolutionFull;
TH2F *resolution;
//TClonesArray *histoarrayCopy = new TClonesArray("TH1D",200);
TClonesArray histoarray("TH1D",200);
TString *ytitle;
TH1F *GetReconstructedEt(TF1 *inputfunc, int nevents, char *name, float lowbound, float highbound, int nbins);
TH1F *GetReconstructedEt(TH1 *input, int nevents, char *name, float lowbound, float highbound, int nbins, bool smooth);
TH1F *GetTrueEt(TF1 *inputfunc, int nevents, char *name, float lowbound, float highbound, int nbins);
Float_t GetChi2(TH1F *reconstructed, TH1F *simulated);
TH1F *RestrictAxes(TH1F *old,float minEt,float maxEt);
TH2F *RestrictAxes(TH1F *old,float minEt,float maxEt);
void Smooth(TH1F *, TF1 *,bool);
//variables I had included in arguments but which are now not used
bool zerolowetbins = false;
int minbin = 4;
bool its = false;
int difftype = 0;
int niter = 3;
bool rescaleguess = false;
//For both training and testing cases the code is:
//0 = MC
//1 = reweighted MC
//2 = exponential
//3 = flat distribution
//4 = Straight line
//5 = Half Gaussian
//6 = double exponential
//7 = Levy function
//8 = Power law
//for had =
//0 = tot et
//1 = had et
//2 = pikp et
//3 = raw et
//testmean is
//the mean for an exponential
//A for a Levy function and for a power function
//testn is the input n for the Levy and Power functions
//Bool_t testFunction = kTRUE;//for checking that the input function is not crazy so I can feel out the sanity of input parameters
Bool_t testFunction = kFALSE;//for checking that the input function is not crazy so I can feel out the sanity of input parameters
void Unfoldpp(int simdataset = 2009,int recodataset = 2009, char *outputfilename = "junk", int had = 3, int trainingcase= 0,int neveUsed = 1e7,float minEt = 0.0, float maxEt = 20.0, int varymean=0, bool unfoldData = true, float testmean = 3.0,float testn = -8,bool smooth = true){

int unfoldcase = trainingcase;
myGaussian->SetParameter(0,1.0);
#ifdef __CINT__
  gSystem->Load("~/alicework/RooUnfold-1.1.1/libRooUnfold");
  //gSystem->Load("libRooUnfold");
#endif
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  char *infilename = NULL;
  char *datainfilename = NULL;
  switch(simdataset){
  case 2009://900 GeV pp
    infilename="rootFiles/LHC11b1a/Et.ESD.sim.LHC11b1a.Run118506.root";
    break;
  case 20111://2.76 TeV pp
    infilename="rootFiles/LHC11b10a/Et.ESD.sim.LHC11b10a.Run146805.root";
    break;
  case 2010://7 TeV pp
    infilename="rootFiles/LHC10e20/Et.ESD.sim.LHC10e20.Run130795.root";
    break;
  case 2012://8 TeV pp
    infilename="rootFiles/LHC12c1b/Et.ESD.sim.LHC12c1b.Run178030.root";
    break;
  case 2013://5.5 TeV pPb
    difftype = -1;
    infilename="rootFiles/LHC13b3/Et.ESD.sim.LHC13b3.Run195483.root";
    break;
  }
  switch(recodataset){
  case 2009://900 GeV pp
    datainfilename="rootFiles/LHC10c/Et.ESD.sim.LHC10c.Run118506.root";
    break;
  case 20111://2.76 TeV pp
    datainfilename="rootFiles/LHC11a/Et.ESD.sim.LHC11a.Run146805.root";
    break;
  case 2010://7 TeV pp
    datainfilename="rootFiles/LHC10e/Et.ESD.sim.LHC10e.Run130795.root";
    break;
  case 2012://8 TeV pp
    datainfilename="rootFiles/LHC12b/Et.ESD.sim.LHC12b.Run178030.root";
    break;
  case 2013://5.5 TeV pPb
    datainfilename="rootFiles/LHC13b/Et.ESD.sim.LHC13b.Run195483.root";
    break;
  }

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
  if(had==3){
    hadStr = new TString("Raw");
    longHadStr = new TString("E_{T}^{#pi,K,p,raw}");
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
  TString *simtail = new TString("");
  if(had==3) simtail = detector;
  sprintf(histoname,"Sim%sEt%s%s",hadStr->Data(),tail->Data(),simtail->Data());
  //sprintf(histoname,"Sim%sEt",hadStr->Data());
  cout<<"Looking for "<<histoname<<endl;
  file->cd();
  hTemp = (TH1F*) out2->FindObject(histoname);
  outfile->cd();
  hSim = (TH1F*) hTemp->Clone(Form("%sCopy",histoname));

  file->cd();
  sprintf(histoname,"Reco%sEtFullAcceptanceITS%s",hadStr->Data(),tail->Data());
  cout<<"Looking for "<<histoname<<endl;
  hTemp = (TH1F*) out2->FindObject(histoname);
  outfile->cd();
  hITS = (TH1F*) hTemp->Clone(Form("%sSim",histoname));

  datafile->cd();
  sprintf(histoname,"Reco%sEtFullAcceptanceITS",hadStr->Data());
  cout<<"Looking for "<<histoname<<endl;
  hTemp = (TH1F*) out2->FindObject(histoname);
  if(unfoldData && !hTemp){ cerr<<"no histogram "<<histoname<<endl; return;}
  outfile->cd();
  hITSData = (TH1F*) hTemp->Clone(Form("%sData",histoname));

  file->cd();
  sprintf(histoname,"Reco%sEtFullAcceptanceTPC%s",hadStr->Data(),tail->Data());
  cout<<"Looking for "<<histoname<<endl;
  //sprintf(histoname,"Reco%sEtFullAcceptanceTPC",hadStr->Data());
  //cout<<"Numerator "<<histoname<<endl;
  hTemp = (TH1F*)  out2->FindObject(histoname);
  outfile->cd();
  cout<<"This histogram has "<<hTemp->GetEntries()<<" entries "<<" and a maximum of "<<hTemp->GetMaximum()<<endl;
  hTPC = (TH1F*) hTemp->Clone(Form("%sSim",histoname));


  datafile->cd();
  sprintf(histoname,"Reco%sEtFullAcceptanceTPC",hadStr->Data());
  cout<<"Looking for "<<histoname<<endl;
  hTemp = (TH1F*) out2->FindObject(histoname);
  if(unfoldData && !hTemp){ cerr<<"no histogram "<<histoname<<endl; return;}
  outfile->cd();
  hTPCData = (TH1F*) hTemp->Clone(Form("%sData",histoname));
  file->cd();


  hSim->Sumw2();
  hITS->Sumw2();
  hTPC->Sumw2();
  hITSData->Sumw2();
  hTPCData->Sumw2();
  int rebin = 1;
  //   if(had==2){
  //     //hSim->Rebin(rebin);
  //     hTPC->Rebin(rebin);
  //     hITS->Rebin(rebin);
  //   }

  outfile->cd();
  //TH1F *hITSClone = (TH1F*)hITS->Clone("hITSClone");
  file->cd();

  //Ex SimTotEtMinusRecoEtFullAcceptanceITS
  //sprintf(histoname,"Sim%sEtMinusReco%sEtFullAcceptance%s",hadStr->Data(),hadStr->Data(),detector->Data());
  sprintf(histoname,"Sim%sEtVsReco%sEtFullAcceptance%s",hadStr->Data(),hadStr->Data(),detector->Data());
  cout<<"Looking for "<<histoname<<endl;
  hTemp2 = (TH2F*)out2->FindObject(histoname);
  outfile->cd();
  resolutionFull = (TH2F*) hTemp2->Clone(Form("%sFull",histoname));
 
   resolutionFull->Rebin2D(rebin,rebin);
   hSim->Rebin(rebin);
   hITS->Rebin(rebin);
   if(hITSData) hITSData->Rebin(rebin);
   hTPC->Rebin(rebin);
   if(hTPCData) hTPCData->Rebin(rebin);
  //cout<<"Rebinning "<<rebin<<"x"<<endl;
  //float minEt = 0.000;
  //float maxEt = 100.00;
  //   hSim->GetXaxis()->SetRange(minbin,maxbin);
  //   hITS->GetXaxis()->SetRange(minbin,maxbin);
  //   hTPC->GetXaxis()->SetRange(minbin,maxbin);
  //   hITSData->GetXaxis()->SetRange(minbin,maxbin);
  //   hTPCData->GetXaxis()->SetRange(minbin,maxbin);
  //   resolution->GetXaxis()->SetRange(minbin,maxbin);
  //   resolution->GetYaxis()->SetRange(minbin,maxbin);

  //cout<<"I run from "<<hSim->GetBinLowEdge(1)<<" to "<<hSim->GetBinLowEdge(hSim->GetNbinsX()+1)<<endl;


  //cout<<"Histo "<<histoname<<endl;
  TCanvas *canvas0 = new TCanvas("canvas0","Resolution",600,600);
  canvas0->SetTopMargin(0.020979);
  canvas0->SetRightMargin(0.0184564);
  canvas0->SetLeftMargin(0.0989933);
  canvas0->SetBottomMargin(0.101399);
  canvas0->SetBorderSize(0);
  canvas0->SetFillColor(0);
  canvas0->SetFillColor(0);
  canvas0->SetBorderMode(0);
  canvas0->SetFrameFillColor(0);
  canvas0->SetFrameBorderMode(0);
//   myGaussian->Draw();
//   return;
  canvas0->SetLogz();
  //cerr<<"239"<<endl;
  if(smooth)resolutionFull->Smooth(5,"R");
  //cerr<<"241"<<endl;
  resolutionFull->Draw("colz");
  //return;
  //if(had==2)  resolution->Rebin2D(rebin,rebin);
  int nbins =  resolutionFull->GetXaxis()->GetNbins();
  if(zerolowetbins){
    for(int j=0;j<minbin;j++){//loop over bins in y=reconstructed et
      for(int i=0;i<nbins;i++){
	resolutionFull->SetBinContent(i,j,0.0);
      }
    }
  }
  for(int i=0;i<nbins;i++){//loop over bins in x
    histoarray[i]=(TH1D*)resolutionFull->ProjectionY(Form("tmp%i",i+1),i+1,i+1);
  }
  float lowrange = 0.0;
  float highrange = 100.0;
  file->cd();
  file->Close();
  outfile->cd();
  //  cerr<<"261"<<endl;
  //================FILLING TRAINING HISTOGRAMS===========================
  int neveUnused = 1e3;
  //int neveUsed = 1e7;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~EXPONENTIAL~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //TF1 *func = new TF1("func","([0]+x*[2])/[1]*exp(-x/[1])",0,100);
  TF1 *func = new TF1("func","([0])/[1]*exp(-x/[1])",0.0,100);
  //TF1 *func = new TF1("func","([0])/[1]*exp(-x/[1])",0,100);
  func->SetParameter(0,1.0);
  func->SetParameter(0,1);
  // cerr<<"271"<<endl;
  //cerr<<" can I get the mean? "<<endl;
  //cerr<<" does histo exist?" <<hSim<<endl;
  // cerr<< hSim->GetMean()<<endl;
  //func->SetParameter(1,1.0/2.23876e-01);
  float mymean = testmean;//hSim->GetMean();
  //cerr<<"273"<<endl;
  if(unfoldData) mymean = testmean;
  // cerr<<"275"<<endl;
  func->SetParameter(1,TMath::Abs((1.0+varymean*0.3)* mymean));
  func->SetParameter(2,1);
  func->SetLineColor(4);
  //   TF1 *funcLong = new TF1("funcLong","([0]+[1]*x+[2]*x*x+[3]*x*x*x+x^[4])/[5]*exp(-(x**[6])/[5])",lowrange,highrange);
  //   funcLong->SetParameter(0,1.00467e-01);
  //   funcLong->SetParameter(1,-2.82339e-01);
  //   funcLong->SetParameter(2,-7.10366e-02);
  //   funcLong->SetParameter(3,1.22634e-02);
  //   funcLong->SetParameter(4,9.25757e-01);
  //   funcLong->SetParameter(5,6.77688e-01);
  //   funcLong->SetParameter(6,6.30298e-01);
  int nevents = neveUnused;//1e7;
  //cerr<<"287"<<endl;
  if(trainingcase==2 || unfoldcase==2) nevents = neveUsed;
  float lowbound = hITS->GetXaxis()->GetBinLowEdge(1);
  float highbound = hITS->GetXaxis()->GetBinLowEdge(nbins+1);
  // cerr<<"289"<<endl;
  cout<<"Maxing histograms with ranges "<<lowbound<<" - "<<highbound<<endl;
  TH1F *hSimExponential = GetTrueEt(func,nevents,Form("testtrue%i",i),lowbound,highbound,hITS->GetNbinsX());
  TH1F *hMeasuredExponential = GetReconstructedEt(func,nevents,Form("testsmeared%i",i),lowbound,highbound,hITS->GetNbinsX(),smooth);

  // cerr<<"293"<<endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~FLAT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TF1 *funcFlat = new TF1("funcFlat","[0]",0,100);
  funcFlat->FixParameter(0,0.01);
  nevents = neveUnused;//1e7;
  if(trainingcase==3 || unfoldcase==3) nevents = neveUsed;
  TH1F *hSimFlat = GetTrueEt(funcFlat,nevents,Form("testtrueflat%i",i),lowbound,highbound,hITS->GetNbinsX());
  TH1F *hMeasuredFlat = GetReconstructedEt(funcFlat,nevents,Form("testsmearedflat%i",i),lowbound,highbound,hITS->GetNbinsX(),smooth);

  //cerr<<"302"<<endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~STRAIGHT LINE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TF1 *funcStraight = new TF1("funcStraight","[0]-[1]*x",0.0,maxEt*2);
  funcStraight->FixParameter(0,1.0);
  funcStraight->FixParameter(1,0.01);
  nevents = neveUnused;//1e7;
  if(trainingcase==4 || unfoldcase==4) nevents = neveUsed;
  TH1F *hSimStraight = GetTrueEt(funcStraight,nevents,Form("testtruestraight%i",i),lowbound,highbound,hITS->GetNbinsX());
  TH1F *hMeasuredStraight = GetReconstructedEt(funcStraight,nevents,Form("testsmearedstraight%i",i),lowbound,highbound,hITS->GetNbinsX(),smooth);


  //cerr<<"313"<<endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~GAUS LINE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TF1 *funcGaus = new TF1("funcGaus","[0]*exp(-x*x/[1]/[1]/TMath::Pi())",0.0,maxEt*2);//[1] is the mean x
  funcGaus->FixParameter(0,1.0);
  funcGaus->FixParameter(1,15);
  nevents = neveUnused;//1e7;
  if(trainingcase==5 || unfoldcase==5) nevents = neveUsed;
  TH1F *hSimGaus = GetTrueEt(funcGaus,nevents,Form("testtruegaus%i",i),lowbound,highbound,hITS->GetNbinsX());
  TH1F *hMeasuredGaus = GetReconstructedEt(funcGaus,nevents,Form("testsmearedgaus%i",i),lowbound,highbound,hITS->GetNbinsX(),smooth);

  //cerr<<"323"<<endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~DOUBLE EXPONENT LINE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TF1 *funcDoubleExponential = new TF1("funcDoubleExponential","[0]/[1]*exp(-x/[1])+[2]/[3]*exp(-x/[3])",0.0,maxEt*2);
  //TF1 *funcDoubleExponential = new TF1("funcDoubleExponential","[0]/[1]*exp(-x/[1])",0,100);
  funcDoubleExponential->SetParameter(0,1.0);
  funcDoubleExponential->SetParameter(1,testmean);
  funcDoubleExponential->FixParameter(2,0.1);
  funcDoubleExponential->FixParameter(3,1.0);
  funcDoubleExponential->SetLineColor(5);
  nevents = neveUnused;//1e7;
  if(trainingcase==6 || unfoldcase==6) nevents = neveUsed;
  TH1F *hSimDoubleExponential = GetTrueEt(funcDoubleExponential,nevents,Form("testtruedoubleexponent%i",i),lowbound,highbound,hITS->GetNbinsX());
  TH1F *hMeasuredDoubleExponential = GetReconstructedEt(funcDoubleExponential,nevents,Form("testsmeareddoubleexponent%i",i),lowbound,highbound,hITS->GetNbinsX(),smooth);


  // cerr<<"388"<<endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~TRUE SIMULATED~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //trying to make something which is not identical to hTPC/hITS
  //if(trainingcase==0 || unfoldcase==0) nevents = neveUsed;
  //TH1F *hMeasuredSimMCSmeared =  GetReconstructedEt(hSim,nevents,Form("testtruemcsmeared%i",i),lowbound,highbound,hITS->GetNbinsX());
  if(trainingcase==0 || unfoldcase==0) TH1F *hMeasuredSimMCSmeared =  hTPC;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~REWEIGHTED SIMULATED~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TH1F *hSimReweighted = (TH1F*) hSim->Clone("hSimReweighted");
  TF1 *fWeight = new TF1("fWeight","1-[0]*([1]-x)+[2]*([3]-x)**2",0.0,2.0*maxEt);
  fWeight->SetParameter(0,0);
  fWeight->SetParameter(1,3);
  fWeight->SetParameter(2,0);
  fWeight->SetParameter(3,3);
  hSimReweighted->Divide(fWeight);
  nevents = neveUnused;//1e7;
  if(trainingcase==1 || unfoldcase==1) nevents = neveUsed;
  TH1F *hMeasReweighted =  GetReconstructedEt(hSimReweighted,nevents,Form("testreweighted%i",i),lowbound,highbound,hITS->GetNbinsX());
  //float chi2 = GetChi2(hITS,test);

  //cerr<<"357"<<endl;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~LEVY~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TF1 *funcLevy = new TF1("funcLevy","[0]*x*(1+x/[1])**[2]",0,100);
  //TF1 *funcLevy = new TF1("funcLevy","[0]*x*([1]+x+[3]*x*x)**[2]",0,100);
  funcLevy->SetParameter(0,1.0);
  float n = testn;//-4;
  float A = 2;//- (1.0+varymean*0.3)* mymean *(n+3)/(n*n-4*n+5)*50;
  //cout<<"Mean "<<mymean<<" A "<<A<<" n "<<n<<endl;
  //funcLevy->SetParameter(1,(1.0+varymean*0.3)* mymean);
  funcLevy->SetParameter(1,-mymean*(testn+3)/2);
  //funcLevy->SetParameter(1,A);
  //funcLevy->SetParameter(1,1.0);
  //funcLevy->SetParameter(2,-5.96110e+00);
  funcLevy->SetParameter(2,n);
  //funcLevy->SetParameter(3,1e-1);
  nevents = neveUnused;//1e7;
  if(trainingcase==7 || unfoldcase==7) nevents = neveUsed;
  //   hSim->Fit(funcLevy,"","",1,100);
  //return;
  TH1F *hSimLevy = GetTrueEt(funcLevy,nevents,Form("testtruelevy%i",i),lowbound,highbound,hITS->GetNbinsX());
  TH1F *hMeasuredLevy = GetReconstructedEt(funcLevy,nevents,Form("testsmearedlevy%i",i),lowbound,highbound,hITS->GetNbinsX(),smooth);

  // cerr<<"379"<<endl;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~POWER~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TF1 *funcPower = new TF1("funcPower","[0]*(1+x/[1])**[2]",0,100);
  funcPower->SetParameter(0,1.0);
  float n = testn;//-100;
  funcPower->SetParameter(1,-mymean*(testn+2));
  //funcPower->SetParameter(1,(1.0+varymean*0.3)* mymean*(-n+2));
  funcPower->SetParameter(2,n);
  nevents = neveUnused;//1e7;
  //cout<<"effective scale for power function "<<funcPower->GetParameter(1)<<endl;
  if(trainingcase==8 || unfoldcase==8) nevents = neveUsed;
  //   hSim->Fit(funcLevy,"","",1,100);
  //   return;
  TH1F *hSimPower = GetTrueEt(funcPower,nevents,Form("testtruepower%i",i),lowbound,highbound,hITS->GetNbinsX());
  //cerr<<__FILE__<<" "<<__LINE__<<endl;
  TH1F *hMeasuredPower = GetReconstructedEt(funcPower,nevents,Form("testsmearedpower%i",i),lowbound,highbound,hITS->GetNbinsX(),smooth);
  // cerr<<__FILE__<<" "<<__LINE__<<endl;

  //================TRAINING===========================
  //  RooUnfoldResponse (const TH1* measured, const TH1* truth, const TH2* response, const char* name= 0, const char* title= 0);  // create from already-filled histograms
  TH1F *hTrainingTruth;
  TH1F *hTrainingReconstructed;
  switch(trainingcase){
  case 0:
    cout<<"Training with pure MC"<<endl;
    hTrainingTruth = hSim;
    hTrainingReconstructed = hMeasuredSimMCSmeared;
    cout<<"Training reconstructed is "<<hMeasuredSimMCSmeared->GetName()<<" entries "<<hMeasuredSimMCSmeared->GetEntries()<<" max "<<hMeasuredSimMCSmeared->GetMaximum()<<endl;
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
  case 7:
    cout<<"Training with Levy distribution"<<endl;
    hTrainingTruth = hSimLevy       ;
    hTrainingReconstructed = hMeasuredLevy       ;
    break;
  case 8:
    cout<<"Training with Power law distribution"<<endl;
    hTrainingTruth = hSimPower       ;
    hTrainingReconstructed = hMeasuredPower       ;
    break;
  }
  cout<<"nbins sim "<<hSim->GetNbinsX()<<" data "<<hITS->GetNbinsX()<<" training truth "<<hTrainingTruth->GetNbinsX()<<" training reconstructed "<<hTrainingReconstructed->GetNbinsX()<<endl;
  cout<<"max sim "<<hSim->GetBinLowEdge(401)<<" data "<<hITS->GetBinLowEdge(401)<<" training truth "<<hTrainingTruth->GetBinLowEdge(401)<<" training reconstructed "<<hTrainingReconstructed->GetBinLowEdge(401)<<endl;
  cout<<"min sim "<<hSim->GetBinLowEdge(1)<<" data "<<hITS->GetBinLowEdge(1)<<" training truth "<<hTrainingTruth->GetBinLowEdge(1)<<" training reconstructed "<<hTrainingReconstructed->GetBinLowEdge(1)<<endl;
  //return;
  hSim = RestrictAxes(hSim,minEt,maxEt);
  hITS = RestrictAxes(hITS,minEt,maxEt);
  hTPC = RestrictAxes(hTPC,minEt,maxEt);
  hITSData = RestrictAxes(hITSData,minEt,maxEt);
  hTPCData = RestrictAxes(hTPCData,minEt,maxEt);
  resolution = RestrictAxes(resolutionFull,minEt,maxEt);
  hTrainingReconstructed = RestrictAxes(hTrainingReconstructed,minEt,maxEt);
  hTrainingTruth = RestrictAxes(hTrainingTruth,minEt,maxEt);
//   if(rescaleguess && unfoldData){
//     cout<<"Rescaling the guess by MC reconstructed/data reconstructed"<<endl;
//     TH1F *hScale = hTrainingReconstructed->Clone("hScale");
//     if(its) hScale->Divide(hITSData);
//     else{hScale->Divide(hTPCData);}
//     //hTrainingReconstructed is now ~MC/data
//     //TH1F *hTrainingTruthCopy = (TH1F*)hTrainingTruth->Clone("hTrainingTruthCopy");
//     hTrainingTruth->Divide(hScale);//This is now MC/(MC/data) and in principle is the best guess at data?
//     //hTrainingTruth->Add(hTrainingTruthCopy);
//     float scale = hTrainingTruth->Integral();
//     cout<<"Integral "<<scale<<endl;
//     hTrainingTruth->Scale(1.0/scale);

//     hTrainingReconstructed = GetReconstructedEt(hTrainingTruth,nevents,Form("testMCreweightedbydata%i",i),lowbound,highbound,hITS->GetNbinsX());
//   }

  RooUnfoldResponse response(hTrainingReconstructed,hTrainingTruth,resolution,"UnfoldResponseFromHistograms","MyRooUnfoldResponse");
  //RooUnfoldResponse response(hTPC,hSimReweighted,resolution,"UnfoldResponseFromHistograms","MyRooUnfoldResponse");
  //RooUnfoldResponse response(hSmeared,hTrue,resolution,"UnfoldResponseFromHistograms","MyRooUnfoldResponse");

  //================NORMALIZING===========================
  //cout<<"Normalizing..."<<endl;

  Float_t binwidth = hSim->GetBinWidth(1);
  Float_t neve = hSim->GetEntries();
  if(neve<=1e-5){
    cerr<<"Warning!! simulation histo not filled!  Stopping!"<<endl;
    return;
  }
  if(binwidth<=1e-5){
    cerr<<"Warning!! binwidth!  Stopping!"<<endl;
    return;
  }
  //true MC
  hSim->Scale(1.0/neve/binwidth);
  hSim->Scale(1.0/hSim->Integral());
  neve = hTPC->GetEntries();
  if(neve<=1e-5){
    cerr<<"Warning!! Histo not filled!  Stopping!"<<endl;
    return;
  }
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
  //Levy MC
  neve = hMeasuredLevy->GetEntries();
  hMeasuredLevy->Scale(1.0/neve/binwidth);
  neve = hSimLevy->GetEntries();
  hSimLevy->Scale(1.0/neve/binwidth);
  if(hITSData){
    neve = hITSData->GetEntries();
    hITSData->Scale(1.0/neve/binwidth);
    hITSData->Scale(1.0/hITSData->Integral());
  }
  if(hITSData){
    neve = hTPCData->GetEntries();
    hTPCData->Scale(1.0/neve/binwidth);
  }
  if(rescaleguess){
    neve = hTrainingReconstructed->GetEntries();
    hTrainingReconstructed->Scale(1.0/neve/binwidth);
    hTPCData->Scale(1.0/hTPCData->Integral());
  }

  //TF1 *funcScale = new TF1("func","([0])/[1]*exp(-x/[1])",0,100);
//   cout<<"integrals plain "<<hSim->Integral("width")<<" "<<hITS->Integral("width")<<" "<<hTPC->Integral("width")<<endl;
//   cout<<"integrals plain "<<hSimReweighted->Integral("width")<<" "<<hMeasReweighted->Integral("width")<<endl;
//   cout<<"integrals plain "<<hSimExponential->Integral("width")<<" "<<hMeasuredExponential->Integral("width")<<endl;
//   cout<<"integrals plain "<<hSimFlat->Integral("width")<<" "<<hMeasuredFlat->Integral("width")<<endl;
//   cout<<"integrals plain "<<hSimStraight->Integral("width")<<" "<<hMeasuredStraight->Integral("width")<<endl;
//   cout<<"integrals plain "<<hSimGaus->Integral("width")<<" "<<hMeasuredGaus->Integral("width")<<endl;
//   cout<<"integrals plain "<<hSimDoubleExponential->Integral("width")<<" "<<hMeasuredDoubleExponential->Integral("width")<<endl;
  //================SETTING TEST CASE========================

  TH1F *hTruth;
  TH1F *hMeasured;
  switch(unfoldcase){
  case 0:
    cout<<"Test case is pure MC"<<endl;
    hTruth = hSim;
    if(its) hMeasured = hITS;
    else{hMeasured = hTPC;}
    break;
  case 1:
    cout<<"Test case is a reweighted MC"<<endl;
    hTruth = hSimReweighted;
    hMeasured = hMeasReweighted;
    break;
  case 2:
    cout<<"Test case is a exponential"<<endl;
    hTruth = hSimExponential;
    hMeasured = hMeasuredExponential;
    break;
  case 3:
    cout<<"Test case is a flat distribution"<<endl;
    hTruth = hSimFlat       ;
    hMeasured = hMeasuredFlat       ;
    break;
  case 4:
    cout<<"Test case is a straight line distribution"<<endl;
    hTruth = hSimStraight       ;
    hMeasured = hMeasuredStraight    ;
    break;
  case 5:
    cout<<"Test case is a half Gaussian distribution"<<endl;
    hTruth = hSimGaus       ;
    hMeasured = hMeasuredGaus    ;
    break;
  case 6:
    cout<<"Test case is a double exponential"<<endl;
    hTruth = hSimDoubleExponential;
    hMeasured = hMeasuredDoubleExponential;
    break;
  case 7:
    cout<<"Test case is a Levy distribution"<<endl;
    hTruth = hSimLevy       ;
    hMeasured = hMeasuredLevy       ;
    break;
  case 8:
    cout<<"Test case is a Power law"<<endl;
    hTruth = hSimPower       ;
    hMeasured = hMeasuredPower       ;
    break;
  }
  if(unfoldData){
    if(its) hMeasured = hITSData;
    else{hMeasured = hTPCData;}
  }
  cout<<"MC truth: "<<hTruth->GetMean()<<" Mean observed "<<hMeasured->GetMean()<<endl;
  hTruth->GetXaxis()->SetTitle(hSim->GetTitle());
  hMeasured->GetXaxis()->SetTitle(hSim->GetTitle());




  //cout<<"Histo "<<histoname<<endl;
  TCanvas *canvasn1 = new TCanvas("canvasn1","Training reconstructed vs data reconstructed",600,600);
  canvasn1->SetTopMargin(0.020979);
  canvasn1->SetRightMargin(0.0184564);
  canvasn1->SetLeftMargin(0.0989933);
  canvasn1->SetBottomMargin(0.101399);
  canvasn1->SetBorderSize(0);
  canvasn1->SetFillColor(0);
  canvasn1->SetFillColor(0);
  canvasn1->SetBorderMode(0);
  canvasn1->SetFrameFillColor(0);
  canvasn1->SetFrameBorderMode(0);
  //canvasn1->Divide(2);
  canvasn1->SetLogy();
  hMeasured->SetMarkerStyle(22);
  hTrainingReconstructed->SetMarkerStyle(26);
  hTrainingReconstructed->SetMarkerColor(2);
  hTrainingReconstructed->SetLineColor(2);
  cout<<"Measured histo is "<<hMeasured->GetName()<<" has entries "<<hMeasured->GetEntries()<<" max "<<hMeasured->GetMaximum()<<endl;
  hMeasured->Draw();
  cout<<"Reconstructed histo is "<<hTrainingReconstructed->GetName()<<" has entries "<<hTrainingReconstructed->GetEntries()<<" max "<<hTrainingReconstructed->GetMaximum()<<endl;
  hTrainingReconstructed->Draw("same");
  TLegend *legendn1 = new TLegend(0.437919,0.800699,0.966443,0.958042);
  legendn1->AddEntry(hMeasured,"Measured");
  legendn1->AddEntry(hTrainingReconstructed,"Simulated reconstructed");
  legendn1->Draw();
  //cout<<"Reconstructed mean "<<hMeasured->GetMean()<<" Simulated reconstructed mean "<<hTrainingReconstructed->GetMean()<<endl;
  //return;
  float matchval = 1.0;
  if(unfoldData){
    int binsim = hTrainingReconstructed->FindBin(matchval+.01);
    float simval = hTrainingReconstructed->GetBinContent(binsim);
    //cerr<<__FILE__<<" "<<__LINE__<<endl;
    int binreco = hMeasured->FindBin(matchval+0.01);
    float recoval = hMeasured->GetBinContent(binreco);
    cout<<"matching at "<<matchval<<" simval "<<simval<<" recoval "<<recoval<<" bins "<<binsim<<" "<<binreco<<endl;
    if(recoval>0.0) hMeasured->Scale(simval/recoval);
    else{cerr<<"Uh-oh!  Can't rescale.  :("<<endl;}
  }
  //return;
  //   for(int i=1;i<=hMeasured->GetNbinsX();i++){
  //     cout<<hMeasured->GetBinContent(i);
  //     if(i!=hMeasured->GetNbinsX())cout<<", ";
  //   }
  //   cout<<endl;
  //hMeasured->Fit(funcPower,"N","",1,100);
  //funcLevy->Draw("same");
  //hMeasured->Fit(funcLevy,"","",1,100);
  //return;
  //return;
//   if(simdataset==2009){
//     cout<<"NOT doing full unfolding because this is 2009.  Exiting."<<endl;
//     cout<<"Mean is "<<hTruth->GetMean()<<" ";
//     GetChi2(hMeasured,hTrainingReconstructed);
//     canvasn1->SaveAs(Form("%s%s.C",outputfilename,"RecoVSSimReco"));
//     cout<<"Saving "<<Form("%s%s.C",outputfilename,"RecoVSSimReco")<<endl;
//     return;
//   }
  if(testFunction) return;
  //================UNFOLDING===========================
  if(zerolowetbins){
    for(int i=0;i<minbin;i++){
      hMeasured->SetBinContent(i,0.0);
    }
  }
  RooUnfoldBayes   unfold (&response, hMeasured, niter);    // OR
  unfold.SetSmoothing(true);

  //================PLOTTING===========================


  TH1D* hReco1= (TH1D*) unfold.Hreco();
  //cout<<"MEAN WITHOUT ANY MASSAGING "<<hReco1->GetMean()<<endl;

  TCanvas *canvas7 = new TCanvas("canvas7","Reconstructed Unfolded E_{T}",600,600);
  canvas7->SetTopMargin(0.020979);
  canvas7->SetRightMargin(0.0184564);
  canvas7->SetLeftMargin(0.0989933);
  canvas7->SetBottomMargin(0.101399);
  canvas7->SetBorderSize(0);
  canvas7->SetFillColor(0);
  canvas7->SetFillColor(0);
  canvas7->SetBorderMode(0);
  canvas7->SetFrameFillColor(0);
  canvas7->SetFrameBorderMode(0);
  canvas7->SetLogy();
  hReco1->Draw();

  TCanvas *canvas1 = new TCanvas("canvas1","Results with fits",600,600);
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

  //rescale just in case
  //   float myscale = hTruth->Integral("width") / hReco1->Integral("width");
  //   cout<<"Rescaling result by "<<myscale<<endl;
  //   hReco1->Scale(myscale);
  hReco1->Fit(func,"NI","",0.0,1.0);
  funcDoubleExponential->SetParameter(1,func->GetParameter(1));
  hReco1->Fit(funcDoubleExponential,"NI","",0.0,1.0);
  float xgoal1 = 1.0;
  float xgoal2 = 0.01;
  if(recodataset ==2010){
    float xgoal1 = 1.0;
    float xgoal2 = 2.0;
  }
  float x1 = hReco1->GetBinCenter(hReco1->FindBin(xgoal1));
  float y1 = hReco1->GetBinContent(hReco1->FindBin(xgoal1));
  float x2 = hReco1->GetBinCenter(hReco1->FindBin(xgoal2));
  float y2 = hReco1->GetBinContent(hReco1->FindBin(xgoal2));
//   float x2 = hReco1->GetBinCenter(1);
//   float y2 = hReco1->GetBinContent(1);
  float myb = -TMath::Log(y1/y2)/(x1-x2);
  float myA = y1/myb/exp(-myb*x1);
  func->SetParameter(0,myA);
  func->SetParameter(1,1.0/myb);
  //cout<<"x1 "<<x1<<" y1 "<<y1<<" x2 "<<x2<<" y2 "<<y2<<" b "<<1/myb<<" A "<<myA/myb<<endl;


  hTruth->SetMarkerStyle(22);
  hReco1->SetMarkerStyle(30);
  hReco1->SetLineColor(2);
  hReco1->SetMarkerColor(2);
  hITS->SetMarkerColor(3);
  hITS->SetLineColor(3);
  hTPC->SetMarkerColor(7);
  hTPC->SetLineColor(7);
  //hTruth->GetXaxis()->SetRange(0,hTruth->FindBin(10.0));
  hReco1->Draw();
  if(!unfoldData) hTruth->Draw("same");
  //hITS->Draw("same");
  //hTPC->Draw("same");
  func->Draw("same");
  funcDoubleExponential->Draw("same");
  hTruth->GetXaxis()->SetRange(0,hTruth->GetNbinsX());
  //cout<<" Means "<<hTruth->GetMean()<<" "<<hReco1->GetMean()<<endl;
//   hTruth->GetXaxis()->SetRange(hTruth->FindBin(1.),hTruth->GetNbinsX());
//   hReco1->GetXaxis()->SetRange(hTruth->FindBin(1.),hReco1->GetNbinsX());
  hTruth->GetXaxis()->SetRange(0,hTruth->FindBin(20.));
  hReco1->GetXaxis()->SetRange(0,hTruth->FindBin(20.));
  //cout<<" Truncated Means "<<hTruth->GetMean()<<" "<<hReco1->GetMean()<<endl;
  hTruth->GetXaxis()->SetRange(0,hTruth->GetNbinsX());
  hReco1->GetXaxis()->SetRange(0,hReco1->GetNbinsX());

  TLegend *leg = new TLegend(0.505034,0.744755,0.605705,0.952797);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  //leg->AddEntry(hSim,"Simulated E_{T}");
  leg->AddEntry(hReco1,"Unfolded result");
  leg->AddEntry(func,"Exponential");
  leg->AddEntry(funcDoubleExponential,"Double exponential");
  hSim->SetLineColor(1);
  hSim->GetYaxis()->SetTitleOffset(1.3);
  hSim->GetYaxis()->SetTitle("1/N_{eve} dN/dE_{T}");
  leg->SetTextSize(0.0437063);
  leg->Draw();
  // canvas1->SaveAs(Form("%s%s.png",outputfilename,"Extrap"));
  //canvas1->SaveAs(Form("%s%s.eps",outputfilename,"Extrap"));
  //canvas1->SaveAs(Form("%s%s.C",outputfilename,"Extrap"));
  //canvas1->cd(2);
  TCanvas *canvas6 = new TCanvas("canvas6","Unfolded data/simulated MC truth",600,600);
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
  if(unfoldData){
    int binsim = hTruth->FindBin(matchval+.01);
    float simval = hTruth->GetBinContent(binsim);
    //cerr<<__FILE__<<" "<<__LINE__<<endl;
    int binreco = hReco1Clone->FindBin(matchval+0.01);
    float recoval = hReco1Clone->GetBinContent(binreco);
    cout<<"matching at "<<matchval<<" simval "<<simval<<" recoval "<<recoval<<" bins "<<binsim<<" "<<binreco<<endl;
    if(recoval>0.0) hReco1Clone->Scale(simval/recoval);
  }
  hReco1Clone->Divide(hTruth);
  hReco1Clone->Draw("");
  hReco1Clone->GetYaxis()->SetTitle("N_{eve}^{reco}/N_{eve}^{true}");
  hReco1Clone->GetXaxis()->SetTitle("E_{T}");
  hReco1Clone->SetMinimum(0.0);
  hReco1Clone->SetMaximum(2.0);
  hReco1Clone->GetXaxis()->SetRange(1,hReco1Clone->FindBin(maxEt));
  //   canvas1->cd(3);

  //   hReco1->Draw("");
  //cout<<"Bin 1: true :"<< hSim->GetBinContent(1)<<"+/-"<<hSim->GetBinError(1)<<" reco:"<<hReco1->GetBinContent(1)<<"+/-"<<hSim->GetBinError(1)<<endl;

  //  unfold.PrintTable(cout,hSim);

  if(trainingcase==1 || unfoldcase==1){
    TCanvas *canvas2 = new TCanvas("WeightingFunction","WeightingFunction",600,600);
    canvas2->SetTopMargin(0.020979);
    canvas2->SetRightMargin(0.0184564);
    canvas2->SetLeftMargin(0.0989933);
    canvas2->SetBottomMargin(0.101399);
    canvas2->SetBorderSize(0);
    canvas2->SetFillColor(0);
    canvas2->SetFillColor(0);
    canvas2->SetBorderMode(0);
    canvas2->SetFrameFillColor(0);
    canvas2->SetFrameBorderMode(0);
    fWeight->Draw();
  }
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
//   TCanvas *canvas4 = new TCanvas("canvas4","canvas4",1200,600);
//   canvas4->SetTopMargin(0.020979);
//   canvas4->SetRightMargin(0.0184564);
//   canvas4->SetLeftMargin(0.0989933);
//   canvas4->SetBottomMargin(0.101399);
//   canvas4->SetBorderSize(0);
//   canvas4->SetFillColor(0);
//   canvas4->SetFillColor(0);
//   canvas4->SetBorderMode(0);
//   canvas4->SetFrameFillColor(0);
//   canvas4->SetFrameBorderMode(0);
//   hTrainingReconstructed->Draw();
//   hTrainingTruth->Draw("same");
//   hTrainingTruth->SetLineColor(1);
//   hTrainingReconstructed->SetLineColor(2);


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

  //Calculating mean/mean error with covariant matrix
  float cutoffLow = hReco1->GetBinLowEdge(hReco1->FindBin(1.0+0.01));
  float cutoffHigh = hReco1->GetBinLowEdge(hReco1->FindBin(30+0.01));
  if(had==1) cutoffHigh = hReco1->GetBinLowEdge(hReco1->FindBin(20+0.01));
  int cutoffbin = hReco1->FindBin(1.01);//Get the bin at 1 GeV, with a little wiggle room to be sure it returns the bin I mean
  int cutoffbinHigh = hReco1->FindBin(30.01);//Get the bin at 1 GeV, with a little wiggle room to be sure it returns the bin I mean
  if(had==1) cutoffbinHigh = hReco1->FindBin(20.01);//Get the bin at 1 GeV, with a little wiggle room to be sure it returns the bin I mean  
  if(had==2) cutoffbinHigh = hReco1->FindBin(17.01);//Get the bin at 1 GeV, with a little wiggle room to be sure it returns the bin I mean
  int nbins = hReco1->GetNbinsX();
  TMatrixD cov = unfold.Ereco();
  TVectorD result = unfold.Vreco();
  TVectorD resultError = unfold.ErecoV();
  //   TMatrixD cov = unfold.GetMeasuredCov();
  //   TVectorD result = unfold.Vmeasured();
  //   TVectorD resultError = unfold.Emeasured();
  float cutoff = hReco1->GetBinLowEdge(cutoffbin);
  float binwidth = hReco1->GetBinWidth(1);
  float mean = 0;
  float meanerr = 0;//this is actually the error squared but I'm just trying a different way of calculating it
  float meanvar = 0;
  float total = 0;
  //truncated values
  float meantrunc = 0;
  float meanerrtrunc = 0;//this is actually the error squared but I'm just trying a different way of calculating it
  float meanvartrunc = 0;
  float totaltrunc = 0;
  for(int i=0;i<cutoffbinHigh;i++){
    total += result(i)*binwidth;
    float energyi = hReco1->GetBinCenter(i+1);
    mean += result(i)*binwidth*energyi;
    meanerr += resultError(i)*resultError(i)*binwidth*energyi*binwidth*energyi;
    if(i>=cutoffbin){
      totaltrunc += result(i)*binwidth;
      meantrunc += result(i)*binwidth*energyi;
      meanerrtrunc += resultError(i)*resultError(i)*binwidth*energyi*binwidth*energyi;
    }
    for(int j=0;j<nbins;j++){
      float energyj = hReco1->GetBinCenter(j+1);
      meanvar += energyi*binwidth  *  energyj*binwidth  * cov(i,j);
      if(i>=cutoffbin && j>=cutoffbin){
	meanvartrunc += energyi*binwidth  *  energyj*binwidth  * cov(i,j);
      }
      //       if(i==j && i<40 && j<40 && cov(i,j)>0.0){
      // 	cout<<"i "<<i<<" bin center "<<hReco1->GetBinCenter(i+1)<<" j "<<j<<" cov "<<cov(i,j)<<" err "<< result(i)*binwidth*energyi  *  result(j)*binwidth*energyj  * cov(i,j);
      // 	if(i==j) cout<<" error "<<resultError(i)<<" rel. err "<< resultError(i)/result(i)*100.0<<"%";
      // 	cout<< endl;
      //       }
    }
    //     cout<<"i "<<i;
    //     cout<<" result "<<result(i);
    //     cout<<" +/- "<<resultError(i);
    //     cout<<endl;
  }
  cout<<"Mean "<<mean<<" total "<<total<<" mean err ";
  if(meanvar>0) cout<<TMath::Sqrt(meanvar);
  cout<<" or ";
  if(meanerr>0)cout<<TMath::Sqrt(meanerr);
  cout<<endl;
  cout<<"Mean trunc";
  if(totaltrunc>0)cout<<meantrunc/totaltrunc;
  cout<<endl;
  TF1 *funcMean = new TF1("funcMean","x*([0])/[1]*exp(-x/[1])",0,100);
  for(int i=0; i<2;i++){funcMean->SetParameter(i,func->GetParameter(i));}
  float expoExtrap = -1;
  if(funcMean->GetParameter(1)!=0) funcMean->Integral(0.0,cutoffLow);
  float expoExtrapTotal = -1;
  if(func->GetParameter(1)!=0) func->Integral(0.0,cutoffLow);
  TF1 *funcDoubleExponentialMean = new TF1("funcDoubleExponentialMean","x*[0]/[1]*exp(-x/[1])+x*[2]/[3]*exp(-x/[3])",0,100);
  for(int i=0; i<4;i++){funcDoubleExponentialMean->SetParameter(i,funcDoubleExponential->GetParameter(i));}
  float dblExpoExtrap = funcDoubleExponentialMean->Integral(0.0,cutoffLow);
  float dblExpoExtrapTotal = funcDoubleExponential->Integral(0.0,cutoffLow);
  float total = hReco1->Integral(1,hReco1->FindBin(cutoffHigh-0.01),"width");
  float truncated = hReco1->Integral(hReco1->FindBin(cutoffLow+0.01),hReco1->FindBin(cutoffHigh-0.01),"width");
  hReco1->GetXaxis()->SetRange(1,hReco1->FindBin(cutoffHigh));
  hReco1->GetXaxis()->SetRange(hReco1->FindBin(cutoffLow),hReco1->FindBin(cutoffHigh));
  cout<<"Truncated Mean "<<hReco1->GetMean()<<" total "<<total<<" mean err ";
  if(meanvar>0) cout<<TMath::Sqrt(meanvar);
  cout<<endl;
  hReco1->GetXaxis()->SetRange(1,hReco1->GetNbinsX());
//   cout<<"Extrapolations:  exponential "<<expoExtrap<<" double exponential "<<dblExpoExtrap<<endl;
//   cout<<"Total: high ";
//   if((totaltrunc+dblExpoExtrapTotal)>0) cout<<(dblExpoExtrap+meantrunc)/(totaltrunc+dblExpoExtrapTotal);
//   cout<<" best ";
//   if((totaltrunc+expoExtrapTotal)>0) cout<<(expoExtrap+meantrunc)/(totaltrunc+expoExtrapTotal);
//   cout<<" low "<<mean<<endl;
//   cout<<"Raw mean "<<hReco1->GetMean();
//   cout<<" Total scaled: low ";
//   if((totaltrunc+dblExpoExtrapTotal)>0) cout<<(dblExpoExtrap+meantrunc)/(totaltrunc+dblExpoExtrapTotal);
//   cout<<" best ";
//   if(totaltrunc+expoExtrapTotal>0) cout<<(expoExtrap+meantrunc)/(totaltrunc+expoExtrapTotal);
//   cout<<" high ";
//   if(total>0) cout<<mean/total;
//   cout<<endl;
  canvas1->cd();
  funcDoubleExponential->SetLineStyle(2);
  func->Draw("same");
  funcDoubleExponential->Draw("same");
  // canvas1->SaveAs(Form("%s.C",outputfilename));
  //canvas1->SaveAs(Form("%s.eps",outputfilename));
  //canvas1->SaveAs(Form("%s.png",outputfilename));



  TCanvas *canvas5 = new TCanvas("RecovsSimReco","RecovsSimReco",600,600);
  canvas5->SetTopMargin(0.020979);
  canvas5->SetRightMargin(0.0184564);
  canvas5->SetLeftMargin(0.0989933);
  canvas5->SetBottomMargin(0.101399);
  canvas5->SetBorderSize(0);
  canvas5->SetFillColor(0);
  canvas5->SetFillColor(0);
  canvas5->SetBorderMode(0);
  canvas5->SetFrameFillColor(0);
  canvas5->SetFrameBorderMode(0);
  canvas5->SetLogy();
  hMeasured->Draw();
  nevents = neveUsed;

  TH1F *hMeasuredSim = GetReconstructedEt((TH1*)hReco1,(int) nevents,(char *)"SmearedTruth",(float) hReco1->GetBinLowEdge(1),(float) hReco1->GetBinLowEdge(hReco1->GetNbinsX()),(int) hReco1->GetNbinsX());//(TH1D*) unfold.Hmeasured();
  float myscale2 = hMeasured->GetBinContent(hMeasured->FindBin(matchval)) / hMeasuredSim->GetBinContent(hMeasuredSim->FindBin(matchval));
  //cout<<"p1 "<<hMeasured->GetBinContent(hMeasured->FindBin(matchval))<<" "<<hMeasured->GetBinContent(hMeasured->FindBin(matchval))<<" p2 "<<
  //  cout<<"Scaling by "<<hMeasured->GetBinContent(hMeasured->FindBin(matchval))<<" / "<<hMeasuredSim->GetBinContent(hMeasuredSim->FindBin(matchval))<<" = "<<myscale2<<endl;
  GetChi2(hMeasured,hMeasuredSim);
  cout<<"Chi^2 of Smeared compared to measured "<<GetChi2(hMeasured,hMeasuredSim)<<endl;
  hMeasuredSim->Scale(myscale2);
  hMeasuredSim->Draw("same");
  hMeasuredSim->SetMarkerStyle(21);
  hMeasuredSim->SetLineColor(2);
  hMeasuredSim->SetMarkerColor(2);
  TLegend *legend5 = new TLegend(0.437919,0.800699,0.966443,0.958042);
  legend5->AddEntry(hMeasured,"Measured");
  legend5->AddEntry(hMeasuredSim,"Simulated Measured");
  legend5->Draw();
  //canvas5->SaveAs(Form("%s%s.png",outputfilename,"Reco"));
  //canvas5->SaveAs(Form("%s%s.eps",outputfilename,"Reco"));
  canvas5->SaveAs(Form("%s%s.C",outputfilename,"Reco"));


  TCanvas *canvas6 = new TCanvas("RecoOverSimReco","RecoOverSimReco",600,600);
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
  canvas6->SetLogy();

  //TH1F *hMeasuredClone = hMeasured->Clone("hMeasuredClone");
  //hMeasuredClone->Divide(hMeasuredSim);
  Float_t avgratio = 0;
  Int_t lowbin = hMeasured->FindBin(0.0);
  Int_t highbin = hMeasured->FindBin(maxEt);
  TH1F *hMeasuredClone = new TH1F("hMeasuredClone","Measured/simulated measured",highbin-lowbin+1,0,maxEt);
  Int_t mylowbin = hMeasured->FindBin(1.0);
  Int_t myhighbin = hMeasured->FindBin(maxEt-2.0);
  Float_t maxratio = 0;
  Float_t minratio = 1e6;
  for(Int_t i = lowbin;i<=highbin;i++){
    Float_t ratio = hMeasured->GetBinContent(i) / hMeasuredSim->GetBinContent(i);
    Float_t ratioerror = 0;
    if(hMeasured->GetBinContent(i)>0 && hMeasuredSim->GetBinContent(i)>0) ratio * TMath::Sqrt(TMath::Power(hMeasured->GetBinError(i)/hMeasured->GetBinContent(i),2)+TMath::Power(hMeasuredSim->GetBinError(i)/hMeasuredSim->GetBinContent(i),2));
    hMeasuredClone->SetBinContent(i,ratio);
    hMeasuredClone->SetBinError(i,ratioerror);
    if(i>=mylowbin && i<=myhighbin){
      avgratio += ratio;
      if(maxratio<ratio) maxratio = ratio;
      if(minratio>ratio) minratio = ratio;
    }
  }
  avgratio = avgratio/((Float_t)(myhighbin-mylowbin+1));
  cout<<"average ratio : "<<avgratio<<" max ratio "<<maxratio<<" min ratio "<<minratio;
  Bool_t goodFit = kFALSE;
  if(avgratio<5 && avgratio>0.1 && maxratio<2 && minratio>0.1){
    goodFit = kTRUE;
    cout<<" GOOD"<<endl;
  }
  else{
    cout<<" BAD"<<endl;
  }
  hMeasuredClone->Draw();
  canvas6->SaveAs(Form("%s%s.C",outputfilename,"RecoOverSim"));
  canvas6->SaveAs(Form("%s%s.png",outputfilename,"RecoOverSim"));


  TString temp = outputfilename;
  TString filename = temp+".root";
  TFile *outfile2 = new TFile(filename.Data(), "RECREATE");
  cout<<"Creating "<<filename.Data()<<endl;
  hReco1->Write();
  hMeasured->Write();
  hMeasuredSim->Write();
  outfile2->Close();

  TString temp = outputfilename;
  TString filename = temp+".root";
  TFile *outfile2 = new TFile(filename.Data(), "RECREATE");
  hReco1->Write();
  hMeasured->Write();
  hMeasuredSim->Write();
  outfile2->Close();
  return;


//   func->Draw();
  //funcDoubleExponential->Draw();
}

TH1F *GetReconstructedEt(TF1 *inputfunc, int nevents, char *name, float lowbound, float highbound, int nbins,bool smooth){
  TH1F *hResult = new TH1F(name,ytitle->Data(),nbins,lowbound,highbound);
  hResult->Sumw2();
  // cerr<<__FILE__<<" "<<__LINE__<<" Creating "<<hResult->GetName()<<" with "<<nevents<<" events"<<endl;
  for(int i=0;i<nevents;i++){
    float et = inputfunc->GetRandom(lowbound,highbound);
    int mybin = resolutionFull->GetXaxis()->FindBin(et);
    if(mybin>0 && mybin<=nbins){//then this is in our range...
      float myres = ((TH1D*)histoarray.At(mybin-1))->GetRandom();
      //float etreco = (1-myres)*et;
      float etreco = myres;
      //cout<<"et true "<<et<<" reco "<<etreco<<endl;
      hResult->Fill(etreco);
    }
    //if(i%100000==0) Smooth(hResult,inputfunc);
  }
  //if(smooth) 
  //cerr<<__FILE__<<" "<<__LINE__<<" "<<name<<endl;
   //Smooth(hResult,inputfunc);
  //float binwidth = hResult->GetBinWidth(1);
  //hResult->Scale(1.0/nevents/binwidth);
  return hResult;
}
TH1F *GetReconstructedEt(TH1 *input, int nevents, char *name, float lowbound, float highbound, int nbins){
  // cerr<<__FILE__<<" "<<__LINE__<<" "<<name<<endl;
  TH1F *hResult = new TH1F(name,ytitle->Data(),nbins,lowbound,highbound);
  hResult->Sumw2();
  for(int i=0;i<nevents;i++){
    float et = input->GetRandom();
    int mybin = resolutionFull->GetXaxis()->FindBin(et);
    if(mybin>0 && mybin<=nbins){//then this is in our range...
      float myres = ((TH1D*)histoarray.At(mybin-1))->GetRandom();
      //float etreco = (1-myres)*et;
      float etreco = myres;
      //cout<<"et true "<<et<<" reco "<<etreco<<endl;
      hResult->Fill(etreco);
    }
  }
  hResult->GetYaxis()->SetTitle("1/N_{eve} dN/dE_{T}");
  hResult->GetXaxis()->SetTitle("E_{T}");
  //float binwidth = hResult->GetBinWidth(1);
  //hResult->Scale(1.0/nevents/binwidth);
  return hResult;
}

TH1F *GetTrueEt(TF1 *inputfunc, int nevents, char *name, float lowbound, float highbound, int nbins){

  //cerr<<__FILE__<<" "<<__LINE__<<" Throwing simulation "<<name<<" with "<<nevents<<" events "<<endl;
  TH1F *hResult = new TH1F(name,ytitle->Data(),nbins,lowbound,highbound);
  hResult->Sumw2();
  //cerr<<"event ";
  for(int i=0;i<nevents;i++){
    //if(i%1000) cerr<<i<<" ";
    float et = inputfunc->GetRandom(lowbound,highbound);
    hResult->Fill(et);
  }
  //cerr<<endl;
  //hResult->SetBinContent(1,hResult->GetBinContent(1)*10);
  float binwidth = hResult->GetBinWidth(1);
  //hResult->Scale(1.0/nevents/binwidth);
  hResult->GetYaxis()->SetTitle("1/N_{eve} dN/dE_{T}");
  hResult->GetXaxis()->SetTitle("E_{T}");
  return hResult;
}


Float_t GetChi2(TH1F *reconstructed, TH1F *simulated){
  Float_t chi2 = 0;
  int nbins = reconstructed->FindBin(10.0);
  int lowbin =  reconstructed->FindBin(5.01);
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
  cout<<"Chi^2/ndf "<<chi2<<"/"<<ndf<<" = "<<chi2/ndf<<endl;
  if(ndf>0)return chi2/ndf;
  else{return 0;}
}

TH1F *RestrictAxes(TH1F *old,float minEt,float maxEt){
  int minbin = old->FindBin(minEt+0.01);
  int maxbin = old->FindBin(maxEt-0.01);
  //cout<<"Restricting range.  Old bin width:"<<old->GetBinWidth(1);
  int totnbins = maxbin-minbin+1;
  TH1F *hNew = new TH1F(Form("%sNew",old->GetName()),old->GetTitle(),totnbins,minEt,maxEt);
  hNew->GetYaxis()->SetTitle(old->GetYaxis()->GetTitle());
  hNew->GetXaxis()->SetTitle(old->GetXaxis()->GetTitle());
  //cout<<old->GetName()<<" New range "<<minEt<<"-"<<maxEt<<" bins "<<minbin<<"-"<<maxbin<<endl;
  for(int i=minbin;i<=maxbin;i++){
    //cout<<hNew->FindBin(old->GetBinCenter(i))<<":"<<old->GetBinContent(i);
    //if(i!=maxbin) cout<<",";
    int oldbin = i;
    float newbin = hNew->FindBin(old->GetBinCenter(oldbin));
    //cout<<"old "<<oldbin<<" new "<<newbin<<" center "<<hNew->GetBinCenter(newbin)<<"="<<old->GetBinCenter(oldbin)<<" content "<<old->GetBinContent(oldbin)<<" +/- "<<old->GetBinError(oldbin)<<endl;
    hNew->SetBinContent(newbin,old->GetBinContent(oldbin));
    hNew->SetBinError(newbin,old->GetBinError(oldbin));
  }
  //cout<<endl;
  hNew->Sumw2();
  //delete old;
  //old = hNew;
  //cout<<" new bin width "<<hNew->GetBinWidth(1)<<endl;
  return hNew;
}
TH2F *RestrictAxes(TH2F *old,float minEt,float maxEt){
  int minbin = old->GetYaxis()->FindBin(minEt+0.01);
  int maxbin = old->GetYaxis()->FindBin(maxEt-0.01);
  //cout<<"Restricting range.  Old bin width:"<<old->GetBinWidth(1);
  int totnbins = maxbin-minbin+1;
  TH2F *hNew = new TH2F(Form("%sNew",old->GetName()),old->GetTitle(),totnbins,minEt,maxEt,totnbins,minEt,maxEt);
  hNew->GetYaxis()->SetTitle(old->GetYaxis()->GetTitle());
  hNew->GetXaxis()->SetTitle(old->GetXaxis()->GetTitle());
  for(int i=minbin;i<=maxbin;i++){
    int newbinx = hNew->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i));
    int oldbinx = i;//hNew->GetXaxis()->FindBin(old->GetXaxis()->GetBinCenter(i));
    for(int j=minbin;j<=maxbin;j++){
      int newbiny = hNew->GetYaxis()->FindBin(old->GetYaxis()->GetBinCenter(j));
      int oldbiny = j;//hNew->GetYaxis()->FindBin(old->GetYaxis()->GetBinCenter(j));
      //int mybinNew = hNew->FindBin(old->GetXaxis()->GetBinCenter(i),old->GetYaxis()->GetBinCenter(j));
      //int mybinOld = old->FindBin(old->GetXaxis()->GetBinCenter(i),old->GetYaxis()->GetBinCenter(j));
      //cout<<"i "<<i<<" j "<<j<<" ";
      //cout<<"Old bin ("<<oldbinx<<","<<oldbiny<<") is becoming ("<<newbinx<<","<<newbiny<<") with content "<<old->GetBinContent(oldbinx,oldbiny)<<" +/- "<<old->GetBinError(oldbinx,oldbiny)<<endl;
      hNew->SetBinContent(newbinx,newbiny,old->GetBinContent(oldbinx,oldbiny));
      hNew->SetBinError(newbinx,newbiny,old->GetBinError(oldbinx,oldbiny));
      //hNew->SetBinError(mybinNew,old->GetBinError(mybinOld));
    }
  }
  hNew->Sumw2();
  return hNew;
  //delete old;
  //old = hNew;
  //cout<<" new bin width "<<hNew->GetBinWidth(1)<<endl;
}

void Smooth(TH1F *histo,TF1 *func){
  cout<<"Smoothing..."<<endl;
  //histo->Smooth(1,"R");
  //cerr<<__FILE__<<" "<<__LINE__<<endl;
  histo->Fit(func,"N","",10,35);
  //cerr<<__FILE__<<" "<<__LINE__<<endl;
  for(int i=histo->FindBin(15.0);i<histo->GetNbinsX();i++){
    float thisval = histo->GetBinContent(i);
    float thiscenter = histo->GetBinCenter(i);
    float thiserr = histo->GetBinError(i);
    float expectation = func->Eval(thiscenter);
    //if the number is more than 10 sigma away from the expectation and we would have expected at least 1 entry in that bin, recalculate it assuming Gaussian smoothing
    if(thisval<1.0) break;
    //cerr<<__FILE__<<" "<<__LINE__<<endl;
    if(expectation>1.0 && TMath::Abs(thiserr)>1e-5 && TMath::Abs(thisval-expectation)/thiserr>5){
      float nSigma = myGaussian->GetRandom();
      //cerr<<__FILE__<<" "<<__LINE__<<endl;
      histo->SetBinContent(i,expectation+nSigma*TMath::Sqrt(expectation));
      //cout<<"Resetting bin "<<i<<" at "<<histo->GetBinCenter(i)<<" from "<<thisval<<" which was "<<TMath::Abs(thisval-expectation)/thiserr<<" sigma away to "<<expectation+nSigma*TMath::Sqrt(expectation)<<" nSigma "<<nSigma<<" away from "<<expectation<<endl;
    }
//     float nextval = histo->GetBinContent(i+1);
//     float lastval = histo->GetBinContent(i-1);
//     float nextcenter = histo->GetBinCenter(i+1);
//     float lastcenter = histo->GetBinCenter(i-1);
//     float nexterr = histo->GetBinError(i+1);
//     float lasterr = histo->GetBinError(i-1);
//     if(nextval>thisval){//We have a problem!  This function should decrease!
//       //impose y=A e(-x/b) and use previous two values to calculate b and A
//       float b = TMath::Log(thisval/lastval)/(lastcenter-thiscenter);
//       float A = thisval/exp(-b*thiscenter);
//       float newval = A*exp(-b*nextcenter);
//       cout<<"A "<<A<<" b "<<" center "<<nextcenter<<" old "<<histo->GetBinContent(i+1)<<" new "<<newval<<endl;
//       //if(newval)histo->SetBinContent(i+1,newval);
//     }
  }
}
