#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TInterpreter.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TF1.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TList.h>
#include <TLegendEntry.h>
#include <TDatabasePDG.h>
#include <TGraph.h>

#include "AliAODRecoDecayHF.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliHFMassFitter.h"
#include "AliNormalizationCounter.h"
#endif

/* $Id$ */ 

enum {kD0toKpi, kDplusKpipi};
enum {kBoth, kParticleOnly, kAntiParticleOnly};
enum {kExpo=0, kLinear, kPol2};
enum {kGaus=0, kDoubleGaus};
enum {kCorr=0, kUnCorr, kNoPid};


// Common variables: to be configured by the user
// const Int_t nPtBins=6;
// Double_t ptlims[nPtBins+1]={1., 2.,4.,6.,8.,12.,24.};
// Int_t rebin[nPtBins]={4,4,6,6,8,8};
// Double_t sigmapt[nPtBins]={ 0.008, 0.010, 0.012, 0.016, 0.018, 0.020 };
const Int_t nPtBins=5;
Double_t ptlims[nPtBins+1]={1.,2.,4.,8.,12.,24.};
Int_t rebin[nPtBins]={4,4,6,8,8};
Double_t sigmapt[nPtBins]={ 0.008, 0.014, 0.019, 0.027, 0.033 };
Bool_t fixPeakSigma = kTRUE;
//
const Int_t nMultbins=7;
Double_t multlims[nMultbins+1]={1.,9.,14.,20.,31.,50.,81.,100.};
// const Int_t nMultbins=1;
// Double_t multlims[nMultbins+1]={0.,500.};
//
Int_t firstUsedBin[nPtBins]={-1,-1,-1,-1,-1}; // -1 uses all bins, >=1 to set the lower bin to be accepted from original histo
//
Bool_t isMC=kFALSE;
Int_t typeb=kExpo;
Int_t types=kGaus;
Int_t optPartAntiPart=kBoth;
Int_t factor4refl=0;
Float_t massRangeForCounting=0.05; // GeV --> it is 3 sigmapt[binpt]
Float_t nSigmaRangeForCounting=3.0; //  3 sigmapt[binpt]
TH2F* hPtMass=0x0;
TString suffix="StdPid";


// Functions
Bool_t LoadDplusHistos(TObjArray* listFiles, TH3F** hPtMassMult, TH2F** hNtrZvtx, TH2F** hNtrZvtxCorr, const char *CutsType, Int_t Option);
Bool_t LoadD0toKpiHistos(TObjArray* listFiles, TH3F** hPtMassMult, TH2F** hNtrZvtx, TH2F** hNtrZvtxCorr, AliNormalizationCounter *counter, const char *CutsType, Int_t Option);
Bool_t CheckNtrVsZvtx(TH2F** hNtrZvtx, TH2F** hNtrZvtxCorr, Int_t nFiles);
TH1F* RebinHisto(TH1F* hOrig, Int_t reb, Int_t firstUse=-1);


void ReadDvsMultiplicity(Int_t analysisType=kD0toKpi,
			 TString fileNameb="AnalysisResults.root",
			 TString fileNamec="",
			 TString fileNamed="",
			 TString fileNamee="",
			 const char *CutsType="",
			 Int_t Option=kCorr)
{
  // gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWGHF -I$ALICE_ROOT/PWGHF/base -I$ALICE_ROOT/PWGHF/vertexingHF -I$ALICE_ROOT/PWG/FLOW/Case -I$ALICE_ROOT/PWG/FLOW/Tasks -g");

  // gInterpreter->ExecuteMacro("$ALICE_ROOT/PWGHF/vertexingHF/macros/LoadLibraries.C");
  gStyle->SetOptTitle(1);

  TString ntrkname="Ntracklets";
   // for(int j=0; j<=nMultbins; j++) multlims[j] *= (68./8.8);
   // ntrkname="Nvzero";

  Int_t nFiles=0;
  TObjArray* listFiles=new TObjArray();
  if(fileNameb!="") { listFiles->AddLast(new TObjString(fileNameb.Data())); nFiles++; }
  if(fileNamec!="") { listFiles->AddLast(new TObjString(fileNamec.Data())); nFiles++; }
  if(fileNamed!="") { listFiles->AddLast(new TObjString(fileNamed.Data())); nFiles++; }
  if(fileNamee!="") { listFiles->AddLast(new TObjString(fileNamee.Data())); nFiles++; }
  if(listFiles->GetEntries()==0){
    printf("Missing file names in input\n");
    return;
  }
  TH3F** hPtMassMult=new TH3F*[1];
  hPtMassMult[0]=0x0;
  TH1F** hmass=new TH1F*[nPtBins*nMultbins];
  for(Int_t i=0; i<nPtBins*nMultbins; i++) hmass[i]=0x0;
  TH2F** hNtrZvtx=new TH2F*[4];
  for(Int_t i=0; i<4; i++) hNtrZvtx[i]=0x0;
  TH2F** hNtrZvtxCorr=new TH2F*[4];
  for(Int_t i=0; i<4; i++) hNtrZvtxCorr[i]=0x0;
  AliNormalizationCounter *counter=0x0;
  TH1F * hNormalization = new TH1F("hNormalization","Events in the norm counter",nMultbins,multlims);

  Float_t massD;
  Bool_t retCode;
  if(analysisType==kD0toKpi) {
    massD=1.86484003067016602e+00;//(Float_t)(TDatabasePDG::Instance()->GetParticle(421)->Mass());
    retCode=LoadD0toKpiHistos(listFiles,hPtMassMult,hNtrZvtx,hNtrZvtxCorr,counter,CutsType,Option);
  }
  else if(analysisType==kDplusKpipi) {
    massD=1.86961996555328369e+00;//(Float_t)(TDatabasePDG::Instance()->GetParticle(411)->Mass());
    retCode=LoadDplusHistos(listFiles,hPtMassMult,hNtrZvtx,hNtrZvtxCorr,CutsType,Option);
  }
  else{
    printf("Wrong analysis type parameter\n");
    return;
  }
  if(!retCode){
    printf("ERROR in reading input files\n");
    return;
  }
  //
  // Check the multiplicity variables
  //
  Bool_t isMult = CheckNtrVsZvtx(hNtrZvtx,hNtrZvtxCorr,nFiles);
  if(!isMult) return;
  //
  //
  printf(" Preparing the mass fits");
  TH1F* hCntSig1[nMultbins];
  TH1F* hCntSig2[nMultbins];
  TH1F* hNDiffCntSig1[nMultbins];
  TH1F* hNDiffCntSig2[nMultbins];
  TH1F* hSignal[nMultbins];
  TH1F* hRelErrSig[nMultbins];
  TH1F* hInvSignif[nMultbins];
  TH1F* hBackground[nMultbins];
  TH1F* hBackgroundNormSigma[nMultbins];
  TH1F* hSignificance[nMultbins];
  TH1F* hMass[nMultbins];
  TH1F* hSigma[nMultbins];
  for(Int_t i=0; i<nMultbins; i++){
    hCntSig1[i]=new TH1F(Form("hCntSig1_%d",i),Form("hCntSig1_%d",i),nPtBins,ptlims);
    hCntSig2[i]=new TH1F(Form("hCntSig2_%d",i),Form("hCntSig2_%d",i),nPtBins,ptlims);
    hNDiffCntSig1[i]=new TH1F(Form("hNDiffCntSig1_%d",i),Form("hNDiffCntSig1_%d",i),nPtBins,ptlims);
    hNDiffCntSig2[i]=new TH1F(Form("hNDiffCntSig2_%d",i),Form("hNDiffCntSig2_%d",i),nPtBins,ptlims);
    hSignal[i]=new TH1F(Form("hSignal_%d",i),Form("hSignal_%d",i),nPtBins,ptlims);
    hRelErrSig[i]=new TH1F(Form("hRelErrSig_%d",i),Form("hRelErrSig_%d",i),nPtBins,ptlims);
    hInvSignif[i]=new TH1F(Form("hInvSignif_%d",i),Form("hInvSignif_%d",i),nPtBins,ptlims);
    hBackground[i]=new TH1F(Form("hBackground_%d",i),Form("hBackground_%d",i),nPtBins,ptlims);
    hBackgroundNormSigma[i]=new TH1F(Form("hBackgroundNormSigma_%d",i),Form("hBackgroundNormSigma_%d",i),nPtBins,ptlims);
    hSignificance[i]=new TH1F(Form("hSignificance_%d",i),Form("hSignificance_%d",i),nPtBins,ptlims);
    hMass[i]=new TH1F(Form("hMass_%d",i),Form("hMass_%d",i),nPtBins,ptlims);
    hSigma[i]=new TH1F(Form("hSigma_%d",i),Form("hSigma_%d",i),nPtBins,ptlims);
  }
  printf(", defined...\n");

  //  std::cout << " htemp :"<<hPtMassMult[0]<<std::endl;
  TH1F* hptaxis = (TH1F*)hPtMassMult[0]->ProjectionZ("hptaxis");
  TH1F* hmassaxis = (TH1F*)hPtMassMult[0]->ProjectionY("hmassaxis");
  TH1F* hmultaxis = (TH1F*)hPtMassMult[0]->ProjectionX("hmultaxis");
  Int_t nMassBins=hmassaxis->GetNbinsX();
  Double_t hmin=hmassaxis->GetBinLowEdge(3);
  Double_t hmax=hmassaxis->GetBinLowEdge(nMassBins-2) + hmassaxis->GetBinWidth(nMassBins-2);
  Int_t iPad=1;
  
  printf("Now initializing the fit functions\n");
  // TF1* funBckStore1=0x0;
  // TF1* funBckStore2=0x0;
  // TF1* funBckStore3=0x0;

  Int_t nPtMultbins = nPtBins*nMultbins;
  AliHFMassFitter** fitter=new AliHFMassFitter*[nPtMultbins];
  Double_t arrchisquare0[nPtBins];
  Double_t arrchisquare1[nPtBins];
  Double_t arrchisquare2[nPtBins];
  Double_t arrchisquare3[nPtBins];
  Double_t arrchisquare4[nPtBins];
  Double_t arrchisquare5[nPtBins];
  for(Int_t i=0; i<nPtBins; i++){
    arrchisquare0[i]=0.;
    arrchisquare1[i]=0.;
    arrchisquare2[i]=0.;
    arrchisquare3[i]=0.;
    arrchisquare4[i]=0.;
    arrchisquare5[i]=0.;
  }
  
  TCanvas* canvas[nMultbins];
  Int_t nx=2, ny=2;
  if(nPtBins>4){
    nx=3;
    ny=2;
  }
  if(nPtBins>6){
    nx=4;
    ny=3;
  }
  if(nPtBins>12){
    nx=5;
    ny=4;
  }
  for(Int_t i=0; i<nMultbins; i++ ){
    canvas[i] = new TCanvas(Form("c%d",i),Form("summary canvas for mult bin %d",i));
    canvas[i]->Divide(nx,ny);
  }
  TCanvas *myCanvas[nPtMultbins];
  
  //
  //
  // Loop on multiplicity bins
  //
  Int_t massBin=0;
  for(Int_t j=0; j<nMultbins; j++){
    Double_t sig,errsig,s,errs,b,errb;
    //    printf(" Studying multiplicity bin %d\n",j);
    Int_t multbinlow = hmultaxis->FindBin(multlims[j]);
    Int_t multbinhigh = hmultaxis->FindBin(multlims[j+1])-1;
    Float_t val = multbinlow + (multbinhigh-multbinlow)/2.;
    Int_t hnbin = hNormalization->FindBin(val);
    Float_t nevents = 0.;
    if(counter) { nevents = counter->GetNEventsForNorm(multbinlow,multbinhigh); std::cout << std::endl<<std::endl<<" Nevents ("<<multbinlow<<","<<multbinhigh<<") ="<<nevents<<std::endl<<std::endl<<std::endl;}
    hNormalization->SetBinContent(hnbin,nevents);
    //
    // Loop on pt bins
    //
    iPad=1;
    for(Int_t iBin=0; iBin<nPtBins; iBin++){
      canvas[j]->cd(iPad++);
      //      printf(" projecting to the mass histogram\n");
      Int_t ptbinlow = hptaxis->FindBin(ptlims[iBin]);
      Int_t ptbinhigh = hptaxis->FindBin(ptlims[iBin+1])-1;
      hmass[massBin] = (TH1F*)hPtMassMult[0]->ProjectionY(Form("hmass_%d_%d",j,iBin),multbinlow,multbinhigh,ptbinlow,ptbinhigh);
      hmass[massBin]->SetTitle( Form("%2.0f<p_{T}<%2.0f GeV/c, %s [%3.0f,%3.0f]",ptlims[iBin],ptlims[iBin+1],ntrkname.Data(),multlims[j],multlims[j+1]) );
      //      std::cout << std::endl<<  Form("%2.0f<p_{T}<%2.0f GeV/c, %s [%3.0f,%3.0f]",ptlims[iBin],ptlims[iBin+1],ntrkname.Data(),multlims[j],multlims[j+1]) << std::endl<< std::endl;
      if(  hmass[massBin]->GetEntries() < 60 ) {
	massBin++;
	continue;
      }
      Int_t origNbins=hmass[massBin]->GetNbinsX(); 
      //      printf(" rebinning the mass histogram\n");
      TH1F* hRebinned=RebinHisto(hmass[massBin],rebin[iBin],firstUsedBin[iBin]);
      hmin=hRebinned->GetBinLowEdge(2);
      hmax=hRebinned->GetBinLowEdge(hRebinned->GetNbinsX());
      //      printf(" define the mass fitter and options \n");
      fitter[massBin] = new AliHFMassFitter(hRebinned,hmin, hmax,1,typeb,types);
      fitter[massBin]->SetRangeFit(1.65,2.10);
      Int_t rebinItem = origNbins/fitter[massBin]->GetBinN();
      fitter[massBin]->SetReflectionSigmaFactor(factor4refl);
      fitter[massBin]->SetInitialGaussianMean(massD);
      fitter[massBin]->SetInitialGaussianSigma(sigmapt[iBin]);
      if(fixPeakSigma) {
	fitter[massBin]->SetFixGaussianMean(massD);
	fitter[massBin]->SetFixGaussianSigma(sigmapt[iBin],kTRUE);
      }
      Bool_t out=fitter[massBin]->MassFitter(0);
      if(!out) {
	fitter[massBin]->GetHistoClone()->Draw();
	massBin++;
	continue;
      }
      //      printf(" getting the fit parameters\n");
      Double_t mass=fitter[massBin]->GetMean();
      Double_t massUnc=fitter[massBin]->GetMeanUncertainty();
      Double_t sigma=fitter[massBin]->GetSigma();
      Double_t sigmaUnc=fitter[massBin]->GetSigmaUncertainty();
      if(j==0) arrchisquare0[iBin]=fitter[massBin]->GetReducedChiSquare();
      else if(j==1) arrchisquare1[iBin]=fitter[massBin]->GetReducedChiSquare();
      else if(j==2) arrchisquare2[iBin]=fitter[massBin]->GetReducedChiSquare();
      else if(j==3) arrchisquare3[iBin]=fitter[massBin]->GetReducedChiSquare();
      else if(j==4) arrchisquare4[iBin]=fitter[massBin]->GetReducedChiSquare();
      else if(j==5) arrchisquare5[iBin]=fitter[massBin]->GetReducedChiSquare();
      TF1* fB1=fitter[massBin]->GetBackgroundFullRangeFunc();
      TF1* fB2=fitter[massBin]->GetBackgroundRecalcFunc();
      //      TF1* fM=fitter[massBin]->GetMassFunc();
      // if(iBin==0 && fB1) funBckStore1=new TF1(*fB1);
      // if(iBin==0 && fB2) funBckStore2=new TF1(*fB2);
      // if(iBin==0 && fM) funBckStore3=new TF1(*fM);

      fitter[massBin]->DrawHere(gPad);
      fitter[massBin]->Signal(3,s,errs);
      fitter[massBin]->Background(3,b,errb);
      fitter[massBin]->Significance(3,sig,errsig);
      Double_t ry=fitter[massBin]->GetRawYield();
      Double_t ery=fitter[massBin]->GetRawYieldError();
      myCanvas[massBin] = new TCanvas(Form("myCanvas_%d_%d",j,iBin),Form("Invariant mass mult bin %d, pt bin %d",j,iBin));
      fitter[massBin]->DrawHere(gPad);
    
      Float_t cntSig1=0.;
      Float_t cntSig2=0.;
      Float_t cntErr=0.;
      massRangeForCounting = nSigmaRangeForCounting*sigmapt[iBin];
      //      std::cout << " pt bin "<< iBin << " mass range = "<< massRangeForCounting<<std::endl;
      Float_t minBinSum=hmassaxis->FindBin(massD-massRangeForCounting);
      Float_t maxBinSum=hmassaxis->FindBin(massD+massRangeForCounting);
      for(Int_t iMB=minBinSum; iMB<=maxBinSum; iMB++){
	Float_t bkg1=fB1 ? fB1->Eval(hmass[massBin]->GetBinCenter(iMB))/rebinItem : 0;
	Float_t bkg2=fB2 ? fB2->Eval(hmass[massBin]->GetBinCenter(iMB))/rebinItem : 0;
	cntSig1+=(hmass[massBin]->GetBinContent(iMB)-bkg1);
	cntSig2+=(hmass[massBin]->GetBinContent(iMB)-bkg2);
	cntErr+=(hmass[massBin]->GetBinContent(iMB));
      }
      hCntSig1[j]->SetBinContent(iBin+1,cntSig1);
      hCntSig1[j]->SetBinError(iBin+1,TMath::Sqrt(cntErr));
      hNDiffCntSig1[j]->SetBinContent(iBin+1,(s-cntSig1)/s);
      hNDiffCntSig1[j]->SetBinError(iBin+1,TMath::Sqrt(cntErr)/s);
      hCntSig2[j]->SetBinContent(iBin+1,cntSig2);
      hNDiffCntSig2[j]->SetBinContent(iBin+1,(s-cntSig2)/s);
      hNDiffCntSig2[j]->SetBinError(iBin+1,TMath::Sqrt(cntErr)/s);
      hCntSig2[j]->SetBinError(iBin+1,TMath::Sqrt(cntErr));
      hSignal[j]->SetBinContent(iBin+1,ry);
      hSignal[j]->SetBinError(iBin+1,ery);
      hRelErrSig[j]->SetBinContent(iBin+1,errs/s);
      hInvSignif[j]->SetBinContent(iBin+1,1/sig);
      hInvSignif[j]->SetBinError(iBin+1,errsig/(sig*sig));
      hBackground[j]->SetBinContent(iBin+1,b); //consider sigma
      hBackground[j]->SetBinError(iBin+1,errb);
      hBackgroundNormSigma[j]->SetBinContent(iBin+1,b/(3*fitter[massBin]->GetSigma())*(3*0.012)); //consider sigma
      hBackgroundNormSigma[j]->SetBinError(iBin+1,errb);
      hSignificance[j]->SetBinContent(iBin+1,sig);
      hSignificance[j]->SetBinError(iBin+1,errsig);
      hMass[j]->SetBinContent(iBin+1,mass);
      hMass[j]->SetBinError(iBin+1,massUnc);
      hSigma[j]->SetBinContent(iBin+1,sigma);
      hSigma[j]->SetBinError(iBin+1,sigmaUnc);

      massBin++;
      delete hRebinned;
    }// end loop on pt bins

    canvas[j]->Update();
    canvas[j]->SaveAs(Form("hMass%s_%d_%d.eps",CutsType,typeb,j));
    //    canvas[j]->SaveAs(Form("hMass%s_%d_%d_MultInt.eps",CutsType,typeb,j));
    
  }// end loop on multiplicity bins


  TCanvas *cpar=new TCanvas("cpar","Fit params",1200,600);
  cpar->Divide(2,1);
  cpar->cd(1);
  for(Int_t imult=0; imult<nMultbins; imult++) {
    hMass[imult]->SetMarkerStyle(20);
    hMass[imult]->GetXaxis()->SetTitle("Pt (GeV/c)");
    hMass[imult]->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
    hMass[imult]->SetMarkerColor(2*imult);
    if(imult==5) hMass[imult]->SetMarkerColor(2*imult-3);
    if(imult==0) {
      hMass[imult]->SetMarkerColor(kBlack);
      hMass[imult]->Draw("PE");
    }
    else  hMass[imult]->Draw("PEsame");
  }
  cpar->cd(2);
  for(Int_t imult=0; imult<nMultbins; imult++) {
    hSigma[imult]->SetMarkerStyle(20);
    //  hSigma[0]->Draw("PE");
    hSigma[imult]->GetXaxis()->SetTitle("Pt (GeV/c)");
    hSigma[imult]->GetYaxis()->SetTitle("Sigma (GeV/c^{2})");
    hSigma[imult]->SetMarkerColor(2*imult);
    if(imult==5) hSigma[imult]->SetMarkerColor(2*imult-3);
    if(imult==0) {
      hSigma[imult]->SetMarkerColor(kBlack);
      hSigma[imult]->Draw("PE");
    }
    else  hSigma[imult]->Draw("PEsame");
  }
  cpar->Update();

  /*
  TCanvas** csig;//= new TCanvas*[nMultbins];
  TCanvas** cDiffS;//=new TCanvas*[nMultbins];
  for(Int_t i=0; i<nMultbins; i++){
    csig[i] =new TCanvas(Form("csig_%d",i),Form("Results, multiplicity bin %d",i),1200,600);
    csig[i]->Divide(3,1);
    csig[i]->cd(1);
    hSignal[i]->SetMarkerStyle(20);
    hSignal[i]->SetMarkerColor(4);
    hSignal[i]->SetLineColor(4);
    hSignal[i]->GetXaxis()->SetTitle("Pt (GeV/c)");
    hSignal[i]->GetYaxis()->SetTitle("Signal");
    hSignal[i]->Draw("P");
    hCntSig1[i]->SetMarkerStyle(26);
    hCntSig1[i]->SetMarkerColor(2);
    hCntSig1[i]->SetLineColor(2);
    hCntSig1[i]->Draw("PSAME");
    hCntSig2[i]->SetMarkerStyle(29);
    hCntSig2[i]->SetMarkerColor(kGray+1);
    hCntSig2[i]->SetLineColor(kGray+1);
    hCntSig2[i]->Draw("PSAME");
    TLegend* leg=new TLegend(0.4,0.7,0.89,0.89);
    leg->SetFillColor(0);
    TLegendEntry* ent=leg->AddEntry(hSignal[i],"From Fit","PL");
    ent->SetTextColor(hSignal[i]->GetMarkerColor());
    ent=leg->AddEntry(hCntSig1[i],"From Counting1","PL");
    ent->SetTextColor(hCntSig1[i]->GetMarkerColor());
    ent=leg->AddEntry(hCntSig2[i],"From Counting2","PL");
    ent->SetTextColor(hCntSig2[i]->GetMarkerColor());
    leg->Draw();
    csig[i]->cd(2);
    hBackground[i]->SetMarkerStyle(20);
    hBackground[i]->Draw("P");
    hBackground[i]->GetXaxis()->SetTitle("Pt (GeV/c)");
    hBackground[i]->GetYaxis()->SetTitle("Background");
    csig[i]->cd(3);
    hSignificance[i]->SetMarkerStyle(20);
    hSignificance[i]->Draw("P");
    hSignificance[i]->GetXaxis()->SetTitle("Pt (GeV/c)");
    hSignificance[i]->GetYaxis()->SetTitle("Significance");
    
    cDiffS[i]=new TCanvas(Form("cDiffS_%d",i),Form("Difference, multiplicity bin %d",i),1200,600);
    cDiffS[i]->Divide(2,1);
    cDiffS[i]->cd(1);
    hRelErrSig[i]->SetMarkerStyle(20); //fullcircle
    hRelErrSig[i]->SetTitleOffset(1.2);  
    hRelErrSig[i]->SetTitle("Rel Error from Fit;p_{t} (GeV/c);Signal Relative Error");
    hRelErrSig[i]->Draw("P");
    hInvSignif[i]->SetMarkerStyle(21); //fullsquare
    hInvSignif[i]->SetMarkerColor(kMagenta+1);
    hInvSignif[i]->SetLineColor(kMagenta+1);
    hInvSignif[i]->Draw("PSAMES");
    TLegend* leg2=new TLegend(0.4,0.7,0.89,0.89);
    leg2->SetFillColor(0);
    TLegendEntry* ent2=leg2->AddEntry(hRelErrSig[i],"From Fit","P");
    ent2->SetTextColor(hRelErrSig[i]->GetMarkerColor());
    ent2=leg2->AddEntry(hInvSignif[i],"1/Significance","PL");
    ent2->SetTextColor(hInvSignif[i]->GetMarkerColor());
    leg2->Draw();

    cDiffS[i]->cd(2);
    hNDiffCntSig1[i]->SetMarkerStyle(26);
    hNDiffCntSig1[i]->SetMarkerColor(2);
    hNDiffCntSig1[i]->SetLineColor(2);
    hNDiffCntSig1[i]->SetTitle("Cmp Fit-Count;p_{t} (GeV/c);(S_{fit}-S_{count})/S_{fit}");
    hNDiffCntSig1[i]->Draw("P");
    hNDiffCntSig2[i]->SetMarkerStyle(29);
    hNDiffCntSig2[i]->SetMarkerColor(kGray+1);
    hNDiffCntSig2[i]->SetLineColor(kGray+1);
    hNDiffCntSig2[i]->Draw("PSAME");
    TLegend* leg1=new TLegend(0.4,0.7,0.89,0.89);
    leg1->SetFillColor(0);
    TLegendEntry* ent1=leg1->AddEntry(hNDiffCntSig1[i],"From Counting1","PL");
    ent1->SetTextColor(hNDiffCntSig1[i]->GetMarkerColor());
    ent1=leg1->AddEntry(hNDiffCntSig2[i],"From Counting2","PL");
    ent1->SetTextColor(hNDiffCntSig2[i]->GetMarkerColor());
    leg1->Draw();
  }
*/

  TGraph* grReducedChiSquare0=new TGraph(nPtBins,ptlims,arrchisquare0);
  grReducedChiSquare0->SetName("grReducedChiSquare0");
  grReducedChiSquare0->SetTitle("Reduced Chi2;p_{t} (GeV/c);#tilde{#chi}^{2}");
  TGraph* grReducedChiSquare1=new TGraph(nPtBins,ptlims,arrchisquare1);
  grReducedChiSquare1->SetName("grReducedChiSquare1");
  grReducedChiSquare1->SetTitle("Reduced Chi2;p_{t} (GeV/c);#tilde{#chi}^{2}");
  TGraph* grReducedChiSquare2=new TGraph(nPtBins,ptlims,arrchisquare2);
  grReducedChiSquare2->SetName("grReducedChiSquare2");
  grReducedChiSquare2->SetTitle("Reduced Chi2;p_{t} (GeV/c);#tilde{#chi}^{2}");
  TGraph* grReducedChiSquare3=new TGraph(nPtBins,ptlims,arrchisquare3);
  grReducedChiSquare3->SetName("grReducedChiSquare3");
  grReducedChiSquare3->SetTitle("Reduced Chi2;p_{t} (GeV/c);#tilde{#chi}^{2}");
  TCanvas *cChi2=new TCanvas("cChi2","reduced chi square",600,600);
  cChi2->cd();
  grReducedChiSquare0->SetMarkerStyle(21);
  grReducedChiSquare0->Draw("AP");
  grReducedChiSquare1->SetMarkerStyle(22);
  grReducedChiSquare1->Draw("Psame");
  grReducedChiSquare2->SetMarkerStyle(23);
  grReducedChiSquare2->Draw("Psame");
  grReducedChiSquare3->SetMarkerStyle(24);
  grReducedChiSquare3->Draw("Psame");

  TCanvas* cbkgNormSigma=new TCanvas("cbkgNormSigma","Background normalized to sigma",400,600);
  cbkgNormSigma->cd();
  for(Int_t i=0; i<nMultbins; i++){
    hBackgroundNormSigma[i]->SetMarkerStyle(20);
    hBackgroundNormSigma[i]->GetXaxis()->SetTitle("Pt (GeV/c)");
    hBackgroundNormSigma[i]->GetYaxis()->SetTitle("Background #times 3 #times 0.012/ (3 #times #sigma)");
    hBackgroundNormSigma[i]->SetMarkerColor(2*i);
    if(i==5) hBackgroundNormSigma[i]->SetMarkerColor(2*i-3);
    if(i==0) { 
      hBackgroundNormSigma[i]->SetMarkerColor(kBlack);
      hBackgroundNormSigma[i]->Draw("PE");
    }
    else  hBackgroundNormSigma[i]->Draw("Psame");
  }


  TString partname="Both";
  if(optPartAntiPart==kParticleOnly) {
    if(analysisType==kD0toKpi) partname="D0";
    if(analysisType==kDplusKpipi) partname="Dplus";
  }
  if(optPartAntiPart==kAntiParticleOnly) {
    if(analysisType==kD0toKpi) partname="D0bar";
    if(analysisType==kDplusKpipi) partname="Dminus";
  }

  TString outfilename = Form("RawYield_Mult_%s_%s",partname.Data(),CutsType);
  //  outfilename += "_MultInt";
  if(fixPeakSigma) outfilename += "_SigmaFixed";
  outfilename += Form("_BCin%1.1fSigma",nSigmaRangeForCounting);
  if(typeb==0) outfilename += "_Expo.root";
  else if(typeb==1) outfilename += "_Linear.root";
  else if(typeb==2) outfilename += "_Pol2.root";

  TFile* outf=new TFile(outfilename,"recreate");
  outf->cd(); 
  hNormalization->Write();
  for(Int_t j=0; j<massBin; j++) hmass[j]->Write();
  for(Int_t j=0; j<nMultbins; j++){
    hMass[j]->Write();
    hSigma[j]->Write();
    hCntSig1[j]->Write();
    hCntSig2[j]->Write();
    hNDiffCntSig1[j]->Write();
    hNDiffCntSig2[j]->Write();
    hSignal[j]->Write();
    hRelErrSig[j]->Write();
    hInvSignif[j]->Write();
    hBackground[j]->Write();
    hBackgroundNormSigma[j]->Write();
    hSignificance[j]->Write();
  }
  grReducedChiSquare0->Write();
  grReducedChiSquare1->Write();
  grReducedChiSquare2->Write();
  grReducedChiSquare3->Write();
  //  hPtMass->Write();
  outf->Close();
  
}

//_____________________________________________________________________________________________
Bool_t LoadD0toKpiHistos(TObjArray* listFiles, TH3F** hPtMassMult, TH2F** hNtrZvtx, TH2F** hNtrZvtxCorr, AliNormalizationCounter *counter, const char *CutsType, Int_t Option)
{
  Int_t nFiles=listFiles->GetEntries();
  TList **hlist=new TList*[nFiles];
  TList **hlistNorm=new TList*[nFiles];
  AliRDHFCutsD0toKpi** dcuts=new AliRDHFCutsD0toKpi*[nFiles];

  Int_t nReadFiles=0;
  for(Int_t iFile=0; iFile<nFiles; iFile++){
    TString fName=((TObjString*)listFiles->At(iFile))->GetString();    
    TFile *f=TFile::Open(fName.Data());
    if(!f){
      printf("ERROR: file %s does not exist\n",fName.Data());
      continue;
    }
    printf("Open File %s\n",f->GetName());

    TString dirname="PWG3_D2H_DMult_D0";
    if(optPartAntiPart==kParticleOnly) dirname+="D0";
    else if(optPartAntiPart==kAntiParticleOnly) dirname+="D0bar";
    dirname += CutsType;
    TDirectory *dir = (TDirectory*)f->Get(dirname);
    if(!dir){
      printf("ERROR: directory %s not found in %s\n",dirname.Data(),fName.Data());
      continue;
    }

    TString listmassname="coutputD0";
    if(optPartAntiPart==kParticleOnly) listmassname+="D0";
    else if(optPartAntiPart==kAntiParticleOnly) listmassname+="D0bar";
    listmassname += CutsType;
    printf("List mass name %s\n",listmassname.Data());
    hlist[nReadFiles]=(TList*)dir->Get(listmassname);

    TString listnorm="coutputNormD0";
    if(optPartAntiPart==kParticleOnly) listnorm+="D0";
    else if(optPartAntiPart==kAntiParticleOnly) listnorm+="D0bar";
    listnorm += CutsType;
    printf("List norm name %s\n",listnorm.Data());
    hlistNorm[nReadFiles]=(TList*)dir->Get(listnorm);
    //    AliNormalizationCounter *tmpcounter = (AliNormalizationCounter*)hlistNorm[nReadFiles]->FindObject("NormCounterCorrMult");
    //    counter->Add(tmpcounter);
    //    delete tmpcounter;
    //    counter = tmpcounter;

    TString cutsobjname="coutputCutsD0";
    if(optPartAntiPart==kParticleOnly) cutsobjname+="D0";
    else if(optPartAntiPart==kAntiParticleOnly) cutsobjname+="D0bar";
    cutsobjname += CutsType;
    printf("Cuts name %s\n",cutsobjname.Data());
    dcuts[nReadFiles]=(AliRDHFCutsD0toKpi*)dir->Get(cutsobjname);
    if(!dcuts[nReadFiles]) {
      printf("ERROR: Cut objects do not match\n");
      return kFALSE;
    }
    /*
    if(nReadFiles>0){
      Bool_t sameCuts=dcuts[nReadFiles]->CompareCuts(dcuts[0]);
      if(!sameCuts){
	printf("ERROR: Cut objects do not match\n");
	return kFALSE;
      }
    }
    */
    nReadFiles++;
  }
  if(nReadFiles<nFiles){
    printf("WARNING: not all requested files have been found\n");
    if (nReadFiles==0) {
      printf("ERROR: Any file/dir found\n");
      return kFALSE;
    }
  }
  //  printf("Cuts type %s, particle/antipart %d\n",CutsType,optPartAntiPart);

  /*
  Int_t nPtBinsCuts=dcuts[0]->GetNPtBins();
  printf("Number of pt bins for cut object = %d\n",nPtBins);
  Float_t *ptlimsCuts=dcuts[0]->GetPtBinLimits();
  ptlimsCuts[nPtBinsCuts]=ptlimsCuts[nPtBinsCuts-1]+4.;
  */

  printf("Get the 3D histogram \n");
  const char *histoname="";
  if(optPartAntiPart==kParticleOnly) histoname= "hPtVsMassvsMultPart";
  else if(optPartAntiPart==kAntiParticleOnly) histoname="hPtVsMassvsMultAntiPart";
  else if(optPartAntiPart==kBoth) histoname="hPtVsMassvsMult";
  if(Option==kUnCorr) histoname="hPtVsMassvsMultUncorr";
  if(Option==kNoPid) histoname="hPtVsMassvsMultNoPid";

  TH3F * htemp;
  for(Int_t iFile=0; iFile<nReadFiles; iFile++){
    printf(" Looking for histo histMass %s for file %d\n",histoname,iFile);
    htemp=(TH3F*)hlist[iFile]->FindObject(Form("%s",histoname));
    //    cout << " htemp :"<<htemp<<endl;
    if(!hPtMassMult[0]){
      hPtMassMult[0]=new TH3F(*htemp);
    }else{
      hPtMassMult[0]->Add(htemp);
    }
    hNtrZvtx[iFile] = (TH2F*)hlist[iFile]->FindObject("hNtrVsZvtx");
    hNtrZvtxCorr[iFile] = (TH2F*)hlist[iFile]->FindObject("hNtrCorrVsZvtx");
  }
  
  //  std::cout<<" hPtMassMult:"<<hPtMassMult[0]<<std::endl;

  TString partname="Both";
  if(optPartAntiPart==kParticleOnly) partname="D0";
  if(optPartAntiPart==kAntiParticleOnly) partname="D0bar";

  TString outfilename = Form("RawYield%s_%s",partname.Data(),CutsType);
  if(fixPeakSigma) outfilename += "_SigmaFixed";
  if(typeb==0) outfilename += "_Expo.root";
  else if(typeb==1) outfilename += "_Linear.root";
  else if(typeb==2) outfilename += "_Pol2.root";
  TFile* outf=new TFile(outfilename,"recreate");
  outf->cd();
  dcuts[0]->Write();
  outf->Close();
  return kTRUE;

}

//_____________________________________________________________________________________________
Bool_t LoadDplusHistos(TObjArray* listFiles, TH3F** hPtMassMult, TH2F** hNtrZvtx, TH2F** hNtrZvtxCorr, const char *CutsType, Int_t Option)
{
Int_t nFiles=listFiles->GetEntries();
  TList **hlist=new TList*[nFiles];
  TList **hlistNorm=new TList*[nFiles];
  AliRDHFCutsDplustoKpipi** dcuts=new AliRDHFCutsDplustoKpipi*[nFiles];

  Int_t nReadFiles=0;
  for(Int_t iFile=0; iFile<nFiles; iFile++){
    TString fName=((TObjString*)listFiles->At(iFile))->GetString();    
    TFile *f=TFile::Open(fName.Data());
    if(!f){
      printf("ERROR: file %s does not exist\n",fName.Data());
      continue;
    }
    printf("Open File %s\n",f->GetName());
 TDirectory *dir = (TDirectory*)f->Get(Form("PWG3_D2H_DMult_Dplus%s",suffix.Data()));
    // TDirectory *dir = (TDirectory*)f->Get("PWG3_D2H_DMult");
    if(!dir){
      printf("ERROR: directory PWG3_D2H_DMult not found in %s\n",fName.Data());
      continue;
    }
    hlist[nReadFiles]=(TList*)dir->Get(Form("coutputDplus%s",suffix.Data()));
    TList *listcut = (TList*)dir->Get(Form("coutputCutsDplus%s",suffix.Data()));
    TList *listNorm = (TList*)dir->Get(Form("coutputNormDplus%s",suffix.Data()));
    dcuts[nReadFiles]=(AliRDHFCutsDplustoKpipi*)listcut->At(0);
     if(nReadFiles>0){
      Bool_t sameCuts=dcuts[nReadFiles]->CompareCuts(dcuts[0]);
      if(!sameCuts){
 
	printf("ERROR: Cut objects do not match\n");
	return kFALSE;
      }
     }



 nReadFiles++;
      }
  if(nReadFiles<nFiles){
    printf("WARNING: not all requested files have been found\n");
    if (nReadFiles==0) {
      printf("ERROR: Any file/dir found\n");
      return kFALSE;
    }
  }


 printf("Get the 3D histogram \n");
  const char *histoname="";
  if(optPartAntiPart==kParticleOnly) histoname= "hPtVsMassvsMultPart";
  else if(optPartAntiPart==kAntiParticleOnly) histoname="hPtVsMassvsMultAntiPart";
  else if(optPartAntiPart==kBoth) histoname="hPtVsMassvsMult";
  if(Option==kUnCorr) histoname="hPtVsMassvsMultUncorr";
  if(Option==kNoPid) histoname="hPtVsMassvsMultNoPid";

  TH3F * htemp;
  for(Int_t iFile=0; iFile<nReadFiles; iFile++){
    printf(" Looking for histo histMass %s for file %d\n",histoname,iFile);
    htemp=(TH3F*)hlist[iFile]->FindObject(Form("%s",histoname));
    //    cout << " htemp :"<<htemp<<endl;
    if(!hPtMassMult[0]){
      hPtMassMult[0]=new TH3F(*htemp);
    }else{
      hPtMassMult[0]->Add(htemp);
    }
    //  (TH2F*)hNtrZvtx[iFile] = (TH2F*)hlist[iFile]->FindObject("hNtrVsZvtx"); 
    //(TH2F*)hNtrZvtxCorr[iFile] = (TH2F*)hlist[iFile]->FindObject("hNtrVsZvtxCorr");
  }
  
  //  cout<<" hPtMassMult:"<<hPtMassMult[0]<<endl;

  TString partname="Both";
  if(optPartAntiPart==kParticleOnly) partname="Dplus";
  if(optPartAntiPart==kAntiParticleOnly) partname="Dminus";

  TString outfilename = Form("RawYield%s_%s",partname.Data(),CutsType);
  if(fixPeakSigma) outfilename += "_SigmaFixed";
  if(typeb==0) outfilename += "_Expo.root";
  else if(typeb==1) outfilename += "_Linear.root";
  else if(typeb==2) outfilename += "_Pol2.root";
  TFile* outf=new TFile(outfilename,"recreate");
  outf->cd();
  dcuts[0]->Write();
  outf->Close();
  return kTRUE;
  
}

//_____________________________________________________________________________________________
TH1F* RebinHisto(TH1F* hOrig, Int_t reb, Int_t firstUse){
  // Rebin histogram, from bin firstUse to lastUse
  // Use all bins if firstUse=-1
  
  Int_t nBinOrig=hOrig->GetNbinsX();
  Int_t firstBinOrig=1;
  Int_t lastBinOrig=nBinOrig;
  Int_t nBinOrigUsed=nBinOrig;
  Int_t nBinFinal=nBinOrig/reb;
  if(firstUse>=1){
    firstBinOrig=firstUse;
    nBinFinal=(nBinOrig-firstUse+1)/reb;
    nBinOrigUsed=nBinFinal*reb;
    lastBinOrig=firstBinOrig+nBinOrigUsed-1;
  }else{
    Int_t exc=nBinOrigUsed%reb;
    if(exc!=0){
      nBinOrigUsed-=exc;
      firstBinOrig+=exc/2;
      lastBinOrig=firstBinOrig+nBinOrigUsed-1;
    }
  }
  
  printf("Rebin from %d bins to %d bins -- Used bins=%d in range %d-%d\n",nBinOrig,nBinFinal,nBinOrigUsed,firstBinOrig,lastBinOrig);
  Float_t lowLim=hOrig->GetXaxis()->GetBinLowEdge(firstBinOrig);
  Float_t hiLim=hOrig->GetXaxis()->GetBinUpEdge(lastBinOrig);
  TH1F* hRebin=new TH1F(Form("%s-rebin",hOrig->GetName()),hOrig->GetTitle(),nBinFinal,lowLim,hiLim);
  Int_t lastSummed=firstBinOrig-1;
  for(Int_t iBin=1;iBin<=nBinFinal; iBin++){
    Float_t sum=0.;
    for(Int_t iOrigBin=0;iOrigBin<reb;iOrigBin++){
      sum+=hOrig->GetBinContent(lastSummed+1);
      lastSummed++;
    }
    hRebin->SetBinContent(iBin,sum);
  }
  return hRebin;
}

//_____________________________________________________________________________________________
Bool_t CheckNtrVsZvtx(TH2F** hNtrackVsVtxZ, TH2F** hNtrackVsVtxZCorr, Int_t nFiles)
{

  TCanvas *cNtrVsZvtx = new TCanvas("cNtrVsZvtx","Ntr Vs Zvtx");
  cNtrVsZvtx->Divide(2,2);
  TProfile **hProf = new TProfile*[nFiles];
  TProfile **hProfCorr = new TProfile*[nFiles];
  for(Int_t i=0; i<nFiles; i++){
    cNtrVsZvtx->cd(i+1);
    cNtrVsZvtx->SetLogz();
    //    hNtrackVsVtxZ[i]->Fit("pol4");
    hNtrackVsVtxZ[i]->Draw("colz");
    hProf[i] = (TProfile*)hNtrackVsVtxZ[i]->ProfileX(Form("%s_%d","hProf",i));
    hProf[i]->SetLineColor(kBlack);
    hProf[i]->Draw("same");
    cNtrVsZvtx->Update();
  }

  TCanvas *cNtrVsZvtxCorr = new TCanvas("cNtrVsZvtxCorr","Ntr Vs Zvtx Corr");
  cNtrVsZvtxCorr->Divide(2,2);
  for(Int_t i=0; i<nFiles; i++){
    cNtrVsZvtxCorr->cd(i+1);
    cNtrVsZvtxCorr->SetLogz();
    //    hNtrackVsVtxZCorr[i]->Fit("pol4");
    hNtrackVsVtxZCorr[i]->Draw("colz");
    hProfCorr[i] = (TProfile*)hNtrackVsVtxZCorr[i]->ProfileX(Form("%s_%d","hProfCorr",i));
    hProfCorr[i]->SetLineColor(kBlack);
    hProfCorr[i]->Draw("same");
    cNtrVsZvtx->Update();
  }

  TH1F *hNtrAxis = (TH1F*)hNtrackVsVtxZ[0]->ProjectionY("hNtrAxis");
  TH1F *hZvtx[nFiles];
  Int_t firstbin = hNtrAxis->FindBin( multlims[0] );
  Int_t lastbin = hNtrAxis->FindBin( multlims[nMultbins] ) -1;
  TCanvas *cZvtx = new TCanvas("cZvtx","Zvtx projections");
  cZvtx->Divide(2,2);
  for(Int_t i=0; i<nFiles; i++){
    cZvtx->cd(i+1);
    hZvtx[i] = (TH1F*)hNtrackVsVtxZ[i]->ProjectionX(Form("hZvtx_%d",i),firstbin,lastbin);
    hZvtx[i]->Draw();
  }
  TH1F *hZvtxCorr[nFiles];
  TCanvas *cZvtxCorr = new TCanvas("cZvtxCorr","Zvtx projections Corr");
  cZvtxCorr->Divide(2,2);
  for(Int_t i=0; i<nFiles; i++){
    cZvtxCorr->cd(i+1);
    hZvtxCorr[i] = (TH1F*)hNtrackVsVtxZCorr[i]->ProjectionX(Form("hZvtxCorr_%d",i),firstbin,lastbin);
    hZvtxCorr[i]->Draw();
  }

  return kTRUE;
}
