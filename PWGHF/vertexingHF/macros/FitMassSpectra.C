#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TInterpreter.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TDatabasePDG.h>
#include <TGraph.h>

#include "AliAODRecoDecayHF.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliHFMassFitter.h"
#include "AliNormalizationCounter.h"
#endif


// MACRO to perform fits to D meson invariant mass spectra
// and store raw yields and cut object into a root output file
// Origin: F. Prino (prino@to.infn.it)
// D0: C. Bianchin (cbianchi@pd.infn.it)



//
enum {kD0toKpi, kDplusKpipi, kDStarD0pi, kDsKKpi};
enum {kBoth, kParticleOnly, kAntiParticleOnly};
enum {kExpo=0, kLinear, kPol2, kThrExpo=5};
enum {kGaus=0, kDoubleGaus};


// Common variables: to be configured by the user
const Int_t nPtBins=6;
Double_t ptlims[nPtBins+1]={2.,3.,4.,5.,6.,8.,12.};
Int_t rebin[nPtBins]={2,4,4,4,4,4};
Int_t firstUsedBin[nPtBins]={-1,-1,-1,-1,-1,-1}; // -1 uses all bins, >=1 to set the lower bin to be accepted from original histo

TString suffix="Loose_SecondSet1236_ForCF08";


//const Int_t nPtBins=7;//6;
//Double_t ptlims[nPtBins+1]={1.,2.,3.,4.,5.,6.,8.,12.};
//Int_t rebin[nPtBins+1]={8,6,10,10,10,10,10,10}; //for looser cuts
//Int_t rebin[nPtBins+1]={10,10,10,10,10,10,10,10}; //for Chiara's cuts
Int_t typeb=kExpo;
Int_t types=kGaus;
Int_t optPartAntiPart=kBoth;
Int_t factor4refl=0;
Double_t minMassForFit=1.7;
Double_t maxMassForFit=2.1;
Float_t massRangeForCounting=0.05; // GeV
TH2F* hPtMass=0x0;
Double_t nEventsForNorm=0.;

//for D0only
Bool_t cutsappliedondistr=kTRUE;
const Int_t nsamples=2;//3;
//Int_t nevents[nsamples]={8.5635859e+07 /*LHC10b+cpart*/,8.9700624e+07/*LHC10dpart*/};
//Int_t nevents[nsamples]={9.0374946e+07 /*LHC10b+c*/,9.7593785e+07/*LHC10d*/};
//Int_t nevents[nsamples]={1.1777812e+08 /*LHC10dnewTPCpid*/,9.7593785e+07/*LHC10d*/};
//Int_t nevents[nsamples]={1.1777812e+08 /*LHC10dnewTPCpid*/,9.7593785e+07/*LHC10d*/,9.0374946e+07 /*LHC10b+c*/};
Int_t nevents[nsamples]={1.18860695e+08 /*LHC10dnewTPCpid*/,9.0374946e+07 /*LHC10b+c*/};
//Int_t nevents[nsamples]={1. /*LHC10dnewTPCpid*/,1 /*LHC10dnewTPCpidrescale*/};
// Functions

Bool_t LoadDplusHistos(TObjArray* listFiles, TH1F** hMass);
Bool_t LoadDsHistos(TObjArray* listFiles, TH1F** hMass);
Bool_t LoadD0toKpiHistos(TObjArray* listFiles, TH1F** hMass);
Bool_t LoadDstarD0piHistos(TObjArray* listFiles, TH1F** hMass);




void FitMassSpectra(Int_t analysisType=kDplusKpipi,
		    TString fileNameb="",
		    TString fileNamec="",
		    TString fileNamed="",
		    TString fileNamee=""
	       ){
  //

  gInterpreter->ExecuteMacro("$ALICE_ROOT/PWGHF/vertexingHF/macros/LoadLibraries.C");
  gStyle->SetOptTitle(0);

  TObjArray* listFiles=new TObjArray();
  if(fileNameb!="") listFiles->AddLast(new TObjString(fileNameb.Data()));
  if(fileNamec!="") listFiles->AddLast(new TObjString(fileNamec.Data()));
  if(fileNamed!="") listFiles->AddLast(new TObjString(fileNamed.Data()));
  if(fileNamee!="") listFiles->AddLast(new TObjString(fileNamee.Data()));
  if(listFiles->GetEntries()==0){
    printf("Missing file names in input\n");
    return;
  }

  TH1F** hmass=new TH1F*[nPtBins];
  for(Int_t i=0;i<nPtBins;i++) hmass[i]=0x0;

  Float_t massD, massD0_fDstar;
  Bool_t retCode;
  if(analysisType==kD0toKpi){
    retCode=LoadD0toKpiHistos(listFiles,hmass);
    massD=TDatabasePDG::Instance()->GetParticle(421)->Mass();
  }
  else if(analysisType==kDplusKpipi){
    retCode=LoadDplusHistos(listFiles,hmass);
    massD=TDatabasePDG::Instance()->GetParticle(411)->Mass();
  }
 else if(analysisType==kDStarD0pi){
    retCode=LoadDstarD0piHistos(listFiles,hmass);
    massD=TDatabasePDG::Instance()->GetParticle(413)->Mass();
    massD0_fDstar=TDatabasePDG::Instance()->GetParticle(421)->Mass();
    massD =massD-massD0_fDstar;
  }
  else if(analysisType==kDsKKpi){
    retCode=LoadDsHistos(listFiles,hmass);
    massD=TDatabasePDG::Instance()->GetParticle(431)->Mass();
  }
  else{
    printf("Wrong analysis type parameter\n");
    return;
  }
  if(!retCode){
    printf("ERROR in reading input files\n");
    return;
  } 
  


  TH1D* hCntSig1=new TH1D("hCntSig1","hCntSig1",nPtBins,ptlims);
  TH1D* hCntSig2=new TH1D("hCntSig2","hCntSig2",nPtBins,ptlims);
  TH1D* hNDiffCntSig1=new TH1D("hNDiffCntSig1","hNDiffCntSig1",nPtBins,ptlims);
  TH1D* hNDiffCntSig2=new TH1D("hNDiffCntSig2","hNDiffCntSig2",nPtBins,ptlims);
  TH1D* hSignal=new TH1D("hSignal","hSignal",nPtBins,ptlims);
  TH1D* hRelErrSig=new TH1D("hRelErrSig","hRelErrSig",nPtBins,ptlims);
  TH1D* hInvSignif=new TH1D("hInvSignif","hInvSignif",nPtBins,ptlims);
  TH1D* hBackground=new TH1D("hBackground","hBackground",nPtBins,ptlims);
  TH1D* hBackgroundNormSigma=new TH1D("hBackgroundNormSigma","hBackgroundNormSigma",nPtBins,ptlims);
  TH1D* hSignificance=new TH1D("hSignificance","hSignificance",nPtBins,ptlims);
  TH1D* hMass=new TH1D("hMass","hMass",nPtBins,ptlims);
  TH1D* hSigma=new TH1D("hSigma","hSigma",nPtBins,ptlims);


  Int_t nMassBins=hmass[0]->GetNbinsX();
  Double_t hmin=TMath::Max(minMassForFit,hmass[0]->GetBinLowEdge(2));
  Double_t hmax=TMath::Min(maxMassForFit,hmass[0]->GetBinLowEdge(nMassBins-2)+hmass[0]->GetBinWidth(nMassBins-2));
  Float_t minBinSum=hmass[0]->FindBin(massD-massRangeForCounting);
  Float_t maxBinSum=hmass[0]->FindBin(massD+massRangeForCounting);
  Int_t iPad=1;

  TF1* funBckStore1=0x0;
  TF1* funBckStore2=0x0;
  TF1* funBckStore3=0x0;

  AliHFMassFitter** fitter=new AliHFMassFitter*[nPtBins];
  Double_t arrchisquare[nPtBins];
  TCanvas* c1= new TCanvas("c1","MassSpectra",1000,800);
  Int_t nx=3, ny=2;
  if(nPtBins>6){
    nx=4;
    ny=3;
  }
  if(nPtBins>12){
    nx=5;
    ny=4;
  }
  c1->Divide(nx,ny);

  Double_t sig,errsig,s,errs,b,errb;
  for(Int_t iBin=0; iBin<nPtBins; iBin++){
    hmass[iBin]->SetName(Form("hInvMass_PtBin%d",iBin));
    hmass[iBin]->SetTitle(Form("%.1f<pt<%.1f",ptlims[iBin],ptlims[iBin+1]));    
    c1->cd(iPad++);
    Int_t origNbins=hmass[iBin]->GetNbinsX();
    TH1F* hRebinned=(TH1F*)AliVertexingHFUtils::RebinHisto(hmass[iBin],rebin[iBin],firstUsedBin[iBin]);
    hRebinned->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
    hRebinned->GetYaxis()->SetTitle(Form("Entries/(%.0f MeV/c^{2})",(hRebinned->GetBinWidth(1)*1000)));
    hRebinned->GetYaxis()->SetTitleOffset(1.1);
    hmin=TMath::Max(minMassForFit,hRebinned->GetBinLowEdge(2));
    hmax=TMath::Min(maxMassForFit,hRebinned->GetBinLowEdge(hRebinned->GetNbinsX()));
    fitter[iBin]=new AliHFMassFitter(hRebinned,hmin, hmax,1,typeb,types);
    rebin[iBin]=origNbins/fitter[iBin]->GetBinN();
    fitter[iBin]->SetReflectionSigmaFactor(factor4refl);
    fitter[iBin]->SetInitialGaussianMean(massD);
    if(analysisType==kDStarD0pi) fitter[iBin]->SetInitialGaussianSigma(0.00065);
    Bool_t out=fitter[iBin]->MassFitter(0);
    if(!out) {
      fitter[iBin]->GetHistoClone()->Draw();
      continue;
    }
    Double_t mass=fitter[iBin]->GetMean();
    Double_t sigma=fitter[iBin]->GetSigma();
    arrchisquare[iBin]=fitter[iBin]->GetReducedChiSquare();
    TF1* fB1=fitter[iBin]->GetBackgroundFullRangeFunc();
    TF1* fB2=fitter[iBin]->GetBackgroundRecalcFunc();
    TF1* fM=fitter[iBin]->GetMassFunc();
    if(iBin==0 && fB1) funBckStore1=new TF1(*fB1);
    if(iBin==0 && fB2) funBckStore2=new TF1(*fB2);
    if(iBin==0 && fM) funBckStore3=new TF1(*fM);

    fitter[iBin]->DrawHere(gPad);    
    fitter[iBin]->Signal(3,s,errs);
    fitter[iBin]->Background(3,b,errb);
    fitter[iBin]->Significance(3,sig,errsig);
    Double_t ry=fitter[iBin]->GetRawYield();
    Double_t ery=fitter[iBin]->GetRawYieldError();
    Float_t cntSig1=0.;
    Float_t cntSig2=0.;
    Float_t cntErr=0.;
    for(Int_t iMB=minBinSum; iMB<=maxBinSum; iMB++){
      Float_t bkg1=fB1 ? fB1->Eval(hmass[iBin]->GetBinCenter(iMB))/rebin[iBin] : 0;
      Float_t bkg2=fB2 ? fB2->Eval(hmass[iBin]->GetBinCenter(iMB))/rebin[iBin] : 0;
      cntSig1+=(hmass[iBin]->GetBinContent(iMB)-bkg1);
      cntSig2+=(hmass[iBin]->GetBinContent(iMB)-bkg2);
      cntErr+=(hmass[iBin]->GetBinContent(iMB));
    }
    hCntSig1->SetBinContent(iBin+1,cntSig1);
    hCntSig1->SetBinError(iBin+1,TMath::Sqrt(cntErr));
    hNDiffCntSig1->SetBinContent(iBin+1,(s-cntSig1)/s);
    hNDiffCntSig1->SetBinError(iBin+1,TMath::Sqrt(cntErr)/s);
    hCntSig2->SetBinContent(iBin+1,cntSig2);
    hNDiffCntSig2->SetBinContent(iBin+1,(s-cntSig2)/s);
    hNDiffCntSig2->SetBinError(iBin+1,TMath::Sqrt(cntErr)/s);
    hCntSig2->SetBinError(iBin+1,TMath::Sqrt(cntErr));
    hSignal->SetBinContent(iBin+1,ry);
    hSignal->SetBinError(iBin+1,ery);
    hRelErrSig->SetBinContent(iBin+1,errs/s);
    hInvSignif->SetBinContent(iBin+1,1/sig);
    hInvSignif->SetBinError(iBin+1,errsig/(sig*sig));
    hBackground->SetBinContent(iBin+1,b); //consider sigma
    hBackground->SetBinError(iBin+1,errb);
    hBackgroundNormSigma->SetBinContent(iBin+1,b/(3*fitter[iBin]->GetSigma())*(3*0.012)); //consider sigma
    hBackgroundNormSigma->SetBinError(iBin+1,errb);
    hSignificance->SetBinContent(iBin+1,sig);
    hSignificance->SetBinError(iBin+1,errsig);
    hMass->SetBinContent(iBin+1,mass);
    hMass->SetBinError(iBin+1,0.0001);
    hSigma->SetBinContent(iBin+1,sigma);
    hSigma->SetBinError(iBin+1,fitter[iBin]->GetSigmaUncertainty());
    
  }

  /*
  c1->cd(1); // is some cases the fitting function of 1st bin get lost
  funBckStore1->Draw("same");
  funBckStore2->Draw("same");
  funBckStore3->Draw("same");
  */

  TCanvas *cpar=new TCanvas("cpar","Fit params",1200,600);
  cpar->Divide(2,1);
  cpar->cd(1);
  hMass->SetMarkerStyle(20);
  hMass->Draw("PE");
  hMass->GetXaxis()->SetTitle("Pt (GeV/c)");
  hMass->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
  cpar->cd(2);
  hSigma->SetMarkerStyle(20);
  hSigma->Draw("PE");
  hSigma->GetXaxis()->SetTitle("Pt (GeV/c)");
  hSigma->GetXaxis()->SetTitle("Sigma (GeV/c^{2})");

  TCanvas* csig=new TCanvas("csig","Results",1200,600);
  csig->Divide(3,1);
  csig->cd(1);
  hSignal->SetMarkerStyle(20);
  hSignal->SetMarkerColor(4);
  hSignal->SetLineColor(4);
  hSignal->GetXaxis()->SetTitle("Pt (GeV/c)");
  hSignal->GetYaxis()->SetTitle("Signal");
  hSignal->Draw("P");
  hCntSig1->SetMarkerStyle(26);
  hCntSig1->SetMarkerColor(2);
  hCntSig1->SetLineColor(2);
  hCntSig1->Draw("PSAME");
  hCntSig2->SetMarkerStyle(29);
  hCntSig2->SetMarkerColor(kGray+1);
  hCntSig2->SetLineColor(kGray+1);
  hCntSig2->Draw("PSAME");
  TLegend* leg=new TLegend(0.4,0.7,0.89,0.89);
  leg->SetFillColor(0);
  TLegendEntry* ent=leg->AddEntry(hSignal,"From Fit","PL");
  ent->SetTextColor(hSignal->GetMarkerColor());
  ent=leg->AddEntry(hCntSig1,"From Counting1","PL");
  ent->SetTextColor(hCntSig1->GetMarkerColor());
  ent=leg->AddEntry(hCntSig2,"From Counting2","PL");
  ent->SetTextColor(hCntSig2->GetMarkerColor());
  leg->Draw();
  csig->cd(2);
  hBackground->SetMarkerStyle(20);
  hBackground->Draw("P");
  hBackground->GetXaxis()->SetTitle("Pt (GeV/c)");
  hBackground->GetYaxis()->SetTitle("Background");
  csig->cd(3);
  hSignificance->SetMarkerStyle(20);
  hSignificance->Draw("P");
  hSignificance->GetXaxis()->SetTitle("Pt (GeV/c)");
  hSignificance->GetYaxis()->SetTitle("Significance");

  TCanvas* cDiffS=new TCanvas("cDiffS","Difference",1200,600);
  cDiffS->Divide(2,1);
  cDiffS->cd(1);
  hRelErrSig->SetMarkerStyle(20); //fullcircle
  hRelErrSig->SetTitleOffset(1.2);  
  hRelErrSig->SetTitle("Rel Error from Fit;p_{t} (GeV/c);Signal Relative Error");
  hRelErrSig->Draw("P");
  hInvSignif->SetMarkerStyle(21); //fullsquare
  hInvSignif->SetMarkerColor(kMagenta+1);
  hInvSignif->SetLineColor(kMagenta+1);
  hInvSignif->Draw("PSAMES");
  TLegend* leg2=new TLegend(0.4,0.7,0.89,0.89);
  leg2->SetFillColor(0);
  TLegendEntry* ent2=leg2->AddEntry(hRelErrSig,"From Fit","P");
  ent2->SetTextColor(hRelErrSig->GetMarkerColor());
  ent2=leg2->AddEntry(hInvSignif,"1/Significance","PL");
  ent2->SetTextColor(hInvSignif->GetMarkerColor());
  leg2->Draw();

  cDiffS->cd(2);
  hNDiffCntSig1->SetMarkerStyle(26);
  hNDiffCntSig1->SetMarkerColor(2);
  hNDiffCntSig1->SetLineColor(2);
  hNDiffCntSig1->SetTitle("Cmp Fit-Count;p_{t} (GeV/c);(S_{fit}-S_{count})/S_{fit}");
  hNDiffCntSig1->Draw("P");
  hNDiffCntSig2->SetMarkerStyle(29);
  hNDiffCntSig2->SetMarkerColor(kGray+1);
  hNDiffCntSig2->SetLineColor(kGray+1);
  hNDiffCntSig2->Draw("PSAME");
  TLegend* leg1=new TLegend(0.4,0.7,0.89,0.89);
  leg1->SetFillColor(0);
  TLegendEntry* ent1=leg1->AddEntry(hNDiffCntSig1,"From Counting1","PL");
  ent1->SetTextColor(hNDiffCntSig1->GetMarkerColor());
  ent1=leg1->AddEntry(hNDiffCntSig2,"From Counting2","PL");
  ent1->SetTextColor(hNDiffCntSig2->GetMarkerColor());
  leg1->Draw();

  TGraph* grReducedChiSquare=new TGraph(nPtBins,ptlims,arrchisquare);
  grReducedChiSquare->SetName("grReducedChiSquare");
  grReducedChiSquare->SetTitle("Reduced Chi2;p_{t} (GeV/c);#tilde{#chi}^{2}");
  TCanvas *cChi2=new TCanvas("cChi2","reduced chi square",600,600);
  cChi2->cd();
  grReducedChiSquare->SetMarkerStyle(21);
  grReducedChiSquare->Draw("AP");

  TCanvas* cbkgNormSigma=new TCanvas("cbkgNormSigma","Background normalized to sigma",400,600);
  cbkgNormSigma->cd();
  hBackgroundNormSigma->SetMarkerStyle(20);
  hBackgroundNormSigma->Draw("P");
  hBackgroundNormSigma->GetXaxis()->SetTitle("Pt (GeV/c)");
  hBackgroundNormSigma->GetYaxis()->SetTitle("Background #times 3 #times 0.012/ (3 #times #sigma)");
  hBackgroundNormSigma->Draw();


  TString partname="Both";
  if(optPartAntiPart==kParticleOnly) {
    if(analysisType==kD0toKpi) partname="D0";
    if(analysisType==kDplusKpipi) partname="Dplus";
    if(analysisType==kDsKKpi) partname="Dsplus";
  }
  if(optPartAntiPart==kAntiParticleOnly) {
    if(analysisType==kD0toKpi) partname="D0bar";
    if(analysisType==kDplusKpipi) partname="Dminus";
    if(analysisType==kDsKKpi) partname="Dsminus";
  }

  printf("Events for norm = %f\n",nEventsForNorm);
  TH1F* hNEvents=new TH1F("hNEvents","",1,0.,1.);
  hNEvents->SetBinContent(1,nEventsForNorm);

  TFile* outf=new TFile(Form("RawYield%s.root",partname.Data()),"update");
  outf->cd();
  hNEvents->Write();
  hMass->Write();
  hSigma->Write();
  hCntSig1->Write();
  hCntSig2->Write();
  hNDiffCntSig1->Write();
  hNDiffCntSig2->Write();
  hSignal->Write();
  hRelErrSig->Write();
  hInvSignif->Write();
  hBackground->Write();
  hBackgroundNormSigma->Write();
  hSignificance->Write();
  grReducedChiSquare->Write();
  hPtMass->Write();
  outf->Close();
}


Bool_t LoadDplusHistos(TObjArray* listFiles, TH1F** hMass){

  Int_t nFiles=listFiles->GetEntries();
  TList **hlist=new TList*[nFiles];
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
    TDirectory *dir = (TDirectory*)f->Get("PWG3_D2H_InvMassDplus");
    if(!dir){
      printf("ERROR: directory PWG3_D2H_InvMassDplus not found in %s\n",fName.Data());
      continue;
    }
    hlist[nReadFiles]=(TList*)dir->Get(Form("coutputDplus%s",suffix.Data()));
    TList *listcut = (TList*)dir->Get(Form("coutputDplusCuts%s",suffix.Data()));
    dcuts[nReadFiles]=(AliRDHFCutsDplustoKpipi*)listcut->At(0);
    if(nReadFiles>0){
      Bool_t sameCuts=dcuts[nReadFiles]->CompareCuts(dcuts[0]);
      if(!sameCuts){
	printf("ERROR: Cut objects do not match\n");
	return kFALSE;
      }
    }
    AliNormalizationCounter* c=(AliNormalizationCounter*)dir->Get(Form("coutputDplusNorm%s",suffix.Data()));
    printf("Events for normalization = %f\n",c->GetNEventsForNorm());
    nEventsForNorm+=c->GetNEventsForNorm();
    nReadFiles++;
  }
  if(nReadFiles<nFiles){
    printf("WARNING: not all requested files have been found\n");
  }

  Int_t nPtBinsCuts=dcuts[0]->GetNPtBins();
  printf("Number of pt bins for cut object = %d\n",nPtBins);
  Float_t *ptlimsCuts=dcuts[0]->GetPtBinLimits();
  ptlimsCuts[nPtBinsCuts]=ptlimsCuts[nPtBinsCuts-1]+4.;

  Int_t iFinBin=0;
  for(Int_t i=0;i<nPtBinsCuts;i++){
    if(ptlimsCuts[i]>=ptlims[iFinBin+1]) iFinBin+=1; 
    if(iFinBin>nPtBins) break;
    if(ptlimsCuts[i]>=ptlims[iFinBin] && 
       ptlimsCuts[i+1]<=ptlims[iFinBin+1]){
      for(Int_t iFile=0; iFile<nReadFiles; iFile++){
	TString histoName;
	if(optPartAntiPart==kBoth) histoName.Form("hMassPt%dTC",i);
	else if(optPartAntiPart==kParticleOnly) histoName.Form("hMassPt%dTCPlus",i);
	else if(optPartAntiPart==kAntiParticleOnly) histoName.Form("hMassPt%dTCMinus",i);
	TH1F* htemp=(TH1F*)hlist[iFile]->FindObject(histoName.Data());
	if(!htemp){
	  printf("ERROR: Histogram %s not found\n",histoName.Data());
	  return kFALSE;
	}
	if(!hMass[iFinBin]){
	  hMass[iFinBin]=new TH1F(*htemp);
	}else{
	  hMass[iFinBin]->Add(htemp);
	}
      }
    }
  }
  TString partname="Both";
  if(optPartAntiPart==kParticleOnly) partname="Dplus";
  if(optPartAntiPart==kAntiParticleOnly) partname="Dminus";

  for(Int_t iFile=0; iFile<nReadFiles; iFile++){
    TH2F* htemp2=(TH2F*)hlist[iFile]->FindObject("hPtVsMassTC");
    if(iFile==0){
      hPtMass=new TH2F(*htemp2);
    }else{
      hPtMass->Add(htemp2);
    }
  }

  TFile* outf=new TFile(Form("RawYield%s.root",partname.Data()),"recreate");
  outf->cd();
  dcuts[0]->Write();
  outf->Close();

  return kTRUE;

}


Bool_t LoadDsHistos(TObjArray* listFiles, TH1F** hMass){

  Int_t nFiles=listFiles->GetEntries();
  TList **hlist=new TList*[nFiles];
  AliRDHFCutsDstoKKpi** dcuts=new AliRDHFCutsDstoKKpi*[nFiles];

  Int_t nReadFiles=0;
  for(Int_t iFile=0; iFile<nFiles; iFile++){
    TString fName=((TObjString*)listFiles->At(iFile))->GetString();    
    TFile *f=TFile::Open(fName.Data());
    if(!f){
      printf("ERROR: file %s does not exist\n",fName.Data());
      continue;
    }
    printf("Open File %s\n",f->GetName()); 
    TDirectory *dir = (TDirectory*)f->Get("PWG3_D2H_InvMassDs");
    if(!dir){
      printf("ERROR: directory PWG3_D2H_InvMassDs not found in %s\n",fName.Data());
      continue;
    }
    hlist[nReadFiles]=(TList*)dir->Get("coutputDs");
    TList *listcut = (TList*)dir->Get("coutputDsCuts");
    dcuts[nReadFiles]=(AliRDHFCutsDstoKKpi*)listcut->At(0);
    cout<< dcuts[nReadFiles]<<endl;
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
  }

  Int_t nPtBinsCuts=dcuts[0]->GetNPtBins();
  printf("Number of pt bins for cut object = %d\n",nPtBins);
  Float_t *ptlimsCuts=dcuts[0]->GetPtBinLimits();
  ptlimsCuts[nPtBinsCuts]=ptlimsCuts[nPtBinsCuts-1]+4.;

  Int_t iFinBin=0;
  for(Int_t i=0;i<nPtBinsCuts;i++){
    if(ptlimsCuts[i]>=ptlims[iFinBin+1]) iFinBin+=1; 
    if(iFinBin>nPtBins) break;
    if(ptlimsCuts[i]>=ptlims[iFinBin] && 
       ptlimsCuts[i+1]<=ptlims[iFinBin+1]){
      for(Int_t iFile=0; iFile<nReadFiles; iFile++){
	TString histoName;
	if(optPartAntiPart==kBoth) histoName.Form("hMassAllPt%dphi",i);
	else if(optPartAntiPart==kParticleOnly){
	  printf("Particle/Antiparticle not yet enabled for Ds");
	  histoName.Form("hMassAllPt%dphi",i);
	}
	else if(optPartAntiPart==kAntiParticleOnly){
	  printf("Particle/Antiparticle not yet enabled for Ds");
	  histoName.Form("hMassAllPt%dphi",i);
	}
	TH1F* htemp=(TH1F*)hlist[iFile]->FindObject(histoName.Data());
	if(!htemp){
	  printf("ERROR: Histogram %s not found\n",histoName.Data());
	  return kFALSE;
	}
	if(!hMass[iFinBin]){
	  hMass[iFinBin]=new TH1F(*htemp);
	}else{
	  hMass[iFinBin]->Add(htemp);
	}
      }
    }
  }
  TString partname="Both";
  if(optPartAntiPart==kParticleOnly) partname="Both";
  if(optPartAntiPart==kAntiParticleOnly) partname="Both";

  for(Int_t iFile=0; iFile<nReadFiles; iFile++){
    TH2F* htemp2=(TH2F*)hlist[iFile]->FindObject("hPtVsMassPhi");
    if(iFile==0){
      hPtMass=new TH2F(*htemp2);
    }else{
      hPtMass->Add(htemp2);
    }
  }

  TFile* outf=new TFile(Form("RawYield%s.root",partname.Data()),"recreate");
  outf->cd();
  dcuts[0]->Write();
  outf->Close();

  return kTRUE;

}


Bool_t LoadD0toKpiHistos(TObjArray* listFiles, TH1F** hMass){
  //
  Int_t nFiles=listFiles->GetEntries();
  TList **hlist=new TList*[nFiles];
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

    TString dirname="PWG3_D2H_D0InvMass";
    if(optPartAntiPart==kParticleOnly) dirname+="D0";
    else if(optPartAntiPart==kAntiParticleOnly) dirname+="D0bar";
    if(cutsappliedondistr) dirname+="C";
    TDirectory *dir = (TDirectory*)f->Get(dirname);
    if(!dir){
      printf("ERROR: directory %s not found in %s\n",dirname.Data(),fName.Data());
      continue;
    }
    TString listmassname="coutputmassD0Mass";
    if(optPartAntiPart==kParticleOnly) listmassname+="D0";
    else if(optPartAntiPart==kAntiParticleOnly) listmassname+="D0bar";
    if(cutsappliedondistr) listmassname+="C";

    hlist[nReadFiles]=(TList*)dir->Get(listmassname);

    TString cutsobjname="cutsD0";
    if(optPartAntiPart==kParticleOnly) cutsobjname+="D0";
    else if(optPartAntiPart==kAntiParticleOnly) cutsobjname+="D0bar";
    if(cutsappliedondistr) cutsobjname+="C";

    dcuts[nReadFiles]=(AliRDHFCutsD0toKpi*)dir->Get(cutsobjname);
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

  Int_t nPtBinsCuts=dcuts[0]->GetNPtBins();
  printf("Number of pt bins for cut object = %d\n",nPtBins);
  Float_t *ptlimsCuts=dcuts[0]->GetPtBinLimits();
  ptlimsCuts[nPtBinsCuts]=ptlimsCuts[nPtBinsCuts-1]+4.;

  Int_t iFinBin=0;
  for(Int_t i=0;i<nPtBinsCuts;i++){
    if(ptlimsCuts[i]>=ptlims[iFinBin+1]) iFinBin+=1; 
    if(iFinBin>nPtBins) break;
    if(ptlimsCuts[i]>=ptlims[iFinBin] && 
       ptlimsCuts[i+1]<=ptlims[iFinBin+1]){
      for(Int_t iFile=0; iFile<nReadFiles; iFile++){
	TH1F* htemp=(TH1F*)hlist[iFile]->FindObject(Form("histMass_%d",i));
	if(!hMass[iFinBin]){
	  hMass[iFinBin]=new TH1F(*htemp);
	}else{
	  hMass[iFinBin]->Add(htemp);
	}
      }
    }
  }
  TString partname="Both";
  if(optPartAntiPart==kParticleOnly) partname="D0";
  if(optPartAntiPart==kAntiParticleOnly) partname="D0bar";

  TFile* outf=new TFile(Form("RawYield%s.root",partname.Data()),"recreate");
  outf->cd();
  dcuts[0]->Write();
  outf->Close();
  return kTRUE;
}

Bool_t LoadDstarD0piHistos(TObjArray* listFiles, TH1F** hMass){
  //
  Int_t nFiles=listFiles->GetEntries();
  TList **hlist=new TList*[nFiles];
  AliRDHFCutsDStartoKpipi** dcuts=new AliRDHFCutsDStartoKpipi*[nFiles];
  Int_t nReadFiles=0;
  for(Int_t iFile=0; iFile<nFiles; iFile++){
    TString fName=((TObjString*)listFiles->At(iFile))->GetString();    
    TFile *f=TFile::Open(fName.Data());
    if(!f){
      printf("ERROR: file %s does not exist\n",fName.Data());
      continue;
    }
    printf("Open File %s\n",f->GetName());
    TString dirname="PWG3_D2H_DStarSpectra";
    TDirectory *dir = (TDirectory*)f->Get(dirname);
    if(!dir){
      printf("ERROR: directory %s not found in %s\n",dirname.Data(),fName.Data());
      continue;
    }
    TString listmassname="DStarPID00";

    hlist[nReadFiles]=(TList*)dir->Get(listmassname);
    TString cutsobjname="DStartoKpipiCuts";
    dcuts[nReadFiles]=(AliRDHFCutsDStartoKpipi*)dir->Get(cutsobjname);
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
  Int_t nPtBinsCuts=dcuts[0]->GetNPtBins();
  printf("Number of pt bins for cut object = %d\n",nPtBins);
  Float_t *ptlimsCuts=dcuts[0]->GetPtBinLimits();
  ptlimsCuts[nPtBinsCuts]=ptlimsCuts[nPtBinsCuts-1]+4.;




  Int_t iFinBin=0;
  for(Int_t i=0;i<nPtBinsCuts;i++){
    if(ptlimsCuts[i]>=ptlims[iFinBin+1]) iFinBin+=1; 
    if(iFinBin>nPtBins) break;
    if(ptlimsCuts[i]>=ptlims[iFinBin] && 
       ptlimsCuts[i+1]<=ptlims[iFinBin+1]){
      for(Int_t iFile=0; iFile<nReadFiles; iFile++){
	TH1F* htemp=(TH1F*)hlist[iFile]->FindObject(Form("histDeltaMass_%d",i));
	if(!hMass[iFinBin]){
	  hMass[iFinBin]=new TH1F(*htemp);
	}else{
	  hMass[iFinBin]->Add(htemp);
	}
      }
    }
  }
  TString partname="Both";
  if(optPartAntiPart==kParticleOnly) partname="DStar";
  if(optPartAntiPart==kAntiParticleOnly) partname="DStarbar";

  TFile* outf=new TFile(Form("RawYield%s.root",partname.Data()),"recreate");
  outf->cd();
  dcuts[0]->Write();
  outf->Close();
  return kTRUE;
}

void CompareFitTypes(TString* paths, TString* legtext,Int_t ncmp=3,TString* filenameYield=0x0){
  //read ncmp RawYield.roots and draw them together
  //set the global variable nevents before running
  //arguments:
  // - paths= vector of ncmp dimension with the paths of the file RawYield.root
  // - legtext= vector of ncmp dimension with the label for the legend
  // - ncmp= number of files to compare (default is 3)
  // - filenameYield= optional ncmp-dimensional array with the difference between the names of the files to be compared (useful if the 2 files are in the same directory but have different names)

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetFrameBorderMode(0);

  if(!filenameYield) filenameYield=new TString[ncmp];

  for(Int_t k=0;k<ncmp;k++){
    if(!filenameYield) filenameYield[k]="RawYield.root";
    filenameYield[k].Prepend(paths[k]);
  }
  
  TCanvas* cSig=new TCanvas("cSig","Results",1200,600);
  cSig->Divide(3,1);
  TCanvas* cBkgN=new TCanvas("cBkgN","Background normalized to sigma",400,600);
  TCanvas* cDiffS=new TCanvas("cDiffS","Difference",1200,600);
  cDiffS->Divide(2,1);
  TCanvas *cChi2=new TCanvas("cChi2","reduced chi square",600,600);

  TLegend* leg1=new TLegend(0.4,0.7,0.89,0.89);
  leg1->SetFillColor(0);
  TLegend* leg2=(TLegend*)leg1->Clone();
  TLegend* leg3=(TLegend*)leg1->Clone();
  TLegend* leg4=new TLegend(0.4,0.6,0.8,0.89);
  leg4->SetFillColor(0);

  for(Int_t iTypes=0;iTypes<ncmp;iTypes++){
    TFile* fin=new TFile(filenameYield[iTypes]);
    if(!fin){
      printf("WARNING: %s not found",filenameYield[iTypes].Data());
      continue;
    }

    TH1F* hSignal=(TH1F*)fin->Get("hSignal");
    TH1F* hBackground=(TH1F*)fin->Get("hBackground");
    TH1F* hBackgroundNormSigma=(TH1F*)fin->Get("hBackgroundNormSigma");
    TH1F* hSignificance=(TH1F*)fin->Get("hSignificance");
    hSignal->SetName(Form("%s%d",hSignal->GetName(),iTypes));
    hBackground->SetName(Form("%s%d",hBackground->GetName(),iTypes));
    hBackgroundNormSigma->SetName(Form("%s%d",hBackgroundNormSigma->GetName(),iTypes));
    hSignificance->SetName(Form("%s%d",hSignificance->GetName(),iTypes));

    hSignal->SetMarkerColor(iTypes+2);
    hSignal->SetLineColor(iTypes+2);
    hBackground->SetMarkerColor(iTypes+2);
    hBackground->SetLineColor(iTypes+2);
    hBackgroundNormSigma->SetMarkerColor(iTypes+2);
    hBackgroundNormSigma->SetLineColor(iTypes+2);
    hSignificance->SetMarkerColor(iTypes+2);
    hSignificance->SetLineColor(iTypes+2);

    TLegendEntry* ent4=leg4->AddEntry(hSignal,Form("%s",legtext[iTypes].Data()),"PL");
    ent4->SetTextColor(hSignal->GetMarkerColor());

    if(ncmp==nsamples){
      printf("Info: Normalizing signal, background and significance to the number of events\n");
      hSignal->Scale(1./nevents[iTypes]);
      hBackground->Scale(1./nevents[iTypes]);
      hBackgroundNormSigma->Scale(1./nevents[iTypes]);
      hSignificance->Scale(1./TMath::Sqrt(nevents[iTypes]));
    }

    if(iTypes==0){
      cSig->cd(1);
      hSignal->DrawClone("P");
      cSig->cd(2);
      hBackground->DrawClone("P");
      cSig->cd(3);
      hSignificance->DrawClone("P");
      cBkgN->cd();
      hBackgroundNormSigma->DrawClone("P");
    } else{
      cSig->cd(1);
      hSignal->DrawClone("Psames");
      cSig->cd(2);
      hBackground->DrawClone("Psames");
      cSig->cd(3);
      hSignificance->DrawClone("Psames");
      cBkgN->cd();
      hBackgroundNormSigma->DrawClone("Psames");
    }

    TH1F* hRelErrSig=(TH1F*)fin->Get("hRelErrSig");
    TH1F* hInvSignif=(TH1F*)fin->Get("hInvSignif");
    hRelErrSig->SetName(Form("%s%d",hRelErrSig->GetName(),iTypes));
    hInvSignif->SetName(Form("%s%d",hInvSignif->GetName(),iTypes));

    hRelErrSig->SetMarkerColor(iTypes+2);
    hRelErrSig->SetLineColor(iTypes+2);
    hInvSignif->SetMarkerColor(iTypes+2);
    hInvSignif->SetLineColor(iTypes+2);

    TLegendEntry* ent1=leg1->AddEntry(hRelErrSig,Form("From Fit (%s)",legtext[iTypes].Data()),"P");
    ent1->SetTextColor(hRelErrSig->GetMarkerColor());
    ent1=leg1->AddEntry(hInvSignif,Form("1/Significance (%s)",legtext[iTypes].Data()),"PL");
    ent1->SetTextColor(hInvSignif->GetMarkerColor());

    cDiffS->cd(1);
    if(iTypes==0){
      hRelErrSig->DrawClone("P");
      hInvSignif->DrawClone();
    } else{
      hRelErrSig->DrawClone("Psames");
      hInvSignif->DrawClone("sames");
    }

    TH1F* hNDiffCntSig1=(TH1F*)fin->Get("hNDiffCntSig1");
    TH1F* hNDiffCntSig2=(TH1F*)fin->Get("hNDiffCntSig2");
    hNDiffCntSig1->SetName(Form("%s%d",hNDiffCntSig1->GetName(),iTypes));
    hNDiffCntSig2->SetName(Form("%s%d",hNDiffCntSig2->GetName(),iTypes));

    hNDiffCntSig1->SetMarkerColor(iTypes+2);
    hNDiffCntSig1->SetLineColor(iTypes+2);
    hNDiffCntSig2->SetMarkerColor(iTypes+2);
    hNDiffCntSig2->SetLineColor(iTypes+2);
    TLegendEntry* ent2=leg2->AddEntry(hNDiffCntSig1,Form("From Counting1 (%s)",legtext[iTypes].Data()),"PL");
    ent2->SetTextColor(hNDiffCntSig1->GetMarkerColor());
    ent2=leg2->AddEntry(hNDiffCntSig2,Form("From Counting2 (%s)",legtext[iTypes].Data()),"PL");
    ent2->SetTextColor(hNDiffCntSig2->GetMarkerColor());

    cDiffS->cd(2);
    if(iTypes==0){
      hNDiffCntSig1->DrawClone();
      hNDiffCntSig2->DrawClone();
    }else{
     hNDiffCntSig1->DrawClone("sames");
     hNDiffCntSig2->DrawClone("sames");
    }

    TGraph* grReducedChiSquare=(TGraph*)fin->Get("grReducedChiSquare");
    grReducedChiSquare->SetName(Form("%s%d",grReducedChiSquare->GetName(),iTypes));

    grReducedChiSquare->SetMarkerColor(iTypes+2);
    grReducedChiSquare->SetLineColor(iTypes+2);
    TLegendEntry* ent3=leg3->AddEntry(grReducedChiSquare,Form("%s",legtext[iTypes].Data()),"PL");
    ent3->SetTextColor(grReducedChiSquare->GetMarkerColor());

    cChi2->cd();
    if(iTypes==0) grReducedChiSquare->DrawClone("AP");
    else grReducedChiSquare->DrawClone("P");
  }

  cSig->cd(1);
  leg4->Draw();

  cDiffS->cd(1);
  leg1->Draw();

  cDiffS->cd(2);
  leg2->Draw();

  cChi2->cd();
  leg3->Draw();

  TFile* fout=new TFile("ComparisonRawYield.root","RECREATE");
  fout->cd();
  cDiffS->Write();
  cChi2->Write();
  fout->Close();
}


