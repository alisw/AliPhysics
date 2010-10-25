#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TInterpreter.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TDatabasePDG.h>

#include "AliAODRecoDecayHF.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliHFMassFitter.h"
#endif


// MACRO to perform fits to D meson invariant mass spectra
// and store raw yields and cut object into a root output file
// Origin: F. Prino (prino@to.infn.it)
// D0: C. Bianchin (cbianchi@pd.infn.it)



//
enum {kD0toKpi, kDplusKpipi};
enum {kBoth, kParticleOnly, kAntiParticleOnly};
enum {kExpo=0, kLinear, kPol2};
enum {kGaus=0, kDoubleGaus};


// Common variables: to be configured by the user
const Int_t nPtBins=6;
Double_t ptlims[nPtBins+1]={2.,3.,4.,5.,6.,8.,12.};
Int_t rebin[nPtBins+1]={2,4,4,4,4,4,4};
// const Int_t nPtBins=7;
// Double_t ptlims[nPtBins+1]={1.,2.,3.,4.,5.,6.,8.,12.};
// Int_t rebin[nPtBins+1]={10,10,10,10,10,10,10,10};
Int_t typeb=kExpo;
Int_t types=kGaus;
Int_t optPartAntiPart=kBoth;
Int_t factor4refl=0;
Float_t massRangeForCounting=0.05; // GeV

//for D0only
Bool_t cutsappliedondistr=kFALSE;

// Functions

Bool_t LoadDplusHistos(TObjArray* listFiles, TH1F** hMass);
Bool_t LoadD0toKpiHistos(TObjArray* listFiles, TH1F** hMass);




void FitMassSpectra(Int_t analysisType=kDplusKpipi,
	       TString fileNameb="",
	       TString fileNamec="",
	       TString fileNamed=""
	       ){
  //

  gInterpreter->ExecuteMacro("$ALICE_ROOT/PWG3/vertexingHF/macros/LoadLibraries.C");
  gStyle->SetOptTitle(0);

  TObjArray* listFiles=new TObjArray();
  if(fileNameb!="") listFiles->AddLast(new TObjString(fileNameb.Data()));
  if(fileNamec!="") listFiles->AddLast(new TObjString(fileNamec.Data()));
  if(fileNamed!="") listFiles->AddLast(new TObjString(fileNamed.Data()));
  if(listFiles->GetEntries()==0){
    printf("Missing file names in input\n");
    return;
  }

  
  TH1F** hmass=new TH1F*[nPtBins];
  for(Int_t i=0;i<nPtBins;i++) hmass[i]=0x0;

  Float_t massD;
  Bool_t retCode;
  if(analysisType==kD0toKpi){
    retCode=LoadD0toKpiHistos(listFiles,hmass);
    massD=TDatabasePDG::Instance()->GetParticle(421)->Mass();
  }
  else if(analysisType==kDplusKpipi){
    retCode=LoadDplusHistos(listFiles,hmass);
    massD=TDatabasePDG::Instance()->GetParticle(411)->Mass();
  }
  else{
    printf("Wronganalysistype parameter\n");
    return;
  }
  if(!retCode){
    printf("ERROR in reading input files\n");
    return;
  } 
  


  TH1F* hCntSig1=new TH1F("hCntSig1","hCntSig1",nPtBins,ptlims);
  TH1F* hCntSig2=new TH1F("hCntSig2","hCntSig2",nPtBins,ptlims);
  TH1F* hSignal=new TH1F("hSignal","hSignal",nPtBins,ptlims);
  TH1F* hBackground=new TH1F("hBackground","hBackground",nPtBins,ptlims);
  TH1F* hSignificance=new TH1F("hSignificance","hSignificance",nPtBins,ptlims);


  Int_t nMassBins=hmass[1]->GetNbinsX();
  Double_t hmin=hmass[1]->GetBinLowEdge(3);
  Double_t hmax=hmass[1]->GetBinLowEdge(nMassBins-2)+hmass[1]->GetBinWidth(nMassBins-2);
  Float_t minBinSum=hmass[1]->FindBin(massD-massRangeForCounting);
  Float_t maxBinSum=hmass[1]->FindBin(massD+massRangeForCounting);
  Int_t iPad=1;

  TF1* funBckStore1=0x0;
  TF1* funBckStore2=0x0;
  TF1* funBckStore3=0x0;

  AliHFMassFitter** fitter=new AliHFMassFitter*[nPtBins];
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
    c1->cd(iPad++);
    Int_t origNbins=hmass[iBin]->GetNbinsX();
    fitter[iBin]=new AliHFMassFitter(hmass[iBin],hmin, hmax,rebin[iBin],typeb,types);
    rebin[iBin]=origNbins/fitter[iBin]->GetBinN();
    fitter[iBin]->SetReflectionSigmaFactor(factor4refl);
    fitter[iBin]->SetInitialGaussianMean(massD);
    Bool_t out=fitter[iBin]->MassFitter(0);
    if(!out) continue;
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
    hCntSig2->SetBinContent(iBin+1,cntSig2);
    hCntSig2->SetBinError(iBin+1,TMath::Sqrt(cntErr));
    hSignal->SetBinContent(iBin+1,s);
    hSignal->SetBinError(iBin+1,errs);
    hBackground->SetBinContent(iBin+1,b);
    hBackground->SetBinError(iBin+1,errb);
    hSignificance->SetBinContent(iBin+1,sig);
    hSignificance->SetBinError(iBin+1,errsig);
    
  }
  c1->cd(1); // is some cases the fitting function of 1st bin get lost
  funBckStore1->Draw("same");
  funBckStore2->Draw("same");
  funBckStore3->Draw("same");

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

  TFile* outf=new TFile("RawYield.root","update");
  outf->cd();
  hCntSig1->Write();
  hCntSig2->Write();
  hSignal->Write();
  hBackground->Write();
  hSignificance->Write();
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
    hlist[nReadFiles]=(TList*)dir->Get("coutputDplus");
    TList *listcut = (TList*)dir->Get("coutputDplusCuts");
    dcuts[nReadFiles]=(AliRDHFCutsDplustoKpipi*)listcut->At(1);
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

  TFile* outf=new TFile("RawYield.root","recreate");
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
    if(optPartAntiPart==kParticleOnly) dirname+="D0";
    else if(optPartAntiPart==kAntiParticleOnly) dirname+="D0bar";
    if(cutsappliedondistr) listmassname+="C";

    hlist[nReadFiles]=(TList*)dir->Get(listmassname);

    TString cutsobjname="cutsD0";
    if(optPartAntiPart==kParticleOnly) dirname+="D0";
    else if(optPartAntiPart==kAntiParticleOnly) dirname+="D0bar";
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

  TFile* outf=new TFile("RawYield.root","recreate");
  outf->cd();
  dcuts[0]->Write();
  outf->Close();
  return kTRUE;
}
