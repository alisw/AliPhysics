#include <TROOT.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TObjArray.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TFile.h>
#include <TString.h>
#include <TLegend.h>
#include <TList.h>

#include "AliDielectronSignalExt.h"
#include "AliDielectronCFdraw.h"
#include "AliDielectron.h"

AliDielectronSignalBase* GetSignalLS(AliDielectronCFdraw &d, Int_t step,
                                     AliDielectronSignalBase::EBackgroundMethod type=AliDielectronSignalBase::kLikeSign);
AliDielectronSignalBase* GetSignalRot(AliDielectronCFdraw &d, Int_t step);
void SetStyle(AliDielectronSignalBase *sig, const char* nameAdd);
void DrawSpectra(AliDielectronSignalBase *sig, const char* cname, TH1  *hEventStat=0x0, Bool_t save=kFALSE);
void FitMCshape(AliDielectronSignalBase *sig);

const char *mcLineShapeFile="$ALICE_ROOT/PWGDQ/dielectron/macros/mcMinv_LHC10f7a.root";
TString addToName="";

//_______________________________________
void PlotDataResults(const char* filenameData, const char* filenameMC="", Bool_t save=kFALSE)
{
  if (!addToName.IsNull()) addToName.Prepend("-");
  AliDielectronCFdraw d(filenameData);
  AliDielectronCFdraw dCorr("corrCont","corrCont");
  TString nameCorr(filenameMC);
  if (!nameCorr.IsNull()) d.SetCFContainers(nameCorr.Data());
  TFile f(filenameData);
  TH1 *hStats=(TH1*)f.Get("hEventStat");
  if (!f.IsOpen() || f.IsZombie() || !hStats) return;
  hStats->SetDirectory(0);
  f.Close();
  
  Int_t stepFirst=0, stepAny=1, stepTOFmix=2;
  
  gStyle->SetOptStat(0);
  //Set common Ranges
  d.SetRangeUser("Leg1_NclsTPC",70.,170.);
  d.SetRangeUser("Leg2_NclsTPC",70.,170.);
  d.SetRangeUser("Leg1_Pt",1.01,100000);
  d.SetRangeUser("Leg2_Pt",1.01,100000);
  d.SetRangeUser("Leg1_Eta",-0.899,0.899);
  d.SetRangeUser("Leg2_Eta",-0.899,0.899);

  d.SetRangeUser("Leg1_TPC_nSigma_Electrons",-3.,2.99);
  d.SetRangeUser("Leg2_TPC_nSigma_Electrons",-3.,2.99);
  d.SetRangeUser("Leg1_TPC_nSigma_Pions",3.51,20); 
  d.SetRangeUser("Leg2_TPC_nSigma_Pions",3.51,20); 
  d.SetRangeUser("Leg1_TPC_nSigma_Protons",3.01,20); 
  d.SetRangeUser("Leg2_TPC_nSigma_Protons",3.01,20);

//   d.SetRangeUser("Pt",0,1000);
  
  d.SetRangeUser("M",0.5,5.);
  //============================
  //SPD first
  //
  
  //--- Like sign subtraction
  AliDielectronSignalBase *sigFirst=GetSignalLS(d,stepFirst);
  SetStyle(sigFirst,"ITS First - Like Sign subtraction");
  DrawSpectra(sigFirst,"cFirst",hStats,save);
  //--- Like sign subtraction Arithmetic mean
  AliDielectronSignalBase *sigFirstArith=GetSignalLS(d,stepFirst,AliDielectronSignalBase::kLikeSignArithm);
  SetStyle(sigFirstArith,"ITS FirstArith - Like Sign subtraction");
  DrawSpectra(sigFirstArith,"cFirstArith",hStats,save);
  
  //============================
  //SPD any
  //
  AliDielectronSignalBase *sigAny=GetSignalLS(d,stepAny);
  SetStyle(sigAny,"ITS Any - Like Sign subtraction");
  DrawSpectra(sigAny,"cAny",hStats,save);
  //--- like sign with arithmetic mean
  AliDielectronSignalBase *sigAnyArith=GetSignalLS(d,stepAny,AliDielectronSignalBase::kLikeSignArithm);
  SetStyle(sigAnyArith,"ITS Any - Like Sign subtraction (Arithm. mean)");
  DrawSpectra(sigAnyArith,"cAnyArith",hStats,save);
  
  if (hStats) delete hStats;
}


//_______________________________________
AliDielectronSignalBase *GetSignalLS(AliDielectronCFdraw &d, Int_t step, AliDielectronSignalBase::EBackgroundMethod type)
{
  //
  // Get Extracted signal from likesign method
  //
  
  TObjArray *arr=new TObjArray;
  arr->SetOwner();
  
  for (Int_t iType=0;iType<3;++iType){
    d.SetRangeUser("PairType",iType,iType);
    arr->AddAt(d.Project("M",step),iType);
  }
  
  AliDielectronSignalExt *sig=new AliDielectronSignalExt;
  sig->SetScaleRawToBackground(3.2,4.9);
  sig->SetIntegralRange(2.92,3.15);
  sig->SetMethod(type);
  sig->Process(arr);
  
  delete arr;
  return sig;
}

//_______________________________________
AliDielectronSignalBase *GetSignalRot(AliDielectronCFdraw &d, Int_t step)
{
  //
  // Get Extracted signal from likesign method
  //
  
  TObjArray *arr=new TObjArray;
  arr->SetOwner();
  
  Int_t iType=AliDielectron::kEv1PM;
  d.SetRangeUser("PairType",iType,iType);
  arr->AddAt(d.Project("M",step),iType);
  
  iType=AliDielectron::kEv1PMRot;
  d.SetRangeUser("PairType",iType,iType);
  arr->AddAt(d.Project("M",step),iType);
  
  AliDielectronSignalExt *sig=new AliDielectronSignalExt;
  sig->SetScaleRawToBackground(3.2,4.9);
  sig->SetIntegralRange(2.93,3.15);
  sig->SetMethod(AliDielectronSignalBase::kRotation);
  sig->Process(arr);
  
  
  delete arr;
  return sig;
}

//_______________________________________
void SetStyle(AliDielectronSignalBase *sig, const char* nameAdd)
{
  //
  //
  //
  TH1 *hUS=sig->GetUnlikeSignHistogram();
  hUS->SetMarkerStyle(20);
  hUS->SetMarkerSize(0.7);
  hUS->SetMarkerColor(kRed);
  hUS->SetLineColor(kRed);
  hUS->SetStats(0);
  hUS->SetTitle(Form("Like sign spectrum %s;M inv. ee;counts per %.4g MeV",
                     nameAdd,1000*hUS->GetXaxis()->GetBinWidth(1)));
  
  TH1* hBackground=sig->GetBackgroundHistogram();
  hBackground->SetMarkerStyle(24);
  hBackground->SetMarkerSize(0.7);
  hBackground->SetStats(0);
  hBackground->SetMarkerColor(kBlue);
  hBackground->SetLineColor(kBlue);
  hBackground->SetTitle(Form("Like sign spectrum %s;M inv. ee;counts per %.4g MeV",
                             nameAdd,1000*hBackground->GetXaxis()->GetBinWidth(1)));
  
  TH1* hSignal=sig->GetSignalHistogram();
  hSignal->SetMarkerStyle(20);
  hSignal->SetMarkerSize(0.7);
  hSignal->SetMarkerColor(kRed);
  hSignal->SetLineColor(kRed);
  hSignal->SetTitle(Form("Like sign spectrum %s;M inv. ee;counts per %.4g MeV",
                         nameAdd,1000*hSignal->GetXaxis()->GetBinWidth(1)));
  
  
}

//_______________________________________
void DrawSpectra(AliDielectronSignalBase *sig, const char* cname, TH1  *hEventStat, Bool_t save)
{
  //
  //
  //
  
  
  FitMCshape(sig);
  
  gStyle->SetOptTitle(0);
  TCanvas *c=(TCanvas*)gROOT->FindObject(cname);
  if (!c) c=new TCanvas(cname,cname,400*1.3,500*1.3);
  c->SetTopMargin(0.04);  c->SetRightMargin(0.04);
  c->SetLeftMargin(0.13); c->SetBottomMargin(0.14);
  c->Clear();
  c->Divide(1,2,0,0);
  c->cd(1);
  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.01);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.001);
  /*
  gPad->SetGridx();
  gPad->SetGridy();
*/
  gPad->SetTickx();
  gPad->SetTicky();
  
  TH1 *hUS=sig->GetUnlikeSignHistogram();
  
  TH1* hBackground=sig->GetBackgroundHistogram();
  
  hUS->GetXaxis()->CenterTitle();
  hUS->GetXaxis()->SetTitleSize(0.07); hUS->GetXaxis()->SetLabelSize(0.059);
//   hUS->GetXaxis()->SetLabelOffset(1.);
  hUS->GetXaxis()->SetTitleOffset(1.1);
  hUS->GetXaxis()->SetNdivisions(510);
  hUS->GetYaxis()->CenterTitle();
  hUS->GetYaxis()->SetTitleSize(0.07); hUS->GetYaxis()->SetLabelSize(0.059);
  hUS->GetYaxis()->SetTitleOffset(1.);
  hUS->GetXaxis()->SetLabelFont(42); hUS->GetYaxis()->SetLabelFont(42);
  hUS->GetXaxis()->SetTitleFont(42); hUS->GetYaxis()->SetTitleFont(42);
  hUS->GetYaxis()->SetNdivisions(510);
  //hUS->GetXaxis()->SetLimits(mMinPlot,mMmaxPlot);
  //hUS->SetLineWidth(.7);  
  hUS->SetMarkerSize(.8);
  hUS->SetMarkerColor(2); hUS->SetMarkerStyle(20); hUS->SetLineColor(2);
  hUS->Draw();
  hUS->SetMaximum(hUS->GetMaximum()*1.3);

  //hBackground->SetLineWidth(.4);  
  hBackground->SetMarkerSize(.7);
  hBackground->SetMarkerColor(4); hBackground->SetMarkerStyle(27); hBackground->SetLineColor(4);

  hBackground->Draw("same");
  
  TLatex *lat=new TLatex;
  lat->SetNDC(kTRUE);
  //  lat->DrawLatex(0.68, 0.67, "ALICE Performance");
  lat->SetTextColor(1);lat->SetTextFont(42);lat->SetTextSize(.045);
  
  Double_t sigN=sig->GetSignal();
  Double_t sigEr=sig->GetSignalError();
  Double_t sigS2B=sig->GetSB();
  Double_t sigS2Ber=sig->GetSBError();
  Double_t sigSignif= sig->GetSignificance();
  Double_t sigSignifEr= sig->GetSignificanceError();
  lat->DrawLatex(0.18, 0.92, Form("S: %3d#pm%4.1f, S/B: %3.1f#pm %4.2f, Signif.: %4.1f#pm%4.2f (%4.2f-%4.2f GeV) ",(int)sigN,sigEr,sigS2B,sigS2Ber,sigSignif,sigSignifEr,
                       hUS->GetBinLowEdge(hUS->FindBin(sig->GetIntegralMin())),
                       hUS->GetBinLowEdge(hUS->FindBin(sig->GetIntegralMax())+1)));
  
  TLegend *leg=new TLegend(0.17,0.72,0.42,0.88);
  leg->SetBorderSize(0); leg->SetFillColor(0); leg->SetTextFont(42);
  leg->SetFillStyle(0); leg->SetMargin(0.25); //separation symbol-text
  leg->SetEntrySeparation(0.15);
  leg->AddEntry(hUS,"OS", "p");
  leg->AddEntry(hBackground, Form("LS*%4.2f",sig->GetScaleFactor()), "p");
  leg->Draw();

  c->cd(2);
  gPad->SetRightMargin(0.01);
  gPad->SetTopMargin(0);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(.15);
  /*
  gPad->SetGridx();
  gPad->SetGridy();
  */
  gPad->SetTickx();
  gPad->SetTicky();
  
  TH1* hSignal=sig->GetSignalHistogram();
  hSignal->GetXaxis()->SetTitle("M_{ee} (GeV/c^{2})");
  hSignal->GetXaxis()->CenterTitle();
  hSignal->GetXaxis()->SetTitleSize(0.06); hSignal->GetXaxis()->SetLabelSize(0.05);
  hSignal->GetXaxis()->SetTitleOffset(1.1);
  hSignal->GetXaxis()->SetNdivisions(510);
  hSignal->GetYaxis()->CenterTitle();
  hSignal->GetYaxis()->SetTitleSize(0.06); hSignal->GetYaxis()->SetLabelSize(0.05);
  hSignal->GetYaxis()->SetTitleOffset(1.1);
  hSignal->GetXaxis()->SetLabelFont(42); hSignal->GetYaxis()->SetLabelFont(42);
  hSignal->GetXaxis()->SetTitleFont(42); hSignal->GetYaxis()->SetTitleFont(42);
  hSignal->GetYaxis()->SetNdivisions(510);
  //hSignal->GetXaxis()->SetLimits(mMinPlot,mMmaxPlot);
  //hSignal->SetLineWidth(.7);  
  hSignal->SetMarkerSize(.8);
  hSignal->SetMarkerColor(2); hSignal->SetMarkerStyle(20); hSignal->SetLineColor(2);

  hSignal->Draw();
  
  TLegend *leg2=new TLegend(0.16,0.8,0.57,0.93);
  leg2->SetBorderSize(0); leg2->SetFillColor(0); leg2->SetTextFont(42);
  leg2->SetFillStyle(0); leg2->SetMargin(0.25); //separation symbol-text
  leg2->SetEntrySeparation(0.15);
  leg2->AddEntry(hSignal, Form("OS-%4.2f*LS",sig->GetScaleFactor()), "p");
  TObject *hMmc=hSignal->GetListOfFunctions()->FindObject("mcMinv");
  if (hMmc) leg2->AddEntry(hMmc, Form("MC (#chi^{2}/dof=%3.1f)",TString(hMmc->GetTitle()).Atof()), "l");
  leg2->Draw();
  
  Double_t beforePhys=0;
  Double_t afterPhys=0;
  Double_t afterV0and=0;
  Double_t afterEventSel=0;
  Double_t afterPileupRej=0;

  if (hEventStat){
    Int_t beforePhysBin=hEventStat->GetXaxis()->FindBin("Before Phys. Sel.");
    Int_t afterPhysBin=hEventStat->GetXaxis()->FindBin("After Phys. Sel.");
    Int_t afterV0andBin=hEventStat->GetXaxis()->FindBin("V0and triggers");
    Int_t afterEventSelBin=hEventStat->GetXaxis()->FindBin("After Event Filter");
    Int_t afterPileupRejBin=hEventStat->GetXaxis()->FindBin("After Pileup rejection");
    
    beforePhys=hEventStat->GetBinContent(beforePhysBin);
    afterPhys=hEventStat->GetBinContent(afterPhysBin);
    afterV0and=hEventStat->GetBinContent(afterV0andBin);
    afterEventSel=hEventStat->GetBinContent(afterEventSelBin);
    afterPileupRej=hEventStat->GetBinContent(afterPileupRejBin);
    
    printf("Mevents: all: %4.1f, PhysSel: %4.1f, V0AND: %4.1f, EventSel: %4.1f, PileupRej: %4.1f\n",(Float_t)beforePhys/1.e+6,(Float_t)afterPhys/1.e+6,(Float_t)afterV0and/1.e+6,(Float_t)afterEventSel/1.e+6,(Float_t)afterPileupRej/1.e+6);
    
    lat->DrawLatex(0.18, 0.2, Form("%4.1f Mevents",(Float_t)afterPhys/1.e+6));
  }
  
  if (save){
//     c->SaveAs(Form("%s%s.eps",cname,addToName.Data()));
    c->SaveAs(Form("%s%s.png",cname,addToName.Data()));
/*
    FILE *out_file;
    if ( (out_file = fopen(Form("sig_%s.txt",cname), "w")) == NULL )
    {   fprintf(stderr, "Cannot open file %s\n", Form("sig_%s.txt",cname)); }
    fprintf(stdout, "Signal file: %s\n", Form("sig_%s.txt",cname));
    fprintf(out_file,"%3d %4.1f  %3.1f %4.2f  %4.1f %4.2f %d\n",(int)sigN,sigEr,sigS2B,sigS2Ber,sigSignif,sigSignifEr,(Int_t)afterPhys);
    fclose(out_file);

    TFile outMinv(Form("Minv_%s.root",cname), "RECREATE");
    hUS->Write();
    hBackground->Write();
    hSignal->Write();
    if (hMmc) hMmc->Write();
    outMinv.Close();*/

  }
  
  return;
}

//_______________________________________
void FitMCshape(AliDielectronSignalBase *sig)
{
  TFile mcFile(mcLineShapeFile);
  if (!mcFile.IsOpen()) {
    printf("mcMinv_LHC10e2 not found!!!\n");
    return;
  }
  TH1D *hMmc = (TH1D*)mcFile.Get("mcMinv");
  if (!hMmc) {
    printf("mcMinv not found!!\n");
    return;
  }
  hMmc->SetDirectory(0);
  hMmc->SetName("mcMinv");
  mcFile.Close();
  
  TH1* hMsub=sig->GetSignalHistogram();
  Double_t mb1=sig->GetIntegralMin();
  Double_t mb2=sig->GetIntegralMax();
  
  Double_t effInt = 0.;
  for(Int_t iBin=hMmc->FindBin(mb1); iBin<hMmc->FindBin(mb2); iBin++) {
    effInt += hMmc->GetBinContent(iBin);
  }
  effInt/=hMmc->Integral();
  printf("MC signal fraction in range %4.2f-%4.2f GeV: %5.3f \n",hMmc->GetBinLowEdge(hMmc->FindBin(mb1)),hMmc->GetBinLowEdge(hMmc->FindBin(mb2)+1),effInt);
 
  Float_t mcScale1=(hMsub->GetXaxis()->GetBinWidth(1)/hMmc->GetXaxis()->GetBinWidth(1))*
    hMsub->Integral(hMsub->FindBin(mb1),hMsub->FindBin(mb2))/
    hMmc->Integral(hMmc->FindBin(mb1),hMmc->FindBin(mb2));
  
  printf("1st guess of MC scale factor: %6.3f\n",mcScale1);
  
  Float_t mcScale=0.;
  Float_t chi2_min=100000.;
  Int_t iMin=0;
  Int_t ndf=0;
  
  for(Int_t i=0; i<20; i++){
    Float_t chi2=0.;
    Float_t scale=(0.4+0.05*(Float_t)i)*mcScale1;
    ndf=0;
    for(Int_t ib=1; ib<=hMsub->GetXaxis()->GetNbins(); ib++){
      Float_t data=(Float_t)hMsub->GetBinContent(ib);
      Float_t err=(Float_t)hMsub->GetBinError(ib);
      Float_t mc=scale*((Float_t)hMmc->GetBinContent(hMmc->FindBin(hMsub->GetBinCenter(ib))));
      if (err>0) {
        chi2 += ((data-mc)*(data-mc))/(err*err);
        ndf++;
      } else {
        //printf("bin %d Err: %6.3f, chi2: %6.1f\n",ib,err,chi2);
      }
    }
      //printf("%d scale factor: %6.3f, chi2: %6.1f\n",i,scale,chi2);
    if(chi2 < chi2_min){
      chi2_min = chi2;
      mcScale = scale;
      iMin=i;
    }
  }
  //Float_t chi2dof=chi2_min/(Float_t)(hMinv->GetXaxis()->GetNbins()-1);
  Float_t chi2dof=chi2_min/((Float_t)(ndf-1));
  printf("MC fit (i=%d): chi2/dof: %6.3f/%d, Scale: %7.4f \n",iMin,chi2_min,(ndf-1),mcScale);
  hMmc->SetTitle(Form("%f",chi2dof));
  
  //mcScale=IntData/IntMC;printf("Int Data, MC: %10.1f %10.1f, MC scale: %6.3f\n",IntData,IntMC,mcScale);
  
  hMmc->Scale(mcScale);
  hMmc->SetOption("sameHISTC");
  hMsub->GetListOfFunctions()->Add(hMmc);
}
