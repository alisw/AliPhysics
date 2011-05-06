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
#include <TF1.h>

#include "AliDielectronSignalExt.h"
#include "AliDielectronCFdraw.h"
#include "AliDielectron.h"
#include "AliDielectronHelper.h"

AliDielectronSignalBase* GetSignalLS(AliDielectronCFdraw &d, Int_t step,
                                     AliDielectronSignalBase::EBackgroundMethod type=AliDielectronSignalBase::kLikeSign);
AliDielectronSignalBase* GetSignalLSY(AliDielectronCFdraw &d, Int_t step,
                                     AliDielectronSignalBase::EBackgroundMethod type=AliDielectronSignalBase::kLikeSign);
AliDielectronSignalBase* GetSignalRot(AliDielectronCFdraw &d, Int_t step);
void SetStyle(AliDielectronSignalBase *sig, const char* nameAdd);
void DrawSpectra(AliDielectronSignalBase *sig, const char* cname, TH1  *hEventStat=0x0, Float_t fNevSel=1., Bool_t save=kTRUE);
Float_t FitMCshape(AliDielectronSignalBase *sig);

const char *mcLineShapeFile="$ALICE_ROOT/PWG3/dielectron/macros/mcMinv_LHC10f7a.root";
 
//_______________________________________
void PlotDataResults(const char* filenameData, const char* filenameMC="", Bool_t save=kTRUE, Bool_t DoPt=kTRUE )
{
  AliDielectronCFdraw d(filenameData);
  AliDielectronCFdraw dCorr("corrCont","corrCont");
  TString nameCorr(filenameMC);
  if (!nameCorr.IsNull()) d.SetCFContainers(nameCorr.Data());

  TFile f(filenameData);
  TH1 *hStats=(TH1*)f.Get("hEventStat");
  hStats->SetDirectory(0);
  TList* listData = (TList*)f.Get("jpsi_QA");
  //THashList *a = (THashList*)listData->FindObject("basicQ+SPDfirst+pt>.6+PID");
  THashList *a = (THashList*)listData->FindObject("basicQ+SPDfirst+pt>.8+PID"); //01.02.11
  //THashList *a = (THashList*)listData->FindObject("FullCuts+SPDany+Pi3.0+P3.0"); //Ionut
  THashList *b = (THashList*)a->FindObject("Event");    
  TH1F *hZvtx= (TH1F*) b->FindObject("VtxZ");
  hZvtx->SetDirectory(0);
  f.Close();
  
  hZvtx->GetXaxis()->SetTitleSize(.05); hZvtx->GetXaxis()->SetLabelSize(.038);
  hZvtx->GetXaxis()->SetTitleOffset(1.1);
  hZvtx->GetXaxis()->SetNdivisions(510);
  //hZvtx->GetYaxis()->CenterTitle();  hZvtx->GetXaxis()->CenterTitle();
  hZvtx->GetYaxis()->SetTitleSize(.05); hZvtx->GetYaxis()->SetLabelSize(.038);
  hZvtx->GetYaxis()->SetTitleOffset(1.3); //1.0
  hZvtx->GetXaxis()->SetLabelFont(42); hZvtx->GetYaxis()->SetLabelFont(42);
  hZvtx->GetXaxis()->SetTitleFont(42); hZvtx->GetYaxis()->SetTitleFont(42);
  hZvtx->GetYaxis()->SetNdivisions(510);
  hZvtx->SetLineWidth(2.5);  
  hZvtx->SetMarkerSize(.8);
  hZvtx->SetMarkerColor(2); hZvtx->SetMarkerStyle(20); hZvtx->SetLineColor(2);
  hZvtx->GetXaxis()->SetTitle("Zvtx (cm)");
  hZvtx->GetYaxis()->SetTitle("Events");
  
  TCanvas *c0=new TCanvas("c0","c0",10,0,450,450);
  c0->SetTopMargin(0.04);  c0->SetRightMargin(0.03);
  c0->SetLeftMargin(0.13); c0->SetBottomMargin(0.14);
  
  TF1* gausFit=new TF1("gausFit", "gaus",-9.,9.);
  gausFit->SetLineColor(4); gausFit->SetLineWidth(.7);// gausFit->Draw("same");
  hZvtx->Fit(gausFit, "QM", "", -9.,9.);
  gStyle->SetOptFit(0);
  hZvtx->Draw();
  
  //Float_t IntGaus=gausFit->Integral(gausFit->GetParameter(1)-2.0*gausFit->GetParameter(2),gausFit->GetParameter(1)+2.0*gausFit->GetParameter(2))/h->GetBinWidth(1);
  Float_t NevSel=hZvtx->Integral();
  Float_t IntGaus=gausFit->Integral(-25.,25.)/hZvtx->GetBinWidth(1);
  printf("NevSel: %6.1f, IntGauss: %6.1f \n",NevSel/1.e+6,IntGaus/1.e+6);
  Float_t fNevSel=NevSel/IntGaus;
  
  TLatex *lat=new TLatex();
  lat->SetNDC(kTRUE);
  lat->SetTextColor(2);lat->SetTextFont(42);lat->SetTextSize(.035);
  lat->DrawLatex(0.16, 0.9, Form("NeventsSel: %5.1f M",NevSel/1.e+6));
  
  TLatex *lat2=new TLatex();
  lat2->SetNDC(kTRUE);
  lat2->SetTextColor(4);lat2->SetTextFont(42);lat2->SetTextSize(.035);
  lat2->DrawLatex(0.7, 0.9, "Gauss fit:");
  lat2->DrawLatex(0.7, 0.86, Form("Nevents: %5.1f M",IntGaus/1.e+6));
  lat2->DrawLatex(0.7, 0.82, Form("Mean: %5.3f cm",gausFit->GetParameter(1)));
  lat2->DrawLatex(0.7, 0.78, Form("Sigma: %4.2f cm",gausFit->GetParameter(2)));
  
  Int_t afterPhysBin=hStats->GetXaxis()->FindBin("After Phys. Sel.");
  Double_t afterPhys=hStats->GetBinContent(afterPhysBin);
  lat2->DrawLatex(0.4, 0.18, Form("NevPhysSel: %5.1f",(Float_t)afterPhys/1.e+6));

  c0->SaveAs("Zvtx_cSel.eps");

  //if (hZvtx) delete hZvtx;

  Int_t stepFirst=0, stepAny=1, stepAnyV0=2, stepAnyV0T=3; //Data ...from 05.02.11 on
  //Int_t stepFirst=0, stepAny=1, stepTOFmix=2; //Data
  //Int_t stepFirst=2, stepAny=4, stepTOFmix=6; d.SetRangeUser("PairType",1,1); //MC only
  //Int_t stepFirst=1, stepAny=3, stepAnyV0T=5; //MC truth
  //Int_t stepFirst=2, stepAny=4, stepAnyV0T=5; //MC, all pairs ...or vice-versa

  gStyle->SetOptStat(0);
  //Set common Ranges
  d.SetRangeUser("Leg1_NclsTPC",70.,170.); //was 90
  d.SetRangeUser("Leg2_NclsTPC",70.,170.);
  d.SetRangeUser("Leg1_Pt",1.01,100000);
  d.SetRangeUser("Leg2_Pt",1.01,100000);
  d.SetRangeUser("Leg1_Eta",-0.899,0.899);
  d.SetRangeUser("Leg2_Eta",-0.899,0.899);

//   d.SetRangeUser("Leg1_TPC_nSigma_Electrons",-3,3);
//   d.SetRangeUser("Leg2_TPC_nSigma_Electrons",-3,3);
  d.SetRangeUser("Leg1_TPC_nSigma_Pions",3.51,20); // 2011-02-09_0027.5195 only
  d.SetRangeUser("Leg2_TPC_nSigma_Pions",3.51,20); 
  d.SetRangeUser("Leg1_TPC_nSigma_Protons",3.01,20); 
  d.SetRangeUser("Leg2_TPC_nSigma_Protons",3.01,20);
  d.SetRangeUser("M",0.5,5.);
  //============================
  //SPD first
  //
  
  //--- Like sign subtraction
  AliDielectronSignalBase *sigFirst=GetSignalLS(d,stepFirst);
  SetStyle(sigFirst,"SPDfirst - Like Sign subtraction");
  DrawSpectra(sigFirst,"cFirst",hStats,fNevSel,save);
  //--- Like sign subtraction Arithmetic mean
  AliDielectronSignalBase *sigFirstArith=GetSignalLS(d,stepFirst,AliDielectronSignalBase::kLikeSignArithm);
  SetStyle(sigFirstArith,"SPDfirstArith - Like Sign subtraction");
  DrawSpectra(sigFirstArith,"cFirstArith",hStats,fNevSel,save);
  //--- Rotation subtraction
  AliDielectronSignalBase *sigFirstRot=GetSignalRot(d,stepFirst);
  SetStyle(sigFirstRot,"SPDfirst - Track rotation subtraction");
  DrawSpectra(sigFirstRot,"cFirstRot",hStats,save);
  
  //============================
  //SPD any
  //
  AliDielectronSignalBase *sigAny=GetSignalLS(d,stepAny);
  SetStyle(sigAny,"ITS Any - Like Sign subtraction");
  DrawSpectra(sigAny,"cAny",hStats,fNevSel,save);
  //--- like sign with arithmetic mean
  AliDielectronSignalBase *sigAnyArith=GetSignalLS(d,stepAny,AliDielectronSignalBase::kLikeSignArithm);
  SetStyle(sigAnyArith,"SPDany - Like Sign subtraction (Arithm. mean)");
  DrawSpectra(sigAnyArith,"cAnyArith",hStats,fNevSel,save);
  
  //--- Rotation subtraction
  AliDielectronSignalBase *sigAnyRot=GetSignalRot(d,stepAny);
  SetStyle(sigAnyRot,"SPDany - Track rotation subtraction");
  DrawSpectra(sigAnyRot,"cAnyRot",hStats,fNevSel,save);

  //AliDielectronSignalBase *sigAnyRot2=GetSignalRot(d,stepAnyV0); //diff. rot. angle
  //SetStyle(sigAnyRot2,"ITS First - Track rotation subtraction");
  //DrawSpectra(sigAnyRot2,"cAnyRot2",hStats,fNevSel,save);

  //AliDielectronSignalBase *sigAnyRot3=GetSignalRot(d,stepAnyV0T);  //diff. rot. angle
  //SetStyle(sigAnyRot3,"ITS First - Track rotation subtraction");
  //DrawSpectra(sigAnyRot3,"cAnyRot3",hStats,fNevSel,save);

  //--- like sign with arithmetic mean, V0 (tender) conversions excl.
  //AliDielectronSignalBase *sigAnyV0tArith=GetSignalLS(d,stepAnyV0T,AliDielectronSignalBase::kLikeSignArithm);
  //SetStyle(sigAnyV0tArith,"ITS Any V0 (Tender) excl. - LS-Arithm.");
  //DrawSpectra(sigAnyV0tArith,"cAnyV0tArith",hStats,fNevSel,save);

  //=============================
  //TOF up to 1.2, parametrisation in TPC ele
  //
/*
  AliDielectronSignalBase *sigTOFmix=GetSignalLS(d,stepTOFmix);
  SetStyle(sigTOFmix,"TOF + TPC - Like Sign subtraction");
  DrawSpectra(sigTOFmix,"cTOFTPC",hStats,fNevSel,save);
  //--- Rotation subtraction
  AliDielectronSignalBase *sigTOFmixRot=GetSignalRot(d,stepTOFmix);
  SetStyle(sigTOFmixRot,"TOF + TPC - Track rotation subtraction");
  DrawSpectra(sigTOFmixRot,"cTOFTPCrot",hStats,fNevSel,save);
*/  

 //================================
 // y bins
  AliDielectronSignalBase *sigY_03_09=GetSignalLSY(d,stepAny,AliDielectronSignalBase::kLikeSignArithm);
  DrawSpectra(sigY_03_09,"cAny_Y03_09",hStats,fNevSel,save);

  d.SetRangeUser("Y",-.299,+.299);
  AliDielectronSignalBase *sigY_03=GetSignalLS(d,stepAny,AliDielectronSignalBase::kLikeSignArithm);
  DrawSpectra(sigY_03,"cAny_Y03",hStats,fNevSel,save);

//  simulator.AliModule::SetDensityFactor(1.07);

  if (DoPt){

  //===============================
  // Pt bins
  //   "0.0, 0.8, 1.4, 2.8, 5., 9.9"
  d.SetRangeUser("Y",-.899,+.899);
  //Char_t *metPt="cAnyArith"; //method for background for Minv in pt bins
  Char_t *metPt="cAnyRot";


  FILE *out_file;
  if (save){      
      if ( (out_file = fopen(Form("sigPt_%s.txt",metPt), "w")) == NULL )
      { fprintf(stderr, "Cannot open file %s\n", Form("sigPt_%s.txt",metPt)); }
  }  
  const Int_t NptBins=6;
  //Float_t PtC[5]={0.5,1.25,2.25,4.,7.5}; //pt bins...need clever way
  //Float_t PtW[5]={0.5,0.25,0.75,1.0,2.5};
  //Float_t PtC[5]={0.5,1.5,2.5,4.,6.5}; //pt bins...23.01.11
  //Float_t PtW[5]={0.5,0.5,0.5,1.0,1.5};
  Float_t PtC[NptBins]={0.5,1.5,2.5,4.,6.,8.5}; //pt bins...01.02.11
  Float_t PtW[NptBins]={0.5,0.5,0.5,1.,1.,1.5};
  

  //TVectorD *vBins=AliDielectronHelper::MakeArbitraryBinning("0.0, 1, 1.5, 3, 5., 9.9");
  //TVectorD *vBins=AliDielectronHelper::MakeArbitraryBinning("0.0, 1, 2, 3, 5, 8"); //26.01  
  TVectorD *vBins=AliDielectronHelper::MakeArbitraryBinning("0.0,1,2,3,5,7,10"); //01.02


  //AliDielectronSignalBase *sig_Pt[5];
  AliDielectronSignalBase *sig_Pt[NptBins];

  //for (Int_t i=0; i<5; ++i){
  for (Int_t i=0; i<NptBins; ++i){
    d.SetRangeUser("Pt",(*vBins)[i]+.001, (*vBins)[i+1]-0.001);
    if ("cAnyArith" == metPt){	
	sig_Pt[i]=GetSignalLS(d,stepAny,AliDielectronSignalBase::kLikeSignArithm);
    } 
    else if ("cAnyRot" == metPt){
	sig_Pt[i]=GetSignalRot(d,stepAny); //rotation
    }
    else {
	cout << " ??? No such method ??? " <<metPt <<endl;
    }
    DrawSpectra(sig_Pt[i],Form("%s_Pt%d",metPt,i),hStats,fNevSel,save);
  }
  
  TH1D *hSigPt=new TH1D("hSigPt","hSigPt",NptBins,vBins->GetMatrixArray());
  TH1D *hSigSign=new TH1D("hSigSign","hSigSign",NptBins,vBins->GetMatrixArray());
  TH1D *hSigSOB=new TH1D("hSigSOB","hSigSOB",NptBins,vBins->GetMatrixArray());

  for (Int_t i=0; i<NptBins; ++i){
      hSigPt->SetBinContent(i+1,sig_Pt[i]->GetSignal());
      hSigPt->SetBinError(i+1,sig_Pt[i]->GetSignalError());
      hSigSign->SetBinContent(i+1,sig_Pt[i]->GetSignificance());
      hSigSign->SetBinError(i+1,sig_Pt[i]->GetSignificanceError());
      hSigSOB->SetBinContent(i+1,sig_Pt[i]->GetSB());
      hSigSOB->SetBinError(i+1,sig_Pt[i]->GetSBError());
      if (save){ 
	  fprintf(out_file,"%4.1f %5.2f %5.1f  %4.1f  %3.1f %4.2f  %4.1f %4.2f \n",PtC[i],PtW[i],hSigPt->GetBinContent(i+1),hSigPt->GetBinError(i+1),hSigSign->GetBinContent(i+1),hSigSign->GetBinError(i+1),hSigSOB->GetBinContent(i+1),hSigSOB->GetBinError(i+1));
      }
  
  }

  if (save){      
      fclose(out_file);      
      fprintf(stdout, " *** Signal file vs. pt: %s\n", Form("sigPt_%s.txt",metPt));
  }

  const Int_t kMarkTyp=20; //marker type
  const Int_t kMarkCol=2; //...and color
  const Float_t kTitSize=0.08; //axis title size
  const Float_t kAsize=0.85*kTitSize; //...and label size
  const Float_t kToffset=0.8;
  const Int_t kFont=42;

  TCanvas *cSigPt=(TCanvas*)gROOT->FindObject("cSigPt");

  if (!cSigPt) cSigPt=new TCanvas("cSigPt","cSigPt",5,5,450,700);
  cSigPt->Clear();
  //cSigPt->Divide(2,2);
  cSigPt->SetTopMargin(0.04);  cSigPt->SetRightMargin(0.04);
  cSigPt->SetLeftMargin(0.13);  cSigPt->SetBottomMargin(0.2);
  cSigPt->Divide(1,3,0,0);
  cSigPt->cd(1);

  hSigPt->SetMarkerStyle(kMarkTyp);
  hSigPt->SetMarkerColor(kMarkCol); hSigPt->SetLineColor(kMarkCol);
  hSigPt->GetYaxis()->SetTitleOffset(kToffset);
  hSigPt->SetTitleSize(kTitSize,"XY");
  hSigPt->SetTitleFont(kFont,"XY");
  hSigPt->SetLabelSize(kAsize,"YX");
  hSigPt->SetLabelFont(kFont,"XY");
  //gS->SetAxisRange(0.,Pt2[Npt-1],"X"); 
  //hSigPt->GSetAxisRange(Y1m,Y1M,"Y"); 
  hSigPt->GetYaxis()->SetTitle("N_{J/#psi}");
  hSigPt->Draw();


  cSigPt->cd(2);
  hSigSign->SetMarkerStyle(kMarkTyp);
  hSigSign->SetMarkerColor(kMarkCol); hSigSign->SetLineColor(kMarkCol);
  hSigSign->GetYaxis()->SetTitleOffset(kToffset);
  hSigSign->SetTitleSize(kTitSize,"XY");
  hSigSign->SetTitleFont(kFont,"XY");
  hSigSign->SetLabelSize(kAsize,"YX");
  hSigSign->SetLabelFont(kFont,"XY");
  hSigSign->GetYaxis()->SetTitle("Significance");
  hSigSign->Draw();

  cSigPt->cd(3);
  hSigSOB->SetMaximum(4.); 
  hSigSOB->SetMinimum(0.); 
  hSigSOB->SetMarkerStyle(kMarkTyp);
  hSigSOB->SetMarkerColor(kMarkCol); hSigSOB->SetLineColor(kMarkCol);
  hSigSOB->GetYaxis()->SetTitleOffset(kToffset);
  hSigSOB->SetTitleSize(kTitSize,"XY");
  hSigSOB->SetTitleFont(kFont,"XY");
  hSigSOB->SetLabelSize(kAsize,"YX");
  hSigSOB->SetLabelFont(kFont,"XY");
  hSigSOB->GetXaxis()->SetTitle("p_{t}^{J/#psi} (GeV/c");
  hSigSOB->GetYaxis()->SetTitle("S/B");
  hSigSOB->Draw();  

  cSigPt->SaveAs(Form("SigPt_%s.eps",metPt));
  
  } //DoPt

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
    TH1D* hh = (TH1D*)d.Project("M",step);
    if(TMath::Abs(hh->GetBinWidth(1)-0.02)<0.0001) {
	//hh->Rebin(2);
    }
    arr->AddAt(hh,iType);
  }
  
  AliDielectronSignalExt *sig=new AliDielectronSignalExt;
  sig->SetScaleRawToBackground(3.2,4.9);
  //sig->SetScaleRawToBackground(3.2,5.0);
  sig->SetIntegralRange(2.92,3.15);
  sig->SetMethod(type);
  sig->Process(arr);
  
  delete arr;
  return sig;
}


//_______________________________________
AliDielectronSignalBase *GetSignalLSY(AliDielectronCFdraw &d, Int_t step, AliDielectronSignalBase::EBackgroundMethod type)
{
  //
  // Get Extracted signal from likesign method
  //

  TObjArray *arr=new TObjArray;
  arr->SetOwner();

  d.SetRangeUser("Y",-.899,-.301);

  for (Int_t iType=0;iType<3;++iType){
    d.SetRangeUser("PairType",iType,iType);
    arr->AddAt(d.Project("M",step),iType);
  }

  d.SetRangeUser("Y",+.301,+.899);

  for (Int_t iType=0;iType<3;++iType){
    d.SetRangeUser("PairType",iType,iType);
    TH1 *h=d.Project("M",step);
    
    TH1 *h2=(TH1*)arr->At(iType);
    h2->Add(h);
    delete h;
  }

  AliDielectronSignalExt *sig=new AliDielectronSignalExt;
  sig->SetScaleRawToBackground(3.2,4.9);
  sig->SetIntegralRange(2.93,3.15);
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
void DrawSpectra(AliDielectronSignalBase *sig, const char* cname, TH1  *hEventStat, Float_t fNevSel, Bool_t save)
{
  //
  //
  //
  Float_t effInt=FitMCshape(sig);
  
  gStyle->SetOptTitle(0);
  TCanvas *c=(TCanvas*)gROOT->FindObject(cname);
  if (!c) c=new TCanvas(cname,cname,400,500);
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
  lat->DrawLatex(0.18, 0.92, Form("S: %5.1f#pm%4.1f, S/B: %3.1f#pm %4.2f, Signif.: %4.1f#pm%4.2f (%4.2f-%4.2f GeV) ",sigN,sigEr,sigS2B,sigS2Ber,sigSignif,sigSignifEr,
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
    
    afterPhys*=fNevSel; //fraction after Zvtx cut

    lat->DrawLatex(0.18, 0.2, Form("%4.1f Mevents",(Float_t)afterPhys/1.e+6));
  }  

  if (save){ 

    FILE *out_file;
    if ( (out_file = fopen(Form("sig_%s.txt",cname), "w")) == NULL )
    {   fprintf(stderr, "Cannot open file %s\n", Form("sig_%s.txt",cname)); }
    fprintf(stdout, "Signal file: %s\n", Form("sig_%s.txt",cname));
    fprintf(out_file,"%3d %4.1f  %3.1f %4.2f  %4.1f %4.2f  %d  %4.2f  %6.3f\n",(int)sigN,sigEr,sigS2B,sigS2Ber,sigSignif,sigSignifEr,(Int_t)afterPhys,sig->GetScaleFactor(),effInt);
    fclose(out_file);

    TFile outMinv(Form("Minv_%s.root",cname), "RECREATE");
    hUS->Write();
    hBackground->Write();
    hSignal->Write();
    if (hMmc) hMmc->Write();
    outMinv.Close();

    c->SaveAs(Form("Minv_%s.eps",cname));
    //c->SaveAs(Form("%s.png",cname));

  }
  
  return;
}

//_______________________________________
Float_t FitMCshape(AliDielectronSignalBase *sig)
{
  TFile mcFile(mcLineShapeFile);
  if (!mcFile.IsOpen()) {
    printf("MC Minv file not found!!!\n");
    return 0;
  }
  //TH1D *hMmc = (TH1D*)mcFile.Get("mcMinv");
  TH1D *hMmc = (TH1D*)mcFile.Get("HistSignal");
  if (!hMmc) {
    printf("mcMinv not found!!\n");
    return 0;
  }
  hMmc->SetDirectory(0);
  hMmc->SetName("mcMinv");
  mcFile.Close();
  
  TH1* hMsub=sig->GetSignalHistogram();
  Double_t mb1=sig->GetIntegralMin();
  Double_t mb2=sig->GetIntegralMax();
  
  Double_t effInt = 0.;
  for(Int_t iBin=hMmc->FindBin(mb1); iBin<=hMmc->FindBin(mb2); iBin++) {
    effInt += hMmc->GetBinContent(iBin);
  }
  effInt/=hMmc->Integral();
  printf("MC signal fraction in range %4.2f-%4.2f GeV: %5.3f \n",hMmc->GetBinLowEdge(hMmc->FindBin(mb1)),hMmc->GetBinLowEdge(hMmc->FindBin(mb2)+1),effInt);
 
  Float_t MminFit=1.2; //min., Minv for fit (and chi2 calc.)
  Float_t mcScale1=(hMsub->GetXaxis()->GetBinWidth(1)/hMmc->GetXaxis()->GetBinWidth(1))*
    hMsub->Integral(hMsub->FindBin(mb1),hMsub->FindBin(mb2))/
    hMmc->Integral(hMmc->FindBin(mb1),hMmc->FindBin(mb2));
  
  printf("1st guess of MC scale factor: %6.3f ...chi2 for Minv>%4.1f\n",mcScale1,MminFit);
  
  Float_t mcScale=0.;
  Float_t chi2_min=100000.;
  Int_t iMin=0;
  Int_t ndf=0;
  
  for(Int_t i=0; i<20; i++){
    Float_t chi2=0.;
    Float_t scale=(0.4+0.05*(Float_t)i)*mcScale1;
    ndf=0;
    for(Int_t ib=1; ib<=hMsub->GetXaxis()->GetNbins(); ib++){
	if (hMsub->GetBinCenter(ib) > MminFit) {
	    Float_t data=(Float_t)hMsub->GetBinContent(ib);
	    Float_t err=(Float_t)hMsub->GetBinError(ib);
	    Float_t mc=scale*((Float_t)hMmc->GetBinContent(hMmc->FindBin(hMsub->GetBinCenter(ib))));
	    if (err>0) {
		chi2 += ((data-mc)*(data-mc))/(err*err);
		ndf++;
	    } else {
		printf("bin %d Err: %6.3f, chi2: %6.1f\n",ib,err,chi2);
	    }
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
  hMmc->SetLineColor(1); hMmc->SetMarkerColor(1);
  hMsub->GetListOfFunctions()->Add(hMmc);

  return effInt;

}
