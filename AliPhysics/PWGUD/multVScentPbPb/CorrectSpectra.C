#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TList.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TF1.h"
#include "TLine.h"
#include "TPaletteAxis.h"
//
//
#endif
#include "/home/shahoian/ALICE/mymacros/SaveCanvas.h"



const char kHStatName[]="hStat";
double kEps = 1e-6;
//
double kdPhiBgTailMin = 0.1; // lower limit of dphi tail to use for bg normalization  
double kdPhiBgTailMax = 0.3; // upper limit of dphi tail to use for bg normalization   
//
double kWDistBgTailMin = 5.; // lower limit of wgh.distance tail to use for bg normalization  
double kWDistBgTailMax = 25.; // upper limit of wgh.distance tail to use for bg normalization

double kdPhiSgCut = 0.06;     // cut in dphi-bent used to extract the signal, extracted from stat histo
double kWDistSgCut = 1.5;        // cut in w.distance used to extract the signal, extracted from stat histo
//
enum { kNormShapeDist,      // normalize bg tails usig weighted distance shape 
       kNormShapeDPhi,      // normalize bg tails usig dPhi-bend shape
       kNormShapes};

enum { kSclWghMean,         // normalize bg tails to data using weighted mean of bin-by-bin ratios 
       kSclIntegral,        // normalize bg tails to data using integral
       kSclTypes};

// histograms to be used for bg nomalization for each of NormShapes used
const char* kNormShapeH[kNormShapes] = {
  "EtaDist",                // Weighted distance vs Eta
  "EtaDPhiS"                // Dphi-bent vs distance
}; 

const char* figDir = "fig161110";
TString  useBgType    = "Inj";
Int_t    useShapeType = kNormShapeDist;    // which distribution to use for bg normalization
Bool_t   useMCLB      = 0;//kFALSE;             // use Comb MC Labels as a template for Bg.
Int_t    useScaleType = kSclIntegral;//kSclWghMean;       // which type of tails normalization to use
Double_t useEtaCut    = 1.;                // cut on eta
Double_t useZvMin     = -5.;                // cut on Z vertex
Double_t useZvMax     =  5.;                // cut on Z vertex
Int_t    useBinGrEta  = -1;                // for bg fits group eta bins
Int_t    useBinGrZv   = -1;                // for bg fits group Zv bins

enum {kBitNormPerEvent=BIT(14)};
  // bins for saved parameters in the hStat histo
enum {kDummyBin,
      kEvTot,       // events read
      kEvProcData,  // events with data mult.object (ESD or reco)
      kEvProcInj,   // events Injected
      kEvProcRot,   // events Rotated
      kEvProcMix,   // events Mixed
      //
      kDPhi,        // dphi window
      kDTht,        // dtheta window
      kNStd,        // N.standard deviations to keep
      kPhiShift,    // bending shift
      kThtS2,       // is dtheta scaled by 1/sin^2
      kPhiOvl,      // overlap params
      kZEtaOvl,     // overlap params
      kNoOvl,       // flag that overlap are suppressed
      //
      kPhiRot,      // rotation phi
      kInjScl,      // injection scaling
      kEtaCut,      // eta cut
      kZVMin,       // min ZVertex to process
      kZVMax,       // max ZVertex to process
      kTrcMin,      // min mult to process
      kTrcMax,      // max mult to process
      //
      kOneUnit=49,  // just 1 to track mergings
      kNWorkers=50, // n workers
      kNStatBins
};


enum {kSigCorr,kMCPrim,kRawDtCut,kSignalEst,kSignalEstMC,kBgEst,k1MBeta,k1MBetaMC,kAlpha,kAlphaMC,kBgMC,kBgRescFc,kDataDist,kBgDist,kBgMCDist, kNHistos, kMCShift=20};

void    CorrectSpectra(const char* flNameData, const char* flNameMC,const char* unique="");
void    PrepareHistos(TList* lst, Bool_t isMC, TObjArray* outArray);
void    ProcessHistos(TObjArray* outArray);
TH1*    NormalizeBg(TH1* dataH, TH1* bgH, double &scl, double &scle);
TObject* FindObject(const char* nameH, const TList* lst, Bool_t normPerEvent=kTRUE);
TList*  LoadList(const char* flName, const char* addPref, const char* nameL="clist");
void    GetRatE(double x,double xe, double y,double ye, double &rat, double &rate);
Int_t   CheckStat(const TList* lst,const char* dtType);
void    Integrate(TH1* hist, double xmn,double xmx, double &val, double& err);
void    CropHisto(TH1* histo, int b00, int b01, int b10=-1, int b11=-1);
void    CropHisto(TH1* histo, double b00, double b01, double b10=-1, double b11=-1);
void    GetRealMinMax(TH1* h, double &vmn, double &vmx);
const char*   HName(const char* prefix,const char* htype);

void ProjDZE(THnSparseF* hrawd, THnSparseF* hgenb, THnSparseF* hmcComb, TObjArray* res, int bstep0=-1,int bstep1=-1);

void PlotResults();
void PlotDNDEta();
void PlotDist();
void PlotAlphaBeta();
void PlotSpecies();

TList *listDt=0, *listMC=0;
TObjArray resArr;
char outStr[1000];
char outTitle[1000];
TString uniqueName="";

void CorrectSpectra(const char* flNameData, const char* flNameMC, const char* uniqueNm)
{
  //
  uniqueName = uniqueNm;
  listDt = LoadList(flNameData,"dt_");
  listMC = LoadList(flNameMC,"mc_");
  //
  resArr.Clear();
  PrepareHistos(listMC, kTRUE,  &resArr);
  PrepareHistos(listDt, kFALSE, &resArr);
  //
  ProcessHistos(&resArr); 
  //
  sprintf(outStr,"CutEta%.1f_Zv%.1f_%.1f_bg%s_Shape_%s_mcLB%d_cutSig%.1f_cutBg%.1f",useEtaCut,useZvMin,useZvMax,useBgType.Data(),
	  useShapeType==kNormShapeDist ? "wdst":"dphi", 
	  useMCLB,
	  useShapeType==kNormShapeDist ? kWDistSgCut:kdPhiSgCut,
	  useShapeType==kNormShapeDist ? kWDistBgTailMin:kdPhiBgTailMin);
  sprintf(outTitle,"%s, |#eta|<%.1f, %.1f<Z_{V}<%.1f,  Bg.:%s, UseMCLB=%d, CutVar:%s, |sig|<%.2f, %.2f<|bg.nrm|<%.2f",
	  uniqueName.Data(),
	  useEtaCut,useZvMin,useZvMax,useBgType.Data(),
	  useMCLB,
	  useShapeType==kNormShapeDist ? "#Delta":"#Delta#varphi-#delta_{#varphi}",
	  useShapeType==kNormShapeDist ? kWDistSgCut:kdPhiSgCut,
	  useShapeType==kNormShapeDist ? kWDistBgTailMin : kdPhiBgTailMin,
	  useShapeType==kNormShapeDist ? kWDistBgTailMax : kdPhiBgTailMax	  
	  );
  PlotResults();
}

//_____________________________________________________________________
void PrepareHistos(TList* lst, Bool_t isMC, TObjArray* outArray)
{
  // params:
  //
  int shift = isMC ? kMCShift : 0;
  //
  const char* xxZvEta;
  if (useShapeType==kNormShapeDist) xxZvEta = "DistZvEta";
  else                              xxZvEta = "DistZvDPhiS";

  printf("PrepareBg  : (%4s) of type %s from %s with shape %d\n",isMC?" MC ":"Data",useBgType.Data(),lst->GetName(),useShapeType);
  if (!CheckStat(lst,useBgType.Data())) {printf("Bg of type %s is absent in list %s\n",useBgType.Data(),lst->GetName()); return;}
  //
  THnSparseF* hD = (THnSparseF*) FindObject(HName("Data",xxZvEta),lst);
  THnSparseF* hB = (THnSparseF*) FindObject(HName(useBgType.Data(),xxZvEta),lst);
  //
  // special feature: use MC Labels bg as a shape instead of generated bg
  if (useMCLB/* && !isMC*/) {
    TString nm  = hB->GetName();   nm  += "_MCLB";
    TString tit = hB->GetTitle();  tit += "_MCLB";
    THnSparseF* hBMCLB = (THnSparseF*) FindObject(HName("Comb",xxZvEta),listMC)->Clone(nm.Data());
    hBMCLB->SetTitle(tit.Data());
    hB = hBMCLB;
  }
  THnSparseF* hBmc = 0;
  if (isMC) hBmc = (THnSparseF*) FindObject(HName("Comb",xxZvEta),lst);
  //
  ProjDZE(hD,hB,hBmc, outArray, useBinGrEta, useBinGrZv);
  //
  if (!isMC) return;
  //
  // prepare MC primary signal histo
  TH2F* mcPrim = (TH2F*)FindObject( "zvEtaPrimMC", lst );
  mcPrim = (TH2F*) mcPrim->Clone("mcTrueSignal");
  CropHisto(mcPrim,-useEtaCut,useEtaCut,useZvMin,useZvMax);
  outArray->AddAtAndExpand(mcPrim, kMCPrim + shift);
  //
}

//_____________________________________________________________________
void ProcessHistos(TObjArray* outArray) 
{
  //
  // build alpha matrix
  TH2* halp = (TH2*)outArray->At(kMCShift + kMCPrim);
  halp = (TH2*) halp->Clone("Alpha");
  halp->SetTitle("#alpha");
  halp->Divide( (TH2*)outArray->At(kMCShift + k1MBeta) );
  halp->Divide( (TH2*)outArray->At(kMCShift + kRawDtCut) );
  halp->SetMinimum(1.5);
  halp->SetMaximum(4.);
  outArray->AddAtAndExpand(halp, kAlpha);
  //
  // build alpha matrix with MC labels bg
  TH2* halpMC = (TH2*)outArray->At(kMCShift + kMCPrim);
  halpMC = (TH2*) halpMC->Clone("AlphaMC");
  halpMC->SetTitle("#alpha MC labels");
  halpMC->Divide( (TH2*)outArray->At(kMCShift + k1MBetaMC) );
  halpMC->Divide( (TH2*)outArray->At(kMCShift + kRawDtCut) );
  halpMC->SetMinimum(1.5);
  halpMC->SetMaximum(4.);
  outArray->AddAtAndExpand(halpMC, kAlphaMC);
  //
  // build corrected signal
  TH2* hsigCorr = (TH2*)outArray->At(kSignalEst);
  hsigCorr = (TH2*) hsigCorr->Clone("SignalEstCorr");
  hsigCorr->SetTitle("Corrected Signal");
  hsigCorr->Multiply( halp );
  outArray->AddAtAndExpand(hsigCorr, kSigCorr);
  //
}

TCanvas *canvFin=0;
TCanvas *canvAlp=0;
TCanvas *canvBet=0;
TCanvas* canvDst=0;
TCanvas* canvSpec=0;

void PlotResults() 
{
  PlotDist();
  PlotAlphaBeta();
  PlotSpecies();
  PlotDNDEta();

}

void PlotDNDEta()
{
  //
  TObjArray *res = &resArr;
  char buff[1000];
  // eta range
  double plEta = useEtaCut*1.1;
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  double mn = 1e6,mx = -1e6;
  canvFin = new TCanvas("canvFin", "canvFin",0,50,700,550);
  canvFin->SetLeftMargin(0.15);
  //  canvFin->ToggleEventStatus();
  //
  // corrected data
  TH1* hsigCorr = ((TH2F*)res->At(kSigCorr))->ProjectionX("DataCorrSignal");
  SetHStyle(hsigCorr,kRed,20,1.0);
  hsigCorr->GetXaxis()->SetRangeUser(-plEta+kEps,plEta-kEps);
  hsigCorr->Scale(1./hsigCorr->GetBinWidth(1));
  hsigCorr->Draw();
  mx = TMath::Max(mx, hsigCorr->GetMaximum());
  mn = TMath::Min(mn, hsigCorr->GetMinimum());
  TF1* pl0 = new TF1("pl0","pol0");
  pl0->SetParameter(0,hsigCorr->GetMinimum());
  hsigCorr->Fit(pl0,"q0","",-0.5,.5);
  double fval = pl0->GetParameter(0);
  double ferr = pl0->GetParError(0);
  char ftres[1000];
  sprintf(ftres,"dN/d#eta_{|#eta|<0.5} = %d #pm %d",int(fval),int(ferr));
  printf("dN/d#eta_{|#eta|<0.5} = %.2f  %.2f\n",fval,ferr);
  TLatex *txfit = new TLatex(-0.2,hsigCorr->GetMinimum()*0.9, ftres);
  txfit->SetTextSize(0.04);
  txfit->Draw();
  //
  // raw data
  TH1* hraw = ((TH2F*)res->At(kRawDtCut))->ProjectionX("DataRaw");
  SetHStyle(hraw,kRed,21,1.0);
  hraw->GetXaxis()->SetRangeUser(-plEta,plEta);
  hraw->Scale(1./hraw->GetBinWidth(1));
  hraw->Draw("same");
  mn = TMath::Min(mn, hraw->GetMinimum());
  mx = TMath::Max(mx, hraw->GetMaximum());
  //  
  // raw data bg sub
  TH1* hraws = ((TH2F*)res->At(kSignalEst))->ProjectionX("DataRawSub");
  SetHStyle(hraws,kRed,23,1.0);
  hraws->GetXaxis()->SetRangeUser(-plEta,plEta);
  hraws->Scale(1./hraw->GetBinWidth(1));
  hraws->Draw("same");
  mn = TMath::Min(mn, hraw->GetMinimum());
  mx = TMath::Max(mx, hraw->GetMaximum());
  //
  // bg
  TH1* hbg = ((TH2F*)res->At(kBgEst))->ProjectionX("BgEst");
  SetHStyle(hbg,kMagenta,22,1.0);
  hbg->GetXaxis()->SetRangeUser(-plEta,plEta);
  hbg->Scale(1./hbg->GetBinWidth(1));
  hbg->Draw("same");
  mn = TMath::Min(mn, hbg->GetMinimum());
  mx = TMath::Max(mx, hbg->GetMaximum());
  //
  // mc part ----------------------------
  // raw data
  TH1* hrawMC = ((TH2F*)res->At(kRawDtCut+kMCShift))->ProjectionX("DataRawMC");
  SetHStyle(hrawMC,kBlue,24,1.0);
  hrawMC->GetXaxis()->SetRangeUser(-plEta,plEta);
  hrawMC->Scale(1./hrawMC->GetBinWidth(1));
  hrawMC->Draw("same");
  mn = TMath::Min(mn, hrawMC->GetMinimum());
  mx = TMath::Max(mx, hrawMC->GetMaximum());
  //  
  // raw data bg sub
  TH1* hrawsMC = ((TH2F*)res->At(kSignalEst+kMCShift))->ProjectionX("DataRawSubMC");
  SetHStyle(hrawsMC,kBlue,26,1.0);
  hrawsMC->GetXaxis()->SetRangeUser(-plEta,plEta);
  hrawsMC->Scale(1./hrawMC->GetBinWidth(1));
  hrawsMC->Draw("same");
  mn = TMath::Min(mn, hrawMC->GetMinimum());
  mx = TMath::Max(mx, hrawMC->GetMaximum());
  //
  // raw data bgMClabels sub
  TH1* hrawsMCLB = ((TH2F*)res->At(kSignalEstMC+kMCShift))->ProjectionX("DataRawSubMCLB");
  SetHStyle(hrawsMCLB,kGreen+2,30,1.0);
  hrawsMCLB->GetXaxis()->SetRangeUser(-plEta,plEta);
  hrawsMCLB->Scale(1./hrawsMCLB->GetBinWidth(1));
  hrawsMCLB->Draw("same");
  mn = TMath::Min(mn, hrawsMCLB->GetMinimum());
  mx = TMath::Max(mx, hrawsMCLB->GetMaximum());
  //
  // bg est
  TH1* hbgMCEst = ((TH2F*)res->At(kBgEst+kMCShift))->ProjectionX("BgEstMC");
  SetHStyle(hbgMCEst,kBlue,26,1.0);
  hbgMCEst->GetXaxis()->SetRangeUser(-plEta,plEta);
  hbgMCEst->Scale(1./hbgMCEst->GetBinWidth(1));
  hbgMCEst->Draw("same");
  mn = TMath::Min(mn, hbgMCEst->GetMinimum());
  mx = TMath::Max(mx, hbgMCEst->GetMaximum());
  //  
  // bg MC
  TH1* hbgMC = ((TH2F*)res->At(kBgMC+kMCShift))->ProjectionX("BgMC");
  SetHStyle(hbgMC,kGreen+2,25,1.0);
  hbgMC->GetXaxis()->SetRangeUser(-plEta,plEta);
  hbgMC->Scale(1./hbgMC->GetBinWidth(1));
  hbgMC->Draw("same");
  mn = TMath::Min(mn, hbgMC->GetMinimum());
  mx = TMath::Max(mx, hbgMC->GetMaximum());
  //
  mn = 0;
  hsigCorr->SetMinimum(mn);
  hsigCorr->SetMaximum(mx + 0.4*(mx-mn));
  gPad->Modified();
  //
  TLegend *legDnDeta = new TLegend(0.15,0.75, 0.45,0.95);
  legDnDeta->SetFillColor(kWhite);
  legDnDeta->SetHeader("Data");
  legDnDeta->AddEntry(hsigCorr,"Corrected","pl");
  legDnDeta->AddEntry(hraw,    "Reconstructed","pl");
  sprintf(buff,"Reconstructed - Bckg.%s.",useBgType.Data());
  legDnDeta->AddEntry(hraws,   buff,"pl");
  sprintf(buff,"Background %s.",useBgType.Data());
  legDnDeta->AddEntry(hbg,     buff,"pl");
  legDnDeta->Draw();
  //
  TLegend *legDnDetaMC = new TLegend(0.60,0.72, 0.95,0.95);
  legDnDetaMC->SetFillColor(kWhite);
  legDnDetaMC->SetHeader("MC");
  legDnDetaMC->AddEntry(hrawMC,    "Reconstructed","pl");
  sprintf(buff,"Reconstructed - Bckg.%s.",useBgType.Data());
  legDnDetaMC->AddEntry(hrawsMC,   buff,"pl");
  sprintf(buff,"Reconstructed - Bckg.%s.","MC.Labels");
  legDnDetaMC->AddEntry(hrawsMCLB,   buff,"pl");
  sprintf(buff,"Background %s.",useBgType.Data());
  legDnDetaMC->AddEntry(hbgMCEst,     buff,"pl");
  sprintf(buff,"Background %s.","MC Labels");
  legDnDetaMC->AddEntry(hbgMC,     buff,"pl");

  legDnDetaMC->Draw();
  //
  gPad->SetGrid(1.1);
  gPad->Modified();
  canvFin->cd();
  AddLabel(outTitle,0.1,0.97, kBlack,0.025);
  sprintf(buff,"%s/%sdNdEta_%s",figDir,uniqueName.Data(),outStr);
  SaveCanvas(canvFin,buff,"cg");
  //
}
//
void PlotAlphaBeta() 
{
  char buff[1000];
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TObjArray* res = &resArr;
  //------------------------------------------------------
  canvBet = new TCanvas("canvBet","canvBet",10,10,1000,400);
  canvBet->Divide(3,1,0.01,0.06);
  canvBet->cd(1);
  TH1* dtBet = (TH1*)res->At(k1MBeta);
  TH1* mcBet = (TH1*)res->At(k1MBeta+kMCShift);
  TH1* mcBetLB = (TH1*)res->At(k1MBetaMC+kMCShift);
  double mn,mx,mnt,mxt;
  GetRealMinMax(dtBet,mn,mx);
  GetRealMinMax(mcBet,mnt,mxt);
  if (mnt<mn) mn = mnt;
  if (mxt>mx) mx = mxt;
  GetRealMinMax(mcBetLB,mnt,mxt);
  if (mnt<mn) mn = mnt;
  if (mxt>mx) mx = mxt;
  //
  dtBet->SetMinimum(mn - 0.05*(mx-mn));
  dtBet->SetMaximum(mx + 0.05*(mx-mn));
  mcBet->SetMinimum(mn - 0.05*(mx-mn));
  mcBet->SetMaximum(mx + 0.05*(mx-mn));
  mcBetLB->SetMinimum(mn - 0.05*(mx-mn));
  mcBetLB->SetMaximum(mx + 0.05*(mx-mn));
  //
  canvBet->cd(1);
  gPad->SetRightMargin(0.15);
  dtBet->Draw("colz");
  AddLabel("#beta Data",0.2,0.95,kBlack,0.05);
  gPad->Modified();
  dtBet->GetYaxis()->SetTitleOffset(1.4);
  TPaletteAxis *p = (TPaletteAxis*)dtBet->FindObject("palette");
  if (p) p->SetX1NDC(0.85);
  canvBet->cd(2);
  gPad->SetRightMargin(0.15);
  mcBet->Draw("colz");
  AddLabel("#beta MC (bckg.estimated)",0.2,0.95,kBlack,0.05);
  gPad->Modified();
  mcBet->GetYaxis()->SetTitleOffset(1.4);
  p = (TPaletteAxis*)mcBet->FindObject("palette");
  if (p) p->SetX1NDC(0.85);
  canvBet->cd(3);
  gPad->SetRightMargin(0.15);
  mcBetLB->Draw("colz");
  AddLabel("#beta MC (bckg.from MC labels)",0.2,0.95,kBlack,0.05);
  gPad->Modified();
  mcBetLB->GetYaxis()->SetTitleOffset(1.4);
  p = (TPaletteAxis*)mcBetLB->FindObject("palette");
  if (p) p->SetX1NDC(0.85);
  //
  canvBet->cd();
  AddLabel(outTitle,0.1,0.97, kBlack, 0.025);
  //
  sprintf(buff,"%s/%sBeta_%s",figDir,uniqueName.Data(),outStr);
  SaveCanvas(canvBet,buff,"cg");

  //------------------------------------------------------
  canvAlp = new TCanvas("canvAlp","canvAlp",10,10,900,400);
  canvAlp->Divide(2,1,0.01,0.06);
  canvAlp->cd(1);
  TH1* dtAlp = (TH1*)res->At(kAlpha);
  TH1* mcAlp = (TH1*)res->At(kAlphaMC);
  GetRealMinMax(dtAlp,mn,mx);
  GetRealMinMax(mcAlp,mnt,mxt);
  if (mnt<mn) mn = mnt;
  if (mxt>mx) mx = mxt;
  dtAlp->SetMinimum(mn - 0.05*(mx-mn));
  dtAlp->SetMaximum(mx + 0.05*(mx-mn));
  mcAlp->SetMinimum(mn - 0.05*(mx-mn));
  mcAlp->SetMaximum(mx + 0.05*(mx-mn));
  //
  canvAlp->cd(1);
  gPad->SetRightMargin(0.15);
  dtAlp->Draw("colz");
  AddLabel("#alpha (bckg.estimated)",0.2,0.95,kBlack,0.05);
  gPad->Modified();
  dtAlp->GetYaxis()->SetTitleOffset(1.4);
  TPaletteAxis *pa = (TPaletteAxis*)dtBet->FindObject("palette");
  if (pa) pa->SetX1NDC(0.85);
  canvAlp->cd(2);
  gPad->SetRightMargin(0.15);
  mcAlp->Draw("colz");
  AddLabel("#alpha (bckg.from MC labels)",0.2,0.95,kBlack,0.05);
  gPad->Modified();
  mcAlp->GetYaxis()->SetTitleOffset(1.4);
  pa = (TPaletteAxis*)mcBet->FindObject("palette");
  if (pa) pa->SetX1NDC(0.85);
  gPad->Modified();
  canvAlp->cd();
  AddLabel(outTitle,0.1,0.97, kBlack, 0.025);
  //
  sprintf(buff,"%s/%sAlpha_%s",figDir,uniqueName.Data(),outStr);
  SaveCanvas(canvAlp,buff,"cg");
  
}

void PlotSpecies() 
{
  char buff[1000];
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //------------------------------------------------------
  TH2F* hSpecPrim  = (TH2F*)FindObject( "pdgPrim", listMC,kFALSE);
  TH2F* hSpecSec   = (TH2F*)FindObject( "pdgSec", listMC,kFALSE);
  TH2F* hSpecPrimP = (TH2F*)FindObject( "pdgPrimPar", listMC,kFALSE);
  TH2F* hSpecSecP  = (TH2F*)FindObject( "pdgSecPar", listMC,kFALSE);
  int nbd = hSpecPrim->GetXaxis()->GetNbins();
  //
  TH1* hSpecPrimAll = hSpecPrim->ProjectionX("specPrimAll",1,nbd+1,"e");
  hSpecPrimAll->Scale(100./hSpecPrimAll->Integral());
  hSpecPrimAll->GetYaxis()->SetTitle("Fraction,%");
  hSpecPrimAll->GetXaxis()->SetLabelSize(0.06);
  hSpecPrimAll->GetXaxis()->LabelsOption("v");
  //
  TH1* hSpecSecAll  = hSpecSec->ProjectionX("specSecAll",1,nbd+1,"e");
  hSpecSecAll->Scale(100./hSpecSecAll->Integral());
  hSpecSecAll->GetYaxis()->SetTitle("Fraction,%");
  hSpecSecAll->GetXaxis()->SetLabelSize(0.05);
  //
  TH1* hSpecPrimPAll = hSpecPrimP->ProjectionX("specPrimPAll",1,nbd+1,"e");
  hSpecPrimPAll->Scale(100./hSpecPrimPAll->Integral());
  hSpecPrimPAll->GetYaxis()->SetTitle("Fraction,%");
  hSpecPrimPAll->GetXaxis()->SetLabelSize(0.06);
  hSpecPrimPAll->GetXaxis()->LabelsOption("v");

  //
  TH1* hSpecSecPAll  = hSpecSecP->ProjectionX("specSecPAll",1,nbd+1,"e");
  hSpecSecPAll->Scale(100./hSpecSecPAll->Integral());
  hSpecSecPAll->GetYaxis()->SetTitle("Fraction,%");
  hSpecSecPAll->GetXaxis()->SetLabelSize(0.05);
  //
  int binCut = hSpecPrim->GetXaxis()->FindBin(kWDistSgCut-kEps);
  TH1* hSpecPrimSel = hSpecPrim->ProjectionX("specPrimSel",1,binCut,"e");
  hSpecPrimSel->Scale(100./hSpecPrimSel->Integral());
  hSpecPrimSel->GetYaxis()->SetTitle("Fraction,%");
  hSpecPrimSel->GetXaxis()->SetLabelSize(0.05);
  //
  TH1* hSpecSecSel  = hSpecSec->ProjectionX("specSecSel",1,binCut,"e");
  hSpecSecSel->Scale(100./hSpecSecSel->Integral());
  hSpecSecSel->GetYaxis()->SetTitle("Fraction,%");
  hSpecSecSel->GetXaxis()->SetLabelSize(0.05);
  //
  TH1* hSpecPrimPSel = hSpecPrimP->ProjectionX("specPrimPSel",1,binCut,"e");
  hSpecPrimPSel->Scale(100./hSpecPrimPSel->Integral());
  hSpecPrimPSel->GetYaxis()->SetTitle("Fraction,%");
  hSpecPrimPSel->GetXaxis()->SetLabelSize(0.05);
  //
  TH1* hSpecSecPSel  = hSpecSecP->ProjectionX("specSecPSel",1,binCut,"e");
  hSpecSecPSel->Scale(100./hSpecSecPSel->Integral());
  hSpecSecPSel->GetYaxis()->SetTitle("Fraction,%");
  hSpecSecPSel->GetXaxis()->SetLabelSize(0.05);


  canvSpec = new TCanvas("canvSpec","canvSpec",10,10,1100,800);
  canvSpec->Divide(1,2,0.01,0.01);
  canvSpec->cd(1);
  hSpecPrimAll->Draw();
  SetHStyle(hSpecPrimAll,kBlue,25,1.1);
  hSpecPrimSel->Draw("same");
  SetHStyle(hSpecPrimSel,kRed,20,1);
  //
  hSpecSecAll->Draw("same");
  SetHStyle(hSpecSecAll,kGreen,32,1.1);
  hSpecSecSel->Draw("same");
  SetHStyle(hSpecSecSel,kBlack,22,1);
  //
  TLegend *legPart = new TLegend(0.8,0.72, 0.999,0.999);
  legPart->SetFillColor(kWhite);
  legPart->SetHeader("Tracklet PDG");
  //
  legPart->AddEntry(hSpecPrimAll,    "Prim., before #Delta cut","pl");
  legPart->AddEntry(hSpecPrimSel,    "Prim., after #Delta cut","pl");
  legPart->AddEntry(hSpecSecAll,     "Sec., before #Delta cut","pl");
  legPart->AddEntry(hSpecSecSel,     "Sec., after #Delta cut","pl");
  //
  legPart->Draw();
  gPad->SetLogy();
  gPad->SetGrid(1,1);
  gPad->Modified();
  //
  canvSpec->cd(2);
  hSpecPrimPAll->Draw();
  SetHStyle(hSpecPrimPAll,kBlue,25,1.1);
  hSpecPrimPSel->Draw("same");
  SetHStyle(hSpecPrimPSel,kRed,20,1);
  //
  hSpecSecPAll->Draw("same");
  SetHStyle(hSpecSecPAll,kGreen,32,1.1);
  hSpecSecPSel->Draw("same");
  SetHStyle(hSpecSecPSel,kBlack,22,1);
  //
  TLegend *legPartP = new TLegend(0.8,0.72, 0.999,0.999);
  legPartP->SetFillColor(kWhite);
  legPartP->SetHeader("Tracklet Parents PDG");
  //
  legPartP->AddEntry(hSpecPrimPAll,    "Prim., before #Delta cut","pl");
  legPartP->AddEntry(hSpecPrimPSel,    "Prim., after #Delta cut","pl");
  legPartP->AddEntry(hSpecSecPAll,     "Sec., before #Delta cut","pl");
  legPartP->AddEntry(hSpecSecPSel,     "Sec., after #Delta cut","pl");
  //
  legPartP->Draw();
  gPad->SetLogy();
  gPad->SetGrid(1,1);
  gPad->Modified();
  //
  canvSpec->cd(1);
  AddLabel(outTitle,0.1,0.97, kBlack, 0.025);
  canvSpec->cd();
  //
  sprintf(buff,"%s/%sSpecies_%s",figDir,uniqueName.Data(),outStr);
  SaveCanvas(canvSpec,buff,"cg");
}

void PlotDist()
{
  TObjArray* res = &resArr;
  char buff[1000];
  canvDst = new TCanvas("canvDst","canvDst",10,10,700,500);
  TH1* mcdst = (TH1*)res->At(kDataDist+kMCShift);
  TH1* mcdstbg = (TH1*)res->At(kBgDist+kMCShift);
  TH1* mcdstbgLB = (TH1*)res->At(kBgMCDist+kMCShift);
  TH1* dtdst = (TH1*)res->At(kDataDist);
  TH1* dtdstbg = (TH1*)res->At(kBgDist);
  //
  TH2* mcDstZN     = (TH2*)FindObject(useShapeType==kNormShapeDist  ?  "TrData_ZvDist":"TrData_ZvDPhiS", listMC );
  TH2* mcDstZSec   = (TH2*)FindObject(useShapeType==kNormShapeDist  ?  "TrSec_ZvDist":"TrSec_ZvDPhiS", listMC );
  TH2* mcDstZCombU = (TH2*)FindObject(useShapeType==kNormShapeDist  ?  "TrCombU_ZvDist":"TrCombU_ZvDPhiS", listMC );
  TH2* mcDstZCombC = (TH2*)FindObject(useShapeType==kNormShapeDist  ?  "TrComb_ZvDist":"TrComb_ZvDPhiS", listMC );
  //
  mcDstZN->GetXaxis()->SetRangeUser(useZvMin+kEps,useZvMax-kEps);
  mcDstZSec->GetXaxis()->SetRangeUser(useZvMin+kEps,useZvMax-kEps);
  mcDstZCombU->GetXaxis()->SetRangeUser(useZvMin+kEps,useZvMax-kEps);
  mcDstZCombC->GetXaxis()->SetRangeUser(useZvMin+kEps,useZvMax-kEps);
  //
  TH1* mcDstN  = mcDstZN->ProjectionY("mcDstN");
  TH1* mcDstSec  = mcDstZSec->ProjectionY("mcDstSec");
  TH1* mcDstCombU = mcDstZCombU->ProjectionY("mcDstCombU");
  TH1* mcDstCombC = mcDstZCombC->ProjectionY("mcDstCombC");
  //
  double scl,sclE;
  mcDstN = NormalizeBg(mcdst,mcDstN,scl,sclE);
  mcDstSec->Scale(scl);
  mcDstCombU->Scale(scl);
  mcDstCombC->Scale(scl);
  mcDstCombC->Add(mcDstCombU,-1);
  //
  dtdst->Draw("");
  gPad->Modified();
  dtdst->GetXaxis()->SetLabelSize(0.03);
  dtdst->GetXaxis()->SetTitleSize(0.03);
  dtdst->GetXaxis()->SetTitleOffset(2);
  dtdstbg->Draw("same");

  mcdst->Draw("same");
  mcDstSec->Draw("same");
  mcdstbgLB->Draw("same");
  mcdstbg->Draw("same");
  mcDstCombC->Draw("same");
  //

  SetHStyle(mcdst,kBlue, 25,0.7);
  SetHStyle(mcdstbgLB,kGreen, 7/*32*/,0.5);
  SetHStyle(mcdstbg,kCyan, 7/*26*/,0.5);
  SetHStyle(mcDstCombC,kGreen+2, 21,0.7);
  SetHStyle(mcDstSec,kBlue+2, 22,0.7);
  //
  SetHStyle(dtdst,kRed, 20,0.7);
  SetHStyle(dtdstbg,kBlue, 34,0.7);
  //
  double vmcTot,vmcTotE;
  double vmcSec,vmcSecE, ratSec,ratSecE;  
  double vmcCmbEst,vmcCmbEstE, ratCmbEst,ratCmbEstE;
  double vmcCmb,vmcCmbE, ratCmb,ratCmbE;
  double vmcCmbC,vmcCmbCE, ratCmbC,ratCmbCE;  
  double cutSgMin,cutSgMax;
  double cutBgMin,cutBgMax;
  if (useShapeType==kNormShapeDist) {
    cutSgMin = 0;
    cutSgMax = kWDistSgCut;
    cutBgMin = kWDistBgTailMin;
    cutBgMax = kWDistBgTailMax;
  }
  else {
    cutSgMin =  0;
    cutSgMax =  kdPhiSgCut;
    cutBgMin =  kdPhiBgTailMin;
    cutBgMax =  kdPhiBgTailMax;
  }  
  Integrate(mcdst,    cutSgMin,cutSgMax, vmcTot,vmcTotE);     
  Integrate(mcDstSec, cutSgMin,cutSgMax, vmcSec,vmcSecE);
  GetRatE(vmcSec,vmcSecE, vmcTot,vmcTotE, ratSec,ratSecE);
  //
  Integrate(mcdstbgLB, cutSgMin,cutSgMax, vmcCmb,vmcCmbE);  
  GetRatE(vmcCmb,vmcCmbE, vmcTot,vmcTotE, ratCmb,ratCmbE);
  //
  Integrate(mcdstbg,  cutSgMin,cutSgMax, vmcCmbEst,vmcCmbEstE); 
  GetRatE(vmcCmbEst,vmcCmbEstE, vmcTot,vmcTotE, ratCmbEst,ratCmbEstE);
  //
  Integrate(mcDstCombC,  cutSgMin,cutSgMax, vmcCmbC,vmcCmbCE);  
  GetRatE(vmcCmbC,vmcCmbCE, vmcTot,vmcTotE, ratCmbC,ratCmbCE);
  //
  double vdtTot,vdtTotE;
  double vdtBg,vdtBgE, ratdtBg,ratdtBgE;  
  //  
  Integrate(dtdst,    cutSgMin,cutSgMax, vdtTot,vdtTotE);     
  Integrate(dtdstbg,  cutSgMin,cutSgMax, vdtBg,vdtBgE);     
  GetRatE( vdtBg,vdtBgE,  vdtTot,vdtTotE, ratdtBg,ratdtBgE);
  //
  //
  double dmn = mcdst->GetMinimum();
  double dmx = mcdst->GetMaximum();
  TLine *ln = new TLine(cutSgMax, dmn, cutSgMax, dmx);
  ln->SetLineColor(kBlack);
  ln->Draw();
  TLine *lnc = new TLine(cutBgMin, dmn, cutBgMin, dmx);
  lnc->SetLineColor(kRed);
  lnc->Draw();
  if (useShapeType==kNormShapeDPhi) {
    ln = new TLine(-cutSgMax, dmn, -cutSgMax, dmx);
    ln->SetLineColor(kBlack);
    ln->Draw();
    //
    lnc = new TLine(-cutBgMin, dmn, -cutBgMin, dmx);
    lnc->SetLineColor(kRed);
    lnc->Draw();
  }
  //
  TLegend *legDstMC1 = new TLegend(0.60,0.72, 0.95,0.95);
  legDstMC1->SetFillColor(kWhite);

  //
  legDstMC1->AddEntry(dtdst,    "Data raw","pl");
  sprintf(buff,"Data Comb. %s. : %.1f%%",useBgType.Data(),ratdtBg*100);
  legDstMC1->AddEntry(dtdstbg,  buff,"pl");
  //


  legDstMC1->AddEntry(mcdst,    "MC raw","pl");
  sprintf(buff,"MC Comb. %s. : %.1f%%",useBgType.Data(),ratCmbEst*100);
  legDstMC1->AddEntry(mcdstbg,  buff,"pl");
  //
  sprintf(buff,"MC Comb. %s. : %.1f%%","MC Lbl.",ratCmb*100);
  legDstMC1->AddEntry(mcdstbgLB,  buff,"pl");

  sprintf(buff,"MC Comb.Uncorr %s. : %.1f%%","MC Lbl.",ratCmbC*100);
  legDstMC1->AddEntry(mcDstCombC,  buff,"pl");

  sprintf(buff,"MC Sec.   : %.1f%%",ratSec*100);
  legDstMC1->AddEntry(mcDstSec,  buff,"pl");

  legDstMC1->Draw();

  if (useShapeType==kNormShapeDist) gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetGrid(1,1);
  gPad->Modified();
  //
  canvDst->cd();
  AddLabel(outTitle,0.1,0.97, kBlack, 0.025);
  //
  sprintf(buff,"%s/%sDst_%s",figDir,uniqueName.Data(),outStr);
  SaveCanvas(canvDst,buff,"cg");
  //
}


//______________________________________________________________________
void CropHisto(TH1* histo, int bx0, int bx1, int by0, int by1)	       
{
  // fill 0 to all bins outside defined range
  TAxis *xax = histo->GetXaxis();
  int nbx = xax->GetNbins(); 
  double vmn=1e16,vmx=-1e16;
  if (histo->InheritsFrom(TH2::Class())) {
    TAxis *yax = histo->GetYaxis();
    int nby = yax->GetNbins();
    for (int ix=nbx+2;ix--;) {
      for (int iy=nby+2;iy--;) {
	if ((ix<bx0||ix>bx1)||(iy<by0||iy>by1)) {
	  histo->SetBinContent(ix,iy,0);
	  histo->SetBinError(ix,iy,0);
	}
	else {
	  double vl = histo->GetBinContent(ix,iy);
	  if (vl<vmn) vmn = vl;
	  if (vl>vmx) vmx = vl;
	}
      }
    }
  }
  else {
    for (int ix=nbx+2;ix--;) {
      if ((ix<bx0||ix>bx1)) {
	histo->SetBinContent(ix,0);
	histo->SetBinError(ix,0);
      }
      else {
	double vl = histo->GetBinContent(ix);
	if (vl<vmn) vmn = vl;
	if (vl>vmx) vmx = vl;
      }
    }
  }
  //
  if (vmn==vmx) {
    vmn = 0.95*vmn;
    vmx = 1.05*vmx;
  }
  histo->SetMaximum(vmx);
  histo->SetMinimum(vmn);
}

//______________________________________________________________________
void CropHisto(TH1* histo, double vx0, double vx1, double vy0, double vy1)	       
{
  // fill 0 to all bins outside defined range
  TAxis *xax = histo->GetXaxis();
  int bx0,bx1,by0=-1,by1=-1;
  bx0 = xax->FindBin(vx0+kEps);
  bx1 = xax->FindBin(vx1-kEps);
  if (histo->InheritsFrom(TH2::Class())) {
    TAxis *yax = histo->GetYaxis();
    by0 = yax->FindBin(vy0+kEps);
    by1 = yax->FindBin(vy1-kEps);
  }
  CropHisto(histo,bx0,bx1,by0,by1);
}

//______________________________________________________________________
TH1* NormalizeBg(TH1* dataH, TH1* bgH, double &scl, double &sclE)
{
  // match generated bg and data tails, calculate normalization, return normalized bg copy
  //
  TAxis* xax = dataH->GetXaxis();
  int nbtot = xax->GetNbins();
  int bgBins[2][2] = {{0}}; // limiting bins for tails integration
  Int_t ntails; // 0 for dphi plot, 1 for weighted dist plot  
  if  (useShapeType == kNormShapeDist) { // only positive tail
    bgBins[0][0] = xax->FindBin(kWDistBgTailMin+kEps); // positive tail min bin
    bgBins[0][1] = xax->FindBin(kWDistBgTailMax-kEps); // positive tail max bin
    ntails = 1;
  }
  else if (useShapeType == kNormShapeDPhi) {         // both tails
    bgBins[0][0] = xax->FindBin( kdPhiBgTailMin+kEps); // positive tail min bin
    bgBins[0][1] = xax->FindBin( kdPhiBgTailMax-kEps); // positive tail max bin
    bgBins[1][0] = xax->FindBin(-kdPhiBgTailMax+kEps); // negative tail min bin
    bgBins[1][1] = xax->FindBin(-kdPhiBgTailMin-kEps); // positive tail max bin    
    ntails = 2;
  }
  else {printf("NormalizeBg: unknown shape type %d\n",useShapeType);exit(1);}
  printf("NormalizeBg: bins for tails: right: %d:%d / left: %d:%d\n",bgBins[0][0],bgBins[0][1],bgBins[1][0],bgBins[1][1]);
  // 
  double meanR=0,meanRE=0,meanRE2=0;
  double meanD=0,meanDE2=0;
  double meanB=0,meanBE2=0;
  double meanRI=0,meanRIE=0;
  for (int itp=0;itp<=ntails;itp++) {
    for (int ib=bgBins[itp][0];ib<=bgBins[itp][1];ib++) {
      if (ib<1||ib>nbtot) continue;
      double vD = dataH->GetBinContent(ib);
      double vB = bgH->GetBinContent(ib);
      double eD = dataH->GetBinError(ib);
      double eB = bgH->GetBinError(ib);
      meanD += vD; meanDE2 += eD*eD;
      meanB += vB; meanBE2 += eB*eB;
      if (vD<=0 || vB<=0 || eD<=0 || eB<=0) continue;
      double rat = vD/vB;
      double ratE2 = rat*rat*(eD*eD/vD/vD + eB*eB/vB/vB); 
      meanR   += rat/ratE2; meanRE2 += 1.0/ratE2;
    }
  }
  //
  if (meanRE2>0) {
    meanR  /= meanRE2;
    meanRE2 = 1./meanRE2;
    meanRE = TMath::Sqrt(meanRE2);
  }
  if (meanDE2>0 && meanBE2>0) {
    meanRI  = meanD/meanB;
    meanRIE =  meanRI*TMath::Sqrt(meanDE2/meanD/meanD + meanBE2/meanB/meanB);
  }
  printf("NormalizeBg: Tails scaling %s wrt %s: Wgh.Mean:%.4f(%.4f) / Integral:%.4f(%.4f)\n",
	 bgH->GetName(),dataH->GetName(), meanR,meanRE, meanRI,meanRIE);
  printf("NormalizeBg: Select scaling type %s\n",useScaleType==kSclWghMean ? "Wgh.Mean":"Integral");
  //
  scl  = useScaleType==kSclWghMean ? meanR  : meanRI;
  sclE = useScaleType==kSclWghMean ? meanRE : meanRIE;
  //
  // rescaled bg
  char buff[1000];
  sprintf(buff,"%s_bgNorm",bgH->GetName());
  bgH = (TH1*)bgH->Clone(buff);
  sprintf(buff,"%s bgNorm%d %.4f+-%.4f",bgH->GetName(),useScaleType,scl,sclE);
  TH1* dumH = (TH1*)bgH->Clone("dummySCL$"); dumH->Reset();
  for (int i=1;i<=nbtot;i++) {
    dumH->SetBinContent(i,scl);
    dumH->SetBinError(i,sclE);
  }
  bgH->Multiply(dumH);
  delete dumH;
  return bgH;
}

//______________________________________________________________________
TObject* FindObject(const char* nameH, const TList* lst, Bool_t normPerEvent)
{
  // get histo, optionally normalizing it per processed event
  if (!lst) {printf("FindObject %s: No list provided\n",nameH); exit(1);}
  int nent = lst->GetEntries();
  TString nm;
  TObject *hst = 0;
  for (int i=nent;i--;) {
    nm = lst->At(i)->GetName();
    if (nm.EndsWith(nameH)) {hst = lst->At(i); break;}
  }
  if (!hst) {printf("FindObject: No %s histo in list %s\n",nameH,lst->GetName()); exit(1);}
  if (!normPerEvent || hst->TestBit(kBitNormPerEvent)) return hst; // already normalized
  TString nameHS = nameH;
  if (nameHS==kHStatName) return hst;                              // never scale stat. histo
  //
  TH1* hstat = (TH1*)FindObject(kHStatName,lst,kFALSE);
  double nrm = hstat->GetBinContent(kEvProcData);
  if (nrm<1) {printf("FindObject: Anomaluous %d number of events processed in list %p\n",int(nrm),lst); exit(1);}
  //
  // account for eventual cut in Z
  TH1* hzv = (TH1*)FindObject("zv",lst,kFALSE);
  int nbz = hzv->GetNbinsX();
  double zvTot = hzv->Integral(1,nbz);
  int zb0 = hzv->FindBin( useZvMin+kEps); if (zb0<1)   zb0 = 1;
  int zb1 = hzv->FindBin( useZvMax-kEps); if (zb1>nbz) zb1 = nbz;
  double zvSel = hzv->Integral(zb0,zb1);
  if (zvTot<1 || zvSel<1) {printf("No statistics: NzvTot: %.1f NzvSel:%.f\n",zvTot,zvSel); return 0;}
  else printf("%f fraction of selected events is used with current Zv cut %.1f:%.1f\n",zvSel/zvTot,useZvMin,useZvMax);
  nrm *= zvSel/zvTot;
  //
  if      (hst->InheritsFrom(TH1::Class())) ((TH1*)hst)->Scale(1./nrm);
  else if (hst->InheritsFrom(THnSparse::Class())) {
    THnSparse* spr = (THnSparse*) hst;
    spr->Sumw2();
    int coord[3] = {0,0,0};
    for (Long64_t i = 0; i < spr->GetNbins(); ++i) {
      // Get the content of the bin from the current histogram
      Double_t v = spr->GetBinContent(i, coord);
      spr->SetBinContent(coord, v/nrm);
      spr->SetBinError(coord,TMath::Sqrt(v)/nrm);
    }    
  }
  //
  hst->SetBit(kBitNormPerEvent);
  return hst;
}

//______________________________________________________________________
TList* LoadList(const char* flName, const char* addPref, const char* nameL)
{
  // load list with histos
  TString nms = flName;
  gSystem->ExpandPathName(nms);
  TFile* fl = TFile::Open(nms.Data());
  if (!fl) {printf("LoadList: No file %s\n",nms.Data()); exit(1);}
  TList* lst = (TList*)fl->Get(nameL);
  if (!lst) {printf("LoadList: No list %s in file %s\n",nameL,nms.Data()); exit(1);}
  lst->SetName(flName);
  //
  int nEnt = lst->GetSize();
  TString nm;
  for (int i=0;i<nEnt;i++) {
    TNamed* ob = (TNamed*)lst->At(i);
    nm = addPref;
    nm += ob->GetName();
    ob->SetName(nm.Data());
  }
  //
  return lst;
}

//____________________________________________________________________________
void GetRatE(double x,double xe, double y,double ye, double &rat, double &rate)
{
  rat = 0; rate = 0;
  if (TMath::Abs(y)<1e-16 || TMath::Abs(x)<1e-16) return;
  rat = x/y;
  rate = rat*TMath::Sqrt( xe*xe/(x*x) + ye*ye/(y*y));
}

//____________________________________________________________________________
void Integrate(TH1* hist, double xmn,double xmx, double &val, double& err)
{
  // integrate 1d histo within given limits
  TAxis* xax = hist->GetXaxis();
  int bmn = xax->FindBin(xmn+kEps); if (bmn<1) bmn = 0; // include 
  int bmx = xax->FindBin(xmx-kEps);
  val = hist->IntegralAndError(bmn,bmx,err);
  // is this histo with symmetric axis ? then integrate also negative half axis
  if (TMath::Abs( xax->GetXmin() + xax->GetXmax() )<1e-6) {
    bmn = xax->FindBin(-xmx+kEps); 
    bmx = xax->FindBin(-xmn-kEps); 
    double errn;
    val += hist->IntegralAndError(bmn,bmx,errn);
    err = TMath::Sqrt(err*err + errn*errn);
  }
}


//____________________________________________________________________________
const char* HName(const char* prefix,const char* htype)
{
  // compose the name of histo in the clist
  static TString strh;
  strh = "Tr"; strh += prefix; strh += "_"; strh += htype;
  return strh.Data();
}

//____________________________________________________________________________
Int_t CheckStat(const TList* lst, const char* dtType)
{
  // check if needed bg was generated
  TH1* hstat = (TH1*)FindObject(kHStatName,lst);
  TString dts = dtType;
  if (dts=="Data") return int( hstat->GetBinContent(kEvProcData) );
  if (dts=="Mix")  return int( hstat->GetBinContent(kEvProcMix) );
  if (dts=="Inj")  return int( hstat->GetBinContent(kEvProcInj) );
  if (dts=="Rot")  return int( hstat->GetBinContent(kEvProcRot) );
  printf("Unknown process %s statistics is checked. Alowed: Data,Mix,Inj,Rot",dtType);
  return 0;
}

//____________________________________________________________________________
void ProjDZE(THnSparseF* hrawd, THnSparseF* hgenb, THnSparseF* hmcComb, TObjArray* res, int bStepEta,int bStepZv)
{
  // project 3d histo of Dist vs Zv vs Eta to Zv vs Eta with all cuts
  int shift = hmcComb ? kMCShift : 0; // is this mc?
  //
  // determine boundaries for zv and eta cuts
  //
  double cutDstMin,cutDstMax;
  double cutBgMin,cutBgMax;
  double cutSgMin,cutSgMax;
  if (useShapeType==kNormShapeDist) {
    cutDstMin = 0;
    cutDstMax = kWDistBgTailMax;
    //
    cutBgMin = kWDistBgTailMin;
    cutBgMax = kWDistBgTailMax;
    //
    cutSgMin = 0;
    cutSgMax = kWDistSgCut;
  }
  else {
    cutDstMin = -kdPhiBgTailMax;
    cutDstMax =  kdPhiBgTailMax;    
    //
    cutBgMin =  kdPhiBgTailMin;
    cutBgMax =  kdPhiBgTailMax;
    //
    cutSgMin =  0;
    cutSgMax =  kdPhiSgCut;
  }
  //
  int bn0[3] = {0}; // 1st bin to count
  int bn1[3] = {0}; // last bin to count
  // eta 
  TAxis* axEta = hrawd->GetAxis(0);
  bn0[0] = axEta->FindBin(-useEtaCut+kEps);
  bn1[0] = axEta->FindBin( useEtaCut-kEps);
  // Zv
  TAxis* axZv = hrawd->GetAxis(1);
  bn0[1] = axZv->FindBin( useZvMin+kEps);
  bn1[1] = axZv->FindBin( useZvMax-kEps);
  // W.dist
  TAxis* axDst = hrawd->GetAxis(2);
  bn0[2] = axDst->FindBin( cutDstMin + kEps);
  bn1[2] = axDst->FindBin( cutDstMax - kEps);
  //
  //
  int nb[3] = { bn1[0]-bn0[0]+1, bn1[1]-bn0[1]+1, bn1[2]-bn0[2]+1};  // number of bins to count
  int nbTot[3] = {axEta->GetNbins(), axZv->GetNbins(), axDst->GetNbins() }; // total bins
  //
  if (bStepEta<1 || bStepEta>nb[0]) bStepEta = nb[0];
  if (bStepZv<1  || bStepEta>nb[1]) bStepZv  = nb[1];
  //
  for (int i=3;i--;) { // set default axis range
    hrawd->GetAxis(i)->SetRange(1,nbTot[i]);
    hgenb->GetAxis(i)->SetRange(1,nbTot[i]);
    if (hmcComb) hmcComb->GetAxis(i)->SetRange(1,nbTot[i]);
  }
  //
  char grpTit[100];
  sprintf(grpTit,"grp_eta%d_zv%d",bStepEta,bStepZv);
  //
  TString pref = hmcComb ? "mc" : "dt";
  //
  // "Data" histo with cut on tails where we look for signal
  hrawd->GetAxis(2)->SetRangeUser(cutDstMin+kEps,cutDstMax-kEps);
  TH2* hRawDtCut = hrawd->Projection(1,0,"e");
  hrawd->GetAxis(2)->SetRange(bn0[2],bn1[2]);
  hRawDtCut->SetName(pref+"_RawWithCut");
  hRawDtCut->SetTitle(pref+" Raw with cut on tracklets");
  res->AddAtAndExpand(hRawDtCut, kRawDtCut+shift);
  //
  // "Data - Est.Bg" histo with cut on tails where we look for signal
  hrawd->GetAxis(2)->SetRangeUser(cutDstMin+kEps,cutDstMax-kEps);
  TH2* hSignalEst = hrawd->Projection(1,0,"e");
  hrawd->GetAxis(2)->SetRange(bn0[2],bn1[2]);
  hSignalEst->SetName(pref+"_SignalCut_"+grpTit);
  hSignalEst->SetTitle(pref+" Signal (raw-bg) with cut on tracklets "+grpTit);
  res->AddAtAndExpand(hSignalEst, kSignalEst+shift);
  //
  // "Data - MC.Bg" histo with cut on tails where we look for signal
  TH2* hSignalEstMC = 0;
  if (hmcComb) {
    hrawd->GetAxis(2)->SetRangeUser(cutDstMin+kEps,cutDstMax-kEps);
    hSignalEstMC  = hrawd->Projection(1,0,"e");
    hrawd->GetAxis(2)->SetRange(bn0[2],bn1[2]);
    hSignalEstMC->SetName(pref+"_SignalCut_bgMCLabels_"+grpTit);
    hSignalEstMC->SetTitle(pref+" Signal (raw-bg_MCLabels) with cut on tracklets "+grpTit);
    res->AddAtAndExpand(hSignalEstMC, kSignalEstMC+shift);
  }
  //
  // Estimated background in the cut range
  hgenb->GetAxis(2)->SetRangeUser(cutDstMin+kEps,cutDstMax-kEps);
  TH2* hBgEst = hgenb->Projection(1,0,"e");
  hgenb->GetAxis(2)->SetRange(bn0[2],bn1[2]);
  hBgEst->SetName(pref+"_BgEst_"+grpTit);
  hBgEst->SetTitle(pref+" Estimated Bg "+grpTit);
  res->AddAtAndExpand(hBgEst, kBgEst+shift);
  //
  // 1-beta for "data" = (Data_cut - Bg_cut) / Data_cut
  TH2* h1mBeta     = hrawd->Projection(1,0,"e"); 
  h1mBeta->Reset();
  h1mBeta->SetName(pref+"_h1mBeta_"+grpTit);
  h1mBeta->SetTitle(pref+" 1-#beta with gen.bg. "+grpTit);
  res->AddAtAndExpand(h1mBeta, k1MBeta+shift);
  //
  // If MC labels info is provided
  TH2* h1mBetaMC = 0;  // 1-beta for MC with bg from labels
  TH2* hBgMC = 0;      // bg from MC labels
  if (hmcComb) {
    hmcComb->GetAxis(2)->SetRangeUser(cutDstMin+kEps,cutDstMax-kEps);
    h1mBetaMC = hmcComb->Projection(1,0,"e");
    h1mBetaMC->SetName(pref+"_h1mBetaMC");
    h1mBetaMC->SetTitle(pref+" 1-#beta with bg. from MC labels");
    h1mBetaMC->Divide(hRawDtCut);
    res->AddAtAndExpand(h1mBetaMC, k1MBetaMC+shift);
    for (int ib0=1;ib0<=nbTot[0];ib0++) 
      for (int ib1=1;ib1<=nbTot[1];ib1++) 
	h1mBetaMC->SetBinContent(ib0,ib1, 1.- h1mBetaMC->GetBinContent(ib0,ib1)); 
    //
    hBgMC = hmcComb->Projection(1,0,"e");
    hBgMC->SetName(pref+"_Bg_MClab");
    hBgMC->SetTitle(pref+" Bg from MC labels");
    res->AddAtAndExpand(hBgMC, kBgMC+shift);
    //
    // finalize estimated signal with bg from MC labels
    hSignalEstMC->Add(hBgMC,-1);
    //
    hmcComb->GetAxis(2)->SetRange(bn0[2],bn1[2]);
  }
  // 
  // rescaling factors for generated bg 
  TH2* hBgRescFc  = hrawd->Projection(1,0,"e"); hBgRescFc->Reset();
  hBgRescFc->SetName(pref+"_hBgRescFactors_"+grpTit);
  hBgRescFc->SetTitle(pref+" Scale.factor for gen.bg. "+grpTit);
  res->AddAtAndExpand(hBgRescFc, kBgRescFc+shift);
  //
  int nbint = bStepEta*bStepZv;
  float nbinstsq = TMath::Sqrt(nbint);
  hrawd->GetAxis(0)->SetRange(bn0[0],bn1[0]);
  hrawd->GetAxis(1)->SetRange(bn0[1],bn1[1]);
  TH1* hDstDt = hrawd->Projection(2,"e");
  hrawd->GetAxis(0)->SetRange(1,nbTot[0]);
  hrawd->GetAxis(1)->SetRange(1,nbTot[1]);
  hDstDt->SetName(pref+"_DistRaw_"+grpTit);   // "Data" projection on tracklet quality (w.distance) axis to check the tails
  hDstDt->SetTitle(pref+" DistRaw "+grpTit); 
  double nrmDst,dumErr = 0;
  Integrate(hDstDt, cutBgMin,cutBgMax, nrmDst, dumErr);
  hDstDt->Scale(1./nrmDst);
  res->AddAtAndExpand(hDstDt, kDataDist+shift);
  //
  hgenb->GetAxis(0)->SetRange(bn0[0],bn1[0]);
  hgenb->GetAxis(1)->SetRange(bn0[1],bn1[1]);
  TH1* hDstBg = hgenb->Projection(2,"e");  hDstBg->Reset();  
  hgenb->GetAxis(0)->SetRange(1,nbTot[0]);
  hgenb->GetAxis(1)->SetRange(1,nbTot[1]);
  hDstBg->SetName(pref+"_DistBgNorm_"+grpTit);  // Gen.Bg projection on tracklet quality (w.distance) axis to check the tails
  hDstBg->SetTitle(pref+" DistBgNorm "+grpTit);  // Gen.Bg projection on tracklet quality (w.distance) axis to check the tails
  res->AddAtAndExpand(hDstBg, kBgDist+shift);
  TH1* hDstBgMC = 0;
  if (hmcComb) {
    hmcComb->GetAxis(0)->SetRange(bn0[0],bn1[0]);
    hmcComb->GetAxis(1)->SetRange(bn0[1],bn1[1]);
    hDstBgMC = hmcComb->Projection(2,"e");
    hmcComb->GetAxis(0)->SetRange(1,nbTot[0]);
    hmcComb->GetAxis(1)->SetRange(1,nbTot[1]);
    hDstBgMC->SetName(pref+"_DistBgMC");
    hDstBgMC->SetTitle(pref+" Bg. Distance from MC labels");
    hDstBgMC->Scale(1./nrmDst);
    res->AddAtAndExpand(hDstBgMC, kBgMCDist+shift);
  }
  //
  // fill 1-beta matrix
  for (int ib0=bn0[0];ib0<=bn1[0];ib0+=bStepEta) { // eta
    hrawd->GetAxis(0)->SetRange(ib0,ib0+bStepEta-1);
    hgenb->GetAxis(0)->SetRange(ib0,ib0+bStepEta-1);
    for (int ib1=bn0[1];ib1<=bn1[1];ib1+=bStepZv) { // zv
      hrawd->GetAxis(1)->SetRange(ib1,ib1+bStepZv-1);
      hgenb->GetAxis(1)->SetRange(ib1,ib1+bStepZv-1);
      //
      TH1D* dstD = hrawd->Projection(2,"e"); // data "qaulity" for given eta:zv bin
      TH1D* dstB = hgenb->Projection(2,"e"); // data "qaulity" for given eta:zv bin
      double scl,sclE;
      TH1* nrmB = NormalizeBg(dstD,dstB,scl,sclE);  // get rescaling factor for bg. from tails comparison
      double bgVal,bgErr;
      double dtVal,dtErr;
      // integral in the range where we look for signal
      Integrate(nrmB, cutSgMin, cutSgMax, bgVal, bgErr);
      Integrate(dstD, cutSgMin, cutSgMax, dtVal, dtErr);
      double beta,betaErr;
      GetRatE(bgVal,bgErr, dtVal, dtErr,beta,betaErr);
      //    betaErr*=nbinstsq; // ??? RS
      for (int i=ib0;i<ib0+bStepEta;i++) {
	for (int j=ib1;j<ib1+bStepZv;j++) {
	  hBgRescFc->SetBinContent(i,j, scl);
	  hBgRescFc->SetBinError(i,j, sclE);
	}
      }
      hDstBg->Add(nrmB);
      delete dstD;
      delete dstB;
      delete nrmB;
      //
    }
  }
  hDstBg->Scale(1./nrmDst);
  //
  // finalize estimated bg and signal matrices
  hBgEst->Multiply(hBgRescFc);
  hSignalEst->Add(hBgEst,-1);
  //
  // finalize 1-beta
  for (int ib0=bn0[0];ib0<=bn1[0];ib0++) { // eta
    for (int ib1=bn0[1];ib1<=bn1[1];ib1++) { // zv
      //      printf("Bin %d %d\n",ib0,ib1);
      double bg  = hBgEst->GetBinContent(ib0,ib1);
      double bgE = hBgEst->GetBinError(ib0,ib1);
      double dt  = hRawDtCut->GetBinContent(ib0,ib1);
      double dtE = hRawDtCut->GetBinError(ib0,ib1);
      double beta,betaE;
      GetRatE(bg,bgE,dt,dtE, beta,betaE );      
      h1mBeta->SetBinContent(ib0,ib1,1.-beta);
      h1mBeta->SetBinError(ib0,ib1,betaE);
      //
    }
  }
  //
  // crop to needed range
  for (int i=shift;i<kNHistos+shift;i++) {
    TH2* hist = (TH2*)res->At(i);
    if (!hist || !hist->InheritsFrom(TH2::Class())) continue;
    CropHisto(hist, bn0[0],bn1[0], bn0[1],bn1[1]);
  }
  //
  h1mBeta->SetMinimum(0.6);
  h1mBeta->SetMaximum(0.85);
  if (hmcComb) {
    h1mBetaMC->SetMinimum(0.6);
    h1mBetaMC->SetMaximum(0.85);
  }
  //
  // restore
  for (int i=3;i--;) { // set default axis range
    hrawd->GetAxis(i)->SetRange(1,nbTot[i]);
    hgenb->GetAxis(i)->SetRange(1,nbTot[i]);
    if (hmcComb) hmcComb->GetAxis(i)->SetRange(1,nbTot[i]);
  }
  //
}


void GetRealMinMax(TH1* histo, double &vmn, double &vmx)
{
  TAxis *xax = histo->GetXaxis();
  int nbx = xax->GetNbins(); 
  vmn=1e6, vmx=-1e6;
  if (histo->InheritsFrom(TH2::Class())) {
    TAxis *yax = histo->GetYaxis();
    int nby = yax->GetNbins();
    for (int ix=nbx+2;ix--;) {
      for (int iy=nby+2;iy--;) {
	double vl = histo->GetBinContent(ix,iy);
	if (vl<kEps) continue;
	if (vl<vmn) vmn = vl;
	if (vl>vmx) vmx = vl;
      }
    }
  }
  //
  else {
    for (int ix=nbx+2;ix--;) {
      double vl = histo->GetBinContent(ix);
      if (vl<vmn) vmn = vl;
      if (vl>vmx) vmx = vl;
    }
  }
  //
}
