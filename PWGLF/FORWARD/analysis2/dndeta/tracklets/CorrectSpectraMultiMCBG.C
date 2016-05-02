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
// #include "TPaletteAxis.h"
#include "TArrayD.h"
#include "TGraphErrors.h"
//
//
#endif
#include "SaveCanvas.h"



const char kHStatName[]="hStat";
double kEps = 1e-6;
//
double kdPhiBgTailMin = 0.1; // lower limit of dphi tail to use for bg normalization  
double kdPhiBgTailMax = 0.3; // upper limit of dphi tail to use for bg normalization   
//
double kWDistBgTailMin = 5.; // lower limit of wgh.distance tail to use for bg normalization  
double kWDistBgTailMax = 25.; // upper limit of wgh.distance tail to use for bg normalization

double kdPhiSgCut=-1;       // cut in dphi-bent used to extract the signal, extracted from stat histo
double kWDistSgCut=-1;      // cut in w.distance used to extract the signal, extracted from stat histo
//
enum { kNormShapeDist,      // normalize bg tails usig weighted distance shape 
       kNormShapeDPhi,      // normalize bg tails usig dPhi-bend shape
       kNormShapes};

enum { kSclWghMean,         // normalize bg tails to data using weighted mean of bin-by-bin ratios 
       kSclIntegral,        // normalize bg tails to data using integral
       kSclTypes};


const char* figDir = "corrFig";//"figMult";
const char* resDir = "corrRes";//"resMult";
Bool_t   use1mBeta    = true;
TString  useBgType    = "Comb"; // "Inj"; // "Inj"; "Comb";
Int_t    useShapeType = kNormShapeDist;    // which distribution to use for bg normalization
// If set 1 one, then 
Int_t    useMCLB      = 2;             // use Comb MC Labels as a template for Bg.
Double_t scaleBG      = 1.3;    // apply scaling to extracted bg for data
Int_t    useScaleType = kSclIntegral;//kSclWghMean;       // which type of tails normalization to use
Bool_t   useZbinWAv  = kFALSE; // kFALSE;
//Bool_t   normToMB =  kFALSE;    // normalize with MB MC truth, rather than for selected events
Bool_t   normToMB =  kTRUE;    // normalize with MB MC truth, rather than for selected events
const double kEtaFitRange = 0.5;

Double_t minAlpha = 0.;
Double_t maxAlpha = 2.5;

enum {kBitNormPerEvent=1<<14};
// bins for saved parameters in the hStat histo
enum {kDummyBin,
      kEvTot0,      // events read
      kEvTot,       // events read after vertex quality selection
      kOneUnit,     // just 1 to track primate merges
      kNWorkers,    // n workers
      //
      kCentVar,     // cetrality var. used
      kDPhi,        // dphi window
      kDTht,        // dtheta window
      kNStd,        // N.standard deviations to keep
      kPhiShift,    // bending shift
      kThtS2,       // is dtheta scaled by 1/sin^2
      kThtCW,       // on top of w.dist cut cut also on 1 sigma dThetaX
      kPhiOvl,      // overlap params
      kZEtaOvl,     // overlap params
      kNoOvl,       // flag that overlap are suppressed
      //
      kPhiRot,      // rotation phi
      kInjScl,      // injection scaling
      kEtaMin,      // eta cut
      kEtaMax,      // eta cut
      kZVMin,       // min ZVertex to process
      kZVMax,       // max ZVertex to process
      //
      kDPiSCut,     // cut on dphi used to extract signal (when WDist is used in analysis, put it equal to kDPhi
      kNStdCut,     // cut on weighted distance (~1) used to extract signal 
      //
      kMCV0Scale,   // scaling value for V0 in MC
      //
      // here we put entries for each mult.bin
      kBinEntries = 50,
      kEvInMltBin = 0,
      kEvProcData,  // events with data mult.object (ESD or reco) passing all selections
      kEvProcInj,   // events Injected, total
      kEvProcRot,   // events Rotated
      kEvProcMix,   // events Mixed
      kEntriesPerBin = 10
};


enum {kSigCorr,kMCPrim,kRawDtCut,kSignalEst,kSignalEstMC,kBgEst,k1MBeta,k1MBetaMC,k1MBetaMCscl,kAlpha,kAlphaMC,kBgMC,kBgRescFc,
      kDataDist,kBgDist,kBgMCDist,kZvDist,kZvMCDistNS,kZvDistCorr,kZvEff,kMCShift=20, kNHistos=kMCShift+kMCShift};

#ifndef __CINT__
#include <TROOT.h>
UShort_t fgDebug = 1;
void _MyPrint(UShort_t lvl, const char* fmt, ...)
{
  if (lvl > fgDebug) return;
  char buf[512];
  va_list ap;
  va_start(ap, fmt);
  vsnprintf(buf, 511, fmt, ap);
  buf[511] = '\0';
  va_end(ap);
  gROOT->IndentLevel();
  Printf("%s",buf);
}
void Incr() { gROOT->IncreaseDirLevel(); }
void Decr() { gROOT->DecreaseDirLevel(); }
struct _MyGuard
{
  _MyGuard(UShort_t lvl, const char* fmt, ...)
    : ok(lvl <= fgDebug)
  {
    if (!ok) return;
    char buf[512];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, 511, fmt, ap);
    buf[511] = '\0';
    va_end(ap);
    gROOT->IndentLevel();
    Printf("%s",buf);
    Incr();
  }
  ~_MyGuard() { if (ok) Decr(); }
  Bool_t ok;
};
#define MyPrint(L,F,...) _MyPrint(L,F, ## __VA_ARGS__)
#define MyGuard(L,F,...) _MyGuard _guard(L,F, ## __VA_ARGS__)
#endif

void PrintAndPause(TCanvas* c, const TString& what, Bool_t wait)
{
  Printf("Wat is %d", wait);
  c->Modified();
  c->Update();
  c->cd();
  c->Print(what);
  if (wait) c->WaitPrimitive();
}

void    CorrectSpectraMultiMCBG(const char* flNameData, const char* flNameMC,const char* unique="",int maxBins=10, Bool_t waitForUser=false, const char* bgType="Comb");
Bool_t  PrepareHistos(int bin, TList* lst, TList* lisMC, Bool_t isMC);
void    ProcessHistos(int bin);
TH1*    NormalizeBg(TH1* dataH, TH1* bgH, double &scl, double &scle);
TObject* FindObject(int bin, const char* nameH, const TList* lst, Bool_t normPerEvent=kTRUE);
TList*  LoadList(const char* flName, const char* addPref, const char* nameL="clist");
void    GetRatE(double x,double xe, double y,double ye, double &rat, double &rate);
Int_t   CheckStat(const TList* lst,const char* dtType);
void    Integrate(TH1* hist, double xmn,double xmx, double &val, double& err);
void    CropHisto(TH1* histo, int b00, int b01, int b10=-1, int b11=-1);
void    CropHisto(TH1* histo, double b00, double b01, double b10=-1, double b11=-1);
void    GetRealMinMax(TH1* h, double &vmn, double &vmx);
const char*   HName(const char* prefix,const char* htype);
TH1*    ProjNorm(TH2* hEtaZ, TH1* hZv,const char* name="_px", Int_t firstbin=0, Int_t lastbin=-1);
TH1*    ProjectWghMean(TH2* hEtaZ, const char* name = "_px", Int_t firstbin = 0, Int_t lastbin = -1, double rejOutliers=6.);
void    CorrectForZV(TH2* hEtaZ, TH1* hZv);
void    KillBadBins(TH2* histo, double mn=-1e50,double mx=1e50);

void PlotResults(Bool_t waitForUser);
void PlotDNDEta(int bin);
void PlotAlphaBeta(int bin);
void PlotSpecies();
void PrintH(TH2* h, Int_t prec=2);
void PrintH(TH1* h, Int_t prec=2);

Int_t nCentBins = -1;
TList *listDt=0, *listMC=0;
TObjArray resArr, resDnDeta;
char outStr[1000];
char outTitle[1000];
TString uniqueName="";
//
TArrayD dNdEta,dNdEtaErr;
TCanvas *canvFin=0;
//
Bool_t creatDnDEtaCMacro = kTRUE;
Bool_t creatAlphaBetaCMacro = kTRUE;
Bool_t creatSpeciesCMacro = kTRUE;

void CorrectSpectraMultiMCBG(const char* flNameData, const char* flNameMC, const char* uniqueNm,int maxBin,Bool_t waitForUser,const char* bgType) 
{
  useBgType = bgType;
  //
  uniqueName = uniqueNm;
  Printf("WaitForUser is %d", waitForUser);
  MyPrint(1,"Loading lists from %s and %s", flNameData, flNameMC);
  {
    MyGuard(1,"Loading data lists");
    listDt = LoadList(flNameData,"dt_");
  }
  {
    MyGuard(1,"loading MC lists");
    listMC = LoadList(flNameMC,"mc_");
  }
  //
  resArr.Clear();
  //
  TH1* hstat = (TH1*)FindObject(-1,kHStatName,listDt,kFALSE);
  MyPrint(2,"Got statistics histogram %p", hstat);
  
  //TH1* hstat = (TH1*)FindObject(-1,kHStatName,listMC,kFALSE);
  //
  int nbstat = hstat->GetNbinsX();
  nCentBins = (nbstat - kBinEntries)/kEntriesPerBin;
  nCentBins = nCentBins>maxBin ? maxBin : nCentBins;
  MyPrint(2,"Will process %d centrality bins", nCentBins);

  if (nCentBins<1) return;
  double myMergeFactor = hstat->GetBinContent(kOneUnit);
  if (myMergeFactor<1) myMergeFactor = 1.;
  MyPrint(2,"Output was merged from %5.1f jobs", myMergeFactor);

  //
  dNdEta.Set(nCentBins);
  dNdEtaErr.Set(nCentBins);
  //
  kdPhiSgCut  = hstat->GetBinContent(kDPiSCut)/myMergeFactor;
  kWDistSgCut = hstat->GetBinContent(kNStdCut)/myMergeFactor;
  MyPrint(2,"Signal cuts used: dPhiS: %f WDist:%f",
	  kdPhiSgCut,kWDistSgCut);

  {
    MyGuard(1,"Now looping over bins");
    for (int ib=0;ib<nCentBins;ib++) {
      {
	MyGuard(1,"Extracting from MC");
	if (!PrepareHistos(ib,listMC,listMC,kTRUE))  continue;//return;
      }
      {
	MyGuard(1,"Extracting from Data");
	if (!PrepareHistos(ib,listDt,listMC,kFALSE)) continue;//return;
      }
      {
	MyGuard(1,"Processing information");
	ProcessHistos(ib);
      }
    }
  }

  //
  sprintf(outStr,"CutEta%.1f_%.1f_Zv%.1f_%.1f_bg%s_Shape_%s_%s_cutSig%.1f_cutBg%.1f",
	  hstat->GetBinContent(kEtaMin)/myMergeFactor,
	  hstat->GetBinContent(kEtaMax)/myMergeFactor,
	  hstat->GetBinContent(kZVMin)/myMergeFactor,
	  hstat->GetBinContent(kZVMax)/myMergeFactor,	 
	  useBgType.Data(),
	  useShapeType==kNormShapeDist ? "wdst":"dphi",
	  (use1mBeta ? "MCBG" : "DTBG"), 
	  useShapeType==kNormShapeDist ? kWDistSgCut:kdPhiSgCut,
	  useShapeType==kNormShapeDist ? kWDistBgTailMin:kdPhiBgTailMin);
  //
  PlotResults(waitForUser);
  //
  printf("Final Results:\n");
  printf("dNdEta:    "); 
  for (int i=nCentBins;i--;) printf("%.2f,",dNdEta[i]); printf("\n");
  printf("dNdEtaErr: "); 
  for (int i=nCentBins;i--;) printf("%.2f,",dNdEtaErr[i]); printf("\n");
  //
  TString rtnm1 = resDir; rtnm1 += "/"; rtnm1 += uniqueName; 
  rtnm1 += "_"; rtnm1 += nCentBins; rtnm1+= "bins_";
  rtnm1 += outStr; rtnm1 += ".root";
  TFile* flRes = TFile::Open(rtnm1.Data(),"recreate");
  flRes->WriteObject(&resDnDeta, "TObjArray", "kSingleKey");
  flRes->WriteObject(&resArr, "TObjArrayAux", "kSingleKey");
  flRes->Close();
  delete flRes;
  printf("Stored result in %s\n",rtnm1.Data());
  //
}

//_____________________________________________________________________
TH1* CopyAdd(TH1* h,
	     const char*     name,
	     const char*     title,
	     TObjArray*      col,
	     Int_t           location,
	     Int_t           shift)
{
  TH1* ret = static_cast<TH1*>(h->Clone(name));
  ret->SetTitle(title);
  if (col) col->AddAtAndExpand(ret, location+shift);
  return ret;
}

//_____________________________________________________________________
Bool_t PrepareHistos(int bin, TList* lst, TList* lstMC, Bool_t isMC)
{
  if (!lst) {
    Warning("PrepareHistos", "No list for %s", isMC ? "simulation" : "data");
    return false;
  }
  // fill standard histos for given bin
  //
  char buffn[500];
  char bufft[500];
  //
  double cutBgMin,cutBgMax;
  double cutSgMin,cutSgMax;
  //
  if (useShapeType==kNormShapeDist) {
    cutBgMin = kWDistBgTailMin;
    cutBgMax = kWDistBgTailMax;
    //
    cutSgMin = 0;
    cutSgMax = kWDistSgCut;
  }
  else {
    cutBgMin =  kdPhiBgTailMin;
    cutBgMax =  kdPhiBgTailMax;
    //
    cutSgMin =  0;
    cutSgMax =  kdPhiSgCut;
  }
  //
  TH1* hstat = (TH1*)FindObject(-1,kHStatName,lst,kFALSE);
  MyPrint(2,"Got statistics histogram %p", hstat);
  
  double nrmBin = hstat->GetBinContent(kBinEntries + kEvProcData+bin*kEntriesPerBin);
  if (nrmBin<1) {
    Warning("PrepareHistos", "Norm is 0 for bin %d",bin);
    return kFALSE;
  }
  MyPrint(0,"Normalization for bin %d: %f", bin, nrmBin);

  const char* zeCut = "ZvEtaCutT";
  TObjArray* res = &resArr;
  int shift = bin*kNHistos + (isMC ? kMCShift : 0);
  //
  // histo for "data" Z vs Eta with default signal cut - scaled by IPz dist
  TH2* tmp2 = (TH2*)FindObject(bin,HName("Data",zeCut),lst);
  if (!tmp2) {
    Warning("PrepareHistos", "Didn't find %s", HName("Data",zeCut));
    return kFALSE;
  }
  TH2* hRawDtCut = (TH2*)CopyAdd(tmp2,
				 Form("bin%d_%s_RawWithCut",bin,isMC?"mc":"dt"),
				 Form("bin%d %s Raw Data with cut on tracklets",
				      bin,isMC ? "mc":"dt"),
				 res,kRawDtCut,shift);
  MyPrint(1,"Measurements %s (clone of %s)", hRawDtCut->GetName(),
	  tmp2->GetName());
  
  // Number of eta, IPz bins 
  int nbEta = hRawDtCut->GetXaxis()->GetNbins();
  int nbZV  = hRawDtCut->GetYaxis()->GetNbins();
  //
  // Zv nonuniformity histo - IPz vs Centrality 
  TH2* hzv2 = (TH2*) FindObject(-1,"zv",lst, kFALSE);
  TH1* hz   = hzv2->ProjectionX(Form("zv%d_%s",bin,isMC ? "mc":"dt"),
				bin+1,bin+1,"e");
  MyPrint(1,"Get histogram of vertex disp for bin %d and scale by %f",
	  bin, nrmBin);
  hz->Scale(1./nrmBin);
  res->AddAtAndExpand(hz, kZvDist+shift);
  
  //
  // "Data - Est.Bg" histo with cut on tails where we look for signal
  TH2* hSignalEst = 
    (TH2*)CopyAdd(hRawDtCut,
		  Form("bin%d_%s_SignalWithCut",bin,isMC ? "mc":"dt"),
		  Form("bin%d %s Signal (raw-bg) with cut on tracklets",
		       bin,isMC ? "mc":"dt"),
		  res,kSignalEst,shift);
  MyPrint(1,"Signal estimate: %s (clone of %s)", hSignalEst->GetName(),
	  hRawDtCut->GetName());
  
  //
  // "Data - MC.Bg" histo with cut on tails where we look for signal
  TH2* hSignalEstMC = 0;
  if (isMC) {
    MyGuard(1, "MC specific code - MC signal estimate");
    hSignalEstMC =
      (TH2*)CopyAdd(hRawDtCut,
		    Form("bin%d_%s_SignalWithCut_bgMCLabels",
			 bin,isMC ? "mc":"dt"),
		    Form("bin%d %s Signal (raw-bg_MCLabels) w/cut on tracklets",
			 bin,isMC ? "mc":"dt"),
	      res,kSignalEstMC,shift);
    MyPrint(1,"MC signal estimate %s (copy of %s)",hSignalEstMC->GetName(),
	    hRawDtCut->GetName());
    // also currently a copy of eta vs IPz
  }
  //

  // Estimated background in the cut range
  TList* localLst = (useBgType.EqualTo("Comb") ? lstMC : lst);
  tmp2        = (TH2*)FindObject(bin,HName(useBgType.Data(),zeCut),localLst);
  TH2* hBgEst = (TH2*)CopyAdd(tmp2,
			      Form("bin%d_%s_BgEst",bin,isMC ? "mc":"dt"),
			      Form("bin%d %s Estimated Bg",bin,isMC?"mc":"dt"),
			      res,kBgEst,shift);
  MyPrint(1,"Bg estimate %s (clone of %s)",hBgEst->GetName(),tmp2->GetName());
  
  TH2* h1mBeta = (TH2*)CopyAdd(hBgEst,
			       Form("bin%d_%s_1mBeta",bin,isMC ? "mc":"dt"),
			       Form("bin%d %s 1-#beta with estimated bg",
				    bin,isMC ? "mc":"dt"),
			       res, k1MBeta, shift);
  h1mBeta->Reset();
  MyPrint(1, "1-beta %s (clone of %s)", h1mBeta->GetName(), hBgEst->GetName());

  //
  // If MC labels info is provided, prepare 1-beta with MC bg
  TH2* h1mBetaMC = 0;  // 1-beta for MC with bg from labels
  TH2* h1mBetaMCscl = 0; // scaled by ad hoc factor
  TH2* hBgMC = 0;      // bg from MC labels
  if (isMC) {
    MyGuard(1, "MC specific code - prepare 1-beta");
    tmp2 = (TH2*)FindObject(bin,HName("Comb",zeCut),lstMC);
    MyPrint(1,"Get MC background from %s", HName("Comb",zeCut));
    if (!tmp2) return kFALSE;
    hBgMC = (TH2*)CopyAdd(tmp2,
			  Form("bin%d_%s_BgMC",bin,isMC ? "mc":"dt"),
			  Form("bin%d %s Bg from MC labels",
			       bin,isMC ? "mc":"dt"),
			  res,kBgMC,shift);
    MyPrint(1,"Set MC background %s (clone of %s)",hBgMC->GetName(),
	    tmp2->GetName()); 
    
    //    printf(">>Here1\n");
    // Now scale by measured eta vs IPz 
    h1mBetaMC = (TH2*)CopyAdd(hBgMC,
			      Form("bin%d_%s_h1mBetaMC",bin,isMC ? "mc":"dt"),
			      Form("bin%d %s 1-#beta with bg. from MC labels",
				   bin,isMC ? "mc":"dt"),
			      res, k1MBetaMC,shift);
    h1mBetaMC->Divide(hRawDtCut);
    MyPrint(1, "1-beta (MC) %s (ratio of %s to %s)",
	    h1mBetaMC->GetName(),
	    hBgMC->GetName(), hRawDtCut->GetName());

    h1mBetaMCscl = (TH2*)CopyAdd(h1mBetaMC,
				 Form("%s_scl",h1mBetaMC->GetName()),
				 h1mBetaMC->GetTitle(),
				 res,k1MBetaMCscl,shift);
    h1mBetaMCscl->Scale(scaleBG);
    MyPrint(1, "1-beta (MC) %s (clone of %s)", h1mBetaMCscl->GetName(),
	    h1mBetaMC->GetName());
    //    printf("<<Here1\n");
    for (int ib0=1;ib0<=nbEta;ib0++) 
      for (int ib1=1;ib1<=nbZV;ib1++) {
	// Do the "1-" part
	Double_t oneMBeta    = h1mBetaMC->GetBinContent(ib0,ib1);
	Double_t oneMBetaScl = h1mBetaMCscl->GetBinContent(ib0,ib1);
	h1mBetaMCscl->SetBinContent(ib0,ib1, 1.- oneMBetaScl);
	h1mBetaMC->SetBinContent(ib0,ib1, 1.- oneMBeta);
      }
    MyPrint(1,"Did 1-X on %s and %s (%f)", h1mBetaMC->GetName(),
	    h1mBetaMCscl->GetName(), scaleBG);
    //
    // Now we have (1-beta)*meassured=no_combinatorics_signal 
    hSignalEstMC->Multiply(h1mBetaMC);
    MyPrint(1,"Multiply %s onto %s", h1mBetaMC->GetName(),
	    hSignalEstMC->GetName());
  }
  //
  // uncut w.distance or dphi distribution for data
  TH1* tmp1 = (TH1*)FindObject(bin,HName("Data",useShapeType==kNormShapeDist?
					   "WDist":"DPhiS"),lst);
  if (!tmp1) return kFALSE;
  TH1* hDstDt = CopyAdd(tmp1,
			Form("bin%d_%s_DistRawData",bin,isMC ? "mc":"dt"),
			Form("bin%d %s Raw Distance for Data",
			     bin,isMC?"mc":"dt"),
			res,kDataDist,shift);
  MyPrint(1,"Delta dist %s (clone of %s)", hDstDt->GetName(), tmp1->GetName());

  // Integrate tail of measured Delta dist 
  double nrmDst,dumErr = 0;
  Integrate(hDstDt, cutBgMin,cutBgMax, nrmDst, dumErr);
  MyPrint(1,"Integral of %s is %f +/- %f", hDstDt->GetName(), nrmDst, dumErr);
  
  if (nrmDst<1e-10) nrmDst = 1.; // RS resque
  // Scale measured Delta distribution by integral
  hDstDt->Scale(1./nrmDst);
  MyPrint(1,"Scaled %s and by integral %f", hDstDt->GetName(),nrmDst);

  //
  // uncut w.distance or dphi distribution for generated bg
  TH1* hDstBg = 0;
  if (!useBgType.IsNull()) {
    tmp1 = (TH1*)FindObject(bin,HName(useBgType.Data(),
				      useShapeType==kNormShapeDist ?
				      "WDist":"DPhiS"),localLst);
    if (!tmp1) return false;
    hDstBg = CopyAdd(tmp1,
		     Form("bin%d_%s_DistRawGenBgNorm",bin,isMC ? "mc":"dt"),
		     Form("bin%d %s Raw Dist. for Gen.Bg. Normalized to data",
			  bin,isMC ? "mc":"dt"),
		     res, kBgDist, shift);
    // Scale bg Delta to mesaured Delta integral 
    hDstBg->Scale(1./nrmDst);
    MyPrint(1,"%s (copy of %s) scaled by %f", hDstBg->GetName(),
	    tmp1->GetName(), nrmDst);
  }
  //
  // fill 1-beta matrix
  double scl,sclE;
  // get rescaling factor for bg. from tails comparison
  hDstBg = NormalizeBg(hDstDt,hDstBg,scl,sclE);
  MyPrint(1, "Scaling between %s and %s: %f +/- %f",
	  hDstDt->GetName(), hDstBg->GetName(), scl, sclE);
  res->AddAtAndExpand(hDstBg, kBgDist+shift);
  double bgVal,bgErr;
  double dtVal,dtErr;
  // integral in the range where we look for signal
  Integrate(hDstBg, cutSgMin, cutSgMax, bgVal, bgErr);
  Integrate(hDstDt, cutSgMin, cutSgMax, dtVal, dtErr);

  double sclb,sclbErr;
  GetRatE(bgVal,bgErr, dtVal, dtErr,sclb,sclbErr);
  MyPrint(1,"Ratio of (%s) %f +/- %f over (%s) %f +/- %f -> %f +/- %f",
	  hDstBg->GetName(), bgVal, bgErr,
	  hDstDt->GetName(), dtVal, dtErr,
	  sclb, sclbErr);
  
  //  hDstBg->Scale(1./nrmDst);
  //
  // finalize estimated bg and signal matrices
  hBgEst->Scale(scl);
  MyPrint(1, "Scaling %s by %f", hBgEst->GetName(), scl);
  //
  if (use1mBeta) {
    MyGuard(1,"Using MC Combinatorials fraction for real data bg.estimate");
    TH2* bet1m = (TH2*)res->At((isMC ? k1MBetaMC : k1MBetaMCscl)
			       + shift + (isMC ? 0 :kMCShift) );
    MyPrint(1,"Multiply %s with %s", hSignalEst->GetName(), bet1m->GetName());
    hSignalEst->Multiply(bet1m);
    // if (isMC) PrintH(hSignalEst);
  }
  else {
    MyPrint(1,"Subtracting %s from %s",hBgEst->GetName(),hSignalEst->GetName());
    hSignalEst->Add(hBgEst,-1);
  }
  //
  // finalize 1-beta
  {
    MyGuard(1,"Finalize 1-beta");
    for (int ib0=1;ib0<=nbEta;ib0++) { // eta
      for (int ib1=1;ib1<=nbZV;ib1++) { // zv
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
    MyPrint(1, "%s calculated as 1 minus ratio of %s to %s",
	    h1mBeta->GetName(), hBgEst->GetName(), hRawDtCut->GetName());
  }
  //
  if (isMC) {
    MyGuard(1,"MC specific stuff");
    // Zv nonuniformity histo: MCGen before cuts
    TH2* hzv2ns = (TH2*) FindObject(-1,"zvMCNoPS",lst, kFALSE);
    MyPrint(1, "MC truth vertex distribution %s", hzv2ns->GetName());
    TH1* hzns   = hzv2ns->ProjectionX(Form("zvMCNS%d",bin),bin+1,bin+1,"e");
    MyPrint(1, "Projection %s", hzns->GetName());
    //
    double myMergeFactor = hstat->GetBinContent(kOneUnit);
    double zmn = hstat->GetBinContent(kZVMin)/myMergeFactor;
    double zmx = hstat->GetBinContent(kZVMax)/myMergeFactor;
    int zbmn = hzns->FindBin(zmn+1e-3);
    int zbmx = hzns->FindBin(zmx-1e-3);
    double nEvTot = hzns->Integral(zbmn,zbmx);
    MyPrint(1,"Scale MC Truth By %9.1f",nEvTot);
    if (nEvTot<1) exit(1);
    hzns->Scale(1./nEvTot);
    /*
    // here there could be some over/under flows, add them to
    // first/last used bins correspondingly
    // 
    double under = hzns->Integral(0,zbmn-1);
    double over  = hzns->Integral(zbmx+1,hzns->GetNbinsX()+1);
    hzns->SetBinContent(zbmn, hzns->GetBinContent(zbmn)+under);
    hzns->SetBinContent(zbmx, hzns->GetBinContent(zbmx)+over);
    printf("Add %e and %e of under/overflows to extreme bins\n",under,over);
    */
    res->AddAtAndExpand(hzns, kZvMCDistNS+shift);
    //
    TH2* mcPrim = 0;
    // prepare MC primary signal histo
    if (normToMB) {
      // don't normalise to NSel
      tmp2   = (TH2*)FindObject(bin,"zvEtaPrimMC", lst, kFALSE );
      mcPrim = (TH2*)CopyAdd(tmp2,
			     Form("bin%d_zvEtaPrimMC",bin),
			     Form("bin%d MC True signal Zv vs Eta",bin),
			     res, kMCPrim, shift);
      mcPrim->Scale(1./nEvTot);
      MyPrint(1,"%s (copy of %s) scaled by %f", mcPrim->GetName(),
	      tmp2->GetName(), nEvTot);
    }
    else {
      // don't normalise to NSel
      tmp2   = (TH2*)FindObject(bin,"zvrecEtaPrimMCSel", lst, kTRUE );
      mcPrim = (TH2*)CopyAdd(tmp2,
			     Form("bin%d_zvEtaPrimMCSel", bin),
			     Form("bin%d MC True signal Zv vs Eta (sel)",bin),
			     res, kMCPrim, shift);
      MyPrint(1,"%s (copy of %s)", mcPrim->GetName(), tmp2->GetName());
    }
  }
  //
  return kTRUE;
}

//_____________________________________________________________________
void ProcessHistos(int bin)
{
  //
  int shift = bin*kNHistos;
  MyPrint(2,"Shift for bin %d is %d", bin, shift);
  //
  TString prefN = "bin"; prefN += bin; prefN+="_";
  TString prefT = "bin"; prefT += bin; prefT+=" ";
  TObjArray* res = &resArr;
  //
  // // build MC truth
  // TH2* mctru2 = (TH2*)res->At(shift + kMCShift + kMCPrim);
  // mctru2 = (TH2*) mctru2->Clone(prefN+"MCtruth2D");
  // mctru2->SetTitle(prefN+"MCtruth2D");
  // TH1F* mctru = mctru2->ProjectionX(prefN+"MCtruth");
  // mctru->SetTitle(prefN+"MCtruth");
  // double binNorm = mctru2->GetXaxis()->GetBinWidth(1)*mctru2->GetYaxis()->GetBinWidth(1);
  // mctru->Scal(1./binNorm);
  // res->AddAtAndExpand(mctru, shift + kMCTruth);
  //
  TH1* hZV = (TH1*)res->At(shift + kZvDist);

  TH1* hZVmcRec = (TH1*)res->At(shift + kMCShift + kZvDist);
  TH1* hZVmcGen = (TH1*)res->At(shift + kMCShift + kZvMCDistNS);
  TH1* mcVtxEff = (TH1*)hZVmcRec->Clone(Form("Eff_%s",hZVmcRec->GetName()));
  MyPrint(1,"Got vertex distributions %s (data) %s %s (MC)",
	  hZV->GetName(), hZVmcRec->GetName(), hZVmcGen->GetName());

  mcVtxEff->Divide(hZVmcRec,hZVmcGen,1,1,"B"); // vertex rec eff
  MyPrint(1,"Calculate vertex efficiency as %s over %s -> %s",
	  hZVmcRec->GetName(), hZVmcGen->GetName(), mcVtxEff->GetName());

  TH1* hZVcorr = (TH1*)hZV->Clone(Form("corr_%s",hZV->GetName()));
  hZVcorr->Divide(mcVtxEff);
  MyPrint(1,"Calculate corrected vertex distribution %s", hZVcorr->GetName());
  //
  res->AddAt(hZVcorr, shift+kZvDistCorr);
  res->AddAt(mcVtxEff, shift+kMCShift+kZvEff);
  MyPrint(2,"Adding vtx corr (%s) and eff (%s) at %d and %d",
	  hZVcorr->GetName(), mcVtxEff->GetName(), shift+kZvDistCorr,
	  shift+kMCShift+kZvEff);
  //
  // build alpha matrix
  TH2* halp = (TH2*)res->At(shift + kMCShift + kMCPrim);
  MyPrint(1,"Start of alpha: %s", halp->GetName());
  halp = (TH2*) halp->Clone(prefN+"Alpha");
  MyPrint(1,"Cloned as %s", halp->GetName());
  halp->SetTitle(prefN+"#alpha");

  TH2* omBeta = (TH2*)res->At(shift + kMCShift + k1MBeta);
  TH2* mcRaw  = (TH2*)res->At(shift + kMCShift + kRawDtCut);
  halp->Divide( omBeta );
  halp->Divide( mcRaw );
  MyPrint(1,"%s scaled by %s and %s", halp->GetName(), omBeta->GetName(),
	  mcRaw->GetName());
  res->AddAtAndExpand(halp, shift + kAlpha);
  MyPrint(2,"Adding alpha %s at %d", halp->GetName(), shift+kAlpha);
  
  //
  KillBadBins(halp,minAlpha,maxAlpha);
  //
  // build alpha matrix with MC labels bg
  TH2* halpMC = (TH2*)res->At(shift + kMCShift + kMCPrim);
  // PrintH(halpMC);
  MyPrint(1,"Start of alpha MC: %s", halpMC->GetName());
  halpMC = (TH2*) halpMC->Clone(prefN + "AlphaMC");
  MyPrint(1,"Cloned as %s", halpMC->GetName());  
  halpMC->SetTitle(prefT + "#alpha MC labels");
  
  TH2* omBetaMC = (TH2*)res->At(shift + kMCShift + k1MBetaMC);
  halpMC->Divide( omBetaMC );
  halpMC->Divide( mcRaw );
  MyPrint(1, "%s scaled by %s and %s", halpMC->GetName(), omBetaMC->GetName(),
	  mcRaw->GetName());
  res->AddAtAndExpand(halpMC, shift + kAlphaMC);
  //
  KillBadBins(halpMC,minAlpha,maxAlpha);
  // PrintH(halpMC);
  //
  // build corrected signal
  TH2* hsigCorr = (TH2*)res->At(shift + kSignalEst);
  MyPrint(1,"Start of corrected signal: %s", hsigCorr->GetName());
  hsigCorr = (TH2*) hsigCorr->Clone(prefN + "SignalEstCorr");
  MyPrint(1,"Cloned as %s", hsigCorr->GetName()); 
  hsigCorr->SetTitle(prefT + "Corrected Signal");

  hsigCorr->Multiply(use1mBeta ? halpMC : halp );
  // PrintH(hsigCorr);
  MyPrint(1,"Multiply by %s", use1mBeta ? halpMC->GetName() : halp->GetName());
  res->AddAtAndExpand(hsigCorr, shift + kSigCorr);
  MyPrint(2,"Add corrected signal %s at %d", hsigCorr->GetName(), shift+kSigCorr);
  //

  TH2* tmpSC = (TH2*)hsigCorr->Clone("tmpSC");
  TH1* useZV = normToMB ? hZVcorr : hZV;
  TH1* hsigCorrX = 0;
  if (useZbinWAv) {
    //    CorrectForZV(tmpSC,hZV);
    CorrectForZV(tmpSC,useZV);
    hsigCorrX = ProjectWghMean(tmpSC,"DataCorrSignalX");
    MyPrint(1,"Corrected X signal from weighted mean: %s",
	    hsigCorrX->GetName());
  }
  else {
    //    hsigCorrX = ProjNorm(tmpSC,hZV,"DataCorrSignalX");
    hsigCorrX = ProjNorm(tmpSC,useZV,"DataCorrSignalX");
    MyPrint(1,"Corrected X signal from normal mean %s", hsigCorrX->GetName());
  }
  //
  delete tmpSC;
  hsigCorrX->Scale(1./hsigCorrX->GetBinWidth(1));
  MyPrint(1,"Scaling corrected X signal %s with bin width",
	  hsigCorrX->GetName());
  
  TF1* pl0 = new TF1("pl0","pol0");
  pl0->SetParameter(0,hsigCorr->GetMinimum());
  hsigCorrX->Fit(pl0,"q0","",-kEtaFitRange,kEtaFitRange);
  MyPrint(1,"Fit straight line to X signal %s", hsigCorrX->GetName());
  double fval = pl0->GetParameter(0);
  double ferr = pl0->GetParError(0);
  delete hsigCorrX;
  dNdEta[bin]    = fval;
  dNdEtaErr[bin] = ferr;
  printf("Bin %d: dN/d#eta_{|#eta|<0.5} = %.2f  %.2f\n",bin, fval,ferr);
  //
}

void PlotResults(Bool_t waitForUser) 
{
  MyGuard(1,"Plottinmg results");
  TString psnm1 = figDir; psnm1 += "/"; psnm1 += uniqueName; 
  psnm1 += "_"; psnm1 += nCentBins; psnm1+= "bins_";
  psnm1 += outStr; psnm1 += ".pdf";
  TString psnm0 = psnm1.Data(); 
  psnm0 += "[";
  TString psnm2 = psnm1.Data(); 
  psnm2 += "]";
  //
  TH1* hstat = (TH1*)FindObject(-1,kHStatName,listDt,kFALSE);
  double myMergeFactor = hstat->GetBinContent(kOneUnit);
  if (myMergeFactor<1) myMergeFactor = 1.;
  //
  TH1* hbins = (TH1*)FindObject(-1,"EvCentr",listDt,kFALSE);
  //
  if (!canvFin) canvFin = new TCanvas("canvFin", "canvFin",0,50,700,1000);
  canvFin->Clear();
  //
  PrintAndPause(canvFin,psnm0.Data(),false);
  //
  canvFin->Divide(1,2);
  canvFin->cd(1);
  gPad->SetLeftMargin(0.15);
  TGraphErrors* grp = new TGraphErrors(nCentBins);
  for (int i=0;i<nCentBins;i++) {
    grp->SetPoint(i,hbins->GetBinCenter(i+1),dNdEta[i]);
    grp->SetPointError(i,hbins->GetBinWidth(i+1)/2,dNdEtaErr[i]);
  }
  grp->SetMarkerStyle(20);
  grp->SetMarkerColor(kRed);
  grp->SetLineColor(kRed);
  grp->SetMinimum(1e-6);
  grp->Draw("ap");
  grp->GetXaxis()->SetTitle("Centrality Variable");
  grp->GetYaxis()->SetTitle("dN/d#eta_{|#eta|<0.5}");
  grp->GetYaxis()->SetTitleOffset(1.8);
  MyPrint(2,"Created graph of central dNch/deta %s", grp->GetName());
  gPad->SetGrid(1,1);
  //
  canvFin->cd(2);
  gPad->SetLeftMargin(0.15);
  hbins->Draw();
  hbins->SetMinimum(1e-6);
  hbins->SetMarkerStyle(20);
  hbins->SetMarkerColor(kRed);
  hbins->SetLineColor(kRed);
  hbins->GetYaxis()->SetTitle("accepted events");
  hbins->GetYaxis()->SetTitleOffset(1.8);
  gPad->SetGrid(1,1);
  MyPrint(2,"Drew %s", hbins->GetName());
  //
  canvFin->cd(0);
  PrintAndPause(canvFin,psnm1.Data(),waitForUser);
  //
  const TArrayD &binArr = *hbins->GetXaxis()->GetXbins();
  //
  {
    MyGuard(1,"Looping over centrality bins");
    for (int i=0;i<nCentBins;i++) {
      //
      sprintf(outTitle,"%s, %d<C_%s<%d, %.1f<#eta<%.1f, %.1f<Z_{V}<%.1f,  Bg.:%s, %s, CutVar:%s, |sig|<%.2f, %.2f<|bg.nrm|<%.2f",
	      uniqueName.Data(),
	      (int)binArr[i],
	      hstat->GetXaxis()->GetBinLabel(kCentVar),
	      (int)binArr[i+1],
	      hstat->GetBinContent(kEtaMin)/myMergeFactor,
	      hstat->GetBinContent(kEtaMax)/myMergeFactor,
	      hstat->GetBinContent(kZVMin)/myMergeFactor,
	      hstat->GetBinContent(kZVMax)/myMergeFactor,	  
	      useBgType.Data(),
	      (use1mBeta ? "MC-BG" : "DT-BG"),
	      useShapeType==kNormShapeDist ? "#Delta":"#Delta#varphi-#delta_{#varphi}",
	      useShapeType==kNormShapeDist ? kWDistSgCut:kdPhiSgCut,
	      useShapeType==kNormShapeDist ? kWDistBgTailMin : kdPhiBgTailMin,
	      useShapeType==kNormShapeDist ? kWDistBgTailMax : kdPhiBgTailMax	  
	      );
      //
      MyPrint(2,"Title: %s", outTitle);
      PlotDNDEta(i);
      PrintAndPause(canvFin,psnm1.Data(),waitForUser);
      PlotAlphaBeta(i);
      PrintAndPause(canvFin,psnm1.Data(),waitForUser);
    }
  }
  MyPrint(2,"Now plotting species");
  PlotSpecies();
  PrintAndPause(canvFin,psnm1.Data(),waitForUser);
  //
  PrintAndPause(canvFin,psnm2.Data(),false);    
}

void PlotDNDEta(int bin)
{
  MyGuard(1,"Plotting dN/deta in bin %d", bin);
  //
  TObjArray *res = &resArr;
  TString prefN = "bin"; prefN += bin; prefN+="_";
  int shift = bin*kNHistos;
  //
  char buff[1000];
  // eta range
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  double mn = 1e6,mx = -1e6;
  if (!canvFin) canvFin = new TCanvas("canvFin", "canvFin",0,50,700,1000);
  canvFin->Clear();
  canvFin->Divide(1,2);
  canvFin->cd(1);
  gPad->SetLeftMargin(0.15);
  //
  TH1* hZV     = (TH1*)res->At(shift + kZvDist);
  TH1* hZVcorr = (TH1*)res->At(shift + kZvDistCorr);
  if (!res->At(shift + kSigCorr)) {
    MyPrint(2,"No histogram at %d", shift+kSigCorr);
    return;
  }
  MyPrint(2,"Got Z vertex and correction histograms: %s %s", hZV->GetName(), hZVcorr->GetName());
  
  // corrected data
  TString nms =  prefN;
  nms += "DataCorrSignal";
  nms += "_"; 
  nms += uniqueName;
  TH2* tmpSCorr = (TH2F*)res->At(shift + kSigCorr)->Clone(Form("%stmpSCorrBin%d",prefN.Data(),bin));
  TH1* hsigCorr = 0;
  TH1* useZV = normToMB ? hZVcorr : hZV;
  MyPrint(2,"Got corrrected signal at %d (%s): %s", shift + kSigCorr,
	  res->At(shift + kSigCorr)->GetName(), tmpSCorr->GetName());
  if (useZbinWAv) {
    //    CorrectForZV(tmpSCorr,hZV);
    CorrectForZV(tmpSCorr,useZV);
    hsigCorr = ProjectWghMean( tmpSCorr, nms.Data());
  }
  else {
    //    hsigCorr = ProjNorm(tmpSCorr,hZV,nms.Data());
    hsigCorr = ProjNorm(tmpSCorr,useZV,nms.Data());
    MyPrint(2,"Corrected X signal from normal mean %s", hsigCorr->GetName());
  }
  //  delete tmpSCorr;
  SetHStyle(hsigCorr,kRed,20,1.0);
  hsigCorr->Scale(1./hsigCorr->GetBinWidth(1));
  hsigCorr->Draw();
  MyPrint(2,"Scaled %s by bin width, and drew", hsigCorr->GetName());
  mx = TMath::Max(mx, hsigCorr->GetMaximum());
  mn = TMath::Min(mn, hsigCorr->GetMinimum());
  resDnDeta.AddAtAndExpand( hsigCorr, bin );
  MyPrint(2,"Added %s at %d", hsigCorr->GetName(), bin);
  hsigCorr->SetDirectory(0);
  resDnDeta.AddAtAndExpand( res->At(shift + kSigCorr), 100+bin );
  MyPrint(2,"Added unnormalized %s at %d", hsigCorr->GetName(), bin);
  resDnDeta.AddAtAndExpand( tmpSCorr, 200+bin );
  tmpSCorr->SetDirectory(0);
  resDnDeta.AddAtAndExpand(hZV , 250+bin );
  MyPrint(2,"Added vertex %s at %d", hZV->GetName(), 250+bin);
  hZV->SetDirectory(0);
  
  //
  // raw data
  TH2* tmpRaw = (TH2F*)res->At(shift+kRawDtCut)->Clone("tmpRaw");
  TH1* hraw = tmpRaw->ProjectionX(prefN+"DataRaw"); 
  delete tmpRaw;
  MyPrint(2,"Got the raw data %s and made projection %s",
	  res->At(shift+kRawDtCut)->GetName(), hraw->GetName());
  SetHStyle(hraw,kRed,21,1.0);
  hraw->Scale(1./hraw->GetBinWidth(1));
  hraw->Draw("same");
  mn = TMath::Min(mn, hraw->GetMinimum());
  mx = TMath::Max(mx, hraw->GetMaximum());
  //  
  // raw data bg sub
  TH2* tmpRawB = (TH2F*)res->At(shift+kSignalEst)->Clone("tmpRawB");
  TH1* hraws = tmpRawB->ProjectionX(prefN+"DataRawSub");
  delete tmpRawB;
  SetHStyle(hraws,kRed,23,1.0);
  MyPrint(2,"Got the raw data %s and made projection %s",
	  res->At(shift+kSignalEst)->GetName(), hraws->GetName());
  hraws->Scale(1./hraw->GetBinWidth(1));
  hraws->Draw("same");
  mn = TMath::Min(mn, hraw->GetMinimum());
  mx = TMath::Max(mx, hraw->GetMaximum());
  //
  // bg
  TH2* tmpBg = (TH2F*)res->At(shift+kBgEst)->Clone("tmpBg");
  TH1* hbg = tmpBg->ProjectionX(prefN+"BgEst"); 
  delete tmpBg;
  MyPrint(2,"Got the bg estimate %s and made projection %s",
	  res->At(shift+kBgEst)->GetName(), hbg->GetName());
  SetHStyle(hbg,kMagenta,22,1.0);
  hbg->Scale(1./hbg->GetBinWidth(1));
  hbg->Draw("same");
  mn = TMath::Min(mn, hbg->GetMinimum());
  mx = TMath::Max(mx, hbg->GetMaximum());
  //
  // mc part ----------------------------
  TH1* hZVMC   = (TH1*)res->At(shift + kZvDist + kMCShift);
  TH1* hZVMCNS = (TH1*)res->At(shift + kZvMCDistNS + kMCShift);
  //
  // raw data
  TH2* tmpRawMC = (TH2F*)res->At(shift+kRawDtCut+kMCShift)->Clone("tmpRawMC");
  TH1* hrawMC = tmpRawMC->ProjectionX(prefN+"DataRawMC"); 
  delete tmpRawMC;
  MyPrint(2,"Got the MC raw %s and made projection %s",
	  res->At(shift+kRawDtCut+kMCShift)->GetName(), hrawMC->GetName());
  SetHStyle(hrawMC,kBlue,24,1.0);
  hrawMC->Scale(1./hrawMC->GetBinWidth(1));
  hrawMC->Draw("same");
  mn = TMath::Min(mn, hrawMC->GetMinimum());
  mx = TMath::Max(mx, hrawMC->GetMaximum());
  //  
  // raw data bg sub
  TH2* tmpRawBMC = (TH2F*)res->At(shift+kSignalEst+kMCShift)->Clone("tmpRawBMC");
  TH1* hrawsMC = tmpRawBMC->ProjectionX(prefN+"DataRawSubMC"); 
  delete tmpRawBMC;
  MyPrint(2,"Got the MC signal estimate %s and made projection %s",
	  res->At(shift+kSignalEst+kMCShift)->GetName(), hrawsMC->GetName());
  SetHStyle(hrawsMC,kBlue,26,1.0);
  hrawsMC->Scale(1./hrawMC->GetBinWidth(1));
  hrawsMC->Draw("same");
  mn = TMath::Min(mn, hrawMC->GetMinimum());
  mx = TMath::Max(mx, hrawMC->GetMaximum());
  //
  // raw data bgMClabels sub
  TH2* tmpRawLBMC = (TH2F*)res->At(shift+kSignalEstMC+kMCShift)->Clone("tmpRawLBMC");
  TH1* hrawsMCLB =  tmpRawLBMC->ProjectionX(prefN+"DataRawSubMCLB");
  delete tmpRawLBMC;
  MyPrint(2,"Got the MC signal estimate %s from labels and made projection %s",
	  res->At(shift+kSignalEstMC+kMCShift)->GetName(), hrawsMCLB->GetName());
  SetHStyle(hrawsMCLB,kGreen+2,30,1.0);
  hrawsMCLB->Scale(1./hrawsMCLB->GetBinWidth(1));
  hrawsMCLB->Draw("same");
  mn = TMath::Min(mn, hrawsMCLB->GetMinimum());
  mx = TMath::Max(mx, hrawsMCLB->GetMaximum());
  //
  // bg est
  TH2* tmpBgMCEst = (TH2F*)res->At(shift+kBgEst+kMCShift)->Clone("tmpRawBbMCEst");
  TH1* hbgMCEst = tmpBgMCEst->ProjectionX(prefN+"BgEstMC");
  MyPrint(2,"Got the MC bg estimate %s and made projection %s",
	  res->At(shift+kBgEst+kMCShift)->GetName(), hbgMCEst->GetName());
  delete tmpBgMCEst;
  SetHStyle(hbgMCEst,kBlue,26,1.0);
  hbgMCEst->Scale(1./hbgMCEst->GetBinWidth(1));
  hbgMCEst->Draw("same");
  mn = TMath::Min(mn, hbgMCEst->GetMinimum());
  mx = TMath::Max(mx, hbgMCEst->GetMaximum());
  //  
  // bg MC
  TH2* tmpBgMC = (TH2F*)res->At(shift+kBgMC+kMCShift)->Clone("tmpBgMC");
  TH1* hbgMC = tmpBgMC->ProjectionX(prefN+"BgMC"); 
  delete tmpBgMC;
  MyPrint(2,"Got the MC bg %s and made projection %s",
	  res->At(shift+kBgMC+kMCShift)->GetName(), hbgMC->GetName());
  SetHStyle(hbgMC,kGreen+2,25,1.0);
  hbgMC->Scale(1./hbgMC->GetBinWidth(1));
  hbgMC->Draw("same");
  mn = TMath::Min(mn, hbgMC->GetMinimum());
  mx = TMath::Max(mx, hbgMC->GetMaximum());
  //
  // mc truth
  TH2* tmpMCTrue = (TH2F*)res->At(shift+kMCPrim+kMCShift)->Clone(Form("%stmpMCTrueBin%d",prefN.Data(),bin));
  TH1* hMCtrue = 0;
  TH1* hzMC = normToMB ? hZVMCNS : hZVMC; //!!!RS
  if (useZbinWAv) {
    CorrectForZV(tmpMCTrue,hzMC);
    hMCtrue = ProjectWghMean(tmpMCTrue, prefN+"MCTruth");
  }
  else {
    hMCtrue = ProjNorm(tmpMCTrue,hzMC,prefN+"MCTruth");
    //    hMCtrue = ProjNorm(tmpMCTrue,hZVMCNS,prefN+"MCTruth");
  }
  //  delete tmpMCTrue;
  MyPrint(2,"Got the MC truth %s and made projection %s",
	  res->At(shift+kMCPrim+kMCShift)->GetName(), hMCtrue->GetName());
  SetHStyle(hMCtrue,kCyan,20,0.8);
  hMCtrue->Scale(1./hMCtrue->GetBinWidth(1));
  hMCtrue->SetLineStyle(1);
  hMCtrue->Draw("same");
  mn = TMath::Min(mn, hMCtrue->GetMinimum());
  mx = TMath::Max(mx, hMCtrue->GetMaximum());
  //
  hMCtrue->SetDirectory(0);
  tmpMCTrue->SetDirectory(0);
  resDnDeta.AddAtAndExpand( hMCtrue, 300+bin);
  resDnDeta.AddAtAndExpand( res->At(shift+kMCPrim+kMCShift), 400+bin );
  resDnDeta.AddAtAndExpand( tmpMCTrue, 500+bin );
  resDnDeta.AddAtAndExpand( hZVMCNS, 550+bin );
  //

  mn = 0;
  hsigCorr->SetMinimum(mn);
  hsigCorr->SetMaximum(mx + 0.4*(mx-mn));
  //
  char ftres[1000];
  sprintf(ftres,"dN/d#eta_{|#eta|<%.2f} = %.2f #pm %.2f",kEtaFitRange,dNdEta[bin],dNdEtaErr[bin]);  
  double xtxt = (hsigCorr->GetXaxis()->GetXmin()+hsigCorr->GetXaxis()->GetXmax())/2;
  TLatex *txfit = new TLatex(xtxt,mn+0.25*(mx-mn), ftres);
  txfit->SetTextSize(0.04);
  txfit->Draw();


  gPad->Modified();
  //
  MyPrint(2,"Making the legend");
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
  //
  legDnDetaMC->Draw();
  //
  gPad->SetGrid(1.1);
  gPad->Modified();
  AddLabel(outTitle,0.1,0.97, kBlack,0.02);
  //
  canvFin->cd();
  //
  //---------------- dsitributions
  canvFin->cd(2);
  //
  MyPrint(2,"Drawing AUX distributions");
  TH1* mcdst = (TH1*)res->At(shift+kDataDist+kMCShift);
  TH1* mcdstbg = (TH1*)res->At(shift+kBgDist+kMCShift);
  TH1* mcdstbgLB = (TH1*)res->At(shift+kBgMCDist+kMCShift);
  TH1* dtdst = (TH1*)res->At(shift+kDataDist);
  TH1* dtdstbg = (TH1*)res->At(shift+kBgDist);
  //
  TH1* mcDstN     = (TH1*)FindObject(bin,HName("Data", useShapeType==kNormShapeDist ? "WDist":"DPhiS"), listMC );
  TH1* mcDstSec   = (TH1*)FindObject(bin,HName("Sec",  useShapeType==kNormShapeDist ? "WDist":"DPhiS"), listMC );
  TH1* mcDstCombU = (TH1*)FindObject(bin,HName("CombU",useShapeType==kNormShapeDist ? "WDist":"DPhiS"), listMC );
  TH1* mcDstCombC = (TH1*)FindObject(bin,HName("Comb", useShapeType==kNormShapeDist ? "WDist":"DPhiS"), listMC );
  //  
  double scl,sclE;
  mcDstN = NormalizeBg(mcdst,mcDstN,scl,sclE);
  mcDstSec->Scale(scl);
  mcDstCombU->Scale(scl);
  mcDstCombC->Scale(scl);
  mcDstCombC->Add(mcDstCombU,-1);

  dtdst->Draw("");
  gPad->Modified();
  dtdst->GetXaxis()->SetLabelSize(0.03);
  dtdst->GetXaxis()->SetTitleSize(0.03);
  dtdst->GetXaxis()->SetTitleOffset(2);
  dtdstbg->Draw("same");
  //
  mcdst->Draw("same");
  mcDstSec->Draw("same");
  if (mcdstbgLB) mcdstbgLB->Draw("same");
  mcdstbg->Draw("same");
  mcDstCombC->Draw("same");
  //

  SetHStyle(mcdst,kBlue, 25,0.7);
  if (mcdstbgLB) SetHStyle(mcdstbgLB,kGreen, 7/*32*/,0.5);
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
  if (mcdstbgLB) {
    Integrate(mcdstbgLB, cutSgMin,cutSgMax, vmcCmb,vmcCmbE);  
    GetRatE(vmcCmb,vmcCmbE, vmcTot,vmcTotE, ratCmb,ratCmbE);
  }
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
  canvFin->cd();
  //
  if (creatDnDEtaCMacro) {
    sprintf(buff,"%s/%s_b%d_dNdEta_%s",figDir,uniqueName.Data(),bin,outStr);
    SaveCanvas(canvFin,buff,"cg");
  }
  //
}
//
void PlotAlphaBeta(int bin) 
{
  MyGuard(2,"Draw alpha and beta");
  char buff[1000];
  int shift = bin*kNHistos;
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TObjArray* res = &resArr;
  //------------------------------------------------------
  if (!canvFin) canvFin = new TCanvas("canvFin","canvFin",10,10,700,1000);
  canvFin->Clear();
  canvFin->Divide(2,3,0.01,0.01);
  canvFin->cd(1);
  TH1* dtBet = (TH1*)res->At(shift + k1MBeta);
  TH1* mcBet = (TH1*)res->At(shift + k1MBeta+kMCShift);
  TH1* mcBetLB = (TH1*)res->At(shift + k1MBetaMC+kMCShift);
  if (!dtBet) return;
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
  canvFin->cd(1);
  gPad->SetRightMargin(0.15);
  dtBet->Draw("colz");
  AddLabel("#beta Data",0.2,0.95,kBlack,0.04);
  gPad->Modified();
  dtBet->GetYaxis()->SetTitleOffset(1.4);

  // TPaletteAxis *p = (TPaletteAxis*)dtBet->FindObject("palette");
  // if (p) p->SetX1NDC(0.85);
  canvFin->cd(2);
  gPad->SetRightMargin(0.15);
  mcBet->Draw("colz");
  AddLabel("#beta MC (bckg.estimated)",0.2,0.95,kBlack,0.04);
  gPad->Modified();
  mcBet->GetYaxis()->SetTitleOffset(1.4);
  // p = (TPaletteAxis*)mcBet->FindObject("palette");
  // if (p) p->SetX1NDC(0.85);
  canvFin->cd(3);
  gPad->SetRightMargin(0.15);
  mcBetLB->Draw("colz");
  AddLabel("#beta MC (bckg.from MC labels)",0.2,0.95,kBlack,0.04);
  gPad->Modified();
  mcBetLB->GetYaxis()->SetTitleOffset(1.4);
  // p = (TPaletteAxis*)mcBetLB->FindObject("palette");
  /// if (p) p->SetX1NDC(0.85);
  //
  //------------------------------------------------------
  TH1* dtAlp = (TH1*)res->At(shift + kAlpha);
  TH1* mcAlp = (TH1*)res->At(shift + kAlphaMC);
  GetRealMinMax(dtAlp,mn,mx);
  GetRealMinMax(mcAlp,mnt,mxt);
  if (mnt<mn) mn = mnt;
  if (mxt>mx) mx = mxt;
  dtAlp->SetMinimum(mn - 0.05*(mx-mn));
  dtAlp->SetMaximum(mx + 0.05*(mx-mn));
  mcAlp->SetMinimum(mn - 0.05*(mx-mn));
  mcAlp->SetMaximum(mx + 0.05*(mx-mn));
  //
  canvFin->cd(4);
  gPad->SetRightMargin(0.15);
  dtAlp->Draw("colz");
  AddLabel("#alpha (bckg.estimated)",0.2,0.95,kBlack,0.04);
  gPad->Modified();
  dtAlp->GetYaxis()->SetTitleOffset(1.4);
  // TPaletteAxis *pa = (TPaletteAxis*)dtBet->FindObject("palette");
  // if (pa) pa->SetX1NDC(0.85);
  canvFin->cd(5);
  gPad->SetRightMargin(0.15);
  mcAlp->Draw("colz");
  AddLabel("#alpha (bckg.from MC labels)",0.2,0.95,kBlack,0.04);
  gPad->Modified();
  mcAlp->GetYaxis()->SetTitleOffset(1.4);
  // pa = (TPaletteAxis*)mcBet->FindObject("palette");
  // if (pa) pa->SetX1NDC(0.85);
  gPad->Modified();
  canvFin->cd(6);
  AddLabel(outTitle,0.1,0.5, kBlack, 0.02);
  //
  if (creatAlphaBetaCMacro) {
    sprintf(buff,"%s/%sAlphaBeta_%s",figDir,uniqueName.Data(),outStr);
    SaveCanvas(canvFin,buff,"cg");
  }
  //
}

void PlotSpecies() 
{
  MyGuard(2,"Plotting species");
  char buff[1000];
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //------------------------------------------------------
  TH2F* hSpecPrim  = (TH2F*)FindObject(-1, "pdgPrim", listMC,kFALSE);
  TH2F* hSpecSec   = (TH2F*)FindObject(-1, "pdgSec", listMC,kFALSE);
  TH2F* hSpecPrimP = (TH2F*)FindObject(-1, "pdgPrimPar", listMC,kFALSE);
  TH2F* hSpecSecP  = (TH2F*)FindObject(-1, "pdgSecPar", listMC,kFALSE);
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
  if (hSpecPrimSel->Integral()<1) return;
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
  //
  if (!canvFin) canvFin = new TCanvas("canvFin","canvFin",10,10,1100,800);
  canvFin->Clear();
  canvFin->Divide(1,2,0.01,0.01);
  canvFin->cd(1);
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
  canvFin->cd(2);
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
  canvFin->cd(1);
  //  AddLabel(outTitle,0.1,0.97, kBlack, 0.02);
  canvFin->cd();
  //
  if (creatSpeciesCMacro) {
    sprintf(buff,"%s/%sSpecies_%s",figDir,uniqueName.Data(),outStr);
    SaveCanvas(canvFin,buff,"cg");
  }
}

//______________________________________________________________________
void CropHisto(TH1* histo, int bx0, int bx1, int by0, int by1)	       
{
  MyGuard(3,"Cropping histogram %s between bins %d,%d x %d,%d", bx0, bx1, by0, by1);
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
  MyGuard(3,"Cropping histogram %s between values %f,%f x %f,%f", vx0, vx1, vy0, vy1);
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
  // match generated bg and data tails, calculate normalization,
  // return normalized bg copy
  //
  MyGuard(1,"Normalize background from %s and %s",
	  dataH->GetName(), bgH->GetName());
  
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
  else {
    MyPrint(2,"NormalizeBg: unknown shape type %d",useShapeType);
    return 0;
  }
  MyPrint(2,"NormalizeBg: bins for tails: right: %d:%d / left: %d:%d",
	  bgBins[0][0],bgBins[0][1],bgBins[1][0],bgBins[1][1]);
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
  MyPrint(1,"NormalizeBg: Tails scaling %s wrt %s: "
	  "Wgh.Mean:%.4f(%.4f) / Integral:%.4f(%.4f)",
	  bgH->GetName(),dataH->GetName(), meanR,meanRE, meanRI,meanRIE);
  MyPrint(2,"NormalizeBg: Select scaling type %s",
	  useScaleType==kSclWghMean ? "Wgh.Mean":"Integral");
  //
  scl  = useScaleType==kSclWghMean ? meanR  : meanRI;
  sclE = useScaleType==kSclWghMean ? meanRE : meanRIE;
  //
  // rescaled bg
  char buff[1000];
  sprintf(buff,"%s_bgNorm",bgH->GetName());
  TH1* tmp = bgH;
  bgH = (TH1*)tmp->Clone(buff);
  sprintf(buff,"%s bgNorm%d %.4f+-%.4f",bgH->GetName(),useScaleType,scl,sclE);
  TH1* dumH = (TH1*)bgH->Clone("dummySCL$"); dumH->Reset();
  for (int i=1;i<=nbtot;i++) {
    dumH->SetBinContent(i,scl);
    dumH->SetBinError(i,sclE);
  }
  bgH->Multiply(dumH);
  MyPrint(1,"Returning normal backg delta distribution %s (copy of %s)",
	  bgH->GetName(), tmp->GetName());
  delete dumH;
  return bgH;
}

//______________________________________________________________________
TObject* FindObject(int bin, const char* nameH, const TList* lst, Bool_t normPerEvent)
{
  // get histo, optionally normalizing it per processed event
  if (!lst) {
    Warning("FindObject","%s: No list provided",nameH);
    return 0;
  }
  
  int nent = lst->GetEntries();
  char buff[200];
  if (bin>=0) sprintf(buff,"b%d_%s",bin,nameH);
  else        sprintf(buff,"%s",nameH);

  MyGuard(1,"Looking for bin %d of %s: %s", bin, nameH, buff);

  TString nm;
  TObject *hst = 0;
  for (int i=nent;i--;) {
    nm = lst->At(i)->GetName();
    if (nm.EndsWith(buff)) {hst = lst->At(i); break;}
  }
  if (!hst) {
    Warning("FindObject",
	    "No bin %d %s histo in list %s\n",bin,nameH,lst->GetName());
    lst->Print();
    return 0;
  }
  // already normalized
  if (!normPerEvent || hst->TestBit(kBitNormPerEvent)) return hst; 

  TString nameHS = nameH;
  // never scale stat. histo
  if (nameHS==kHStatName) return hst;

  // 
  TH1* hstat = (TH1*)FindObject(-1,kHStatName,lst,kFALSE);
  // MyPrint(2,"Find statistics histogram: %p", hstat);
  
  //double nrm = hstat->GetBinContent(kBinEntries +
  //                                  kEvProcInj+bin*kEntriesPerBin);
  double nrm = hstat->GetBinContent(kBinEntries + kEvProcData+
				    bin*kEntriesPerBin); // HACK
  if (nrm<1) {
    Warning("FindObject",
	    "%s: Anomaluous %d number of events processed in bin %d of list %p",
	    buff,int(nrm),bin,lst);
    return 0;
  }
  MyPrint(2,"Normalization: %f", nrm);
  
  //
  if      (hst->InheritsFrom(TH1::Class())) {
    MyPrint(0,"Scaling histogram %s by %f", hst->GetName(), nrm);
    ((TH1*)hst)->Scale(1./nrm);
  }
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
  MyPrint(2,"Flagging %s a scaled", hst->GetName());
  hst->SetBit(kBitNormPerEvent);
  return hst;
}

//______________________________________________________________________
TList* LoadList(const char* flName, const char* addPref, const char* nameL)
{
  MyPrint(2,"Loading lists from %s (prefix=%s, name=%s)", flName, addPref, nameL);
  
  // load list with histos
  TString nms = flName;
  gSystem->ExpandPathName(nms);

  TFile* fl = TFile::Open(nms.Data());
  if (!fl) {
    Error("LoadList","No file %s\n",nms.Data());
    return 0;
  }

  TList* lst = (TList*)fl->Get(nameL);
  if (!lst) {
    Error("LoadList","No list %s in file %s\n",nameL,nms.Data());
    return 0;
  }

  MyPrint(2,"Setting name to %s", flName);
  lst->SetName(flName);
  int nEnt = lst->GetSize();
  TString nm;
  for (int i=0;i<nEnt;i++) {
    TNamed* ob = (TNamed*)lst->At(i);
    nm = addPref;
    nm += ob->GetName();
    MyPrint(2," Renaming %40s to %s", ob->GetName(), nm.Data());
    ob->SetName(nm.Data());
  }
  //
  return lst;
}

//____________________________________________________________________________
void GetRatE(double x,double xe, double y,double ye, double &rat, double &rate)
{
  MyGuard(10,"Calculating ratio of %f +/- %f over %f +/- %f", x, xe, y, ye);
  rat = 0; rate = 0;
  if (TMath::Abs(y)<1e-16 || TMath::Abs(x)<1e-16) return;
  rat = x/y;
  rate = rat*TMath::Sqrt( xe*xe/(x*x) + ye*ye/(y*y));
  MyPrint(10,"Result: %f +/- %f", rat, rate);
}

//____________________________________________________________________________
void Integrate(TH1* hist, double xmn,double xmx, double &val, double& err)
{
  MyGuard(1,"Integrate %s over [%f,%f]", hist->GetName(), xmn, xmx);
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
  MyPrint(1,"Integral: %f +/- %f", val, err);
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
  MyGuard(3,"Check input list %s for %s", lst->GetName(), dtType);
  // check if needed bg was generated
  TH1* hstat = (TH1*)FindObject(-1,kHStatName,lst);
  TString dts = dtType;
  if (dts=="Data") return int( hstat->GetBinContent(kEvProcData) );
  if (dts=="Mix")  return int( hstat->GetBinContent(kEvProcMix) );
  if (dts=="Inj")  return int( hstat->GetBinContent(kEvProcInj) );
  if (dts=="Rot")  return int( hstat->GetBinContent(kEvProcRot) );
  MyPrint(3,"Unknown process %s statistics is checked. Alowed: Data,Mix,Inj,Rot",dtType);
  return 0;
}


void GetRealMinMax(TH1* histo, double &vmn, double &vmx)
{
  MyGuard(3,"Get real min and max of %s", histo->GetName());
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

//____________________________________________________________________________
void CorrectForZV(TH2* hEtaZ, TH1* hZv)
{
  MyGuard(2,"Correcting %s for vertex efficiency using %s",
	  hEtaZ->GetName(), hZv->GetName());
  // correct for the nonuniformity of Zvertices along Z
  int nbEta = hEtaZ->GetNbinsX();
  int nbZv  = hEtaZ->GetNbinsY();
  //
  //  printf("CorrectForZV: %s %s\n",hEtaZ->GetName(), hZv->GetName());
  TAxis* zax = hEtaZ->GetYaxis();
  for (int ibz=1;ibz<=nbZv;ibz++) {
    double z = zax->GetBinCenter(ibz);
    int izbh = hZv->FindBin(z);
    double scl  = hZv->GetBinContent(izbh);
    double scle = hZv->GetBinError(izbh);
    MyPrint(2,"Vertex correction for z=%f: %f +/- %f", z, scl, scle);
    //scle = 0; 
    if (scl<1e-9) {
      MyPrint(2,"Did not find ZV weight for Z=%f, stop",z);
      continue;
    }
    for (int ibe=1;ibe<=nbEta;ibe++) {
      double val = hEtaZ->GetBinContent(ibe,ibz);
      double err = hEtaZ->GetBinError(ibe,ibz);
      //      err = 0;
      double rat = val/scl;
      double ratE = (val>0 && scl>0) ? rat*TMath::Sqrt((err*err)/(val*val) + (scle*scle)/(scl*scl) ) : 0.0;
      //      double eta = hEtaZ->GetXaxis()->GetBinCenter(ibe);
      //      printf("Zc: %+.5e Eta: %+.5e  Val: %+.5e Scl:%+.5e ->%.5e\n",z,eta,val,scl,val/scl);
      hEtaZ->SetBinContent(ibe,ibz, rat);
      hEtaZ->SetBinError(ibe,ibz, ratE);
    }
  }
}

//____________________________________________________________________________
TH1* ProjNorm(TH2* hEtaZ, TH1* hZv,const char* name, Int_t firstbin, Int_t lastbin)
{
  MyGuard(1,"Projection of %s - %s", hEtaZ->GetName(), hEtaZ->GetTitle());
  // correct for the nonuniformity of Zvertices along Z
  int nbEta = hEtaZ->GetNbinsX();
  int nbZv  = hEtaZ->GetNbinsY();
  if (firstbin<1) firstbin = 1;
  if (lastbin<firstbin) lastbin = nbZv;
  //
  TH1* prj = hEtaZ->ProjectionX(name,firstbin,lastbin,"e");
  prj->Reset();
  MyPrint(2,"Made projection and reset it");
  double *arrV = new double[nbZv];
  double *arrE = new double[nbZv];
  double *arZV = new double[nbZv];
  double *arZE = new double[nbZv];
  double *relE = new double[nbZv];
  int    *ind  = new int[nbZv];
  //
  TAxis* zax = hEtaZ->GetYaxis();
  for (int ibe=1;ibe<=nbEta;ibe++) {
    //
    int cnt = 0;
    for (int ibz=firstbin;ibz<=lastbin;ibz++) {
      double val = hEtaZ->GetBinContent(ibe,ibz);
      double err = hEtaZ->GetBinError(ibe,ibz);
      if (val<1e-9) continue;
      //
      double z = zax->GetBinCenter(ibz);
      int izbh = hZv->FindBin(z);
      //
      arrV[cnt] = val;
      arrE[cnt] = err;
      relE[cnt] = err/val;
      arZV[cnt] = hZv->GetBinContent(izbh);
      arZE[cnt] = hZv->GetBinError(izbh);
      cnt++;
    }
    // suppress bins increasing the error
    TMath::Sort(cnt,relE,ind,kFALSE);
    double av=0,ave=0,vsm=0,vsme=0,zsm=0,zsme=0;    
    double rat = 1e99;
    int iu=0;
    for (iu=0;iu<cnt;iu++) {
      int j = ind[iu];
      double rt = TMath::Sqrt(vsme+arrE[j]*arrE[j])/(vsm+ arrV[j]);
      // printf("%2d(%2d/%2d) New:%.e Old:%.e  | %e %e %e\n",ibe,j,iu,rt,rat,arrV[j],arrE[j],arrE[j]/arrV[j]);
      MyPrint(10,"  %10s - bin %3d,%3d %9.3f+/-%9.3f (%5.2f%%) [%9.3f+/-%9.3f] "
	      "-> %9.3f (%9.3g)",
	      hEtaZ->GetName(), ibe, j, arrV[j], arrE[j],
	      100*relE[j], arZV[j], arZE[j], rt, rat);
      if (rt>rat) continue;//break;
      rat = rt;
      vsm += arrV[j];
      vsme+= arrE[j]*arrE[j];
      zsm += arZV[j];
      zsme+= arZE[j]*arZE[j];
      //      printf("%2d(%2d/%2d) New:%.e Old:%.e  | %e %e %e | %e %e | %e %e\n",ibe,j,iu,rt,rat,arrV[j],arrE[j],arrE[j]/arrV[j], 
      //	     vsm,TMath::Sqrt(vsme), zsm, TMath::Sqrt(zsme) );
      //
    }
    if (iu==0) continue; // no contribution
    av = vsm/zsm;
    zsme /= zsm*zsm;
    vsme /= vsm*vsm;
    ave = av*TMath::Sqrt( vsme + zsme );
    MyPrint(10,
	    "%10s - bin %3d (%+6.3f) count=%2d  n=%9.3f+/-%9.3f d=%9.3f+/-%9.3f"
	    " -> %9.3f +/- %9.3f", hEtaZ->GetName(), ibe,
	    hEtaZ->GetXaxis()->GetBinCenter(ibe), iu, vsm, TMath::Sqrt(vsme),
	    zsm, TMath::Sqrt(zsme), av, ave);
    //    printf("->%d %.4e(%.4e) %2d bins out of %2d\n",ibe, av,ave,iu,cnt);
    prj->SetBinContent(ibe, av);
    prj->SetBinError(ibe, ave);
  }
  //
  delete[] arrV;
  delete[] arrE;
  delete[] arZV;
  delete[] arZE;
  delete[] relE;
  delete[] ind;
  //
  return prj;
}

//____________________________________________________________________________
TH1* ProjectWghMean(TH2* hEtaZ, const char* name, Int_t firstbin, Int_t lastbin, double rejOutliers)
{
  MyGuard(2,"Calculating weighted mean projection of %s from %d to %d",
	  hEtaZ, firstbin, lastbin);
  // for each X bin calculated the weighted average over Y bins
  int nbEta = hEtaZ->GetNbinsX();
  int nbZv  = hEtaZ->GetNbinsY();
  if (firstbin<1) firstbin = 1;
  if (lastbin<firstbin) lastbin = nbZv;  
  //
  MyPrint(2,"HISTO: %s %s %p | bins: %d %d",hEtaZ->GetName(),hEtaZ->GetTitle(), hEtaZ, firstbin, lastbin);
  TArrayD val(nbZv),err(nbZv),prc(nbZv);
  TArrayI ind(nbZv);
  //
  TH1* prj = hEtaZ->ProjectionX(name,firstbin,lastbin,"e");
  prj->Reset();
  MyPrint(2,"Made projection, and reset it");
  for (int ibe=1;ibe<=nbEta;ibe++) {
    int cnt = 0;
    double valAv = 0;
    double errAv = 0;
    //
    for (int ibz=firstbin;ibz<=lastbin;ibz++) {
      val[cnt] = hEtaZ->GetBinContent(ibe,ibz);
      err[cnt] = hEtaZ->GetBinError(ibe,ibz);
      if (!err[cnt]>0) continue;
      valAv += val[cnt]/(err[cnt]*err[cnt]);
      errAv += 1./(err[cnt]*err[cnt]);
      // if (errAv>0) printf("%2d| (%+e %+e) %+e %+e\n",cnt,val[cnt],err[cnt],valAv/errAv,1./sqrt(errAv));
      cnt++;
    }
    if (cnt<1) continue;
    valAv /= errAv;
    errAv = 1./TMath::Sqrt(errAv);
    //    printf("#%3d Full W.Mean for %3d bins: %+e +- %e\n",ibe,cnt,valAv,errAv);
    //
    int cntSkip=-1;
    int useCnt = 0;
    double valFin=0,errFin=0;
    while(cntSkip<cnt-1) {
      // get truncated mean
      TMath::Sort(cnt,val.GetArray(),ind.GetArray(), kFALSE); // sort in increasing value
      if (cntSkip>-1) err[ind[cntSkip]] = -1;
      valAv = 0;
      errAv = 0;
      for (int ic=0;ic<cnt;ic++) {
	double erb = err[ind[ic]];
	double vlb = val[ind[ic]];
	if (!(erb>0)) continue;
	valAv += vlb/(erb*erb);
	errAv += 1./(erb*erb);
	useCnt++;
      }
      if (useCnt<1 || !(errAv>0)) break;
      valAv /= errAv;
      errAv = 1./TMath::Sqrt(errAv);
      //      printf("#%3d Truncated0 W.Mean for %d out of %d bins: %+e +- %e\n",ibe,useCnt,cnt,valAv,errAv);
      //
      for (int ic=0;ic<cnt;ic++) prc[ic] = err[ic]>0 ? TMath::Abs(valAv-val[ic])/err[ic] : 1e9;
      // take most precize values    
      TMath::Sort(cnt,prc.GetArray(),ind.GetArray(), kFALSE); // sort in increasing error
      useCnt = cnt>10 ? cnt/2 : cnt*0.7;
      if (useCnt<1) useCnt = 1;
      valAv = 0;
      errAv = 0;
      for (int ic=0;ic<useCnt;ic++) {
	double erb = err[ind[ic]];
	double vlb = val[ind[ic]];
	if (!(erb>0)) continue;
	valAv += vlb/(erb*erb);
	errAv += 1./(erb*erb);
      }
      if (!(errAv>0)) break;
      valAv /= errAv;
      errAv = 1./TMath::Sqrt(errAv);
      //      printf("#%3d Truncated W.Mean for %d out of %d bins: %+e +- %e\n",ibe,useCnt,cnt,valAv,errAv);
      valFin=0;
      errFin=0;
      useCnt = 0;
      for (int ic=0;ic<cnt;ic++) {
	double erb = err[ic];
	double vlb = val[ic];
	if (!(erb>0)) continue;
	double dev = TMath::Abs(valAv-vlb)/erb;
	if (rejOutliers>0 && dev>rejOutliers) continue;
	//
	valFin += vlb/(erb*erb);
	errFin += 1./(erb*erb);
	useCnt++;
      }
      if (useCnt<1 || !(errFin)>0) {
	cntSkip++;
	continue;
      }
      valFin /= errFin;
      errFin = 1./TMath::Sqrt(errFin); 
      break;
    }
    // estimate error
    if (!(errFin>0)) {
      printf("failed to evaluate eta bin %d\n",ibe);
      continue;
    }
    //
    //    printf("#%3d Truncated W.Mean for %d out of %d selected bins: %+e +- %e\n",ibe,useCnt,cnt,valFin,errFin);
    //
    prj->SetBinContent(ibe,valFin);
    prj->SetBinError(ibe,errFin);
  }
  //
  return prj;
}

void KillBadBins(TH2* histo, double mn,double mx)
{
  MyGuard(2,"Killing bad (<%f or >%f) bins of %s",
	  mn, mx, histo->GetName());
  int nbx = histo->GetNbinsX();
  int nby = histo->GetNbinsY();
  for (int ix=1;ix<=nbx;ix++) {
    for (int iy=1;iy<=nby;iy++) {
      double vl = histo->GetBinContent(ix,iy);
      if (vl>mx || vl<mn) {
	histo->SetBinContent(ix,iy,0);
	histo->SetBinError(ix,iy,0);
      }
    }
  }
  histo->SetMinimum(mn);
  histo->SetMaximum(mx);
}

void PrintH(TH2* h, Int_t prec)
{
  Printf("Content of %s - %s", h->GetName(), h->GetTitle());
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    printf("%3d: ", i);
    for (Int_t j = 1; j <= h->GetNbinsY(); j++) {
      printf("%.*f+/-%.*f ",
	     prec, h->GetBinContent(i,j),
	     prec, h->GetBinError(i,j));
    }
    printf("\n");
  }
}
//____________________________________________________________________
void PrintH(TH1* h, Int_t prec)
{
  Printf("Content of %s - %s", h->GetName(), h->GetTitle());
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    if (h->GetBinContent(i) <= 1e-6) continue;
    printf("%3d (%+6.3f): %.*f+/-%.*f\n", i,
	   h->GetXaxis()->GetBinCenter(i),
	   prec, h->GetBinContent(i),
	   prec, h->GetBinError(i));
  }
}
