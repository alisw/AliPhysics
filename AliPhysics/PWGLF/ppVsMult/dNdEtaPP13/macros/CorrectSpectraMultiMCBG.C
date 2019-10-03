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
TString  useBgType    = "";
Int_t    useShapeType = kNormShapeDist;    // which distribution to use for bg normalization
Int_t    useMCLB      = 2;             // use Comb MC Labels as a template for Bg.
Double_t scaleBG      = 1.3;    // apply scaling to extracted bg for data
Int_t    useScaleType = kSclIntegral;//kSclWghMean;       // which type of tails normalization to use
Bool_t   useZbinWAv  = kFALSE;//kFALSE;
Bool_t   normToMB =  kFALSE;    // normalize with MB MC truth, rather than for selected events FIXME: WHAT IS THIS FOR???? only used to decide if one has to scale mc prim, and no other histo. Conflicting with ForceMB?
//Bool_t   normToMB =  kTRUE;    // normalize with MB MC truth, rather than for selected events
const double kEtaFitRange = 1.;//0.5;

Double_t minAlpha = 0.;
//Double_t maxAlpha = 2.5;
Double_t maxAlpha = 2500;

enum {kBitNormPerEvent=BIT(14)};
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
    kDataDist,kBgDist,kBgMCDist,kZvDist,kZvMCDistNS,kMCShift=20, kNHistos=kMCShift+kMCShift};

void    CorrectSpectraMultiMCBG(const char* flNameData, const char* flNameMC, const char *nameL = "clist", const char* unique="", Bool_t forceMCMB = kFALSE , int maxBins=20);
Bool_t  PrepareHistos(int bin, TList* lst, Bool_t isMC);
Bool_t  PrepareHistosHack(int bin, TList* lst, Bool_t isMC, Bool_t forceMCMB);
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

void PlotResults(Bool_t forceMCMB = kFALSE);
void PlotDNDEta(int bin, Bool_t forceMCMB = kFALSE);
void PlotAlphaBeta(int bin);
void PlotSpecies();

Int_t zbinMin = 0;
Int_t zbinMax = -1;

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
Bool_t creatDnDEtaCMacro = kFALSE;
Bool_t creatAlphaBetaCMacro = kFALSE;
Bool_t creatSpeciesCMacro = kFALSE;

void CorrectSpectraMultiMCBG(const char* flNameData, const char* flNameMC, const char *nameL, const char* uniqueNm, Bool_t forceMCMB ,int maxBin)
{
    
    printf("--->>> CorrectSpectraMultiMCBG(\"%s\", \"%s\", \"%s\", \"%s\", %d, %d)\n", flNameData, flNameMC, nameL, uniqueNm, forceMCMB, maxBin);
    
    //
    //gROOT->LoadMacro("SaveCanvas.C+");
    uniqueName = uniqueNm;
    listDt = LoadList(flNameData,"dt_", nameL);
    if (!listDt) return;
    if (forceMCMB && 0) { // R+HACK FIXME: MF is skipping this, as roberto what this meams
      TString str = nameL;
      if (str.EndsWith("_L")) str = "clist_MB_L";
      else if (str.EndsWith("_LL")) str = "clist_MB_LL";
      else if (str.EndsWith("_R")) str = "clist_MB_R";
      else if (str.EndsWith("_RR")) str = "clist_MB_RR";
      else if (str.EndsWith("_W")) str = "clist_MB_W";
      else if (str.EndsWith("_WW")) str = "clist_MB_WW";
      else str = "clist_MB";
      printf("--->>> STR = %s\n", str.Data());
      listMC = LoadList(flNameMC,"mc_", str.Data()); // R+HACK
    }
    else
      listMC = LoadList(flNameMC,"mc_", nameL);
    if (!listMC) return;
    //
    resArr.Clear();
    //
    TH1* hstat = (TH1*)FindObject(-1,kHStatName,listDt,kFALSE);
    if (!hstat) return;
    //TH1* hstat = (TH1*)FindObject(-1,kHStatName,listMC,kFALSE);
    //
    int nbstat = hstat->GetNbinsX();
    nCentBins = (nbstat - kBinEntries)/kEntriesPerBin;
    nCentBins = nCentBins>maxBin ? maxBin : nCentBins;
    printf("%d bins will be processed\n",nCentBins);
    if (nCentBins<1) return;
    double myMergeFactor = hstat->GetBinContent(kOneUnit);
    printf("Detected %f mergings\n",myMergeFactor);
    if (myMergeFactor<1) myMergeFactor = 1.;
    //
    dNdEta.Set(nCentBins);
    dNdEtaErr.Set(nCentBins);
    //
    kdPhiSgCut  = hstat->GetBinContent(kDPiSCut)/myMergeFactor;
    kWDistSgCut = hstat->GetBinContent(kNStdCut)/myMergeFactor;
    printf("Signal cuts used: dPhiS: %f WDist:%f\n",kdPhiSgCut,kWDistSgCut);
    //
    for (int ib=0;ib<nCentBins;ib++) {
      if (!PrepareHistosHack(ib,listMC,kTRUE,forceMCMB))  continue; // R+HACK
      if (!PrepareHistosHack(ib,listDt,kFALSE,forceMCMB)) continue;
      // // FIXME: should this be Hack or not hack?
      // if (!PrepareHistos(ib,listMC,kTRUE))  continue; // 
      // if (!PrepareHistos(ib,listDt,kFALSE)) continue;
      ProcessHistos(ib);
      //
    }
    //
    sprintf(outStr,"CutEta%.1f_%.1f_Zv%.1f_%.1f_bg%s_Shape_%s_mcLB%d_cutSig%.1f_cutBg%.1f",
            hstat->GetBinContent(kEtaMin)/myMergeFactor,
            hstat->GetBinContent(kEtaMax)/myMergeFactor,
            hstat->GetBinContent(kZVMin)/myMergeFactor,
            hstat->GetBinContent(kZVMax)/myMergeFactor,
            useBgType.Data(),
            useShapeType==kNormShapeDist ? "wdst":"dphi",
            useMCLB,
            useShapeType==kNormShapeDist ? kWDistSgCut:kdPhiSgCut,
            useShapeType==kNormShapeDist ? kWDistBgTailMin:kdPhiBgTailMin);
    //
    PlotResults(forceMCMB);
    //
    printf("Final Results:\n");
    printf("dNdEta:    ");
    for (int i=nCentBins;i--;) printf("%.2f,",dNdEta[i]); printf("\n");
    printf("dNdEtaErr: ");
    for (int i=nCentBins;i--;) printf("%.2f,",dNdEtaErr[i]); printf("\n");
    //
    printf("--->>> resDir = %s\n", resDir);
    TString rtnm1 = resDir; rtnm1 += "/"; rtnm1 += uniqueName;
    rtnm1 += "_"; rtnm1 += nCentBins; rtnm1+= "bins_";
    rtnm1 += outStr; rtnm1 += ".root";
    TFile* flRes = TFile::Open(rtnm1.Data(),"recreate");
    flRes->WriteObject(&resDnDeta, "TObjArray", "kSingleKey");
    flRes->WriteObject(&resArr, "TObjArrayAux", "kSingleKey");
    flRes->Close();
    delete flRes;
    printf("Stored result in %s\n",rtnm1.Data());
    // R+HACK
    printf("--->>> resDir = %s\n", resDir);
    rtnm1 = resDir; rtnm1 += "/"; rtnm1 += uniqueName;
    rtnm1 += ".root";
    flRes = TFile::Open(rtnm1.Data(),"recreate");
    flRes->WriteObject(&resDnDeta, "TObjArray", "kSingleKey");
    flRes->WriteObject(&resArr, "TObjArrayAux", "kSingleKey");
    flRes->Close();
    delete flRes;
    printf("Stored result in %s\n",rtnm1.Data());
    //
}

//_____________________________________________________________________
Bool_t PrepareHistos(int bin, TList* lst, Bool_t isMC)
{

  printf("--->>> PrepareHistos(%d, %p, %d)\n", bin, lst, isMC);
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
    double nrmBin = hstat->GetBinContent(kBinEntries + kEvProcData+bin*kEntriesPerBin);
    if (nrmBin<1) {
        printf("Norm is 0 for bin %d\n",bin);
        //    exit(1);
        return kFALSE;
    }
    
    const char* zeCut = "ZvEtaCutT";
    TObjArray* res = &resArr;
    int shift = bin*kNHistos + (isMC ? kMCShift : 0);
    //
    // histo for "data" Z vs Eta with default signal cut
    TH2* hRawDtCut = (TH2*) FindObject(bin,HName("Data",zeCut),lst);

    if (!hRawDtCut) return kFALSE;
    sprintf(buffn,"bin%d_%s_RawWithCut",bin,isMC ? "mc":"dt");
    sprintf(bufft,"bin%d %s Raw Data with cut on tracklets",bin,isMC ? "mc":"dt");
    hRawDtCut = (TH2*)hRawDtCut->Clone(buffn);
    hRawDtCut->SetTitle(bufft);
    res->AddAtAndExpand(hRawDtCut, kRawDtCut+shift);
    //
    int nbEta = hRawDtCut->GetXaxis()->GetNbins();
    int nbZV  = hRawDtCut->GetYaxis()->GetNbins();
    //
    // Zv nonuniformity histo
    TH2* hzv2 = (TH2*) FindObject(-1,"zv",lst, kFALSE);
    TH1* hz   = hzv2->ProjectionX(Form("zv%d_%s",bin,isMC ? "mc":"dt"),bin+1,bin+1,"e");
    hz->Scale(1./nrmBin);
    res->AddAtAndExpand(hz, kZvDist+shift);
    //
    // "Data - Est.Bg" histo with cut on tails where we look for signal
    sprintf(buffn,"bin%d_%s_SignalWithCut",bin,isMC ? "mc":"dt");
    sprintf(bufft,"bin%d %s Signal (raw-bg) with cut on tracklets",bin,isMC ? "mc":"dt");
    TH2* hSignalEst = (TH2F*)hRawDtCut->Clone(buffn);
    hSignalEst->SetTitle(bufft);
    res->AddAtAndExpand(hSignalEst, kSignalEst+shift);
    //
    // "Data - MC.Bg" histo with cut on tails where we look for signal
    TH2* hSignalEstMC = 0;
    if (isMC) {
        sprintf(buffn,"bin%d_%s_SignalWithCut_bgMCLabels",bin,isMC ? "mc":"dt");
        sprintf(bufft,"bin%d %s Signal (raw-bg_MCLabels) with cut on tracklets",bin,isMC ? "mc":"dt");
        hSignalEstMC = (TH2F*)hRawDtCut->Clone(buffn);
        hSignalEstMC->SetTitle(bufft);
        res->AddAtAndExpand(hSignalEstMC, kSignalEstMC+shift);
    }
    //
    // Estimated background in the cut range
    TH2* hBgEst = 0;//(TH2*) FindObject(bin,HName(useBgType.Data(),zeCut),lst);
    // special feature: use MC Labels bg as a shape instead of generated bg
    if (useMCLB/* && !isMC*/) {
        sprintf(buffn,"bin%d_%s_BgEst",bin,isMC ? "mc":"dt");
        sprintf(bufft,"bin%d %s Estimated Bg",bin,isMC ? "mc":"dt");
        TString nm  = buffn;   nm  += "_MCLB";
        TString tit = bufft;   tit += "_MCLB";
        TH2* hBMCLB = (TH2*) FindObject(bin,HName("Comb",zeCut),listMC);
        if (!hBMCLB) return kFALSE;
        hBMCLB = (TH2F*)hBMCLB->Clone(nm.Data());
        hBMCLB->SetTitle(tit.Data());
        delete hBgEst;
        hBgEst = hBMCLB;
        res->AddAtAndExpand(hBgEst,kBgEst +shift);
    }
    //
    // 1-beta for "data" = (Data_cut - Bg_cut) / Data_cut
    sprintf(buffn,"bin%d_%s_1mBeta",bin,isMC ? "mc":"dt");
    sprintf(bufft,"bin%d %s 1-#beta with estimated bg",bin,isMC ? "mc":"dt");
    TH2* h1mBeta = (TH2*)hBgEst->Clone(buffn);
    h1mBeta->SetTitle(bufft);
    h1mBeta->Reset();
    res->AddAtAndExpand(h1mBeta, k1MBeta+shift);
    //
    // If MC labels info is provided, prepare 1-beta with MC bg
    TH2* h1mBetaMC = 0;  // 1-beta for MC with bg from labels
    TH2* h1mBetaMCscl = 0; // scaled by ad hoc factor
    TH2* hBgMC = 0;      // bg from MC labels
    if (isMC) {
        sprintf(buffn,"bin%d_%s_BgMC",bin,isMC ? "mc":"dt");
        sprintf(bufft,"bin%d %s Bg from MC labels",bin,isMC ? "mc":"dt");
        hBgMC = (TH2*) FindObject(bin,HName("Comb",zeCut),listMC);
        if (!hBgMC) return kFALSE;
        hBgMC = (TH2F*)hBgMC->Clone(buffn);
        hBgMC->SetTitle(bufft);
        res->AddAtAndExpand(hBgMC, kBgMC+shift);
        //
        sprintf(buffn,"bin%d_%s_h1mBetaMC",bin,isMC ? "mc":"dt");
        sprintf(bufft,"bin%d %s 1-#beta with bg. from MC labels",bin,isMC ? "mc":"dt");
        h1mBetaMC = (TH2F*) hBgMC->Clone(buffn);
        h1mBetaMC->SetTitle(bufft);
        //    printf(">>Here1\n");
        h1mBetaMC->Divide(hRawDtCut);
        h1mBetaMCscl = (TH2*)h1mBetaMC->Clone(Form("%s_scl",h1mBetaMC->GetName()));
        //    printf("<<Here1\n");
        for (int ib0=1;ib0<=nbEta;ib0++)
            for (int ib1=1;ib1<=nbZV;ib1++) {
                h1mBetaMCscl->SetBinContent(ib0,ib1, 1.- scaleBG*h1mBetaMC->GetBinContent(ib0,ib1));
                h1mBetaMC->SetBinContent(ib0,ib1, 1.- h1mBetaMC->GetBinContent(ib0,ib1));
            }
        //

        res->AddAtAndExpand(h1mBetaMC, k1MBetaMC+shift);
        res->AddAtAndExpand(h1mBetaMCscl, k1MBetaMCscl+shift);
        hSignalEstMC->Multiply(h1mBetaMC);
        //    hSignalEstMC->Add(hBgMC,-1);
    }
    //
    // uncut w.distance or dphi distribution for data
    TH1* hDstDt = (TH1*) FindObject(bin,HName("Data",useShapeType==kNormShapeDist ? "WDist":"DPhiS"),lst);
    if (!hDstDt) return kFALSE;
    sprintf(buffn,"bin%d_%s_DistRawData",bin,isMC ? "mc":"dt");
    sprintf(bufft,"bin%d %s Raw Distance for Data",bin,isMC ? "mc":"dt");
    hDstDt = (TH1*) hDstDt->Clone(buffn);
    hDstDt->SetTitle(bufft);
    double nrmDst,dumErr = 0;
    Integrate(hDstDt, cutBgMin,cutBgMax, nrmDst, dumErr);
    if (nrmDst<1e-10) nrmDst = 1.; // RS resque
    hDstDt->Scale(1./nrmDst);
    res->AddAtAndExpand(hDstDt, kDataDist+shift);
    //
    // uncut w.distance or dphi distribution for generated bg
    TH1* hDstBg = (TH1*) 0;//FindObject(bin,HName(useBgType.Data(),useShapeType==kNormShapeDist ? "WDist":"DPhiS"),lst);
    if (!hDstBg) {} //return kFALSE;
    else {
        sprintf(buffn,"bin%d_%s_DistRawGenBgNorm",bin,isMC ? "mc":"dt");
        sprintf(bufft,"bin%d %s Raw Distance for Gen.Bg. Normalized to data",bin,isMC ? "mc":"dt");
        hDstBg = (TH1*) hDstBg->Clone(buffn);
        hDstBg->SetTitle(bufft);
        hDstBg->Scale(1./nrmDst);
    }
    //  res->AddAtAndExpand(hDstBg, kBgDist+shift);
    //
    // uncut w.distance or dphi distribution for comb. MC labels
    TH1* hDstBgMC = 0;
    if (isMC) {
        hDstBgMC = (TH1*) FindObject(bin,HName("Comb",useShapeType==kNormShapeDist ? "WDist":"DPhiS"),lst);
        if (!hDstBgMC) return kFALSE;
        sprintf(buffn,"bin%d_%s_DistBgMC",bin,isMC ? "mc":"dt");
        sprintf(bufft,"bin%d %s Bg. Distance from MC labels",bin,isMC ? "mc":"dt");
        hDstBgMC = (TH1*) hDstBgMC->Clone(buffn);
        hDstBgMC->SetTitle(bufft);
        hDstBgMC->Scale(1./nrmDst);
        res->AddAtAndExpand(hDstBgMC, kBgMCDist+shift);
    }
    if (!hDstBg && useMCLB) {
        sprintf(buffn,"bin%d_%s_DistRawGenBgNorm",bin,isMC ? "mc":"dt");
        sprintf(bufft,"bin%d %s Raw Distance for Gen.Bg. Normalized to data",bin,isMC ? "mc":"dt");
        hDstBgMC = (TH1*)res->At(kBgMCDist+bin*kNHistos + kMCShift);
        printf("Subst. bin %d hDstBg by hDstBgMC @ %p\n",bin,hDstBgMC);
        hDstBg = (TH1*)hDstBgMC->Clone(buffn);
        hDstBg->SetTitle(bufft);
    }
    //
    // fill 1-beta matrix
    double scl,sclE;
    hDstBg = NormalizeBg(hDstDt,hDstBg,scl,sclE);  // get rescaling factor for bg. from tails comparison
    res->AddAtAndExpand(hDstBg, kBgDist+shift);
    double bgVal,bgErr;
    double dtVal,dtErr;
    // integral in the range where we look for signal
    Integrate(hDstBg, cutSgMin, cutSgMax, bgVal, bgErr);
    Integrate(hDstDt, cutSgMin, cutSgMax, dtVal, dtErr);
    double sclb,sclbErr;
    GetRatE(bgVal,bgErr, dtVal, dtErr,sclb,sclbErr);
    //  hDstBg->Scale(1./nrmDst);
    //
    // finalize estimated bg and signal matrices
    hBgEst->Scale(scl);
    //
    if (useMCLB==2) {
        printf("Using MC Combinatorials fraction for real data bg.estimate\n");
        TH2* bet1m = (TH2*)res->At( (isMC ? k1MBetaMC : k1MBetaMCscl) + shift + (isMC ? 0 :kMCShift) );
        //TH2* bet1m = (TH2*)res->At( k1MBetaMC  + shift + (isMC ? 0 :kMCShift) );
        //    printf("Here %d (%d): %p\n",isMC,(isMC ? k1MBetaMC : k1MBetaMCscl) + shift + (isMC ? 0 :kMCShift), bet1m); bet1m->Print();
        hSignalEst->Multiply(bet1m);
    }
    else {
        hSignalEst->Add(hBgEst,-1);
    }
    //
    // finalize 1-beta
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
    //
    if (isMC) {
        // Zv nonuniformity histo: MCGen before cuts
        TH2* hzv2ns = (TH2*) FindObject(-1,"zvMCNoPS",lst, kFALSE);
        TH1* hzns   = hzv2ns->ProjectionX(Form("zvMCNS%d",bin),bin+1,bin+1,"e");
        //
        double myMergeFactor = hstat->GetBinContent(kOneUnit);
        double zmn = hstat->GetBinContent(kZVMin)/myMergeFactor;
        double zmx = hstat->GetBinContent(kZVMax)/myMergeFactor;
        int zbmn = hzns->FindBin(zmn+1e-3);
        int zbmx = hzns->FindBin(zmx-1e-3);
        double nEvTot = hzns->Integral(zbmn,zbmx);
        printf("Scale MC Truth By %9.1f\n",nEvTot);
        if (nEvTot<1) exit(1);
        hzns->Scale(1./nEvTot);
        /*
         // here there could be some over/under flows, add them to first/last used bins correspondingly
         double under = hzns->Integral(0,zbmn-1);
         double over  = hzns->Integral(zbmx+1,hzns->GetNbinsX()+1);
         hzns->SetBinContent(zbmn, hzns->GetBinContent(zbmn)+under);
         hzns->SetBinContent(zbmx, hzns->GetBinContent(zbmx)+over);
         printf("Add %e and %e of under/overflows to extreme bins\n",under,over);
         */
        res->AddAtAndExpand(hzns, kZvMCDistNS+shift);
        //
        TH2F* mcPrim = 0;
        sprintf(buffn,"bin%d_zvEtaPrimMC",bin);
        sprintf(bufft,"bin%d MC True signal Zv vs Eta",bin);
        //
        // prepare MC primary signal histo
        if (normToMB) {
            mcPrim = (TH2F*)FindObject(bin,"zvEtaPrimMC", lst, kFALSE ); // don't normalise to NSel
            mcPrim = (TH2F*) mcPrim->Clone(buffn);
            mcPrim->Scale(1./nEvTot);
        }
        else {
            mcPrim = (TH2F*)FindObject(bin,"zvrecEtaPrimMCSel", lst, kTRUE ); // don't normalise to NSel
        }
        mcPrim->SetTitle(bufft);
        res->AddAtAndExpand(mcPrim, kMCPrim + shift);
        //
    }
    //

    return kTRUE;
}

//_____________________________________________________________________
Bool_t PrepareHistosHack(int bin, TList* lst, Bool_t isMC, Bool_t forceMCMB)
{

  // MF 2015-12-08 This was created by Roberto. If I understand it correctly, if forces the usage of the MB bin for corrections, by setting the bin_ variables
  int bin_ = bin;
  int bin__ = bin;
  if (forceMCMB)
    bin_ = 0;
  if (forceMCMB && isMC)
    bin__ = 0;
  
  printf("--->>> PrepareHistosHack(%d, %p, %d, %d)\n", bin, lst, isMC, forceMCMB);
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
    double nrmBin = hstat->GetBinContent(kBinEntries + kEvProcData+bin__*kEntriesPerBin);
    if (nrmBin<1) {
        printf("Norm is 0 for bin %d\n",bin);
        //    exit(1);
        return kFALSE;
    }
    
    const char* zeCut = "ZvEtaCutT";
    TObjArray* res = &resArr;
    int shift = bin*kNHistos + (isMC ? kMCShift : 0);
    int shift_ = bin_*kNHistos + (isMC ? kMCShift : 0);
    int shift__ = bin__*kNHistos + (isMC ? kMCShift : 0);
    //
    // histo for "data" Z vs Eta with default signal cut
    TH2* hRawDtCut = (TH2*) FindObject(bin__,HName("Data",zeCut),lst);
    if (!hRawDtCut) return kFALSE;
    sprintf(buffn,"bin%d_%s_RawWithCut",bin,isMC ? "mc":"dt");
    sprintf(bufft,"bin%d %s Raw Data with cut on tracklets",bin,isMC ? "mc":"dt");
    hRawDtCut = (TH2*)hRawDtCut->Clone(buffn);
    hRawDtCut->SetTitle(bufft);
    res->AddAtAndExpand(hRawDtCut, kRawDtCut+shift);
    //
    int nbEta = hRawDtCut->GetXaxis()->GetNbins();
    int nbZV  = hRawDtCut->GetYaxis()->GetNbins();
    //
    // Zv nonuniformity histo
    TH2* hzv2 = (TH2*) FindObject(-1,"zv",lst, kFALSE);
    TH1* hz   = hzv2->ProjectionX(Form("zv%d_%s",bin,isMC ? "mc":"dt"),bin__+1,bin__+1,"e");
    hz->Scale(1./nrmBin);
    res->AddAtAndExpand(hz, kZvDist+shift);
    //
    // "Data - Est.Bg" histo with cut on tails where we look for signal
    sprintf(buffn,"bin%d_%s_SignalWithCut",bin,isMC ? "mc":"dt");
    sprintf(bufft,"bin%d %s Signal (raw-bg) with cut on tracklets",bin,isMC ? "mc":"dt");
    TH2* hSignalEst = (TH2F*)hRawDtCut->Clone(buffn);
    hSignalEst->SetTitle(bufft);
    res->AddAtAndExpand(hSignalEst, kSignalEst+shift);
    //
    // "Data - MC.Bg" histo with cut on tails where we look for signal
    TH2* hSignalEstMC = 0;
    if (isMC) {
        sprintf(buffn,"bin%d_%s_SignalWithCut_bgMCLabels",bin,isMC ? "mc":"dt");
        sprintf(bufft,"bin%d %s Signal (raw-bg_MCLabels) with cut on tracklets",bin,isMC ? "mc":"dt");
        hSignalEstMC = (TH2F*)hRawDtCut->Clone(buffn);
        hSignalEstMC->SetTitle(bufft);
        res->AddAtAndExpand(hSignalEstMC, kSignalEstMC+shift);
    }
    //
    // Estimated background in the cut range
    TH2* hBgEst = 0;//(TH2*) FindObject(bin,HName(useBgType.Data(),zeCut),lst);
    // special feature: use MC Labels bg as a shape instead of generated bg
    if (useMCLB/* && !isMC*/) {
        sprintf(buffn,"bin%d_%s_BgEst",bin,isMC ? "mc":"dt");
        sprintf(bufft,"bin%d %s Estimated Bg",bin,isMC ? "mc":"dt");
        TString nm  = buffn;   nm  += "_MCLB";
        TString tit = bufft;   tit += "_MCLB";
        TH2* hBMCLB = (TH2*) FindObject(bin_,HName("Comb",zeCut),listMC);
        if (!hBMCLB) return kFALSE;
        hBMCLB = (TH2F*)hBMCLB->Clone(nm.Data());
        hBMCLB->SetTitle(tit.Data());
        delete hBgEst;
        hBgEst = hBMCLB;
        res->AddAtAndExpand(hBgEst,kBgEst +shift);
    }
    
    //
    // 1-beta for "data" = (Data_cut - Bg_cut) / Data_cut
    sprintf(buffn,"bin%d_%s_1mBeta",bin,isMC ? "mc":"dt");
    sprintf(bufft,"bin%d %s 1-#beta with estimated bg",bin,isMC ? "mc":"dt");
    TH2* h1mBeta = (TH2*)hBgEst->Clone(buffn);
    h1mBeta->SetTitle(bufft);
    h1mBeta->Reset();
    res->AddAtAndExpand(h1mBeta, k1MBeta+shift);
    //
    // If MC labels info is provided, prepare 1-beta with MC bg
    TH2* h1mBetaMC = 0;  // 1-beta for MC with bg from labels
    TH2* h1mBetaMCscl = 0; // scaled by ad hoc factor
    TH2* hBgMC = 0;      // bg from MC labels
    if (isMC) {
        sprintf(buffn,"bin%d_%s_BgMC",bin,isMC ? "mc":"dt");
        sprintf(bufft,"bin%d %s Bg from MC labels",bin,isMC ? "mc":"dt");
        hBgMC = (TH2*) FindObject(bin__,HName("Comb",zeCut),listMC);
        if (!hBgMC) return kFALSE;
        hBgMC = (TH2F*)hBgMC->Clone(buffn);
        hBgMC->SetTitle(bufft);
        res->AddAtAndExpand(hBgMC, kBgMC+shift);
        //
        sprintf(buffn,"bin%d_%s_h1mBetaMC",bin,isMC ? "mc":"dt");
        sprintf(bufft,"bin%d %s 1-#beta with bg. from MC labels",bin,isMC ? "mc":"dt");
        h1mBetaMC = (TH2F*) hBgMC->Clone(buffn);
        h1mBetaMC->SetTitle(bufft);
        //    printf(">>Here1\n");
        h1mBetaMC->Divide(hRawDtCut);
        h1mBetaMCscl = (TH2*)h1mBetaMC->Clone(Form("%s_scl",h1mBetaMC->GetName()));
        //    printf("<<Here1\n");
        for (int ib0=1;ib0<=nbEta;ib0++)
            for (int ib1=1;ib1<=nbZV;ib1++) {
                h1mBetaMCscl->SetBinContent(ib0,ib1, 1.- scaleBG*h1mBetaMC->GetBinContent(ib0,ib1));
                h1mBetaMC->SetBinContent(ib0,ib1, 1.- h1mBetaMC->GetBinContent(ib0,ib1));
            }
        //

        res->AddAtAndExpand(h1mBetaMC, k1MBetaMC+shift);
        res->AddAtAndExpand(h1mBetaMCscl, k1MBetaMCscl+shift);
        hSignalEstMC->Multiply(h1mBetaMC);
        //    hSignalEstMC->Add(hBgMC,-1);
    }
    //
    // uncut w.distance or dphi distribution for data
    TH1* hDstDt = (TH1*) FindObject(bin__,HName("Data",useShapeType==kNormShapeDist ? "WDist":"DPhiS"),lst);
    if (!hDstDt) return kFALSE;
    sprintf(buffn,"bin%d_%s_DistRawData",bin,isMC ? "mc":"dt");
    sprintf(bufft,"bin%d %s Raw Distance for Data",bin,isMC ? "mc":"dt");
    hDstDt = (TH1*) hDstDt->Clone(buffn);
    hDstDt->SetTitle(bufft);
    double nrmDst,dumErr = 0;
    Integrate(hDstDt, cutBgMin,cutBgMax, nrmDst, dumErr);
    if (nrmDst<1e-10) nrmDst = 1.; // RS resque
    hDstDt->Scale(1./nrmDst);
    res->AddAtAndExpand(hDstDt, kDataDist+shift);
    //
    // uncut w.distance or dphi distribution for generated bg
    TH1* hDstBg = (TH1*) 0;//FindObject(bin,HName(useBgType.Data(),useShapeType==kNormShapeDist ? "WDist":"DPhiS"),lst);
    if (!hDstBg) {} //return kFALSE;
    else {
        sprintf(buffn,"bin%d_%s_DistRawGenBgNorm",bin,isMC ? "mc":"dt");
        sprintf(bufft,"bin%d %s Raw Distance for Gen.Bg. Normalized to data",bin,isMC ? "mc":"dt");
        hDstBg = (TH1*) hDstBg->Clone(buffn);
        hDstBg->SetTitle(bufft);
        hDstBg->Scale(1./nrmDst);
    }
    //  res->AddAtAndExpand(hDstBg, kBgDist+shift);
    //
    // uncut w.distance or dphi distribution for comb. MC labels
    TH1* hDstBgMC = 0;
    if (isMC) {
        hDstBgMC = (TH1*) FindObject(bin__,HName("Comb",useShapeType==kNormShapeDist ? "WDist":"DPhiS"),lst);
        if (!hDstBgMC) return kFALSE;
        sprintf(buffn,"bin%d_%s_DistBgMC",bin,isMC ? "mc":"dt");
        sprintf(bufft,"bin%d %s Bg. Distance from MC labels",bin,isMC ? "mc":"dt");
        hDstBgMC = (TH1*) hDstBgMC->Clone(buffn);
        hDstBgMC->SetTitle(bufft);
        hDstBgMC->Scale(1./nrmDst);
        res->AddAtAndExpand(hDstBgMC, kBgMCDist+shift);
    }
    if (!hDstBg && useMCLB) {
        sprintf(buffn,"bin%d_%s_DistRawGenBgNorm",bin,isMC ? "mc":"dt");
        sprintf(bufft,"bin%d %s Raw Distance for Gen.Bg. Normalized to data",bin,isMC ? "mc":"dt");
        hDstBgMC = (TH1*)res->At(kBgMCDist+bin_*kNHistos + kMCShift);
        printf("Subst. bin %d hDstBg by hDstBgMC @ %p\n",bin,hDstBgMC);
        hDstBg = (TH1*)hDstBgMC->Clone(buffn);
        hDstBg->SetTitle(bufft);
    }
    //
    // fill 1-beta matrix
    double scl,sclE;
    hDstBg = NormalizeBg(hDstDt,hDstBg,scl,sclE);  // get rescaling factor for bg. from tails comparison
    res->AddAtAndExpand(hDstBg, kBgDist+shift);
    double bgVal,bgErr;
    double dtVal,dtErr;
    // integral in the range where we look for signal
    Integrate(hDstBg, cutSgMin, cutSgMax, bgVal, bgErr);
    Integrate(hDstDt, cutSgMin, cutSgMax, dtVal, dtErr);
    double sclb,sclbErr;
    GetRatE(bgVal,bgErr, dtVal, dtErr,sclb,sclbErr);
    //  hDstBg->Scale(1./nrmDst);
    //
    // finalize estimated bg and signal matrices
    
    hBgEst->Scale(scl);
    //
    if (useMCLB==2) {
        printf("Using MC Combinatorials fraction for real data bg.estimate\n");
        TH2* bet1m = (TH2*)res->At( (isMC ? k1MBetaMC : k1MBetaMCscl) + shift_ + (isMC ? 0 :kMCShift) );
        //TH2* bet1m = (TH2*)res->At( k1MBetaMC  + shift + (isMC ? 0 :kMCShift) );
        //    printf("Here %d (%d): %p\n",isMC,(isMC ? k1MBetaMC : k1MBetaMCscl) + shift + (isMC ? 0 :kMCShift), bet1m); bet1m->Print();
        hSignalEst->Multiply(bet1m);
    }
    else {
        hSignalEst->Add(hBgEst,-1);
    }
    //
    // finalize 1-beta
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

    //
    if (isMC) {
        // Zv nonuniformity histo: MCGen before cuts
        TH2* hzv2ns = (TH2*) FindObject(-1,"zvMCNoPS",lst, kFALSE);
        TH1* hzns   = hzv2ns->ProjectionX(Form("zvMCNS%d",bin),bin__+1,bin__+1,"e");
        //
        double myMergeFactor = hstat->GetBinContent(kOneUnit);
        double zmn = hstat->GetBinContent(kZVMin)/myMergeFactor;
        double zmx = hstat->GetBinContent(kZVMax)/myMergeFactor;
        int zbmn = hzns->FindBin(zmn+1e-3);
        int zbmx = hzns->FindBin(zmx-1e-3);
        double nEvTot = hzns->Integral(zbmn,zbmx);
        printf("Scale MC Truth By %9.1f\n",nEvTot);
        if (nEvTot<1) exit(1);
        hzns->Scale(1./nEvTot);
        /*
         // here there could be some over/under flows, add them to first/last used bins correspondingly
         double under = hzns->Integral(0,zbmn-1);
         double over  = hzns->Integral(zbmx+1,hzns->GetNbinsX()+1);
         hzns->SetBinContent(zbmn, hzns->GetBinContent(zbmn)+under);
         hzns->SetBinContent(zbmx, hzns->GetBinContent(zbmx)+over);
         printf("Add %e and %e of under/overflows to extreme bins\n",under,over);
         */
        res->AddAtAndExpand(hzns, kZvMCDistNS+shift);
        //
        TH2F* mcPrim = 0;
        sprintf(buffn,"bin%d_zvEtaPrimMC",bin);
        sprintf(bufft,"bin%d MC True signal Zv vs Eta",bin);
        //
        // prepare MC primary signal histo
        if (normToMB) {
            mcPrim = (TH2F*)FindObject(bin__,"zvEtaPrimMC", lst, kFALSE ); // don't normalise to NSel
            mcPrim = (TH2F*) mcPrim->Clone(buffn);
            mcPrim->Scale(1./nEvTot);
        }
        else {
            mcPrim = (TH2F*)FindObject(bin__,"zvrecEtaPrimMCSel", lst, kTRUE ); // don't normalise to NSel
        }
        mcPrim->SetTitle(bufft);
        res->AddAtAndExpand(mcPrim, kMCPrim + shift);
        //
    }
    //
    return kTRUE;
}

//_____________________________________________________________________
void ProcessHistos(int bin)
{
    printf("--->>> ProcessHistos(%d)\n", bin);
    
    //
    int shift = bin*kNHistos;
    //
    TString prefN = "bin"; prefN += Form("%02d", bin); prefN+="_";
    TString prefT = "bin"; prefT += Form("%02d", bin); prefT+=" ";
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

    // build alpha matrix
    TH2* halp = (TH2*)res->At(shift + kMCShift + kMCPrim);
    printf("XXXXX %p XXXXX\n", halp);
    halp->Print();

    halp = (TH2*) halp->Clone(prefN+"Alpha");
    halp->SetTitle(prefN+"#alpha");
    printf("Halp: %p %d %d\n",halp,halp->GetNbinsX(),halp->GetNbinsY()); halp->Print();
    TH2* xxx = (TH2*)res->At(shift + kMCShift + k1MBeta);
    printf("other: %p %d %d\n",xxx,xxx->GetNbinsX(),xxx->GetNbinsY()); xxx->Print();
    halp->Divide( (TH2*)res->At(shift + kMCShift + k1MBeta) );
    halp->Divide( (TH2*)res->At(shift + kMCShift + kRawDtCut) );
    res->AddAtAndExpand(halp, shift + kAlpha);
    //

    KillBadBins(halp,minAlpha,maxAlpha);
    //

    // build alpha matrix with MC labels bg
    TH2* halpMC = (TH2*)res->At(shift + kMCShift + kMCPrim);
    halpMC = (TH2*) halpMC->Clone(prefN + "AlphaMC");
    halpMC->SetTitle(prefT + "#alpha MC labels");
    halpMC->Divide( (TH2*)res->At(shift + kMCShift + k1MBetaMC) );
    halpMC->Divide( (TH2*)res->At(shift + kMCShift + kRawDtCut) );
    res->AddAtAndExpand(halpMC, shift + kAlphaMC);
    //
    KillBadBins(halpMC,minAlpha,maxAlpha);
    //
    // build corrected signal
    TH2* hsigCorr = (TH2*)res->At(shift + kSignalEst);
    hsigCorr = (TH2*) hsigCorr->Clone(prefN + "SignalEstCorr");
    hsigCorr->SetTitle(prefT + "Corrected Signal");
    hsigCorr->Multiply( useMCLB==2 ? halpMC : halp );
    // hsigCorr->Multiply( halp );
    res->AddAtAndExpand(hsigCorr, shift + kSigCorr);
    //
    
    TH2* tmpSC = (TH2*)hsigCorr->Clone("tmpSC");
    TH1* hsigCorrX = 0;
    if (useZbinWAv) {
        CorrectForZV(tmpSC,hZV);
        hsigCorrX = ProjectWghMean(tmpSC,"DataCorrSignalX");
    }
    else {
      hsigCorrX = ProjNorm(tmpSC,hZV,"DataCorrSignalX",zbinMin,zbinMax); //R+HACK
    }
    //
    delete tmpSC;
    hsigCorrX->Scale(1./hsigCorrX->GetBinWidth(1));
    TF1* pl0 = new TF1("pl0","pol0");
    pl0->SetParameter(0,hsigCorr->GetMinimum());
    hsigCorrX->Fit(pl0,"q0","",-kEtaFitRange,kEtaFitRange);
    double fval = pl0->GetParameter(0);
    double ferr = pl0->GetParError(0);
    delete hsigCorrX;
    dNdEta[bin]    = fval;
    dNdEtaErr[bin] = ferr;
    printf("Bin %d: dN/d#eta_{|#eta|<0.5} = %.2f  %.2f\n",bin, fval,ferr);
    //
}

void PlotResults(Bool_t forceMCMB)
{
    
    printf("--->>> PlotResults()\n");
    
    TString psnm1 = figDir; psnm1 += "/"; psnm1 += uniqueName;
    psnm1 += "_"; psnm1 += nCentBins; psnm1+= "bins_";
    psnm1 += outStr; psnm1 += ".ps";
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
    canvFin->Print(psnm0.Data());
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
    gPad->SetGrid(1,1);
    //    resArr->AddAtAndExpand
    //
    grp->SaveAs("result.C");
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
    //
    canvFin->cd(0);
    canvFin->Print(psnm1.Data());
    //
    const TArrayD &binArr = *hbins->GetXaxis()->GetXbins();
    //
    for (int i=0;i<nCentBins;i++) {
        //
        sprintf(outTitle,"%s, %d<C_%s<%d, %.1f<#eta<%.1f, %.1f<Z_{V}<%.1f,  Bg.:%s, UseMCLB=%d, CutVar:%s, |sig|<%.2f, %.2f<|bg.nrm|<%.2f",
                uniqueName.Data(),
                (int)binArr[i],
                hstat->GetXaxis()->GetBinLabel(kCentVar),
                (int)binArr[i+1],
                hstat->GetBinContent(kEtaMin)/myMergeFactor,
                hstat->GetBinContent(kEtaMax)/myMergeFactor,
                hstat->GetBinContent(kZVMin)/myMergeFactor,
                hstat->GetBinContent(kZVMax)/myMergeFactor,
                useBgType.Data(),
                useMCLB,
                useShapeType==kNormShapeDist ? "#Delta":"#Delta#varphi-#delta_{#varphi}",
                useShapeType==kNormShapeDist ? kWDistSgCut:kdPhiSgCut,
                useShapeType==kNormShapeDist ? kWDistBgTailMin : kdPhiBgTailMin,
                useShapeType==kNormShapeDist ? kWDistBgTailMax : kdPhiBgTailMax
                );
        //
        PlotDNDEta(i, forceMCMB);
        canvFin->Print(psnm1.Data());
        PlotAlphaBeta(i);
        canvFin->Print(psnm1.Data());
    }
    PlotSpecies();
    canvFin->Print(psnm1.Data());
    //
    canvFin->Print(psnm2.Data());
}

void PlotDNDEta(int bin, Bool_t forceMCMB)
{

  int bin_ = bin;
  if (forceMCMB) bin_ = 0;
  
    printf("--->>> PlotDNDEta(%d)\n", bin);
    //
    TObjArray *res = &resArr;
    TString prefN = "bin"; prefN += Form("%02d", bin); prefN+="_";
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
    TH1* hZV = (TH1*)res->At(shift + kZvDist);
    if (!res->At(shift + kSigCorr)) return;
    // corrected data
    TString nms =  prefN;
    nms += "DataCorrSignal";
    nms += "_";
    nms += uniqueName;
    TH2* tmpSCorr = (TH2F*)res->At(shift + kSigCorr)->Clone(Form("%stmpSCorrBin%d",prefN.Data(),bin));
    TH1* hsigCorr = 0;
    if (useZbinWAv) {
        CorrectForZV(tmpSCorr,hZV);
        hsigCorr = ProjectWghMean( tmpSCorr, nms.Data());
    }
    else {
      hsigCorr = ProjNorm(tmpSCorr,hZV,nms.Data(),zbinMin,zbinMax); //R+HACK
    }
    //  delete tmpSCorr;
    SetHStyle(hsigCorr,kRed,20,1.0);
    hsigCorr->Scale(1./hsigCorr->GetBinWidth(1));
    hsigCorr->Draw();
    mx = TMath::Max(mx, hsigCorr->GetMaximum());
    mn = TMath::Min(mn, hsigCorr->GetMinimum());
    resDnDeta.AddAtAndExpand( hsigCorr, bin );
    hsigCorr->SetDirectory(0);
    resDnDeta.AddAtAndExpand( res->At(shift + kSigCorr), 100+bin );
    resDnDeta.AddAtAndExpand( tmpSCorr, 200+bin );
    tmpSCorr->SetDirectory(0);
    resDnDeta.AddAtAndExpand(hZV , 250+bin );
    hZV->SetDirectory(0);
    
    //
    // raw data
    TH2* tmpRaw = (TH2F*)res->At(shift+kRawDtCut)->Clone("tmpRaw");
    TH1* hraw = tmpRaw->ProjectionX(prefN+"DataRaw");
    delete tmpRaw;
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
    hraws->Scale(1./hraw->GetBinWidth(1));
    hraws->Draw("same");
    mn = TMath::Min(mn, hraw->GetMinimum());
    mx = TMath::Max(mx, hraw->GetMaximum());
    //
    // bg
    TH2* tmpBg = (TH2F*)res->At(shift+kBgEst)->Clone("tmpBg");
    TH1* hbg = tmpBg->ProjectionX(prefN+"BgEst");
    delete tmpBg;
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
    SetHStyle(hrawsMCLB,kGreen+2,30,1.0);
    hrawsMCLB->Scale(1./hrawsMCLB->GetBinWidth(1));
    hrawsMCLB->Draw("same");
    mn = TMath::Min(mn, hrawsMCLB->GetMinimum());
    mx = TMath::Max(mx, hrawsMCLB->GetMaximum());
    //
    // bg est
    TH2* tmpBgMCEst = (TH2F*)res->At(shift+kBgEst+kMCShift)->Clone("tmpRawBbMCEst");
    TH1* hbgMCEst = tmpBgMCEst->ProjectionX(prefN+"BgEstMC");
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
    SetHStyle(hbgMC,kGreen+2,25,1.0);
    hbgMC->Scale(1./hbgMC->GetBinWidth(1));
    hbgMC->Draw("same");
    mn = TMath::Min(mn, hbgMC->GetMinimum());
    mx = TMath::Max(mx, hbgMC->GetMaximum());
    //
    // mc truth
    printf("=====+++++=====+++++=====+++++=====+++++=====+++++=====+++++\n");
    res->At(shift+kMCPrim+kMCShift)->Print();
    printf("=====+++++=====+++++=====+++++=====+++++=====+++++=====+++++\n");
    TH2* tmpMCTrue = (TH2F*)res->At(shift+kMCPrim+kMCShift)->Clone(Form("%stmpMCTrueBin%d",prefN.Data(),bin));
    TH1* hMCtrue = 0;
    if (useZbinWAv) {
        CorrectForZV(tmpMCTrue,hZVMCNS);
        hMCtrue = ProjectWghMean(tmpMCTrue, prefN+"MCTruth");
    }
    else {
      hMCtrue = ProjNorm(tmpMCTrue,hZVMCNS,prefN+"MCTruth");
    }
    //  delete tmpMCTrue;
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
    TH1* mcdst = (TH1*)res->At(shift+kDataDist+kMCShift);
    TH1* mcdstbg = (TH1*)res->At(shift+kBgDist+kMCShift);
    TH1* mcdstbgLB = (TH1*)res->At(shift+kBgMCDist+kMCShift);
    TH1* dtdst = (TH1*)res->At(shift+kDataDist);
    TH1* dtdstbg = (TH1*)res->At(shift+kBgDist);
    //
    TH1* mcDstN     = (TH1*)FindObject(bin_,HName("Data", useShapeType==kNormShapeDist ? "WDist":"DPhiS"), listMC );
    TH1* mcDstSec   = (TH1*)FindObject(bin_,HName("Sec",  useShapeType==kNormShapeDist ? "WDist":"DPhiS"), listMC );
    TH1* mcDstCombU = (TH1*)FindObject(bin_,HName("CombU",useShapeType==kNormShapeDist ? "WDist":"DPhiS"), listMC );
    TH1* mcDstCombC = (TH1*)FindObject(bin_,HName("Comb", useShapeType==kNormShapeDist ? "WDist":"DPhiS"), listMC );
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
    
    printf("--->>> PlotAplhaBeta(%d)\n", bin);
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
    TPaletteAxis *p = (TPaletteAxis*)dtBet->FindObject("palette");
    if (p) p->SetX1NDC(0.85);
    canvFin->cd(2);
    gPad->SetRightMargin(0.15);
    mcBet->Draw("colz");
    AddLabel("#beta MC (bckg.estimated)",0.2,0.95,kBlack,0.04);
    gPad->Modified();
    mcBet->GetYaxis()->SetTitleOffset(1.4);
    p = (TPaletteAxis*)mcBet->FindObject("palette");
    if (p) p->SetX1NDC(0.85);
    canvFin->cd(3);
    gPad->SetRightMargin(0.15);
    mcBetLB->Draw("colz");
    AddLabel("#beta MC (bckg.from MC labels)",0.2,0.95,kBlack,0.04);
    gPad->Modified();
    mcBetLB->GetYaxis()->SetTitleOffset(1.4);
    p = (TPaletteAxis*)mcBetLB->FindObject("palette");
    if (p) p->SetX1NDC(0.85);
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
    TPaletteAxis *pa = (TPaletteAxis*)dtBet->FindObject("palette");
    if (pa) pa->SetX1NDC(0.85);
    canvFin->cd(5);
    gPad->SetRightMargin(0.15);
    mcAlp->Draw("colz");
    AddLabel("#alpha (bckg.from MC labels)",0.2,0.95,kBlack,0.04);
    gPad->Modified();
    mcAlp->GetYaxis()->SetTitleOffset(1.4);
    pa = (TPaletteAxis*)mcBet->FindObject("palette");
    if (pa) pa->SetX1NDC(0.85);
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
    printf("--->>> PlotSpecies()\n");
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
    
    printf("--->>> CropHisto(\"%p\", %d, %d, %d, %d)\n", histo, bx0, bx1, by0, by1);
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
TObject* FindObject(int bin, const char* nameH, const TList* lst, Bool_t normPerEvent)
{

  //  printf("--->>> FindObject: %d, %s, %p, %d\n", bin, nameH, lst, normPerEvent);
  
    // get histo, optionally normalizing it per processed event
    if (!lst) {printf("FindObject %s: No list provided\n",nameH); return NULL;}
    int nent = lst->GetEntries();
    //
    char buff[200];
    if (bin>=0) sprintf(buff,"b%d_%s",bin,nameH);
    else  sprintf(buff,"%s",nameH);
    TString nm;
    TObject *hst = 0;
    for (int i=nent;i--;) {
        nm = lst->At(i)->GetName();
        if (nm.EndsWith(buff)) {hst = lst->At(i); break;}
    }
    if (!hst) {printf("FindObject: No bin %d %s histo in list %s\n",bin,nameH,lst->GetName()); return NULL;}
    if (!normPerEvent || hst->TestBit(kBitNormPerEvent)) return hst; // already normalized
    TString nameHS = nameH;
    if (nameHS==kHStatName) return hst;                              // never scale stat. histo
    //
    TH1* hstat = (TH1*)FindObject(-1,kHStatName,lst,kFALSE);
    //double nrm = hstat->GetBinContent(kBinEntries + kEvProcInj+bin*kEntriesPerBin);
    double nrm = hstat->GetBinContent(kBinEntries + kEvProcData+bin*kEntriesPerBin); // HACK
    if (nrm<1) {printf("FindObject %s: Anomaluous %d number of events processed in bin %d of list %p\n",
                       buff,int(nrm),bin,lst); return 0;}
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
    if (!fl) {printf("LoadList: No file %s\n",nms.Data()); return NULL;}
    TList* lst = (TList*)fl->Get(nameL);
    if (!lst) {printf("LoadList: No list %s in file %s\n",nameL,nms.Data()); return NULL;}
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
    TH1* hstat = (TH1*)FindObject(-1,kHStatName,lst);
    TString dts = dtType;
    if (dts=="Data") return int( hstat->GetBinContent(kEvProcData) );
    if (dts=="Mix")  return int( hstat->GetBinContent(kEvProcMix) );
    if (dts=="Inj")  return int( hstat->GetBinContent(kEvProcInj) );
    if (dts=="Rot")  return int( hstat->GetBinContent(kEvProcRot) );
    printf("Unknown process %s statistics is checked. Alowed: Data,Mix,Inj,Rot",dtType);
    return 0;
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

//____________________________________________________________________________
void CorrectForZV(TH2* hEtaZ, TH1* hZv)
{
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
        //scle = 0;
        if (scl<1e-9) {printf("Did not find ZV weight for Z=%f, stop\n",z); exit(1);}
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
    // correct for the nonuniformity of Zvertices along Z
    int nbEta = hEtaZ->GetNbinsX();
    int nbZv  = hEtaZ->GetNbinsY();
    if (firstbin<1) firstbin = 1;
    if (lastbin<firstbin) lastbin = nbZv;
    //
    TH1* prj = hEtaZ->ProjectionX(name,firstbin,lastbin,"e");
    prj->Reset();
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
    // for each X bin calculated the weighted average over Y bins
    int nbEta = hEtaZ->GetNbinsX();
    int nbZv  = hEtaZ->GetNbinsY();
    if (firstbin<1) firstbin = 1;
    if (lastbin<firstbin) lastbin = nbZv;  
    //
    printf("HISTO: %s %s %p | bins: %d %d\n",hEtaZ->GetName(),hEtaZ->GetTitle(), hEtaZ, firstbin, lastbin);
    TArrayD val(nbZv),err(nbZv),prc(nbZv);
    TArrayI ind(nbZv);
    //
    TH1* prj = hEtaZ->ProjectionX(name,firstbin,lastbin,"e");
    prj->Reset();
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
