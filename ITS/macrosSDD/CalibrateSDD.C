#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TProfile.h>
#include <TSpline.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TFile.h>
#include <TGrid.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <TLatex.h>
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBStorage.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSCorrMapSDD.h"
#include "AliITSCorrMap1DSDD.h"
#include "AliITSDriftSpeedArraySDD.h"
#include "AliITSDriftSpeedSDD.h"
#include "AliITSresponseSDD.h"
#endif

/*  $Id$    */

/*
  // master equations
  DeltaX_L = -deltaX + deltaV_L*(T_L - T0_true - deltaT0) - V_Ltrue*deltaT0 - corrMapL
  DeltaX_R =  deltaX + deltaV_R*(T_R - T0_true - deltaT0) - V_Rtrue*deltaT0 - corrMapR
  //
  
  where 
  corrMapL(R) : Vdrift vs X nonuniformity correction maps
  V_L(R)true: true drift speed
  deltaV_L(R) : error in drift speed
  V_Ltrue = V_Lassumed - deltaV_L
  V_Rtrue = V_Rassumed - deltaV_R
  //
  DeltaX :  measured residuals (track extrapolation - measurement) for left (L) 
  and right (R) sides (right side multiplied by -1!!!)
  deltaX : misalignment of module (Xtrue = Xassumed - deltaX)
  T_L(R)    : measured drift time (w/o T0 subtraction)
  T0_true   : real T0 
  deltaT0   : error in T0
*/

//------------------- Details of CDB objects to create ---------------
int firstRun = -1;
int lastRun  = -1;
TString cdbComment = "";
//--------------------------------------------------------------------

const Int_t kSDDMin=240,kSDDMax=499,kNSDD=kSDDMax-kSDDMin+1;
enum {kSXvsX,kSXvsZ,kSTDvsX,kSTDvsZ,kSDXvsXclean,kSVDvxX,kSCorrMapX,kSXvsXCorr,kSXvsZCorr,kSDVvsZ,kSDVvsZOrig,kSGrDxx,kSGrTDx,kSGrTXCorr, kNStore};
const Int_t maxVZOrd    = 4;            // max order of polinomial for VD vs Z correction (normal sensors)
const Double_t chi2BadVZ = 2.9;         // treat as bad module if chi2/ndf is worse
//
enum {kDxX,kDxZ,kTDX,kTDZ,kNQATypes};
// Prefixes for residuals histos prepared by PrepSDDResidMap.C from MP2 output
// For QA histograms the name should be changed to "hpSDDResXvsXD","hpSDDResXvsZ","hpSDDDrTimevsXD","hpSDDDrTimevsZ" respectively
const Char_t* kQAHistoNames[] = {    // 
  "Xdrift_Mod_",                     // DX vs Xdrift 
  "ZAnode_Mod_",                     // DX vs local Z
  "TD_vs_Xdrift_Mod_",               // TDrift-T0 vs Xdrift
  "TD_vs_ZAnode_Mod_"                // TDrift-T0 vs local Z
};
const Char_t* kLRNames[] = {"left","right"};  // identifiers for left/right sides in histo names
//
double minBEntryX = 10;        // min entries in the bin to account (vs Xdrift)
double minBEntryZ = 20;        // min entries in the bin to account (vs Z)
double skipDXEdge = 1000.;     // microns to skip from each side of XDrift when smoothing
double edgeSmearMinErr = 60.;        // minimal error for extrapolation
double edgeSmearMaxErr = 600.;       // maximal error for extrapolation
double wDXEdge    = 5000.;     // width in microns for XDrift edges to average
double wVDEdge    = 2500;      // width in microns for VDrift averaging region on the edges
double threshClean= 0.3;       // clean edge bins with occupancy smaller than in the neighbourhood
int    nVDEdge    = 20;        // min width in n bins for VDrift averaging region on the edges
int    smoothLevel= 5;         // smoothing level
int    minBinsEdge = 5;        // min number of X>0 good bins to fit the edge
Bool_t useSplieForR2V = kFALSE;      // if true, extract Vd profile using Derivateve of spline (may be inadequate for sharp turns)
//
Bool_t forceT0CorrTo0   = kTRUE;
Bool_t forceRPhiCorrTo0 = kTRUE;
//
//
AliITSsegmentationSDD sddSeg;
//
//----------------- working variables --------------------------------
Int_t currSDDId=-1,currSDDSide=-1,currMod=-1;
TProfile*  currHDXvsX = 0;                // raw DX vs X histo from QA
TProfile*  currHDXvsZ = 0;                // raw DX vs Z histo from QA
TProfile*  currHTDvsX = 0;                // raw DriftTime-T0 vs X histo from QA
TProfile*  currHTDvsZ = 0;                // raw DriftTime-T0 vs Z histo from QA
TH1*       currHDXvsXclean = 0;           // smoothed DX vs X histo
TH1*       currHVDvxX = 0;                // extracted VDrift vs X profile
TH1*       currHCorrMapX = 0;             // extracted correction map
TH1*       currHDXvsXCorr = 0;            // DX vs X histo after application of the map
TH1*       currHDXvsZCorrMT0 = 0;         // DX vs X histo after application of the map and t0 correction
TH1*       currHDVvsZ = 0;                // correction for the VDrift vs Z
TH1*       currHDVvsZOrig = 0;            // correction for the VDrift vs Z (before error processing)
TGraphErrors* currGrDxx = 0;              // DX final vs XDrift
TGraphErrors* currGrTDx = 0;              // DX final vs TDrift 
TGraphErrors* currGrTXCorr = 0;           // TDrift vs XDritf
TCollection* qaHistos = 0;                // input QA collection
TObjArray procHistos;                     // processed histos buffer
AliITSresponseSDD* sddResp=0;             // SDD response, updated
TObjArray *vdarrayOld = 0;                // olf VDrifts array
TObjArray* vDriftArr=0;                   // VDrifts array, updated
TObjArray* corrMaps=0;                    // holder for correction maps
//
TH1*       resOffsDXraw[2]={0};           // Dx vs X offset at Xdrift=0 raw
TH1*       resOffsDXAMap[2]={0};          // Dx vs X offset (mean) after correction by map
TH1*       resOffsDX[2]={0};              // Dx vs X offset at Xdrift=0 after correction by map (for each side)
TH1*       resVDCorr[2]={0};              // VDrift correction  (for each side)
TH1*       resVDMean[2]={0};              // average VDrift (for each side)
TH1*       resVDCorrZ[2]={0};             // mean VDrift vs Z corrections (for each side)
TH1*       resT0Corr = 0;                 // correction to modules T0
TH1*       resXLocCorr = 0;               // correction to module location
TCanvas*   sddCanv = 0;                   // report canvas
//
Bool_t     sddRespUpdT0OK = kFALSE;       // flag that SDDresponse object T0 was updated
Bool_t     sddRespUpdVDOK = kFALSE;       // flag that SDDresponse object VDrift correction was updated
Bool_t     sddVDriftUpdOK = kFALSE;       // flag that SDD VDrift object was updated
//
TString pathSDDRespOld    = "";
TString pathSDDVDriftOld  = "";
TString pathSDDCorrMapOld = "";
//
//--------------------------------------------------------------------
void      Process(const char* pathSDDResp=0,const char* pathSDDVDrift=0, const char* pathSDDCorrMap=0);
void      Process(TCollection* qa, const char* pathSDDResp=0,const char* pathSDDVDrift=0, const char* pathSDDCorrMap=0);
Bool_t    ProcSDDSideDXvsX();
Bool_t    ProcSDDSideDXvsZ();
TProfile* GetQAHisto(Int_t qaType);
TProfile* H2Profile(TH2* h2);
TH1*      H2ProfileAsTH1(TH2* h2);
TH1*      ProfileAsTH1(TProfile* prof, const char* addName);
TH1*      CleanProfile(TProfile* profR);
TH1*      Vdrift2Resid(TH1* vd);
TH1*      Resid2Vdrift(TH1* res);
void      RedoProfileErrors(TH1* profH1,TProfile* prof);
void      SetModuleID(Int_t mdID,Int_t side=-1);
void      PrepareModuleHistos(Int_t mdID, Int_t side, Bool_t vsZ);
void      CheckResultDX();
void      CalcDXCorrections();
double    GetVOld(double z);
Int_t     GetStoreID(int type, int imd=-1,int side=-1);
void      CleanPrev();
void      StoreCurrent();
double    ftVdZ(double *x, double *par);
void      PlotReport(const char* psname="repSDDQA.ps");
TH1*      GetPadBaseHisto(TPad* pad);
Bool_t    PlotHisto(TH1* h, Option_t* opt="", int mrkStyle=20,int mrkCol=kBlack, double mrkSize=1.);
TLatex*   AddPadLabel(const char*txt,float x=0.1,float y=0.9,int color=kBlack,float size=0.04);
void      GetHistoMinMaxInRange(TH1* h, double &mn,double &mx);
double    edgeLow(double* x, double *par);
double    edgeHigh(double* x, double *par);
void      CureEdges(TH1* prof);
//
Bool_t    LoadSDDVDrift(TString& path, TObjArray *&arr);
Bool_t    LoadSDDResponse(TString& path, AliITSresponseSDD *&resp);
Bool_t    LoadSDDCorrMap(TString& path, TObjArray *&maps);
AliCDBEntry* GetCDBEntry(const char* path);
//
void                UpdateSDDResponse(AliITSresponseSDD *resp, Bool_t t0=kTRUE, Bool_t vdrift=kTRUE);
void                UpdateSDDVDrift(AliITSDriftSpeedArraySDD* vdArr, int imd, int side);
TObjArray*          UpdateSDDVDrift();
TObjArray*          CreateCorrMaps();
AliITSresponseSDD*  UpdateSDDResponse(Bool_t t0=kTRUE, Bool_t vdrift=kTRUE);
AliITSCorrMap1DSDD* CreateCorrMap(TH1* mapHisto, int imd, int side, AliITSCorrMap1DSDD* updateMap=0);
void      PrepCDBObj(TObject *obj,const char* path,int firstrun=0,int lastrun=999999999,const char* comment="");

//-------------------------------------------------------------------

//_________________________________________________________________________
void Process(TCollection* qa, const char* pathSDDResp, const char* pathSDDVDrift, const char* pathSDDCorrMap)
{
  // process all
  qaHistos = qa;
  Process(pathSDDResp,pathSDDVDrift,pathSDDCorrMap);
}

//_________________________________________________________________________
void Process(const char* pathSDDResp, const char* pathSDDVDrift, const char* pathSDDCorrMap)
{
  // process all
  procHistos.Clear();
  pathSDDRespOld     = pathSDDResp;
  pathSDDVDriftOld   = pathSDDVDrift;
  pathSDDCorrMapOld  = pathSDDCorrMap;
  sddRespUpdT0OK = sddRespUpdVDOK = sddVDriftUpdOK = kFALSE;
  //
  if (pathSDDVDriftOld.IsNull()) printf("Attention: Old VDrift is missing!\n");
  //
  for (int imd=kSDDMin;imd<=kSDDMax;imd++) {
    for (int ix=0;ix<2;ix++) {
      CleanPrev();            // clean data from previous module
      printf("Processing %d/%d\n",imd,ix);
      PrepareModuleHistos(imd,ix, kTRUE);
      ProcSDDSideDXvsX();     // process DX vs X histos to get corr.maps, deltaT0, deltaX0, deltaV_mean
      StoreCurrent();
    }
    CalcDXCorrections();      // calculate delta's
    //
    for (int ix=0;ix<2;ix++) {
      CleanPrev();            // clean data from previous module
      PrepareModuleHistos(imd,ix, kFALSE);
      ProcSDDSideDXvsZ();     // process deltaV vs anode profile
      StoreCurrent();         
    }
  }
  //
  corrMaps = CreateCorrMaps();   // create correction maps
  PrepCDBObj(corrMaps,"ITS/Calib/MapsTimeSDD",firstRun,lastRun,cdbComment.Data());
  //
  if (!pathSDDVDriftOld.IsNull()) {
    vDriftArr = UpdateSDDVDrift();
    PrepCDBObj(vDriftArr,"ITS/Calib/DriftSpeedSDD",firstRun,lastRun,cdbComment.Data());
  }
  //
  if (!pathSDDRespOld.IsNull()) {
    sddResp = UpdateSDDResponse();
    PrepCDBObj(sddResp,"ITS/Calib/RespSDD",firstRun,lastRun,cdbComment.Data());
  }
  //
  PlotReport();
}

//_________________________________________________________________________
Bool_t ProcSDDSideDXvsX()
{
  // process one side of the SDD module
  currHDXvsXCorr   = ProfileAsTH1(currHDXvsX,"corrCheck");
  RedoProfileErrors(currHDXvsXCorr,currHDXvsX);
  //
  if ( (currHDXvsXclean=CleanProfile(currHDXvsX)) ) {
    //
    delete currHDXvsXCorr;
    currHDXvsXCorr = (TH1*)currHDXvsXclean->Clone( Form("%s_corrCheck",currHDXvsX->GetName()) );
    //
    currHVDvxX      = Resid2Vdrift(currHDXvsXclean);
    currHCorrMapX   = Vdrift2Resid(currHVDvxX);
    //
    currHDXvsXCorr->Add(currHCorrMapX,-1);
  }
  // check results
  CheckResultDX();
  //
  return kTRUE;
}

//_________________________________________________________________________________
Bool_t ProcSDDSideDXvsZ()
{
  // extract correction for Vdrift vs Z from DXvsZ profile and mean_drift_time vs Z
  //
  static TF1* fitP0 = new TF1("fitP0","pol0",-4000,4000);
  double chi2s[maxVZOrd+1] = {0};
  static TF1* fitVvsZs[maxVZOrd+1] = {0};
  static Bool_t iniDone = kFALSE;
  double zrange = sddSeg.Dz()/2 - 1.;
  if (!iniDone) {
    iniDone = kTRUE;
    for (int iord=0;iord<=maxVZOrd;iord++) {
      fitVvsZs[iord] = new TF1("fitVvsZ",ftVdZ, -zrange, zrange ,iord+2);
      fitVvsZs[iord]->FixParameter(0, iord+0.1);
    }
  }
  //
  int nb = currHDXvsZ->GetNbinsX();  
  if (currHDXvsZ->GetEntries()<minBEntryZ*nb) return kFALSE;
  //
  currHDVvsZ = ProfileAsTH1(currHDXvsZ,"_vzcorr"); 
  currHDVvsZ->Reset();
  //
  currHDXvsZCorrMT0 = ProfileAsTH1(currHDXvsZ,"_zcorrMapT0"); 
  currHDXvsZCorrMT0->Reset();
  //
  int nbUse=0, ib0 = currHDVvsZ->FindBin(-zrange), ib1 = currHDVvsZ->FindBin( zrange);
  double *statZB = new double[nb+1];   // entries per Z bin
  memset(statZB,0,sizeof(double)*(nb+1)); 
  double vmean=0, vmeanE=0, norm=0;
  //
  for (int ib=ib0;ib<=ib1;ib++) {
    double ne = currHTDvsZ->GetBinEntries(ib);
    if (ne<1) continue;
    double dx = currHDXvsZ->GetBinContent(ib);   // raw residual X vs Z
    double dxe= currHDXvsZ->GetBinError(ib);     
    double t  = currHTDvsZ->GetBinContent(ib);   // mean assumed (TDrift - T0)
    double te = currHTDvsZ->GetBinError(ib);
    double vCorr = resVDCorr[currSDDSide]->GetBinContent(currMod+1);
    dx -= vCorr*t;          // subtract mean V correction
    dxe = TMath::Sqrt(dxe*dxe + vCorr*vCorr*te*te);
    currHDXvsZCorrMT0->SetBinContent(ib,dx);
    currHDXvsZCorrMT0->SetBinError(ib,dxe);
  }  
  currHDXvsZCorrMT0->Fit(fitP0,"q0","");
  double pedestal = fitP0->GetParameter(0);
  //
  for (int ib=ib0;ib<=ib1;ib++) {
    double ne = currHTDvsZ->GetBinEntries(ib);
    if (ne<minBEntryZ) continue;
    double dx = currHDXvsZCorrMT0->GetBinContent(ib);   // raw residual X vs Z
    double dxe= currHDXvsZCorrMT0->GetBinError(ib);     
    double t  = currHTDvsZ->GetBinContent(ib);   // mean assumed (TDrift - T0)
    double te = currHTDvsZ->GetBinError(ib);
    if (t<1) continue;
    // corrections
    // account for the corrections leading the mean DX to 0
    dx -= pedestal;
    double v = dx/t;
    double ve = TMath::Sqrt(dxe*dxe + v*v*te*te)/t;
    //
    nbUse++;
    vmeanE += 1./(ve*ve); 
    vmean  += v/(ve*ve); 
    //    printf("%d %f %f | %f %f | -> %f %f\n",ib,dx,dxe,t,te,v,ve);
    int ibDest = currSDDSide==0 ? ib : nb+1-ib;
    statZB[ibDest] = ne;
    currHDVvsZ->SetBinContent(ibDest,v);
    currHDVvsZ->SetBinError(ibDest,ve);
    //
  }
  //
  if (nbUse<maxVZOrd+1) {delete statZB; return kFALSE;}
  if (vmeanE>0) { vmean /= vmeanE; vmeanE = 1./TMath::Sqrt(vmeanE); }
  // fill empty and bad bins by mean value with large error
  vmeanE *= nbUse;
  /*
  for (int ib=ib0;ib<=ib1;ib++) {
    //if (statZB[ib]<minBEntryZ) currHDVvsZ->SetBinContent(ib,vmean);
    if (currHDVvsZ->GetBinError(ib)>vmeanE || statZB[ib]<minBEntryZ) {
      currHDVvsZ->SetBinError(ib,vmeanE*2);
      currHDVvsZ->SetBinContent(ib,vmean);
    }
  }
  */
  // save original
  currHDVvsZOrig = (TH1*) currHDVvsZ->Clone(Form("%s_orig",currHDVvsZ->GetName()));
  // find leftmost good bin
  int bL,bR;
  double errL=0,errR=0, valL=0,valR=0;
  for (bL=ib0;bL<=ib1;bL++) {
    if (/*currHDVvsZ->GetBinError(bL)>vmeanE ||*/ statZB[bL]<minBEntryZ) continue;
    errL = currHDVvsZ->GetBinError(bL);
    valL = currHDVvsZ->GetBinContent(bL);
    break;
  }
  for (bR=ib1;bR>=ib0;bR--) {
    if (/*currHDVvsZ->GetBinError(bR)>vmeanE ||*/ statZB[bR]<minBEntryZ) continue;
    errR = currHDVvsZ->GetBinError(bR);
    valR = currHDVvsZ->GetBinContent(bR);
    break;
  }
  // change bad bins 
  for (int ib=bL-1;ib>=ib0;ib--) {
    double err = errL + (vmeanE-errL)*float(ib-ib0)/(bL-ib0+1);
    double val = vmean + (valL-vmean)*float(ib-ib0)/(bL-ib0+1);
    currHDVvsZ->SetBinError(ib,err);
    currHDVvsZ->SetBinContent(ib,val);//currHDVvsZ->GetBinContent(ib+1));
  }
  for (int ib=bR+1;ib<=ib1;ib++) {
    double err = errR + (vmeanE-errR)*float(ib1-ib)/(ib1-bR+1);
    double val = vmean + (valR-vmean)*float(ib1-ib)/(ib1-bR+1);
    currHDVvsZ->SetBinError(ib,err);
    currHDVvsZ->SetBinContent(ib,val);//currHDVvsZ->GetBinContent(ib-1));
  }
  //
  for (int iord=0;iord<=maxVZOrd;iord++) {
    currHDVvsZ->Fit(fitVvsZs[iord],"q0");
    chi2s[iord] = fitVvsZs[iord]->GetChisquare();
    int ndf = fitVvsZs[iord]->GetNDF();
    if (ndf>0) chi2s[iord] /= ndf;
  }
  //
  // analyse chi2's
  int bestOrd = 0;
  double bestChi2 = 9999;
  for (int iOrd=0;iOrd<=maxVZOrd;iOrd++) {
    if (chi2s[iOrd]<bestChi2-0.5) {bestChi2 = chi2s[iOrd]; bestOrd = iOrd;}
  }
  TF1* fitVvsZ = fitVvsZs[bestOrd];
  currHDVvsZ->Fit(fitVvsZ,"q");
  //  
  // extract mean correction neglecting the constant term
  vmean = 0;
  double freePar = 0;//fitVvsZ->GetParameter(1);
  for (int ib=1;ib<=nb;ib++) {
    if (statZB[ib]<1) continue;
    vmean += statZB[ib]*(fitVvsZ->Eval(currHDVvsZ->GetBinCenter(ib)) - freePar); // account only Z dependent part
    norm  += statZB[ib];
  }
  //
  if (!resVDCorrZ[currSDDSide]) {
    resVDCorrZ[currSDDSide] = new TH1F(Form("VDcorrvsZ%d",currSDDSide),Form("mean VDrift correction vs Z, side %d",currSDDSide),kNSDD,kSDDMin-0.5,kSDDMax+0.5);
    resVDCorrZ[currSDDSide]->SetMarkerColor(2+currSDDSide);
    resVDCorrZ[currSDDSide]->SetMarkerStyle(20+4*currSDDSide);    
  }
  //
  //  if (currSDDSide) vmean = -vmean;
  resVDCorrZ[currSDDSide]->SetBinContent(currMod+1, norm>0 ? vmean/norm : 0);
  delete statZB;
  return kTRUE;
}

//_________________________________________________________________________
void PrepareModuleHistos(Int_t mdID, Int_t side, Bool_t vsX)
{
  // retrieve QA histos, if needed, convert to TH1
  SetModuleID(mdID, side);
  if (vsX) {
    currHDXvsX = GetQAHisto(kDxX);
    currHTDvsX = GetQAHisto(kTDX);
  }
  else {
    currHTDvsZ = GetQAHisto(kTDZ);
    currHDXvsZ = GetQAHisto(kDxZ);
  }
}

//_________________________________________________________________________
void SetModuleID(Int_t mdID,Int_t side)
{
  // set current module ID,side
  currSDDId   = mdID;
  currSDDSide = side;
  currMod     = currSDDId - kSDDMin;
}

//_________________________________________________________________________
TProfile* GetQAHisto(Int_t qaType)
{
  // retrieve needed QA histo
  if (!qaHistos) {printf("QA Histos collection is not set\n"); exit(1);}
  //
  if (currSDDId<kSDDMin || currSDDId>kSDDMax || currSDDSide<0 || currSDDSide>2) 
    {printf("Illegal SDD module ID/Side: %d/%d\n",currSDDId,currSDDSide); exit(1);}
  //
  if (qaType<0 || qaType>=kNQATypes) 
    {printf("Illegal QA Histo Type %d\n",qaType); exit(1);}
  //
  const char* hname = Form("%s%d_%s",kQAHistoNames[qaType],currSDDId,kLRNames[currSDDSide]);
  TH1* h = (TH1*) qaHistos->FindObject(hname);
  if (!h) {printf("Did not find histo %s in QA output collection\n",hname); exit(1);}
  //
  if (h->InheritsFrom(TH2::Class())) h = H2Profile((TH2*)h); // make sure it is not TH2 histo
  //
  return (TProfile*)h;
}

//_________________________________________________________________________
TH1* H2ProfileAsTH1(TH2* h2)
{
  // extract profile as TH1 histo
  TProfile* prf = h2->ProfileX();
  TH1* prof = ProfileAsTH1(prf, "_profH1");
  delete prf;
  return prof;
}

//_________________________________________________________________________
TH1* ProfileAsTH1(TProfile* prf, const char* addName)
{
  // convert profile to TH1 histo
  TString nm = prf->GetName(); nm += addName;
  TAxis* xax = prf->GetXaxis();
  TH1* prof = new TH1F(nm.Data(),nm.Data(),prf->GetNbinsX(),xax->GetXmin(),xax->GetXmax());
  for (int i=1;i<=prof->GetNbinsX();i++) {
    prof->SetBinContent(i, prf->GetBinContent(i));
    prof->SetBinError(i, prf->GetBinError(i));
  }
  return prof;
}

//_________________________________________________________________________
TProfile* H2Profile(TH2* h2)
{
  // extract profile 
  TString nm = h2->GetName(); nm += "_prof";
  TProfile* prf = h2->ProfileX(nm.Data());
  return prf;
}

//_________________________________________________________________________________
TH1* CleanProfile(TProfile* profR)
{
  // clean profile copy
  int ib0 = profR->FindBin(1), ib1 = profR->FindBin(sddSeg.Dx()-1);
  int nbn = profR->GetNbinsX();
  int nSkip = skipDXEdge/profR->GetBinWidth(1);
  int nMean = wDXEdge/profR->GetBinWidth(1);
  //
  TH1* prof = ProfileAsTH1(profR,"_Clean");
  RedoProfileErrors(prof,profR);
  //
  CureEdges(prof);
  //
  //  for (int i=1;i<ib0;i++) {prof->SetBinContent(i,0); prof->SetBinError(i,0);}
  //  for (int i=ib1+1;i<=nbn;i++) {prof->SetBinContent(i,0); prof->SetBinError(i,0);}
  // get mean occupancy skipping first and last nSkip occupied bins
  //
  int nbCntL=0, nbCntR=0,nskipL=0,nskipR=0;
  double smL=0,smR=0;
  for (int ib=ib0;ib<=ib1 && nbCntL<nMean;ib++) { // skipping left, get average of good bins
    double vl = profR->GetBinEntries(ib);
    if (vl<minBEntryX) continue;
    if (nskipL<nSkip) {nskipL++; continue;}
    smL += vl;
    nbCntL++;
  }
  if (nbCntL<1) return 0;
  smL /= nbCntL;
  //
  for (int ib=ib1+1;ib>ib0 && nbCntR<nMean;ib--) { // skipping right
    double vl = profR->GetBinEntries(ib);
    if (vl<minBEntryX) continue;
    if (nskipR<nSkip) {nskipR++; continue;}
    smR += vl;
    nbCntR++;
  }
  //
  if (nbCntR<1) return 0;
  smR /= nbCntR;
  //
  prof->GetXaxis()->SetRange(ib0,ib1);
  prof->Smooth(smoothLevel,"r");
  prof->GetXaxis()->SetRange(1,nbn);
  // kill bins with small statistics, from left and right
  for (int ib=1;ib<ib1;ib++)   if (profR->GetBinEntries(ib)<threshClean*smL) {prof->SetBinContent(ib,0);prof->SetBinError(ib,0);} else break;
  for (int ib=ib1;ib>ib0;ib--) if (profR->GetBinEntries(ib)<threshClean*smR) {prof->SetBinContent(ib,0);prof->SetBinError(ib,0);} else break;
  //
  return prof;
}


//_________________________________________________________________________________
TH1* Resid2Vdrift(TH1* res)
{
  // extract Vdrift profile from residuals profile
  //
  TSpline3 spl(res);
  TString nm = res->GetName(); nm += "_vd";
  TH1* vd = (TH1*) res->Clone(nm.Data());
  vd->Reset();
  int nb = vd->GetNbinsX();
  //
  int nbcnt=0, nbfl = wVDEdge/vd->GetBinWidth(1);
  if (nbfl<nVDEdge) nbfl=nVDEdge;
  //
  double vAv=0,eAv=0;                           // average vdrift / v0
  for (int i=1;i<=nb;i++) {
    double deriv = useSplieForR2V ? spl.Derivative(vd->GetBinCenter(i)) : (res->GetBinContent(i)-res->GetBinContent(i-1))/res->GetBinWidth(i);
    double v = TMath::Abs(deriv-1)>1e-6 ? 1./(1-deriv) : 1.;
    vd->SetBinContent(i, v);
    double err2 = 0; // set relative error
    int nbd = 1;
    for (int i1=i-1;i1<=i+1;i1++) {
      if (i1<1 || i1>nb) continue;
      double err = res->GetBinError(i);
      if (err<1e-9) err = 1e4;
      err2 += err*err;
      nbd++;
    }
    if (nbd) err2/=nbd;
    if (err2>0 && v>0 && v<3) {vAv += v/err2; eAv += 1./err2;}
    err2 = TMath::Sqrt(err2)/res->GetBinWidth(1);
    vd->SetBinError(i, err2 );
    nbcnt++;
  }
  vAv = eAv>0 ? vAv/eAv : 1.0;
  //  printf("mean V = %f in %d bins\n",vAv,nbcnt);
  // cure anomalous bins
  //  for (int i=1;i<=nb;i++) if (vd->GetBinContent(i)<0.1) vd->SetBinContent(i,vAv);
  //
  int ib = 0;
  for (ib=1;ib<nb;ib++) { // detect up to which bin the left tail is unstable
    double x = vd->GetBinCenter(ib);
    if (x<0 || res->GetBinError(ib)<1e-9) {vd->SetBinContent(ib,1); vd->SetBinError(ib,0); continue;}
    double vl = vd->GetBinContent(ib);
    if (TMath::Abs(vl-vAv)<0.2) break;
    vd->SetBinContent(ib,1);
    vd->SetBinError(ib,0);    
  }
  //
  double vL=0,eL=0,vR=0,eR=0;
  // get the mean of leftmost stable vdrift
  int lastCounted = 0;
  nbcnt = 0;
  for (int i=ib;i<=nb && nbcnt<nbfl;i++) {
    double berr = vd->GetBinError(i);
    if (berr<1e-9 || berr>0.5 || TMath::Abs(vd->GetBinContent(i)-vAv)>0.2) continue;
    vL += vd->GetBinContent(i)/berr/berr;
    eL += 1./berr/berr;
    nbcnt++;
    lastCounted = i;
  }
  vL = eL>0 ? vL/eL : vAv;
  //printf("VLeft: %f, in %d bins (%d %d)\n",vL,nbcnt,ib,lastCounted);
  for (int i=1;i<=ib;i++) vd->SetBinContent(i, vL);
  // for safety check if there are no outliers in first 3 "stable" bins
  for (int i=ib+1;i<ib+3;i++) if (  vd->GetBinError(i)<1e-9 || TMath::Abs(vd->GetBinContent(i)-vL)>0.2) vd->SetBinContent(i, vL);
  vd->SetBinContent(vd->FindBin(1),1.); // no correction at t=0 !!
  //
  double lmax = sddSeg.Dx()-1;
  for (ib=nb+1;ib--;) { // detect up to which bin the right tail is unstable
    double x = vd->GetBinCenter(ib);
    if (x>=lmax || res->GetBinError(ib)<1e-9) {vd->SetBinContent(ib,1);vd->SetBinError(ib,0);  continue;}
    if (TMath::Abs(vd->GetBinContent(ib)-vAv)<0.4) break;
    vd->SetBinContent(ib,vAv);
    vd->SetBinError(ib,0);    
  }
  nbcnt= 0;
  lastCounted = 0;
  for (int i=ib;i>=1 && nbcnt<nbfl;i--) {
    double berr = vd->GetBinError(i);
    if (berr<1e-9 || berr>0.5 || TMath::Abs(vd->GetBinContent(i)-vAv)>0.4) continue;
    vR += vd->GetBinContent(i)/berr/berr;
    eR += 1./berr/berr;
    nbcnt++;
    lastCounted = i;
  }
  vR = eR>0 ? vR/eR : vAv;
  // printf("VRight: %f, in %d bins (%d %d)\n",vR,nbcnt,lastCounted,ib);
  for (int i=ib;i<=nb;i++) vd->SetBinContent(i, vR);
  // for safety check if there are no outliers in first 3 "stable bins
  for (int i=(lastCounted+ib)/2;i<ib;i++) 
    if (  vd->GetBinError(i)<1e-9 || (vd->GetBinError(i)>0.3 && TMath::Abs(vd->GetBinContent(i)-vR)>0.2)) vd->SetBinContent(i, vR);
  // fit the empty bins on the right
  //  vd->Fit(fitP1,"+","",vd->GetBinCenter(ib-nbfl),vd->GetBinCenter(nb));
  //for (int i=ib;i<=nb;i++) vd->SetBinContent(i, fitP1->Eval(vd->GetBinCenter(i)));
  //
  return vd;
}


//_________________________________________________________________________________
TH1* Vdrift2Resid(TH1* vd)
{
  // convert Vdrift profile to residuals profile (final correction map)
  //
  TString nm = vd->GetName(); nm += "_vd2res";
  TH1* res = (TH1*) vd->Clone(nm.Data());
  res->Reset();
  int ib0 = res->FindBin(1);
  int ib1 = res->FindBin(sddSeg.Dx()-1);
  double resv=0;
  int lastB=0;
  for (lastB=ib1+1;lastB--;) if (vd->GetBinError(lastB)>1e-9) break;   // find last good bin
  double lastX = vd->GetBinCenter(lastB);
  // extend by 1mm
  lastX += 1000;
  if (lastX>sddSeg.Dx()) lastX = sddSeg.Dx();
  lastB = vd->FindBin(lastX);
  // 
  // 1st iteration : estimate correction at max Xdrift
  for (int i=ib0;i<=lastB;i++) {
    double dx = res->GetBinWidth(i);
    double v = vd->GetBinContent(i);
    resv += dx*(1.-1./v);    
  }
  double vcorr = (resv)/lastX;
  //
  // 2nd iteration : create new residuals forcing them to be 0 at Xdrift=0 and maxXDrift
  resv = res->GetBinWidth(ib0)*vcorr;
  for (int i=ib0;i<=lastB;i++) {
    double dx = res->GetBinWidth(i);
    double v = vd->GetBinContent(i);
    resv += dx*(1.-1./v - vcorr);
    res->SetBinContent(i, resv);
  }
  //
  return res;
}

//__________________________________________________________________________________
void CheckResultDX()
{
  // check mean residuals before and after correction
  const double kFOffs = 0.05; // skip edges
  static TF1* fitP1 = new TF1("fitP1","pol1",-5000,40000);
  static TF1* fitP0 = new TF1("fitP0","pol0",-5000,40000);
  //
  if (currHDXvsX->GetEntries()<minBEntryX*currHDXvsX->GetNbinsX()) {
    if (currHCorrMapX) currHCorrMapX->Reset();
    return;
  }
  //
  // vdrift correction 
  int b0 = currHDXvsXCorr->FindBin(1),b1 = currHDXvsXCorr->FindBin(sddSeg.Dx()-1);
  int nb = 0;
  for (int ib=b0;ib<b1;ib++) if (currHTDvsX->GetBinEntries(ib)>=minBEntryX) nb++;
  currHDXvsX->Fit(fitP0,"q","");
  double offsRaw = fitP0->GetParameter(0);
  currHDXvsXCorr->Fit(fitP0,"q0","");
  double offsAMap = fitP0->GetParameter(0);

  //
  currGrDxx    = new TGraphErrors(nb);               // residual vs Xdrift
  currGrTDx    = new TGraphErrors(nb);               // residual vs tdrift
  currGrTXCorr = new TGraphErrors(nb);               // Xdrift vs tdrift
  double tmin = 1e6, tmax = -1e6, xmin = 1e6, xmax = -1e6;
  int ip = 0;
  for (int i=b0;i<b1;i++) {
    if (currHTDvsX->GetBinEntries(i)<minBEntryX) continue;
    double t = currHTDvsX->GetBinContent(i);
    double x = currHDXvsXCorr->GetBinCenter(i);
    if (tmin>t) tmin = t;
    if (tmax<t) tmax = t;
    if (xmin>x) xmin = x;
    if (xmax<x) xmax = x;    
    currGrDxx->SetPoint(ip,x,currHDXvsXCorr->GetBinContent(i));
    currGrDxx->SetPointError(ip,currHDXvsXCorr->GetBinWidth(i),currHDXvsXCorr->GetBinError(i));
    //
    currGrTDx->SetPoint(ip,t, currHDXvsXCorr->GetBinContent(i));
    currGrTDx->SetPointError(ip, currHTDvsX->GetBinError(i), currHDXvsXCorr->GetBinError(i));    
    //
    currGrTXCorr->SetPoint(ip,t, currHTDvsX->GetBinCenter(i));
    currGrTXCorr->SetPointError(ip,currHTDvsX->GetBinError(i), currHTDvsX->GetBinWidth(i));    
    //
    ip++;
  }
  double del = tmax-tmin;
  tmin += kFOffs*del;
  tmax -= kFOffs*del;  
  del = xmax - xmin;
  xmin += kFOffs*del;
  xmax -= kFOffs*del;
  //
  currGrDxx->Fit(fitP1,"q","",xmin, xmax);
  double offs = fitP1->GetParameter(0);           // offset of correction line at Xdrift=0
  //
  //  printf("Fitting VD correction in the range %.1f:%.1f\n",tmin,tmax);
  currGrTDx->Fit(fitP1,"q","",tmin,tmax);
  double vcor = fitP1->GetParameter(1);
  //
  currGrTXCorr->Fit(fitP1,"q","",tmin,tmax);
  double vav = fitP1->GetParameter(1);
  //
  // store results
  if (!resOffsDXraw[currSDDSide]) {
    resOffsDXraw[currSDDSide] = new TH1F(Form("OffsRaw%d",currSDDSide),Form("DX Raw Offset, side %d",currSDDSide),kNSDD,kSDDMin-0.5,kSDDMax+0.5);
    resOffsDXraw[currSDDSide]->SetMarkerColor(2+currSDDSide);
    resOffsDXraw[currSDDSide]->SetMarkerStyle(20+4*currSDDSide);
  }
  //
  if (!resOffsDXAMap[currSDDSide]) {
    resOffsDXAMap[currSDDSide] = new TH1F(Form("OffsAMap%d",currSDDSide),Form("DX Offset after Map corr. Mean, side %d",currSDDSide),kNSDD,kSDDMin-0.5,kSDDMax+0.5);
    resOffsDXAMap[currSDDSide]->SetMarkerColor(3+currSDDSide);
    resOffsDXAMap[currSDDSide]->SetMarkerStyle(20+4*currSDDSide);
  }
  //
  if (!resOffsDX[currSDDSide]) {
    resOffsDX[currSDDSide] = new TH1F(Form("Offs%d",currSDDSide),Form("DX Offset at TD=0 after Map corr., side %d",currSDDSide),kNSDD,kSDDMin-0.5,kSDDMax+0.5);
    resOffsDX[currSDDSide]->SetMarkerColor(2+currSDDSide);
    resOffsDX[currSDDSide]->SetMarkerStyle(20+4*currSDDSide);    
  }
  //
  if (!resVDCorr[currSDDSide]) {
    resVDCorr[currSDDSide] = new TH1F(Form("VDcorr%d",currSDDSide),Form("VDrift correction, side %d",currSDDSide),kNSDD,kSDDMin-0.5,kSDDMax+0.5);
    resVDCorr[currSDDSide]->SetMarkerColor(2+currSDDSide);
    resVDCorr[currSDDSide]->SetMarkerStyle(20+4*currSDDSide);    
  }
  //
  if (!resVDMean[currSDDSide]) {
    resVDMean[currSDDSide] = new TH1F(Form("VDmean%d",currSDDSide),Form("VDrift mean, side %d",currSDDSide),kNSDD,kSDDMin-0.5,kSDDMax+0.5);
    resVDMean[currSDDSide]->SetMarkerColor(2+currSDDSide);
    resVDMean[currSDDSide]->SetMarkerStyle(20+4*currSDDSide);    
  }
  //
  resOffsDXraw[currSDDSide]->SetBinContent(currMod+1, offsRaw);
  resOffsDXAMap[currSDDSide]->SetBinContent(currMod+1, offsAMap);
  resOffsDX[currSDDSide]->SetBinContent(currMod+1, offs);
  resVDCorr[currSDDSide]->SetBinContent(currMod+1, vcor);
  resVDMean[currSDDSide]->SetBinContent(currMod+1, vav);
  //
}

//__________________________________________________________________________
void CalcDXCorrections()
{
  // estimate time0 and alignment correction for the whole module
  if (!resT0Corr) {
    resT0Corr = new TH1F("T0Corr","T0 Correction",kNSDD,kSDDMin-0.5,kSDDMax+0.5);
    resT0Corr->SetMarkerColor(2);
    resT0Corr->SetMarkerStyle(20);    
  }
  //
  if (!resXLocCorr) {
    resXLocCorr = new TH1F("XLocCorr","XLoc Correction",kNSDD,kSDDMin-0.5,kSDDMax+0.5);
    resXLocCorr->SetMarkerColor(2);
    resXLocCorr->SetMarkerStyle(20);    
  }
  //
  if (!resVDMean[0] || !resVDMean[1]) return;
  if (!resOffsDX[0] || !resOffsDX[1]) return;
  double vL = resVDMean[0]->GetBinContent(currMod+1);   // true mean VL
  double vR = resVDMean[1]->GetBinContent(currMod+1);   // true mean VR
  double offsL = resOffsDX[0]->GetBinContent(currMod+1);
  double offsR = resOffsDX[1]->GetBinContent(currMod+1);
  //
  double vsum=0,t0Corr=0,xlCorr=0;
  if (vL>1 && vR>1) { // both sides available
    vsum = vL + vR;
    t0Corr = -(offsL+offsR)/vsum;
    xlCorr = -(offsL*vR - offsR*vL)/vsum;
  }
  /*
  else if (vL>1) t0Corr = -offsL/vL; // only one side is available
  else if (vR>1) t0Corr = -offsR/vR;
  */
  else if (vL>1) xlCorr = -offsL; // only one side is available
  else if (vR>1) xlCorr =  offsR;
  //
  //  printf("SDD%d VL:%f VR:%f offsL:%+f offsR:%+f  dT:%+f dX:%+f\n",currSDDId, vL,vR, offsL,offsR, t0Corr,xlCorr);
  resT0Corr->SetBinContent(currMod+1, t0Corr);        // T0 correction
  resXLocCorr->SetBinContent(currMod+1, xlCorr);      // X alignment correction
  //
  double addMap[2]={0,0};
  Bool_t redoMaps = kFALSE;
  //
  if (forceT0CorrTo0) { // T0 correction was forced to be 0, attribute this to map
    addMap[0] -= vL*t0Corr;
    addMap[1] -= vR*t0Corr;
    redoMaps = kTRUE;
  }
  if (forceRPhiCorrTo0) { // alignment correction was forced to be 0, attribute this to map
    addMap[0] -=  xlCorr;
    addMap[1] -= -xlCorr;
    redoMaps = kTRUE;
  }
  //
  if (redoMaps) {
    for (int ix=0;ix<2;ix++) {
      TH1* map  = (TH1*)procHistos.At( GetStoreID(kSCorrMapX, currSDDId, ix) );
      TH1* mapc = (TH1*)procHistos.At( GetStoreID(kSXvsXCorr, currSDDId, ix) );
      if (!map || !mapc) continue;
      int ib0 = map->FindBin(1);
      int ib1 = map->FindBin(sddSeg.Dx()-1);
      for (int ib=ib0+1;ib<ib1;ib++) {
	map->AddBinContent(ib,  addMap[ix]);
	mapc->AddBinContent(ib, -addMap[ix]);
      }
    }
  }

  //
}

//______________________________________________________________
Int_t GetStoreID(int type, int imd,int side)
{
  // entry of the histo/graph of type in the procHistos array
  //
  if (imd<0)  imd  = currSDDId;
  if (side<0) side = currSDDSide;
  if (type<0||type>=kNStore || imd<kSDDMin || imd>kSDDMax || side<0 || side>1) {
    printf("Wrong object requested: type: %d, Mod:%d/%d\n",type,imd,side);
    exit(1);
  }
  return (2*(imd-kSDDMin)+side)*kNStore + type;
}

//______________________________________________________________
void CleanPrev()
{
  // clean "current" objects from last event
  currHDXvsX = 0;
  currHDXvsZ = 0;
  currHTDvsX = 0;
  currHTDvsZ = 0;
  currHDXvsXclean = 0;
  currHVDvxX = 0;
  currHCorrMapX = 0;
  currHDXvsXCorr = 0;
  currHDVvsZ = 0;
  currGrDxx = 0;
  currGrTDx = 0;
  currGrTXCorr = 0;
  //
}

//______________________________________________________________
void StoreCurrent()
{
  // store "current" objects in procHistos
  if (currHDXvsX)       procHistos.AddAtAndExpand(currHDXvsX,      GetStoreID(kSXvsX));
  if (currHDXvsZ)       procHistos.AddAtAndExpand(currHDXvsZ,      GetStoreID(kSXvsZ));
  if (currHTDvsX)       procHistos.AddAtAndExpand(currHTDvsX,      GetStoreID(kSTDvsX));
  if (currHTDvsZ)       procHistos.AddAtAndExpand(currHTDvsZ,      GetStoreID(kSTDvsZ));
  if (currHDXvsXclean)  procHistos.AddAtAndExpand(currHDXvsXclean, GetStoreID(kSDXvsXclean));
  if (currHVDvxX)       procHistos.AddAtAndExpand(currHVDvxX,      GetStoreID(kSVDvxX));
  if (currHCorrMapX)    procHistos.AddAtAndExpand(currHCorrMapX,   GetStoreID(kSCorrMapX));
  if (currHDXvsXCorr)   procHistos.AddAtAndExpand(currHDXvsXCorr,  GetStoreID(kSXvsXCorr));
  if (currHDXvsZCorrMT0) procHistos.AddAtAndExpand(currHDXvsZCorrMT0,  GetStoreID(kSXvsZCorr));
  if (currHDVvsZ)       procHistos.AddAtAndExpand(currHDVvsZ,      GetStoreID(kSDVvsZ));
  if (currHDVvsZOrig)   procHistos.AddAtAndExpand(currHDVvsZOrig,  GetStoreID(kSDVvsZOrig));
  if (currGrDxx)        procHistos.AddAtAndExpand(currGrDxx,       GetStoreID(kSGrDxx));
  if (currGrTDx)        procHistos.AddAtAndExpand(currGrTDx,       GetStoreID(kSGrTDx));
  if (currGrTXCorr)     procHistos.AddAtAndExpand(currGrTXCorr,    GetStoreID(kSGrTXCorr));
  //
}

//_________________________________________________________________________________
TObjArray* CreateCorrMaps()
{
  // create correction maps for all modules
  printf("Creating correction maps (update %s)\n",pathSDDCorrMapOld.Data());
  TObjArray *dest = new TObjArray(2*kNSDD);
  TObjArray* update = 0;
  if (!pathSDDCorrMapOld.IsNull() && !LoadSDDCorrMap(pathSDDCorrMapOld,update)) {
    printf("The update of correction map was requested but the source %s is not found\n",pathSDDCorrMapOld.Data());
    exit(1);
  }
  //
  dest->Clear();
  AliITSCorrMap1DSDD *updMap = 0;
  for (int imd=kSDDMin;imd<=kSDDMax;imd++) {
    for (int side=0;side<2;side++) {
      TH1* mph = (TH1*)procHistos.At( GetStoreID(kSCorrMapX,imd,side) );
      //if (!mph) printf("Correction map is missing for module %d/%d\n",imd,side);
      if (update) updMap = (AliITSCorrMap1DSDD*)update->At(2*(imd-kSDDMin) + side);
      AliITSCorrMap1DSDD* mp = CreateCorrMap(mph,imd,side, updMap);
      dest->AddAtAndExpand(mp, 2*(imd-kSDDMin) + side);
    }
  }
  //
  return dest;
}

//_________________________________________________________________________________
AliITSCorrMap1DSDD* CreateCorrMap(TH1* mapHisto, int imd, int side, AliITSCorrMap1DSDD* updateMap)
{
  // create or update correction map from histo
  int nbCorr = 1, nbOld = 0;
  int b0=0,b1=0;
  if (mapHisto) {
    b0 = mapHisto->FindBin(1);
    b1 = mapHisto->FindBin(sddSeg.Dx()-1);
  }
  nbCorr = b1-b0+1;
  AliITSCorrMap1DSDD* mpCorr = 0;
  //
  if (updateMap && nbCorr>1) {
    if (mapHisto) {
      TSpline3 spl(mapHisto);
      nbOld = updateMap->GetNBinsDrift();
      double dx = sddSeg.Dx()/nbOld;
      for (int ip=0;ip<nbOld;ip++) {
	double x = dx*(0.5+ip);
	updateMap->SetCellContent(0,ip,updateMap->GetCellContent(0,ip)-spl.Eval(x));
      }
    }
    mpCorr = updateMap;
  }
  else {
    mpCorr = new AliITSCorrMap1DSDD(Form("DriftTimeMap_%d_%d",imd,side),nbCorr);
    if (side==0) mpCorr->SetInversionBit(); // !!! left side should return correction*-1
    if (mapHisto) for (int ib=b0;ib<=b1;ib++) mpCorr->SetCellContent(0,ib-b0,-mapHisto->GetBinContent(ib));
  }
  //
  return mpCorr;
}

//_________________________________________________________________________________
TObjArray* UpdateSDDVDrift()
{
  // retrieve SDD VDrift object and update it
  if (!vdarrayOld && !LoadSDDVDrift(pathSDDVDriftOld,vdarrayOld)) return 0;
  TObjArray *vdarrayUp = new TObjArray(2*kNSDD);
  //
  for (int imd=kSDDMin;imd<=kSDDMax;imd++) {
    for (int side=0;side<2;side++) {
      int iad = 2*(imd-kSDDMin)+side;
      AliITSDriftSpeedArraySDD* drarr = (AliITSDriftSpeedArraySDD*) vdarrayOld->At( iad );
      AliITSDriftSpeedArraySDD* drarrUp = new AliITSDriftSpeedArraySDD();
      AliITSDriftSpeedSDD* vOr = drarr->GetDriftSpeedObject(0);
      AliITSDriftSpeedSDD* vUp = new AliITSDriftSpeedSDD(*vOr); 
      drarrUp->AddDriftSpeed(vUp);
      vdarrayUp->AddAt(drarrUp, iad);
      UpdateSDDVDrift(drarrUp, imd, side);
    }
  }
  //
  sddVDriftUpdOK = kTRUE;
  return vdarrayUp;
}


//_________________________________________________________________________________
void UpdateSDDVDrift(AliITSDriftSpeedArraySDD* vdArr, int imd, int side)
{
  // update vdrift vs anode in the object
  AliITSDriftSpeedSDD* ds;
  if (!vdArr || !(ds=vdArr->GetDriftSpeedObject(0))) {printf("No VDrift object for  module %d/%d\n",imd,side); exit(1);}
  TH1* vdh = (TH1*)procHistos.At( GetStoreID(kSDVvsZ,imd,side) );
  if (!vdh) {
    //printf("VDrift vs Z correction is not processed for module %d/%d\n",imd,side); 
    return;
  }
  TF1* fp = vdh->GetFunction("fitVvsZ");
  if (!fp)  {printf("VDrift vs Z correction fit is missing SDD%d/%d\n",imd,side); return;}
  //
  int ord = (int)fp->GetParameter(0); // 1st param is the order of poly
  int ordOld = ds->GetDegreeofPoly();
  if (ord>ordOld) ds->SetDegreeofPoly(ord);
  for (int ip=0;ip<ord+1;ip++) {       // don't store offset (par[1])
    double par = ds->GetDriftSpeedParameter(ip) - fp->GetParameter(ip+1);
    ds->SetDriftSpeedParameter(ip, par);
  }
  //
}

//_________________________________________________________________________________
AliITSresponseSDD* UpdateSDDResponse(Bool_t t0, Bool_t vdrift)
{
  // retrieve RespSDD object and update it
  AliITSresponseSDD* resp = 0;
  if (!LoadSDDResponse(pathSDDRespOld, resp)) return 0;
  UpdateSDDResponse(resp, t0, vdrift);
  sddRespUpdT0OK = t0;
  sddRespUpdVDOK = vdrift;
  //
  return resp;
}

//_________________________________________________________________________________
void UpdateSDDResponse(AliITSresponseSDD *resp, Bool_t t0, Bool_t vdrift)
{
  // update the map with extracted values
  printf("Updating RespSDD object: T0:%s VDrift:%s\n",t0?"ON":"OFF",vdrift?"ON":"OFF");
  //
  if (t0 && !resT0Corr) 
    {printf("T0 update is requested but corrections were not processed"); exit(1);}
  if (vdrift && !(resVDCorr[0] && resVDCorr[1]))
    {printf("VDrift update is requested but corrections were not processed"); exit(1);}
  //
  for (int imd=kSDDMin;imd<=kSDDMax;imd++) {
    if (t0 && !forceT0CorrTo0) resp->SetModuleTimeZero(imd, resp->GetTimeZero(imd) - resT0Corr->GetBinContent(imd-kSDDMin+1));
    if (vdrift) {
      for (int ix=0;ix<2;ix++) {
	double vdZ = sddVDriftUpdOK&&resVDCorrZ[ix] ? resVDCorrZ[ix]->GetBinContent(imd-kSDDMin+1) : 0; // contribution from DXvsZ correction
	double vdX = resVDCorr[ix]->GetBinContent(imd-kSDDMin+1); // contribution from DXvsX correction
	resp->SetDeltaVDrift(imd, resp->GetDeltaVDrift(imd,ix) - (vdX-vdZ), ix);
      }
    }
  }
  //
}

//___________________________________________________________________
double GetVOld(double z)
{
  // return VDrift assumed in reconstruction
  if (!vdarrayOld && !LoadSDDVDrift(pathSDDVDriftOld,vdarrayOld)) return 0;
  AliITSDriftSpeedArraySDD* drarr = (AliITSDriftSpeedArraySDD*) vdarrayOld->At( 2*currMod + currSDDSide);
  float anode = sddSeg.GetAnodeFromLocal( currSDDSide==0 ? 1.:-1. ,z*1e-4);
  double v = drarr->GetDriftSpeed(0, anode);
  return v;
}

//___________________________________________________________________
double ftVdZ(double *x, double *par)
{
  // function to fit the vdrift dependence on Z
  //
  // convert Z to anode
  double z = x[0];
  double ian = (z/sddSeg.Dz() + 0.5);
  if (ian<0) ian = 0.;
  else if (ian>1) ian = 1.;
  ian *= sddSeg.NpzHalf();
  //
  int ord = int(par[0]);
  double v = par[ord+1];
  for (int i=ord;i--;) v = par[i+1]+ian*v;
  return v;
}

//________________________________________________________________________________________________________
Bool_t LoadSDDVDrift(TString& path, TObjArray *&arr)
{
  // load VDrift object
  if (path.IsNull()) return kFALSE;
  printf("Loading SDD VDrift from %s\n",path.Data());
  //
  AliCDBEntry *entry = 0;
  delete arr;
  arr = 0;
  while(1) {
    if (path.BeginsWith("path: ")) { // must load from OCDB
      entry = GetCDBEntry(path.Data());
      if (!entry) break;
      arr = (TObjArray*) entry->GetObject();
      entry->SetObject(NULL);
      entry->SetOwner(kTRUE);
      break;
    }
    //
    if (gSystem->AccessPathName(path.Data())) break;
    TFile* precf = TFile::Open(path.Data());
    if (precf->FindKey("TObjArray")) arr = (TObjArray*)precf->Get("TObjArray");
    else if (precf->FindKey("AliCDBEntry") && (entry=(AliCDBEntry*)precf->Get("AliCDBEntry"))) {
      arr = (TObjArray*) entry->GetObject();
      if (arr && arr->InheritsFrom(TObjArray::Class())) entry->SetObject(NULL);
      else arr = 0;
      entry->SetObject(NULL);
      entry->SetOwner(kTRUE);
      delete entry;
    }
    //
    precf->Close();
    delete precf;
    break;
  } 
  //
  if (!arr) {printf("Failed to load SDD vdrift from %s\n",path.Data()); return kFALSE;}
  arr->SetOwner(kTRUE);
  return kTRUE;
}

//________________________________________________________________________________________________________
Bool_t LoadSDDResponse(TString& path, AliITSresponseSDD *&resp)
{
  // load SDD response
  if (path.IsNull()) return kFALSE;
  printf("Loading SDD response from %s\n",path.Data());
  //
  AliCDBEntry *entry = 0;
  delete resp;
  resp = 0;
  //
  while(1) {
    if (path.BeginsWith("path: ")) { // must load from OCDB
      entry = GetCDBEntry(path.Data());
      if (!entry) break;
      resp = (AliITSresponseSDD*) entry->GetObject();
      entry->SetObject(NULL);
      entry->SetOwner(kTRUE);
      break;
    }
    //
    if (gSystem->AccessPathName(path.Data())) break;
    TFile* precf = TFile::Open(path.Data());
    if (precf->FindKey("AliITSresponseSDD")) resp = (AliITSresponseSDD*)precf->Get("AliITSresponseSDD");
    else if (precf->FindKey("AliCDBEntry") && (entry=(AliCDBEntry*)precf->Get("AliCDBEntry"))) {
      resp = (AliITSresponseSDD*) entry->GetObject();
      if (resp && resp->InheritsFrom(AliITSresponseSDD::Class())) entry->SetObject(NULL);
      else resp = 0;
      entry->SetObject(NULL);
      entry->SetOwner(kTRUE);
      delete entry;
    }
    //
    precf->Close();
    delete precf;
    break;
  } 
  //
  if (!resp) {printf("Error: Failed to load SDD response from %s\n",path.Data()); return kFALSE;}
  return kTRUE;
}

//________________________________________________________________________________________________________
Bool_t LoadSDDCorrMap(TString& path, TObjArray *&maps)
{
  // Load SDD correction map
  //
  if (path.IsNull()) return kFALSE;
  printf("Loading SDD Correction Maps from %s\n",path.Data());
  //
  AliCDBEntry *entry = 0;
  delete maps;
  maps = 0;
  while(1) {
    if (path.BeginsWith("path: ")) { // must load from OCDB
      entry = GetCDBEntry(path.Data());
      if (!entry) break;
      maps = (TObjArray*) entry->GetObject();
      entry->SetObject(NULL);
      entry->SetOwner(kTRUE);
      break;
    }
    //
    if (gSystem->AccessPathName(path.Data())) break;
    TFile* precf = TFile::Open(path.Data());
    if (precf->FindKey("TObjArray")) maps = (TObjArray*)precf->Get("TObjArray");
    else if (precf->FindKey("AliCDBEntry") && (entry=(AliCDBEntry*)precf->Get("AliCDBEntry"))) {
      maps = (TObjArray*) entry->GetObject();
      if (maps && maps->InheritsFrom(TObjArray::Class())) entry->SetObject(NULL);
      else maps = 0;
      entry->SetObject(NULL);
      entry->SetOwner(kTRUE);
      delete entry;
    }
    //
    precf->Close();
    delete precf;
    break;
  } 
  //
  if (!maps) {printf("Failed to load SDD Correction Map from %s\n",path.Data()); return kFALSE;}
  
  return kTRUE;
}

//_______________________________________________________________________________________
AliCDBEntry* GetCDBEntry(const char* path)
{
  // return object from the OCDB
  AliCDBEntry *entry = 0;
  printf("Loading object %s\n",path);
  AliCDBManager* man = AliCDBManager::Instance();
  AliCDBId* cdbId = AliCDBId::MakeFromString(path);
  if (!cdbId) {
    printf("Failed to create cdbId\n");
    return 0;
  }
  //
  AliCDBStorage* stor = man->GetDefaultStorage();
  if (!stor && !man->GetRaw()) man->SetDefaultStorage("raw://");
  if (man->GetRaw()) man->SetRun(cdbId->GetFirstRun());
  if (stor) {
    TString tp = stor->GetType();
    if (tp.Contains("alien",TString::kIgnoreCase) && !gGrid) TGrid::Connect("alien:"); 
  } 
  entry = man->Get(cdbId->GetPath(),cdbId->GetFirstRun(),cdbId->GetVersion(),cdbId->GetSubVersion());
  //  entry = man->Get( *cdbId );
  man->ClearCache();
  //
  delete cdbId;
  return entry;
  //
}
//

//_______________________________________________________________________________________
Bool_t PlotHisto(TH1* h, Option_t* opt, int mrkStyle,int mrkCol, double mrkSize)
{
  const double kOffsH = 0.15;
  if (!h) return kFALSE;
  TString opts = opt; opts.ToLower();
  if (opts.Contains("p")) {
    h->SetMarkerStyle(mrkStyle);
    h->SetMarkerColor(mrkCol);
    h->SetMarkerSize(mrkSize);
  }
  h->SetLineColor(mrkCol);
  h->Draw(opt);
  //
  h->SetMinimum(); h->SetMaximum();
  double hmn=h->GetMinimum(),hmx=h->GetMaximum(); // new histo min/max
  //
  TH1* hbase = GetPadBaseHisto((TPad*)gPad); if (!hbase) return 0;
  double smn = hbase->GetMinimum(),smx = hbase->GetMaximum(); // base set min/max?
  hbase->SetMinimum();  hbase->SetMaximum();
  double omn = hbase->GetMinimum(),omx = hbase->GetMaximum(); // base real min max
  if (smn<omn && smx>omx) { // min/max for bas histo was set by hand: extract original min/max
    omx = (smn*kOffsH+smx*(1+kOffsH))/(1+2*kOffsH);
    omn = (smn-kOffsH*omx)/(1+kOffsH);
  }
  if (hmn<omn) omn = hmn;
  if (hmx>omx) omx = hmx;
  double del = omx-omn;
  hbase->SetMinimum( omn - kOffsH*del );
  hbase->SetMaximum( omx + kOffsH*del );
  gPad->Update();
  return kTRUE;
}

//_______________________________________________________________________________________
void GetHistoMinMaxInRange(TH1* h, double &mn,double &mx)
{
  // compute min/max of histo in the user range
  mn = 1e50;
  mx =-1e50;
  int b0 = h->GetXaxis()->GetFirst(), b1 = h->GetXaxis()->GetLast();
  for (int i=b0;i<=b1;i++) {
    double e = h->GetBinError(i); 
    if (TMath::Abs(e)<1e-9) continue;
    double v = h->GetBinContent(i);
    if (mn>v-e) mn = v-e;
    if (mx<v+e) mx = v+e;
  }
}

//_______________________________________________________________________________________
void PlotReport(const char* psname)
{
  // report results
  sddCanv = new TCanvas("sddCanv","sddCanv",700,900);
  //
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //
  TString psnm1 = psname;
  if (psnm1.IsNull()) psnm1 = "sddQAreport.ps";
  TString psnm0 = psnm1 + "["; 
  TString psnm2 = psnm1 + "]";
  sddCanv->Print(psnm0.Data());
  //
  // mean corrections per module/side
  sddCanv->Clear();
  sddCanv->Divide(2,3);
  int cntPad = 0;
  //
  for (int ix=0;ix<2;ix++) { // mean residuals before/after correction
    sddCanv->cd(++cntPad);
    PlotHisto(resOffsDXraw[ix],"p"     ,7,kBlack,1);
    PlotHisto(resOffsDX[ix]   ,"p same",7,kRed  ,1);    
    AddPadLabel(Form("<#DeltaX> %s : Raw",ix?"Right":"Left"), 0.1,0.93,kBlack,0.05);
    AddPadLabel("After Map", 0.5,0.93,kRed,0.05);
  }
  //
  for (int ix=0;ix<2;ix++) { // mean residuals before/after correction
    sddCanv->cd(++cntPad);
    PlotHisto(resVDCorr[ix] ,"p"     ,7,kBlack,1);
    PlotHisto(resVDCorrZ[ix],"p same",7,kRed  ,1);    
    AddPadLabel(Form("<#DeltaV> %s : from #DeltaX vs X",ix?"Right":"Left"), 0.1,0.93,kBlack,0.05);
    AddPadLabel("from #DeltaX vs Z", 0.6,0.93,kRed,0.05);
  }
  //
  sddCanv->cd(++cntPad);
  PlotHisto(resT0Corr,"p",7,kBlue,1);
  AddPadLabel("T0 Correction", 0.3,0.93,kRed,0.05);
  if (forceT0CorrTo0) AddPadLabel("Forced to 0 by transferring to maps", 0.15,0.88,kRed,0.05);
  //
  sddCanv->cd(++cntPad);
  PlotHisto(resXLocCorr,"p",7,kBlack,1);
  AddPadLabel("#DeltaX Correction", 0.3,0.93,kRed,0.05);
  if (forceRPhiCorrTo0) AddPadLabel("Forced to 0 by transferring to maps", 0.15,0.88,kRed,0.05);
  //
  sddCanv->cd();
  sddCanv->Print(psnm1.Data());
  //
  //-------------------------------------------------------------------
  TH1* hsdd = 0;
  //
  cntPad = 999;
  int nModPerPage = 3;
  int nRowPerMod = 2;
  Bool_t saved = kFALSE;
  int ib0=1,ib1=999;
  //
  for (int imd=kSDDMin;imd<=kSDDMax;imd++) {
    if (cntPad>=2*nModPerPage*nRowPerMod) {
      sddCanv->cd();
      if (imd!=kSDDMin) sddCanv->Print(psnm1.Data());
      sddCanv->Clear();
      sddCanv->Divide(2,nModPerPage*nRowPerMod);
      cntPad = 0;
      saved = kTRUE;
    }
    for (int ix=0;ix<2;ix++) {
      sddCanv->cd(++cntPad);
      hsdd = (TH1*)procHistos.At(GetStoreID(kSXvsX,imd,ix)); // raw residuals
      if (hsdd) {
	ib0 = hsdd->FindBin(1);
	ib1 = hsdd->FindBin(sddSeg.Dx()-1);
	hsdd->GetXaxis()->SetRange(ib0,ib1);
      }
      PlotHisto(hsdd," ",7,kBlack,1);
      hsdd = (TH1*)procHistos.At(GetStoreID(kSCorrMapX,imd,ix)); // map
      if (hsdd) hsdd->GetXaxis()->SetRange(ib0,ib1);
      PlotHisto(hsdd,"same",7,kRed,1);
      hsdd = (TH1*)procHistos.At(GetStoreID(kSXvsXCorr,imd,ix)); // map
      if (hsdd) hsdd->GetXaxis()->SetRange(ib0,ib1);
      PlotHisto(hsdd,"same",7,kBlue,1);
      //
      AddPadLabel(Form("<#DeltaX> %d %s: Raw",imd,ix?"Right":"Left"), 0.1,0.93,kBlack,0.07);
      AddPadLabel("Map", 0.4,0.93,kRed,0.07);
      AddPadLabel("Corr.Map", 0.5,0.93,kBlue,0.07);
      //
      AddPadLabel(Form("#DeltaV:%+.4f | #DeltaT0:%+5.0f | #DeltaX:%+4.0f",
		       resVDCorr[ix] ? resVDCorr[ix]->GetBinContent(imd-kSDDMin+1):0,
		       resT0Corr     ?     resT0Corr->GetBinContent(imd-kSDDMin+1):0,
		       resXLocCorr   ?   resXLocCorr->GetBinContent(imd-kSDDMin+1):0),
		  0.5, 0.15, kRed, 0.07);
      //
      saved = kFALSE;
    }
    //
    for (int ix=0;ix<2;ix++) {
      sddCanv->cd(++cntPad);
      hsdd = (TH1*)procHistos.At(GetStoreID(kSDVvsZ,imd,ix)); // correction Vd vs Z
      PlotHisto(hsdd," ",7,kBlack,1);
      hsdd = (TH1*)procHistos.At(GetStoreID(kSDVvsZOrig,imd,ix)); // correction Vd vs Z
      PlotHisto(hsdd,"same",24,kBlue,1);
      AddPadLabel(Form("<#DeltaV> vs Z %d %s | Stat:%.2e",imd,ix?"Right":"Left",
		       ((TH1*)procHistos.At(GetStoreID(kSXvsX,imd,ix)))->GetEntries()), 0.1,0.93,kBlack,0.07);
      //
      AddPadLabel(Form("<#DeltaV>:%+.4f",resVDCorrZ[ix] ? resVDCorrZ[ix]->GetBinContent(imd-kSDDMin+1):0), 0.5, 0.15, kRed, 0.07);
      //
      saved = kFALSE;
    }
    //
  }
  //
  sddCanv->cd();
  if (!saved) sddCanv->Print(psnm1.Data());
  sddCanv->Print(psnm2.Data());
}

//__________________________________
TH1* GetPadBaseHisto(TPad* pad)
{
  if (!pad) pad = (TPad*)gPad;
  if (!pad) return 0;
  TList* lst = pad->GetListOfPrimitives();
  int size = lst->GetSize();
  TH1* hst=0;
  for (int i=0;i<size;i++) {
    TObject* obj = lst->At(i);
    if (!obj) continue;
    if (obj->InheritsFrom("TH1")) {hst = (TH1*)obj; break;}
  }
  return hst;
}

//__________________________________
TLatex* AddPadLabel(const char*txt,float x,float y,int color,float size)
{
  TLatex* lt = new TLatex(x,y,txt); 
  lt->SetNDC(); 
  lt->SetTextColor(color);
  lt->SetTextSize(size);
  lt->Draw();
  return lt;
}

//__________________________________
void SetCDBObjData(int firstrun,int lastrun,const char* comment)
{
  // change range and comment of the objects to store
  firstRun = firstrun;
  lastRun  = lastrun;
  cdbComment = comment;  
}

//__________________________________
void PrepCDBObj(TObject *obj,const char* path,int firstrun,int lastrun,const char* comment)
{
  if (firstrun<0) firstrun = 0;
  //
  AliCDBManager* man = AliCDBManager::Instance();
  man->UnsetDefaultStorage();
  man->SetDefaultStorage("local://");
  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetResponsible("Ruben Shahoyan");
  md->SetComment(comment);
  AliCDBId id(path,firstrun,lastrun<=0 ? (AliCDBRunRange::Infinity()) : lastrun);
  //AliCDBStorage* st = man->GetStorage("local//.");
  man->Put(obj,id,md); 
  //
}

//__________________________________________________________
double edgeLow(double* x, double *par)
{
  // Low TDrift edge:
  // residuals assuming linear dependence of "natural" residual vs Xtrue and smeared
  // by the finite track resolution
  double x0    = par[0];
  double sigma = par[1];
  double offs  = par[2];
  double slop  = par[3];
  //
  if (sigma<1) return 0;
  if (x0<-sigma) return 0; 
  //
  double xex = x[0];
  xex -= x0;
  //
  double arg = xex/sigma;
  arg *= arg/2;
  double res = arg<50 ? slop*sigma*TMath::Exp(-arg)/TMath::Sqrt(2*TMath::Pi()) : 0;
  double erftrm = 1.+TMath::Erf(xex/sigma/TMath::Sqrt(2));
  //printf("%+e %+e %+e\n",x[0],x0,erftrm);
  if (xex<0 && TMath::Abs(erftrm)<1e-10) res = -xex*slop;
  else res /= erftrm/2.;
  res += (offs + (slop-1.)*xex);
  return res;
  //
}

//__________________________________________________________
double edgeHigh(double* x, double *par)
{
  // High TDrift edge
  // residual assuming linear dependence of "natural" residual vs Xtrue and smeared
  // by the finite track resolution
  double x0    = par[0];
  double sigma = par[1];
  double offs  = par[2];
  double slop  = par[3];
  double tailCorr  = par[4];
  //
  if (sigma<1) return 0;
  if (x0<-sigma) return 0; 
  //
  double xex = (x0 - x[0])*tailCorr;
  //
  double arg = xex/sigma;
  arg *= arg/2;
  double res = arg<50 ? slop*sigma*TMath::Exp(-arg)/TMath::Sqrt(2*TMath::Pi()) : 0;
  double erftrm = 1.+TMath::Erf(xex/sigma/TMath::Sqrt(2));
  if (xex<0 && TMath::Abs(erftrm)<1e-10) res = xex*slop;
  else res /= -erftrm/2.;
  res += (offs + (slop-1.)*xex);
  return res;
  //
}

//_________________________________________________________________________
void RedoProfileErrors(TH1* profH, TProfile* prof)
{
  // cure errors of low.stat bins
  int nbCnt = 0, nbn = prof->GetNbinsX();
  double meanStat = 0, meanSpread = 0, wghStat = 0;
  for (int i=1;i<=nbn;i++) {
    double stat = prof->GetBinEntries(i);
    if (stat>0) {meanStat+=stat; nbCnt++;}
  }
  if (nbCnt>0) meanStat/= nbCnt;   // mean occupancy
  //
  for (int i=1;i<=nbn;i++) {
    double stat = prof->GetBinEntries(i);
    if (stat<meanStat/2) continue;
    meanSpread += prof->GetBinError(i)*TMath::Sqrt(stat)*stat;
    wghStat += stat;
  }
  if (wghStat) meanSpread /= wghStat;             // mean spread
  //
  for (int i=1;i<=nbn;i++) {                       // assign error acording to occupancy
    double stat = prof->GetBinEntries(i);
    if (stat>meanStat/2 || stat<1) continue;
    profH->SetBinError(i, meanSpread/TMath::Sqrt(stat));
  }
}

//_________________________________________________________________________
void CureEdges(TH1* prof)
{
  // cure edges of the profile histo
  const double kMaxChi2 = 20.;
  const double kSlpDf = 0.05;
  static TF1* fitEdgeLow  = new TF1("fitEdgeLow" ,edgeLow , -5000, sddSeg.Dx()+5000,4);
  static TF1* fitEdgeHigh = new TF1("fitEdgeHigh",edgeHigh, -5000, sddSeg.Dx()+5000,5);
  //
  int ndf,ib0,ib1,nbn = prof->GetNbinsX();
  double sigma,offs,slp,chi2,x0;
  //
  // LowT edge
  // find 1st non-empty bin
  for (ib0=1;ib0<=nbn;ib0++) if (prof->GetBinError(ib0)>1e-9) break;
  x0 = prof->GetBinCenter(ib0);
  ib1 = prof->FindBin(wDXEdge); 
  if (ib1-ib0<minBinsEdge) ib1 = ib0+minBinsEdge;
  //
  fitEdgeLow->SetParameters(100,100,0,1);
  fitEdgeLow->SetParLimits(0,0, sddSeg.Dx());
  fitEdgeLow->SetParLimits(1,edgeSmearMinErr, edgeSmearMaxErr);
  fitEdgeLow->SetParLimits(3,1.-kSlpDf, 1.+kSlpDf);
  //
  prof->Fit(fitEdgeLow,"q","",prof->GetBinLowEdge(ib0)+1, prof->GetBinCenter(ib1+1)-1);
  chi2 = fitEdgeLow->GetChisquare();
  ndf = fitEdgeLow->GetNDF();
  if (ndf>0) chi2 /= ndf;
  //
  x0    = fitEdgeLow->GetParameter(0);
  sigma = fitEdgeLow->GetParameter(1);
  offs  = fitEdgeLow->GetParameter(2);
  slp   = fitEdgeLow->GetParameter(3);
  if ( chi2<kMaxChi2) { 
    x0 += 3*sigma;
    ib1 = prof->FindBin(x0);
    for (int i=ib0;i<=ib1;i++) {
      if (prof->GetBinError(i)<1e-9) continue;
      double xb = prof->GetBinCenter(i);
      double polval = offs+(slp-1.)*(xb-x0);
      if (xb>0) prof->AddBinContent(i, polval -  fitEdgeLow->Eval( xb ) );
      else prof->SetBinContent(i,polval);
    }
  }
  //
  // find last non-empty bin
  for (ib1=nbn;ib1>=1;ib1--) if (prof->GetBinError(ib1)>1e-9) break;
  x0 = prof->GetBinCenter(ib1);
  ib0 = prof->FindBin(sddSeg.Dx() - wDXEdge); 
  if (ib1-ib0<minBinsEdge) ib0 = ib1-minBinsEdge;
  //
  fitEdgeHigh->SetParameters(prof->GetBinCenter(ib0)+wDXEdge-100,100,0,1,1.);
  fitEdgeHigh->SetParLimits(0,0, sddSeg.Dx()+150);
  fitEdgeHigh->SetParLimits(1,edgeSmearMinErr, edgeSmearMaxErr);
  fitEdgeHigh->SetParLimits(3,1.-kSlpDf, 1.+kSlpDf);
  fitEdgeHigh->SetParLimits(4,0.3, 3.);
  prof->Fit(fitEdgeHigh,"q+","",prof->GetBinLowEdge(ib0)+1, prof->GetBinCenter(ib1+1)+1);
  //
  chi2 = fitEdgeHigh->GetChisquare();
  ndf = fitEdgeHigh->GetNDF();
  if (ndf>0) chi2 /= ndf;
  //
  x0    = fitEdgeHigh->GetParameter(0);
  sigma = fitEdgeHigh->GetParameter(1);
  offs  = fitEdgeHigh->GetParameter(2);
  slp   = fitEdgeHigh->GetParameter(3);
  if ( chi2<kMaxChi2 ) {
    x0 -= 3*sigma;
    ib0 = prof->FindBin(x0);
    for (int i=ib0;i<=ib1;i++) {
      if (prof->GetBinError(i)<1e-9) continue;
      double xb = prof->GetBinCenter(i);
      double polval = offs+(slp-1.)*(xb-x0);
      if (xb<sddSeg.Dx()) prof->AddBinContent(i, polval -  fitEdgeHigh->Eval( xb ) );
      else prof->SetBinContent(i,polval);
    }
  }
  //
}
