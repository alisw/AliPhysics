#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TH3.h"
#include "TH3D.h"
#include "TNtuple.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"

#include "AliHFSystErr.h"
#include <Riostream.h>
#endif

//------------------------------------------------------------------------------------------
//
//  MACRO TO COMPUTE THE AVERAGE D0, D+ and D*+ meson RAA
//
//   Different options considered as weights (relative stat, relative stat+uncorr-syst,
//   relative stat+global-syst, absolute stat) as option of the macro
//
//  Output: merged RAA root file, plus a tex file with the results
//
//  Usage:
//    1. First set the average pt bining (variables nbins and ptbinlimits on top of the macro)
//    2. Load the libraries and compile the macro
//    3. Macro parameters:
//       a: D0 RAA file,  b: D0 pp reference file
//       d: D+ RAA file,  d: D+ pp reference file
//       e: D*+ RAA file, f: D*+ pp reference file
//       g: output filename
//       h: average option (weights)
//       i: centrality range,       j: centrality estimator
//       k: flag in case the pp reference has the split separated uncertainties
//       l,m,n: flag in case D0/D+/D*+ pp reference has some extrapolated bin
//
//
//  Author: Z. Conesa del Valle
//
//------------------------------------------------------------------------------------------


/////
// Compilation instructions
//
//  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/PWGHF/base -I$ALICE_PHYSICS/PWGHF/vertexingHF -I$ALICE_PHYSICS/PWG/FLOW/Base -I$ALICE_PHYSICS/PWG/FLOW/Tasks -I$ALICE_PHYSICS/PWG -g"); 
// .x $ALICE_PHYSICS/PWGHF/vertexingHF/macros/LoadLibraries.C
// .L $ALICE_PHYSICS/PWGHF/vertexingHF/macros/AverageDmesonRaa.C+
//
/////

//const Int_t nbins=10;
//Double_t ptbinlimits[nbins+1]={1.,2.,3.,4.,5.,6.,7.,8.,12.,16.,24.};
//Double_t ptbinlimits[nbins+1]={1.,2.,3.,4.,5.,6.,8.,12.,16.,24.,36.};
const Int_t nbins=7;
Double_t ptbinlimits[nbins+1]={1.,2.,4.,6.,8.,12.,16.,24.};
Bool_t useExtrapPPref=kFALSE;
Bool_t isDebug=true;

//const Int_t ptmaxPPRefData[3] = { 16., 24., 24. };

enum AverageOption{ kRelativeStatUnc=0, kRelativeStatUncorrWoPidSyst, kRelativeStatUncorrWPidSyst, kRelativeStatRawYieldSyst, kRelativeStatGlobalSyst, kAbsoluteStatUnc };

enum centrality{ kpp, k07half, kpPb0100, k010, k1020, k020, k2040, k2030, k3040, k4050, k3050, k5060, k4060, k6080, k4080, k5080, k80100, kpPb020, kpPb2040, kpPb4060, kpPb60100 };
enum centestimator{ kV0M, kV0A, kZNA, kCL1 };

//____________________________________________________________
Int_t FindGraphBin(TGraphAsymmErrors *gr, Double_t pt)
{
    Int_t istart =0;
    Int_t npoints = gr->GetN();
    for(Int_t i=0; i<=npoints; i++){
        Double_t x=0.,y=0.;
        gr->GetPoint(i,x,y);
        if ( TMath::Abs ( x - pt ) < 0.4 ) {
            istart = i;
            break;
        }
    }
    return istart;
}
//____________________________________________________________
void FindGraphRelativeUnc(TGraphAsymmErrors *gr, Double_t pt, Double_t &uncLow, Double_t &uncHigh)
{
    Int_t istart =-1;
    Int_t npoints = gr->GetN();
    Double_t x=0.,y=0.;
    for(Int_t i=0; i<=npoints; i++){
        gr->GetPoint(i,x,y);
        if ( TMath::Abs ( x - pt ) < 0.4 ) {
            istart = i;
            break;
        }
    }
    
    uncLow=0.; uncHigh=0.;
    
    if(istart>-1){
        uncLow = gr->GetErrorYlow(istart)/y;
        uncHigh = gr->GetErrorYhigh(istart)/y;
        //    cout<<"      pt="<<pt<<" -"<<uncLow<<" +"<<uncHigh<<endl;
    }
}

//____________________________________________________________
Double_t GetWeight(Int_t averageoption, Double_t pt,
                   TH1D* hRaa, Double_t raaSystLow, Double_t raaSystHigh,
		   Double_t ppSystRawYield,
                   Double_t ppSystRawYieldCutVar, Double_t ppSystRawYieldCutVarPid,
		   Double_t ABSystRawYield,
                   Double_t ABSystRawYieldCutVar, Double_t ABSystRawYieldCutVarPid)
{
    Double_t weight=1.0;
    Int_t hbin = hRaa->FindBin(pt);
    
    Double_t stat = hRaa->GetBinError(hbin);
    Double_t relativeStat = stat / hRaa->GetBinContent(hbin);
    Double_t weightStat=0.;
    Double_t relativeSyst = 0.;
    if(averageoption==kRelativeStatUnc) {
        weightStat = relativeStat;
    }
    else if(averageoption==kAbsoluteStatUnc) {
        weightStat = stat;
    }
    else if(averageoption==kRelativeStatUncorrWoPidSyst) {
        weightStat = relativeStat;
        relativeSyst = TMath::Sqrt( ppSystRawYieldCutVar*ppSystRawYieldCutVar + ABSystRawYieldCutVar*ABSystRawYieldCutVar );
    }
    else if(averageoption==kRelativeStatUncorrWPidSyst) {
        weightStat = relativeStat;
        relativeSyst = TMath::Sqrt( ppSystRawYieldCutVarPid*ppSystRawYieldCutVarPid + ABSystRawYieldCutVar*ABSystRawYieldCutVarPid );
    }
    else if(averageoption==kRelativeStatRawYieldSyst){
        weightStat = relativeStat;
        relativeSyst = TMath::Sqrt( ppSystRawYield*ppSystRawYield + ABSystRawYield*ABSystRawYield );      
    }
    else if(averageoption==kRelativeStatGlobalSyst) {
        weightStat = relativeStat;
        relativeSyst = raaSystHigh>raaSystLow ? raaSystHigh : raaSystLow;
    }

    //  weight = TMath::Sqrt( relativeStat*relativeStat + relativeSyst*relativeSyst );
    weight = TMath::Sqrt( weightStat*weightStat + relativeSyst*relativeSyst );
    // cout<< endl<<" rel stat "<< relativeStat<<" rel syst "<< relativeSyst<<" weight="<<(1.0/(weight*weight))<<endl;
    
    return (1.0/(weight*weight));
}


//____________________________________________________________
void AverageDmesonRaa( const char* fD0Raa="",    const char* fD0ppRef="",
                      const char* fDplusRaa="", const char* fDplusppRef="",
                      const char* fDstarRaa="", const char* fDstarppRef="",
                      const char* outfile="", Int_t averageOption=kRelativeStatUnc, Int_t cc=kpPb0100, Int_t ccestimator=kV0M,
                      Bool_t isReadAllPPUnc=false, Bool_t isPPRefExtrapD0=false, Bool_t isPPRefExtrapDplus=false, Bool_t isPPRefExtrapDstar=false)
{
    
    FILE *resultFile;
    TString foutname = "Average";
    if(fD0Raa) foutname += "Dzero";
    if(fDplusRaa) foutname += "Dplus";
    if(fDstarRaa) foutname += "Dstar";
    if(averageOption==kRelativeStatUnc)  foutname+= "_RelStatUncWeight";
    else if(averageOption==kAbsoluteStatUnc)  foutname+= "_AbsStatUncWeight";
    else if(averageOption==kRelativeStatUncorrWoPidSyst) foutname+= "_RelStatUncorrWeight";
    else if(averageOption==kRelativeStatRawYieldSyst) foutname+= "_RelStatRawYieldSystWeight";
    if(!useExtrapPPref) foutname+= "_NoExtrapBins";
    TDatime d;
    TString ndate = Form("%02d%02d%04d",d.GetDay(),d.GetMonth(),d.GetYear());
    resultFile = fopen( Form("%s_result_%s.txt",foutname.Data(),ndate.Data()),"w");
    fprintf(resultFile,"Ptmin (GeV/c)   Ptmax (GeV/c)   Raa(Daverage)   +-(stat)    +(syst)  - (syst) \n\n");
    
    
    // Define here all needed histograms/graphs to be retrieved
    TH1D *hDmesonRaa[3];
    TH1I *hCombinedReferenceFlag[3];
    TGraphAsymmErrors *gDataSystematicsPP[3], *gDataSystematicsAB[3];
    TGraphAsymmErrors *gScalingUncPP[3];
    TGraphAsymmErrors *gRABFeedDownSystematicsElossHypothesis[3];
    TGraphAsymmErrors *gRAB_GlobalSystematics[3];
    TH1D *hDmesonPPRef[3], *hDmesonPPYieldExtrUnc[3], *hDmesonPPCutVarUnc[3], *hDmesonPPIDPUnc[3], *hDmesonPPMCPtUnc[3];
    
    // Define here all output histograms/graphs
    TH1D *hDmesonAverageRAB = new TH1D("hDmesonAverageRAB","D meson average Raa ; p_{T} (GeV/c)",nbins,ptbinlimits);
    TGraphAsymmErrors *gRABNorm;
    TGraphAsymmErrors *gRAB_DmesonAverage_GlobalSystematics = new TGraphAsymmErrors(0);
    TGraphAsymmErrors *gRAB_DmesonAverage_FeedDownSystematicsElossHypothesis = new TGraphAsymmErrors(0);
    TGraphAsymmErrors *gRAB_DmesonAverage_ScalingSystematicsPP = new TGraphAsymmErrors(0);
    TGraphAsymmErrors *gRAB_DmesonAverage_DataSystematicsPP = new TGraphAsymmErrors(0);
    TGraphAsymmErrors *gRAB_DmesonAverage_DataSystematicsAB = new TGraphAsymmErrors(0);
    TGraphAsymmErrors *gRAB_DmesonAverage_TrackingSystematicsPP = new TGraphAsymmErrors(0);
    TGraphAsymmErrors *gRAB_DmesonAverage_TrackingSystematicsAB = new TGraphAsymmErrors(0);
    gRAB_DmesonAverage_GlobalSystematics->SetNameTitle("gRAB_DmesonAverage_GlobalSystematics","DmesonAverage GlobalSystematics");
    gRAB_DmesonAverage_FeedDownSystematicsElossHypothesis->SetNameTitle("gRAB_DmesonAverage_FeedDownSystematicsElossHypothesis","DmesonAverage FeedDownSystematicsElossHypothesis");
    gRAB_DmesonAverage_ScalingSystematicsPP->SetNameTitle("gRAB_DmesonAverage_ScalingSystematicsPP","DmesonAverage Scaling uncertainty PP");
    gRAB_DmesonAverage_DataSystematicsPP->SetNameTitle("gRAB_DmesonAverage_DataSystematicsPP","DmesonRaaAverage DataSystematicsPP - tracking uncertainty PP");
    gRAB_DmesonAverage_DataSystematicsAB->SetNameTitle("gRAB_DmesonAverage_DataSystematicsAB","DmesonRaaAverage DataSystematicsAB - tracking uncertainty AB");
    gRAB_DmesonAverage_TrackingSystematicsPP->SetNameTitle("gRAB_DmesonAverage_TrackingSystematicsPP","DmesonRaaAverage tracking uncertainty PP");
    gRAB_DmesonAverage_TrackingSystematicsAB->SetNameTitle("gRAB_DmesonAverage_TrackingSystematicsAB","DmesonRaaAverage tracking uncertainty AB");
    
    
    const char *filenamesRaa[3] = { fD0Raa, fDplusRaa, fDstarRaa };
    const char *filenamesReference[3] = { fD0ppRef, fDplusppRef, fDstarppRef };
    const char *filenamesSuffixes[3] = { "Dzero", "Dplus", "Dstar" };
    Bool_t isDmeson[3] = { true, true, true };
    
    if(strcmp(filenamesRaa[0],"")==0) { cout<<" Dzero not set, error"<<endl; return; }
    
    //
    // Get Raa file histos and graphs
    //
    AliHFSystErr *ppSyst[3];
    AliHFSystErr *ABSyst[3];
    Double_t ppTracking[3][nbins], ppSystRawYield[3][nbins], ppSystCutVar[3][nbins], ppSystPid[3][nbins];
    Double_t ABTracking[3][nbins], ABSystRawYield[3][nbins], ABSystCutVar[3][nbins], ABSystPid[3][nbins];
    for(Int_t j=0; j<3; j++) {
        if(strcmp(filenamesRaa[j],"")==0)  { isDmeson[j]=false; continue; }
        cout<<" Reading file "<<filenamesRaa[j]<<"..."<<endl;
        TFile *fDRaa = TFile::Open(filenamesRaa[j],"read");
        if(!fDRaa){ cout<<" Error on file !!!!!"<<filenamesRaa[j]<<endl; return; }
        hDmesonRaa[j] = (TH1D*)fDRaa->Get("hRABvsPt");
        hDmesonRaa[j]->SetName(Form("%s_%s",hDmesonRaa[j]->GetName(),filenamesSuffixes[j]));
        
        //    cout<< hDmesonRaa[j]<< " bins="<< hDmesonRaa[j]->GetNbinsX()<<", value bin 1 ="<<hDmesonRaa[j]->GetBinContent(1)<<endl;
        
        if(j==0) gRABNorm = (TGraphAsymmErrors*)fDRaa->Get("gRAB_Norm");
        //
        gDataSystematicsPP[j] = (TGraphAsymmErrors*)fDRaa->Get("gRAB_DataSystematicsPP");
        gDataSystematicsPP[j]->SetName(Form("%s_%s",gDataSystematicsPP[j]->GetName(),filenamesSuffixes[j]));
        gDataSystematicsAB[j] = (TGraphAsymmErrors*)fDRaa->Get("gRAB_DataSystematicsAB");
        gDataSystematicsAB[j]->SetName(Form("%s_%s",gDataSystematicsAB[j]->GetName(),filenamesSuffixes[j]));
        gRABFeedDownSystematicsElossHypothesis[j] = (TGraphAsymmErrors*)fDRaa->Get("gRAB_FeedDownSystematicsElossHypothesis");
        gRABFeedDownSystematicsElossHypothesis[j]->SetName(Form("%s_%s",gRABFeedDownSystematicsElossHypothesis[j]->GetName(),filenamesSuffixes[j]));
        gRAB_GlobalSystematics[j] = (TGraphAsymmErrors*)fDRaa->Get("gRAB_GlobalSystematics");
        gRAB_GlobalSystematics[j]->SetName(Form("%s_%s",gRAB_GlobalSystematics[j]->GetName(),filenamesSuffixes[j]));
	Bool_t shouldDelete=kFALSE;
	if(fDRaa->Get("AliHFSystErrPP")){
	  ppSyst[j]=(AliHFSystErr*)fDRaa->Get("AliHFSystErrPP");
	  printf("AliHFSystErr object for meson %d in pp (%s) read from HFPtSpectrumRaa output file\n",j,ppSyst[j]->GetTitle());
	}else{   
	  printf("Create instance of AliHFSystErr for meson %d in pp \n",j);
	  ppSyst[j] = new AliHFSystErr(Form("ppSyst_%d",j),Form("ppSyst_%d",j));
	  ppSyst[j]->SetIsPass4Analysis(kTRUE);
	  ppSyst[j]->Init(j+1);
	  shouldDelete=kTRUE;
	}
        for(Int_t ipt=0; ipt<nbins; ipt++) {
	  Double_t ptval = ptbinlimits[ipt] + (ptbinlimits[ipt+1]-ptbinlimits[ipt])/2.;
	  ppTracking[j][ipt]=0.; ppSystRawYield[j][ipt]=0; ppSystCutVar[j][ipt]=0; ppSystPid[j][ipt]=0.;
	  ppTracking[j][ipt]    =ppSyst[j]->GetTrackingEffErr(ptval);
	  ppSystRawYield[j][ipt]=ppSyst[j]->GetRawYieldErr(ptval);
	  ppSystCutVar[j][ipt]  =ppSyst[j]->GetCutsEffErr(ptval);
	  ppSystPid[j][ipt]     =ppSyst[j]->GetPIDEffErr(ptval);
        }
	if(shouldDelete)  delete ppSyst[j];
	shouldDelete=kFALSE;
	if(fDRaa->Get("AliHFSystErrAA")){
	  ABSyst[j]=(AliHFSystErr*)fDRaa->Get("AliHFSystErrAA");
	  printf("AliHFSystErr object for meson %d in AA (%s) read from HFPtSpectrumRaa output file\n",j,ppSyst[j]->GetTitle());
	}else{
	  printf("Create instance of AliHFSystErr for meson %d in AA \n",j);
	  ABSyst[j] = new AliHFSystErr(Form("ABSyst_%d",j),Form("ABSyst_%d",j));
	  ABSyst[j]->SetCollisionType(1); // PbPb by default
	  if ( cc == k010 ) ABSyst[j]->SetCentrality("010");
	  else if ( cc == k1020 ) ABSyst[j]->SetCentrality("1020");
	  else if ( cc == k2040 || cc == k2030 || cc == k3040 ) {
            ABSyst[j]->SetCentrality("2040");
            ABSyst[j]->SetIsPbPb2010EnergyScan(true);
	  }
	  else if ( cc == k3050 ) ABSyst[j]->SetCentrality("3050");
	  else if ( cc == k4060 || cc == k4050 || cc == k5060 ) ABSyst[j]->SetCentrality("4060");
	  else if ( cc == k6080 || cc == k5080 ) ABSyst[j]->SetCentrality("6080");
	  else if ( cc == k4080 ) ABSyst[j]->SetCentrality("4080");
	  // Going to pPb systematics
	  else if ( cc == kpPb0100 || cc == kpPb020 || cc == kpPb2040 || cc == kpPb4060 || cc == kpPb60100 ) {
            ABSyst[j]->SetCollisionType(2);
            ABSyst[j]->SetRunNumber(16);
	    if(ccestimator==kV0A) {
	      if(cc == kpPb020) ABSyst[j]->SetCentrality("020V0A");
	      else if(cc == kpPb2040) ABSyst[j]->SetCentrality("2040V0A");
	      else if(cc == kpPb4060) ABSyst[j]->SetCentrality("4060V0A");
	      else if(cc == kpPb60100) ABSyst[j]->SetCentrality("60100V0A");
            } else if (ccestimator==kZNA) {
	      if(cc == kpPb020) ABSyst[j]->SetCentrality("020ZNA");
	      else if(cc == kpPb2040) ABSyst[j]->SetCentrality("2040ZNA");
	      else if(cc == kpPb4060) ABSyst[j]->SetCentrality("4060ZNA");
	      else if(cc == kpPb60100) ABSyst[j]->SetCentrality("60100ZNA");
            } else if (ccestimator==kCL1) {
	      if(cc == kpPb020) ABSyst[j]->SetCentrality("020CL1");
	      else if(cc == kpPb2040) ABSyst[j]->SetCentrality("2040CL1");
	      else if(cc == kpPb4060) ABSyst[j]->SetCentrality("4060CL1");
	      else if(cc == kpPb60100) ABSyst[j]->SetCentrality("60100CL1");
            } else {
	      if(!(cc == kpPb0100)) {
		cout <<" Error on the pPb options"<<endl;
		return;
	      }
            }
	  }
	  ABSyst[j]->Init(j+1);
 	  shouldDelete=kTRUE;
	}
	for(Int_t ipt=0; ipt<nbins; ipt++) {
	  Double_t ptval = ptbinlimits[ipt] + (ptbinlimits[ipt+1]-ptbinlimits[ipt])/2.;
	  ABTracking[j][ipt]=0.; ABSystRawYield[j][ipt]=0; ABSystCutVar[j][ipt]=0; ABSystPid[j][ipt]=0.;
	  ABTracking[j][ipt]    =ABSyst[j]->GetTrackingEffErr(ptval);
	  ABSystRawYield[j][ipt]=ABSyst[j]->GetRawYieldErr(ptval);
	  ABSystCutVar[j][ipt]  =ABSyst[j]->GetCutsEffErr(ptval);
	  ABSystPid[j][ipt]     =ABSyst[j]->GetPIDEffErr(ptval);
        }
        if(shouldDelete)  delete ABSyst[j];    
    }

    //
    // Get pp-reference file histos and graphs
    //
    const char *pprefhgnames[6] = { "fhScaledData","gScaledDataSystExtrap",
        "fhScaledSystRebinYieldExtraction","fhScaledSystRebinCutVariation","fhScaledSystRebinPIDUnc","fhScaledSystRebinMCPt"};
    Bool_t isPPRefExtrap[3] = { isPPRefExtrapD0, isPPRefExtrapDplus, isPPRefExtrapDstar };
    for(Int_t j=0; j<3; j++) {
        if(strcmp(filenamesReference[j],"")==0)  { isDmeson[j]=false; continue; }
        cout<<" Reading file "<<filenamesReference[j]<<"..."<<endl;
        TFile *fRef = TFile::Open(filenamesReference[j],"read");
        gScalingUncPP[j] = (TGraphAsymmErrors*)fRef->Get(pprefhgnames[1]);
        gScalingUncPP[j]->SetName(Form("%s_%s",gScalingUncPP[j]->GetName(),filenamesSuffixes[j]));
        if(isPPRefExtrap[j]) {
            hCombinedReferenceFlag[j] = (TH1I*)fRef->Get("hCombinedReferenceFlag");
            hCombinedReferenceFlag[j]->SetName(Form("%s_%s",hCombinedReferenceFlag[j]->GetName(),filenamesSuffixes[j]));
        }
        if(isReadAllPPUnc){
            const char*hname="fhScaledData";
            if(isPPRefExtrap[j]) hname="hReference";
            hDmesonPPRef[j] = (TH1D*)fRef->Get(hname);
            hDmesonPPRef[j]->SetName(Form("%s_%s",hDmesonPPRef[j]->GetName(),filenamesSuffixes[j]));
            hDmesonPPYieldExtrUnc[j] = (TH1D*)fRef->Get(pprefhgnames[2]);
            hDmesonPPYieldExtrUnc[j]->SetName(Form("%s_%s",hDmesonPPYieldExtrUnc[j]->GetName(),filenamesSuffixes[j]));
            hDmesonPPCutVarUnc[j] = (TH1D*)fRef->Get(pprefhgnames[3]);
            hDmesonPPCutVarUnc[j]->SetName(Form("%s_%s",hDmesonPPCutVarUnc[j]->GetName(),filenamesSuffixes[j]));
            hDmesonPPIDPUnc[j] = (TH1D*)fRef->Get(pprefhgnames[4]);
            hDmesonPPIDPUnc[j]->SetName(Form("%s_%s",hDmesonPPIDPUnc[j]->GetName(),filenamesSuffixes[j]));
            hDmesonPPMCPtUnc[j] = (TH1D*)fRef->Get(pprefhgnames[5]);
            hDmesonPPMCPtUnc[j]->SetName(Form("%s_%s",hDmesonPPMCPtUnc[j]->GetName(),filenamesSuffixes[j]));
        }
    }
    
    
    //
    // Loop per pt bin
    //
    for(Int_t ipt=0; ipt<nbins; ipt++) {
        
        cout<<" Calculation for pt bin ("<<ptbinlimits[ipt]<<","<<ptbinlimits[ipt+1]<<")"<<endl;
        
        Double_t ptval = ptbinlimits[ipt] + (ptbinlimits[ipt+1]-ptbinlimits[ipt])/2.;
        Double_t RaaDmeson[3]={0.,0.,0.};
        Double_t RaaDmesonStat[3]={0.,0.,0.};
        Double_t RaaDmesonSystLow[3]={0.,0.,0.};
        Double_t RaaDmesonSystHigh[3]={0.,0.,0.};
        Double_t weight[3]={0.,0.,0.};
        Double_t ppSystLow[3]={0.,0.,0.};
        Double_t ppSystHigh[3]={0.,0.,0.};
        Double_t ppSystUncorrLow[3]={0.,0.,0.};
        Double_t ppSystUncorrHigh[3]={0.,0.,0.};
        //    Double_t ppTracking[3]={0.,0.,0.};
        Double_t ScalingLow[3]={0.,0.,0.};
        Double_t ScalingHigh[3]={0.,0.,0.};
	Double_t ppSystRawYieldOnly[3]={0.,0.,0.};
	Double_t ppSystRawYieldCutVar[3]={0.,0.,0.};
        Double_t ppSystRawYieldCutVarPid[3]={0.,0.,0.};
        Double_t ABSystLow[3]={0.,0.,0.};
        Double_t ABSystHigh[3]={0.,0.,0.};
        Double_t ABSystUncorrLow[3]={0.,0.,0.};
        Double_t ABSystUncorrHigh[3]={0.,0.,0.};
        Double_t ABSystRawYieldOnly[3]={0.,0.,0.};
        Double_t ABSystRawYieldCutVar[3]={0.,0.,0.};
        Double_t ABSystRawYieldCutVarPid[3]={0.,0.,0.};
        //    Double_t ABTracking[3]={0.,0.,0.};
        Double_t RabFdElossLow[3]={0.,0.,0.};
        Double_t RabFdElossHigh[3]={0.,0.,0.};
        Double_t RabGlobalLow[3]={0.,0.,0.};
        Double_t RabGlobalHigh[3]={0.,0.,0.};
        
        Double_t average=0., averageStat=0.;
        Double_t weightTot=0.;
        Double_t ppTrackingAv=0., ABTrackingAv=0.;
        Double_t ppDataSystAvLow=0., ppDataSystAvHigh=0.;
        Double_t ABDataSystAvLow=0., ABDataSystAvHigh=0.;
        Double_t scalingLowAv=0., scalingHighAv=0.;
        Double_t raaSystUncorrLow=0., raaSystUncorrHigh=0.;
        Double_t raabeautyLow=0., raabeautyHigh=0.;
        
        Int_t histoBin=-1;
        
        // Get tracking uncertainties and raw yield and cut-variation and pid-systematics
        if(isDebug) cout<<" Retrieving tracking + rawyield systematics"<<endl;
        for(Int_t j=0; j<3; j++) {
            if(!isDmeson[j]) continue;
	    ppSystRawYieldOnly[j] = ppSystRawYield[j][ipt];
            ppSystRawYieldCutVar[j] = TMath::Sqrt( ppSystRawYield[j][ipt]*ppSystRawYield[j][ipt]
                                                  + ppSystCutVar[j][ipt]*ppSystCutVar[j][ipt] );
            ppSystRawYieldCutVarPid[j] = TMath::Sqrt( ppSystRawYield[j][ipt]*ppSystRawYield[j][ipt]
                                                     + ppSystCutVar[j][ipt]*ppSystCutVar[j][ipt]
                                                     + ppSystPid[j][ipt]*ppSystPid[j][ipt] );
	    ABSystRawYieldOnly[j] = ABSystRawYield[j][ipt];
            ABSystRawYieldCutVar[j] = TMath::Sqrt( ABSystRawYield[j][ipt]*ABSystRawYield[j][ipt]
                                                  + ABSystCutVar[j][ipt]*ABSystCutVar[j][ipt] );
            ABSystRawYieldCutVarPid[j] = TMath::Sqrt( ABSystRawYield[j][ipt]*ABSystRawYield[j][ipt]
                                                     + ABSystCutVar[j][ipt]*ABSystCutVar[j][ipt]
                                                     + ABSystPid[j][ipt]*ABSystPid[j][ipt] );
            if(isDebug) cout<<" j="<<j<<" pt="<< ptval<<" ppref unc RY+CV="<<ppSystRawYieldCutVar[j]<<" RY+CV+PID="<<ppSystRawYieldCutVarPid[j]<<endl;
            if(isDebug) cout<<" j="<<j<<" pt="<< ptval<<" AB unc RY+CV="<<ABSystRawYieldCutVar[j]<<" RY+CV+PID="<<ABSystRawYieldCutVarPid[j]<<endl;
        }
        
        if(isReadAllPPUnc){
            if(isDebug) cout<<" Retrieving all pp reference systematics from the rebinned file"<<endl;
            for(Int_t j=0; j<3; j++) {
                if(!isDmeson[j]) continue;
                Int_t ibin = hDmesonPPRef[j]->FindBin(ptval);
                Double_t ppval = hDmesonPPRef[j]->GetBinContent(ibin);
                Double_t rawyield = hDmesonPPYieldExtrUnc[j]->GetBinContent( hDmesonPPYieldExtrUnc[j]->FindBin(ptval) )/ppval;
                Double_t cutvar = hDmesonPPCutVarUnc[j]->GetBinContent( hDmesonPPCutVarUnc[j]->FindBin(ptval) )/ppval;
                Double_t pid = hDmesonPPIDPUnc[j]->GetBinContent(  hDmesonPPIDPUnc[j]->FindBin(ptval) )/ppval;
                ppSystRawYieldCutVar[j] = TMath::Sqrt( rawyield*rawyield + cutvar*cutvar );
                ppSystRawYieldCutVarPid[j] = TMath::Sqrt( rawyield*rawyield + cutvar*cutvar + pid*pid );
                if(isDebug) cout<<"redo j="<<j<<" pt="<< ptval<<" ppref unc RY+CV="<<ppSystRawYieldCutVar[j]<<" RY+CV+PID="<<ppSystRawYieldCutVarPid[j]<<endl;
            }
        }
        // Check for the pp reference systematics for extrapolated pt bins
        for(Int_t j=0; j<3; j++) {
            if(isPPRefExtrap[j]){
                //      if(ptval>ptmaxPPRefData[j]) {
                Int_t ippbin = hCombinedReferenceFlag[j]->FindBin(ptval);
                Bool_t flag = hCombinedReferenceFlag[j]->GetBinContent(ippbin);
                if(!flag) continue;
                //	cout<<" pp ref flag="<<flag<<" >>> ";
                // Get pp reference relative systematics
                Double_t ppSystTotLow=0., ppSystTotHigh=0.;
                FindGraphRelativeUnc(gDataSystematicsPP[j],ptval,ppSystTotLow,ppSystTotHigh);
                ppSystRawYieldCutVar[j] = ppSystTotLow > ppSystTotHigh ? ppSystTotLow : ppSystTotHigh ;
                ppSystRawYieldCutVarPid[j] = ppSystRawYieldCutVar[j];
		ppSystRawYieldOnly[j] = ppSystRawYieldCutVar[j];
	    }
        }
        //
        // Loop per meson to get the Raa values and uncertainties for the given pt bin
        //
        if(isDebug) cout<<" Retrieving all Raa values and uncertainties"<<endl;
        for(Int_t j=0; j<3; j++) {
            // Get value, stat unc and weight
            if(!isDmeson[j]) continue;
            if(!hDmesonRaa[j]) continue;
            histoBin = hDmesonRaa[j]->FindBin(ptval);
            RaaDmeson[j] = hDmesonRaa[j]->GetBinContent( histoBin );
            if (RaaDmeson[j]<=0) continue;
            RaaDmesonStat[j] = hDmesonRaa[j]->GetBinError( histoBin );
            // Get global systematics
            FindGraphRelativeUnc(gRAB_GlobalSystematics[j],ptval,RaaDmesonSystLow[j],RaaDmesonSystHigh[j]);
            // Evaluate the weight
            weight[j] = GetWeight(averageOption,ptval,hDmesonRaa[j],
                                  RaaDmesonSystLow[j],RaaDmesonSystHigh[j],
                                  ppSystRawYieldOnly[j],ppSystRawYieldCutVar[j],ppSystRawYieldCutVarPid[j],
				  ABSystRawYieldOnly[j],ABSystRawYieldCutVar[j],ABSystRawYieldCutVarPid[j]);
            cout<<" raa "<<filenamesSuffixes[j]<<" meson  = "<<RaaDmeson[j]<<" +-"<<RaaDmesonStat[j]<<"(stat) -> (weight="<<weight[j]<<") ,";
            // Get pp reference relative systematics
            FindGraphRelativeUnc(gDataSystematicsPP[j],ptval,ppSystLow[j],ppSystHigh[j]);
            // Get pp-extrapolation relative uncertainty
            FindGraphRelativeUnc(gScalingUncPP[j],ptval,ScalingHigh[j],ScalingLow[j]); // exchanging low-high bc has oposite influence on Raa
            if(isPPRefExtrap[j]){
                Int_t ippbin = hCombinedReferenceFlag[j]->FindBin(ptval);
                Bool_t flag = hCombinedReferenceFlag[j]->GetBinContent(ippbin);
                if(isDebug) cout<< " bin="<<j<<" pp ref flag on? "<<flag;
                if(flag){ ScalingHigh[j]=0.; ScalingLow[j]=0.; ppTracking[j][ipt]=0.; }
		if(flag && !useExtrapPPref){ 
		  weight[j] =0;
		  cout<<"weight set to 0";
		}
            }
            // Get pp reference systematics minus tracking systematics minus extrapolation uncertainties
            ppSystUncorrLow[j] = TMath::Sqrt( ppSystLow[j]*ppSystLow[j]
                                             - ScalingLow[j]*ScalingLow[j]
                                             - ppTracking[j][ipt]*ppTracking[j][ipt] );
            ppSystUncorrHigh[j] = TMath::Sqrt( ppSystHigh[j]*ppSystHigh[j]
                                              - ScalingHigh[j]*ScalingHigh[j]
                                              - ppTracking[j][ipt]*ppTracking[j][ipt] );
            // Get AB relative systematics
            FindGraphRelativeUnc(gDataSystematicsAB[j],ptval,ABSystLow[j],ABSystHigh[j]);
            // Get AB relative systematics minus tracking systematics
            ABSystUncorrLow[j] = TMath::Sqrt( ABSystLow[j]*ABSystLow[j]
                                             - ABTracking[j][ipt]*ABTracking[j][ipt] );
            ABSystUncorrHigh[j] = TMath::Sqrt( ABSystHigh[j]*ABSystHigh[j]
                                              - ABTracking[j][ipt]*ABTracking[j][ipt] );
            // Get Feed-Down and Eloss relative uncertainties on the Raa
            FindGraphRelativeUnc(gRABFeedDownSystematicsElossHypothesis[j],ptval,RabFdElossLow[j],RabFdElossHigh[j]);
            //
            // Check with global Raa uncertainties
            FindGraphRelativeUnc(gRAB_GlobalSystematics[j],ptval,RabGlobalLow[j],RabGlobalHigh[j]);
            Double_t testLow = TMath::Sqrt( RabFdElossLow[j]*RabFdElossLow[j] + ABSystLow[j]*ABSystLow[j]
                                           + ppSystLow[j]*ppSystLow[j] + ScalingLow[j]*ScalingLow[j] );
            Double_t testHigh = TMath::Sqrt( RabFdElossHigh[j]*RabFdElossHigh[j] + ABSystHigh[j]*ABSystHigh[j]
                                            + ppSystHigh[j]*ppSystHigh[j] + ScalingHigh[j]*ScalingHigh[j] );
            if (TMath::Abs( testLow - RabGlobalLow[j] ) > 0.015) {
                cout << endl<<" >>>> Error on the global Raa uncertainties low : test-sum = "<< testLow<<", global = "<< RabGlobalLow[j]<<" ppref="<<ppSystLow[j]<<endl;
            }
            if (TMath::Abs( testHigh - RabGlobalHigh[j] ) > 0.015) {
                cout << endl<<" >>>> Error on the global Raa uncertainties high : test-sum = "<< testHigh<<", global = "<< RabGlobalHigh[j]<<" ppref="<<ppSystHigh[j]<<endl<<endl;
            }
            //
            histoBin = -1;
	    cout<<endl;
        }
        cout<<endl;
        
        //
        // Evaluate Dmeson average
        //
        if(isDebug) cout<<" Evaluating the average"<<endl;
        for(Int_t j=0; j<3; j++){
            if(!isDmeson[j]) continue;
            if( !(RaaDmeson[j]>0.) ) continue;
            weightTot += weight[j];
            // weighted average
            average += RaaDmeson[j]*weight[j];
            // stat absolute uncertainty (uncorrelated) : sum in quadrature
            averageStat += (RaaDmesonStat[j]*weight[j])*(RaaDmesonStat[j]*weight[j]);
            // pp tracking relative uncertainty (correlated) : linear sum
            ppTrackingAv += ppTracking[j][ipt]*RaaDmeson[j]*weight[j];
            // AB tracking relative uncertainty (correlated) : linear sum
            ABTrackingAv += ABTracking[j][ipt]*RaaDmeson[j]*weight[j];
            // pp scaling relative uncertainty (correlated) : linear sum
            scalingLowAv += ScalingLow[j]*RaaDmeson[j]*weight[j];
            scalingHighAv += ScalingHigh[j]*RaaDmeson[j]*weight[j];
            // Get pp and AB relative uncorrelated systematics : sum in quadrature
            raaSystUncorrLow += (ppSystUncorrLow[j]*RaaDmeson[j]*weight[j])*(ppSystUncorrLow[j]*RaaDmeson[j]*weight[j])
            + (ABSystUncorrLow[j]*RaaDmeson[j]*weight[j])*(ABSystUncorrLow[j]*RaaDmeson[j]*weight[j]);
            raaSystUncorrHigh += (ppSystUncorrHigh[j]*RaaDmeson[j]*weight[j])*(ppSystUncorrHigh[j]*RaaDmeson[j]*weight[j])
            + (ABSystUncorrHigh[j]*RaaDmeson[j]*weight[j])*(ABSystUncorrHigh[j]*RaaDmeson[j]*weight[j]);
            ppDataSystAvLow += (ppSystUncorrLow[j]*RaaDmeson[j]*weight[j])*(ppSystUncorrLow[j]*RaaDmeson[j]*weight[j]);
            ABDataSystAvLow += (ABSystUncorrLow[j]*RaaDmeson[j]*weight[j])*(ABSystUncorrLow[j]*RaaDmeson[j]*weight[j]);
            ppDataSystAvHigh += (ppSystUncorrHigh[j]*RaaDmeson[j]*weight[j])*(ppSystUncorrHigh[j]*RaaDmeson[j]*weight[j]);
            ABDataSystAvHigh += (ABSystUncorrHigh[j]*RaaDmeson[j]*weight[j])*(ABSystUncorrHigh[j]*RaaDmeson[j]*weight[j]);
            // Beauty uncertainties: evaluate Raa average for the upper / lower bands
            raabeautyLow += (1-RabFdElossLow[j])*RaaDmeson[j]*weight[j];
            raabeautyHigh += (1+RabFdElossHigh[j])*RaaDmeson[j]*weight[j];
        }
        average /= weightTot;
        averageStat = TMath::Sqrt(averageStat)/weightTot;
        ppTrackingAv /= weightTot;
        ABTrackingAv /= weightTot;
        scalingLowAv /= weightTot;
        scalingHighAv /= weightTot;
        raaSystUncorrLow = TMath::Sqrt(raaSystUncorrLow)/weightTot;
        raaSystUncorrHigh = TMath::Sqrt(raaSystUncorrHigh)/weightTot;
        ppDataSystAvLow = TMath::Sqrt(ppDataSystAvLow)/weightTot;
        ppDataSystAvHigh = TMath::Sqrt(ppDataSystAvHigh)/weightTot;
        ABDataSystAvLow = TMath::Sqrt(ABDataSystAvLow)/weightTot;
        ABDataSystAvHigh = TMath::Sqrt(ABDataSystAvHigh)/weightTot;
        
        // finalization beauty uncertainties
        raabeautyLow /= weightTot;
        raabeautyHigh /= weightTot;
        Double_t RaaBeauty[3] = { average, raabeautyLow, raabeautyHigh };
        Double_t beautyUncLow = average-TMath::MinElement(3,RaaBeauty);
        Double_t beautyUncHigh = TMath::MaxElement(3,RaaBeauty)-average;
        
        // finalization global uncertainties
        Double_t totalUncLow = TMath::Sqrt( raaSystUncorrLow*raaSystUncorrLow
                                           + ppTrackingAv*ppTrackingAv + ABTrackingAv*ABTrackingAv
                                           + scalingLowAv*scalingLowAv
                                           + beautyUncLow*beautyUncLow );
        Double_t totalUncHigh = TMath::Sqrt( raaSystUncorrHigh*raaSystUncorrHigh
                                            + ppTrackingAv*ppTrackingAv + ABTrackingAv*ABTrackingAv
                                            + scalingHighAv*scalingHighAv
                                            + beautyUncHigh*beautyUncHigh );
        if(isDebug) cout<<" Raa="<<average<<" +-"<<averageStat<<"(stat) +"<<ppDataSystAvHigh<<" -"<<ppDataSystAvLow<<" (pp-data) +-"<<ppTrackingAv<<" (pp-track) +"<<ABDataSystAvHigh<<" -"<<ABDataSystAvLow<<" (ab-data) +-"<<ABTrackingAv<<" (ab-track) +"<<scalingHighAv<<" -"<<scalingLowAv<<" (scal) +"<<beautyUncHigh<<" -"<<beautyUncLow<<" (fd)"<<endl;
        //
        // Fill output histos/graphs
        //
        histoBin = hDmesonAverageRAB->FindBin(ptval);
        hDmesonAverageRAB->SetBinContent(histoBin,average);
        hDmesonAverageRAB->SetBinError(histoBin,averageStat);
        Double_t ept = hDmesonAverageRAB->GetBinWidth(histoBin)/2.;
        ept=0.3;
        gRAB_DmesonAverage_GlobalSystematics->SetPoint(ipt,ptval,average);
        gRAB_DmesonAverage_GlobalSystematics->SetPointError(ipt,ept,ept,totalUncLow,totalUncHigh);
        ept=0.2;
        gRAB_DmesonAverage_FeedDownSystematicsElossHypothesis->SetPoint(ipt,ptval,average);
        gRAB_DmesonAverage_FeedDownSystematicsElossHypothesis->SetPointError(ipt,ept,ept,beautyUncLow,beautyUncHigh);
        //
        ept=0.1;
        gRAB_DmesonAverage_TrackingSystematicsPP->SetPoint(ipt,ptval,average);
        gRAB_DmesonAverage_TrackingSystematicsPP->SetPointError(ipt,ept,ept,ppTrackingAv,ppTrackingAv);
        gRAB_DmesonAverage_TrackingSystematicsAB->SetPoint(ipt,ptval,average);
        gRAB_DmesonAverage_TrackingSystematicsAB->SetPointError(ipt,ept,ept,ABTrackingAv,ABTrackingAv);
        gRAB_DmesonAverage_ScalingSystematicsPP->SetPoint(ipt,ptval,average);
        gRAB_DmesonAverage_ScalingSystematicsPP->SetPointError(ipt,ept,ept,scalingLowAv,scalingHighAv);
        gRAB_DmesonAverage_DataSystematicsPP->SetPoint(ipt,ptval,average);
        gRAB_DmesonAverage_DataSystematicsPP->SetPointError(ipt,ept,ept,ppDataSystAvLow,ppDataSystAvHigh);
        gRAB_DmesonAverage_DataSystematicsAB->SetPoint(ipt,ptval,average);
        gRAB_DmesonAverage_DataSystematicsAB->SetPointError(ipt,ept,ept,ABDataSystAvLow,ABDataSystAvHigh);
        histoBin = -1;
        //
        // Printout
        cout<< " pt min (GeV/c),  pt max (GeV/c), Raa(Daverage), +- (stat), + (syst) , - (syst) "<<endl;
        cout<< ptbinlimits[ipt] <<"  "<< ptbinlimits[ipt+1]<< "  "<< average<< "  "<< averageStat<< "  "<< totalUncHigh<<"  "<<totalUncLow<<endl;
        fprintf(resultFile,"%02.0f   %02.0f   %5.3f   %5.3f   %5.3f   %5.3f\n",ptbinlimits[ipt],ptbinlimits[ipt+1],average,averageStat,totalUncHigh,totalUncLow);
    } // end loop on pt bins
    
    fclose(resultFile);
    
    //
    // Now can start drawing
    //
    TH2F* hempty=new TH2F("hempty"," ; p_{T} (GeV/c} ; Nucl. modif. fact.",100,0.,ptbinlimits[nbins],100,0.,2.);
    hempty->SetStats(0);

    TCanvas *cAvCheck = new TCanvas("cAvCheck","Average Dmeson check");
    hempty->Draw();
    hDmesonAverageRAB->SetLineColor(kBlack);
    hDmesonAverageRAB->SetMarkerStyle(20);
    hDmesonAverageRAB->SetMarkerColor(kBlack);
    hDmesonAverageRAB->Draw("esame");
    for(Int_t j=0; j<3; j++) {
        if(!isDmeson[j]) continue;
        hDmesonRaa[j]->SetLineColor(kBlack);
        hDmesonRaa[j]->SetMarkerColor(2+j);
        hDmesonRaa[j]->SetMarkerStyle(21+j);
        gRAB_GlobalSystematics[j]->SetFillStyle(0);
        gRAB_GlobalSystematics[j]->SetLineColor(2+j);
        gRAB_GlobalSystematics[j]->Draw("2");
        hDmesonRaa[j]->Draw("esame");
    }
    gRAB_DmesonAverage_GlobalSystematics->SetFillStyle(0);
    gRAB_DmesonAverage_GlobalSystematics->SetLineColor(kBlack);
    gRAB_DmesonAverage_GlobalSystematics->Draw("2");
    hDmesonAverageRAB->Draw("esame");
    cAvCheck->Update();
    cAvCheck->SaveAs(Form("%s_result_%s.gif",foutname.Data(),ndate.Data()));

    TCanvas *cAv = new TCanvas("cAv","Average Dmeson");
    hDmesonAverageRAB->Draw("e");
    gRAB_DmesonAverage_FeedDownSystematicsElossHypothesis->SetFillStyle(1001);
    gRAB_DmesonAverage_FeedDownSystematicsElossHypothesis->SetFillColor(kMagenta-7);
    gRAB_DmesonAverage_FeedDownSystematicsElossHypothesis->Draw("2");
    gRAB_DmesonAverage_GlobalSystematics->Draw("2");
    hDmesonAverageRAB->Draw("esame");
    cAv->Update();
    
    //
    // Now can start saving the output
    //
    TFile *fout = new TFile(Form("HFPtSpectrumRaa_%s_%s.root",foutname.Data(),ndate.Data()),"recreate");
    hDmesonAverageRAB->Write();
    gRABNorm->Write();
    gRAB_DmesonAverage_GlobalSystematics->Write();
    gRAB_DmesonAverage_FeedDownSystematicsElossHypothesis->Write();
    gRAB_DmesonAverage_TrackingSystematicsPP->Write();
    gRAB_DmesonAverage_TrackingSystematicsAB->Write();
    gRAB_DmesonAverage_ScalingSystematicsPP->Write();
    gRAB_DmesonAverage_DataSystematicsPP->Write();
    gRAB_DmesonAverage_DataSystematicsAB->Write();
    fout->Write();
    
}
