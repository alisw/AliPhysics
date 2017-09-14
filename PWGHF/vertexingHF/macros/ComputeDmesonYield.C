#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TFile.h"
#include "TNtuple.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TMarker.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "AliHFSystErr.h"
#endif

// input files
TString filnamPPref="~/alice/Charm/PbYield/2011/Final/ppref/D0Kpi_276TeV_FONLLExtrapolationAndExtrapolation_from7TeVAlicedata_combinedFD_141211_010312_160712_rebin_1_2_3_4_5_6_8_12_16_24.root";

TString filnamSpectrumNb="~/alice/Charm/PbYield/2011/Final/010/Dzero/HFPtSpectrum_D0toKpi_010_Nb_160615.root";
TString filnamSpectrumFc="~/alice/Charm/PbYield/2011/Final/010/Dzero/HFPtSpectrum_D0toKpi_010_fc_160615.root";
TString filnamRaaNb="~/alice/Charm/PbYield/2011/Final/010/Dzero/HFPtSpectrumRaa_D0toKpi_010_Nb_160615.root";
TString filnamRaaFc="~/alice/Charm/PbYield/2011/Final/010/Dzero/HFPtSpectrumRaa_D0toKpi_010_fc_160615.root";

//configuration
enum dspec{kDzero,kDplus,kDstar,kDs};
Int_t mesonSpecie=kDzero;
TString centrality="0-10";
const Int_t nPtBins=9;
Double_t binlim[nPtBins+1]={1.,2.,3.,4.,5.,6.,8.,12.,16.,24.};
Int_t method=2;    // 1=fc; 2=Nb --> default=2
Int_t optErrFD=1;  // 0=from histos, not combined; 
                   // 1= from ntuple with Rbc hypo, not combining Nb and fc; 
                   // 2= from ntuple with Rbc hypo, combining Nb and fc
                   // --> default=2; change to 1 only to produce files for D-meson ratios
Float_t centHypoFdOverPr=2.;
Float_t lowHypoFdOverPr=1.;
Float_t highHypoFdOverPr=3.;
Double_t normToCsec=1.; // put here the trigger cross section in ub; if ==1 per-event yield are computed
TString collSyst="Pb-Pb";


// Graphical styles
Bool_t draw[nPtBins]={1,1,0,1,0,1,0,0,1};
Int_t colors[nPtBins]={kGray+2,kMagenta+1,kMagenta,kBlue,kCyan,kGreen+2,kYellow+2,kOrange+1,kRed+1};
Int_t lstyle[nPtBins]={9,10,3,5,7,1,3,6,2};


Bool_t PbPbDataSyst(AliHFSystErr *syst, TH1D* heff, Double_t pt, Double_t &dataSystUp, Double_t &dataSystDown);


void ComputeDmesonYield(){

  TString mesName="Dzero";
  Int_t mesCode=1;
  Double_t brat=0.0393;
  TString mesSymb="D^{0}";
  Float_t ptForExtrap=36.;
  if(mesonSpecie==kDplus){
    mesName="Dplus";
    mesCode=2;
    brat=0.0946;
    mesSymb="D^{+}";
    ptForExtrap=24.;
 }else if(mesonSpecie==kDstar){
    mesName="Dstar";
    mesCode=3;
    brat=0.0393*0.677;
    mesSymb="D^{*+}";
    ptForExtrap=24.;
  }else if(mesonSpecie==kDs){
    mesName="Ds";
    mesCode=4;
    brat=0.0227;
    mesSymb="D_{s}^{+}";
    if(collSyst=="Pb-Pb"){
     centHypoFdOverPr=1.;
     lowHypoFdOverPr=1./3.;
     highHypoFdOverPr=3.;
    }
    ptForExtrap=12.;
  }

  TString centralityNoDash=centrality.Data();
  centralityNoDash.ReplaceAll("-","");
  TString filnamChi;
  TString filnamCnt;
  TString filnamCntS;
  Int_t optCombineFDEloss=1; // 0=enevlope, 1= sum in quadrature
  Bool_t correctForBR=kTRUE;
  Bool_t showRcbSystNorm=kTRUE;

  if(method==1){
    filnamChi=filnamSpectrumFc.Data();
    filnamCnt=filnamRaaFc.Data();
    filnamCntS=filnamRaaNb.Data();
  }
  else if(method==2){
    filnamChi=filnamSpectrumNb.Data();
    filnamCnt=filnamRaaNb.Data();
    filnamCntS=filnamRaaFc.Data();
  }

  Int_t colorSystFD=kGray+1;
  Int_t colorSystRb=kOrange;
  Int_t colorSystDatappC=kBlue+1;//kMagenta+1;
  Int_t colorSystDataaaC=kRed+1;//kMagenta+1;
  Int_t colorppC=kBlue+1;
  Int_t coloraaC=kRed+1;
  Int_t linestppC=1;
  Int_t linestaaC=1;
  Int_t linewippC=2;
  Int_t linewiaaC=2;
  Int_t markerppC=21;
  Int_t markeraaC=23;
  Int_t msizppC=1.2;
  Int_t msizaaC=1.2;
  Float_t sizesystdata=0.4;
  Float_t sizesystfd=0.3;
  Float_t sizesystrb=0.2;
  Float_t sizesysttot=0.15;


  // pp reference

  TFile *filPP=new TFile(filnamPPref.Data());
  TH1D *hppRef = (TH1D*)filPP->Get("hReference");
  if(!hppRef) hppRef = (TH1D*)filPP->Get("fhScaledData");
  TH1D *hppRefSystData = (TH1D*)filPP->Get("fhScaledSystData"); 

  Float_t relsysterrLowDatapp[nPtBins],relsysterrHiDatapp[nPtBins];
  TGraphAsymmErrors* gppRefSyst=(TGraphAsymmErrors*)filPP->Get("gReferenceSyst");
  if(gppRefSyst){
    for(Int_t i=0; i<gppRefSyst->GetN(); i++){
      Double_t x,y;
      gppRefSyst->GetPoint(i,x,y);
      Int_t hBin=TMath::BinarySearch(nPtBins,binlim,x);
      if(x>ptForExtrap){
	relsysterrLowDatapp[hBin]=gppRefSyst->GetErrorYlow(i)/y;
	relsysterrHiDatapp[hBin]=gppRefSyst->GetErrorYhigh(i)/y;
      }
    }
  }
  Float_t relstaterrpp[nPtBins];
  for(Int_t ib=0; ib<nPtBins; ib++){
    Int_t hBin=hppRef->FindBin(0.5*(binlim[ib]+binlim[ib+1]));
    relstaterrpp[ib]=hppRef->GetBinError(hBin)/hppRef->GetBinContent(hBin);
    printf("-- Relative stat errors PP: Bin %d(%f)--- histobin=%d Err %f\n",ib,0.5*(binlim[ib]+binlim[ib+1]),hBin,relstaterrpp[ib]);
    Int_t hBin2=hppRefSystData->FindBin(0.5*(binlim[ib]+binlim[ib+1]));
    if(hppRefSystData->GetBinCenter(hBin2)<ptForExtrap){
      Float_t relsysterrDatapp=hppRefSystData->GetBinError(hBin2)/hppRefSystData->GetBinContent(hBin2);
      relsysterrLowDatapp[ib]=relsysterrDatapp;
      relsysterrHiDatapp[ib]=relsysterrDatapp;
    }
    printf(" ---- Check SYST err DATA PP Bin %d  CS=%f Err+%f -%f\n",ib, hppRef->GetBinContent(hBin),relsysterrHiDatapp[ib],relsysterrLowDatapp[ib]);
  }

  Float_t relsysterrLowEnScalpp[nPtBins];
  Float_t relsysterrHiEnScalpp[nPtBins];  
  for(Int_t ib=0; ib<nPtBins; ib++){
    relsysterrLowEnScalpp[ib]=0.;
    relsysterrHiEnScalpp[ib]=0.;
  }
  TGraphAsymmErrors* gsystppEnSc=(TGraphAsymmErrors*)filPP->Get("gScaledDataSystExtrap");
  for(Int_t i=0; i<gsystppEnSc->GetN(); i++){
    Double_t x,y;
    gsystppEnSc->GetPoint(i,x,y);
    Int_t hBin=TMath::BinarySearch(nPtBins,binlim,x);
    if(hBin>=0 && hBin<nPtBins && x<binlim[nPtBins]){
      if(x<ptForExtrap){
	relsysterrLowEnScalpp[hBin]=gsystppEnSc->GetErrorYlow(i)/y;
	relsysterrHiEnScalpp[hBin]=gsystppEnSc->GetErrorYhigh(i)/y;
      }
    }
  }

  TGraphAsymmErrors* gsystppFD=(TGraphAsymmErrors*)filPP->Get("gScaledDataSystFeedDown");
  Float_t relsysterrLowFDpp[nPtBins];
  Float_t relsysterrHiFDpp[nPtBins];  
  for(Int_t i=0; i<nPtBins; i++){ 
    relsysterrHiFDpp[i]=0.;
    relsysterrLowFDpp[i]=0.;
  }
  for(Int_t i=0; i<gsystppFD->GetN(); i++){
    Double_t x,y;
    gsystppFD->GetPoint(i,x,y);
    Int_t hBin=TMath::BinarySearch(nPtBins,binlim,x);
    if(hBin>=0 && hBin<nPtBins && x<binlim[nPtBins]){
      relsysterrLowFDpp[hBin]=gsystppFD->GetErrorYlow(i)/y;
      relsysterrHiFDpp[hBin]=gsystppFD->GetErrorYhigh(i)/y;
      printf(" ---- Check syst err FD pp Bin %d  Err+%f -%f\n",hBin,relsysterrHiFDpp[hBin],relsysterrLowFDpp[hBin]);
      //      printf("%d %f %f\n",hBin,relsysterrLowFDpp[hBin],relsysterrHiFDpp[hBin]);
    }
  }

  // A-A events
							  
  printf("--- %s events ---\n",collSyst.Data());

  TFile *filChi= new TFile(filnamChi.Data());
  TH1D *hSpC = (TH1D*)filChi->Get("hRECpt");
  Float_t relstaterrPbPb[nPtBins];
  for(Int_t ib=0; ib<nPtBins; ib++){
    Int_t hBin=hSpC->FindBin(0.5*(binlim[ib]+binlim[ib+1]));
    relstaterrPbPb[ib]=hSpC->GetBinError(hBin)/hSpC->GetBinContent(hBin);
    printf("-- Relative stat errors AA from yield: Bin %d(%d) pt=%f Err %f\n",ib,hBin,hSpC->GetBinCenter(hBin),relstaterrPbPb[ib]);
  }

  AliHFSystErr *systematicsABcent =(AliHFSystErr*)filChi->Get("AliHFSystErr");
  if(!systematicsABcent){
    printf("AliHFSysstErr generated on the fly\n");
    systematicsABcent=new AliHFSystErr("AliHFSystErr","on the fly");
    if(collSyst=="p-Pb"){
      systematicsABcent->SetCollisionType(2);
      systematicsABcent->SetRunNumber(16);
    }else{
      systematicsABcent->SetCollisionType(1);
    }
    systematicsABcent->SetCentrality(centralityNoDash.Data());
    systematicsABcent->Init(mesCode);
  }else{
    printf("AliHFSystErr read from HFPtSpectrum output\n");
  }
  TH1D* heffC=(TH1D*)filChi->Get("hDirectEffpt");


  TGraphAsymmErrors* gsystaaFDc=(TGraphAsymmErrors*)filChi->Get("gSigmaCorrConservative"); // FD due to scales
  Float_t relsysterrLowFDaa[nPtBins];
  Float_t relsysterrHiFDaa[nPtBins];  
  for(Int_t i=0; i<nPtBins; i++){ 
    relsysterrLowFDaa[i]=0.;
    relsysterrHiFDaa[i]=0.;
  }
  for(Int_t i=0; i<gsystaaFDc->GetN(); i++){
    Double_t x,y;
    gsystaaFDc->GetPoint(i,x,y);
    Int_t hBin=TMath::BinarySearch(nPtBins,binlim,x);
    if(hBin>=0 && hBin<nPtBins && x<binlim[nPtBins]){
      relsysterrLowFDaa[hBin]=gsystaaFDc->GetErrorYlow(i)/y;
      relsysterrHiFDaa[hBin]=gsystaaFDc->GetErrorYhigh(i)/y;
      printf(" ---- Check syst err FD AA Bin %d  Err+%f -%f\n",hBin,relsysterrHiFDaa[hBin],relsysterrLowFDaa[hBin]);
    }
  }

  TFile *filCnt=new TFile(filnamCnt.Data());
  Float_t relstaterrPbPb2[nPtBins];
  TH1D* hraaCcheck2=(TH1D*)filCnt->Get("hRABvsPt");
  for(Int_t ib=0; ib<nPtBins; ib++){
    Int_t hBin=hraaCcheck2->FindBin(0.5*(binlim[ib]+binlim[ib+1]));
    Double_t aux=hraaCcheck2->GetBinError(hBin)/hraaCcheck2->GetBinContent(hBin);
    relstaterrPbPb2[ib]=TMath::Sqrt(aux*aux-relstaterrpp[ib]*relstaterrpp[ib]);
    
    printf("-- Relative stat errors AA from RAA-PP: Bin %d(%f)---%d Err %f\n",ib,0.5*(binlim[ib]+binlim[ib+1]),hBin,relstaterrPbPb2[ib]);
  }
  AliHFSystErr *systematicsPP =(AliHFSystErr*)filCnt->Get("AliHFSystErrPP");

  TNtuple* ntC=(TNtuple*)filCnt->Get("ntupleRAB");
  Float_t pt,TAB,sigmaPP,invyieldAB,RABCharm,RABBeauty;
  Float_t invyieldABFDHigh,invyieldABFDLow,fprompt;
  ntC->SetBranchAddress("pt",&pt);
  ntC->SetBranchAddress("TAB",&TAB);
  ntC->SetBranchAddress("sigmaPP",&sigmaPP);
  ntC->SetBranchAddress("invyieldAB",&invyieldAB);
  ntC->SetBranchAddress("RABCharm",&RABCharm);
  ntC->SetBranchAddress("RABBeauty",&RABBeauty);
  ntC->SetBranchAddress("invyieldABFDHigh",&invyieldABFDHigh);
  ntC->SetBranchAddress("invyieldABFDLow",&invyieldABFDLow);
  ntC->SetBranchAddress("fc",&fprompt);

  Float_t raac[nPtBins],invypp[nPtBins],invyPbPb[nPtBins],invyPbPbLo[nPtBins],invyPbPbHi[nPtBins],minval[nPtBins];
  Float_t raacLowHyp[nPtBins],invyPbPbLowHyp[nPtBins],minvalLowHyp[nPtBins];
  Float_t raacHigHyp[nPtBins],invyPbPbHigHyp[nPtBins],minvalHigHyp[nPtBins];
  Float_t fPromptCent[nPtBins],fPromptHigHyp[nPtBins],fPromptLowHyp[nPtBins];
  Float_t invyPbPbLoSingleSyst[nPtBins],invyPbPbHiSingleSyst[nPtBins];

  for(Int_t ib=0; ib<nPtBins; ib++){
    minval[ib]=9999.;
    invypp[ib]=-999.;
    invyPbPb[ib]=-999.;
    invyPbPbLo[ib]=-999.;
    invyPbPbHi[ib]=-999.;
    invyPbPbLoSingleSyst[ib]=-999.;
    invyPbPbHiSingleSyst[ib]=-999.;
    raac[ib]=-999.;
    raacLowHyp[ib]=-999.;
    invyPbPbLowHyp[ib]=-999.;
    minvalLowHyp[ib]=9999.;
    raacHigHyp[ib]=-999.;
    invyPbPbHigHyp[ib]=-999.;
    minvalHigHyp[ib]=9999.;
    fPromptCent[ib]=-9999.;
    fPromptHigHyp[ib]=-9999.;
    fPromptLowHyp[ib]=-9999.;
  }
  
  TGraph** gcb=new TGraph*[nPtBins];
  TGraph** gcrbc=new TGraph*[nPtBins];
  for(Int_t ib=0; ib<nPtBins; ib++){
    gcb[ib]=new TGraph(0);
    gcrbc[ib]=new TGraph(0);
  }
  
  for(Int_t ie=0; ie<ntC->GetEntries(); ie++){
    ntC->GetEvent(ie);
    if(correctForBR){
      invyieldAB/=brat;
      invyieldABFDLow/=brat;
      invyieldABFDHigh/=brat;
      sigmaPP/=brat;
    }
    if(pt>binlim[nPtBins]) continue;
    Int_t theBin=TMath::BinarySearch(nPtBins,binlim,(Double_t)pt);
    if(theBin<0 || theBin>=nPtBins) continue;
    Float_t rFdPr=RABBeauty/RABCharm;
    Float_t dist=TMath::Abs(rFdPr-centHypoFdOverPr);
    if(dist<minval[theBin]){
      raac[theBin]=RABCharm;
      minval[theBin]=dist;
      invyPbPb[theBin]=invyieldAB*normToCsec;
      invyPbPbLo[theBin]=invyieldABFDLow*normToCsec;
      invyPbPbHi[theBin]=invyieldABFDHigh*normToCsec;
      invypp[theBin]=TAB*sigmaPP;
      if(collSyst=="p-Pb" && TMath::Abs(normToCsec-1.)>0.001) invypp[theBin]=208.*sigmaPP*1e6;
      fPromptCent[theBin]=fprompt;
    }
    Float_t distLowHyp=TMath::Abs(rFdPr-lowHypoFdOverPr);
    if(distLowHyp<minvalLowHyp[theBin]){
      // LowHyp -> lower Raa(feeddown) -> less Fd to be subtracted -> higher fprompt  -> higher prompt yield -> higher Raa(prompt)
      raacLowHyp[theBin]=RABCharm;
      invyPbPbLowHyp[theBin]=invyieldAB*normToCsec;
      invyPbPbHiSingleSyst[theBin]=invyieldABFDHigh*normToCsec;
      minvalLowHyp[theBin]=distLowHyp;
      fPromptLowHyp[theBin]=fprompt;
    }
    Float_t distHigHyp=TMath::Abs(rFdPr-highHypoFdOverPr);
    if(distHigHyp<minvalHigHyp[theBin]){
      // HigHyp -> higher Raa(feeddown) -> more Fd to be subtracted -> lower fprompt  -> lower prompt yield -> lower Raa(prompt)
      raacHigHyp[theBin]=RABCharm;
      invyPbPbHigHyp[theBin]=invyieldAB*normToCsec;
      invyPbPbLoSingleSyst[theBin]=invyieldABFDLow*normToCsec;
      minvalHigHyp[theBin]=distHigHyp;
      fPromptHigHyp[theBin]=fprompt;
    }
    if(theBin<nPtBins && theBin>=0 && pt<binlim[nPtBins]){
      gcb[theBin]->SetPoint(gcb[theBin]->GetN(),RABBeauty,RABCharm);
      gcrbc[theBin]->SetPoint(gcrbc[theBin]->GetN(),rFdPr,RABCharm);
    }
  }

  TH1D* hcheck[nPtBins];
  for(Int_t i=0; i<nPtBins; i++){
    hcheck[i]=(TH1D*)filCnt->Get(Form("hRCharmVsElossHypo_%d",i+1));
  }

  TMarker** mC=new TMarker*[nPtBins];
  for(Int_t ib=0; ib<nPtBins; ib++){
    mC[ib]=new TMarker(1.,raac[ib],20);
    mC[ib]->SetMarkerSize(1.2);
    if(showRcbSystNorm){
      Double_t auxx,auxy;
      for(Int_t ip=0; ip<gcrbc[ib]->GetN(); ip++){
	gcrbc[ib]->GetPoint(ip,auxx,auxy);
	auxy/=raac[ib];
	gcrbc[ib]->SetPoint(ip,auxx,(auxy-1)*100.);
      }
    }
  }

  TFile *filCntS=new TFile(filnamCntS.Data());
  TNtuple* ntCS=(TNtuple*)filCntS->Get("ntupleRAB");
  ntCS->SetBranchAddress("pt",&pt);
  ntCS->SetBranchAddress("TAB",&TAB);
  ntCS->SetBranchAddress("sigmaPP",&sigmaPP);
  ntCS->SetBranchAddress("invyieldAB",&invyieldAB);
  ntCS->SetBranchAddress("RABCharm",&RABCharm);
  ntCS->SetBranchAddress("RABBeauty",&RABBeauty);
  ntCS->SetBranchAddress("invyieldABFDHigh",&invyieldABFDHigh);
  ntCS->SetBranchAddress("invyieldABFDLow",&invyieldABFDLow);
  ntCS->SetBranchAddress("fc",&fprompt);

  Float_t invyPbPbS[nPtBins],invyPbPbLoS[nPtBins],invyPbPbHiS[nPtBins],minvalS[nPtBins];
  Float_t invyPbPbLoSingleSystS[nPtBins],invyPbPbHiSingleSystS[nPtBins],minvalLowHypS[nPtBins],minvalHigHypS[nPtBins];
  Float_t fPromptCentS[nPtBins],fPromptLowHypS[nPtBins],fPromptHigHypS[nPtBins];
  for(Int_t ib=0; ib<nPtBins; ib++){
    minvalS[ib]=9999.;
    minvalLowHypS[ib]=9999.;
    minvalHigHypS[ib]=9999.;
    invyPbPbS[ib]=-999.;
    invyPbPbLoS[ib]=-999.;
    invyPbPbHiS[ib]=-999.;
    invyPbPbLoSingleSystS[ib]=-999.;
    invyPbPbHiSingleSystS[ib]=-999.;
    fPromptLowHypS[ib]=-999.;
    fPromptHigHypS[ib]=-999.;
  }
  for(Int_t ie=0; ie<ntCS->GetEntries(); ie++){
    ntCS->GetEvent(ie);
    if(correctForBR){
      invyieldAB/=brat;
      invyieldABFDLow/=brat;
      invyieldABFDHigh/=brat;
      sigmaPP/=brat;
    }
    if(pt>binlim[nPtBins]) continue;
    Int_t theBin=TMath::BinarySearch(nPtBins,binlim,(Double_t)pt);
    if(theBin<0 || theBin>=nPtBins) continue;
    Float_t rFdPr=RABBeauty/RABCharm;
    Float_t dist=TMath::Abs(rFdPr-centHypoFdOverPr);
    if(dist<minvalS[theBin]){
      minvalS[theBin]=dist;
      invyPbPbS[theBin]=invyieldAB*normToCsec;
      invyPbPbLoS[theBin]=invyieldABFDLow*normToCsec;
      invyPbPbHiS[theBin]=invyieldABFDHigh*normToCsec;
      fPromptCentS[theBin]=fprompt;
    }
    Float_t distLowHyp=TMath::Abs(rFdPr-lowHypoFdOverPr);
    if(distLowHyp<minvalLowHypS[theBin]){
      // LowHyp -> lower Raa(feeddown) -> less Fd to be subtracted -> higher fprompt  -> higher prompt yield -> higher Raa(prompt)
      invyPbPbHiSingleSystS[theBin]=invyieldABFDHigh*normToCsec;
      minvalLowHypS[theBin]=distLowHyp;
      fPromptLowHypS[theBin]=fprompt;
    }
    Float_t distHigHyp=TMath::Abs(rFdPr-highHypoFdOverPr);
    if(distHigHyp<minvalHigHypS[theBin]){
      // HigHyp -> higher Raa(feeddown) -> more Fd to be subtracted -> lower fprompt  -> lower prompt yield -> lower Raa(prompt)
      invyPbPbLoSingleSystS[theBin]=invyieldABFDLow*normToCsec;
      minvalHigHypS[theBin]=distHigHyp;
      fPromptHigHypS[theBin]=fprompt;
    }
  }

  TH1F* hfPromptCent=new TH1F("hfPromptCent"," ; p_{T} (Gev/c) ; f_{prompt}",nPtBins,binlim);
  TH1F* hfPromptMinNb=new TH1F("hfPromptMinNb"," ; p_{T} (Gev/c) ; f_{prompt}",nPtBins,binlim);
  TH1F* hfPromptMinfc=new TH1F("hfPromptMinfc"," ; p_{T} (Gev/c) ; f_{prompt}",nPtBins,binlim);
  TH1F* hfPromptMaxNb=new TH1F("hfPromptMaxNb"," ; p_{T} (Gev/c) ; f_{prompt}",nPtBins,binlim);
  TH1F* hfPromptMaxfc=new TH1F("hfPromptMaxfc"," ; p_{T} (Gev/c) ; f_{prompt}",nPtBins,binlim);
  hfPromptCent->SetStats(0);
  hfPromptMinNb->SetStats(0);
  hfPromptMaxNb->SetStats(0);
  hfPromptMinfc->SetStats(0);
  hfPromptMaxfc->SetStats(0);

  for(Int_t ib=0; ib<nPtBins; ib++){
    printf("Bin %d\n",ib);
    printf("   fprompt central=%f ---  Nb fpromptmin=%f fpromptmax=%f ---  fc fpromptmin=%f fpromptmax=%f\n",fPromptCent[ib],fPromptHigHyp[ib],fPromptLowHyp[ib],fPromptHigHypS[ib],fPromptLowHypS[ib]);
    
    //add error from FONLL scale in quadrature
    Double_t relerrFDscaleHigh = (invyPbPbHi[ib]-invyPbPb[ib])/invyPbPb[ib];
    Double_t relerrFDscaleLow = (invyPbPb[ib]-invyPbPbLo[ib])/invyPbPb[ib];

    Double_t relerrElossHigh = (fPromptLowHyp[ib]-fPromptCent[ib])/fPromptCent[ib];
    Double_t relerrElossLow = (fPromptCent[ib]-fPromptHigHyp[ib])/fPromptCent[ib];
    
    Double_t toterrFDhigh = TMath::Sqrt(relerrFDscaleHigh*relerrFDscaleHigh+relerrElossHigh*relerrElossHigh)*fPromptCent[ib];
    Double_t toterrFDlow = TMath::Sqrt(relerrFDscaleLow*relerrFDscaleLow+relerrElossLow*relerrElossLow)*fPromptCent[ib];

    Double_t relerrFDscaleHighS = (invyPbPbHiS[ib]-invyPbPbS[ib])/invyPbPbS[ib];
    Double_t relerrFDscaleLowS = (invyPbPbS[ib]-invyPbPbLoS[ib])/invyPbPbS[ib];
    
    Double_t relerrElossHighS = (fPromptLowHypS[ib]-fPromptCentS[ib])/fPromptCentS[ib];
    Double_t relerrElossLowS = (fPromptCentS[ib]-fPromptHigHypS[ib])/fPromptCentS[ib];
    
    Double_t toterrFDhighS = TMath::Sqrt(relerrFDscaleHighS*relerrFDscaleHighS+relerrElossHighS*relerrElossHighS)*fPromptCentS[ib];
    Double_t toterrFDlowS = TMath::Sqrt(relerrFDscaleLowS*relerrFDscaleLowS+relerrElossLowS*relerrElossLowS)*fPromptCentS[ib];

    hfPromptCent->SetBinContent(ib+1,fPromptCent[ib]);
    hfPromptMinNb->SetBinContent(ib+1,fPromptCent[ib]-toterrFDlow);
    hfPromptMaxNb->SetBinContent(ib+1,fPromptCent[ib]+toterrFDhigh);
    hfPromptMinfc->SetBinContent(ib+1,fPromptCentS[ib]-toterrFDlowS);
    hfPromptMaxfc->SetBinContent(ib+1,fPromptCentS[ib]+toterrFDhighS);
    
    printf("Bin %d\n",ib);
    printf("   fprompt central=%f ---  Nb fpromptmin=%f fpromptmax=%f ---  fc fpromptmin=%f fpromptmax=%f\n",fPromptCent[ib],hfPromptMinNb->GetBinContent(ib+1),hfPromptMaxNb->GetBinContent(ib+1),hfPromptMinfc->GetBinContent(ib+1),hfPromptMaxfc->GetBinContent(ib+1));
  }

  TH1F* hppC=new TH1F("hppC",Form("pp reference for %s%% CC",centrality.Data()),nPtBins,binlim);
  TGraphAsymmErrors *gppCsystdata=new TGraphAsymmErrors(0);
  gppCsystdata->SetName("gppCsystdata");  
  gppCsystdata->SetTitle(Form("Data Syst. Err. pp, scaled to %s%% CC",centrality.Data()));  
  TGraphAsymmErrors *gppCsystFD=new TGraphAsymmErrors(0);
  gppCsystFD->SetName("gppCsystFD");
  gppCsystFD->SetTitle(Form("B feed-down Syst. Err. pp, scaled to  %s%% CC",centrality.Data()));
  TH1F* hAAC=new TH1F("hAAC",Form("PbPb %s%% CC",centrality.Data()),nPtBins,binlim);
  TGraphAsymmErrors *gaaCsystdata=new TGraphAsymmErrors(0);
  gaaCsystdata->SetName("gaaCsystdata");  
  gaaCsystdata->SetTitle(Form("Data Syst. Err. PbPb  %s%% CC",centrality.Data()));
  TGraphAsymmErrors *gaaCsystFD=new TGraphAsymmErrors(0);
  gaaCsystFD->SetName("gaaCsystFD");  
  gaaCsystFD->SetTitle(Form("B feed-down Syst. Err. PbPb %s%% CC",centrality.Data()));
  TGraphAsymmErrors *gaaCsystRb=new TGraphAsymmErrors(0);
  gaaCsystRb->SetName("gaaCsystRb");  
  gaaCsystRb->SetTitle(Form("Raa(B) Syst. Err. PbPb %s%% CC",centrality.Data()));
  TGraphAsymmErrors *gaaCsystB=new TGraphAsymmErrors(0);
  gaaCsystB->SetName("gaaCsystB");  
  gaaCsystB->SetTitle(Form("B Syst. Err. PbPb %s%% CC",centrality.Data()));
  TGraphAsymmErrors *gaaCsystTot=new TGraphAsymmErrors(0);
  gaaCsystTot->SetName("gaaCsystTot");  
  gaaCsystTot->SetTitle(Form("Tot Syst. Err. PbPb %s%% CC",centrality.Data()));
  TH1F* hRAAC=new TH1F("hRAAC","",nPtBins,binlim);
  TGraphAsymmErrors *graaC=new TGraphAsymmErrors(0);

  for(Int_t ib=0; ib<nPtBins; ib++){
    printf("Bin %d\n",ib);
    printf("raa(centHyp)=%f   raa(lowHyp)=%f   raa(highHyp)=%f  RelDiff=%f\n",
	     raac[ib],raacLowHyp[ib],raacHigHyp[ib],(raacHigHyp[ib]-raacLowHyp[ib])/raac[ib]);
    if(raac[ib]>0.){
      Float_t relstaterrRaa=TMath::Sqrt(relstaterrpp[ib]*relstaterrpp[ib]+relstaterrPbPb2[ib]*relstaterrPbPb2[ib]);
      hRAAC->SetBinContent(ib+1,raac[ib]);
      hRAAC->SetBinError(ib+1,relstaterrRaa*raac[ib]);
      Int_t nP=graaC->GetN();
      graaC->SetPoint(nP,hRAAC->GetBinCenter(ib+1),raac[ib]);
      graaC->SetPointEXlow(nP,hRAAC->GetBinCenter(ib+1)-hRAAC->GetBinLowEdge(ib+1));
      graaC->SetPointEXhigh(nP,hRAAC->GetBinLowEdge(ib+2)-hRAAC->GetBinCenter(ib+1));
      graaC->SetPointEYlow(nP,raac[ib]-raacHigHyp[ib]);
      graaC->SetPointEYhigh(nP,raacLowHyp[ib]-raac[ib]);      
    }
    if(invypp[ib]>0.){
      hppC->SetBinContent(ib+1,invypp[ib]);
      hppC->SetBinError(ib+1,relstaterrpp[ib]*invypp[ib]);
      Int_t nP=gppCsystdata->GetN();
      gppCsystdata->SetPoint(nP,hppC->GetBinCenter(ib+1),invypp[ib]);
      gppCsystdata->SetPointEXlow(nP,sizesystdata);
      gppCsystdata->SetPointEXhigh(nP,sizesystdata);
      Double_t edatl=relsysterrLowDatapp[ib]*invypp[ib];
      Double_t edath=relsysterrHiDatapp[ib]*invypp[ib];
      Double_t eenscl=relsysterrLowEnScalpp[ib]*invypp[ib];
      Double_t eensch=relsysterrHiEnScalpp[ib]*invypp[ib];
      Double_t elow=TMath::Sqrt(edatl*edatl+eenscl*eenscl);
      Double_t ehig=TMath::Sqrt(edath*edath+eensch*eensch);
      gppCsystdata->SetPointEYlow(nP,elow);
      gppCsystdata->SetPointEYhigh(nP,ehig);
      nP=gppCsystFD->GetN();
      gppCsystFD->SetPoint(nP,hppC->GetBinCenter(ib+1),invypp[ib]);
      gppCsystFD->SetPointEXlow(nP,sizesystfd);
      gppCsystFD->SetPointEXhigh(nP,sizesystfd);
      gppCsystFD->SetPointEYlow(nP,relsysterrLowFDpp[ib]*invypp[ib]);
      gppCsystFD->SetPointEYhigh(nP,relsysterrHiFDpp[ib]*invypp[ib]);
    }
    if(invyPbPb[ib]>0.){
      printf("yielddAA(centHyp)=%f   yielddAA(lowHyp)=%f   yielddAA(highHyp)=%f  RelDiff=%f\n",
	     invyPbPb[ib],invyPbPbLowHyp[ib],invyPbPbHigHyp[ib],(invyPbPbHigHyp[ib]-invyPbPbLowHyp[ib])/invyPbPb[ib]);
      hAAC->SetBinContent(ib+1,invyPbPb[ib]);
      hAAC->SetBinError(ib+1,relstaterrPbPb2[ib]*invyPbPb[ib]);
      Int_t nP=gaaCsystdata->GetN();
      gaaCsystdata->SetPoint(nP,hAAC->GetBinCenter(ib+1),invyPbPb[ib]);
      gaaCsystdata->SetPointEXlow(nP,sizesystdata);
      gaaCsystdata->SetPointEXhigh(nP,sizesystdata);
      Double_t systundatalow,systundatahigh;
      Bool_t isOk=PbPbDataSyst(systematicsABcent,heffC,hAAC->GetBinCenter(ib+1),systundatahigh,systundatalow);
      printf("Check PID syst: %f %f %f\n",systematicsABcent->GetTotalSystErr(hAAC->GetBinCenter(ib+1)),systundatahigh,systundatalow);
      systundatalow*=invyPbPb[ib];
      systundatahigh*=invyPbPb[ib];
      gaaCsystdata->SetPointEYlow(nP,systundatalow);
      gaaCsystdata->SetPointEYhigh(nP,systundatahigh);

      nP=gaaCsystRb->GetN();
      gaaCsystRb->SetPoint(nP,hAAC->GetBinCenter(ib+1),invyPbPb[ib]);
      gaaCsystRb->SetPointEXlow(nP,sizesystrb);
      gaaCsystRb->SetPointEXhigh(nP,sizesystrb);
      gaaCsystRb->SetPointEYlow(nP,invyPbPb[ib]-invyPbPbHigHyp[ib]);
      gaaCsystRb->SetPointEYhigh(nP,invyPbPbLowHyp[ib]-invyPbPb[ib]);      
      nP=gaaCsystFD->GetN();
      gaaCsystFD->SetPoint(nP,hAAC->GetBinCenter(ib+1),invyPbPb[ib]);
      gaaCsystFD->SetPointEXlow(nP,sizesystfd);
      gaaCsystFD->SetPointEXhigh(nP,sizesystfd);
      gaaCsystB->SetPoint(nP,hAAC->GetBinCenter(ib+1),invyPbPb[ib]);
      gaaCsystB->SetPointEXlow(nP,sizesystfd);
      gaaCsystB->SetPointEXhigh(nP,sizesystfd);
      gaaCsystTot->SetPoint(nP,hAAC->GetBinCenter(ib+1),invyPbPb[ib]);
      gaaCsystTot->SetPointEXlow(nP,sizesysttot);
      gaaCsystTot->SetPointEXhigh(nP,sizesysttot);
      if(optErrFD==0){
	gaaCsystFD->SetPointEYlow(nP,relsysterrLowFDaa[ib]*invyPbPb[ib]);
	gaaCsystFD->SetPointEYhigh(nP,relsysterrHiFDaa[ib]*invyPbPb[ib]);
      }else{	
	Double_t minyield,maxyield;
	Double_t minyieldsing,maxyieldsing;
	if(optErrFD==1){
	  minyield=invyPbPbLo[ib];
	  maxyield=invyPbPbHi[ib];
	  minyieldsing=invyPbPbLoSingleSyst[ib];
	  maxyieldsing=invyPbPbHiSingleSyst[ib];
	}else{
	  minyield=TMath::Min(invyPbPbLo[ib],invyPbPbLoS[ib]);
	  maxyield=TMath::Max(invyPbPbHi[ib],invyPbPbHiS[ib]);
	  minyieldsing=TMath::Min(invyPbPbLoSingleSyst[ib],invyPbPbLoSingleSystS[ib]);
	  maxyieldsing=TMath::Max(invyPbPbHiSingleSyst[ib],invyPbPbHiSingleSystS[ib]);	  
	}
	gaaCsystFD->SetPointEYlow(nP,invyPbPb[ib]-minyield);
	gaaCsystFD->SetPointEYhigh(nP,maxyield-invyPbPb[ib]);
	printf("Relative syst from FD (FONLL scales+masses)=-%f +%f\n",gaaCsystFD->GetErrorYlow(nP)/invyPbPb[ib],gaaCsystFD->GetErrorYhigh(nP)/invyPbPb[ib]);
	Double_t systunBlow=0.;
	Double_t systunBhigh=0.;
	if(optCombineFDEloss==0){ // envelope
	  systunBlow=invyPbPb[ib]-minyieldsing;
	  systunBhigh=maxyieldsing-invyPbPb[ib];
	}else{ // sum in quadrature
	  Double_t e1l=gaaCsystFD->GetErrorYlow(nP);
	  Double_t e1h=gaaCsystFD->GetErrorYhigh(nP);
	  Double_t e2l=gaaCsystRb->GetErrorYlow(nP);
	  Double_t e2h=gaaCsystRb->GetErrorYhigh(nP);
	  systunBlow=TMath::Sqrt(e1l*e1l+e2l*e2l);
	  systunBhigh=TMath::Sqrt(e1h*e1h+e2h*e2h);
	}
	gaaCsystB->SetPointEYlow(nP,systunBlow);
	gaaCsystB->SetPointEYhigh(nP,systunBhigh);
	
	Double_t totSystlow=TMath::Sqrt(systunBlow*systunBlow+systundatalow*systundatalow);
	Double_t totSysthigh=TMath::Sqrt(systunBhigh*systunBhigh+systundatahigh*systundatahigh);
	gaaCsystTot->SetPointEYlow(nP,totSystlow);
	gaaCsystTot->SetPointEYhigh(nP,totSysthigh);
      }
    }
  }

  TH1F* hcheckRAAC=(TH1F*)hAAC->Clone("hcheckRAAC");
  hcheckRAAC->Divide(hppC);


  // Plots

  gStyle->SetOptTitle(0);

  TCanvas* c1=new TCanvas("c1","RcVsRb");
  TLegend* leg= new TLegend(0.19,0.18,0.55,0.41);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(43);
  leg->SetTextSize(28);
  leg->SetMargin(0.32);
  TLegendEntry* ent;
  Bool_t first=kTRUE;

  for(Int_t ib=0; ib<nPtBins; ib++){
    if(draw[ib] && gcb[ib]->GetN()>0){
      gcb[ib]->SetLineColor(colors[ib]);
      gcb[ib]->SetLineWidth(3);
      gcb[ib]->SetLineStyle(lstyle[ib]);
      if(first){
	gcb[ib]->Draw("AL");
	gcb[ib]->GetXaxis()->SetLimits(0.,1.);
	gcb[ib]->GetXaxis()->SetTitle("#it{R}_{AA} feed-down");
	gcb[ib]->GetYaxis()->SetTitle("#it{R}_{AA} prompt");
	gcb[ib]->SetMinimum(0.);
	gcb[ib]->SetMaximum(1.);
	first=kFALSE;
      }else{
	gcb[ib]->Draw("lsame");
      }
      ent=leg->AddEntry(gcb[ib],Form("%d<#it{p}_{T}<%d GeV/#it{c}",(Int_t)binlim[ib],(Int_t)binlim[ib+1]),"L");
    }
  }
  leg->Draw();

  first=kTRUE;

  TCanvas* c2=new TCanvas("c2x","RcVsRcb",700,700);
  c2->SetTickx();
  c2->SetTicky();
  // c2->SetGridx();
  // c2->SetGridy();
  c2->SetLeftMargin(0.14);
  c2->SetBottomMargin(0.13);
  c2->SetTopMargin(0.035);
  c2->SetRightMargin(0.045);
  for(Int_t ib=0; ib<nPtBins; ib++){
    if(draw[ib] && gcrbc[ib]->GetN()>0){
      gcrbc[ib]->SetLineColor(colors[ib]);
      gcrbc[ib]->SetLineWidth(3);
      gcrbc[ib]->SetLineStyle(lstyle[ib]);
      if(first){
	gcrbc[ib]->Draw("AL");
	gcrbc[ib]->GetXaxis()->SetLimits(lowHypoFdOverPr,highHypoFdOverPr);
	gcrbc[ib]->GetXaxis()->SetTitle("Hypothesis on (#it{R}_{AA} feed-down)/(#it{R}_{AA} prompt)");
	gcrbc[ib]->GetYaxis()->SetTitleOffset(1.2);
	gcrbc[ib]->GetXaxis()->SetTitleOffset(1.2);
	gcrbc[ib]->GetYaxis()->SetTitleFont(43);
	gcrbc[ib]->GetXaxis()->SetTitleFont(43);
	gcrbc[ib]->GetYaxis()->SetTitleSize(30);
	gcrbc[ib]->GetXaxis()->SetTitleSize(30);
	gcrbc[ib]->GetYaxis()->SetLabelFont(43);
	gcrbc[ib]->GetXaxis()->SetLabelFont(43);
	gcrbc[ib]->GetYaxis()->SetLabelSize(28);
	gcrbc[ib]->GetXaxis()->SetLabelSize(28);
	if(!showRcbSystNorm){
	  gcrbc[ib]->GetYaxis()->SetTitle("#it{R}_{AA} prompt");
	  gcrbc[ib]->SetMinimum(0.);
	  gcrbc[ib]->SetMaximum(0.8);
	}else{
	  gcrbc[ib]->GetYaxis()->SetTitle("Relative variation of #it{R}_{AA} prompt (%)");
	  gcrbc[ib]->SetMinimum(-20);
	  gcrbc[ib]->SetMaximum(20);
	}
	first=kFALSE;	
      }else{
	gcrbc[ib]->Draw("lsame");
      }
      if(!showRcbSystNorm){
	mC[ib]->SetMarkerColor(colors[ib]);
	mC[ib]->Draw("same");
      }
    }
  }
//   TBox* b=new TBox(0.5,0.,2.,1.);
//   b->SetFillColor(kYellow);      
//   b->SetFillStyle(3003);
  //  b->Draw("same");    
  leg->Draw();
  TLatex* tali0=new TLatex(0.18,0.89,"ALICE");
  tali0->SetNDC();
  tali0->SetTextFont(43);
  tali0->SetTextSize(28);
  tali0->Draw();
  TLatex* tpbpb=new TLatex(0.41,0.89,Form("%s, %s%%",collSyst.Data(),centrality.Data()));
  tpbpb->SetNDC();
  tpbpb->SetTextFont(43);
  tpbpb->SetTextSize(28);
  tpbpb->Draw();
  // TLatex* tdmes=new TLatex(0.63,0.85,Form("%s meson",mesSymb.Data()));
  TLatex* tdmes=new TLatex(0.3,0.43,Form("%s meson",mesSymb.Data()));
  tdmes->SetNDC();
  tdmes->SetTextFont(43);
  tdmes->SetTextSize(28);
  tdmes->Draw();
 
  TLine* lin0=new TLine(lowHypoFdOverPr,0.,highHypoFdOverPr,0.);
  lin0->SetLineStyle(2);
  lin0->SetLineColor(kGray+1);
  lin0->Draw();
  if(method==2 && optErrFD==2){
    c2->SaveAs(Form("%s-RcVsRcb_method%d_optErrFD%d_br%d.eps",mesName.Data(),method,optErrFD,correctForBR));
  }


  TCanvas* cfp=new TCanvas("cfp","fprompt",700,700);
  gPad->SetTickx();
  gPad->SetTicky();
  hfPromptCent->SetMinimum(0);
  hfPromptCent->SetMaximum(1.05);
  hfPromptCent->SetLineWidth(3);
  hfPromptCent->Draw();
  hfPromptMinNb->SetLineWidth(2);
  hfPromptMaxNb->SetLineWidth(2);
  hfPromptMinNb->SetLineStyle(2);
  hfPromptMaxNb->SetLineStyle(2);
  hfPromptMinNb->Draw("same");
  hfPromptMaxNb->Draw("same");
  hfPromptMinfc->SetLineWidth(2);
  hfPromptMaxfc->SetLineWidth(2);
  hfPromptMinfc->SetLineStyle(2);
  hfPromptMaxfc->SetLineStyle(2);
  hfPromptMinfc->SetLineColor(2);
  hfPromptMaxfc->SetLineColor(2);
  if(optErrFD==2){
    hfPromptMinfc->Draw("same");
    hfPromptMaxfc->Draw("same");
    TLatex* t1=new TLatex(0.17,0.15,"fc");
    t1->SetNDC();
    t1->SetTextColor(2);
    t1->Draw();
  }
  TLatex* t2=new TLatex(0.17,0.2,"Nb");
  t2->SetNDC();
  t2->Draw();
  cfp->SaveAs(Form("fprompt-%s.eps",mesName.Data()));

  hppC->SetMarkerStyle(markerppC);
  hppC->SetMarkerSize(msizppC);
  hppC->SetMarkerColor(colorppC);
  hppC->SetLineColor(colorppC);
  hppC->SetLineWidth(linewippC);
  hppC->SetLineStyle(linestppC);
  gppCsystdata->SetLineColor(colorSystDatappC);
  gppCsystdata->SetLineWidth(2);
  gppCsystdata->SetFillStyle(0);
  gppCsystdata->SetFillColor(colorSystDatappC);
  gppCsystFD->SetLineColor(colorSystFD);
  gppCsystFD->SetFillColor(colorSystFD);
  hAAC->SetMarkerStyle(markeraaC);
  hAAC->SetMarkerSize(msizaaC);
  hAAC->SetMarkerColor(coloraaC);
  hAAC->SetLineColor(coloraaC);
  hAAC->SetLineWidth(linewiaaC);
  hAAC->SetLineStyle(linestaaC);
  gaaCsystdata->SetLineColor(colorSystDataaaC);
  gaaCsystdata->SetLineWidth(2);
  gaaCsystdata->SetFillStyle(0);

  gaaCsystB->SetFillColor(colorSystRb);
  gaaCsystRb->SetFillColor(colorSystRb);
  gaaCsystFD->SetFillColor(colorSystFD);

  gaaCsystTot->SetLineColor(1);
  gaaCsystTot->SetLineWidth(3);
  gaaCsystTot->SetFillStyle(0);
  

  Float_t ymax=4.7*hppC->GetMaximum();
  Float_t ymin=999.;
  Float_t y16C=hAAC->GetBinContent(hAAC->FindBin(14.));
  if(y16C>0. && y16C<ymin) ymin=y16C;
  Float_t y12C=hAAC->GetBinContent(hAAC->FindBin(10.));
  if(y12C>0. && y12C<ymin) ymin=y12C;
  ymin*=0.255;
  ymax=19.6;
  ymin=7E-7;
  if(collSyst=="p-Pb" && TMath::Abs(normToCsec-1.)>0.001){
    ymax=110000.;
    ymin=1.1;
  }

  TH2F *hempty=new TH2F("hempty","",100,0.,binlim[nPtBins]*1.02,100,ymin,ymax);
  hempty->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hempty->GetYaxis()->SetTitle("d#it{N}/d#it{p}_{T}_{ }|_{ |#it{y}|<0.5} (1/GeV/#it{c})");
  hempty->GetYaxis()->SetTitleSize(0.05);
  hempty->GetXaxis()->SetTitleSize(0.05);
  hempty->GetYaxis()->SetTitleOffset(1.3);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->SetStats(0);
  TLatex* tdec =new TLatex(0.18,0.89,Form("%s",mesSymb.Data()));
  tdec->SetNDC();
  tdec->SetTextSize(0.038);
  tdec->SetTextFont(42);

  gStyle->SetMarkerSize(1.8);
  
  TCanvas* c3=new TCanvas("c3","Yield1Pad",750,800);
  c3->SetLeftMargin(0.15);
  c3->SetRightMargin(0.07);
  c3->SetBottomMargin(0.12);
  c3->SetTopMargin(0.07);
  c3->SetLogy();
  c3->SetTicky();
  hempty->Draw();
  gppCsystFD->Draw("e2same");
  gppCsystdata->DrawClone("e2same");
  gaaCsystFD->Draw("E2same");
  gaaCsystRb->Draw("E2same");
  gaaCsystdata->DrawClone("e2same");
  hppC->DrawClone("Psame");
  hAAC->DrawClone("Psame");
  //  TLegend* legC= new TLegend(0.59,0.77,0.89,0.91);
  TLegend* legC= new TLegend(0.53,0.55,0.89,0.69);
  //TLegend* legC= new TLegend(0.59,0.37,0.89,0.51);
  legC->SetHeader(Form("Centrality %s",centrality.Data()));
  legC->SetFillColor(0);
  legC->SetTextFont(42);
  ent=legC->AddEntry(hppC,"p-p rescaled reference","PL");
  ent=legC->AddEntry("","(#pm 6% norm. unc. not shown)","");
  //ent->SetTextColor(hppC->GetLineColor());
  ent=legC->AddEntry(hAAC,collSyst.Data(),"PL");
  //ent->SetTextColor(hAAC->GetLineColor());
  legC->Draw();
  //TLegend* legSy= new TLegend(0.18,0.16,0.45,0.32);
  //TLegend* legSy= new TLegend(0.62,0.57,0.89,0.75);
  TLegend* legSy= new TLegend(0.566,0.715,0.89,0.86);
  legSy->SetFillColor(0);
  legSy->SetLineColor(kGray+2);
  legSy->SetTextFont(42);
  ent=legSy->AddEntry(gppCsystdata,"Syst. unc. from Data","F");
  ent=legSy->AddEntry(gppCsystFD,"Syst. unc. from FONLL feed-down corr.","F");
  ent=legSy->AddEntry(gaaCsystRb,"Syst. unc. from #it{R}_{AA} feed-down","F");
  legSy->Draw();
  tdec->Draw();
  c3->SaveAs(Form("%s-Yields_1pad_method%d_optErrFD%d_br%d.eps",mesName.Data(),method,optErrFD,correctForBR));



  TH2F *hemptyr=new TH2F("hemptyr","",100,0.,binlim[nPtBins]*1.02,100,0.,1.7);
  hemptyr->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hemptyr->GetYaxis()->SetTitle(Form("#it{R}_{AA} (prompt %s)",mesSymb.Data()));
  hemptyr->GetYaxis()->SetTitleOffset(1.2);
  hemptyr->SetStats(0);

  TCanvas* c5=new TCanvas("c5","RAAcheck");
  c5->SetGridy();
  hRAAC->SetLineWidth(2);
  hRAAC->SetMarkerStyle(20);
  hRAAC->SetMarkerColor(1);
  hRAAC->SetLineColor(1);
  hRAAC->SetMinimum(0.);
  hRAAC->SetMaximum(1.2);
  hcheckRAAC->SetMarkerStyle(24);
  hcheckRAAC->SetMarkerSize(1.2);
  hcheckRAAC->SetMarkerColor(2);
  hraaCcheck2->SetMarkerStyle(25);
  hraaCcheck2->SetMarkerColor(6);
  hraaCcheck2->SetMarkerSize(1.5);

  graaC->SetFillColor(kGray);
  hemptyr->Draw();
  graaC->Draw("E2same");
  hRAAC->Draw("PEsame");
  hraaCcheck2->Draw("PSAME");
  hcheckRAAC->Draw("PSAME");
  TLegend* legr=new TLegend(0.7,0.7,0.89,0.89);
  ent=legr->AddEntry(hRAAC,Form("%s%%",centrality.Data()),"PL");
  legr->SetFillColor(0);
  legr->SetTextFont(42);
  legr->Draw();
  c5->Update();
  //  c5->SaveAs(Form("%s-RAA_check_method%d.eps",mesName.Data(),method));


  // TCanvas* c6=new TCanvas("c6","Rccheck");
  // c6->Divide(3,2);
  // for(Int_t ib=0; ib<nPtBins; ib++){
  //   c6->cd(ib+1);
  //   if(gcrbc[ib]->GetN()>0){
  //     hcheck[ib]->Draw("col");
  //     gcrbc[ib]->SetLineColor(colors[ib]);
  //     gcrbc[ib]->Draw("lsame");
  //   }
  // }


  TString type="Yield";
  if(TMath::Abs(normToCsec-1)>0.001) type="CrossSec";
  TFile* outfil=new TFile(Form("%s%s_method%d_fd%d_br%d.root",mesName.Data(),type.Data(),method,optErrFD,correctForBR),"recreate");
  hAAC->Write();
  gppCsystFD->Write();
  gppCsystdata->Write();
  gaaCsystFD->Write();
  gaaCsystRb->Write();
  gaaCsystB->Write();
  gaaCsystdata->Write();
  gaaCsystTot->Write();
  hppC->Write();
  hfPromptCent->Write();
  hfPromptMinNb->Write();
  hfPromptMaxNb->Write();
  hfPromptMinfc->Write();
  hfPromptMaxfc->Write();
  systematicsABcent->SetName("AliHFSystErrAA");
  systematicsABcent->Write();
  if(systematicsPP) systematicsPP->Write();
  outfil->Close();
}


//____________________________________________________________
Bool_t PbPbDataSyst(AliHFSystErr *syst, TH1D* heff, Double_t pt, Double_t &dataSystUp, Double_t &dataSystDown)
{

  Double_t err = syst->GetTotalSystErr(pt)*syst->GetTotalSystErr(pt);
  Int_t theBin=heff->FindBin(pt);
  Double_t errRel=heff->GetBinError(theBin)/heff->GetBinContent(theBin);
 
  err += (errRel*errRel);

  Double_t errDown = err ;
  Double_t errUp = err ;


  dataSystUp = TMath::Sqrt(errUp);
  dataSystDown = TMath::Sqrt(errDown);

  printf("Pt %f  Bin %d  Eff %f  RelErrEff %f  TotSyst +%f -%f\n",pt,theBin,heff->GetBinContent(theBin),errRel,dataSystUp,dataSystDown);

  return kTRUE;
}


