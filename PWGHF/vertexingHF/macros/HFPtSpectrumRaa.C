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

/* $Id$ */
///////////////////////////////////////////////////////////////////////////////
//
// Macro to compute the Raa, taking as inputs the output of the corrected yields
//  and the pp reference
//
// R_AB =  [ ( dsigma/dpt )_AB / sigma_AB ]  / <TAB> *  ( dsigma/dpt )_pp
//
// 
// Parameters: 
//      1. ppfile = reference pp in the same pt binning
//	2. ABfile = corrected AB yields 
// 	3. outfile = output file name
//	4. decay = decay as in HFSystErr class
//      5. sigmaABCINT1B = cross section for normalization (**USE THE SAME AS ON 2.**)
//	6. fdMethod = feed-down subtraction method (kNb, kfc)
//      7. cc = centrality class
//      8. Energy = colliding energy (k276,k55)
//      9. MinHypo  = minimum energy loss hypothesis (Default 1./3.)
//      10. MaxHypo = maximum energy loss hypothesis (Default 3.0)
//      11. MaxRb = maximum Raa(b) hypothesis (Default 6.0, won't do anything)
//      12. RbHypo : flag to decide whether the Eloss hypothesis is Rb or Rb/Rc
//      13. CentralHypo = central energy loss hypothesis, DEFAULT TO 1.0
//      14. isRaavsEP = flag to compute the Raa IN/OUT of plane, divides the reference by 2.0
//      15. ScaledAndExtrapRef: flag to tag scaled+reference pp-scaled-data
//
//  Complains to : Zaida Conesa del Valle
//
///////////////////////////////////////////////////////////////////////////////

enum centrality{ kpp, k07half, kpPb0100, k010, k1020, k020, k2040, k2030, k3040, k4050, k3050, k5060, k4060, k6080, k4080, k5080, k80100,kpPb010, kpPb020, kpPb2040, kpPb4060, kpPb60100 };
enum centestimator{ kV0M, kV0A, kZNA, kCL1 };
enum energy{ k276, k5dot023, k55 };
enum BFDSubtrMethod { kfc, kNb };
enum RaavsEP {kPhiIntegrated, kInPlane, kOutOfPlane};
enum rapidity{ kdefault, k08to04, k07to04, k04to01, k01to01, k01to04, k04to07, k04to08, k01to05 };
enum particularity{ kTopological, kLowPt, kPP7TeVPass4 };


Bool_t printout = false;
Double_t ptprintout = 1.5;
Double_t NormPPUnc = 0.035;
Double_t NormABUnc = 0.037;
Bool_t elossFDQuadSum = true;

//____________________________________________________________
Bool_t PbPbDataSyst(AliHFSystErr *syst, Double_t pt, Int_t cc, Double_t &dataSystUp, Double_t &dataSystDown);

//____________________________________________________________
Double_t ExtractFDSyst(Double_t total, Double_t fd) {
  // total^2 = data^2 + fd^2
  Double_t data2 = total*total - fd*fd ;
  return TMath::Sqrt( data2 );
}

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


//
//
// R_AB =  [ ( dsigma/dpt )_AB / sigma_AB ]  / <TAB> *  ( dsigma/dpt )_pp
//
//
//____________________________________________________________
void HFPtSpectrumRaa(const char *ppfile="HFPtSpectrum_D0Kpi_method2_rebinnedth_230311_newsigma.root",
		     const char *ABfile="HFPtSpectrum_D0Kpi_PbPbcuts_method2_rebinnedth_230311_newsigma.root",
		     const char *outfile="HFPtSpectrumRaa.root",
		     Int_t decay=1,
		     Double_t sigmaABCINT1B=54.e9,
		     Int_t fdMethod = kNb, Int_t cc=kpp, Int_t Energy=k276,
		     Double_t MinHypo=1./3., Double_t MaxHypo=3.0, Double_t MaxRb=6.0,
		     Bool_t isRbHypo=false, Double_t CentralHypo = 1.0,
		     Int_t ccestimator = kV0M,
		     Bool_t isUseTaaForRaa=true, const char *shadRbcFile="", Int_t nSigmaShad=3.0,
		     Int_t isRaavsEP=kPhiIntegrated, Bool_t isScaledAndExtrapRef=kFALSE,
		     Int_t rapiditySlice=kdefault, Int_t analysisSpeciality=kTopological)
{

  gROOT->Macro("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/LoadLibraries.C");

  //
  // Defining the TAB values for the given centrality class
  //
  Double_t Tab = 1., TabSyst = 0., A=207.2, B=207.2;
  if ( Energy!=k276 && Energy!=k5dot023) {
    printf("\n The Tab values for this cms energy have not yet been implemented, please do it ! \n");
    return;
  }
  if ( cc == kpp ){
    Tab = 1.;
  }
  // Values from Alberica's twiki:
  //   https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CentStudies
  if( (ccestimator == kV0M) && (Energy==k276) ) {
    if ( cc == k07half ) {
      Tab = 24.81; TabSyst = 0.8037;
    } else if ( cc == k010 ) {
      Tab = 23.48; TabSyst = 0.97;
    } else if ( cc == k1020 ) {
      Tab = 14.4318; TabSyst = 0.5733;
    } else if ( cc == k020 ) {
      Tab = 18.93; TabSyst = 0.74;
    } else if ( cc == k2040 ) {
      Tab = 6.86; TabSyst = 0.28;
    } else if ( cc == k2030 ) {
      Tab = 8.73769; TabSyst = 0.370219;
    } else if ( cc == k3040 ) {
      Tab = 5.02755; TabSyst = 0.22099;
    } else if ( cc == k4050 ) {
      Tab = 2.68327; TabSyst = 0.137073;
    } else if ( cc == k3050 ) {
      Tab = 3.87011; TabSyst = 0.183847;
    } else if ( cc == k4060 ) {
      Tab = 2.00;  TabSyst= 0.11;
    } else if ( cc == k4080 ) {
      Tab = 1.20451; TabSyst = 0.071843;
    } else if ( cc == k5060 ) {
      Tab = 1.32884; TabSyst = 0.0929536;
    } else if ( cc == k6080 ) {
      Tab = 0.419; TabSyst = 0.033;
    } else if ( cc == k5080 ) {
      Tab = 0.719; TabSyst = 0.054;
    } else if ( cc == k80100 ){
      Tab = 0.0690; TabSyst = 0.0062;
    }
  }
  if( (ccestimator == kV0M) && (Energy==k5dot023) ) {
      if ( cc == k3050 ) {
          Tab = 3.76; TabSyst = 0.13;
      }
  }


  // pPb Glauber (A. Toia)
  // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PACentStudies#Glauber_Calculations_with_sigma
  if( cc == kpPb0100 ){
    Tab = 0.098334; TabSyst = 0.0070679;
    A=207.2; B=1.;
  }
  else if( ccestimator == kV0A ){
    if ( cc == kpPb020 ) {
      Tab = 0.183; TabSyst = 0.006245;
    } else if ( cc == kpPb2040 ) {
      Tab = 0.134; TabSyst = 0.004899;
    } else if ( cc == kpPb4060 ) {
      Tab = 0.092; TabSyst = 0.004796;
    } else if ( cc == kpPb60100 ) {
      Tab = 0.041; TabSyst = 0.008832;
    }
  }
  else if( ccestimator == kZNA ){
    if ( cc== kpPb010 ){
      Tab = 0.17; TabSyst = 0.01275;
   else if ( cc == kpPb020 ) {
      Tab = 0.164; TabSyst = 0.010724;
    } else if ( cc == kpPb2040 ) {
      Tab = 0.137; TabSyst = 0.005099;
    } else if ( cc == kpPb4060 ) {
      Tab = 0.1011; TabSyst = 0.006;
    } else if ( cc == kpPb60100 ) {
      Tab = 0.0459; TabSyst = 0.003162;
    }
  }
  else if( ccestimator == kCL1 ){
    if ( cc == kpPb020 ) {
      Tab = 0.19; TabSyst = 0.007;
    } else if ( cc == kpPb2040 ) {
      Tab = 0.136; TabSyst = 0.005;
    } else if ( cc == kpPb4060 ) {
      Tab = 0.088; TabSyst = 0.005;
    } else if ( cc == kpPb60100 ) {
      Tab = 0.0369; TabSyst = 0.0085;
    }
  }

  //
  // Reading the pp file 
  //
  TFile * ppf = new TFile(ppfile,"read");
  TH1D * hSigmaPP;
  TGraphAsymmErrors * gSigmaPPSyst;
  TGraphAsymmErrors * gSigmaPPSystData = (TGraphAsymmErrors*)ppf->Get("gScaledDataSystData");
  TGraphAsymmErrors * gSigmaPPSystTheory = (TGraphAsymmErrors*)ppf->Get("gScaledDataSystExtrap");
  TGraphAsymmErrors * gSigmaPPSystFeedDown = (TGraphAsymmErrors*)ppf->Get("gScaledDataSystFeedDown");
  TH1I * hCombinedReferenceFlag=0x0;
  TGraphAsymmErrors * gReferenceFdSyst=0x0;
  if(isScaledAndExtrapRef){
    hCombinedReferenceFlag = (TH1I*)ppf->Get("hCombinedReferenceFlag");
    hSigmaPP = (TH1D*)ppf->Get("hReference");
    gSigmaPPSyst = (TGraphAsymmErrors*)ppf->Get("gReferenceSyst");
    gReferenceFdSyst = (TGraphAsymmErrors*)ppf->Get("gReferenceFdSyst");
  } else { 
    hSigmaPP = (TH1D*)ppf->Get("fhScaledData");
    gSigmaPPSyst = (TGraphAsymmErrors*)ppf->Get("gScaledData");
  }
  Double_t scalePPRefToMatchRapidityBin = 1.0;

  
  // Call the systematics uncertainty class for a given decay
  AliHFSystErr *systematicsPP = 0x0;
  if(ppf->Get("AliHFSystErr")){
    systematicsPP = (AliHFSystErr*)ppf->Get("AliHFSystErr");
    printf("   --> AliHFSystErr object for pp reference read from file. Object title = %s\n",systematicsPP->GetTitle());
  }else{
    systematicsPP = new AliHFSystErr();
    if(analysisSpeciality==kLowPt){
      systematicsPP->SetIsLowPtAnalysis(true);
    }
    if(analysisSpeciality==kPP7TeVPass4){
      systematicsPP->SetIsPass4Analysis(true);
    }
    systematicsPP->Init(decay);
    printf("   --> AliHFSystErr object for pp reference created based on macro arguments. Object title = %s\n",systematicsPP->GetTitle());
  }
  //
  // Reading the AB collisions file
  // 
  TFile * ABf = new TFile(ABfile,"read");
  TH1D *hSigmaAB = (TH1D*)ABf->Get("histoSigmaCorr");
  //  TH2D *hSigmaABRcb = (TH2D*)ABf->Get("histoSigmaCorrRcb");
  //  TGraphAsymmErrors * gSigmaABSyst = (TGraphAsymmErrors*)ABf->Get("gSigmaCorr");
  TGraphAsymmErrors * gSigmaABSystFeedDown = (TGraphAsymmErrors*)ABf->Get("gSigmaCorrConservative");
  TNtuple * nSigmaAB = (TNtuple*)ABf->Get("fnSigma");
  //
  TH1D *hMassAB = (TH1D*)ABf->Get("hRECpt");
  TH1D *hDirectEffptAB = (TH1D*)ABf->Get("hDirectEffpt");
  TH1D *histofcAB = (TH1D*)ABf->Get("histofc");
  //
  TH1D* fhStatUncEffcSigmaAB = (TH1D*)ABf->Get("fhStatUncEffcSigma");
  TH1D* fhStatUncEffbSigmaAB = (TH1D*)ABf->Get("fhStatUncEffbSigma");
  TH1D* fhStatUncEffcFDAB = (TH1D*)ABf->Get("fhStatUncEffcFD");
  TH1D* fhStatUncEffbFDAB = (TH1D*)ABf->Get("fhStatUncEffbFD");
  //
  TH1D* fhStatUncEffcSigmaAB_Raa = (TH1D*)fhStatUncEffcSigmaAB->Clone("fhStatUncEffcSigmaAB_Raa");
  TH1D* fhStatUncEffbSigmaAB_Raa = (TH1D*)fhStatUncEffbSigmaAB->Clone("fhStatUncEffbSigmaAB_Raa");
  TH1D* fhStatUncEffcFDAB_Raa = (TH1D*)fhStatUncEffcFDAB->Clone("fhStatUncEffcFDAB_Raa");
  TH1D* fhStatUncEffbFDAB_Raa = (TH1D*)fhStatUncEffbFDAB->Clone("fhStatUncEffbFDAB_Raa");
  fhStatUncEffcSigmaAB_Raa->Reset();
  fhStatUncEffbSigmaAB_Raa->Reset();
  fhStatUncEffcFDAB_Raa->Reset();
  fhStatUncEffbFDAB_Raa->Reset();
  fhStatUncEffcSigmaAB_Raa->SetName("fhStatUncEffcSigmaAB_Raa");
  fhStatUncEffbSigmaAB_Raa->SetName("fhStatUncEffbSigmaAB_Raa");
  fhStatUncEffcFDAB_Raa->SetName("fhStatUncEffcFDAB_Raa");
  fhStatUncEffbFDAB_Raa->SetName("fhStatUncEffvFDAB_Raa");


  //
  // Call the systematics uncertainty class for a given decay
  AliHFSystErr *systematicsAB = 0x0;
  if(ABf->Get("AliHFSystErr")){
    systematicsAB=(AliHFSystErr*)ABf->Get("AliHFSystErr");
    printf("   --> AliHFSystErr object for A-A read from HFPtSpectrum file. Object title = %s\n",systematicsAB->GetTitle());
  }else{
    systematicsAB = new AliHFSystErr();
    systematicsAB->SetCollisionType(1);
    systematicsAB->SetRunNumber(2016);// check this
    if(Energy==k276){
      if ( cc == k07half ) systematicsAB->SetCentrality("07half");
      else if ( cc == k010 ) systematicsAB->SetCentrality("010");
      else if ( cc == k1020 ) systematicsAB->SetCentrality("1020");
      else if ( cc == k020 ) systematicsAB->SetCentrality("020");
      else if ( cc == k2040 || cc == k2030 || cc == k3040 ) {
	systematicsAB->SetCentrality("2040");
	systematicsAB->SetIsPbPb2010EnergyScan(true);
      }
      else if ( cc == k4060 || cc == k4050 || cc == k5060 ) systematicsAB->SetCentrality("4060");
      else if ( cc == k6080 || cc == k5080 ) systematicsAB->SetCentrality("6080");
      else if ( cc == k4080 ) systematicsAB->SetCentrality("4080");
      else if ( cc == k3050 ) {
	if (isRaavsEP == kPhiIntegrated) systematicsAB->SetCentrality("4080");
	else if (isRaavsEP == kInPlane) systematicsAB->SetCentrality("3050InPlane");
	else if (isRaavsEP == kOutOfPlane) systematicsAB->SetCentrality("3050OutOfPlane");
      }
    } else if (Energy==k5dot023){
      if ( cc == k3050 ){
	systematicsAB->SetRunNumber(15);
	systematicsAB->SetCentrality("3050");
      }
    }
    //
    else if ( cc == kpPb0100 || cc == kpPb010 || cc == kpPb020 || cc == kpPb2040 || cc == kpPb4060 || cc == kpPb60100 ) {
      systematicsAB->SetCollisionType(2);
      // Rapidity slices
      if(rapiditySlice!=kdefault){
	systematicsAB->SetIspPb2011RapidityScan(true);
	TString rapidity="";
	switch(rapiditySlice) {
	case k08to04: rapidity="0804"; scalePPRefToMatchRapidityBin=(0.093+0.280)/1.0; break;
	case k07to04: rapidity="0804"; scalePPRefToMatchRapidityBin=0.280/1.0; break;
	case k04to01: rapidity="0401"; scalePPRefToMatchRapidityBin=0.284/1.0; break;
	case k01to01: rapidity="0101"; scalePPRefToMatchRapidityBin=0.191/1.0; break;
	case k01to04: rapidity="0104"; scalePPRefToMatchRapidityBin=0.288/1.0; break;
	case k04to07: rapidity="0408"; scalePPRefToMatchRapidityBin=0.288/1.0; break;
	case k04to08: rapidity="0408"; scalePPRefToMatchRapidityBin=(0.288+0.096)/1.0; break;
	case k01to05: rapidity="0401"; scalePPRefToMatchRapidityBin=0.4; break;
	}
	systematicsAB->SetRapidity(rapidity);
      }
      // Centrality slices
      if(ccestimator==kV0A) {
	if(cc == kpPb020) systematicsAB->SetCentrality("020V0A");
	else if(cc == kpPb2040) systematicsAB->SetCentrality("2040V0A");
	else if(cc == kpPb4060) systematicsAB->SetCentrality("4060V0A");
	else if(cc == kpPb60100) systematicsAB->SetCentrality("60100V0A");
      } else if (ccestimator==kZNA) {
        if(cc == kpPb010) systematicsAB->SetCentrality("010ZNA");
	else if(cc == kpPb020) systematicsAB->SetCentrality("020ZNA");
	else if(cc == kpPb2040) systematicsAB->SetCentrality("2040ZNA");
	else if(cc == kpPb4060) systematicsAB->SetCentrality("4060ZNA");
	else if(cc == kpPb60100) systematicsAB->SetCentrality("60100ZNA");
      } else if (ccestimator==kCL1) {
	if(cc == kpPb020) systematicsAB->SetCentrality("020CL1");
	else if(cc == kpPb2040) systematicsAB->SetCentrality("2040CL1");
	else if(cc == kpPb4060) systematicsAB->SetCentrality("4060CL1");
	else if(cc == kpPb60100) systematicsAB->SetCentrality("60100CL1");
      }else {
	if(!(cc == kpPb0100)) {
	  cout <<" Error on the pPb options"<<endl;
	  return;
	}
      }
    }
		}
    else { 
      cout << " Systematics not yet implemented " << endl;
      return;
    }
    if(analysisSpeciality==kLowPt){
      systematicsAB->SetIsLowPtAnalysis(true);
    }
		else if(analysisSpeciality==kBDT){
      systematicsAB->SetIsBDTAnalysis(true);
    }
    //
    systematicsAB->Init(decay);
    printf("   --> AliHFSystErr object for A-A created based on macro arguments. Object title = %s\n",systematicsAB->GetTitle());
  }
  //
  Int_t entries = nSigmaAB->GetEntries();
  Float_t pt=0., signal=0., Rb=0., Rcb=0., fcAB=0., yieldAB=0., sigmaAB=0., statUncSigmaAB=0., sigmaABMin=0.,sigmaABMax=0.;
  nSigmaAB->SetBranchAddress("pt",&pt);
  nSigmaAB->SetBranchAddress("Signal",&signal);
  if (fdMethod==kNb) nSigmaAB->SetBranchAddress("Rb",&Rb);
  else if (fdMethod==kfc) nSigmaAB->SetBranchAddress("Rcb",&Rcb);
  nSigmaAB->SetBranchAddress("fc",&fcAB);
  nSigmaAB->SetBranchAddress("Yield",&yieldAB);
  nSigmaAB->SetBranchAddress("Sigma",&sigmaAB);
  nSigmaAB->SetBranchAddress("SigmaStatUnc",&statUncSigmaAB);
  nSigmaAB->SetBranchAddress("SigmaMax",&sigmaABMax);
  nSigmaAB->SetBranchAddress("SigmaMin",&sigmaABMin);
	
  
  // define the binning
  Int_t nbins = hSigmaAB->GetNbinsX();
  Double_t binwidth = hSigmaAB->GetBinWidth(1);
  Double_t *limits = new Double_t[nbins+1];
  Double_t *binwidths = new Double_t[nbins];
  Double_t xlow=0.;
  for (Int_t i=1; i<=nbins; i++) {
    binwidth = hSigmaAB->GetBinWidth(i);
    xlow = hSigmaAB->GetBinLowEdge(i);
    limits[i-1] = xlow;
    binwidths[i-1] = binwidth;
  }
  limits[nbins] = xlow + binwidth;


  //
  // Read the shadowing file if given as input
  //
  Double_t centralRbcShad[nbins+1], minRbcShad[nbins+1], maxRbcShad[nbins+1];
  for(Int_t i=0; i<=nbins; i++) { centralRbcShad[i]=1.0; minRbcShad[i]=6.0; maxRbcShad[i]=0.0; }
  Bool_t isShadHypothesis = false;
  if( strcmp(shadRbcFile,"")!=0 ) {
    isShadHypothesis = true;
    cout<<endl<<">>  Beware, using the shadowing prediction file with an "<<nSigmaShad<<"*sigma <<"<<endl<<endl;
    TFile *fshad = new TFile(shadRbcFile,"read");
    if(!fshad){ cout <<" >> Shadowing file not properly opened!!!"<<endl<<endl; return;}
    // TH1D *hRbcShadCentral = (TH1D*)fshad->Get("hDfromBoverPromptD_Shadowing_central");
    // TH1D *hRbcShadMin = (TH1D*)fshad->Get("hDfromBoverPromptD_Shadowing_upper");
    // TH1D *hRbcShadMax = (TH1D*)fshad->Get("hDfromBoverPromptD_Shadowing_lower");
    TH1D *hRbcShadCentral = (TH1D*)fshad->Get("hDfromBoverDfromc_L0");
    TH1D *hRbcShadMin = (TH1D*)fshad->Get("hDfromBoverDfromc_L0");
    TH1D *hRbcShadMax = (TH1D*)fshad->Get("hDfromBoverDfromc_L1");
    if(!hRbcShadCentral || !hRbcShadMin || !hRbcShadMax) {
      cout<< endl <<">> Shadowing input histograms are not ok !! "<<endl<<endl;
      return;
    }
    //       nSigmaShad
    //    nSigmaShad
    for(Int_t i=1; i<=nbins; i++) {
      Double_t xpt = hSigmaAB->GetBinCenter(i);
      if(xpt>24) xpt = 20;
      centralRbcShad[i] = hRbcShadCentral->GetBinContent( hRbcShadCentral->FindBin(xpt) );
      Double_t minValue0 = hRbcShadMin->GetBinContent( hRbcShadMin->FindBin(xpt) );
      Double_t maxValue0 = hRbcShadMax->GetBinContent( hRbcShadMax->FindBin(xpt) );
      Double_t arrayEl[3] = {minValue0,maxValue0, centralRbcShad[i]};
      Double_t minValue = TMath::MinElement(3,arrayEl);
      Double_t maxValue = TMath::MaxElement(3,arrayEl);
      cout<<">> Shadowing pt="<<xpt<<"  central="<<centralRbcShad[i]<<"  min="<<minValue<<"  max="<<maxValue<<endl;
      if(minValue>centralRbcShad[i]){ minValue = centralRbcShad[i]; }
      if(maxValue<centralRbcShad[i]){ maxValue = centralRbcShad[i]; }
      minRbcShad[i] = centralRbcShad[i] - nSigmaShad*(centralRbcShad[i] - minValue);
      maxRbcShad[i] = centralRbcShad[i] + nSigmaShad*(maxValue - centralRbcShad[i]);
      cout<<">> Shadowing hypothesis pt="<<xpt<<"  central="<<centralRbcShad[i]<<"  min="<<minRbcShad[i]<<"  max="<<maxRbcShad[i]<<endl;
    }
  }

  //
  // define the bins correspondence bw histos/files/graphs
  //
  //
  TH2D * hRABvsRcb = new TH2D("hRABvsRcb"," R_{AB}(c) vs Rcb Eloss hypothesis; p_{T} [GeV/c] ;  R_{AB}(c) ; Rcb Eloss hypothesis  ",nbins,limits,800,0.,4.);
  TH2D * hRABvsRb = new TH2D("hRABvsRb"," R_{AB}(c) vs Rb Eloss hypothesis; p_{T} [GeV/c] ;  R_{AB}(c) ; Rb Eloss hypothesis ",nbins,limits,800,0.,4.);
  //  TH2D * hRABBeautyvsRCharm = new TH2D("hRABBeautyvsRCharm"," R_{AB}(c) vs Rb Eloss hypothesis; p_{T} [GeV/c] ;  R_{AB}(b) ;  R_{AB}(c) ",nbins,limits,800,0.,4.);
  Int_t nbinsHypo=800;//200;
  Double_t *limitsHypo = new Double_t[nbinsHypo+1];
  for(Int_t i=1; i<=nbinsHypo+1; i++) limitsHypo[i-1]= i*4./800.;
  TH3D * hRABCharmVsRBeautyVsPt = new TH3D("hRABCharmVsRBeautyVsPt"," R_{AB}(c) vs Rb vs p_{T} Eloss hypothesis; p_{T} [GeV/c] ;  R_{AB}(b) ;  R_{AB}(c) ",nbins,limits,nbinsHypo,limitsHypo,nbinsHypo,limitsHypo);
  TH2D *hRCharmVsRBeauty[nbins+1];
  for(Int_t i=0; i<=nbins; i++) hRCharmVsRBeauty[i] = new TH2D(Form("hRCharmVsRBeauty_%i",i),Form("RAB(c) vs RAB(b) for pt bin %i ; R_{AB}(b) ;  R_{AB}(c)",i),nbinsHypo,limitsHypo,nbinsHypo,limitsHypo);
  TH2D *hRCharmVsElossHypo[nbins+1];
  for(Int_t i=0; i<=nbins; i++) hRCharmVsElossHypo[i] = new TH2D(Form("hRCharmVsElossHypo_%i",i),Form("RAB(c) vs ElossHypo for pt bin %i ; Eloss Hypothesis (c/b) ;  R_{AB}(c)",i),nbinsHypo,limitsHypo,nbinsHypo,limitsHypo);
  //
  TH1D *hRABEloss00= new TH1D("hRABEloss00","hRABEloss00",nbins,limits);
  TH1D *hRABEloss05= new TH1D("hRABEloss05","hRABEloss05",nbins,limits);
  TH1D *hRABEloss10= new TH1D("hRABEloss10","hRABEloss10",nbins,limits);
  TH1D *hRABEloss15= new TH1D("hRABEloss15","hRABEloss15",nbins,limits);
  TH1D *hRABEloss20= new TH1D("hRABEloss20","hRABEloss20",nbins,limits);
  //
  TH2D * hRABvsRbFDlow = new TH2D("hRABvsRbFDlow"," R_{AB}(c) vs Rb Eloss hypothesis (FD low); p_{T} [GeV/c] ; Rb Eloss hypothesis ; R_{AB}(c) ",nbins,limits,800,0.,4.);
  TH2D * hRABvsRbFDhigh = new TH2D("hRABvsRbFDhigh"," R_{AB}(c) vs Rb Eloss hypothesis (FD high); p_{T} [GeV/c] ; Rb Eloss hypothesis ; R_{AB}(c) ",nbins,limits,800,0.,4.);
  //
  TH1D * hRABvsRbFDhigh_proj = new TH1D("hRABvsRbFDhigh_proj","hRABvsRbFDhigh_proj",nbins,limits);
  TH1D * hRABvsRbFDlow_proj = new TH1D("hRABvsRbFDlow_proj","hRABvsRbFDlow_proj",nbins,limits);
  //
  TNtuple *ntupleRAB=0x0 ;
  if (fdMethod==kNb) {
    ntupleRAB = new TNtuple("ntupleRAB","ntupleRAB (Nb)","pt:TAB:sigmaPP:sigmaAB:invyieldAB:invyieldABFDHigh:invyieldABFDLow:RABCharm:RABCharmFDHigh:RABCharmFDLow:RABBeauty:fc",100000);
  } else if (fdMethod==kfc) {
    ntupleRAB = new TNtuple("ntupleRAB","ntupleRAB (fc)","pt:TAB:sigmaPP:sigmaAB:invyieldAB:invyieldABFDHigh:invyieldABFDLow:Rcb:RABCharm:RABCharmFDHigh:RABCharmFDLow:RABBeauty:RABBeautyFDHigh:RABBeautyFDLow:fc",100000);
  }
  if(!ntupleRAB) printf("ERROR: Wrong method option");

  TH1D * hYieldABvsPt = new TH1D("hYieldABvsPt"," Yield_{AB}(c) vs p_{T} (no Eloss hypothesis); p_{T} [GeV/c] ; Yield_{charm} ",nbins,limits);
  TH1D * hRABvsPt = new TH1D("hRABvsPt"," R_{AB}(c) vs p_{T} (no Eloss hypothesis); p_{T} [GeV/c] ; R_{charm} ",nbins,limits);
  TH1D * hRABvsPt_DataSystematics = new TH1D("hRABvsPt_DataSystematics"," Systematics of R_{AB} (c) vs p_{T} (no Eloss hypothesis); p_{T} [GeV/c] ; R_{charm} ",nbins,limits);
  TGraphAsymmErrors *gRAB_ElossHypothesis = new TGraphAsymmErrors(nbins+1);
  gRAB_ElossHypothesis->SetNameTitle("gRAB_ElossHypothesis","RAB Eloss systematics");
  TGraphAsymmErrors *gRAB_FeedDownSystematics = new TGraphAsymmErrors(nbins+1);
  gRAB_FeedDownSystematics->SetNameTitle("gRAB_FeedDownSystematics","RAB Feed-Down systematics");
  TGraphAsymmErrors *gRAB_fcFeedDownOnly = new TGraphAsymmErrors(nbins+1);
  gRAB_fcFeedDownOnly->SetNameTitle("gRAB_fcFeedDownOnly","RAB fc Feed-Down Only");
  TGraphAsymmErrors *gRAB_FeedDownSystematicsElossHypothesis = new TGraphAsymmErrors(nbins+1);
  gRAB_FeedDownSystematicsElossHypothesis->SetNameTitle("gRAB_FeedDownSystematicsElossHypothesis","RAB Feed-Down systematics considering Eloss hypothesis");
  TGraphAsymmErrors *gRAB_DataSystematics = new TGraphAsymmErrors(nbins+1);
  gRAB_DataSystematics->SetNameTitle("gRAB_DataSystematics","RAB Measurement (no FD, no Eloss) systematics");
  TGraphAsymmErrors *gRAB_DataSystematicsPP = new TGraphAsymmErrors(nbins+1);
  gRAB_DataSystematicsPP->SetNameTitle("gRAB_DataSystematicsPP","RAB Measurement PP meas. systematics (data+scaling)");
  TGraphAsymmErrors *gRAB_DataSystematicsAB = new TGraphAsymmErrors(nbins+1);
  gRAB_DataSystematicsAB->SetNameTitle("gRAB_DataSystematicsAB","RAB Measurement AB (no FD, no Eloss, no PP data) systematics");
  TGraphAsymmErrors *gRAB_GlobalSystematics = new TGraphAsymmErrors(nbins+1);
  gRAB_GlobalSystematics->SetNameTitle("gRAB_GlobalSystematics","RAB Measurement global (data, FD, Eloss) systematics");
  Double_t ElossMax[nbins+1], ElossMin[nbins+1];
  for(Int_t i=0; i<=nbins; i++) { ElossMax[i]=0.; ElossMin[i]=6.; }
  Double_t fcElossMax[nbins+1], fcElossMin[nbins+1];
  for(Int_t i=0; i<=nbins; i++) { fcElossMax[i]=0.; fcElossMin[i]=6.; }
  Double_t FDElossMax[nbins+1], FDElossMin[nbins+1];
  for(Int_t i=0; i<=nbins; i++) { FDElossMax[i]=0.; FDElossMin[i]=6.; }

  TGraphAsymmErrors *gRAB_Norm = new TGraphAsymmErrors(1);
  gRAB_Norm->SetNameTitle("gRAB_Norm","RAB Normalization systematics (pp norm + Tab)");
  Double_t normUnc = TMath::Sqrt ( NormPPUnc*NormPPUnc + (TabSyst/Tab)*(TabSyst/Tab) );
  if(!isUseTaaForRaa) normUnc = TMath::Sqrt ( NormPPUnc*NormPPUnc + NormABUnc*NormABUnc );
  gRAB_Norm->SetPoint(1,0.5,1.);
  gRAB_Norm->SetPointError(1,0.25,0.25,normUnc,normUnc);

  //
  // R_AB =  ( dN/dpt )_AB  / <Ncoll_AB> *  ( dN/dpt )_pp ; <Ncoll> = <Tab> * sigma_NN^inel
  // R_AB =  [ ( dsigma/dpt )_AB / sigma_AB ]  / <TAB> *  ( dsigma/dpt )_pp
  //
  Int_t istartPPfd=0, istartPPsyst=0, istartABfd=0, istartPPextr=0;
  Double_t yPPh=0., yPPl=0., yABh=0., yABl=0.;
  Double_t RaaCharm =0., RaaBeauty=0.;
  Double_t RaaCharmFDhigh = 0., RaaCharmFDlow = 0.;
  Double_t RaaBeautyFDhigh = 0., RaaBeautyFDlow = 0.;
  Double_t systUp=0., systLow=0., systPPUp=0., systPPLow=0., systABUp=0., systABLow=0.;
  //
  //
  // Search the central value of the energy loss hypothesis Rb = Rc (bin)
  //
  Double_t ElossCentral[nbins+1];
  for(Int_t i=0; i<=nbins; i++) { ElossCentral[i]=0.; }
  //
  for(Int_t ientry=0; ientry<=entries; ientry++){

    nSigmaAB->GetEntry(ientry);
    //    cout << " pt="<< pt<<" sigma-AB="<<sigmaAB<<endl;
    if ( !(sigmaAB>0.) ) continue;
    //if(decay==2 && pt<2.) continue;

    // Compute RAB and the statistical uncertainty
    Int_t hppbin = hSigmaPP->FindBin( pt );
    Int_t hABbin = hSigmaAB->FindBin( pt );
    Double_t sigmapp = hSigmaPP->GetBinContent( hppbin );
    sigmapp *= scalePPRefToMatchRapidityBin; // scale to the proper rapidity bin width
    //    cout << " pt="<< pt<<", sigma-pp="<< sigmapp<<endl;
    if (isRaavsEP>0.) sigmapp = 0.5*sigmapp;
    if ( !(sigmapp>0.) ) continue;

    RaaCharm =  ( sigmaAB / sigmaABCINT1B ) / ((Tab*1e3) * sigmapp *1e-12 ) ;
    if(!isUseTaaForRaa) { 
      RaaCharm =  ( sigmaAB ) / ( (A*B) * sigmapp ) ;
    }

    if (fdMethod==kNb) {
      RaaBeauty = Rb ; 
    }
    else if (fdMethod==kfc) {
      RaaBeauty = ( RaaCharm / Rcb ) ;
    }

    Double_t ElossHypo = 0.;
    if (fdMethod==kfc) { ElossHypo = 1. / Rcb; }
    else  { ElossHypo = 1. / (RaaCharm / RaaBeauty) ; }
    if(isRbHypo) ElossHypo = RaaBeauty;

    // If using shadowing hypothesis, change the central hypothesis too
    if(isShadHypothesis) CentralHypo = centralRbcShad[hABbin];

    //    cout <<" pt "<< pt << " Raa charm " << RaaCharm << " Raa beauty " << RaaBeauty << " eloss hypo "<< ElossHypo<<endl; 
    //
    // Find the bin for the central Eloss hypo
    //
    if( TMath::Abs( ElossHypo - CentralHypo ) < 0.075 ){
      Double_t DeltaIni = TMath::Abs( ElossCentral[ hABbin ] - CentralHypo );
      Double_t DeltaV = TMath::Abs( ElossHypo - CentralHypo );
      //      cout << " pt " << pt << " ECentral " << ElossCentral[ hABbin ] << " Ehypo "<< ElossHypo ;
      if ( DeltaV < DeltaIni ) ElossCentral[ hABbin ] = ElossHypo;
      //      cout << " final ECentral " << ElossCentral[ hABbin ] << endl;
    }
  }
  //
  // Calculation of the Raa and its uncertainties
  //
  for(Int_t ientry=0; ientry<entries; ientry++){

    nSigmaAB->GetEntry(ientry);
    if ( !(sigmaAB>0.) ) continue;
    //    if ( pt<2 || pt>16) continue;


    // Compute RAB and the statistical uncertainty
    Int_t hppbin = hSigmaPP->FindBin( pt );
    Double_t sigmapp = hSigmaPP->GetBinContent( hppbin );
    if (isRaavsEP>0.) sigmapp = 0.5*sigmapp;
    sigmapp *= scalePPRefToMatchRapidityBin; // scale to the proper rapidity bin width
    if ( !(sigmapp>0.) ) continue;

    RaaCharm =  ( sigmaAB / sigmaABCINT1B ) / ((Tab*1e3) * sigmapp *1e-12 );
    if(!isUseTaaForRaa) {
      RaaCharm =  ( sigmaAB ) / ( (A*B) * sigmapp ) ;
    }

    // Flag to know if it is an scaled or extrapolated point of the pp reference
    Bool_t isExtrapolatedBin = kFALSE;
    if(isScaledAndExtrapRef && hCombinedReferenceFlag) isExtrapolatedBin = hCombinedReferenceFlag->GetBinContent( hppbin );
    istartPPsyst = -1;
    istartPPsyst = FindGraphBin(gSigmaPPSyst,pt);

    //
    // FONLL Feed-Down systematics
    //
    istartPPfd = -1;
    if(!isExtrapolatedBin) istartPPfd = FindGraphBin(gSigmaPPSystFeedDown,pt);
    istartABfd = -1;
    istartABfd = FindGraphBin(gSigmaABSystFeedDown,pt);

    //      cout << " Starting bin for pp is "<< istartPPfd <<", for AA is "<<istartABfd << endl;
    if(isExtrapolatedBin){
      if(gReferenceFdSyst){
	Int_t ibinfd = FindGraphBin(gReferenceFdSyst,pt);
	yPPh = gReferenceFdSyst->GetErrorYhigh(ibinfd);
	yPPl = gReferenceFdSyst->GetErrorYlow(ibinfd);
      }
    } else { 
      yPPh = gSigmaPPSystFeedDown->GetErrorYhigh(istartPPfd);
      yPPl = gSigmaPPSystFeedDown->GetErrorYlow(istartPPfd);
    }
    if (isRaavsEP>0.) {
      yPPh = yPPh*0.5;
      yPPl = yPPl*0.5;
    }
    yPPh *= scalePPRefToMatchRapidityBin; // scale to the proper rapidity bin width
    yPPl *= scalePPRefToMatchRapidityBin; // scale to the proper rapidity bin width

    yABh = gSigmaABSystFeedDown->GetErrorYhigh(istartABfd);
    yABl = gSigmaABSystFeedDown->GetErrorYlow(istartABfd);


    RaaCharmFDhigh = ( sigmaABMax / sigmaABCINT1B ) / ((Tab*1e3) * (sigmapp+yPPh) *1e-12 ) ;
    RaaCharmFDlow =  ( sigmaABMin / sigmaABCINT1B ) / ((Tab*1e3) * (sigmapp-yPPl) *1e-12 ) ;
    if(printout && TMath::Abs(ptprintout-pt)<0.1 ) cout << endl<<" pt "<< pt << " Raa " << RaaCharm <<" high "<< RaaCharmFDhigh << " low "<< RaaCharmFDlow<<endl;
    if(!isUseTaaForRaa) {
      RaaCharmFDhigh = ( sigmaABMax ) / ( (A*B)* (sigmapp+yPPh) ) ;
      RaaCharmFDlow =  ( sigmaABMin ) / ( (A*B)* (sigmapp-yPPl) ) ;
    }


    if (fdMethod==kNb) {
      RaaBeauty = Rb ; 
      RaaBeautyFDlow = Rb ;
      RaaBeautyFDhigh = Rb ;
      ntupleRAB->Fill( pt, Tab*1e3, sigmapp*1e-12, sigmaAB*1e-12, sigmaAB/sigmaABCINT1B,
		       sigmaABMax / sigmaABCINT1B, sigmaABMin / sigmaABCINT1B,
		       RaaCharm, RaaCharmFDhigh, RaaCharmFDlow, RaaBeauty, fcAB );
    }
    else if (fdMethod==kfc) {
      RaaBeauty = ( RaaCharm / Rcb ) ;
      RaaBeautyFDlow = ( RaaCharmFDlow / Rcb ) ;
      RaaBeautyFDhigh = ( RaaCharmFDhigh / Rcb ) ;
      hRABvsRcb->Fill( pt, RaaCharm, RaaBeauty );
      ntupleRAB->Fill( pt, Tab*1e3, sigmapp*1e-12, sigmaAB*1e-12, sigmaAB/sigmaABCINT1B,
		       sigmaABMax / sigmaABCINT1B, sigmaABMin / sigmaABCINT1B,
		       Rcb, RaaCharm, RaaCharmFDhigh, RaaCharmFDlow, RaaBeauty, RaaBeautyFDhigh, RaaBeautyFDlow, fcAB );
    }
    hRABvsRb->Fill( pt, RaaCharm, RaaBeauty );
    hRABvsRbFDlow->Fill( pt, RaaCharmFDlow, RaaBeautyFDlow );
    hRABvsRbFDhigh->Fill( pt, RaaCharmFDhigh, RaaBeautyFDhigh );
    if(printout && TMath::Abs(ptprintout-pt)<0.1) cout << " pt "<< pt << " Rb " << RaaBeauty <<" high "<< RaaBeautyFDhigh << " low "<< RaaBeautyFDlow <<endl;

    hRABCharmVsRBeautyVsPt->Fill( pt, RaaBeauty, RaaCharm );
    Int_t ptbin = hRABvsPt->FindBin( pt );
    hRCharmVsRBeauty[ptbin]->Fill( RaaBeauty, RaaCharm );
    hRCharmVsRBeauty[ptbin]->Fill( RaaBeautyFDlow, RaaCharmFDlow );
    hRCharmVsRBeauty[ptbin]->Fill( RaaBeautyFDhigh, RaaCharmFDhigh );


    if (fdMethod==kfc) {
      if( TMath::Abs(Rcb-0.015)<0.009 ) hRABEloss00->Fill(pt,RaaCharm);
      if( TMath::Abs(Rcb-0.5)<0.009 ) hRABEloss05->Fill(pt,RaaCharm);
      if( TMath::Abs(Rcb-1.0)<0.009 ) {
	hRABEloss10->Fill(pt,RaaCharm);
	hRABvsRbFDhigh_proj->Fill(pt,RaaCharmFDhigh);
	hRABvsRbFDlow_proj->Fill(pt,RaaCharmFDlow);
      }
      if( TMath::Abs(Rcb-1.5)<0.009 ) hRABEloss15->Fill(pt,RaaCharm);
      if( TMath::Abs(Rcb-2.0)<0.009 ) hRABEloss20->Fill(pt,RaaCharm);
    }
    else if (fdMethod==kNb) {
      if( TMath::Abs(RaaBeauty-0.015)<0.009 ) hRABEloss00->Fill(pt,RaaCharm);
      if( TMath::Abs(RaaBeauty-0.5)<0.009 ) hRABEloss05->Fill(pt,RaaCharm);
      if( TMath::Abs(RaaBeauty-1.0)<0.009 ) {
	hRABEloss10->Fill(pt,RaaCharm);
	hRABvsRbFDhigh_proj->Fill(pt,RaaCharmFDhigh);
	hRABvsRbFDlow_proj->Fill(pt,RaaCharmFDlow);
      }
      if( TMath::Abs(RaaBeauty-1.5)<0.009 ) hRABEloss15->Fill(pt,RaaCharm);
      if( TMath::Abs(RaaBeauty-2.0)<0.009 ) hRABEloss20->Fill(pt,RaaCharm);
    }


    Int_t hABbin = hMassAB->FindBin( pt );
    if(isShadHypothesis) CentralHypo = centralRbcShad[hABbin];

    if(printout && TMath::Abs(ptprintout-pt)<0.1)
    if ( fdMethod==kNb && TMath::Abs(Rb -CentralHypo)< 0.05) {
      cout << " pt "<< pt <<", at bin "<<hABbin<<endl;
      cout<<" entries "<<entries<<", i="<<ientry<<", pt="<<pt<<", Rb="<<Rb<<", Tab="<<Tab<<", sigmaAB="<<sigmaAB<<", sigmapp="<<sigmapp<<", Raacharm="<<RaaCharm<<", RaaBeauty="<<RaaBeauty<<endl;
      cout << "  AB  basis: mass "<< hMassAB->GetBinContent(hABbin)<<", eff "<< hDirectEffptAB->GetBinContent(hABbin)<<endl;
      cout<<"   FD low,  err low AB "<< (sigmaAB-sigmaABMin)<<"  err low PP "<< yPPl<<" Raacharm="<<RaaCharmFDlow<<", RaaBeauty="<<RaaBeautyFDlow<<endl;
      cout<<"   FD high, err high AB "<< (sigmaABMax-sigmaAB)<<"  err high PP "<< yPPh<<" Raacharm="<<RaaCharmFDhigh<<", RaaBeauty="<<RaaBeautyFDhigh<<endl;
    }
    if(printout && TMath::Abs(ptprintout-pt)<0.1)
    if ( fdMethod==kfc) if(TMath::Abs(Rcb -CentralHypo)< 0.05 ){
      cout << " pt "<< pt <<", at bin "<<hABbin<<endl;
      cout<<" entries "<<entries<<", i="<<ientry<<", pt="<<pt<<", Rcb="<<Rcb<<", Tab="<<Tab<<", sigmaAB="<<sigmaAB<<", sigmapp="<<sigmapp<<", Raacharm="<<RaaCharm<<", RaaBeauty="<<RaaBeauty<<endl;
      cout << "  AB  basis: mass "<< hMassAB->GetBinContent(hABbin)<<", eff "<< hDirectEffptAB->GetBinContent(hABbin)<<", fc "<<histofcAB->GetBinContent(hABbin)<< endl;       
      cout<<"   FD low,  err low AB "<< (sigmaAB-sigmaABMin)<<"  err low PP "<< yPPl<<" Raacharm="<<RaaCharmFDlow<<", RaaBeauty="<<RaaBeautyFDlow<<endl;
      cout<<"   FD high, err high AB "<< (sigmaABMax-sigmaAB)<<"  err high PP "<< yPPh<<" Raacharm="<<RaaCharmFDhigh<<", RaaBeauty="<<RaaBeautyFDhigh<<endl;
    }


    //
    // Fill in the global properties ?
    //
    Double_t ElossHypo = 0.;
    if (fdMethod==kfc) { ElossHypo = 1./ Rcb; }
    else  { ElossHypo = 1. / (RaaCharm / RaaBeauty); }
    if(isRbHypo) ElossHypo = RaaBeauty;
    hRCharmVsElossHypo[ptbin]->Fill( ElossHypo, RaaCharm );

   // If using shadowing hypothesis, change the limit hypothesis too
    if(isShadHypothesis) {
      MinHypo = minRbcShad[ hABbin ];
      MaxHypo = maxRbcShad[ hABbin ];
    }

    //    cout <<" pt "<< pt << " Raa charm " << RaaCharm << " Raa beauty " << RaaBeauty << " eloss hypo "<< ElossHypo
    if(ientry==0) cout<<" pt"<< pt<< " ElossCentral "<< ElossCentral[hABbin] << " min-hypo "<<MinHypo << " max-hypo "<<MaxHypo<<endl;

    //
    // Fill in histos charm (null Eloss hypothesis)
    //
    Double_t minFdSyst = 0., maxFdSyst = 0.;
    if ( ElossHypo == ElossCentral[ hABbin ] ) {

      //
      // Data stat uncertainty
      //
      Double_t sigmappStat = hSigmaPP->GetBinError( hppbin );
      if (isRaavsEP>0.) sigmappStat = sigmappStat*0.5;
      sigmappStat *= scalePPRefToMatchRapidityBin; // scale to the proper rapidity bin width
      Int_t hRABbin = hRABvsPt->FindBin( pt );
      Double_t stat = RaaCharm * TMath::Sqrt( (statUncSigmaAB/sigmaAB)*(statUncSigmaAB/sigmaAB) + 
					      (sigmappStat/sigmapp)*(sigmappStat/sigmapp) ) ;
      if ( RaaCharm==0 ) stat =0.;
      if ( RaaCharm>0 ) {
	hRABvsPt->SetBinContent( hRABbin, RaaCharm );
	hRABvsPt->SetBinError( hRABbin, stat );
	hYieldABvsPt->SetBinContent( hRABbin, sigmaAB/sigmaABCINT1B );
	hYieldABvsPt->SetBinError( hRABbin, statUncSigmaAB/sigmaABCINT1B );
	
	cout << "pt="<< pt<< " Raa " << RaaCharm << " stat unc. "<< stat <<
	  " sigma-pp "<< sigmapp <<" sigma-AB "<< sigmaAB<<endl;
	if(printout && TMath::Abs(ptprintout-pt)<0.1) {
	  cout << " Raa " << RaaCharm << " stat unc. "<< stat << " is "<< stat/RaaCharm * 100. <<
	    "%, stat-pp "<< sigmappStat/sigmapp*100. <<"% stat-AB "<< statUncSigmaAB/sigmaAB*100.<<"%"<<endl;
	}

	Double_t errstatEff = fhStatUncEffcSigmaAB->GetBinError( hRABbin );
	fhStatUncEffcSigmaAB_Raa->SetBinError( hRABbin, errstatEff*RaaCharm );
	errstatEff = fhStatUncEffbSigmaAB->GetBinError( hRABbin );
	fhStatUncEffbSigmaAB_Raa->SetBinError( hRABbin, errstatEff*RaaCharm );
	errstatEff = fhStatUncEffcFDAB->GetBinError( hRABbin );
	fhStatUncEffcFDAB_Raa->SetBinError( hRABbin, errstatEff*RaaCharm );
	errstatEff = fhStatUncEffbFDAB->GetBinError( hRABbin );
	fhStatUncEffbFDAB_Raa->SetBinError( hRABbin, errstatEff*RaaCharm );
      }


      //
      //
      // Data systematics (sigma syst-but FD + extrap) syst
      //
      //
      // Data syst: a) Syst in p-p 
      //
      Double_t ptwidth = hSigmaAB->GetBinWidth(hABbin) / 2. ;
      istartPPextr = -1;
      if(!isExtrapolatedBin) istartPPextr = FindGraphBin(gSigmaPPSystTheory,pt);

      Double_t dataPPUp=0., dataPPLow=0.;
      if(isExtrapolatedBin) {
	dataPPUp = gSigmaPPSyst->GetErrorYhigh(istartPPsyst);
	dataPPLow = gSigmaPPSyst->GetErrorYlow(istartPPsyst);
	systPPUp = dataPPUp;
	systPPLow = dataPPLow;
      } else { 
	dataPPUp = ExtractFDSyst( gSigmaPPSystData->GetErrorYhigh(istartPPextr), gSigmaPPSystFeedDown->GetErrorYhigh(istartPPfd) );
	dataPPLow = ExtractFDSyst( gSigmaPPSystData->GetErrorYlow(istartPPextr), gSigmaPPSystFeedDown->GetErrorYlow(istartPPfd) );
	systPPUp = TMath::Sqrt( dataPPUp*dataPPUp + gSigmaPPSystTheory->GetErrorYhigh(istartPPextr)*gSigmaPPSystTheory->GetErrorYhigh(istartPPextr) );
	systPPLow = TMath::Sqrt( dataPPLow*dataPPLow + gSigmaPPSystTheory->GetErrorYlow(istartPPextr)*gSigmaPPSystTheory->GetErrorYlow(istartPPextr) );
      }
      if (isRaavsEP>0.) {
	dataPPUp = dataPPUp*0.5;
	dataPPLow = dataPPLow*0.5;
	if(isExtrapolatedBin) {
	  systPPUp = dataPPUp;
	  systPPLow = dataPPLow;
	} else {  
	  systPPUp = TMath::Sqrt( dataPPUp*dataPPUp + 0.5*gSigmaPPSystTheory->GetErrorYhigh(istartPPextr)*0.5*gSigmaPPSystTheory->GetErrorYhigh(istartPPextr) );
	  systPPLow = TMath::Sqrt( dataPPLow*dataPPLow + 0.5*gSigmaPPSystTheory->GetErrorYlow(istartPPextr)*0.5*gSigmaPPSystTheory->GetErrorYlow(istartPPextr) );
	}
      }
      systPPUp *= scalePPRefToMatchRapidityBin;  // scale to the proper rapidity bin width
      systPPLow *= scalePPRefToMatchRapidityBin; // scale to the proper rapidity bin width


      if(printout && TMath::Abs(ptprintout-pt)<0.1) {
	cout << " pt : "<< pt<<" Syst-pp-data "<< dataPPUp/sigmapp << "%, ";
	if(!isExtrapolatedBin){
	  if (isRaavsEP>0.) cout <<" extr unc + "<< 0.5*gSigmaPPSystTheory->GetErrorYhigh(istartPPextr)/sigmapp <<" - "<< 0.5*gSigmaPPSystTheory->GetErrorYlow(istartPPextr)/sigmapp <<" %";
	  else cout <<" extr unc + "<< (gSigmaPPSystTheory->GetErrorYhigh(istartPPextr)*scalePPRefToMatchRapidityBin)/sigmapp <<" - "<< (gSigmaPPSystTheory->GetErrorYlow(istartPPextr)*scalePPRefToMatchRapidityBin)/sigmapp <<" %";
	}
	cout << endl;
      }

      //
      // Data syst: b) Syst in PbPb
      //
      Double_t dataSystUp=0., dataSystDown=0.;
      Bool_t PbPbDataSystOk = PbPbDataSyst(systematicsAB,pt,cc,dataSystUp,dataSystDown);
      if (!PbPbDataSystOk) { cout <<" There is some issue with the PbPb data systematics, please check and rerun"<<endl; return; }
      systABUp = sigmaAB * TMath::Sqrt( dataSystUp*dataSystUp + 
					(hDirectEffptAB->GetBinError(hABbin)/hDirectEffptAB->GetBinContent(hABbin))*(hDirectEffptAB->GetBinError(hABbin)/hDirectEffptAB->GetBinContent(hABbin)) );

      systABLow = sigmaAB * TMath::Sqrt( dataSystDown*dataSystDown + 
					(hDirectEffptAB->GetBinError(hABbin)/hDirectEffptAB->GetBinContent(hABbin))*(hDirectEffptAB->GetBinError(hABbin)/hDirectEffptAB->GetBinContent(hABbin)) );
      //
      // Data syst : c) combine pp & PbPb
      //
      systLow = sigmapp>0. ? 
	RaaCharm * TMath::Sqrt( (systABLow/sigmaAB)*(systABLow/sigmaAB) + (systPPUp/sigmapp)*(systPPUp/sigmapp) )
	: 0.;

      systUp = sigmapp>0. ? 
	RaaCharm * TMath::Sqrt( (systABUp/sigmaAB)*(systABUp/sigmaAB) + (systPPLow/sigmapp)*(systPPLow/sigmapp) )
	: 0.;
      if ( RaaCharm==0 ) { systPPUp =0.; systPPLow =0.; }
      
      //      if(printout) 
	cout << " Syst-pp-up "<< systPPUp/sigmapp <<"%, syst-pp-low "<< systPPLow/sigmapp <<"%, syst-AB-up "<<systABUp/sigmaAB<<"%, syst-AB-low "<<systABLow/sigmaAB<<"%, tot-syst-up "<<systUp/RaaCharm<<"%, tot-syst-low "<<systLow/RaaCharm<<"%"<<endl;

      if ( RaaCharm>0 ) {
	hRABvsPt_DataSystematics->SetBinContent( hRABbin, RaaCharm );
	hRABvsPt_DataSystematics->SetBinError( hRABbin, systUp );
	gRAB_DataSystematics->SetPoint( hABbin, pt, RaaCharm ); // i, x, y
	gRAB_DataSystematics->SetPointError( hABbin, ptwidth, ptwidth, systLow, systUp );
	gRAB_DataSystematics->SetPointEXlow(hABbin, 0.4); gRAB_DataSystematics->SetPointEXhigh(hABbin,0.4);
	gRAB_DataSystematicsPP->SetPoint( hABbin, pt, RaaCharm ); // i, x, y
	gRAB_DataSystematicsPP->SetPointError( hABbin, ptwidth, ptwidth, RaaCharm *(systPPUp/sigmapp), RaaCharm *systPPLow/sigmapp );
	gRAB_DataSystematicsAB->SetPoint( hABbin, pt, RaaCharm ); // i, x, y
	gRAB_DataSystematicsAB->SetPointError( hABbin, ptwidth, ptwidth, RaaCharm *systABLow/sigmaAB, RaaCharm *systABUp/sigmaAB );
      }

      //
      // Feed-down Systematics
      //
      Double_t FDL=0., FDH=0.;
      if ( RaaCharmFDhigh > RaaCharmFDlow ){
	FDH = RaaCharmFDhigh; FDL = RaaCharmFDlow;
      } else {
	FDL = RaaCharmFDhigh; FDH = RaaCharmFDlow;
      } 
      
      if(printout && TMath::Abs(ptprintout-pt)<0.1) cout<<" Raa "<<RaaCharm<<", Raa-fd-low "<<RaaCharmFDlow <<", Raa-fd-high "<<RaaCharmFDhigh <<endl;
      maxFdSyst = TMath::Abs(FDH - RaaCharm);
      minFdSyst = TMath::Abs(RaaCharm - FDL);
      if ( RaaCharm>0 ) {
	gRAB_FeedDownSystematics->SetPoint( hABbin, pt, RaaCharm ); // i, x, y
	gRAB_FeedDownSystematics->SetPointError( hABbin, 0.3, 0.3, minFdSyst, maxFdSyst ); // i, x, y
	gRAB_fcFeedDownOnly->SetPoint( hABbin, pt,fcAB );
	gRAB_fcFeedDownOnly->SetPointError(hABbin, 0.3, 0.3, fcAB-(sigmaABMin/sigmaAB*fcAB), (sigmaABMax/sigmaAB*fcAB)-fcAB );
      }
      
      //      if(printout) { 
	cout<<" FD syst  +"<< maxFdSyst/RaaCharm <<" - "<<minFdSyst/RaaCharm<<endl;
	cout<<"  fc = "<<fcAB<<", ("<< sigmaABMax/sigmaAB * fcAB <<","<< sigmaABMin/sigmaAB * fcAB <<")"<<endl;
	//      }

      //
      // Filling part of the Eloss scenarii information
      //
      if(RaaCharm>0 ) {
	gRAB_ElossHypothesis->SetPoint( hABbin, pt, RaaCharm ); // i, x, y
	gRAB_ElossHypothesis->SetPointEXlow( hABbin, ptwidth);
	gRAB_ElossHypothesis->SetPointEXhigh( hABbin, ptwidth);
	gRAB_FeedDownSystematicsElossHypothesis->SetPoint( hABbin, pt, RaaCharm ); // i, x, y
	gRAB_FeedDownSystematicsElossHypothesis->SetPointEXlow( hABbin, ptwidth);
	gRAB_FeedDownSystematicsElossHypothesis->SetPointEXhigh( hABbin, ptwidth);
	gRAB_GlobalSystematics->SetPoint( hABbin, pt, RaaCharm ); // i, x, y
	gRAB_GlobalSystematics->SetPointEXlow(hABbin,0.4); gRAB_GlobalSystematics->SetPointEXhigh(hABbin,0.4);
      }
    }

    //
    // Filling Eloss scenarii information
    //
    //    trick in case not fine enough Rb hypothesis to cope with the min/max range
    //    if( RaaCharm>0 && ( (ElossHypo >= MinHypo && ElossHypo <=MaxHypo) || ElossHypo == ElossCentral[ hABbin ] ) && RaaBeauty<=MaxRb ) {
    //      by default better not use it, to monitor when this happens (could affect results)
    if( RaaCharm>0 && ElossHypo >= MinHypo && ElossHypo <=MaxHypo && RaaBeauty<=MaxRb ) {

      Double_t Ehigh =  ElossMax[ hABbin ] ;
      Double_t Elow =  ElossMin[ hABbin ] ;
      if ( RaaCharm > Ehigh ) ElossMax[ hABbin ] = RaaCharm ;
      if ( RaaCharm < Elow ) ElossMin[ hABbin ] = RaaCharm ;
      if(printout && TMath::Abs(ptprintout-pt)<0.1) {
	cout<<" Hypothesis " << ElossHypo << " sigma-AB "<< sigmaAB <<", Raa "<< RaaCharm <<", Raa Eloss max "<< ElossMax[hABbin] <<" Raa Eloss min "<< ElossMin[hABbin] << " Rb="<< RaaBeauty <<endl;
	cout<<"  Rb="<< RaaBeauty <<" max "<< RaaBeautyFDhigh <<" min "<< RaaBeautyFDlow <<endl;
      }
      Double_t fcEhigh =  fcElossMax[ hABbin ] ;
      Double_t fcElow =  fcElossMin[ hABbin ] ;
      if ( fcAB > fcEhigh ) fcElossMax[ hABbin ] = fcAB ;
      if ( fcAB < fcElow ) fcElossMin[ hABbin ] = fcAB ;
      Double_t FDEhigh = FDElossMax[ hABbin ];
      Double_t FDEmin = FDElossMin[ hABbin ];
      Double_t RFDhigh = RaaCharmFDhigh>RaaCharmFDlow ? RaaCharmFDhigh : RaaCharmFDlow;
      Double_t RFDlow = RaaCharmFDlow<RaaCharmFDhigh ? RaaCharmFDlow : RaaCharmFDhigh;
      if ( RFDhigh > FDEhigh ) FDElossMax[ hABbin ] = RFDhigh ;
      if ( RFDlow < FDEmin ) FDElossMin[ hABbin ] = RFDlow ;
      if(printout && TMath::Abs(ptprintout-pt)<0.1) 
	cout<<" Hypothesis " << ElossHypo << " sigma-AB "<< sigmaAB <<", Raa FD-max Eloss max "<< FDElossMax[hABbin] <<" Raa FD-min Eloss min "<< FDElossMin[hABbin] <<endl;
    }


  }


  // Finish filling the y-uncertainties of the Eloss scenarii 
  for (Int_t ibin=0; ibin<=nbins; ibin++){
    Double_t ipt=0., value =0.;
    gRAB_ElossHypothesis->GetPoint(ibin,ipt,value);
    if(ipt<=0) continue;
    //
    // Uncertainty on Raa due to the Eloss hypothesis
    Double_t elossYhigh = TMath::Abs( ElossMax[ibin] - value );
    Double_t elossYlow = TMath::Abs( value - ElossMin[ibin] );
    gRAB_ElossHypothesis->SetPointEYhigh(ibin, elossYhigh );
    gRAB_ElossHypothesis->SetPointEYlow(ibin, elossYlow );
    gRAB_ElossHypothesis->SetPointEXhigh(ibin, 0.2);
    gRAB_ElossHypothesis->SetPointEXlow(ibin, 0.2);
    cout << " pt "<< ipt << " Raa "<< value <<" max "<< ElossMax[ibin] << " min " <<ElossMin[ibin] <<endl;
    cout<<" Eloss syst  +"<< elossYhigh <<" - "<< elossYlow <<endl; 
    //    cout << " fc max "<< fcElossMax[ibin] << " fc min " <<fcElossMin[ibin] <<endl;   
    //
    // Uncertainty on Raa due to the FD unc & Eloss hypothesis
    Double_t fdElossEYhigh = TMath::Abs( FDElossMax[ibin] - value );
    Double_t fdElossEYlow = TMath::Abs( value - FDElossMin[ibin] );
    if(elossFDQuadSum){
      Double_t fdEYhigh = gRAB_FeedDownSystematics->GetErrorYhigh(ibin);
      fdElossEYhigh = TMath::Sqrt( elossYhigh*elossYhigh + fdEYhigh*fdEYhigh );
      Double_t fdEYlow = gRAB_FeedDownSystematics->GetErrorYlow(ibin);
      fdElossEYlow = TMath::Sqrt( elossYlow*elossYlow + fdEYlow*fdEYlow );
    }
    gRAB_FeedDownSystematicsElossHypothesis->SetPointEYhigh(ibin, fdElossEYhigh );
    gRAB_FeedDownSystematicsElossHypothesis->SetPointEYlow(ibin, fdElossEYlow );
    gRAB_FeedDownSystematicsElossHypothesis->SetPointEXhigh(ibin, 0.25);
    gRAB_FeedDownSystematicsElossHypothesis->SetPointEXlow(ibin, 0.25);
    cout<<" FD & Eloss syst  +"<< fdElossEYhigh <<" - "<< fdElossEYlow 
	<<" = + "<< fdElossEYhigh/value <<" - "<< fdElossEYlow/value <<" %" <<endl; 
    //
    // All uncertainty on Raa (FD unc & Eloss + data)
    Double_t systdatal = gRAB_DataSystematics->GetErrorYlow(ibin);
    Double_t systdatah = gRAB_DataSystematics->GetErrorYhigh(ibin);
    Double_t systgbhUnc = TMath::Sqrt( systdatah*systdatah + fdElossEYhigh*fdElossEYhigh );
    Double_t systgblUnc = TMath::Sqrt( systdatal*systdatal + fdElossEYlow*fdElossEYlow );
    gRAB_GlobalSystematics->SetPointEYhigh(ibin,systgbhUnc);
    gRAB_GlobalSystematics->SetPointEYlow(ibin,systgblUnc);
    cout<<" Data syst  +"<< systdatah <<" - "<<  systdatal <<" = + "<< systdatah/value <<" - " <<  systdatal/value << " % "<<endl; 
    cout<<" Global syst  +"<< systgbhUnc <<" - "<<  systgblUnc << " = + "<< systgbhUnc/value <<" - "<<  systgblUnc/value << " %" <<endl; 
    //
  }
    
  cout<<endl<<"  Calculation finished, now drawing"<<endl<<endl;


  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);


  TCanvas *cRABvsRb = new TCanvas("RABvsRb","RAB vs Rb");
  hRABvsRb->Draw("colz");
  cRABvsRb->Update();

//  TCanvas *cRABvsRbvsPt = new TCanvas("cRABvsRbvsPt","RAB vs Rb vs pt");
//  hRABCharmVsRBeautyVsPt->Draw("lego3z");
//  cRABvsRbvsPt->Update();

    
  cout<< "    Drawing feed-down contribution"<<endl;
  TCanvas *cRABvsRbFDl = new TCanvas("RABvsRbFDl","RAB vs Rb (FD low)");
  hRABvsRbFDlow->Draw("cont4z");
  cRABvsRbFDl->Update();
  TCanvas *cRABvsRbFDh = new TCanvas("RABvsRbFDh","RAB vs Rb (FD high)");
  hRABvsRbFDhigh->Draw("cont4z");
  cRABvsRbFDh->Update();

  TCanvas * cSigmaABptEloss = new TCanvas("cSigmaABptEloss","SigmaAB vs pt, Eloss hypothesis");
  TH1D *hSigmaABEloss00= new TH1D("hSigmaABEloss00","hSigmaABEloss00",nbins,limits);
  TH1D *hSigmaABEloss05= new TH1D("hSigmaABEloss05","hSigmaABEloss05",nbins,limits);
  TH1D *hSigmaABEloss10= new TH1D("hSigmaABEloss10","hSigmaABEloss10",nbins,limits);
  TH1D *hSigmaABEloss15= new TH1D("hSigmaABEloss15","hSigmaABEloss15",nbins,limits);
  TH1D *hSigmaABEloss20= new TH1D("hSigmaABEloss20","hSigmaABEloss20",nbins,limits);

  delete [] limits;
  delete [] binwidths;

  for (Int_t i=0; i<=nSigmaAB->GetEntriesFast(); i++) {
    nSigmaAB->GetEntry(i);
    if (fdMethod==kfc) {
      if( TMath::Abs(Rcb-0.015)<0.009 ) hSigmaABEloss00->Fill(pt,sigmaAB);
      if( TMath::Abs(Rcb-0.5)<0.009 ) hSigmaABEloss05->Fill(pt,sigmaAB);
      if( TMath::Abs(Rcb-1.0)<0.009 ) hSigmaABEloss10->Fill(pt,sigmaAB);
      if( TMath::Abs(Rcb-1.5)<0.009 ) hSigmaABEloss15->Fill(pt,sigmaAB);
      if( TMath::Abs(Rcb-2.0)<0.009 ) hSigmaABEloss20->Fill(pt,sigmaAB);
    }
    else if (fdMethod==kNb) {
      if( TMath::Abs(Rb-0.015)<0.009 ) hSigmaABEloss00->Fill(pt,sigmaAB);
      if( TMath::Abs(Rb-0.5)<0.009 ) hSigmaABEloss05->Fill(pt,sigmaAB);
      if( TMath::Abs(Rb-1.0)<0.009 ) hSigmaABEloss10->Fill(pt,sigmaAB);
      if( TMath::Abs(Rb-1.5)<0.009 ) hSigmaABEloss15->Fill(pt,sigmaAB);
      if( TMath::Abs(Rb-2.0)<0.009 ) hSigmaABEloss20->Fill(pt,sigmaAB);
    }
  }
  hSigmaABEloss00->SetLineColor(2);
  hSigmaABEloss05->SetLineColor(3);
  hSigmaABEloss10->SetLineColor(4);
  hSigmaABEloss15->SetLineColor(kMagenta+1);
  hSigmaABEloss20->SetLineColor(kGreen+2);
  hSigmaABEloss00->SetMarkerStyle(22);
  hSigmaABEloss05->SetMarkerStyle(26);
  hSigmaABEloss10->SetMarkerStyle(20);
  hSigmaABEloss15->SetMarkerStyle(25);
  hSigmaABEloss20->SetMarkerStyle(21);
  if (fdMethod==kNb) {
    hSigmaABEloss05->Draw("ph");
    hSigmaABEloss10->Draw("phsame");
    hSigmaABEloss15->Draw("phsame");
    hSigmaABEloss20->Draw("phsame");
  }
  else {
    hSigmaABEloss20->Draw("p");
    hSigmaABEloss00->Draw("phsame");
    hSigmaABEloss05->Draw("phsame");
    hSigmaABEloss10->Draw("phsame");
    hSigmaABEloss15->Draw("phsame");
    hSigmaABEloss20->Draw("phsame");
  }
  TLegend *legrcb = new TLegend(0.8,0.8,0.95,0.9);
  legrcb->SetFillColor(0);
  legrcb->AddEntry(hSigmaABEloss00,"Rc/b=0.0","lp");
  legrcb->AddEntry(hSigmaABEloss05,"Rc/b=0.5","lp");
  legrcb->AddEntry(hSigmaABEloss10,"Rc/b=1.0","lp");
  legrcb->AddEntry(hSigmaABEloss15,"Rc/b=1.5","lp");
  legrcb->AddEntry(hSigmaABEloss20,"Rc/b=2.0","lp");
  legrcb->Draw();
  cSigmaABptEloss->Update();


  TCanvas * cRABptEloss = new TCanvas("cRABptEloss","RAB vs pt, Eloss hypothesis");
  hRABEloss00->SetLineColor(2);
  hRABEloss05->SetLineColor(3);
  hRABEloss10->SetLineColor(4);
  hRABEloss15->SetLineColor(kMagenta+1);
  hRABEloss20->SetLineColor(kGreen+2);
  hRABEloss00->SetMarkerStyle(22);
  hRABEloss05->SetMarkerStyle(26);
  hRABEloss10->SetMarkerStyle(20);
  hRABEloss15->SetMarkerStyle(25);
  hRABEloss20->SetMarkerStyle(21);
  if (fdMethod==kNb) {
    hRABEloss05->Draw("ph");
    hRABEloss10->Draw("phsame");
    hRABEloss15->Draw("phsame");
    hRABEloss20->Draw("phsame");
  }
  else {
    hRABEloss20->Draw("p");
    hRABEloss00->Draw("phsame");
    hRABEloss05->Draw("phsame");
    hRABEloss10->Draw("phsame");
    hRABEloss15->Draw("phsame");
    hRABEloss20->Draw("phsame");
  }
  legrcb = new TLegend(0.8,0.8,0.95,0.9);
  legrcb->SetFillColor(0);
  if (fdMethod==kfc) {
    legrcb->AddEntry(hRABEloss00,"Rc/b=0.0","lp");
    legrcb->AddEntry(hRABEloss05,"Rc/b=0.5","lp");
    legrcb->AddEntry(hRABEloss10,"Rc/b=1.0","lp");
    legrcb->AddEntry(hRABEloss15,"Rc/b=0.5","lp");
    legrcb->AddEntry(hRABEloss20,"Rc/b=2.0","lp");
  }
  else if (fdMethod==kNb) {
    legrcb->AddEntry(hRABEloss00,"Rb=0.0","lp");
    legrcb->AddEntry(hRABEloss05,"Rb=0.5","lp");
    legrcb->AddEntry(hRABEloss10,"Rb=1.0","lp");
    legrcb->AddEntry(hRABEloss15,"Rb=0.5","lp");
    legrcb->AddEntry(hRABEloss20,"Rb=2.0","lp");
  }	
  legrcb->Draw();
  cRABptEloss->Update();

    
  cout<< "    Drawing summary results"<<endl;
  TCanvas * cRABpt = new TCanvas("cRABpt","RAB vs pt, no hypothesis");
  hRABEloss10->Draw("");
  cRABpt->Update();

  TCanvas * cRABptFDUnc = new TCanvas("cRABptFDUnc","RAB vs pt, FD Uncertainties");
  hRABvsRbFDlow_proj->Draw("");
  hRABEloss10->Draw("phsame");
  hRABvsRbFDhigh_proj->SetLineColor(kMagenta+1);
  hRABvsRbFDhigh_proj->Draw("same");
  hRABvsRbFDlow_proj->SetLineColor(kGreen+2);
  hRABvsRbFDlow_proj->Draw("same");
  legrcb = new TLegend(0.8,0.8,0.95,0.9);
  legrcb->SetFillColor(0);
  legrcb->AddEntry(hRABEloss10,"FD Central","lp");
  legrcb->AddEntry(hRABvsRbFDhigh_proj,"FD Upper unc.","l");
  legrcb->AddEntry(hRABvsRbFDlow_proj,"FD Lower unc.","l");
  legrcb->Draw();
  cRABptFDUnc->Update();

  TCanvas *RaaPlot = new TCanvas("RaaPlot","RAB vs pt, plot all");
  RaaPlot->SetTopMargin(0.085);
  RaaPlot->SetBottomMargin(0.1);
  RaaPlot->SetTickx();
  RaaPlot->SetTicky();
  TH2D *hRaaCanvas = new TH2D("hRaaCanvas"," R_{AB}(c) vs p_{T} (no Eloss hypothesis); p_{t} [GeV/c] ; R_{AA} prompt D",40,0.,40.,100,0.,3.0);
  hRaaCanvas->GetXaxis()->SetTitleSize(0.05);
  hRaaCanvas->GetXaxis()->SetTitleOffset(0.9);
  hRaaCanvas->GetYaxis()->SetTitleSize(0.05);
  hRaaCanvas->GetYaxis()->SetTitleOffset(0.9);
  hRaaCanvas->Draw();
  gRAB_Norm->SetFillStyle(1001);
  gRAB_Norm->SetFillColor(kGray+2);
  gRAB_Norm->Draw("2");
  TLine *line = new TLine(0.0172415,1.0,40.,1.0);
  line->SetLineStyle(2);
  line->Draw();
  hRABvsPt->SetMarkerColor(kBlue);
  hRABvsPt->SetMarkerColor(kBlue);
  hRABvsPt->SetMarkerStyle(21);
  hRABvsPt->SetMarkerSize(1.1);
  hRABvsPt->SetLineWidth(2);
  hRABvsPt->Draw("psame");
  gRAB_DataSystematics->SetLineColor(kBlue);
  gRAB_DataSystematics->SetLineWidth(3);
  gRAB_DataSystematics->SetLineWidth(2);
  gRAB_DataSystematics->SetFillColor(kRed);
  gRAB_DataSystematics->SetFillStyle(0);
  gRAB_DataSystematics->Draw("2");
  gRAB_FeedDownSystematics->SetFillColor(kViolet+1);
  gRAB_FeedDownSystematics->SetFillStyle(1001);
  gRAB_FeedDownSystematics->Draw("2");
  gRAB_ElossHypothesis->SetLineColor(kMagenta-7);
  gRAB_ElossHypothesis->SetFillColor(kMagenta-7);
  gRAB_ElossHypothesis->SetFillStyle(1001);
  gRAB_ElossHypothesis->Draw("2");
  hRABvsPt->Draw("psame");
  gRAB_DataSystematics->Draw("2");
  legrcb = new TLegend(0.5517241,0.6504237,0.8520115,0.8728814,NULL,"brNDC");
  legrcb->SetBorderSize(0);
  legrcb->SetTextSize(0.03389831);
  legrcb->SetLineColor(1);
  legrcb->SetLineStyle(1);
  legrcb->SetLineWidth(1);
  legrcb->SetFillColor(0);
  legrcb->SetFillStyle(1001);
  if(cc==k020) legrcb->AddEntry(hRABvsPt,"R_{AA} 0-20% CC","pe");
  else if(cc==k4080) legrcb->AddEntry(hRABvsPt,"R_{AA} 40-80% CC","pe");
  else legrcb->AddEntry(hRABvsPt,"R_{AA} and stat. unc.","pe");
  legrcb->AddEntry(gRAB_DataSystematics,"Syst. from data","f");
  legrcb->AddEntry(gRAB_ElossHypothesis,"Syst. from R_{AA}(B)","f");
  legrcb->AddEntry(gRAB_FeedDownSystematics,"Syst. from B feed-down","f");
  legrcb->Draw();
  TLatex* tc;
  TString system = "Pb-Pb   #sqrt{s_{NN}}=2.76 TeV";
  if( cc==kpPb0100 || cc==kpPb020 || cc==kpPb2040 || cc==kpPb4060 || cc==kpPb60100 ) system = "p-Pb   #sqrt{s_{NN}}=5.023 TeV";
  if(decay==1) tc =new TLatex(0.18,0.82,Form("D^{0},  %s ",system.Data()));
  else if(decay==2) tc =new TLatex(0.18,0.82,Form("D^{+},  %s ",system.Data()));
  else if(decay==3) tc =new TLatex(0.18,0.82,Form("D^{*+},  %s ",system.Data()));
  else if(decay==4) tc =new TLatex(0.18,0.82,Form("D_{s}^{+},  %s ",system.Data()));
  else  tc =new TLatex(0.18,0.82,Form("any (?) D meson,  %s ",system.Data()));
  tc->SetNDC();
  tc->SetTextSize(0.038);
  tc->SetTextFont(42);
  tc->Draw();
  RaaPlot->Update();


  TCanvas *RaaPlotFDEloss = new TCanvas("RaaPlotFDEloss","RAB vs pt, plot FD & ElossUnc");
  RaaPlotFDEloss->SetTopMargin(0.085);
  RaaPlotFDEloss->SetBottomMargin(0.1);
  hRaaCanvas->Draw();
  line->Draw();
  hRABvsPt->Draw("psame");
  gRAB_FeedDownSystematics->SetFillColor(kViolet+1);
  gRAB_FeedDownSystematics->SetFillStyle(1001);
  gRAB_FeedDownSystematics->Draw("2");
  gRAB_ElossHypothesis->SetLineColor(kMagenta-7);
  gRAB_ElossHypothesis->SetFillColor(kMagenta-7);
  gRAB_ElossHypothesis->SetFillStyle(1001);
  gRAB_ElossHypothesis->Draw("2");
  gRAB_FeedDownSystematicsElossHypothesis->SetLineColor(kBlack);
  gRAB_FeedDownSystematicsElossHypothesis->SetFillStyle(0);
  gRAB_FeedDownSystematicsElossHypothesis->SetFillColor(kViolet+1);
  gRAB_FeedDownSystematicsElossHypothesis->Draw("2");
  hRABvsPt->Draw("psame");
  legrcb = new TLegend(0.6,0.6,0.9,0.9);
  legrcb->SetBorderSize(0);
  legrcb->SetTextSize(0.03389831);
  legrcb->SetLineColor(1);
  legrcb->SetLineStyle(1);
  legrcb->SetLineWidth(1);
  legrcb->SetFillColor(0);
  legrcb->SetFillStyle(1001);
  legrcb->AddEntry(hRABvsPt,"R_{PbPb} and stat. unc.","pe");
  legrcb->AddEntry(gRAB_ElossHypothesis,"Energy loss syst.","f");
  legrcb->AddEntry(gRAB_FeedDownSystematics,"Feed down syst.","f");
  legrcb->AddEntry(gRAB_FeedDownSystematicsElossHypothesis,"Feed down & Eloss syst.","f");
  legrcb->Draw();
  RaaPlotFDEloss->Update();
  

  TCanvas *RaaPlotGlob = new TCanvas("RaaPlotGlob","RAB vs pt, plot Global unc");
  RaaPlotGlob->SetTopMargin(0.085);
  RaaPlotGlob->SetBottomMargin(0.1);
  RaaPlotGlob->SetTickx();
  RaaPlotGlob->SetTicky();
  hRaaCanvas->Draw();
  line->Draw();
  hRABvsPt->Draw("psame");
  gRAB_DataSystematics->Draw("2");
  gRAB_FeedDownSystematicsElossHypothesis->Draw("2");
  gRAB_GlobalSystematics->SetLineColor(kRed);
  gRAB_GlobalSystematics->SetLineWidth(2);
  gRAB_GlobalSystematics->SetFillColor(kRed);
  gRAB_GlobalSystematics->SetFillStyle(3002);
  gRAB_GlobalSystematics->Draw("2");
  hRABvsPt->Draw("psame");
  legrcb = new TLegend(0.6,0.6,0.9,0.9);
  legrcb->SetBorderSize(0);
  legrcb->SetTextSize(0.03389831);
  legrcb->SetLineColor(1);
  legrcb->SetLineStyle(1);
  legrcb->SetLineWidth(1);
  legrcb->SetFillColor(0);
  legrcb->SetFillStyle(1001);
  legrcb->AddEntry(hRABvsPt,"R_{PbPb} and stat. unc.","pe");
  legrcb->AddEntry(gRAB_DataSystematics,"Data syst.","f");
  legrcb->AddEntry(gRAB_FeedDownSystematicsElossHypothesis,"Feed down & Eloss syst.","f");
  legrcb->AddEntry(gRAB_GlobalSystematics,"Global syst.","f");
  legrcb->Draw();
  RaaPlotGlob->Update();


  
  TCanvas *RaaPlotSimple = new TCanvas("RaaPlotSimple","RAB vs pt, plot Simple unc");
  RaaPlotSimple->SetTopMargin(0.085);
  RaaPlotSimple->SetBottomMargin(0.1);
  RaaPlotSimple->SetTickx();
  RaaPlotSimple->SetTicky();
  hRaaCanvas->Draw();
  line->Draw();
  hRABvsPt->Draw("psame");
  gRAB_GlobalSystematics->SetLineColor(kBlue);
  gRAB_GlobalSystematics->SetLineWidth(2);
  gRAB_GlobalSystematics->SetFillStyle(0);
  gRAB_GlobalSystematics->Draw("2");
  gRAB_Norm->Draw("2");
  hRABvsPt->Draw("psame");
  legrcb = new TLegend(0.5991379,0.6949153,0.8534483,0.8559322,NULL,"brNDC");
  legrcb->SetBorderSize(0);
  legrcb->SetTextSize(0.03389831);
  legrcb->SetLineColor(1);
  legrcb->SetLineStyle(1);
  legrcb->SetLineWidth(1);
  legrcb->SetFillColor(0);
  legrcb->SetFillStyle(1001); 
  if(cc==k020) legrcb->AddEntry(hRABvsPt,"R_{AA} 0-20% CC","pe");
  else if(cc==k4080) legrcb->AddEntry(hRABvsPt,"R_{AA} 40-80% CC","pe");
  else legrcb->AddEntry(hRABvsPt,"R_{AA} and stat. unc.","pe");
  legrcb->AddEntry(gRAB_GlobalSystematics,"Systematics","f");
  legrcb->Draw();
  tc->Draw();
  RaaPlotSimple->Update();


  TCanvas *c = new TCanvas("c","");
  systematicsAB->DrawErrors();
  c->Update();

  TCanvas *cStatUnc = new TCanvas("cStatUnc","stat unc");
  cStatUnc->Divide(2,2);
  cStatUnc->cd(1);
  fhStatUncEffcSigmaAB_Raa->Draw("e");
  cStatUnc->cd(2);
  fhStatUncEffbSigmaAB_Raa->Draw("e");
  cStatUnc->cd(3);
  fhStatUncEffcFDAB_Raa->Draw("e");
  cStatUnc->cd(4);
  fhStatUncEffbFDAB_Raa->Draw("e");
  cStatUnc->Update();

  //
  // Write the information to a file
  //
  cout<<endl<< "  Save results in the output file"<<endl<<endl;
  TFile * out = new TFile(outfile,"recreate");

  ntupleRAB->Write();
  hRABvsRcb->Write();
  hRABvsRb->Write();
  hRABCharmVsRBeautyVsPt->Write();
  for(Int_t j=0; j<=nbins; j++) hRCharmVsRBeauty[j]->Write();
  //  for(Int_t j=0; j<=nbins; j++) hRCharmVsElossHypo[j]->Write();
  hRABvsPt->Write();
  hRABvsPt_DataSystematics->Write();
  gRAB_ElossHypothesis->Write();
  gRAB_FeedDownSystematics->Write();
  gRAB_fcFeedDownOnly->Write();
  gRAB_DataSystematics->Write();
  gRAB_DataSystematicsPP->Write();
  gSigmaPPSystTheory->Write();
  gRAB_DataSystematicsAB->Write();
  gRAB_Norm->Write();
  gRAB_FeedDownSystematicsElossHypothesis->Write();
  gRAB_GlobalSystematics->Write();
  if(isScaledAndExtrapRef && hCombinedReferenceFlag) hCombinedReferenceFlag->Write();
  systematicsPP->SetName("AliHFSystErrPP");
  systematicsPP->Write();
  systematicsAB->SetName("AliHFSystErrAA");
  systematicsAB->Write();
  out->Write();

}

//____________________________________________________________
Bool_t PbPbDataSyst(AliHFSystErr *syst, Double_t pt, Int_t cc, Double_t &dataSystUp, Double_t &dataSystDown)
{

  Double_t err=0., errUp=1., errDown=1.;
  Bool_t isOk=false;
  Double_t pidunc=0.;

  err = syst->GetTotalSystErr(pt)*syst->GetTotalSystErr(pt);
  errUp = err ;
  errDown = err ;

  // Apply an asymetric PID uncertainty for 2010 PbPb data only
  if( syst->GetRunNumber()==10 && syst->GetCollisionType()==1 ) {
    if( cc==k07half || cc==k020 || cc==k010 || cc==k1020 || cc==k2040 ) {
      if(pt<6) pidunc = 0.15;
      else pidunc = 0.05;
      errUp = err + pidunc*pidunc - syst->GetPIDEffErr(pt)*syst->GetPIDEffErr(pt);
      isOk = true;
    }
    else if ( cc==k3050 || cc==k4080 || cc==k4060 || cc==k6080 ){
      if(pt<3.1) pidunc = 0.10;
      else pidunc = 0.05;
      errUp = err + pidunc*pidunc - syst->GetPIDEffErr(pt)*syst->GetPIDEffErr(pt);
      isOk = true;
    }
  }
  else { isOk = true; }

  dataSystUp = TMath::Sqrt(errUp);
  dataSystDown = TMath::Sqrt(errDown);

  return isOk;
}
