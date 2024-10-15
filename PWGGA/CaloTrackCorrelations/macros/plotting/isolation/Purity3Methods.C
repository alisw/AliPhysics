///
/// \file Purity3Methods.C
/// \ingroup CaloTrackCorrMacrosPlotting
/// \brief Isolated photon purity calculation with 3 methods
///
/// Example macro to calculate the isolated photon purity via 3 methods
///    - ABCD ALICE
///    - ABCD ATLAS
///    - Template
///    Similar to the one used for the calculation in pp and Pb-Pb collisions in paper  https://alice-publications.web.cern.ch/node/11329
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

// ROOT includes
#include <TSystem.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TChain.h>
#include <TGraph.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TPDF.h>
#include <TPostScript.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TMath.h>
#include <TF1.h>
#include <TLegend.h>
#include <TPad.h>
#include <THistPainter.h>
#include <TCanvas.h>
#include <TDirectoryFile.h>
#include <TFractionFitter.h>
#include "TCustomBinning.h"

/// Triggered data samples, minimum bias or EMCal triggered L1-gamma or sum
TString trigger [] = {"MB","L1","L1MB"};

/// Available periods, calculation separate
TString period  [] = {"LHC17q","LHC17p"};

/// Final spectra and purity pT binning
const Int_t nBins = 19;
Double_t binPt[] = { 5., 6., 7., 8.,9.,10.,11.,12.,14.,16.,18.,20.,25.,30.,40.,60.,80.,100.,140.,180};

/// Apply pT dependent cut on shower shape selection
///
TH1F* M02PtDependentCut
(
 TH3F*   h3source,
 Float_t m02Min   = 0.1,
 Float_t m02Max   = 0.3,
 Bool_t  ptDepMin = kFALSE,
 Bool_t  ptDepMax = kTRUE,
 Float_t gapMax   = 0.0,
 Float_t gapMin   = 0.1,
 Float_t maxWidth = 0.,
 TH2F*   hsource  = 0
 )
{
  if ( !hsource && h3source )
    hsource = (TH2F*)h3source->Project3D("yx");
  TH2F* temp = (TH2F*) hsource->Clone("temp"); temp->Reset();
  
//  printf("Dep cut setting: Max %2.2f, Min %2.2f, "
//         "dep min %d, max %d, gapMax %2.2f, gapMin %2.2f \n",
//         m02Min,m02Max,
//         ptDepMin,ptDepMax,
//         gapMax,gapMin);
  
  for(int i=1 ; i<(hsource->GetNbinsX())+1 ; ++i)
  {
    double pt = hsource->GetXaxis()->GetBinCenter(i);
    
    Float_t max = m02Max;
    if ( ptDepMax )
    {
      max = 0.6-0.016*pt+gapMax+maxWidth;
      if ( m02Max > max ) max = m02Max;
    }
    
    Float_t min = m02Min;
    if ( ptDepMin )
    {
      min = 0.6-0.016*pt+gapMin;
      if ( m02Min > min ) min = m02Min;
    }
    
    //printf("\t pt %2.2f, min %2.4f < M02 <  max %2.4f\n",pt,min,max);
    
    for(int ii=1; ii<(hsource->GetNbinsY()); ++ii)
    {
      if ( min < hsource->GetYaxis()->GetBinCenter(ii)  &&  hsource->GetYaxis()->GetBinCenter(ii) < max )
      {
//        printf("\t selected center %2.4f width/2 %2.4f\n",
//               hsource->GetYaxis()->GetBinCenter(ii),
//               hsource->GetYaxis()->GetBinWidth(ii)/2);
        temp ->SetBinContent(i,ii,hsource->GetBinContent(i,ii));
        temp ->SetBinError(i,ii,hsource->GetBinError(i,ii));
      }
    }
  }
  
  TH1F* htarget = (TH1F*)temp->ProjectionX("htarget",1,-1) ;//Clone("htarget");
  return htarget;
}


/// Method to calculate the error of a fraction
///
static Double_t GetFractionError
(Float_t num   , Float_t den,
 Float_t numErr, Float_t denErr)
{
  if ( num == 0 || den == 0 ) return 0.;
  
  Double_t err = num/den * TMath::Sqrt( ( numErr * numErr ) / ( num * num ) +
                                        ( denErr * denErr ) / ( den * den )   );
  
  //printf("\t num %e den %e numErr %e denErr %e, num/den %f, err %e\n",num,den,numErr,denErr,num/den,err);

  return err;
}

/// Method to scale a histogram bin by its size.
/// Used after histogram rebinning.
///
static void ScaleBinBySize(TH1F* h)
{
  for(Int_t ibin = 1; ibin <= h->GetNbinsX();ibin++)
  {
    Double_t width   = h->GetBinWidth(ibin);
    Double_t content = h->GetBinContent(ibin);
    //Double_t center  = h->GetBinCenter(ibin);
    Double_t error   = h->GetBinError(ibin);
    
    //if     (center > 1.0) width/=3;
    //else if(center > 0.5) width/=2.25;

    //printf("bin %d, width %f, content %e\n",ibin,width,content);
    h->SetBinContent(ibin,content/width);
    h->SetBinError  (ibin,  error/width);
  }
}

/// Temporal containers for ABCD cluster histograms used in ProjectABCD
///
TH1F * hAA = 0;
TH1F * hBB = 0;
TH1F * hCC = 0;
TH1F * hDD = 0;

/// For a TH3 histogram containing candidate cluster pt vs shower shape vs isolation momentum
/// Obtain the pT spectra for the ABCD regions put them temporary in histograms hAA, hBB, hCC and hDD
/// Optionally correct the cluste spectra by the trigger efficiency for EMCal triggered data
///
void ProjectABCD(const Float_t sumSigMin  , const Float_t sumSigMax, const Float_t sumBkgMin, const Float_t sumBkgMax,
                 const Float_t  shSigMin  , const Float_t  shSigMax, const Float_t  ssBkgMin, const Float_t  ssBkgMax,
                 const Int_t    shPtDep, const Float_t cone,
                 const TString tagIsoCone , TFile* fileEffTrig, TH3F * h33)
{
  //printf("******* %s\n",tagIsoCone.Data());
  Float_t ssBkgGap   = ssBkgMin-shSigMax;
  Float_t ssBkgWidth = ssBkgMax-ssBkgMin;

//  printf("Project:\n");
//  printf("\t Sig: M02 %1.2f-%1.2f; Iso %1.2f-%1.2f\n", shSigMin, shSigMax, sumSigMin, sumSigMax);
//  printf("\t Bkg: M02 %1.2f-%1.2f; Iso %1.2f-%1.2f\n", ssBkgMin, ssBkgMax, sumBkgMin, sumBkgMax);

  // A, isolated, photon shape
  //printf(">>> A\n");
  h33->SetAxisRange(sumSigMin, sumSigMax-0.001,"Z"); // Sum pt Cone
  if ( !shPtDep )
  {
    h33->SetAxisRange(shSigMin , shSigMax-0.001 ,"Y"); // M02
    hAA = (TH1F*) h33->Project3D("x");
  }
  else  if ( shPtDep == 1 )
  {
    h33->SetAxisRange(shSigMin , 70,"Y"); // M02, reopen
    hAA = M02PtDependentCut(h33,shSigMin,shSigMax,0,1);
  }
 
  hAA->SetName(Form("A_%s",tagIsoCone.Data()));
  //hA[iprod] = (TH1F*) hA[iprod]->Rebin(nBins,"",binPt);
  TH1F* tmpA = hAA;
  hAA = (TH1F*) tmpA->Rebin(nBins,"",binPt);
  tmpA->Delete();
  
  // B, isolated, bkg shape
  //printf(">>> B\n");
  h33->SetAxisRange(sumSigMin, sumSigMax-0.001,"Z"); // Sum pt Cone
  if ( !shPtDep )
  {
    h33->SetAxisRange(ssBkgMin , ssBkgMax-0.001,"Y"); // M02
    hBB = (TH1F*) h33->Project3D("x");
  }
  else  if ( shPtDep == 1 )
  {
    h33->SetAxisRange(0 , ssBkgMax-0.001,"Y"); // M02, reopen
    hBB = M02PtDependentCut(h33, ssBkgMin, ssBkgMax, 1, 1, ssBkgGap, ssBkgGap, ssBkgWidth);
  }

  hBB->SetName(Form("B_%s",tagIsoCone.Data()));
  //hBB = (TH1F*) hBB->Rebin(nBins,"",binPt);
  TH1F* tmpB = hBB;
  hBB = (TH1F*) tmpB->Rebin(nBins,"",binPt);
  tmpB->Delete();
  
  // C, not isolated, photon shape
  //printf(">>> C\n");
  h33->SetAxisRange(sumBkgMin, sumBkgMax-0.001,"Z"); // Sum pt Cone
  if ( !shPtDep )
  {
    h33->SetAxisRange(shSigMin , shSigMax-0.001,"Y"); // M02
    hCC = (TH1F*) h33->Project3D("x");
  }
  else  if ( shPtDep == 1 )
  {
    h33->SetAxisRange(shSigMin , 70,"Y"); // M02, reopen
    hCC = M02PtDependentCut(h33,shSigMin, shSigMax,0,1);
  }

  hCC->SetName(Form("C_%s",tagIsoCone.Data()));
  //hCC = (TH1F*) hCC->Rebin(nBins,"",binPt);
  TH1F* tmpC = hCC;
  hCC = (TH1F*) tmpC->Rebin(nBins,"",binPt);
  tmpC->Delete();
  
  // D, not isolated, bkg shape
  //printf(">>> D\n");
  h33->SetAxisRange(sumBkgMin, sumBkgMax-0.001,"Z"); // Sum pt Cone
  if ( !shPtDep )
  {
    h33->SetAxisRange(ssBkgMin , ssBkgMax-0.001,"Y"); // M02
    hDD = (TH1F*) h33->Project3D("x");
  }
  else  if ( shPtDep == 1 )
  {
    h33->SetAxisRange(0 , ssBkgMax-0.001,"Y"); // M02, reopen
    hDD = M02PtDependentCut(h33, ssBkgMin, ssBkgMax, 1, 1, ssBkgGap, ssBkgGap, ssBkgWidth);
  }

  
  hDD->SetName(Form("D_%s",tagIsoCone.Data()));
  //hDD = (TH1F*) hDD->Rebin(nBins,"",binPt);
  TH1F* tmpD = hDD;
  hDD = (TH1F*) tmpD->Rebin(nBins,"",binPt);
  tmpD->Delete();
  
  // Apply trigger efficiency
  //
  if ( fileEffTrig )
  {
    TString endOfHistoName = Form("_R%1.2f_L1",cone); // low below 16 GeV, high above
    //TString endOfHistoName = Form("_R%1.2f_L1_JJlow",cone);
    //TString endOfHistoName = Form("_R%1.2f_L1_JJhigh",cone);
    TH1F* hEffA = (TH1F*) fileEffTrig->Get("hIsoNarrow"+endOfHistoName);
    TH1F* hEffB = (TH1F*) fileEffTrig->Get("hIsoWide"+endOfHistoName);
    TH1F* hEffC = (TH1F*) fileEffTrig->Get("hNoIsoNarrow"+endOfHistoName);
    TH1F* hEffD = (TH1F*) fileEffTrig->Get("hNoIsoWide"+endOfHistoName);
    
    if(!hEffA || !hEffB || !hEffC || !hEffD)
      printf("Eff histo %p, %p, %p, %p\n", hEffA, hEffB, hEffC, hEffD);
    
    hAA->Divide(hEffA);
    hBB->Divide(hEffB);
    hCC->Divide(hEffC);
    hDD->Divide(hEffD);
  }
  
}

/// Calculate the purity using ATLAS ABCD method
/// Done per bin.
///
void GetPurityLeakInBin( Double_t nA, Double_t nAe, Double_t nB, Double_t nBe, Double_t nC, Double_t nCe, Double_t nD, Double_t nDe,
                         Double_t cB, Double_t cBe, Double_t cC, Double_t cCe, Double_t cD, Double_t cDe,
                         Float_t rBkg, Float_t rBkge, Float_t & purity, Float_t & purityE)
{
  if ( nA == 0 ) return;

  //cB = 0.5; cD = 0.2; cC = 0.1;
  Double_t a = cB*cC*rBkg - cD;
  Double_t b = nD  + cD*nA - nB*cC*rBkg - cB * nC * rBkg;
  Double_t c = nB * nC * rBkg  - nD * nA ;
  
  Double_t a1E =  cB*cC*rBkg * TMath::Sqrt( cBe*cBe/(cB*cB) + cCe*cCe/(cC*cC) + rBkge*rBkge/(rBkg*rBkg) ) ;
  Double_t aE  =  TMath::Sqrt( cDe*cDe + a1E*a1E );
  
  Double_t b1E =  cD * nA        * TMath::Sqrt( cDe*cDe/(cD*cD) + nAe*nAe/(nA*nA) );
  Double_t b2E =  cC * nB * rBkg * TMath::Sqrt( cCe*cCe/(cC*cC) + rBkge*rBkge/(rBkg*rBkg) + nBe*nBe/(nB*nB) ) ;
  Double_t b3E =  cB * nC * rBkg * TMath::Sqrt( cBe*cBe/(cB*cB) + rBkge*rBkge/(rBkg*rBkg) + nCe*nCe/(nC*nC) ) ;
  Double_t bE  =  TMath::Sqrt( nDe*nDe + b1E*b1E + b2E*b2E + b3E*b3E );

  Double_t c1E =  nC * nB * rBkg * TMath::Sqrt( nCe*nCe/(nC*nC) + rBkge*rBkge/(rBkg*rBkg) + nBe*nBe/(nB*nB) ) ;
  Double_t c2E =  nD * nA        * TMath::Sqrt( nDe*nDe/(nD*nD) + nAe*nAe/(nA*nA) );
  Double_t cE  =  TMath::Sqrt(  c1E*c1E + c2E*c2E );
  
  //printf("\t \t a %e (%e) b %e (%e) c %e (%e)\n", a,aE, b, bE, c, cE);

  if ( a == 0  || b*b - 4*a*c < 0 )
  {
    if ( nA > 0 && nD > 0 && nB > 0 )
    {
      //printf("*** Second order term null or negative sqrt !!!! back to no leak\n");
      purity = 1 - rBkg*(nC/nA)/(nD/nB) ;
    }
    else
    {
      //printf("*** Null nA,B,D,C = {%e, %e, %e, %e}!!!\n",nA,nB,nC,nD);
      purity = 0 ;
    }
    
    purityE = purity*TMath::Sqrt( rBkge*rBkge/(rBkg*rBkg) + nAe*nAe / (nA*nA) + nBe*nBe / (nB*nB) + nCe*nCe / (nC*nC) + nDe*nDe / (nD*nD) );
    return;
  }
  
  Float_t xPlus  = (-b+TMath::Sqrt(b*b - 4*a*c))/(2*a);
  Float_t xMinus = (-b-TMath::Sqrt(b*b - 4*a*c))/(2*a);
  
  Double_t acE = a * c * TMath::Sqrt(aE*aE/(a*a)+cE*cE/(c*c));
  Double_t bbE = b*bE ;
  Double_t bbMinus4acE = TMath::Sqrt( bbE*bbE + acE*acE );
  
  //printf("\t \t ac %e (%e) bb %e (%e), bbMinus4acE %e (%e) \n",4 * a * c, acE,b*b, bbE, b*b - 4*a*c, bbMinus4acE);

  Double_t sqrtE = TMath::Sqrt(b*b - 4*a*c) * 1./2. * bbMinus4acE /  (b*b - 4*a*c);
  Double_t numE = TMath::Sqrt( bE*bE + sqrtE*sqrtE );
  Double_t xPlusE  = xPlus  * TMath::Sqrt(numE*numE /( (-b+TMath::Sqrt(b*b - 4*a*c)) * (-b+TMath::Sqrt(b*b - 4*a*c)) ) + aE*aE / (a*a));
  Double_t xMinusE = xMinus * TMath::Sqrt(numE*numE /( (-b-TMath::Sqrt(b*b - 4*a*c)) * (-b-TMath::Sqrt(b*b - 4*a*c)) ) + aE*aE / (a*a));

  //printf("\t \t nAsig (+) %e (%e) (-) %e (%e)\n", xPlus, xPlusE, xMinus, xMinusE);
  purity = xPlus / nA ;

  //purityE = purity * TMath::Sqrt(xPlusE*xPlusE/(xPlus*xPlus) + nAe*nAe / (nA*nA));
  purityE = purity * 0.01;

  
//  Float_t x = (1 - (nC/nA)/(nD/nB)) * nA ;
//  //purity = 1 - (nC/nA)/(nD/nB) ;
//  //purity = 1 - nB*(nC/nA)/nD ;
//  Float_t pCAD = (nC/nA) / nD;
//  //purity = 1 - nB*pCAD ;
//
//  //cB=0;
//  //rBkg = 1;
//  //purity = ((nA - nB*pCAD) / (1-cB*pCAD)) / nA;
//  Float_t newX = (nA*nD - nB*nC) / (nD - cB * nC);
//  purity = (nA*nD - nB * nC * rBkg) / (nD - cB * nC * rBkg) / nA;
//  printf("\t \t x %e new x %e\n",x, newX);
}

/// Calculate the MC alpha factor used in ABCD methods
///
TH1F* GetAlphaFactor(TH1F* hA, TH1F* hB, TH1F* hC, TH1F* hD, TString tag)
{
  if ( !hA || !hB || !hC || !hD ) return 0x0;
  
  TH1F* hACRatio = (TH1F*) hA->Clone();
  TH1F* hBDRatio = (TH1F*) hB->Clone();
  
  hACRatio->Divide(hACRatio,hC,1,1,"");
  hBDRatio->Divide(hBDRatio,hD,1,1,"");
  
  TH1F* hAlpha = (TH1F*) hACRatio->Clone(Form("Alpha_%s",tag.Data()));
  hAlpha->Divide(hAlpha,hBDRatio,1,1,"");
  
  hACRatio->Delete();
  hBDRatio->Delete();
//
//  delete hACRatio;
//  delete hBDRatio;
  
  return hAlpha;
}

///  Main method. Input parameters: the period, and trigger index,
///  the signal to background scale,  the type of signal shower shape selection (0 fix max limit, 1 pT dependent limit)
///  the background isolation momentum and shower shape ranges, (provided as array  indeces since set in different arrays in the code),
///  the type of the MC input file index ( in case several provided, here for different cross talk), the cone value index, in case several provided
///  the A region selection ranges, an string for analysis done with a different global cut (here time cut)
///  and general bools to save the output in a file or do some plotting or printing.
void  Purity3Methods(Int_t iper = 0, Int_t itrig = 1,
             Float_t scaleGJ     = 1.0,
             TString triggerCluster = "MuonCalo",
             Int_t  shPtDep     = 0,

             Int_t   isoBkgMinBin =  4,
             Int_t   isoBkgMaxBin = 16,
             Int_t   m02BkgBinMin = 2,
             Int_t   m02BkgBinMax = 12,
             Int_t   ixTalk       = 1,
             Int_t   icone        = 1,
             
             Float_t isoSigMax   = 1.5, // 1.5, 8.
             Float_t ssSigMin    = 0.1,
             Float_t ssSigMax    = 0.3,
             TString  timeCut    = "",//"_Time10",
             Bool_t  kSave = 1, Bool_t kPlot = 0, Bool_t kPrint = 0)
{
  gStyle->SetOptTitle(1);
  gStyle->SetTitleOffset(2.0,"Y");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(000000);
  gStyle->SetPadRightMargin(0.03);
  //gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  
  // Isolation R values
  //
  Float_t cone[] = {0.4, 0.2};
  
  // Shower shape background limits
  //
  Float_t ssBkgArray[] = {0.30, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0};
  
  Float_t ssBkgMin = ssBkgArray[m02BkgBinMin];
  Float_t ssBkgMax = ssBkgArray[m02BkgBinMax];
  
  if ( ssSigMax > ssBkgMin )
  {
    printf("Skip ssSigMax %f < ssBkgMin %f\n",ssSigMax,ssBkgMin);
    return;
  }
  
  // Isolation momentum background limits
  //
  Float_t isoSigMin =-70.0;
  Float_t isoBkgArray  [] = {
     1.5,  2.0,  2.5,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0,
    12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0,
    55.0, 65.0, 75.0};
  
  Float_t isoBkgMin = isoBkgArray[isoBkgMinBin];
  Float_t isoBkgMax = isoBkgArray[isoBkgMaxBin];

  if ( isoSigMax > isoBkgMin || isoBkgMin >  isoBkgMax )
  {
    printf("Skip isoSigMax %1.2f, isoBkgMin %1.2f, isoBkgMax %1.2f\n",isoSigMax,isoBkgMin,isoBkgMax);
    return;
  }
  
  // Additional options, set the strings for output/input file names
  //
  TString strShShPtDep = "";
  if(shPtDep==1) strShShPtDep += "PtDep";

  if(timeCut!="") strShShPtDep+=timeCut;
  
  TString xTalk[] = {"","_xTalkLess","_xTalkOff"};
  strShShPtDep += xTalk[ixTalk];
 
  TString tagName = "";
    
  tagName+=Form("GJx%2.1f",scaleGJ);
  
  TString outTag = Form("M02%s_S%0.2f-%0.2f_IsoBkg_%1.1f_%1.1f/%s",
                        strShShPtDep.Data() , ssSigMin, ssSigMax,
                        isoBkgMin, isoBkgMax, tagName.Data());
  printf("%s\n",outTag.Data());

  TString outFileName = Form("%s/OutputHisto/%s_%s_Cluster%s_R%1.1f_Iso%1.1f_M02Bkg%1.2f_%1.2f.root",
                             outTag.Data(),period[iper].Data(),trigger[itrig].Data(), triggerCluster.Data(),
                             cone[icone], isoSigMax,ssBkgMin,ssBkgMax);


  // ---- Define shower shape rebin ranges ----
  
  TCustomBinning ssBinning;
  ssBinning.SetMinimum(0.1);
  ssBinning.AddStep(0.22,0.02);
  ssBinning.AddStep(0.30,0.01);
  ssBinning.AddStep(0.50,0.02);
  ssBinning.AddStep(1.00,0.05);
  ssBinning.AddStep(2.00,0.10);
  
  TArrayD binsArray;
  ssBinning.CreateBinEdges(binsArray);

  TCustomBinning ssBinning2;
  ssBinning2.SetMinimum(0.1);
  ssBinning2.AddStep(0.22,0.02);
  ssBinning2.AddStep(0.30,0.01);
  ssBinning2.AddStep(0.50,0.05);
  ssBinning2.AddStep(1.00,0.10);
  ssBinning2.AddStep(2.00,0.20);
  
  TArrayD binsArray2;
  ssBinning2.CreateBinEdges(binsArray2);
    
  // ---- Define/Init Histograms ------

  TH1F *hSig  [nBins];
  TH1F *hSigB [nBins];
  TH1F *hBkg  [nBins];
  TH1F *hBkgW [nBins];
  TH1F *hJJsig[nBins];
  TH1F *hJJbkg[nBins];
  TH1F *hData [nBins];
  TH1F *hRat  [nBins];
  TH1F *hTot  [nBins];
  
  Double_t purityPt   [nBins] = {0.};
  Double_t purityPtErr[nBins] = {0.};
  Double_t chi2       [nBins] = {0.};


  TString processline = "";
  if ( kSave )
  {
    processline = Form(".! mkdir -pv %s/OutputHisto",outTag.Data()) ;
    gROOT->ProcessLine(processline.Data());
  }
  
  if ( kPlot )
  {
    processline = Form(".! mkdir -pv %s/figures",outTag.Data()) ;
    gROOT->ProcessLine(processline.Data());
  }
  
  TH1F* hChi2 ;
  TH1F* hPurityTemp ;
  TH1F* hPurityABCD ;
  TH1F* hPurityAtlas;
  
  // Get MC files, usually one signal MC GJ and 2 bkg jet-jet MC JJlow and JJhigh
  //
  TString mcFileName = "AnalysisResults";
  TFile * fgj  = new TFile(Form("MC/GJ%s/%s/%s.root"    , xTalk[ixTalk].Data(), period[iper].Data(), mcFileName.Data()));
  TFile * fjjl = new TFile(Form("MC/JJlow%s/%s/%s.root" , xTalk[ixTalk].Data(), period[iper].Data(), mcFileName.Data()));
  TFile * fjjh = new TFile(Form("MC/JJhigh%s/%s/%s.root", xTalk[ixTalk].Data(), period[iper].Data(), mcFileName.Data()));
  
  if ( !fgj || !fjjl || !fjjh )
  {
    printf("One MC file not found: gj %p, jjlow %p, jjhigh %p\n",fgj, fjjl, fjjh);
    return;
  }
  
  TFile *fDataL1 = new TFile(Form("../Data/%s_%s_L1%s.root", period[iper].Data(), triggerCluster.Data(),timeCut.Data()));
  TFile *fDataMB = new TFile(Form("../Data/%s_%s_MB%s.root", period[iper].Data(), triggerCluster.Data(),timeCut.Data()));

  if ( !fDataL1 || !fDataMB )
  {
    printf("One Data file not found: L1 %p, MB %p\n", fDataL1, fDataMB);
    return;
  }
  
  // Get  TH3 histogram (pT cluster, sigma, pT iso)  in data and MC
  // Inside analysis output file one per R value
  //
  TString histoName = Form("AnaIsolPhoton_R%1.2f_hPtM02SumPtCone",cone[icone]);
  TH3F * h3GJ   = (TH3F *) fgj  ->Get(histoName);
  if(!h3GJ)
  {
    printf("GJ histo not found!\n");
    return;
  }
    
  // In case of signal/bkg scale variation, scale the signal histogram here
  //
  if ( scaleGJ > 0 ) h3GJ->Scale(scaleGJ);

  TH3F * h3JJl = (TH3F *) fjjl   ->Get(histoName);
  if(!h3JJl)
  {
    printf("h3JJl histo not found!\n");
    return;
  }
  
  TH3F * h3JJh  = 0;
  if ( fjjh ) h3JJh = (TH3F *) fjjh->Get(histoName);
  
  //
  // Remove photon fragmentation from MC, if not suppressed already at the analysis level
  //
  TString histoNameFrag = Form("AnaIsolPhoton_R%1.2f_hPtM02SumPtCone_MCPhotonFrag",cone[icone]);
  h3JJl->Add((TH3F *) fjjl->Get(histoNameFrag),-1);
  if ( fjjh ) h3JJh->Add((TH3F *) fjjh ->Get(histoNameFrag),-1);

  TH3F * h3DataL1 = (TH3F *) fDataL1->Get(histoName);
  TH3F * h3DataMB = (TH3F *) fDataMB->Get(histoName);
  h3DataL1->Sumw2();
  h3DataMB->Sumw2();
    
  // Prepare titles for histograms and temporal strings
  //
  TString tag = Form("Iso%1.1f_Cone%1.1f", isoSigMax, cone[icone]);
  
  TString title1 = Form("%s, %s, cluster-%s, #it{R} = %1.1f", period[iper].Data(),trigger[itrig].Data(), triggerCluster.Data(), cone[icone]);
  TString title2 = Form("Sig: #it{p}_{T}^{iso} < %1.1f GeV/#it{c} - Bkg: %1.1f<#it{p}_{T}^{iso} < %1.1f GeV/#it{c}",
                        isoSigMax, isoBkgMin, isoBkgMax);

  printf("%s - S iso<%1.1f  %1.1f<M02<%1.1f - B %1.1f<iso<%1.1f - %1.1f<M02<%1.1f - GJx%1.1f\n", title1.Data(),
         isoSigMax, ssSigMin, ssSigMax, isoBkgMin, isoBkgMax, ssBkgMin, ssBkgMax, scaleGJ);
  
  //----------------------------
  // Shower shape TEMPLATE FIT
  //----------------------------
  
  Float_t chi2Cut   = 500;//cut to accept fit results
  Bool_t  bApplyBkgWeight = 1;
  Int_t nBinsFit = nBins;
  if(trigger[itrig] == "MB") nBinsFit = 15;
  
  // Loop on each cluster pT bin and do the calculation
  //
  for(int i = 0; i < nBinsFit; i++)
  {
    hSig  [i] = 0x0;
    hSigB [i] = 0x0;
    hBkg  [i] = 0x0;
    hBkgW [i] = 0x0;
    hJJsig[i] = 0x0;
    hJJbkg[i] = 0x0;
    hData [i] = 0x0;
    hRat  [i] = 0x0;
    hTot  [i] = 0x0;
    
    //printf("Bin %d, %1.f<pT<%1.f\n ",i, binPt[i],binPt[i+1]);
    h3GJ    ->SetAxisRange(binPt[i],binPt[i+1]-0.001,"X");
    h3DataL1->SetAxisRange(binPt[i],binPt[i+1]-0.001,"X");
    h3DataMB->SetAxisRange(binPt[i],binPt[i+1]-0.001,"X");
    h3JJl   ->SetAxisRange(binPt[i],binPt[i+1]-0.001,"X");
    if ( h3JJh ) h3JJh->SetAxisRange(binPt[i],binPt[i+1]-0.001,"X");

    // Isolated signal and data
    h3GJ    ->SetAxisRange(isoSigMin,isoSigMax-0.001,"Z");
    h3DataL1->SetAxisRange(isoSigMin,isoSigMax-0.001,"Z");
    h3DataMB->SetAxisRange(isoSigMin,isoSigMax-0.001,"Z");
    
    hSig[i] = (TH1F*) h3GJ->Project3D("y");
    hSig[i]->SetName(Form("hMCSignal_Ebin%d_%s",i,tag.Data()));
    
    TH1F* hDataMB = (TH1F*) h3DataMB->Project3D("y");
    hDataMB->SetName(Form("hIsoDataMB_Ebin%d_%s",i,tag.Data()));
    
    TH1F* hDataL1 = (TH1F*) h3DataL1->Project3D("y");
    hDataL1->SetName(Form("hIsoDataL1_Ebin%d_%s",i,tag.Data()));
    
    // Not isolated signal and data
    h3GJ    ->SetAxisRange(isoBkgMin,isoBkgMax-0.001,"Z");
    h3DataL1->SetAxisRange(isoBkgMin,isoBkgMax-0.001,"Z");
    h3DataMB->SetAxisRange(isoBkgMin,isoBkgMax-0.001,"Z");
    
    hSigB[i] = (TH1F*) h3GJ->Project3D("y");
    hSigB[i]->SetName(Form("hMCSignalNoIso_Ebin%d_%s",i,tag.Data()));
    
    TH1F* hBkgMB = (TH1F*) h3DataMB->Project3D("y");
    hBkgMB->SetName(Form("hNoIsoDataMB_Ebin%d_%s",i,tag.Data()));
    
    TH1F* hBkgL1 = (TH1F*) h3DataL1->Project3D("y");
    hBkgL1->SetName(Form("hNoIsoDataL1_Ebin%d_%s",i,tag.Data()));
    
    if(trigger[itrig] == "L1")
    {
      hData[i] = hDataL1;
      hBkg [i] = hBkgL1;
    }
    else
    {
      hData[i] = hDataMB;
      hBkg [i] = hBkgMB ;
      
      // Add here MB and L1
      if ( trigger[itrig] == "L1MB" )
      {
        hData[i]->Add(hDataL1);
        hBkg [i]->Add(hBkgL1);
      }
    }
    
    hJJsig[i] = 0;
    hJJbkg[i] = 0;
    // Extract the isolated and not isolated bkg shower shape distributions
    // A different MC depending the cluster pT
    //
    if ( bApplyBkgWeight )
    {
     if ( h3JJh && binPt[i] > 16)
      {
        h3JJh->SetAxisRange(isoSigMin,isoSigMax-0.001,"Z");
        hJJsig[i] = (TH1F*) h3JJh->Project3D("y");
        hJJsig[i]->SetName(Form("hMCBkgIsoSig_Ebin%d_%s",i,tag.Data()));
        
        h3JJh->SetAxisRange(isoBkgMin,isoBkgMax-0.001,"Z");
        hJJbkg[i] = (TH1F*) h3JJh->Project3D("y");
        hJJbkg[i]->SetName(Form("hMCBkgIsoBkg_Ebin%d_%s",i,tag.Data()));
      }
      else
      {
        h3JJl->SetAxisRange(isoSigMin,isoSigMax-0.001,"Z");
        hJJsig[i] = (TH1F*) h3JJl->Project3D("y");
        hJJsig[i]->SetName(Form("hMCBkgIsoSig_Ebin%d_%s",i,tag.Data()));
        
        h3JJl->SetAxisRange(isoBkgMin,isoBkgMax-0.001,"Z");
        hJJbkg[i] = (TH1F*) h3JJl->Project3D("y");
        hJJbkg[i]->SetName(Form("hMCBkgIsoBkg_Ebin%d_%s",i,tag.Data()));
      }
    }
    
    //
    // REBIN HISTOGRAMS
    //
    //printf("Before Rebin, pT %1.0f \n",binPt[i]);
    if(binPt[i]<50)
    {
      hSig [i] = (TH1F*) hSig [i]->Rebin(binsArray.GetSize()-1,Form("Sig__%s_r",hSig [i]->GetName()),binsArray.GetArray());
      hSigB[i] = (TH1F*) hSigB[i]->Rebin(binsArray.GetSize()-1,Form("SigB_%s_r",hSigB[i]->GetName()),binsArray.GetArray());
      hBkg [i] = (TH1F*) hBkg [i]->Rebin(binsArray.GetSize()-1,Form("Bkg__%s_r",hBkg [i]->GetName()),binsArray.GetArray());
      hData[i] = (TH1F*) hData[i]->Rebin(binsArray.GetSize()-1,Form("Data_%s_r",hData[i]->GetName()),binsArray.GetArray());
      if(hJJbkg[i]) hJJbkg [i] = (TH1F*) hJJbkg [i]->Rebin(binsArray.GetSize()-1,Form("%sJJbkg_r",hJJbkg[i]->GetName()),binsArray.GetArray());
      if(hJJsig[i]) hJJsig [i] = (TH1F*) hJJsig [i]->Rebin(binsArray.GetSize()-1,Form("%sJJsig_r",hJJsig[i]->GetName()),binsArray.GetArray());
    }
    else
    {
      //printf("New rb\n");
      hSig [i] = (TH1F*) hSig [i]->Rebin(binsArray2.GetSize()-1,Form("Sig__%s_r",hSig [i]->GetName()),binsArray2.GetArray());
      hSigB[i] = (TH1F*) hSigB[i]->Rebin(binsArray2.GetSize()-1,Form("SigB_%s_r",hSigB[i]->GetName()),binsArray2.GetArray());
      hBkg [i] = (TH1F*) hBkg [i]->Rebin(binsArray2.GetSize()-1,Form("Bkg__%s_r",hBkg [i]->GetName()),binsArray2.GetArray());
      hData[i] = (TH1F*) hData[i]->Rebin(binsArray2.GetSize()-1,Form("Data_%s_r",hData[i]->GetName()),binsArray2.GetArray());
      if(hJJbkg[i]) hJJbkg [i] = (TH1F*) hJJbkg [i]->Rebin(binsArray2.GetSize()-1,Form("%sJJbkg_r",hJJbkg[i]->GetName()),binsArray2.GetArray());
      if(hJJsig[i]) hJJsig [i] = (TH1F*) hJJsig [i]->Rebin(binsArray2.GetSize()-1,Form("%sJJsig_r",hJJsig[i]->GetName()),binsArray2.GetArray());
    }
    
    //printf("After Rebin \n");
    ScaleBinBySize(hSig [i]);
    ScaleBinBySize(hBkg [i]);
    ScaleBinBySize(hData[i]);
    
    ScaleBinBySize(hSigB[i]);
    ScaleBinBySize(hJJsig[i]);
    ScaleBinBySize(hJJbkg[i]);
    
    Int_t rb = 1;
    //if ( binPt[i] >= 60 ) rb = 2;
    if(rb>1)
    {
      hSig [i]->Rebin(rb);
      hBkg [i]->Rebin(rb);
      hData[i]->Rebin(rb);
      hSigB [i]->Rebin(rb);
      if(hJJbkg[i]) hJJbkg[i]->Rebin(rb);
      if(hJJsig[i]) hJJsig[i]->Rebin(rb);
    }
    //
    // END REBINS
    //
    
    //
    // SCALE TO DATA
    //
    Float_t fitMin = hData[i]->GetXaxis()->FindBin(0.10);
    Float_t fitMax = hData[i]->GetXaxis()->FindBin(1.99);

    Double_t ratio_sig = 0.;
    Double_t ratio_bkg = 0.;
    Double_t nDat =  hData[i]->Integral(fitMin,fitMax);
    if ( nDat > 200) // with 100 many fit errors in last bins
    {
      ratio_sig = hSig[i]->Integral(fitMin,fitMax) / nDat;
      ratio_bkg = hBkg[i]->Integral(fitMin,fitMax) / nDat;
    }
    else
    {
      purityPt   [i] = 0.0; //to avoid NaN problem when subtracting pi0 contamination
      purityPtErr[i] = 0.;
      continue;
    }
    
    // Comment if Alwina's way
    if(ratio_bkg <= 0 || ratio_sig <= 0)
    {
      purityPt   [i] = 0.0; //to avoid NaN problem when subtracting pi0 contamination
      purityPtErr[i] = 0.;
      
      continue;
    }
    
    hSig [i]->Scale(1./ratio_sig);
    hSigB[i]->Scale(1./ratio_sig);
    hBkg [i]->Scale(1./ratio_bkg);
    
    // ----------------------------------------------------
    // Apply weight to bkg based on MC
    // If not like Alwina: Do it in 2 iterations first without correction and
    // second removing a first/second/etc estimation of the signal in the background regions (like ABCD ATLAS)
    // This first estimation stored in output file before.
    // ----------------------------------------------------
    if ( bApplyBkgWeight )
    {
      // Like Alwina:
      //
      // ADD iso GJ to JJ bkg??
//      if ( scaleGJ > 0.1 )
//      {
//        printf("Scale bkg template");
//        hJJbkg[i]->Add(hSig[i],scaleGJ);
//        //hJJsig[i]->Add(hSig[i],scaleGJ);
//        //hJJbkg[i]->Add(hSigB[i],scaleGJ);
//      }
      
      // Comment if not like Alwina
      //
      if ( scaleGJ > 0.01 )
      {
        // Iteration 0
        TString tagInName = "GJx0.0";
        // Iteration > 0
        //TString tagInName = Form("GJx%2.1f_",scaleGJ);

        TString inTag = Form("M02%s_S%0.2f-%0.2f_IsoBkg_%1.1f_%1.1f/%s",
                             strShShPtDep.Data() , ssSigMin, ssSigMax,
                             isoBkgMin, isoBkgMax, tagInName.Data());

        TFile * inFile = TFile::Open(Form("%s/OutputHisto/%s_%s_Cluster%s_R%1.1f_Iso%1.1f_M02Bkg%1.2f_%1.2f.root",
                                          inTag.Data(),period[iper].Data(),trigger[itrig].Data(), triggerCluster.Data(),
                                          cone[icone], isoSigMax,ssBkgMin,ssBkgMax));
        if(!inFile)
        {
          printf("Make sure you have done scale 0!!!!\n");
          return;
        }
        
        TH1F* hSigBscaled = (TH1F*) inFile->Get(Form("GJMCSigBkg_PtBin%d",i));
        
        if(hSigBscaled)
        {
          //hSigBscaled->Scale(scaleGJ);
          //hBkg[i]->Add(hSigBscaled,-1);
          for(Int_t ix = 1; ix < hSigBscaled->GetNbinsX(); ix++)
          {
            // Avoid fitting problems
            if ( hBkg[i]->GetBinContent(ix) > 0 && hBkg[i]->GetBinContent(ix) > hSigBscaled->GetBinContent(ix) )
              hBkg[i]->SetBinContent(ix,hBkg[i]->GetBinContent(ix)-hSigBscaled->GetBinContent(ix));
          }
        }
        else printf("\t Previous iteration histo bin %d not found!\n",i);
        
        inFile->Close();
        //inFile->Delete();
      }
      
      hBkgW[i] = (TH1F*) hJJsig[i]->Clone(Form("hBkgW_Ebin%d_%s",i,tag.Data()));
      hBkgW[i]->Divide(hJJbkg[i]);
      //hBkgW[i]->SetYTitle(Form("JJ-iso #times #it{R}_{AA} / (JJ-no iso #times #it{R}_{AA} + GJ-iso #times %1.1f)",scaleGJ));
      hBkgW[i]->SetYTitle(Form("JJ-iso #times #it{R}_{AA} / JJ-no iso #times #it{R}_{AA}"));
      
      // Rescale to data
      hBkg[i]->Multiply(hBkgW[i]);
      
      //Coment if alwina
      //hSig [i]->Scale(1./ratio_sig);
      //hSigB[i]->Scale(1./ratio_sig);
      
      ratio_bkg = hBkg[i]->Integral(fitMin,fitMax) / nDat;

      hBkg[i]->Scale(1./ratio_bkg);
    } // Apply weight
      
    // ----------------------------------------------------
    // Fit
    // ----------------------------------------------------

    //printf("Fit\n");

    TObjArray *templates = new TObjArray(2);
    templates->Add(hSig[i]);
    templates->Add(hBkg[i]);
    
    TFractionFitter *fitter = new TFractionFitter(hData[i],templates, "Q");
    fitter->SetRangeX(fitMin,fitMax);
    //fitter->SetRangeX(hData[i]->GetXaxis()->FindBin(0.1),hData[i]->GetXaxis()->FindBin(ssBkgMax-0.001));
    //fitter->ReleaseRangeX();
    fitter->Constrain(0,0.,1.);
    //fitter->Constrain(1,0.,1.);
    
    Int_t status = fitter->Fit();

    Double_t sig_alpha    = 0.;
    Double_t bkg_alpha    = 0.;
    Double_t sig_alphaerr = 0.;
    Double_t bkg_alphaerr = 0.;
    
    fitter->GetResult(0,sig_alpha,sig_alphaerr);
    fitter->GetResult(1,bkg_alpha,bkg_alphaerr);
    //bkg_alpha = 1-sig_alpha;
    chi2[i] = -1;
    
    if ( fitter->GetNDF() > 0  &&  fitter->GetNDF() < 1e8 )
    {
      if ( fitter->GetChisquare() > 0  &&  fitter->GetChisquare() < 1e8 )
      {
          chi2[i] = fitter->GetChisquare()/fitter->GetNDF();
      }
    }
    //hSig[i] = (TH1F*) fitter->GetMCPrediction(0);
    //hBkg[i] = (TH1F*) fitter->GetMCPrediction(1);
    hSig[i]->Scale(sig_alpha);
    hSigB[i]->Scale(sig_alpha);
    hBkg[i]->Scale(bkg_alpha);
    
    hTot[i] = (TH1F*) fitter->GetPlot();
    hRat[i] = (TH1F*) hData[i]->Clone();
    hRat[i]->Divide(hRat[i],hTot[i],1,1,"");
    hRat[i]->SetYTitle("Data / Fit");
    if ( kPrint )
    {
      cout << " sig   = " << sig_alpha << " " << sig_alphaerr ;
      cout << " bkg   = " << bkg_alpha << " " << bkg_alphaerr ;
      cout << " ratio = " << ratio_sig << "|" << ratio_bkg ;
      cout << " chi2  = " << fitter->GetChisquare();
      cout << " ndf   = " << fitter->GetNDF() ;
      cout << " chi2/ndf = " << chi2[i] << endl;
    }
    
    Float_t shSigMax = ssSigMax;
    if ( shPtDep == 1 )
    {
      Float_t pt = binPt[i]+(binPt[i+1]-binPt[i]) / 2.;
      shSigMax = 0.6-0.016*pt;
      if(shSigMax < ssSigMax) shSigMax = ssSigMax;
      //printf("pT %2.2f m02 old %2.2f new %2.2f\n",pt,ssSigMax,shSigMax);
    }
//    if ( shPtDep == 2 )
//    {
//      Float_t pt = binPt[i]+(binPt[i+1]-binPt[i]) / 2.;
//      shSigMax = TMath::Exp(0.19-0.08*pt);
//      if(shSigMax < ssSigMax) shSigMax = ssSigMax;
//      //printf("pT %2.2f m02 old %2.2f new %2.2f\n",pt,ssSigMax,shSigMax);
//    }
    
    Int_t fitSigMax = hData[i]->GetXaxis()->FindBin(shSigMax-0.0001);
    Int_t fitSigMin = hData[i]->GetXaxis()->FindBin(0.1);
    
    //printf("i %d, %2.2f < M02 < %2.2f-%2.2f, bin range %d-%d\n",i,ssSigMin,ssSigMax,shSigMax,fitSigMin,fitSigMax);
    
    //printf("Norm %f\n",hBkg[i]->Integral(fitSigMin,fitSigMax));
    Float_t den = hBkg[i]->Integral(fitSigMin,fitSigMax);
    if ( den < 200 )
    {
      purityPt   [i] = 0.0; //to avoid NaN problem when subtracting pi0 contamination
      purityPtErr[i] = 0.;
      continue;
    }
    
    //printf("Den %f\n",den);
    
    Double_t sigBkg = hSig[i]->Integral(fitSigMin,fitSigMax) / den;
    
    //              Double_t sigBkg_err2 = pow(hSig[i]->Integral(fitSigMin,fitSigMax) * sig_alphaerr,2) / pow(hBkg[i]->Integral(fitSigMin,fitSigMax),2) +
    //              pow(hSig[i]->Integral(fitSigMin,fitSigMax) / pow(hBkg[i]->Integral(fitSigMin,fitSigMax),2),2) * pow(hBkg[i]->Integral(fitSigMin,fitSigMax)*bkg_alphaerr,2);
    //              Double_t sigBkg_err = TMath::Sqrt(sigBkg_err2);
    //
    
    Double_t sigBkg_err = GetFractionError(hSig[i]            ->Integral(fitSigMin,fitSigMax),
                                           hBkg[i]            ->Integral(fitSigMin,fitSigMax),
                                           TMath::Sqrt(hSig[i]->Integral(fitSigMin,fitSigMax)),
                                           TMath::Sqrt(hBkg[i]->Integral(fitSigMin,fitSigMax)));
    
    Double_t purity = sigBkg/(1+sigBkg);
    //Double_t purity_err = TMath::Sqrt(pow(sigBkg_err,2)/pow((1+sigBkg),2));
    Double_t purity_err = GetFractionError(sigBkg, 1+sigBkg, sigBkg_err, sigBkg_err);
    
    if ( kPrint )
    {
      cout << "# sig =" << hSig[i]->Integral(fitSigMin,fitSigMax) << endl;
      cout << "# bkg =" << hBkg[i]->Integral(fitSigMin,fitSigMax) << endl;
      cout << "S/B =" << sigBkg << endl;
      cout << "purity =" << purity << endl;
      cout << "purity error =" << purity_err << endl;
    }
    
    if (purity > 0.01 && chi2[i] < chi2Cut)
    {
      purityPt   [i] = purity;
      purityPtErr[i] = purity_err;
    }
    else
    {
      //cout << "to avoid NaN problem when subtracting pi0 contamination"<<endl;
      purityPt   [i] = 0.0; //to avoid NaN problem when subtracting pi0 contamination
      purityPtErr[i] = 0.;
    }
    
  }
  
  hChi2 = new TH1F("Chi2Template_"+tag,title1+" "+title2,nBins,binPt);
  hPurityTemp = new TH1F("PurityTemplate_"+tag,title1+" "+title2,nBins,binPt);
  for(int i=0;i<nBinsFit;i++)
  {
    hPurityTemp->SetBinContent(i+1,purityPt   [i]);
    hPurityTemp->SetBinError  (i+1,purityPtErr[i]);
    hChi2      ->SetBinContent(i+1,chi2       [i]);
    hChi2      ->SetBinError  (i+1, 0.001);

    if(kPrint)
      cout << "i = " << i << " purity = " << purityPt[i] << " error " <<purityPtErr[i]<<  " chi2 = " << chi2[i] << endl;
  }

  //printf("Plot\n");
  if ( kPlot )
  {
    Int_t firstPt = 2;
    Int_t lastPt = 2+15;
    if(trigger[itrig] == "MB") lastPt = 15;

    TCanvas *c = new TCanvas("c"+tag,"c"+tag,4*800,4*600);
    c->Divide(4,4);

    for(int i = firstPt; i < lastPt; i++) {
      c->cd(i+1-firstPt);

      gPad->SetLogx();
      //gPad->SetLogy();
      //printf("%d hData %p Sig %p Bkg %p SigB%p Tot%p\n",i,hData[i],hSig[i],hBkg[i],hSigB[i],hTot[i]);
      
      if ( hData[i] )
      {
        hData[i]->SetTitle(Form("%1.0f<#it{p}_{T}<%1.0f GeV/#it{c}, P = %2.2f #pm %2.2f, #Chi^{2}/ndf %2.2f",
                                binPt[i], binPt[i+1],purityPt[i],purityPtErr[i], chi2[i]));
        
        hData[i]->GetXaxis()->SetMoreLogLabels(kTRUE);
        hData[i]->SetLineColor(1);
        hData[i]->SetMarkerColor(1);
        hData[i]->Draw();
      }
      
      if ( hSig[i] )
      {
        hSig[i]->SetLineColor(4);
        hSig[i]->SetMarkerColor(4);
        hSig[i]->Draw("same");
        hSig[i]->SetTitle(Form("GJ template, %1.0f<#it{p}_{T}<%1.0f GeV/#it{c}", binPt[i], binPt[i+1]));
      }
      
      if ( hSigB[i] )
      {
        hSigB[i]->SetLineColor(6);
        hSigB[i]->SetMarkerColor(6);
        hSigB[i]->Draw("same");
      }
      
      if ( hBkg[i] )
      {
        hBkg[i]->SetLineColor(2);
        hBkg[i]->SetMarkerColor(2);
        hBkg[i]->SetTitle(Form("Bkg data, %1.0f<#it{p}_{T}<%1.0f GeV/#it{c}", binPt[i], binPt[i+1]));
        hBkg[i]->SetFillStyle(3002);
        hBkg[i]->SetFillColor(kAzure-9);
        hBkg[i]->Draw("same");
      }
      
      if (hTot[i])
      {
        hTot[i]->SetLineColor(kGreen-2);
        hTot[i]->SetMarkerColor(kGreen-2);
        hTot[i]->Draw("same");
      }
    }
    c->cd(16);
    TLegend * legSh = new TLegend(0.,0.,1,1);
    legSh->SetTextSize(0.07);
    legSh->SetFillColor(kWhite);
    legSh->SetLineColor(0);
    legSh->SetHeader(title2);
    legSh->AddEntry("",title1,"");
    if(hData[firstPt])legSh->AddEntry(hData[firstPt],"Data","PL");
    if(hSig[firstPt])legSh->AddEntry(hSig[firstPt],"MC signal","PL");
    if(hSigB[firstPt])legSh->AddEntry(hSigB[firstPt],"MC signal, no iso","PL");
    if(hBkg[firstPt])legSh->AddEntry(hBkg[firstPt],"Background x MC weight","PL");
    if(hTot[firstPt])legSh->AddEntry(hTot[firstPt],"Signal+bkg","PL");
    legSh->Draw();

    c->Print(Form("%s/figures/TemplateM02_%s_%s_Cluster%s_R%1.1f_Iso%1.1f.pdf",
                  outTag.Data(), period[iper].Data(), trigger[itrig].Data(), triggerCluster.Data(), cone[icone], isoSigMax));
    

    TCanvas *cRat = new TCanvas("cRat"+tag,"cRat"+tag,4*800,4*600);
    cRat->Divide(4,4);

    for(int i = firstPt; i < lastPt; i++) {
      cRat->cd(i+1-firstPt);

      gPad->SetLogx();
      gPad->SetGridy();
      //gPad->SetLogy();
      //printf("hData %p Sig %p Bkg %p\n",hData[i],hSig[i],hBkg[i]);
      if(!hRat[i]) continue;
//      hRat[i]->SetTitle(Form("%1.0f<#it{p}_{T}<%1.0f GeV/#it{c}, P = %2.2f #pm %2.2f, #Chi^{2}/ndf %2.2f",
//                             binPt[i], binPt[i+1],purityPt[i],purityPtErr[i], chi2[i]));
      
      hRat[i]->SetTitle(Form("%1.0f<#it{p}_{T}<%1.0f GeV/#it{c}, P = %2.2f #pm %2.2f",
                             binPt[i], binPt[i+1],purityPt[i],purityPtErr[i]));


      hRat[i]->GetXaxis()->SetMoreLogLabels(kTRUE);
      hRat[i]->SetLineColor(1);
      hRat[i]->SetMarkerColor(1);
      hRat[i]->SetMaximum(1.6);
      hRat[i]->SetMinimum(0.6);
      hRat[i]->SetAxisRange(0.1,2);
      hRat[i]->Draw();
      hRat[i]->Fit("pol0","Q","",0.1,2);
      TF1* fifu = hRat[i]->GetFunction("pol0");
      if(fifu)
      {
        TLegend * legF = new TLegend(0.2,0.7,0.93,0.8);
        legF->SetTextSize(0.06);
        legF->SetFillColor(kWhite);
        legF->SetLineColor(0);
        
        legF->AddEntry(fifu,Form("fit %1.3f#pm%1.3f Chi2/ndf %2.2f",
                            fifu->GetParameter(0),
                            fifu->GetParError(0),
                            fifu->GetChisquare()/fifu->GetNDF()),"L");
        legF->Draw();
      }
    }

    cRat->Print(Form("%s/figures/TemplateM02_DataOverFit_%s_%s_Cluster%s_R%1.1f_Iso%1.1f.pdf",
                  outTag.Data(), period[iper].Data(), trigger[itrig].Data(), triggerCluster.Data(), cone[icone], isoSigMax));


    TCanvas *cBkg = new TCanvas("cBkg"+tag,"cBkg"+tag,4*800,4*600);
    cBkg->Divide(4,4);

    for(int i = firstPt; i < lastPt; i++) {
      cBkg->cd(i+1-firstPt);

      gPad->SetLogx();
      //gPad->SetLogy();
      gPad->SetGridy();
      if(!hBkgW[i]) continue;
      
      hBkgW[i]->SetAxisRange(0.1,2);
      //hBkgW[i]->SetMaximum(2.5);
      //hBkgW[i]->SetMinimum(0);
      hBkgW[i]->SetTitle(Form("%1.0f < #it{p}_{T}< %1.0f GeV/#it{c}",binPt[i],binPt[i+1]));
      hBkgW[i]->Draw();
    }

    cBkg->Print(Form("%s/figures/TemplateM02_BkgMCWeight_%s_%s_Cluster%s_R%1.1f_Iso%1.1f.pdf",
                  outTag.Data(), period[iper].Data(), trigger[itrig].Data(), triggerCluster.Data(), cone[icone], isoSigMax));

    TCanvas *cChi2 = new TCanvas("cChi2"+tag,"cChi2"+tag,800,600);
    // gStyle->SetOptStat(0000);
    gPad->SetLogx();
    gPad->SetGridy();

    cChi2->SetRightMargin(0.02);
    cChi2->SetLeftMargin(1.2);
    //cChi2->SetTopMargin(0.02);
    hChi2->GetYaxis()->SetTitleOffset(1.1);
    hChi2->GetXaxis()->SetTitleOffset(1.2);
    hChi2->GetXaxis()->SetMoreLogLabels(kTRUE);

    hChi2->SetMarkerStyle(20);
    hChi2->SetMarkerColor(2);
    hChi2->SetLineColor(2);
    hChi2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hChi2->GetYaxis()->SetTitle("#Chi^{2} / ndf");
    hChi2->SetTitle(title1+" - "+title2);
    hChi2->SetMaximum(500);
    hChi2->SetMinimum(0);
    hChi2->Draw("PE");

//    cChi2->SaveAs(Form("%s/figures/TemplateM02_Chi2_%s_%s_R%1.1f_Iso%1.1f_M02Bkg%1.2f_%1.2f.pdf",
//                             outTag.Data(), period[iper].Data(), trigger[itrig].Data(),
//                             cone[icone], isoSigMax,ssBkgMin,ssBkgMax));

    cChi2->SaveAs(Form("%s/figures/TemplateM02_Chi2_%s_%s_Cluster%s_R%1.1f_Iso%1.1f.pdf",
                       outTag.Data(), period[iper].Data(), trigger[itrig].Data(), triggerCluster.Data(),
                             cone[icone], isoSigMax));
    delete cChi2;
  }
  
  //printf("TEMPLATES\n");
  //
  // TEMPLATE END
  //
  
  //
  // ABCD - ALICE-ATLAS
  //
  //printf("ABCD\n");
  // Open pT range, restricted in template
  //
  h3GJ    ->SetAxisRange(binPt[0],binPt[nBins-1]-0.001,"X");
  h3DataL1->SetAxisRange(binPt[0],binPt[nBins-1]-0.001,"X");
  h3DataMB->SetAxisRange(binPt[0],binPt[nBins-1]-0.001,"X");
  h3JJl   ->SetAxisRange(binPt[0],binPt[nBins-1]-0.001,"X");
  if(h3JJh) h3JJh->SetAxisRange(binPt[0],binPt[nBins-1]-0.001,"X");

  // Project ABCD
  //
  ProjectABCD(isoSigMin, isoSigMax, isoBkgMin, isoBkgMax,
               ssSigMin,  ssSigMax,  ssBkgMin,  ssBkgMax, shPtDep,
              cone[icone], tag+"_L1", 0x0, h3DataL1);
  TH1F* hADataL1 = hAA;
  TH1F* hBDataL1 = hBB;
  TH1F* hCDataL1 = hCC;
  TH1F* hDDataL1 = hDD;
  
  ProjectABCD(isoSigMin, isoSigMax, isoBkgMin, isoBkgMax,
               ssSigMin,  ssSigMax,  ssBkgMin,  ssBkgMax, shPtDep,
              cone[icone], tag+"_MB", 0x0, h3DataMB);
  TH1F* hADataMB = hAA;
  TH1F* hBDataMB = hBB;
  TH1F* hCDataMB = hCC;
  TH1F* hDDataMB = hDD;
  
  TH1F* hAData = 0;
  TH1F* hBData = 0;
  TH1F* hCData = 0;
  TH1F* hDData = 0;
  
  if      ( trigger[itrig] == "MB" )
  {
    hAData = hADataMB;
    hBData = hBDataMB;
    hCData = hCDataMB;
    hDData = hDDataMB;
  }
  else if ( trigger[itrig] == "L1" )
  {
    hAData = hADataL1;
    hBData = hBDataL1;
    hCData = hCDataL1;
    hDData = hDDataL1;
  }
  else
  {
    hAData = (TH1F*) hADataL1->Clone();
    hBData = (TH1F*) hBDataL1->Clone();
    hCData = (TH1F*) hCDataL1->Clone();
    hDData = (TH1F*) hDDataL1->Clone();
    
    hAData->Add(hADataMB);
    hBData->Add(hBDataMB);
    hCData->Add(hCDataMB);
    hDData->Add(hDDataMB);
    //printf("L1MB!\n");
  }
  
  ProjectABCD(isoSigMin, isoSigMax, isoBkgMin, isoBkgMax,
               ssSigMin,  ssSigMax,  ssBkgMin,  ssBkgMax, shPtDep,
              cone[icone], tag+"_JJlow", 0x0, h3JJl);
  TH1F* hAJJlow = hAA;
  TH1F* hBJJlow = hBB;
  TH1F* hCJJlow = hCC;
  TH1F* hDJJlow = hDD;
  
  ProjectABCD(isoSigMin, isoSigMax, isoBkgMin, isoBkgMax,
               ssSigMin,  ssSigMax,  ssBkgMin,  ssBkgMax, shPtDep,
              cone[icone], tag+"_JJhigh", 0x0, h3JJh);
  TH1F* hAJJhigh = hAA;
  TH1F* hBJJhigh = hBB;
  TH1F* hCJJhigh = hCC;
  TH1F* hDJJhigh = hDD;
  
  TH1F* hAJJ = (TH1F*) hAJJhigh->Clone();
  TH1F* hBJJ = (TH1F*) hBJJhigh->Clone();
  TH1F* hCJJ = (TH1F*) hCJJhigh->Clone();
  TH1F* hDJJ = (TH1F*) hDJJhigh->Clone();
  
  // Replace the bins of the final histogram depending the trigger or the MC as a function of a cluster pT threshold
  // Here only the MC
  for(Int_t ix = 0; ix < hAJJ->GetNbinsX(); ix++)
  {
    Float_t center = hAJJ->GetBinCenter(ix);
    if ( center > 16) continue;
    
    hAJJ->SetBinContent(ix,hAJJlow->GetBinContent(ix));
    hBJJ->SetBinContent(ix,hBJJlow->GetBinContent(ix));
    hCJJ->SetBinContent(ix,hCJJlow->GetBinContent(ix));
    hDJJ->SetBinContent(ix,hDJJlow->GetBinContent(ix));
    
    hAJJ->SetBinError(ix,hAJJlow->GetBinError(ix));
    hBJJ->SetBinError(ix,hBJJlow->GetBinError(ix));
    hCJJ->SetBinError(ix,hCJJlow->GetBinError(ix));
    hDJJ->SetBinError(ix,hDJJlow->GetBinError(ix));
  }
  
  //printf("Project high GJ\n");

  ProjectABCD(isoSigMin, isoSigMax, isoBkgMin, isoBkgMax,
               ssSigMin,  ssSigMax,  ssBkgMin,  ssBkgMax, shPtDep,
              cone[icone], tag+"_MB", 0x0, h3GJ);
  TH1F* hAGJ = hAA;
  TH1F* hBGJ = hBB;
  TH1F* hCGJ = hCC;
  TH1F* hDGJ = hDD;
  
  //
  // ABCD ATLAS
  //
  //printf("ATLAS\n");

  // Get BKG correction factor
  //
  TH1F* hAlphaATLAS  = GetAlphaFactor(hAJJ    , hBJJ    , hCJJ    , hDJJ    , tag+"_ATLAS"       );
  TH1F* hAlphaATLASl = GetAlphaFactor(hAJJlow , hBJJlow , hCJJlow , hDJJlow , tag+"_ATLAS_JJlow" );
  TH1F* hAlphaATLASh = GetAlphaFactor(hAJJhigh, hBJJhigh, hCJJhigh, hDJJhigh, tag+"_ATLAS_JJhigh");
  
  // Calculate how much signal leaks to the bkg regions
  //
  //printf("ATLAS leak\n");

  TH1F* hcB = (TH1F*) hBGJ->Clone("cB"+tag);
  TH1F* hcC = (TH1F*) hCGJ->Clone("cC"+tag);
  TH1F* hcD = (TH1F*) hDGJ->Clone("cD"+tag);
  
  hcB->Divide(hcB,hAGJ,1,1,"");
  hcC->Divide(hcC,hAGJ,1,1,"");
  hcD->Divide(hcD,hAGJ,1,1,"");
  
  if ( kPlot )
  {
    TCanvas *cGJSigFrac = new TCanvas("cGJSigFrac"+tag,"cGJSigFrac"+tag,800,600);
    // gStyle->SetOptStat(0000);
    gPad->SetLogx();
    gPad->SetGridy();
    
    cGJSigFrac->SetRightMargin(0.02);
    cGJSigFrac->SetLeftMargin(1.2);
    //cAtlasAlpha->SetTopMargin(0.02);
    hcB->GetYaxis()->SetTitleOffset(1.1);
    hcB->GetXaxis()->SetTitleOffset(1.2);
    hcB->GetXaxis()->SetMoreLogLabels(kTRUE);
    
    hcB->SetMarkerStyle(20);
    hcB->SetMarkerColor(2);
    hcB->SetLineColor(2);
    hcB->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hcB->GetYaxis()->SetTitle("GJ N(B,C,D) / N(A)");
    hcB->SetTitle(title1+" - "+title2);
    hcB->SetMaximum(2);
    hcB->SetMinimum(0);
    
    hcC->SetMarkerStyle(20);
    hcC->SetMarkerColor(1);
    hcC->SetLineColor(1);
    hcC->SetTitle(title1+" - "+title2);

    hcD->SetMarkerStyle(24);
    hcD->SetMarkerColor(4);
    hcD->SetLineColor(4);
    hcD->SetTitle(title1+" - "+title2);

    hcB->Draw("PE");
    hcC->Draw("PEsame");
    hcD->Draw("PEsame");
    
    TLegend * legP = new TLegend(0.6,0.6,0.85,0.8);
    legP->SetTextSize(0.04);
    legP->SetFillColor(kWhite);
    legP->SetLineColor(0);
    //legP->SetHeader(title2);
    //legP->AddEntry("",title1,"");
    legP->AddEntry(hcB,"GJ B / A","PL");
    legP->AddEntry(hcC,"GJ C / A","PL");
    legP->AddEntry(hcD,"GJ D / A","PL");
    legP->Draw();
    
    cGJSigFrac->SaveAs(Form("%s/figures/GJ_SignalFractionPerRegion_%s_R%1.1f_Iso%1.1f_M02Bkg%1.2f_%1.2f.pdf",
                    outTag.Data(), period[iper].Data(),
                            cone[icone], isoSigMax,ssBkgMin,ssBkgMax));
    delete cGJSigFrac;
    
    hcB->GetYaxis()->SetTitle("GJ N(B) / N(A)");
    hcC->GetYaxis()->SetTitle("GJ N(C) / N(A)");
    hcD->GetYaxis()->SetTitle("GJ N(D) / N(A)");
    
    TCanvas *cAtlasAlpha = new TCanvas("cAtlasAlpha"+tag,"cAtlasAlpha"+tag,800,600);
    // gStyle->SetOptStat(0000);
    gPad->SetLogx();
    gPad->SetGridy();
    
    cAtlasAlpha->SetRightMargin(0.02);
    cAtlasAlpha->SetLeftMargin(1.2);
    //cAtlasAlpha->SetTopMargin(0.02);
    hAlphaATLAS->GetYaxis()->SetTitleOffset(1.1);
    hAlphaATLAS->GetXaxis()->SetTitleOffset(1.2);
    hAlphaATLAS->GetXaxis()->SetMoreLogLabels(kTRUE);
    
    hAlphaATLAS->SetMarkerStyle(20);
    hAlphaATLAS->SetMarkerColor(2);
    hAlphaATLAS->SetLineColor(2);
    hAlphaATLAS->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hAlphaATLAS->GetYaxis()->SetTitle("GJ N(B,C,D) / N(A)");
    hAlphaATLAS->SetTitle(title1+" - "+title2);
    hAlphaATLAS->SetMaximum(3);
    hAlphaATLAS->SetMinimum(0);
    hAlphaATLAS->SetYTitle("(JJ(A)#times JJ(D)) / (JJ(B)#times JJ(C))");
    hAlphaATLAS->Draw("PE");
    
    hAlphaATLASl->SetMarkerStyle(24);
    hAlphaATLASl->SetMarkerColor(4);
    hAlphaATLASl->SetLineColor(4);
    hAlphaATLASl->Draw("same");

    if(hAlphaATLASh)
    {
      hAlphaATLASh->SetMarkerStyle(24);
      hAlphaATLASh->SetMarkerColor(2);
      hAlphaATLASh->SetLineColor(2);
      hAlphaATLASh->Draw("same");
    }
    
    TLegend * legA = new TLegend(0.6,0.6,0.9,0.89);
    legA->SetTextSize(0.04);
    legA->SetFillColor(kWhite);
    legA->SetLineColor(0);
    //legP->SetHeader(title2);
    //legP->AddEntry("",title1,"");
    legA->AddEntry(hAlphaATLAS,"JJ Combined","PL");
    legA->AddEntry(hAlphaATLASl,"JJ low","PL");
    if(hAlphaATLASh)
      legA->AddEntry(hAlphaATLASh,"JJ high","PL");
    legA->Draw();
    
    cAtlasAlpha->SaveAs(Form("%s/figures/AlphaMC_ATLAS_%s_R%1.1f_Iso%1.1f_M02Bkg%1.2f_%1.2f.pdf",
                             outTag.Data(), period[iper].Data(),
                             cone[icone], isoSigMax,ssBkgMin,ssBkgMax));
    delete cAtlasAlpha;
  }
  
  //printf("ATLAS purity\n");
  
  // Create the purity histogram with the proper binning
  //
  hPurityAtlas = (TH1F*) hAGJ->Clone(Form("PurityAtlas_%s%s",tag.Data(),period[iper].Data()));
  hPurityAtlas->SetYTitle("#it{Purity}_{ATLAS}");
  
  // So far data driven contamination, get purity 1-C
  //
  for(Int_t ix = 0; ix < hPurityAtlas->GetNbinsX(); ix++)
  {
    Double_t nA  = hAData->GetBinContent(ix);
    Double_t nAe = hAData->GetBinError  (ix);
    Double_t nB  = hBData->GetBinContent(ix);
    Double_t nBe = hBData->GetBinError  (ix);
    Double_t nC  = hCData->GetBinContent(ix);
    Double_t nCe = hCData->GetBinError  (ix);
    Double_t nD  = hDData->GetBinContent(ix);
    Double_t nDe = hDData->GetBinError  (ix);
    
    Double_t cB  = hcB->GetBinContent(ix)*scaleGJ;
    Double_t cBe = hcB->GetBinError  (ix);
    Double_t cC  = hcC->GetBinContent(ix)*scaleGJ;
    Double_t cCe = hcC->GetBinError  (ix);
    Double_t cD  = hcD->GetBinContent(ix)*scaleGJ;
    Double_t cDe = hcD->GetBinError  (ix);
    
    Double_t rBkg  = hAlphaATLAS->GetBinContent(ix);
    Double_t rBkgE = hAlphaATLAS->GetBinError  (ix);
    if ( kPrint )
    {
      printf("ATLAS: Bin %d pt %1.2f\n",ix,hAlphaATLAS->GetBinCenter(ix));
    }
    
    //          cB*=0.1;
    //          cD*=0.1;
    //          printf("\t   nA %1.3e (%1.3e) sqrt(nA) %e nB  %1.3e (%1.3e) sqrt(nB) %e nC  %1.3e (%1.3e) sqrt(nC) %e nD %1.3e (%1.3e) sqrt(nD) %e\n",
    //                 nA, nAe, TMath::Sqrt(nA),
    //                 nB, nBe, TMath::Sqrt(nB),
    //                 nC, nCe, TMath::Sqrt(nC),
    //                 nD, nDe, TMath::Sqrt(nD));
    //          printf("\t   nA %1.3e (%1.3e) nB  %1.3e (%1.3e) nC  %1.3e (%1.3e) nD %1.3e (%1.3e) \n",
    //                 nA, nAe,
    //                 nB, nBe,
    //                 nC, nCe,
    //                 nD, nDe);
    if ( kPrint )
    {
      printf("\t rBkg %1.3f (%1.3f) cB  %1.3f (%1.3f) cC %1.3f (%1.3f) cD  %1.3f (%1.3f) \n",rBkg,rBkgE, cB, cBe, cC, cCe, cD, cDe);
    }
    Float_t purity = 0, purityE = 0;
    GetPurityLeakInBin( nA, nAe, nB, nBe, nC, nCe, nD, nDe,
                       cB, cBe, cC, cCe, cD, cDe,
                       rBkg, rBkgE, purity, purityE);
    
    if ( kPrint )
    {
      printf("\t purity %1.3f (%1.3f) \n", purity, purityE);
    }
    
    hPurityAtlas->SetBinContent(ix, purity );
    hPurityAtlas->SetBinError  (ix, purityE);
  }
  
  TH1F* hsbA = (TH1F*) hAGJ->Clone("sbA"+tag);
  TH1F* hsbB = (TH1F*) hBGJ->Clone("sbB"+tag);
  TH1F* hsbC = (TH1F*) hCGJ->Clone("sbC"+tag);
  TH1F* hsbD = (TH1F*) hDGJ->Clone("sbD"+tag);
  
  hsbA->Divide(hsbA,hAJJ,1,1,"");
  hsbB->Divide(hsbB,hBJJ,1,1,"");
  hsbC->Divide(hsbC,hCJJ,1,1,"");
  hsbD->Divide(hsbD,hCJJ,1,1,"");
  
  if ( kPlot )
  {
    TCanvas *cSB = new TCanvas("cSB"+tag,"cSB"+tag,800,600);
    // gStyle->SetOptStat(0000);
    gPad->SetLogy();
    gPad->SetLogx();
    gPad->SetGridy();
    
    cSB->SetRightMargin(0.02);
    cSB->SetLeftMargin(1.2);
    //cAtlasAlpha->SetTopMargin(0.02);
    hsbB->GetYaxis()->SetTitleOffset(1.1);
    hsbB->GetXaxis()->SetTitleOffset(1.2);
    hsbB->GetXaxis()->SetMoreLogLabels(kTRUE);
    
    hsbB->SetMarkerStyle(20);
    hsbB->SetMarkerColor(2);
    hsbB->SetLineColor(2);
    hsbB->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hsbB->GetYaxis()->SetTitle("GJ / JJ");
    hsbB->SetTitle(title1+" - "+title2);
    hsbB->SetMaximum(3);
    hsbB->SetMinimum(1e-3);
    
    hsbC->SetMarkerStyle(20);
    hsbC->SetMarkerColor(1);
    hsbC->SetLineColor(1);
    hsbC->SetTitle(title1+" - "+title2);
    hsbC->GetYaxis()->SetTitle("GJ / JJ");

    hsbD->SetMarkerStyle(24);
    hsbD->SetMarkerColor(4);
    hsbD->SetLineColor(4);
    hsbD->SetTitle(title1+" - "+title2);
    hsbD->GetYaxis()->SetTitle("GJ / JJ");

    hsbA->SetMarkerStyle(24);
    hsbA->SetMarkerColor(6);
    hsbA->SetLineColor(6);
    hsbA->SetTitle(title1+" - "+title2);
    hsbA->GetYaxis()->SetTitle("GJ / JJ");

    hsbB->Draw("PE");
    hsbC->Draw("PEsame");
    hsbD->Draw("PEsame");
    hsbA->Draw("PEsame");

    TLegend * legSB = new TLegend(0.8,0.65,0.9,0.85);
    legSB->SetTextSize(0.04);
    legSB->SetFillColor(kWhite);
    legSB->SetLineColor(0);
    //legSB->SetHeader(title2);
    //legSB->AddEntry("",title1,"");
    legSB->AddEntry(hsbA,"A","P");
    legSB->AddEntry(hsbB,"B","P");
    legSB->AddEntry(hsbC,"C","P");
    legSB->AddEntry(hsbD,"D","P");
    legSB->Draw();
    
    cSB->SaveAs(Form("%s/figures/SigOverBkg_PerRegion_%s_R%1.1f_Iso%1.1f_M02Bkg%1.2f_%1.2f.pdf",
                            outTag.Data(), period[iper].Data(),
                            cone[icone], isoSigMax,ssBkgMin,ssBkgMax));
    delete cSB;
  }
  
  // ********************
  // ABCD ALICE
  // ********************

  if (scaleGJ > 0)
  {
    hBJJ->Add(hBGJ);
    hCJJ->Add(hCGJ);
    hDJJ->Add(hDGJ);
    
    hBJJlow->Add(hBGJ);
    hCJJlow->Add(hCGJ);
    hDJJlow->Add(hDGJ);
  }
  
  // Get BKG correction factor
  //
  TH1F* hAlphaALICE  = GetAlphaFactor(hAJJ    , hBJJ    , hCJJ    , hDJJ    , tag+"_ALICE");
  TH1F* hAlphaALICEl = GetAlphaFactor(hAJJlow , hBJJlow , hCJJlow , hDJJlow , tag+"_ALICE_JJlow");
  TH1F* hAlphaALICEh = GetAlphaFactor(hAJJhigh, hBJJhigh, hCJJhigh, hDJJhigh, tag+"_ALICE_JJhigh");
  
  if ( kPlot )
  {
    TCanvas *cAliceAlpha = new TCanvas("cAliceAlpha"+tag,"cAliceAlpha"+tag,2*800,1*600);
    // gStyle->SetOptStat(0000);
    
    cAliceAlpha->SetRightMargin(0.02);
    cAliceAlpha->SetLeftMargin(1.2);
    //cAliceAlpha->SetTopMargin(0.02);
    
    cAliceAlpha->Divide(2,1);

    cAliceAlpha->cd(1);
    
    gPad->SetLogx();
    gPad->SetGridy();
    
    hAlphaALICE->GetYaxis()->SetTitleOffset(1.1);
    hAlphaALICE->GetXaxis()->SetTitleOffset(1.2);
    hAlphaALICE->GetXaxis()->SetMoreLogLabels(kTRUE);
    
    hAlphaALICE->SetMarkerStyle(20);
    hAlphaALICE->SetMarkerColor(1);
    hAlphaALICE->SetLineColor(1);
    hAlphaALICE->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hAlphaALICE->GetYaxis()->SetTitle("GJ N(B,C,D) / N(A)");
    hAlphaALICE->SetTitle(title1+" - "+title2);
    hAlphaALICE->SetMaximum(3);
    hAlphaALICE->SetMinimum(0);
    hAlphaALICE->SetYTitle(Form("( JJ(A) #times (JJ*#it{R}_{AA}+GJ*%1.1f)(D) ) / ( (JJ*#it{R}_{AA}+GJ*%1.1f)(B) #times (JJ*#it{R}_{AA}+GJ*%1.1f)(C))",
                                scaleGJ,scaleGJ,scaleGJ));
    hAlphaALICE->Draw("PE");
    
    hAlphaALICEl->SetMarkerStyle(24);
    hAlphaALICEl->SetMarkerColor(4);
    hAlphaALICEl->SetLineColor(4);
    hAlphaALICEl->Draw("same");

    if(hAlphaALICEh)
    {
      hAlphaALICEh->SetMarkerStyle(24);
      hAlphaALICEh->SetMarkerColor(2);
      hAlphaALICEh->SetLineColor(2);
      hAlphaALICEh->Draw("same");
    }
    
    TLegend * legA = new TLegend(0.6,0.6,0.9,0.89);
    legA->SetTextSize(0.04);
    legA->SetFillColor(kWhite);
    legA->SetLineColor(0);
    //legP->SetHeader(title2);
    //legP->AddEntry("",title1,"");
    legA->AddEntry(hAlphaALICE,"JJ Combined","PL");
    legA->AddEntry(hAlphaALICEl,"JJ low","PL");
    if(hAlphaALICEh)
      legA->AddEntry(hAlphaALICEh,"JJ high","PL");
    legA->Draw();
    
    cAliceAlpha->cd(2);
    
    gPad->SetLogx();
    gPad->SetGridy();
    
    TH1F * hRatioCToL = (TH1F*) hAlphaALICE->Clone("AlphaALICERatioCToL");
    hRatioCToL->Divide(hRatioCToL, hAlphaALICEl,1,1,"");
    hRatioCToL->SetMaximum(1.5);
    hRatioCToL->SetMinimum(0.5);
    hRatioCToL->SetYTitle("Ratio #alpha_{MC} JJ x / JJ low");
    hRatioCToL->Draw();

    if(hAlphaALICEh)
    {
      TH1F * hRatioHToL = (TH1F*) hAlphaALICEh->Clone("AlphaALICERatioHToL");
      hRatioHToL->Divide(hRatioHToL, hAlphaALICEl,1,1,"");
      hRatioHToL->Draw("same");
    }
    
    cAliceAlpha->SaveAs(Form("%s/figures/AlphaMC_ALICE_%s_R%1.1f_Iso%1.1f_M02Bkg%1.2f_%1.2f.pdf",
                             outTag.Data(), period[iper].Data(),
                             cone[icone], isoSigMax,ssBkgMin,ssBkgMax));
    delete cAliceAlpha;
  }
  
  // Contamination DD
  //
  TH1F* hCARatio = (TH1F*) hCData->Clone();
  TH1F* hDBRatio = (TH1F*) hDData->Clone();
  
  hCARatio->Divide(hCARatio,hAData,1,1,"");
  hDBRatio->Divide(hDBRatio,hBData,1,1,"");
  
  TH1F* hContamDD = (TH1F*) hCARatio->Clone(Form("ContamDD%s",tag.Data()));
  hContamDD->Divide(hContamDD,hDBRatio,1,1,"");
  
  hCARatio->Delete();
  hDBRatio->Delete();
  
  hPurityABCD = (TH1F*) hContamDD->Clone("PurityALICE_"+tag);
  hPurityABCD->Multiply(hAlphaALICE);
  
  for(Int_t ix = 0; ix < hPurityABCD->GetNbinsX(); ix++)
  {
    hPurityABCD->SetBinContent(ix,1.-hPurityABCD->GetBinContent(ix));
    if ( kPrint )
    {
      printf("ALICE: Bin %d pt %1.2f, purity %2.2f\n",ix, hPurityABCD->GetBinCenter(ix),  hPurityABCD->GetBinContent(ix));
    }
  }
  
  //
  // ABCD - ALICE-ATLAS
  //
  
  for(Int_t ix = 1; ix < hPurityTemp->GetNbinsX(); ix++)
  {
    if ( kPrint )
    {
      printf("Bin %d pt %1.2f, purity: ",ix, hPurityTemp->GetBinCenter(ix));
      printf(" Templ: %2.2f (chi2/ndf %1.0f)", hPurityTemp->GetBinContent(ix),hChi2->GetBinContent(ix));
      printf(" ALICE: %2.2f ", hPurityABCD->GetBinContent(ix));
      printf(" ATLAS: %2.2f \n", hPurityAtlas->GetBinContent(ix));
    }
  }
  
  if(kPlot)
  {
    TCanvas *cPurity = new TCanvas("cPurity"+tag,"cPurity"+tag,800,600);
    // gStyle->SetOptStat(0000);
    gPad->SetLogx();
    gPad->SetGridy();
    
    cPurity->SetRightMargin(0.02);
    cPurity->SetLeftMargin(1.2);
    //cPurity->SetTopMargin(0.02);
    hPurityTemp->GetYaxis()->SetTitleOffset(1.1);
    hPurityTemp->GetXaxis()->SetTitleOffset(1.2);
    hPurityTemp->GetXaxis()->SetMoreLogLabels(kTRUE);
    
    hPurityTemp->SetMarkerStyle(20);
    hPurityTemp->SetMarkerColor(2);
    hPurityTemp->SetLineColor(2);
    hPurityTemp->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hPurityTemp->GetYaxis()->SetTitle("Purity");
    hPurityTemp->SetTitle(title1+" - "+title2);
    hPurityTemp->SetMaximum(0.99);
    hPurityTemp->SetMinimum(0);
    
    hPurityABCD->SetMarkerStyle(20);
    hPurityABCD->SetMarkerColor(1);
    hPurityABCD->SetLineColor(1);

    hPurityAtlas->SetMarkerStyle(24);
    hPurityAtlas->SetMarkerColor(4);
    hPurityAtlas->SetLineColor(4);

    hPurityTemp->Draw("PE");
    hPurityABCD->Draw("PEsame");
    hPurityAtlas->Draw("PEsame");
    
    TLegend * legP = new TLegend(0.6,0.12,0.85,0.3);
    legP->SetTextSize(0.04);
    legP->SetFillColor(kWhite);
    legP->SetLineColor(0);
    //legP->SetHeader(title2);
    //legP->AddEntry("",title1,"");
    legP->AddEntry(hPurityAtlas,"ABCD ATLAS","PL");
    legP->AddEntry(hPurityABCD,"ABCD ALICE","PL");
    legP->AddEntry(hPurityTemp,"Template","PL");
    legP->Draw();
    
    cPurity->SaveAs(Form("%s/figures/Purity_%s_%s_Cluster%s_R%1.1f_Iso%1.1f_M02Bkg%1.2f_%1.2f.pdf",
                    outTag.Data(), period[iper].Data(), trigger[itrig].Data(), triggerCluster.Data(),
                         cone[icone], isoSigMax, ssBkgMin, ssBkgMax));
    delete cPurity;
  }
      
  if ( kSave )
  {
    printf("Write %s\n",outFileName.Data());

    TFile * outFile = new TFile(outFileName,"recreate");
    
    hPurityTemp ->Write("Purity_Template");
    hPurityABCD ->Write("Purity_ABCD_ALICE");
    hPurityAtlas->Write("Purity_ABCD_ATLAS");
    hAlphaATLAS ->Write("AlphaMC_ATLAS");
    hAlphaATLASl ->Write("AlphaMC_ATLAS_JJlow");
    if ( hAlphaATLASh )
    {
      hAlphaATLASh ->Write("AlphaMC_ATLAS_JJhigh");
    }
    
    hAlphaALICEl ->Write("AlphaMC_ALICE_JJlow");
    
    if ( hAlphaALICEh )
    {
      hAlphaALICEh ->Write("AlphaMC_ALICE_JJhigh");
    }
    hContamDD   ->Write("ContamDD_ALICE");
    hChi2       ->Write("Chi2_Template");

    hcB         ->Write("cB");
    hcC         ->Write("cC");
    hcD         ->Write("cD");
    
    hsbA         ->Write("sbA");
    hsbB         ->Write("sbB");
    hsbC         ->Write("sbC");
    hsbD         ->Write("sbD");
    
    for(int i = 0; i < nBinsFit; i++)
    {
      if ( chi2[i] <= 0.001 ) continue;

      if ( hData[i] ) hData[i]->Write(Form("IsoData_PtBin%d",i));
      if ( hBkg[i]  ) hBkg [i]->Write(Form("BkgData_PtBin%d",i));
      if ( hBkgW[i] ) hBkgW[i]->Write(Form("BkgDataWeight_PtBin%d",i));
      if ( hSig[i]  ) hSig [i]->Write(Form("GJMCSig_PtBin%d",i));
      if ( hSigB[i] ) hSigB[i]->Write(Form("GJMCSigBkg_PtBin%d",i));
      if ( hTot[i]  ) hTot [i]->Write(Form("Total_PtBin%d",i));
      if ( hRat[i]  ) hRat [i]->Write(Form("DataFitRat_PtBin%d",i));
    }
    
    outFile->Close();
  }
  
  fgj    ->Close();
  fjjl   ->Close();
  fDataL1->Close();
  fDataMB->Close();
}

