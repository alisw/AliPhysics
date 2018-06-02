///
/// \file InvMassFit.C
/// \ingroup CaloTrackCorrMacrosQA
/// \brief Fit invariant mass distributions
///
/// Macro using as input the 2D histograms mass vs pT of AliAnaPi0
/// For a given set of pT bins invariant mass plots are fitted and mass vs pT 
/// and width vs pT and neutral meson spectra plots are obtained.
/// Also pure MC histograms checking the origin of the clusters and generated spectra
/// is used to plot efficiencies and acceptance
///
/// Based on old macros
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TString.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPad.h>
#include <TFile.h>
#include <TLegend.h>
#include <TObject.h>
#include <TDirectoryFile.h>
#include <TGraphErrors.h>
#include <TList.h>
#include <TArrayD.h>
#include <TGaxis.h>

#endif

// Settings
Bool_t  kMix        = kFALSE;              /// use mixed event to constrain combinatorial background
Int_t   kPolN       = 1;                   /// polinomyal type for residual background under the peak
Bool_t  kSumw2      = kTRUE;               /// Apply Root method Sumw2()
Float_t kNPairCut   = 20;                  /// Minimum number of cluster pairs in pi0 or eta window 
TString kHistoStartName = "AnaPi0_";       /// Common starting name in histograms
TString kProdName   = "LHC18c3_NystromOn"; /// Input production directory name where input file is located
TString kFileName   = "AnalysisResults";   /// Name of file with input histograms
TString kPlotFormat = "eps";               /// Automatic plots format
TString kCalorimeter= "EMCAL";             /// Calorimeter, EMCAL, DCAL (PHOS)
TString kParticle   = "Pi0";               /// Particle searched: "Pi0", "Eta"
Bool_t  kTrunMixFunc= kTRUE;               /// Use a truncated function to get the mixed event scale

// Initializations
Int_t   nEvt  = 0;
TFile * fil   = 0;
TFile * fout  = 0;
TList * lis   = 0;
TF1   * fitfun= 0;
TDirectoryFile * direc =0;
// Mixing
TF1 *tPolPi0 = 0;
TF1 *tPolEta = 0;

//                        - SM                       SM-Sector       - Side                 Side
Int_t modColorIndex[]={1 , 2, 2, 3, 3, 4, 4, 7, 7, 6, 6, 2, 3, 4, 7, 6, 2, 2, 3, 3, 4, 4, 6, 6};
Int_t modStyleIndex[]={20,24,25,24,25,24,25,24,25,24,25,21,21,21,21,21,22,26,22,26,22,26,22,26};

///
/// Open the file and the list and the number of analyzed events
/// 
//-----------------------------------------------------------------------------
Bool_t GetFileAndEvents( TString dirName , TString listName )
{
  fil = new TFile(Form("%s/%s.root",kProdName.Data(),kFileName.Data()),"read");
  
  printf("Open Input File: %s/%s.root\n",kProdName.Data(),kFileName.Data());
  
  if ( !fil ) return kFALSE;
  
  direc  = (TDirectoryFile*) fil->Get(dirName);
  
  if ( !direc && dirName != "" ) return kFALSE;
  
  //printf("dir %p, list %s\n",dir,listName.Data());
  
  if(direc)
    lis = (TList*) direc ->Get(listName);
  else
    lis = (TList*) fil->Get(listName);
  
  if ( !lis && listName != "") return kFALSE;
  
  if(!lis)
    nEvt = ((TH1F*) fil->Get("hNEvents"))->GetEntries();
  else
    nEvt = ((TH1F*) lis->FindObject("hNEvents"))->GetEntries();
  
  printf("nEvt = %d\n",nEvt);
  //nEvt = 1;
  
  return kTRUE;
}

///
/// Gaussian plus polinomial order 3 function
///
//-----------------------------------------------------------------------------
Double_t pi0massP3(Double_t *x, Double_t *par)
{
  Double_t gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                      (2*par[2]*par[2]) );
  Double_t back = par[3] + par[4]*x[0] + par[5]*x[0]*x[0] + par[6]*x[0]*x[0]*x[0];
  return gaus+back;
}

///
/// Gaussian plus polinomial order 2 function
///
//-----------------------------------------------------------------------------
Double_t pi0massP2(Double_t *x, Double_t *par)
{
  Double_t gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                      (2*par[2]*par[2]) );
  Double_t back = par[3] + par[4]*x[0] + par[5]*x[0]*x[0];
  return gaus+back;
}

///
/// Gaussian plus polinomial order 1 function
///
//-----------------------------------------------------------------------------
Double_t pi0massP1(Double_t *x, Double_t *par)
{
  Double_t gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                      (2*par[2]*par[2]) );
  Double_t back = par[3] + par[4]*x[0];
  return gaus+back;
}

///
/// Gaussian plus constant function
///
//-----------------------------------------------------------------------------
Double_t pi0massP0(Double_t *x, Double_t *par)
{
  Double_t gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                      (2*par[2]*par[2]) );
  Double_t back = 0;//par[3];
  return gaus+back;
}

///
/// Truncated function in eta region
///
//-----------------------------------------------------------------------------
Double_t truncatedPolEta(Double_t *x, Double_t *par)
{
  if ((x[0] > 0.425 && x[0] < 0.65) ) {
    TF1::RejectPoint();
    return 0;
  }
  
  return par[0] + par[1]*x[0];// + x[0]*x[0]*par[2]  + x[0]*x[0]*x[0]*par[3];//+ x[0]*x[0]*x[0]*x[0]*par[4];
}

///
/// Truncated function in pi0 region
///
//-----------------------------------------------------------------------------
Double_t truncatedPolPi0(Double_t *x, Double_t *par)
{
  if ((x[0] > 0.07 && x[0] < 0.2)) {
    TF1::RejectPoint();
    return 0;
  }
  
  return par[0] + par[1]*x[0];// + x[0]*x[0]*par[2] + x[0]*x[0]*x[0]*par[3] ;//+ x[0]*x[0]*x[0]*x[0]*par[4];
}

///
/// Crystal ball function
///
//-----------------------------------------------------------------------------
Double_t CrystalBall(Double_t *x, Double_t *par)
{
  Double_t N = par[0];
  Double_t width = par[1];
  
  Double_t mean = par[2];
  Double_t sigma = par[3];
  Double_t alpha = par[4];
  Double_t n = par[5];
  
  Double_t A = pow(n/fabs(alpha),n)*exp(-pow(alpha,2)/2);
  Double_t B = n/fabs(alpha) - fabs(alpha);
  
  if ((x[0]-mean)/sigma>-alpha)
    return N*width*TMath::Gaus(x[0],mean,sigma,1);
  else
    return N/(sqrt(TMath::TwoPi())*sigma)*width*A*pow(B-(x[0]-mean)/sigma,-n);
}

///
/// Initialize the fitting function
///
//-----------------------------------------------------------------------------
void SetFitFun()
{
  //Fitting function
  //if(mix) kPolN = 0;
  
  if(kParticle == "Pi0")
  {
    if ( kPolN == 0)
      fitfun = new TF1("fitfun",pi0massP0,0.100,0.250,4);
    else if      (kPolN == 1)
      fitfun = new TF1("fitfun",pi0massP1,0.100,0.250,5);
    else if (kPolN == 2)
      fitfun = new TF1("fitfun",pi0massP2,0.100,0.250,6);
    else if (kPolN == 3)
      fitfun = new TF1("fitfun",pi0massP3,0.100,0.250,7);
    else
    {
      printf("*** <<< Set Crystal Ball!!! >>> ***\n");
      fitfun = new TF1("fitfun",CrystalBall,0.100,0.250,6);
    }
    
    if(kPolN < 4)
    {
      fitfun->SetParLimits(0,  kNPairCut/5,kNPairCut*1.e4);
      fitfun->SetParLimits(1,  0.105,0.185);
      fitfun->SetParLimits(2,  0.001,0.040);
    }
    else
    {
      fitfun->SetParLimits(0,  kNPairCut/5,kNPairCut*1.e6);
      fitfun->SetParLimits(2,  0.105,0.185);
      fitfun->SetParLimits(1,  0.001,0.040);
      fitfun->SetParLimits(3,  0.001,0.040);
      fitfun->SetParLimits(4,  0,10);
      fitfun->SetParLimits(5,  1,1.e+6);
    }
    
  } // Pi0
  else  // Eta
  {
    if ( kPolN == 0)
      fitfun = new TF1("fitfun",pi0massP0,0.400,0.650,4);
    else if      (kPolN == 1)
      fitfun = new TF1("fitfun",pi0massP1,0.400,0.650,5);
    else if (kPolN == 2)
      fitfun = new TF1("fitfun",pi0massP2,0.400,0.650,6);
    else if (kPolN == 3)
      fitfun = new TF1("fitfun",pi0massP3,0.400,0.650,7);
    if(kPolN < 4){
      fitfun->SetParLimits(0,  kNPairCut/10,1.e+6);
      fitfun->SetParLimits(1,  0.20,0.80);
      fitfun->SetParLimits(2,  0.001,0.06);
    }
  }
  
  fitfun->SetLineColor(kRed);
  fitfun->SetLineWidth(2);
  
  fitfun->SetParName(0,"A");
  fitfun->SetParName(1,"m_{0}");
  fitfun->SetParName(2,"#sigma");
  //  fitfun->SetParName(3,"a_{0}");
  //  fitfun->SetParName(4,"a_{1}");
  
  if (kPolN > 1)
  {
    fitfun->SetParName(5,"a_{2}");
    if (kPolN > 2) fitfun->SetParName(6,"a_{3}"); 
  }  
  
  //Mix fit func
  tPolPi0 = new TF1("truncatedPolPi0", truncatedPolPi0, 0.02, 0.25, 2);
  tPolEta = new TF1("truncatedPolEta", truncatedPolEta, 0.2, 1, 2);
}


//-----------------------
/// Efficiency calculation called in ProjectAndFit()
//-----------------------
void Efficiency(Int_t nPt, TArrayD xPt, TArrayD exPt,
                TArrayD mesonPt, TArrayD mesonPtBC, TString hname)
{
  TGraphErrors * gPrim    = (TGraphErrors*) fout->Get("Primary");
  TGraphErrors * gPrimAcc = (TGraphErrors*) fout->Get("PrimaryInAcceptance");
  if ( !gPrim || !gPrimAcc ) return ;
    
  TArrayD  effPt        ;  effPt        .Set(nPt); 
  TArrayD  effBCPt      ;  effBCPt      .Set(nPt); 
  TArrayD  effPtAcc     ;  effPtAcc     .Set(nPt); 
  TArrayD  effBCPtAcc   ;  effBCPtAcc   .Set(nPt);  
  TArrayD  effPtErr     ;  effPtErr     .Set(nPt); 
  TArrayD  effBCPtErr   ;  effBCPtErr   .Set(nPt); 
  TArrayD  effPtAccErr  ;  effPtAccErr  .Set(nPt); 
  TArrayD  effBCPtAccErr;  effBCPtAccErr.Set(nPt);  
  
  for(Int_t ibin = 0; ibin < nPt; ibin++ )
  {
    //printf("Bin %d \n",ibin);
    //printf("- Fit Reco %2.3e, Prim %2.3e, PrimAcc %2.3e\n",
    //       mesonPt[ibin],gPrim->GetX()[ibin],gPrimAcc->GetX()[ibin]);
    
    if ( gPrim->GetX()[ibin] > 0 && mesonPt[ibin] > 0 )
    {
      effPt[ibin] = mesonPt[ibin] / gPrim->GetX()[ibin];
      effPtErr[ibin] = effPt[ibin] * TMath::Sqrt(1./(mesonPt[ibin]       * mesonPt[ibin])   +
                                                 1./(gPrim->GetX()[ibin] * gPrim->GetX()[ibin]));
      //printf("\t EffxAcc %f, err %f\n",effPt[ibin],effPtErr[ibin]);
    }
    else { effPt[ibin] = 0; effPtErr[ibin] = 0; }
    
    if ( gPrimAcc->GetX()[ibin] > 0 && mesonPt[ibin] > 0 )
    {
      effPtAcc[ibin] = mesonPt[ibin] / gPrimAcc->GetX()[ibin];
      effPtAccErr[ibin] = effPtAcc[ibin] * TMath::Sqrt(1./(mesonPt[ibin]          * mesonPt[ibin])   +
                                                       1./(gPrimAcc->GetX()[ibin] * gPrimAcc->GetX()[ibin]));
      //printf("\t Eff %f, err %f\n",effPtAcc[ibin],effPtAccErr[ibin]);
    }
    else { effPtAcc[ibin] = 0; effPtAccErr[ibin] = 0; }
    
    if ( kMix )
    {
      //printf("- BC  Reco %2.3e, Prim %2.3e, PrimAcc %2.3e\n",
      // mesonPtBC[ibin],gPrim->GetX()[ibin],gPrimAcc->GetX()[ibin]);
      if ( gPrim->GetX()[ibin] > 0 && mesonPtBC[ibin] > 0 )
      {
        effBCPt[ibin] = mesonPtBC[ibin] / gPrim->GetX()[ibin];
        effBCPtErr[ibin] = effBCPt[ibin] * TMath::Sqrt(1./(mesonPtBC[ibin]     * mesonPtBC[ibin])   +
                                                       1./(gPrim->GetX()[ibin] * gPrim->GetX()[ibin]));
        //printf("\t EffxAcc %f, err %f\n",effBCPt[ibin],effBCPtErr[ibin]);
      }
      else { effBCPt[ibin] = 0; effBCPtErr[ibin] = 0; }
      
      if ( gPrimAcc->GetX()[ibin] > 0 && mesonPtBC[ibin] > 0 )
      {
        effBCPtAcc[ibin] = mesonPtBC[ibin] / gPrimAcc->GetX()[ibin];
        effBCPtAccErr[ibin] = effBCPtAcc[ibin] * TMath::Sqrt(1./(mesonPtBC[ibin]     * mesonPtBC[ibin])   +
                                                             1./(gPrimAcc->GetX()[ibin] * gPrimAcc->GetX()[ibin]));
        //printf("\t EffxAcc %f, err %f\n",effBCPtAcc[ibin],effBCPtAccErr[ibin]);
      }
      else { effBCPtAcc[ibin] = 0; effBCPtAccErr[ibin] = 0; }
      
    } // Bin counting
    
  } // pt bin loop
  
  TGraphErrors *gEff = new TGraphErrors(nPt,xPt.GetArray(),effPt.GetArray(),exPt.GetArray(),effPtErr.GetArray());
  gEff->SetName(Form("EfficiencyxAcceptance_%s",hname.Data()));
  gEff->GetHistogram()->SetYTitle("#epsilon_{Reco #times PID #times Acc}");
  gEff->Write();
  
  TGraphErrors *gEffAcc = new TGraphErrors(nPt,xPt.GetArray(),effPtAcc.GetArray(),exPt.GetArray(),effPtAccErr.GetArray());
  gEffAcc->SetName(Form("Efficiency_%s",hname.Data()));
  gEffAcc->GetHistogram()->SetYTitle("#epsilon_{Reco #times PID}");
  gEffAcc->Write();    
  
  TGraphErrors *gBCEff    = 0;
  TGraphErrors *gBCEffAcc = 0;
  if(kMix)
  {
    gBCEff = new TGraphErrors(nPt,xPt.GetArray(),effBCPt.GetArray(),exPt.GetArray(),effBCPtErr.GetArray());
    gBCEff->SetName(Form("EfficiencyxAcceptance_BC_%s",hname.Data()));
    gBCEff->GetHistogram()->SetYTitle("#epsilon_{Reco #times PID #times Acc}");
    gBCEff->Write();
    
    gBCEffAcc = new TGraphErrors(nPt,xPt.GetArray(),effBCPtAcc.GetArray(),exPt.GetArray(),effBCPtAccErr.GetArray());
    gBCEffAcc->SetName(Form("Efficiency_BC_%s",hname.Data()));
    gBCEffAcc->GetHistogram()->SetYTitle("#epsilon_{Reco #times PID}");
    gBCEffAcc->Write();   
  }
  
  // Plot efficiencies
  //
  TCanvas *cEff = new TCanvas(Form("cEff_%s",hname.Data()),
                              Form("Efficiency Graphs for %s",hname.Data()),2*600,600);
  cEff->Divide(2, 1);
  
  cEff->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  //gPad->SetLogy();
  
  //gEff->SetMaximum(1);
  //gEff->SetMinimum(1e-8);
  
  gEff->Draw("AP");
  
  if(kMix)
  {
    gBCEff ->Draw("P");
    
    TLegend * legend =  new TLegend(0.5,0.7,0.9,0.9);
    legend->AddEntry(gEff  ,"From fit","P");
    legend->AddEntry(gBCEff,"From bin counting","P");
    legend->Draw();
  }
  
  cEff->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  //gPad->SetLogy();
  
  //gEff->SetMaximum(1);
  //gEff->SetMinimum(1e-8);
  
  gEffAcc->Draw("AP");
  
  if(kMix)
  {
    gBCEffAcc->Draw("P");
    
    TLegend * legend =  new TLegend(0.5,0.7,0.9,0.9);
    legend->AddEntry(gEffAcc  ,"From fit","P");
    legend->AddEntry(gBCEffAcc,"From bin counting","P");
    legend->Draw();
  }
  
  cEff->Print(Form("IMfigures/%s_%s_Efficiency_%s_%s.%s",
                   kProdName.Data(),kCalorimeter.Data(),kParticle.Data(),hname.Data(),
                   /*kFileName.Data(),*/ kPlotFormat.Data()));   
}

//-------------------------------------
/// Do energy projections, fits and plotting
/// of invariant mass distributions.
//-------------------------------------
void ProjectAndFit
(
 const Int_t nPt, TString comment, TString leg, TString hname,   
 TArrayD xPtLimits, TArrayD xPt, TArrayD exPt, 
 TH2F * hRe, TH2F * hMi
 )
{
  if ( !hRe ) 
  {
    printf("ProjectAndFit() - Null inv mass 2D histo for %s, skip!\n",hname.Data());
    return;
  }
  
  //gROOT->Macro("$MACROS/style.C");//Set different root style parameters
  gStyle->SetOptTitle(1);
  gStyle->SetTitleOffset(2,"Y");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(000000);
  
  Double_t mmin = 0;
  Double_t mmax = 1;
  TString  particleN = "";

  if(kParticle=="Pi0")
  {
    mmin = 0.00;
    mmax = 0.45;
    particleN = " #pi^{0}";
  }
  else // eta
  {
    mmin = 0.25;
    mmax = 0.95;
    particleN = " #eta";
  }
  
  TLegend * pLegendIM[nPt];
  TH1D* hIM           [nPt];
  TH1D* hMix          [nPt];
  TH1D* hMixCorrected [nPt];
  TH1D* hRatio        [nPt];
  TH1D* hSignal       [nPt];
  for(Int_t ipt = 0; ipt< nPt; ipt++)
  {
    hIM          [ipt] = 0 ;
    hMix         [ipt] = 0 ;
    hMixCorrected[ipt] = 0 ;
    hRatio       [ipt] = 0 ;
    hSignal      [ipt] = 0 ;
  }

  TArrayD  mesonPt   ;  mesonPt   .Set(nPt); 
  TArrayD  mesonPtBC ;  mesonPtBC .Set(nPt); 
  TArrayD  mesonMass ;  mesonMass .Set(nPt); 
  TArrayD  mesonWidth;  mesonWidth.Set(nPt);  
  TArrayD emesonPt   ; emesonPt   .Set(nPt); 
  TArrayD emesonPtBC ; emesonPtBC .Set(nPt); 
  TArrayD emesonMass ; emesonMass .Set(nPt); 
  TArrayD emesonWidth; emesonWidth.Set(nPt); 
  
  Int_t rebin  = 1 ;
  Int_t rebin2 = 2*rebin;
  
  Int_t col = TMath::Ceil(TMath::Sqrt(nPt));

  //hRe->Scale(1./nEvt);
  
  // Mix correction factor
  TF1 *line3[nPt];// = new TF1("line3", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", 0.0, 1);
  
  TF1 * fitFunction = 0;
  
  TCanvas  * cIMModi = new TCanvas(Form("c%s", hname.Data()), Form("%s", leg.Data()), 1200, 1200) ;
  cIMModi->Divide(col, col);
  
  for(Int_t i = 0; i < nPt; i++)
  {
    //if(icomb >10 && i > 8) continue; 
    
    cIMModi->cd(i+1) ; 
    //gPad->SetLogy();
    
    Double_t ptMin = xPtLimits[i];
    Double_t ptMax = xPtLimits[i+1];
    //printf("Bin %d (%2.1f, %2.1f)\n",i, ptMin,ptMax);
    
    hIM[i] = hRe->ProjectionY(Form("IM_%s_PtBin%d",hname.Data(),i),
                                            hRe->GetXaxis()->FindBin(ptMin),
                                            hRe->GetXaxis()->FindBin(ptMax));
    hIM[i]->SetTitle(Form("%2.1f < p_{T, #gamma#gamma} < %2.1f GeV/c",ptMin,ptMax));
    if(i < nPt-3)
      hIM[i]->Rebin(rebin);
    else
      hIM[i]->Rebin(rebin2);
    
    if ( kSumw2 )
      hIM[i]->Sumw2();
    
    hIM[i]->SetXTitle("M_{#gamma,#gamma} (GeV/c^{2})");
    //hIM[i]->SetLineColor(modColorIndex);
    //printf("mmin %f, mmax %f\n",mmin,mmax);
    hIM[i]->SetAxisRange(mmin,mmax,"X");
    hIM[i]->SetLineWidth(2);
    hIM[i]->SetLineColor(4);
    if ( !kMix ) hIM[i]->SetMinimum( 0.1);
    else         hIM[i]->SetMinimum(-0.1);
    Double_t mStep =  hIM[i]->GetBinWidth(1);
    
    hSignal[i] = (TH1D*) hIM[i]->Clone();
        
    //--------------------------------------------------
    //   Mix
    //--------------------------------------------------
    if ( kMix && hMi )
    {
      hMix[i] = (TH1D*) hMi->ProjectionY(Form("MiMass_PtBin%d_%s",i,hname.Data()),
                                                       hMi->GetXaxis()->FindBin(ptMin), 
                                                       hMi->GetXaxis()->FindBin(ptMax));
      hMix[i]->SetLineColor(1);  
      hMix[i]->SetTitle(Form("%2.2f < p_{T, #gamma#gamma} < %2.2f GeV/c",ptMin,ptMax));
      if(i < nPt-3)
        hMix[i]->Rebin(rebin);
      else
        hMix[i]->Rebin(rebin2);
      
      //         Double_t sptmin = 0.2;
      //         Double_t sptmax = 0.45;
      
      //         Double_t scalereal = hIM->Integral(hIM->FindBin(sptmin),
      //                                                             hIM->FindBin(sptmax));
      //         Double_t scalemix  = hMix->Integral(hMix->FindBin(sptmin),
      //                                                             hMix->FindBin(sptmax));
      //         if(scalemix > 0)
      //           hMix->Scale(scalereal/scalemix);
      //         else{ 
      //           //printf("could not scale:  re %f  mi %f\n",scalereal,scalemix);
      //           continue;
      //         }
      
      //      hMix[i]->SetLineWidth(1);  
      
      if ( kSumw2 ) 
        hMix[i]->Sumw2();
      
      // --- Ratio ---  
      hRatio[i] = (TH1D*)hIM[i]->Clone(Form("RatioRealMix_PtBin%d_%s",i,hname.Data()));
      hRatio[i]->SetAxisRange(mmin,mmax,"X");
      hRatio[i]->Divide(hMix[i]);
      Double_t p0 =0;
      Double_t p1 =0;
      Double_t p2 =0;
      Double_t p3 =0;
      
      // --- Subtract from signal ---
      hSignal[i] ->SetName(Form("Signal_PtBin%d_%s",i,hname.Data()));

      if ( hMix[i]->GetEntries() > 100 )
      {
        hSignal[i] ->SetLineColor(kViolet);

        if ( kTrunMixFunc )
        {
          if ( kParticle == "Pi0" )
          {
            hRatio[i]->Fit("truncatedPolPi0", "NQRL","",0.11,0.4);
            p0=  tPolPi0->GetParameter(0);           
            p1=  tPolPi0->GetParameter(1);
            p2=  tPolPi0->GetParameter(2);
            p3=  tPolPi0->GetParameter(3);
          } // pi0
          else  // eta
          {
            hRatio[i]->SetAxisRange(mmin,mmax);
            hRatio[i]->Fit("truncatedPolEta", "NQRL","",0.25,0.95);
            p0 =  tPolEta->GetParameter(0);      
            p1 =  tPolEta->GetParameter(1);
            p2 =  tPolEta->GetParameter(2);
            p3 =  tPolEta->GetParameter(3);
          }
          
          line3[i] = new TF1("line3", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", mmin, mmax);
          //line3[i] = new TF1("line3", "[0]+[1]*x", mmin, mmax);
          line3[i]->SetParameters(p0,p1,p2,p3);
          
          //printf("Correct\n");
          hMixCorrected[i] = (TH1D*) hMix[i]->Clone(Form("MixCorrected_PtBin%d_%s",i,hname.Data()));
          for(Int_t j = 0; j< hMix[i]->GetNbinsX(); j++)
          {
            Double_t x          = hMix[i]->GetBinCenter(j);
            Double_t correction = line3[i]->Eval(x);
            Double_t corrected  = hMix[i]->GetBinContent(j)*correction;
            Double_t ecorrected = hMix[i]->GetBinError(j)  *correction;
            //printf("bin %d, center %f, mix %f, corrected %f\n",i,x,hMix[i]->GetBinContent(j),corrected);
            hMixCorrected[i]->SetBinContent(j, corrected);
            hMixCorrected[i]->SetBinError  (j,ecorrected);
          }

          // Subtract
          hSignal[i] ->Add(hMixCorrected[i],-1);
        }
        else
        {
          Float_t scale = hRatio[i]->GetBinContent(hRatio[i]->FindBin(0.35));
          hMix[i]->Scale(scale);
          // Subtract
          hSignal[i] ->Add(hMix[i],-1);
        }
      } // enough mixed entries
    } // mixed event
    
    
    // ---- Fit subtracted signal
    Int_t nMax = 0;
    if(kParticle=="Pi0")
    {
      nMax= hSignal[i]->Integral(hIM[i]->FindBin(0.1),
                                 hIM[i]->FindBin(0.2));
    }
    else
    {
      nMax = hSignal[i]->Integral(hIM[i]->FindBin(0.4),
                                  hIM[i]->FindBin(0.65));
    }
    
    if(nMax > kNPairCut)
    {
      fitfun->SetParLimits(0,nMax/100,nMax*100);
      
      if(kParticle=="Pi0")
      {
        if(kPolN < 4 )fitfun->SetParameters(nMax/5,0.135,20,0);
        else         fitfun->SetParameters(nMax/5,20,0.135,20,20,nMax/5);
        
        if(i < nPt-4)
          hSignal[i]->Fit("fitfun","QR","",0.11,0.3);
        else
          hSignal[i]->Fit("fitfun","QR","",0.11,0.3);
      }// pi0
      else // eta
      {
        
        if ( kPolN < 4 ) fitfun->SetParameters(nMax/5,0.54,40,0);
        else             fitfun->SetParameters(nMax/5,40,0.54,40,40,nMax/5);
        //if(i<4)
        hSignal[i]->Fit("fitfun","QR","",0.4,0.7);
        //else if(i < 15)
        //  hSignal[i]->Fit("fitfun","QRL","",0.3,0.8);
        //else 
        //hSignal[i]->Fit("fitfun","QRL","",0.3,0.8);
      }
    }
    
    // Init results arrays
     mesonPt   .SetAt(-1,i);
    emesonPt   .SetAt(-1,i);    
     mesonPtBC .SetAt(-1,i);
    emesonPtBC .SetAt(-1,i);
     mesonMass .SetAt(-1,i);
    emesonMass .SetAt(-1,i);
     mesonWidth.SetAt(-1,i);
    emesonWidth.SetAt(-1,i);
    
    fitFunction = (TF1*) hSignal[i]->GetFunction("fitfun");
    
    if ( fitFunction )
    {
      //            printf("ipt %d: Chi/NDF %f, Chi %f, NDF %d\n",i,
      //                   fitFunction->GetChisquare()/fitFunction->GetNDF(),
      //                   fitFunction->GetChisquare(),
      //                   fitFunction->GetNDF());
      
      if(fitFunction->GetChisquare()/fitFunction->GetNDF()<1000)
      {
        Double_t A    = fitFunction->GetParameter(0);
        Double_t mean = fitFunction->GetParameter(1);
        Double_t sigm = fitFunction->GetParameter(2);
        //              Double_t a0   = fitFunction->GetParameter(3);
        //              Double_t a1   = fitFunction->GetParameter(4);
        //              Double_t a2   = 0;
        //              if (kPolN == 2)
        //                a2 = fitFunction->GetParameter(5);
        
        Double_t eA    = fitFunction->GetParError(0);
        Double_t emean = fitFunction->GetParError(1);
        Double_t esigm = fitFunction->GetParError(2);
        //              Double_t ea0   = fitFunction->GetParError(3);
        //              Double_t ea1   = fitFunction->GetParError(4);
        //              Double_t ea2   = 0;
        //              if (kPolN == 2)
        //                ea2 = fitFunction->GetParError(5);
        pLegendIM[i] = new TLegend(0.55,0.65,0.95,0.93);
        pLegendIM[i]->SetTextSize(0.035);
        pLegendIM[i]->SetFillColor(10);
        pLegendIM[i]->SetBorderSize(1);
        pLegendIM[i]->SetHeader(Form("#pi %s",leg.Data()));
        pLegendIM[i]->AddEntry("",Form("A = %2.2f #pm %2.2f ",A,eA),"");
        pLegendIM[i]->AddEntry("",Form("#mu = %2.3f #pm %2.3f ",mean,emean),"");
        pLegendIM[i]->AddEntry("",Form("#sigma = %2.3f #pm %2.3f ",sigm,esigm),"");
        //pLegendIM[i]->AddEntry("",Form("p_{0} = %2.1f #pm %2.1f ",a0,ea0),"");
        //pLegendIM[i]->AddEntry("",Form("p_{1} = %2.1f #pm %2.1f ",a1,ea1),"");
        //if(kPolN==2)pLegendIM[i]->AddEntry("",Form("p_{2} = %2.1f#pm %2.1f ",a2,ea2),"");
        
        Double_t  counts = A*sigm / mStep * TMath::Sqrt(TMath::TwoPi());
        Double_t eCounts =
        TMath::Power(eA/A,2) +
        TMath::Power(esigm/sigm,2);
        eCounts = TMath::Sqrt(eCounts)* counts;
        eCounts = TMath::Min(eCounts, TMath::Sqrt(counts));
        
         counts/=(nEvt*(exPt.At(i)*2));
        eCounts/=(nEvt*(exPt.At(i)*2));
        
         mesonPt.SetAt( counts,i);
        emesonPt.SetAt(eCounts,i);
        
        //cout << "N(pi0) fit  = " <<  counts << "+-"<<  eCounts << " mstep " << mStep<<endl;
        
         mesonMass .SetAt( mean*1000.,i);
        emesonMass .SetAt(emean*1000.,i);
        
         mesonWidth.SetAt( sigm*1000.,i);
        emesonWidth.SetAt(esigm*1000.,i);
      } //Good fit
     
    }
    
    Double_t mass  = mesonMass [i];
    Double_t width = mesonWidth[i];
    Double_t mMinBin =   mass - 2 * width;
    Double_t mMaxBin =   mass + 2 * width;
    if(mass <= 0 || width <= 0)
    {
      if(kParticle=="Pi0")
      {
        mMinBin = 0.115; mMaxBin = 0.3 ;
      }
      else{
        mMinBin = 0.4;  mMaxBin = 0.9 ; 
      }
      
    }
    
    //
    // Bin counting instead of fit 
    //
    if ( kMix && hSignal[i] )
    {
      Double_t  countsBin = hSignal[i]->Integral(hSignal[i]->FindBin(mMinBin),  hSignal[i]->FindBin(mMaxBin)) ;
      Double_t eCountsBin = TMath::Sqrt(countsBin);
      
       countsBin/=(nEvt*(exPt.At(i)*2));
      eCountsBin/=(nEvt*(exPt.At(i)*2));
      
       mesonPtBC.SetAt( countsBin,i);
      emesonPtBC.SetAt(eCountsBin,i);
      //printf("pt bin %d - N(pi0) BC = %2.3e +- %2.4e\n",i, countsBin, eCountsBin);
      
      //printf("Bin %d, [%1.1f,%1.1f]\n",i,xPt[i],xPt[i+1]);
    }
    
    if(nMax > kNPairCut)
    {
      if( kMix )
      {
        hIM[i]->Draw("HE");
        pLegendIM[i]->AddEntry(hIM    [i],"Raw pairs"   ,"L");

        if ( hMix[i]->GetEntries() > 100 )
        {
          if ( kTrunMixFunc )
          {
            hMixCorrected[i]->Draw("same");
            pLegendIM[i]->AddEntry(hMixCorrected[i],"Mixed pairs" ,"L");
            
          }
          else            
          {
            hMix[i]->Draw("same");
            pLegendIM[i]->AddEntry(hMix[i],"Mixed pairs" ,"L");
          }
        }
        
        hSignal[i] -> Draw("same");
        pLegendIM[i]->AddEntry(hSignal[i],"Signal pairs","L");
      }
      else 
      {
        pLegendIM[i]->AddEntry(hSignal[i],"Raw pairs","L");
        hSignal[i] -> Draw("HE");
      }
    
      // Plot mass from pairs originated from pi0/eta
      if ( hname.Contains("All") )
      {
        TH1F* hMesonPairOrigin = (TH1F*) fout->Get(Form("IM_MCpTReco_PtBin%d",i));
        if ( hMesonPairOrigin )
        {
          hMesonPairOrigin->SetLineColor(8);
          hMesonPairOrigin->Draw("same");
          pLegendIM[i]->AddEntry(hMesonPairOrigin,Form("%s origin",particleN.Data()),"L");
        }
      } // MC origin
      
//      if ( fitFunction ) 
//        pLegendIM[i]->AddEntry(fitFunction,"Fit","L");

      pLegendIM[i]->Draw();
    }
    
  } //pT bins
  
  cIMModi->Print(Form("IMfigures/%s_%s_Mgg_%s_%s.%s",
                      kProdName.Data(),kCalorimeter.Data(),kParticle.Data(),/*kFileName.Data(),*/
                      hname.Data(),kPlotFormat.Data()));
  
  //xxxx Real / Mixxxxx
  if(kMix)
  {
    //printf("Do real/mix\n");
    
    TCanvas  * cRat = new TCanvas(Form("Ratio_%s\n",hname.Data()), Form("Ratio %s\n", leg.Data()), 1200, 1200) ;
    cRat->Divide(col, col);
    
    for(Int_t i = 0; i < nPt; i++)
    {
      cRat->cd(i+1) ;
      //gPad->SetGridy();
      //gPad->SetLog();
      if(!hRatio[i]) 
      {
        //printf("No ratio in pt bin %d continue\n",i);
        continue;
      }
      
      Double_t nMax =0;
      if(kParticle=="Pi0")
        nMax = hRatio[i]->Integral(hRatio[i]->FindBin(0.005),
                                   hRatio[i]->FindBin(0.20));
      else
        nMax = hRatio[i]->Integral(hRatio[i]->FindBin(0.4),
                                   hRatio[i]->FindBin(0.6));
      
      hRatio[i]->SetMaximum(nMax/4);
      hRatio[i]->SetMinimum(1e-6);
      hRatio[i]->SetAxisRange(mmin,mmax,"X");
      
      hRatio[i]->Draw();
      
      //if(line3[i]) line3[i]->Draw("same");
    }
    
    cRat->Print(Form("IMfigures/%s_%s_MggRatio_%s_%s.%s",
                     kProdName.Data(),kCalorimeter.Data(),kParticle.Data(),/*kFileName.Data(),*/ 
                     hname.Data(),kPlotFormat.Data()));
  }    
  
  //------------------------------
  // Fit parameters 
  //------------------------------
  // Put fit results in TGraphErrors
  
  TGraphErrors* gPt   =  new TGraphErrors(nPt,xPt.GetArray(),mesonPt.GetArray(),exPt.GetArray(),emesonPt.GetArray());
  gPt->SetName(Form("gPt_%s",hname.Data()));
  gPt->GetHistogram()->SetTitle(Form("p_{T} of reconstructed %s, %s ",particleN.Data(), comment.Data()));
  gPt->GetHistogram()->SetXTitle("p_{T} (GeV/c)    ");
  gPt->GetHistogram()->SetYTitle(Form("dN_{%s}/dp_{T} (GeV/c)^{-1}  / N_{events}  ",particleN.Data()));
  gPt->GetHistogram()->SetAxisRange(0.,30);
  gPt->GetHistogram()->SetTitleOffset(1.5,"Y");
  gPt->SetMarkerStyle(20);
  gPt->SetMarkerSize(1);
  gPt->SetMarkerColor(1);
  //   gPt->GetHistogram()->SetMaximum(1e8);
  //   gPt->GetHistogram()->SetMinimum(1e-8);
  
  TGraphErrors* gPtBC =  new TGraphErrors(nPt,xPt.GetArray(),mesonPtBC.GetArray(),exPt.GetArray(),emesonPtBC.GetArray());
  gPtBC->SetName(Form("gPtBC_%s",hname.Data()));
  gPtBC->GetHistogram()->SetTitle(Form("p_{T} of reconstructed %s, %s ",particleN.Data(), comment.Data()));
  gPtBC->GetHistogram()->SetXTitle("p_{T} (GeV/c)    ");
  gPtBC->GetHistogram()->SetYTitle(Form("dN_{%s}/dp_{T} (GeV/c)^{-1}  / N_{events}  ",particleN.Data()));
  gPtBC->GetHistogram()->SetAxisRange(0.,30);
  gPtBC->GetHistogram()->SetTitleOffset(1.5,"Y");
  gPtBC->SetMarkerStyle(24);
  gPtBC->SetMarkerSize(1);
  gPtBC->SetMarkerColor(1);
  //   gPtBC->GetHistogram()->SetMaximum(1e8);
  //   gPtBC->GetHistogram()->SetMinimum(1e-8);
  
  TGraphErrors* gMass =  new TGraphErrors(nPt,xPt.GetArray(),mesonMass.GetArray(),exPt.GetArray(),emesonMass.GetArray());
  gMass->SetName(Form("gMass_%s",hname.Data()));
  gMass->GetHistogram()->SetTitle(Form("mass vs p_{T} of reconstructed %s, %s ",particleN.Data(), comment.Data()));
  gMass->GetHistogram()->SetXTitle("p_{T} (GeV/c)    ");
  gMass->GetHistogram()->SetYTitle(Form("mass_{%s} (MeV/c)^{2}    ",particleN.Data()));
  //gMass->GetHistogram()->SetAxisRange(0.,30);
  gMass->GetHistogram()->SetTitleOffset(1.5,"Y");
  gMass->SetMarkerStyle(20);
  gMass->SetMarkerSize(1);
  gMass->SetMarkerColor(1);
  
  TGraphErrors* gWidth=  new TGraphErrors(nPt,xPt.GetArray(),mesonWidth.GetArray(),exPt.GetArray(),emesonWidth.GetArray());
  gWidth->SetName(Form("gWidth_%s",hname.Data()));
  gWidth->GetHistogram()->SetTitle(Form("Width vs p_{T} of reconstructed %s, %s ",particleN.Data(), comment.Data()));
  gWidth->GetHistogram()->SetXTitle("p_{T} (GeV/c)    ");
  gWidth->GetHistogram()->SetYTitle(Form("#sigma_{%s} (MeV/c)^{2}    ",particleN.Data()));
  //gWidth->GetHistogram()->SetAxisRange(0.,30);
  gWidth->GetHistogram()->SetTitleOffset(1.5,"Y");
  gWidth->SetMarkerStyle(20);
  gWidth->SetMarkerSize(1);
  gWidth->SetMarkerColor(1);
  
  // Plot the fitted results
  
  TCanvas *cFitGraph = new TCanvas(Form("cFitGraph_%s",hname.Data()),
                                   Form("Fit Graphs for %s",hname.Data()),600,600);
  cFitGraph->Divide(2, 2);

  // Mass  
  cFitGraph->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  
  if ( kParticle=="Pi0" )
  {
    gMass->SetMaximum(160);
    gMass->SetMinimum(100);
  }
  else
  {
    gMass->SetMaximum(660);
    gMass->SetMinimum(420);
  }
  
  gMass->Draw("AP");
  
  cFitGraph->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  
  if(kParticle=="Pi0")
  {
    gWidth->SetMaximum(16);
    gWidth->SetMinimum(6);
  }
  else
  {
    gWidth->SetMaximum(60);
    gWidth->SetMinimum(0);
  }
  
  gWidth->Draw("AP");
 
  cFitGraph->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  
  gPt->SetMaximum(1);
  gPt->SetMinimum(1e-8);
 
  gPt->Draw("AP");
    
  if(kMix)
  {
    gPtBC ->Draw("P");
    
    TLegend * legend =  new TLegend(0.5,0.7,0.9,0.9);
    legend->AddEntry(gPt  ,"From fit","P");
    legend->AddEntry(gPtBC,"From bin counting","P");
    legend->Draw();
  }
  
  cFitGraph->Print(Form("IMfigures/%s_%s_MassWidthPtSpectra_%s_%s.%s",
                        kProdName.Data(),kCalorimeter.Data(),kParticle.Data(),hname.Data(),
                        /*kFileName.Data(),*/ kPlotFormat.Data())); 
  
  //-----------------------
  // Write results to file
  //-----------------------
  hRe->Write();
  if ( hMi ) hMi->Write();
  
  for(Int_t ipt = 0; ipt < nPt; ipt++)
  {
    //if(hIM[ipt]) { hIM[ipt]->Scale(1./nEvt); hIM[ipt]->Write(); }
    if(hIM[ipt]) { hIM[ipt]->Write(); }
    if(kMix)
    {
      //        if(hMix         [ipt]) { hMix         [ipt]->Scale(1./nEvt); hMix         [ipt]->Write();}
      //        if(hMixCorrected[ipt]) { hMixCorrected[ipt]->Scale(1./nEvt); hMixCorrected[ipt]->Write();}
      //        if(hRatio       [ipt]) { hRatio       [ipt]->Scale(1./nEvt); hRatio       [ipt]->Write();}
      //        if(hSignal      [ipt]) { hSignal      [ipt]->Scale(1./nEvt); hSignal      [ipt]->Write();}
      if(hMix         [ipt]) { hMix         [ipt]->Write();}
      if(hMixCorrected[ipt]) { hMixCorrected[ipt]->Write();}
      if(hRatio       [ipt]) { hRatio       [ipt]->Write();}
      if(hSignal      [ipt]) { hSignal      [ipt]->Write();}
    }
  }
  
  gPt   ->Write();
  gMass ->Write();
  gWidth->Write();  
  if(kMix) gPtBC->Write();

  // Do some final efficiency calculations if MC   
  Efficiency(nPt, xPt, exPt,mesonPt,mesonPtBC,hname);
}

///
/// Plot combinations of graphs with fit results
/// SM/Sector/Side/SM group by SM/Sector/Side/SM group comparisons
///
//-----------------------------------------------------------------------------
void PlotGraphs(TString opt, Int_t first, Int_t last)
{
  const Int_t nCombi = last-first+1;
  
  TGraphErrors* gPt   [nCombi];
  TGraphErrors* gPtBC [nCombi];
  TGraphErrors* gMass [nCombi];
  TGraphErrors* gWidth[nCombi];
  
  Float_t xmin = 0.7;
  Float_t ymin = 0.3;
  Float_t xmax = 0.9;
  Float_t ymax = 0.9;
  
  if(opt.Contains("Group"))
  {
    xmin = 0.3;
    ymin = 0.7;
  }
  
  TLegend * legend =  new TLegend(xmin,ymin,xmax,ymax);
  legend->SetTextSize(0.05);
  
  TString particleN = " #pi^{0}";
  if(kParticle=="Eta") particleN = " #eta";
  
  // recover the graphs from output file
  //printf("%s\n",opt.Data());
  for(Int_t icomb = 0; icomb < nCombi; icomb ++)
  {
    gPt   [icomb] = (TGraphErrors*) fout->Get(Form("gPt_%s%d"   ,opt.Data(),icomb+first));
    gPtBC [icomb] = (TGraphErrors*) fout->Get(Form("gPtBC_%s%d" ,opt.Data(),icomb+first));
    gMass [icomb] = (TGraphErrors*) fout->Get(Form("gMass_%s%d" ,opt.Data(),icomb+first));
    gWidth[icomb] = (TGraphErrors*) fout->Get(Form("gWidth_%s%d",opt.Data(),icomb+first));
    //printf("\t %d %p %p %p %p\n",icomb,gPt[icomb],gPtBC[icomb],gMass[icomb],gWidth[icomb]);
    
    gPt   [icomb]->SetMarkerStyle(modStyleIndex[icomb]);
    gPt   [icomb]->SetMarkerColor(modColorIndex[icomb]);
    gPt   [icomb]->SetLineColor  (modColorIndex[icomb]);    
    
    gMass[icomb]->SetMarkerStyle(modStyleIndex[icomb]);
    gMass[icomb]->SetMarkerColor(modColorIndex[icomb]);
    gMass[icomb]->SetLineColor  (modColorIndex[icomb]);   
    
    gWidth[icomb]->SetMarkerStyle(modStyleIndex[icomb]);
    gWidth[icomb]->SetMarkerColor(modColorIndex[icomb]);
    gWidth[icomb]->SetLineColor  (modColorIndex[icomb]);
        
    gPt   [icomb]->GetHistogram()->SetTitle(Form("p_{T} of %s from fit, %s "   ,particleN.Data(), opt.Data()));
    gWidth[icomb]->GetHistogram()->SetTitle(Form("%s mass #sigma vs p_{T}, %s ",particleN.Data(), opt.Data()));
    gMass [icomb]->GetHistogram()->SetTitle(Form("%s mass #mu vs p_{T}, %s "   ,particleN.Data(), opt.Data()));
    
    if ( kMix )
    {
      gPtBC [icomb]->SetMarkerStyle(modStyleIndex[icomb]);
      gPtBC [icomb]->SetMarkerColor(modColorIndex[icomb]);
      gPtBC [icomb]->SetLineColor  (modColorIndex[icomb]);   
      gPtBC [icomb]->GetHistogram()->SetTitle(Form("p_{T} of %s bin counted, %s ",particleN.Data(), opt.Data()));
    }
    
    if ( !opt.Contains("Group") )
      legend->AddEntry(gPt[icomb],Form("%s %d",opt.Data(),icomb),"P");
    else
    {
      if(icomb == 0) legend->AddEntry(gPt[icomb],"SM0+4+5+6+7+8+9","P");
      if(icomb == 1) legend->AddEntry(gPt[icomb],"SM1+2","P");
      if(icomb == 2) legend->AddEntry(gPt[icomb],"SM3+7","P");
    }
  }
  
 
  
  TCanvas *cFitGraph = new TCanvas(Form("cFitGraph_%sCombinations",opt.Data()),
                                   Form("Fit Graphs for %s combinations",opt.Data()),600,600);
  cFitGraph->Divide(2, 2);
  
  // Mass  
  cFitGraph->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  
//  gMass[0]->SetMaximum(1);
//  gMass[0]->SetMinimum(1e-8);
  
  gMass[0]->Draw("AP");
  for(Int_t icomb = 1; icomb < nCombi; icomb ++)
    gMass[icomb]->Draw("P");
  
  legend->Draw();

  cFitGraph->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  
//  gWidth[0]->SetMaximum(1);
//  gWidth[0]->SetMinimum(1e-8);
  
  gWidth[0]->Draw("AP");
  for(Int_t icomb = 1; icomb < nCombi; icomb ++)
    gWidth[icomb]->Draw("P");

  legend->Draw();

  cFitGraph->cd(3);

  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  
//  gPt[0]->SetMaximum(1);
//  gPt[0]->SetMinimum(1e-8);
  
  gPt[0]->Draw("AP");
  for(Int_t icomb = 1; icomb < nCombi; icomb ++)
    gPt[icomb]->Draw("P");
  
  legend->Draw();

  if(kMix)
  {
    cFitGraph->cd(4);

//    gPtBC[0]->SetMaximum(1);
//    gPtBC[0]->SetMinimum(1e-8);
    
    gPtBC[0]->Draw("P");
    for(Int_t icomb = 1; icomb < nCombi; icomb ++)
      gPtBC[icomb]->Draw("P");
    legend->Draw();
  }
  
  cFitGraph->Print(Form("IMfigures/%s_%s_MassWidthPtSpectra_%s_%sCombinations.%s",
                        kProdName.Data(),kCalorimeter.Data(),kParticle.Data(),opt.Data(),
                        /*kFileName.Data(),*/kPlotFormat.Data()));
}

///
/// Get primary generated spectra and detector acceptance
/// Plot them
///
//-----------------------------------------------------------------------------
void PrimarySpectra(Int_t nPt, TArrayD xPtLimits, TArrayD xPt, TArrayD exPt)
{
  TArrayD   primMesonPt   ;   primMesonPt   .Set(nPt); 
  TArrayD  eprimMesonPt   ;  eprimMesonPt   .Set(nPt); 
  TArrayD   primMesonPtAcc;   primMesonPtAcc.Set(nPt); 
  TArrayD  eprimMesonPtAcc;  eprimMesonPtAcc.Set(nPt); 
  TArrayD   primMesonAcc  ;   primMesonAcc.Set(nPt); 
  TArrayD  eprimMesonAcc  ;  eprimMesonAcc.Set(nPt); 
  
  TH2D* h2PrimMesonPtY    = (TH2D*) fil->Get( Form("%s_hPrim%sRapidity"   ,kHistoStartName.Data(), kParticle.Data() ) );
  TH2D* h2PrimMesonPtYAcc = (TH2D*) fil->Get( Form("%s_hPrim%sAccRapidity",kHistoStartName.Data(), kParticle.Data() ) );
  
  if ( !h2PrimMesonPtY ) return ;
  
  TH1D* hY = h2PrimMesonPtY->ProjectionY("Rapidity",-1,-1);
  Int_t binYmin = hY->FindBin(+0.65);
  Int_t binYmax = hY->FindBin(-0.65);
  TH1D* hPrimMesonPt = (TH1D*) h2PrimMesonPtY->ProjectionX("PrimaryPt" ,binYmin   ,binYmax   );
  hPrimMesonPt->Write();
  
  TH1D* hYAcc = h2PrimMesonPtY->ProjectionY("RapidityAcc",-1,-1);
  binYmin = hYAcc->FindBin(+0.65);
  binYmax = hYAcc->FindBin(-0.65);
  TH1D* hPrimMesonPtAcc = (TH1D*) h2PrimMesonPtYAcc->ProjectionX("PrimaryPtAccepted" ,binYmin   ,binYmax   );
  hPrimMesonPtAcc->Write();
  
  for(Int_t ibin = 0; ibin < nPt; ibin++ )
  {
    Double_t ptMin = xPtLimits[ibin];
    Double_t ptMax = xPtLimits[ibin+1];
    
     primMesonPt[ibin] = hPrimMesonPt->Integral(hPrimMesonPt->FindBin(ptMin), hPrimMesonPt->FindBin(ptMax));
    eprimMesonPt[ibin] = TMath::Sqrt(primMesonPt[ibin]);
    
     primMesonPtAcc[ibin] = hPrimMesonPt->Integral(hPrimMesonPtAcc->FindBin(ptMin), hPrimMesonPtAcc->FindBin(ptMax));
    eprimMesonPtAcc[ibin] = TMath::Sqrt(primMesonPtAcc[ibin]);
    
    if(primMesonPtAcc[ibin] && primMesonPt[ibin])
    {
       primMesonAcc[ibin] = primMesonPtAcc[ibin] / primMesonPt[ibin] ;
      eprimMesonAcc[ibin] =  primMesonAcc[ibin] * TMath::Sqrt(1./(primMesonPtAcc[ibin] * primMesonPtAcc[ibin]) + 
                                                              1./ primMesonPt   [ibin] * primMesonPt   [ibin]) ;
    }
    else
    {
       primMesonAcc[ibin] = 0;
      eprimMesonAcc[ibin] = 0;
    }
    
    //printf("Bin %d, [%1.1f,%1.1f], content num %f, content den %f times evt %e, bin size %f\n",ibin,xPt[ibin],xPt[ibin+1],mesonPtBC[0][ibin],denBin,nEvt,exPt[ibin]*2);
    
     primMesonPt[ibin]/=(nEvt*(exPt[ibin]*2));
    eprimMesonPt[ibin]/=(nEvt*(exPt[ibin]*2));
    
    //        // Scale biased production
    //        if(kFileName.Contains("pp_7TeV_Pi0"))
    //        {
    //          primMeson[ibin]*=1.e6;
    //          primMeson[ibin]*=360./100.;
    //          
    //          //eprimMeson[ibin]*=1.e6*1.e6;
    //          eprimMeson[ibin]*=360./100.;
    //          
    //          mesonPtBC [0][ibin] /= (0.62 +xPt[ibin]*0.0017 );
    //          emesonPtBC[0][ibin] /= (0.62 +xPt[ibin]*0.0017 );
    //          
    //          primMeson [ibin] /= (0.56 +xPt[ibin]*0.0096 );
    //          eprimMeson[ibin] /= (0.56 +xPt[ibin]*0.0096 );
    //        }
  } // pT bin loop
  
  TGraphErrors * gPrim =  new TGraphErrors(nPt,xPt.GetArray(),primMesonPt.GetArray(),exPt.GetArray(),eprimMesonPt.GetArray());
  gPrim->SetName("Primary");
  gPrim->GetHistogram()->SetYTitle("dN/d p_{T}");
  gPrim->GetHistogram()->SetTitle("|Y| < 0.65");
  gPrim->Write();
  
  TGraphErrors * gPrimAcc =  new TGraphErrors(nPt,xPt.GetArray(),primMesonPtAcc.GetArray(),exPt.GetArray(),eprimMesonPtAcc.GetArray());
  gPrimAcc->SetName("PrimaryInAcceptance");
  gPrimAcc->GetHistogram()->SetYTitle("dN/d p_{T}");
  gPrimAcc->GetHistogram()->SetTitle(Form("|Y| < 0.65 in %s",kCalorimeter.Data()));
  gPrimAcc->Write();    

  // Acceptance
  TH1D* hAcc = (TH1D*) hPrimMesonPtAcc->Clone("hAcceptance");
  hAcc->Divide(hPrimMesonPt);
  
  TGraphErrors * gAcc =  new TGraphErrors(nPt,xPt.GetArray(),primMesonAcc.GetArray(),exPt.GetArray(),eprimMesonAcc.GetArray());
  gAcc->SetName("Acceptance");
  gAcc->GetHistogram()->SetYTitle("Acceptance");
  gAcc->GetHistogram()->SetTitle(Form("Acceptance for |Y| < 0.65 in %s",kCalorimeter.Data()));
  gAcc->Write();    
 
  // Plot spectra and acceptance
  //
  TCanvas *cAcc = new TCanvas(Form("cAcceptance"),
                              Form("Primary generated p_{T} spectra and acceptance"),2*600,600);
  cAcc->Divide(2, 1);
  
  cAcc->cd(1);
  //gPad->SetGridx();
  //gPad->SetGridy();
  gPad->SetLogy();
  
  //gPrim->SetMaximum(1);
  //Prim->SetMinimum(1e-8);
  
  gPrim   ->Draw("AP");
  gPrimAcc->Draw("P");
  
  TLegend * legendS =  new TLegend(0.4,0.7,0.9,0.9);
  legendS->AddEntry(hAcc,"|Y| < 0.65","P");
  legendS->AddEntry(gAcc,Form("Both in %s",kCalorimeter.Data()),"P");
  legendS->Draw();
  
  cAcc->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  //gPad->SetLogy();
  
  //gAcc->SetMaximum(1);
  //gAcc->SetMinimum(1e-8);
  
  hAcc->Draw("HE");
  gAcc->Draw("P");
  
  TLegend * legendA =  new TLegend(0.5,0.15,0.9,0.3);
  legendA->AddEntry(hAcc,"Histo","P");
  legendA->AddEntry(gAcc,"Graph","P");
  legendA->Draw();
  
  cAcc->Print(Form("IMfigures/%s_%s_PrimarySpectraAcceptance_%s.%s",
                   kProdName.Data(),kCalorimeter.Data(),kParticle.Data(),
                   /*kFileName.Data(),*/ kPlotFormat.Data()));   
}
  
///
/// Main method
///
/// \param prodname   : name of directory with histogram file
/// \param filename   : histogram file name
/// \param histoDir   : TDirectoryFile folder name
/// \param histoList  : TList folder name
/// \param calorimeter: "EMCAL","DCAL"
/// \param particle   : "Pi0","Eta", define fitting and plotting ranges for particle
/// \param mixed      : bool, use mixed event to constrain combinatorial background
/// \param truncmix   : bool, 1: use truncated function to constrain backgroun; 0: use factor from fixed bin
/// \param pol        : int, polinomyal type for residual background under the peak
/// \param ncomb      : total number of SM combinations (Single SM, same side 2 SM, same sector 2 SM)
/// \param drawAllCombi: Plot also many SM combinations
/// \param nPairMin   : minimum number of entries un the pi0 or eta peak integral to do the fits and plotting. Careful in MC scaled productions.
/// \param plotFormat : define the type of figures: eps, pdf, etc.
///
//-----------------------------------------------------------------------------
void InvMassFit
(TString prodname     = "LHC18c3_NystromOn",//"LHC17l3b_fast", 
 TString filename     = "AnalysisResults",
 TString histoDir     = "Pi0IM_GammaTrackCorr_EMCAL", 
 TString histoList    = "default", 
 TString calorimeter  = "EMCAL",
 TString particle     = "Pi0",
 Bool_t  mixed        = 1,
 Bool_t  truncmix     = 1,
 Int_t   pol          = 1,
 Int_t   ncomb        = 1,
 Bool_t  drawAllCombi = 0, 
 Float_t nPairMin     = 20, 
 TString plotFormat   = "eps")
{
  kPolN       = pol;
  kMix        = mixed;
  kTrunMixFunc= truncmix;
  kParticle   = particle;
  kProdName   = prodname;
  kFileName   = filename;
  kPlotFormat = plotFormat;
  kCalorimeter= calorimeter;
  kNPairCut = nPairMin;
  //  kNPairCut/=nEvt;
  
  if ( !drawAllCombi ) ncomb = 1;
  
  kHistoStartName = "AnaPi0_Calo0";
  if(kCalorimeter=="DCAL") 
    kHistoStartName = "AnaPi0_Calo1";
  
  //---------------------------------------
  // Get input file
  //---------------------------------------
  Int_t ok = GetFileAndEvents(histoDir, histoList);
  
  printf("Settings: prodname %s, filename %s, histoDir %s, histoList %s, particle %s, calorimeter %s, \n"
         " \t mix %d, kPolN %d, n combi %d, n pairs cut %2.2e\n",
         kProdName.Data(),kFileName.Data(),histoDir.Data(),histoList.Data(), kParticle.Data(), 
         kCalorimeter.Data(), kMix, kPolN, ncomb, kNPairCut);
  
  if ( !ok )
  {
    printf("Could not recover file <%p> or dir <%p> or list <%p>\n",fil,direc,lis);
    return ; 
  }
  
  //---------------------------------------
  // Set-up output file with processed histograms
  //---------------------------------------
  
  // Open output file
  fout = new TFile(Form("IMfigures/%s_%s_MassWidthPtHistograms_%s.root",
                        kProdName.Data(),kCalorimeter.Data(),kParticle.Data() /*,kFileName.Data()*/), "recreate");
  printf("Open Output File: IMfigures/%s_%s_MassWidthPtHistograms_%s.root\n",
         kProdName.Data(),kCalorimeter.Data(),kParticle.Data());//,kFileName.Data());
  
  //---------------------------------------
  // Set the pt bins and total range
  //---------------------------------------

  const Int_t nPtEta = 8;
  Double_t xPtLimitsEta[] = {2.,4.,6.,8.,10.,14.,18.,24.};
  
  const Int_t nPtPi0 = 12;
  Double_t xPtLimitsPi0[] = {2.,3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,18.};
  
  Int_t nPt = nPtPi0;
  if(kParticle == "Eta") nPt = nPtEta;
  
  TArrayD xPtLimits; xPtLimits.Set(nPt+1); 
  TArrayD xPt      ; xPt      .Set(nPt);
  TArrayD exPt     ; exPt     .Set(nPt);
  
  for(Int_t i = 0; i < nPt; i++)
  {
    if(kParticle == "Pi0")
    {
      xPtLimits.SetAt(xPtLimitsPi0[i],i);
      exPt     .SetAt((xPtLimitsPi0[i+1]-xPtLimitsPi0[i])/2,i);
      xPt      .SetAt(xPtLimitsPi0[i]+exPt[i],i);
    }
    else
    {
      xPtLimits.SetAt(xPtLimitsEta[i],i);
      exPt     .SetAt((xPtLimitsEta[i+1]-xPtLimitsEta[i])/2,i);
      xPt      .SetAt(xPtLimitsEta[i]+exPt[i],i);
    }
    
    //printf("%s pT bin %d, pT %2.2f, dpT %2.2f, limit %2.2f\n",kParticle.Data(),i,xPt[i],exPt[i],xPtLimits[i]);
  }
  // Extra entry
  if ( kParticle == "Pi0" ) xPtLimits.SetAt(xPtLimitsPi0[nPt],nPt);
  else                      xPtLimits.SetAt(xPtLimitsEta[nPt],nPt);
  
  TString comment = "";
  TString leg     = "";
  TString hname   = "";   
 
  //---------------------------------------
  // Fitting function and histograms initialization
  //---------------------------------------
  SetFitFun();
  
  TH2F  * hRe ;
  TH2F  * hMi ;
  
  //=======================================
  // Get histograms, project per pT bin and fit
  // Do it for different kind of histograms: 
  //     Total pairs, MC pairs, pairs per SM combinations
  //=======================================
  
  //---------------------------------------
  // MC input
  //---------------------------------------
  
  {
    // Get the generated primary pi0 spectra, for efficiency calculation 
    PrimarySpectra(nPt,xPtLimits, xPt, exPt);
    
    // Reconstructed pT of real meson pairs
    //hRe = (TH2F *) fil->Get( Form("%s_hMCOrgMass_2", kHistoStartName.Data()) ); // Pi0
    //hRe = (TH2F *) fil->Get( Form("%s_hMCOrgMass_3", kHistoStartName.Data()) ); // Eta
    hRe = (TH2F *) fil->Get( Form("%s_hMC%sMassPtRec", kHistoStartName.Data(), kParticle.Data()) );
    comment = "MC pT reconstructed";
    leg     = "MC pT reconstructed";
    hname   = "MCpTReco" ; 
    ProjectAndFit(nPt,comment,leg,hname,xPtLimits, xPt, exPt, hRe, 0x0 );    
    
    // Reconstructed pT of real eta pairs
    hRe = (TH2F *) fil->Get( Form("%s_hMC%sMassPtTrue", kHistoStartName.Data(), kParticle.Data()) );
    comment = "MC pT generated";
    leg     = "MC pT generated";
    hname   = "MCpTGener" ; 
    ProjectAndFit(nPt,comment,leg,hname,xPtLimits, xPt, exPt, hRe, 0x0 );    
  } // MC treatment
  
  //---------------------------------------
  // All SM
  //---------------------------------------
  
  {
    if ( !lis )
    {
      hRe = (TH2F *) fil->Get(Form("%s_hRe_cen0_pidbit0_asy0_dist1", kHistoStartName.Data()));
      hMi = (TH2F *) fil->Get(Form("%s_hMi_cen0_pidbit0_asy0_dist1", kHistoStartName.Data()));
    }
    else
    {
      hRe = (TH2F *) lis->FindObject(Form("%s_hRe_cen0_pidbit0_asy0_dist1", kHistoStartName.Data()));
      hMi = (TH2F *) lis->FindObject(Form("%s_hMi_cen0_pidbit0_asy0_dist1", kHistoStartName.Data()));
    }
    
    printf("histo Re %p - Mix %p All SM\n",hRe,hMi);
    
    if( !hMi ) kMix = kFALSE; 
    
    comment = "All SM";
    leg     = "All SM";
    hname   = "AllSM";
    ProjectAndFit(nPt,comment,leg,hname,xPtLimits, xPt, exPt, hRe, hMi );    
  }
  
  //-------------------------------------------
  // Histograms for cluster pairs in different 
  // combinations of SM
  //-------------------------------------------
  
  //
  // Pairs in same SM
  // 
  if ( drawAllCombi )
  {
    Int_t imod;
    Int_t firstMod = 0;
    Int_t lastMod  = 11;   
    Int_t firstSide = 0;
    Int_t lastSide  = 9;   
    Int_t firstSector = 0;
    Int_t lastSector  = 5;
    if ( kCalorimeter=="DCAL" )
    {
      firstMod  = 12;
      lastMod   = 19;
      firstSide = 10;
      lastSide  = 15;   
      firstSector = 6;
      lastSector  = 9;
    }
    
    // Initialize SM grouping histograms
    // 0 : SM 0+4+5+6+8+9
    // 1 : SM 1+2 
    // 2 : SM 3+7 
    TH2F  * hReSMGroups[3] ;
    TH2F  * hMiSMGroups[3] ;
    
    for(Int_t imodgroup = 0; imodgroup<=2; imodgroup++)
    {
      hReSMGroups[imodgroup] = hMiSMGroups[imodgroup] = 0x0;
    }
    
    for(imod = firstMod; imod<=lastMod; imod++)
    {
      if(!lis)
      {
        hRe = (TH2F *) fil->Get(Form("%s_hReMod_%d", kHistoStartName.Data(), imod));
        hMi = (TH2F *) fil->Get(Form("%s_hMiMod_%d", kHistoStartName.Data(), imod));
      }
      else
      {
        hRe = (TH2F *) lis->FindObject(Form("%s_hReMod_%d", kHistoStartName.Data(), imod));
        hMi = (TH2F *) lis->FindObject(Form("%s_hMiMod_%d", kHistoStartName.Data(), imod));
      }
      
      printf("histo mod %d Re %p - Mix %p\n",imod, hRe,hMi);
      
      comment = Form("both clusters in SM %d",imod);
      leg     = Form("SM %d",imod);
      hname   = Form("SM%d" ,imod);
      
      ProjectAndFit(nPt,comment,leg,hname,xPtLimits, xPt, exPt, hRe, hMi );    
      
      if      ( imod == 0 )
      {
        hReSMGroups[0] = (TH2F*) hRe->Clone("h2D_Re_IM_SM045689");
        if(kMix) hMiSMGroups[0] = (TH2F*) hMi->Clone("h2D_Mi_IM_SM045689");
      }
      else if ( imod == 1 )
      {
        hReSMGroups[1] = (TH2F*) hRe->Clone("h2D_Re_IM_SM12");
        if(kMix) hMiSMGroups[1] = (TH2F*) hMi->Clone("h2D_Mi_IM_SM12");
      }
      else if ( imod == 3 )
      {
        hReSMGroups[2] = (TH2F*) hRe->Clone("h2D_Re_IM_SM37");
        if(kMix) hMiSMGroups[2] = (TH2F*) hMi->Clone("h2D_Mi_IM_SM37");
      }     
      else if ( imod == 2 )
      {
        hReSMGroups[1]->Add(hRe);
        if(kMix) hMiSMGroups[1]->Add(hMi);
      }     
      else if ( imod == 7 )
      {
        hReSMGroups[2]->Add(hRe);
        if(kMix) hMiSMGroups[2]->Add(hMi);
      }
      else if ( imod < 10 )
      {
        hReSMGroups[0]->Add(hRe);
        if(kMix) hMiSMGroups[0]->Add(hMi);
      }
    }
    
    // Group SMs in with particular behaviors
    for(Int_t imodgroup = 0; imodgroup<=2; imodgroup++)
    {
      printf("histo modgroup %d Re %p - Mix %p\n",
             imod, hReSMGroups[imodgroup], hMiSMGroups[imodgroup]);
      
      comment = Form("both clusters in same SM, group %d",imodgroup);
      leg     = Form("SMGroup %d",imodgroup);
      hname   = Form("SMGroup%d" ,imodgroup);
      
      if ( imodgroup != 0 )
      {
        hReSMGroups[imodgroup]->Scale(1./2.);
        if ( kMix ) hMiSMGroups[imodgroup]->Scale(1./2.);
      }
      else
      {
        hReSMGroups[imodgroup]->Scale(1./6.);
        if ( kMix ) hMiSMGroups[imodgroup]->Scale(1./6.);
      }
      
      ProjectAndFit(nPt, comment, leg, hname, xPtLimits, xPt, exPt, 
                    hReSMGroups[imodgroup], hMiSMGroups[imodgroup] );    
    }
    
    //
    // Pairs in same sector
    //
    Int_t isector;
    for(isector = firstSector; isector<=lastSector; isector++)
    { 
      if(!lis)
      {
        hRe = (TH2F *) fil->Get(Form("%s_hReSameSectorEMCALMod_%d", kHistoStartName.Data(), isector));
        hMi = (TH2F *) fil->Get(Form("%s_hMiSameSectorEMCALMod_%d", kHistoStartName.Data(), isector));
      }     
      else 
      {
        hRe = (TH2F *) lis->FindObject(Form("%s_hReSameSectorEMCALMod_%d", kHistoStartName.Data(), isector));
        hMi = (TH2F *) lis->FindObject(Form("%s_hMiSameSectorEMCALMod_%d", kHistoStartName.Data(), isector));
      }
      
      printf("histo sector %d Re %p - Mix %p\n",isector, hRe,hMi);
      
      comment = Form("both clusters in Sector %d",isector);
      leg     = Form("Sector %d",isector);
      hname   = Form("Sector%d" ,isector);
      
      ProjectAndFit(nPt,comment,leg,hname,xPtLimits, xPt, exPt, hRe, hMi );    
    }
    
    // Pairs in same side
    //
    Int_t iside;
    for(iside = firstSide; iside <= lastSide; iside++)
    { 
      if(!lis)
      {
        hRe = (TH2F *) fil->Get(Form("%s_hReSameSideEMCALMod_%d", kHistoStartName.Data(), iside));
        hMi = (TH2F *) fil->Get(Form("%s_hMiSameSideEMCALMod_%d", kHistoStartName.Data(), iside));
      }
      else
      {
        hRe = (TH2F *) lis->FindObject(Form("%s_hReSameSideEMCALMod_%d", kHistoStartName.Data(), iside));
        hMi = (TH2F *) lis->FindObject(Form("%s_hMiSameSideEMCALMod_%d", kHistoStartName.Data(), iside));
      }
      
      printf("histo side %d Re %p - Mix %p\n",iside, hRe,hMi);
      
      comment = Form("both clusters in Side %d",iside);
      leg     = Form("Side %d",iside);
      hname   = Form("Side%d" ,iside);      
      
      ProjectAndFit(nPt,comment,leg,hname,xPtLimits, xPt, exPt, hRe, hMi );    
    }
    
    //------------------------------
    // Fit parameters final plotting
    //------------------------------
    
    PlotGraphs("SM"     ,firstMod   ,lastMod   );
    PlotGraphs("SMGroup",0          ,2         );
    PlotGraphs("Side"   ,firstSide  ,lastSide  );
    PlotGraphs("Sector" ,firstSector,lastSector);
    
  } // multiple combinations
  
  fout->Close();  
}
