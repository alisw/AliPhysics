///
/// \file InvMassFit.C
/// \ingroup CaloTrackCorrMacrosQA
/// \brief Fit invariant mass distributions
///
/// Macro using as input the 2D histograms mass vs pT of AliAnaPi0
/// For a given set of pT bins invariant mass plots are fitted and mass vs pT 
/// and width vs pT and neutral meson spectra plots are obtained
///
/// Based on old macros, to be properly cleaned
///
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TString.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TFile.h"
#include "TLegend.h"
#include "TObject.h"
#include "TDirectoryFile.h"
#include "TGraphErrors.h"
#include "TList.h"
#include <TGaxis.h>

#endif

//-----------------------------------------------------------------------------
Double_t pi0massP0(Double_t *x, Double_t *par);
Double_t pi0massP1(Double_t *x, Double_t *par);
Double_t pi0massP2(Double_t *x, Double_t *par);
Double_t pi0massP3(Double_t *x, Double_t *par);
Double_t truncatedPolPi0(Double_t *x, Double_t *par);
Double_t truncatedPolEta(Double_t *x, Double_t *par);
Double_t CrystalBall(Double_t *x, Double_t *par);

Bool_t  GetFileAndEvents( TString prodname, TString filename, TString dirName, TString listName );
void    SetFitFun();

// Settings
Bool_t  mix   = kFALSE;   /// use mixed event to constrain combinatorial background
TString part  = "Eta";    /// define fitting and plotting ranges for particle
Int_t   polN  = 1;        /// polinomyal type for residual background under the peak
Bool_t  sumw2 = kTRUE;    /// Apply Root method Sumw2()
Bool_t  drawAllCombi = 0; /// Plot also many SM combinations
Float_t nPairCut = 20;    /// Minimum number of cluster pairs in pi0 or eta window 

// Initializations
Int_t nEvt = 0;
TFile * fil   = 0;
TList * lis   = 0;
TDirectoryFile * direc =0;
TF1 *fitfun   = 0;

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
/// \param pol        : int, polinomyal type for residual background under the peak
/// \param ncomb      : total number of SM combinations (Single SM, same side 2 SM, same sector 2 SM)
/// \param nPairMin   : minimum number of entries un the pi0 or eta peak integral to do the fits and plotting. Careful in MC scaled productions.
/// \param fileFormat : define the type of figures: eps, pdf, etc.
///
//-----------------------------------------------------------------------------
void InvMassFit
(TString prodname     = "LHC18c3_NystromOn",//"LHC17l3b_fast", 
 TString filename     = "AnalysisResults",
 TString histoDir     = "Pi0IM_GammaTrackCorr_EMCAL", 
 TString histoList    = "default", 
 TString calorimeter  = "DCAL",
 TString particle     = "Pi0",
 Bool_t  mixed        = 0,
 Int_t   pol          = 1,
 Int_t   ncomb        = 1,
 Float_t nPairMin     = 20, 
 TString fileFormat   = "eps")
{
  part = particle;
  polN = pol;
  mix  = mixed;
  
  Int_t ok = GetFileAndEvents(prodname,filename,histoDir, histoList);

  nPairCut = nPairMin;
//  nPairCut/=nEvt;
  
  printf("Settings: prodname %s, filename %s, histoDir %s, histoList %s, particle %s, calorimeter %s, \n"
         " \t mix %d, polN %d, n combi %d, n pairs cut %2.2e\n",
         prodname.Data(),filename.Data(),histoDir.Data(),histoList.Data(), part.Data(), calorimeter.Data(), 
         mix, polN, ncomb,nPairCut);
  
  if ( !ok )
  {
    printf("Could not recover file <%p> or dir <%p> or list <%p>\n",fil,direc,lis);
    return ; 
  }
  
  //---------------------------------------
  // Set-up fitting and histogram plotting
  //---------------------------------------
  const Int_t ncombFix = 100;
  const Int_t nFixBins = 50; // to initialize arrays size
  
  const Int_t nPtEta = 8;
  Double_t xPtLimitsEta[] = {2.,4.,6.,8.,10.,14.,18.,24.};
  
  const Int_t nPtPi0 = 12;
  Double_t xPtLimitsPi0[] = {2.,3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,18.};

  Double_t xPtLimits[nFixBins]; 
  Double_t xPt      [nFixBins] ;
  Double_t exPt     [nFixBins] ;
  
  Int_t nPt = nPtPi0;
  if(part == "Eta") nPt = nPtEta;
  
  for(Int_t i = 0; i < nPt+1; i++)
  {
    if(part == "Pi0")
    {
      xPtLimits[i] = xPtLimitsPi0[i];
      exPt[i] = (xPtLimitsPi0[i+1]-xPtLimitsPi0[i])/2;
      xPt [i] = xPtLimitsPi0[i]+exPt[i];
    }
    else
    {
      xPtLimits[i] = xPtLimitsEta[i];
      exPt[i] = (xPtLimitsEta[i+1]-xPtLimitsEta[i])/2;
      xPt [i] = xPtLimitsEta[i]+exPt[i];
    }
    //printf("%s pT bin %d, pT %f, dpT %f\n",part.Data(),i,xPt[i],exPt[i]);
  }
    
  Double_t mmin = 0;
  Double_t mmax = 1;
  
  if(part=="Pi0")
  {
    mmin = 0.00;
    mmax = 0.45;
  }
  else // eta
  {
    mmin = 0.25;
    mmax = 0.95;
  }
  
  Int_t rebin  = 1 ;
  Int_t rebin2 = 2*rebin;
  
  //                        - SM                       SM-Sector       - Side                 Side
  Int_t modColorIndex[]={1 , 2, 2, 3, 3, 4, 4, 7, 7, 6, 6, 2, 3, 4, 7, 6, 2, 2, 3, 3, 4, 4, 6, 6};
  Int_t modStyleIndex[]={20,24,25,24,25,24,25,24,25,24,25,21,21,21,21,21,22,26,22,26,22,26,22,26};
  
  //Mix fit func
//  TF1 *tPolPi0 = new TF1("truncatedPolPi0", truncatedPolPi0, 0.02, 0.25, 2);
//  TF1 *tPolEta = new TF1("truncatedPolEta", truncatedPolEta, 0.2, 1, 2);
  
  // Fitting function initialization
  SetFitFun();
  
  //---------------------------------------
  // Get the histograms from file 
  //---------------------------------------
  
  TH2F  * hRe    [ncombFix];
  TH2F  * hPu    ;
  TH2F  * hMi    [ncombFix];
  
  TString comment[ncombFix];
  TString leg    [ncombFix];
  TString hname  [ncombFix];

  Int_t icalo = 0;
  if(calorimeter=="DCAL") icalo = 1;
  
  if ( !lis )
  {
    hRe[0] = (TH2F *) fil->Get(Form("AnaPi0_Calo%d_hRe_cen0_pidbit0_asy0_dist1",icalo));
    hMi[0] = (TH2F *) fil->Get(Form("AnaPi0_Calo%d_hMi_cen0_pidbit0_asy0_dist1",icalo));
  }
  else
  {
    hRe[0] = (TH2F *) lis->FindObject(Form("AnaPi0_Calo%d_hRe_cen0_pidbit0_asy0_dist1",icalo));
    hMi[0] = (TH2F *) lis->FindObject(Form("AnaPi0_Calo%d_hMi_cen0_pidbit0_asy0_dist1",icalo));
  }

  //hPu     = (TH2F *) fil->Get("AnaPi0_TM1_MCOrgMass_2"                );
  
  printf("histo Re %p - Mix %p\n",hRe[0],hMi[0]);
  
  comment[0] = "All SM";
  leg    [0] = "All SM";
  hname  [0] = "AllSM";
  
  //hRe[0]->Scale(1./nEvt);
  
  // Histograms for cluster pairs in different combinations of SM
  //
  
  // Pairs in same SM
  // 
  if(drawAllCombi)
  {
    Int_t imod;
    Int_t firstMod = 0;
    Int_t lastMod = 11;
    if(calorimeter=="DCAL")
    {
      Int_t firstMod = 12;
      Int_t lastMod  = 19;
    }
    
    for(imod = firstMod; imod<lastMod; imod++)
    {
      printf("imod+1 %d\n",imod+1);
      if(!lis)
      {
        hRe[imod+1] = (TH2F *) fil->Get(Form("AnaPi0_Calo%d_hReMod_%d",icalo,imod));
        hMi[imod+1] = (TH2F *) fil->Get(Form("AnaPi0_Calo%d_hMiMod_%d",icalo,imod));
      }
      else
      {
        hRe[imod+1] = (TH2F *) lis->FindObject(Form("AnaPi0_Calo%d_hReMod_%d",icalo,imod));
        hMi[imod+1] = (TH2F *) lis->FindObject(Form("AnaPi0_Calo%d_hMiMod_%d",icalo,imod));
      }
      
      comment[imod+1] = Form("both clusters in SM %d",imod);
      leg    [imod+1] = Form("SM %d",imod);
      hname  [imod+1] = Form("SM%d" ,imod);
      
      printf("histo mod %d Re %p - Mix %p\n",imod, hRe[imod],hMi[imod]);
    }
    
    // Pairs in same sector
    //
    Int_t isector;
    for(isector = 0; isector<5; isector++)
    {
      printf("imod+1+isector %d\n",imod+1+isector);
 
      if(!lis)
      {
        hRe[imod+1+isector] = (TH2F *) fil->Get(Form("AnaPi0_Calo%d_hReSameSectorEMCALMod_%d",icalo,isector));
        hMi[imod+1+isector] = (TH2F *) fil->Get(Form("AnaPi0_Calo%d_hMiSameSectorEMCALMod_%d",icalo,isector));
      }     
      else 
      {
        hRe[imod+1+isector] = (TH2F *) lis->FindObject(Form("AnaPi0_Calo%d_hReSameSectorEMCALMod_%d",icalo,isector));
        hMi[imod+1+isector] = (TH2F *) lis->FindObject(Form("AnaPi0_Calo%d_hMiSameSectorEMCALMod_%d",icalo,isector));
      }
      
      comment[imod+1+isector] = Form("both clusters in Sector %d",isector);
      leg    [imod+1+isector] = Form("Sector %d",isector);
      hname  [imod+1+isector] = Form("Sector%d" ,isector);
      
      printf("histo sector %d Re %p - Mix %p\n",isector, hRe[imod+1+isector],hMi[imod+1+isector]);
    }
    
    // Pairs in same side
    //
    Int_t iside;
    for(iside = 0; iside < 8; iside++)
    {
      printf("imod+1+isector+iside %d\n",imod+1+isector+iside);
 
      if(!lis)
      {
        hRe[imod+1+isector+iside] = (TH2F *) fil->Get(Form("AnaPi0_Calo%d_hReSameSideEMCALMod_%d",icalo,iside));
        hMi[imod+1+isector+iside] = (TH2F *) fil->Get(Form("AnaPi0_Calo%d_hMiSameSideEMCALMod_%d",icalo,iside));
      }
      else
      {
        hRe[imod+1+isector+iside] = (TH2F *) lis->FindObject(Form("AnaPi0_Calo%d_hReSameSideEMCALMod_%d",icalo,iside));
        hMi[imod+1+isector+iside] = (TH2F *) lis->FindObject(Form("AnaPi0_Calo%d_hMiSameSideEMCALMod_%d",icalo,iside));
      }
      comment[imod+1+isector+iside] = Form("both clusters in Side %d",iside);
      leg    [imod+1+isector+iside] = Form("Side %d",iside);
      hname  [imod+1+isector+iside] = Form("Side%d" ,iside);
      
      printf("histo side %d Re %p - Mix %p\n",isector, hRe[imod+1+isector+iside],hMi[imod+1+isector+iside]);
    }
    
    for(Int_t icomb = 0; icomb<ncomb; icomb++)
      printf("icomb %d p %p ; %p ; %s, %s, %s\n",icomb,hRe[icomb],hMi[icomb],comment[icomb].Data(),leg[icomb].Data(),hname[icomb].Data());
    
  }
  
  //---------------------------------------
  // Init some histograms array size
  //---------------------------------------
  
  TLegend * pLegendIM[nFixBins];

  TH1D* hIM           [ncombFix][nFixBins];
  TH1D* hIMPu         [nFixBins];
  TH1D* hMix          [ncombFix][nFixBins];
  TH1D* hMixCorrected [ncombFix][nFixBins];
  TH1D* hRatio        [ncombFix][nFixBins];
  TH1D* hSignal       [ncombFix][nFixBins];
  
  Double_t mesonPt    [ncombFix][nFixBins];
  Double_t mesonPtBC  [ncombFix][nFixBins];
  Double_t mesonPtBCPu       [nFixBins];
  Double_t mesonMass  [ncombFix][nFixBins];
  Double_t mesonWidth [ncombFix][nFixBins];
  Double_t mesonEffPt [ncombFix][nFixBins];
  Double_t primMeson[nFixBins];

  Double_t emesonPt   [ncombFix][nFixBins];
  Double_t emesonPtBC [ncombFix][nFixBins];
  Double_t emesonPtBCPu      [nFixBins];
  Double_t emesonMass [ncombFix][nFixBins];
  Double_t emesonWidth[ncombFix][nFixBins];
  Double_t emesonEffPt[ncombFix][nFixBins];
  Double_t eprimMeson[nFixBins];

  // Mix correction factor
  //TF1 *line3[ncombFix][nFixBins];// = new TF1("line3", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", 0.0, 1);
  
  
  for(Int_t icomb = 0; icomb< ncomb ; icomb++)
  {
    for(Int_t ipt = 0; ipt< nPt; ipt++)
    {
      hIM          [icomb][ipt] = 0 ;
      hIMPu               [ipt] = 0 ;
      hMix         [icomb][ipt] = 0 ;
      hMixCorrected[icomb][ipt] = 0 ;
      hRatio       [icomb][ipt] = 0 ;
      hSignal      [icomb][ipt] = 0 ;
      
      mesonPt      [icomb][ipt] = 0 ;
      mesonPtBC    [icomb][ipt] = 0 ;
      mesonPtBCPu          [ipt] = 0 ;
      mesonMass    [icomb][ipt] = 0 ;
      mesonWidth   [icomb][ipt] = 0 ;
      mesonEffPt   [icomb][ipt] = 0 ;

      emesonPt     [icomb][ipt] = 0 ;
      emesonPtBC   [icomb][ipt] = 0 ;
      emesonPtBCPu        [ipt] = 0 ;
      emesonMass   [icomb][ipt] = 0 ;
      emesonWidth  [icomb][ipt] = 0 ;
      emesonEffPt  [icomb][ipt] = 0 ;

      //line3        [icomb][ipt] = 0;
    }
  }
  
  //-------------------------------------
  // Do energy projections, fits and plotting
  //-------------------------------------

  
  // Invariant mass
  //
  //gROOT->Macro("$MACROS/style.C");//Set different root style parameters
  gStyle->SetOptTitle(1);
  gStyle->SetTitleOffset(2,"Y");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(000000);
  
  Int_t col = TMath::Ceil(TMath::Sqrt(nPt));

  for(Int_t icomb = 0; icomb< ncomb; icomb++)
  {
    TCanvas  * cIMModi = new TCanvas(Form("Combination_%d\n", icomb), Form("Combination_%d\n", icomb), 1200, 1200) ;
    cIMModi->Divide(col, col);
    
    for(Int_t i = 0; i < nPt; i++)
    {
      if(icomb >10 && i > 8) continue; 

      cIMModi->cd(i+1) ; 
      //gPad->SetLogy();
      Double_t ptMin = xPtLimits[i];
      Double_t ptMax = xPtLimits[i+1];
      printf("Bin %d (%f, %f)\n",i, ptMin,ptMax);
      
      hIM[icomb][i] = hRe[icomb]->ProjectionY(Form("IM_Comb%d_PtBin%d",icomb,i),
                                             hRe[icomb]->GetXaxis()->FindBin(ptMin),
                                             hRe[icomb]->GetXaxis()->FindBin(ptMax));
      hIM[icomb][i]->SetTitle(Form("%2.1f < p_{T, #gamma#gamma} < %2.1f GeV/c",ptMin,ptMax));
      if(i < nPt-3)
        hIM[icomb][i]->Rebin(rebin);
      else
        hIM[icomb][i]->Rebin(rebin2);
            
      if(sumw2)
        hIM[icomb][i]->Sumw2();
      
      hIM[icomb][i]->SetXTitle("M_{#gamma,#gamma} (GeV/c^{2})");
      hIM[icomb][i]->SetLineColor(modColorIndex[icomb]);
      //printf("mmin %f, mmax %f\n",mmin,mmax);
      hIM[icomb][i]->SetAxisRange(mmin,mmax,"X");
      hIM[icomb][i] ->SetLineWidth(2);
      hIM[icomb][i] ->SetLineColor(4);
      //hIM[icomb][i]->Draw();
      
      
//      if(filename.Contains("imu") && icomb==0)
//      {
//        hIMPu[i] = hPu->ProjectionY(Form("IMPu_Comb%d_PtBin%d",icomb,i),
//                                    hPu->GetXaxis()->FindBin(ptMin),
//                                    hPu->GetXaxis()->FindBin(ptMax));
//        hIMPu[i]->SetTitle(Form("%2.1f < p_{T, #gamma#gamma} < %2.1f GeV/c",ptMin,ptMax));
//        if(i < nPt-3)
//          hIMPu[i]->Rebin(rebin);
//        else
//          hIMPu[i]->Rebin(rebin2);
//        
//        if(sumw2)
//          hIMPu[icomb][i]->Sumw2();
//        else
//          hIMPu[icomb][i]->Scale(1.e8);
//        
//        hIMPu[i]->SetXTitle("M_{#gamma,#gamma} (GeV/c^{2})");
//        hIMPu[i]->SetLineColor(modColorIndex[icomb]);
//        //printf("mmin %f, mmax %f\n",mmin,mmax);
//        hIMPu[i]->SetAxisRange(mmin,mmax,"X");
//        hIMPu[i]->SetLineWidth(2);
//        hIMPu[i]->SetLineColor(kOrange-3);
//      }
      
      Double_t nMax     = 0;
      if(part=="Pi0")
      {
        nMax= hIM[icomb][i]->Integral(hIM[icomb][i]->FindBin(0.09),
                                      hIM[icomb][i]->FindBin(0.25));
      }
      else
      {
        nMax = hIM[icomb][i]->Integral(hIM[icomb][i]->FindBin(0.4),
                                       hIM[icomb][i]->FindBin(0.65));
      }
      
      //printf("icombi %d, bin %d, nMax %f \n",icomb,i,nMax);
      
      fitfun->SetParLimits(0,nMax/100,nMax*100);

      if(nMax > nPairCut && !mix)
      {
        if(part=="Pi0")
        {
          fitfun->SetParameters(nMax/4,0.135,0.02,1.,0.);
          if(i<9)
            hIM[icomb][i]->Fit("fitfun","QR","",0.09,0.25);
          else
            hIM[icomb][i]->Fit("fitfun","QR","",0.10,0.25);
        } // pi0
        else  //eta
        {
          fitfun->SetParameters(nMax/4,0.54,0.04,1.,0.);
          if(i<4)
            hIM[icomb][i]->Fit("fitfun","QR","",0.4,0.7);
          else if(i < 15)
            hIM[icomb][i]->Fit("fitfun","QR","",0.3,0.8);
          else 
            hIM[icomb][i]->Fit("fitfun","QR","",0.3,0.8);
        }
      }
      
      if(part=="Pi0")
      {
        if(i < nPt -3)
          hIM[icomb][i]->SetMaximum(nMax/5.);
        else
          hIM[icomb][i]->SetMaximum(nMax/2.);
        
      }
      else
      {
        if(i < nPt -3)
          hIM[icomb][i]->SetMaximum(nMax/20.);
        else
          hIM[icomb][i]->SetMaximum(nMax/8.);
        
      }
      
      Double_t mStep =  hIM[icomb][i]->GetBinWidth(i);
      hIM[icomb][i]->SetMinimum(0.1);
      hIM[icomb][i]->SetAxisRange(mmin,mmax,"X");
      
      pLegendIM[i] = 0;
      
      if ( hIM[icomb][i]->GetFunction("fitfun") && !mix )
      {
        if(hIM[icomb][i]->GetFunction("fitfun")->GetChisquare()/hIM[icomb][i]->GetFunction("fitfun")->GetNDF()<100)
        {
          Double_t A    = hIM[icomb][i]->GetFunction("fitfun")->GetParameter(0);
          Double_t mean = hIM[icomb][i]->GetFunction("fitfun")->GetParameter(1);
          Double_t sigm = hIM[icomb][i]->GetFunction("fitfun")->GetParameter(2);
//          Double_t a0   = hIM[icomb][i]->GetFunction("fitfun")->GetParameter(3);
//          Double_t a1   = hIM[icomb][i]->GetFunction("fitfun")->GetParameter(4);
//          Double_t a2   = 0;
//          if (polN == 2)
//            a2 = hIM[icomb][i]->GetFunction("fitfun")->GetParameter(5);
          
          Double_t eA    = hIM[icomb][i]->GetFunction("fitfun")->GetParError(0);
          Double_t emean = hIM[icomb][i]->GetFunction("fitfun")->GetParError(1);
          Double_t esigm = hIM[icomb][i]->GetFunction("fitfun")->GetParError(2);
//          Double_t ea0   = hIM[icomb][i]->GetFunction("fitfun")->GetParError(3);
//          Double_t ea1   = hIM[icomb][i]->GetFunction("fitfun")->GetParError(4);
//          Double_t ea2   = 0;   
//          if (polN == 2)
//            ea2 = hIM[icomb][i]->GetFunction("fitfun")->GetParError(5);
          
          pLegendIM[i] = new TLegend(0.55,0.65,0.95,0.93);
          pLegendIM[i]->SetTextSize(0.035);
          pLegendIM[i]->SetFillColor(10);
          pLegendIM[i]->SetBorderSize(1);
          pLegendIM[i]->SetHeader(Form("#pi %s",leg[icomb].Data()));
          pLegendIM[i]->AddEntry("",Form("A = %2.2f #pm %2.2f ",A,eA),"");
          pLegendIM[i]->AddEntry("",Form("#mu = %2.3f #pm %2.3f ",mean,emean),"");
          pLegendIM[i]->AddEntry("",Form("#sigma = %2.3f #pm %2.3f ",sigm,esigm),"");
          //pLegendIM[i]->AddEntry("",Form("p_{0} = %2.1f #pm %2.1f ",a0,ea0),"");
          //pLegendIM[i]->AddEntry("",Form("p_{1} = %2.1f #pm %2.1f ",a1,ea1),"");
          //if(polN==2)pLegendIM[i]->AddEntry("",Form("p_{2} = %2.1f#pm %2.1f ",a2,ea2),"");
          
          mesonPt[icomb] [i] = A*sigm / mStep * TMath::Sqrt(TMath::TwoPi());
          Double_t ePi0 =
          TMath::Power(eA/A,2) +
          TMath::Power(esigm/sigm,2);
          emesonPt[icomb][i]= TMath::Sqrt(ePi0)* mesonPt[icomb][i];
          
          emesonPt[icomb][i] = TMath::Min(emesonPt[icomb][i], TMath::Sqrt(mesonPt[icomb] [i]));
          
          mesonPt [icomb] [i]/=(nEvt*(exPt[i]*2));
          emesonPt[icomb] [i]/=(nEvt*(exPt[i]*2));
          
          cout << "N(pi0) fit  = " <<  mesonPt[icomb][i] << "+-"<<  emesonPt[icomb][i] << " mstep " << mStep<<endl;
          
          mesonMass [icomb][i] = mean*1000.;
          emesonMass[icomb][i] = emean*1000.;
          
          mesonWidth [icomb][i] = sigm*1000.;
          emesonWidth[icomb][i] = esigm*1000.;
          
        } //Good fit
        else{
          
          mesonPt[icomb][i] =-1;
          emesonPt[icomb][i] = -1;
          
          mesonMass[icomb][i] = -1;
          emesonMass[icomb][i] = -1;
          
          mesonWidth[icomb][i] = -1;
          emesonWidth[icomb][i] = -1;
        }
      }
      else{
        mesonPt[icomb][i] =-1;
        emesonPt[icomb][i] = -1;
        
        mesonMass[icomb][i] = -1;
        emesonMass[icomb][i] = -1;
        
        mesonWidth[icomb][i] = -1;
        emesonWidth[icomb][i] = -1;
      }
      

      //--------------------------------------------------
      //   Mix
      //--------------------------------------------------
      if(mix)
      {
        hMix[icomb][i] = (TH1D*) hMi[icomb]->ProjectionY(Form("MiMass_PtBin%d_Combi%d",i,icomb),
                                                         hMi[icomb]->GetXaxis()->FindBin(ptMin), 
                                                         hMi[icomb]->GetXaxis()->FindBin(ptMax));
        hMix[icomb][i]->SetLineColor(1);  
        hMix[icomb][i]->SetTitle(Form("%2.2f < p_{T, #gamma#gamma} < %2.2f GeV/c",ptMin,ptMax));
        if(i < nPt-3)
          hMix[icomb][i]->Rebin(rebin);
        else
          hMix[icomb][i]->Rebin(rebin2);
        

        
        //         Double_t sptmin = 0.2;
        //         Double_t sptmax = 0.45;
        
        //         Double_t scalereal = hIM[icomb]->Integral(hIM[icomb]->FindBin(sptmin),
        //                                                             hIM[icomb]->FindBin(sptmax));
        //         Double_t scalemix  = hMix[icomb]->Integral(hMix[icomb]->FindBin(sptmin),
        //                                                             hMix[icomb]->FindBin(sptmax));
        //         if(scalemix > 0)
        //           hMix[icomb]->Scale(scalereal/scalemix);
        //         else{ 
        //           //printf("could not scale:  re %f  mi %f\n",scalereal,scalemix);
        //           continue;
        //         }
        
        //      hMix[icomb][i]->SetLineWidth(1);  

        if(!sumw2) 
          hMix[icomb][i]->Sumw2();

        // --- Ratio ---  
        hRatio[icomb][i] = (TH1D*)hIM[icomb][i]->Clone(Form("RatioRealMix_PtBin%d_comb%d",i,icomb));
        hRatio[icomb][i]->SetAxisRange(mmin,mmax,"X");
        hRatio[icomb][i]->Divide(hMix[icomb][i]);
        Double_t p0 =0;
        Double_t p1 =0;
        Double_t p2 =0;
        Double_t p3 =0;

        if(hRatio[icomb][i]->GetEntries()>100)
        {
//          if(part=="Pi0")
//          {
//            hRatio[icomb][i]->Fit("truncatedPolPi0", "NQRL","",0.11,0.4);
//            p0=  tPolPi0->GetParameter(0);           
//            p1=  tPolPi0->GetParameter(1);
//            p2=  tPolPi0->GetParameter(2);
//            p3=  tPolPi0->GetParameter(3);
//          } // pi0
//          else  // eta
//          {
//            hRatio[icomb][i]->SetAxisRange(mmin,mmax);
//            hRatio[icomb][i]->Fit("truncatedPolEta", "NQRL","",0.25,0.95);
//            p0 =  tPolEta->GetParameter(0);      
//            p1 =  tPolEta->GetParameter(1);
//            p2 =  tPolEta->GetParameter(2);
//            p3 =  tPolEta->GetParameter(3);
//          }
//
//          //hRatio[icomb][i]->Draw("same");
//          
//          //line3[icomb][i]= new TF1("line3", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", mmin, mmax);
//          line3[icomb][i]= new TF1("line3", "[0]+[1]*x", mmin, mmax);
//          line3[icomb][i]->SetParameters(p0,p1,p2,p3);
//
//          if(hMix[icomb][i])
//          {
//            //printf("Correct\n");
//            hMixCorrected[icomb][i] = (TH1D*)hMix[icomb][i]->Clone(Form("MixCorrected_PtBin%d_comb%d",i,icomb));
//            for(Int_t j = 0; j< hMix[icomb][i]->GetNbinsX(); j++){
//              Double_t x          = hMix[icomb][i]->GetBinCenter(j);
//              Double_t correction = line3[icomb][i]->Eval(x);
//              Double_t corrected  = hMix[icomb][i]->GetBinContent(j)*correction;
//              Double_t ecorrected = hMix[icomb][i]->GetBinError(j)  *correction;
//              //printf("bin %d, center %f, mix %f, corrected %f\n",i,x,hMix[icomb][i]->GetBinContent(j),corrected);
//              hMixCorrected[icomb][i]->SetBinContent(j, corrected);
//              hMixCorrected[icomb][i]->SetBinError  (j,ecorrected);
//            }
//          }

          Float_t scale = hRatio[icomb][i]->GetBinContent(hRatio[icomb][i]->FindBin(0.9));
          hMix[icomb][i]->Scale(scale);
          
          // --- Subtracted ---
          //printf("Signal\n");
          hSignal[icomb][i] = (TH1D*)hIM[icomb][i]->Clone(Form("Signal_PtBin%d_comb%d",i,icomb));
          //hSignal[icomb][i] ->Add(hMixCorrected[icomb][i],-1);
          hSignal[icomb][i] ->Add(hMix[icomb][i],-1);
          hSignal[icomb][i] ->SetLineColor(kViolet);
         
          // ---- Fit subtracted signal
          if(part=="Pi0")
          {
            nMax= hSignal[icomb][i]->Integral(hIM[icomb][i]->FindBin(0.1),
                                              hIM[icomb][i]->FindBin(0.2));
          }
          else
          {
            nMax = hSignal[icomb][i]->Integral(hIM[icomb][i]->FindBin(0.4),
                                               hIM[icomb][i]->FindBin(0.65));
          }
          
          if(nMax > nPairCut)
          {
            fitfun->SetParLimits(0,nMax/100,nMax*100);
            
            if(part=="Pi0"){
              if(polN < 4 )fitfun->SetParameters(nMax/5,0.135,20,0);
              else         fitfun->SetParameters(nMax/5,20,0.135,20,20,nMax/5);

              if(i < nPt-4)
                hSignal[icomb][i]->Fit("fitfun","QR","",0.11,0.3);
              else
                hSignal[icomb][i]->Fit("fitfun","QR","",0.11,0.3);
            }// pi0
            else // eta
            {
              
              if(polN < 4 )fitfun->SetParameters(nMax/5,0.54,40,0);
              else         fitfun->SetParameters(nMax/5,40,0.54,40,40,nMax/5);
              //if(i<4)
              hSignal[icomb][i]->Fit("fitfun","QR","",0.4,0.7);
              //else if(i < 15)
              //  hSignal[icomb][i]->Fit("fitfun","QRL","",0.3,0.8);
              //else 
              //hSignal[icomb][i]->Fit("fitfun","QRL","",0.3,0.8);
            }
          }
          
          if(hSignal[icomb][i]->GetFunction("fitfun"))
          {
//            printf("ipt %d: Chi/NDF %f, Chi %f, NDF %d\n",i,
//                   hSignal[icomb][i]->GetFunction("fitfun")->GetChisquare()/hSignal[icomb][i]->GetFunction("fitfun")->GetNDF(),
//                   hSignal[icomb][i]->GetFunction("fitfun")->GetChisquare(),
//                   hSignal[icomb][i]->GetFunction("fitfun")->GetNDF());
            if(hSignal[icomb][i]->GetFunction("fitfun")->GetChisquare()/hSignal[icomb][i]->GetFunction("fitfun")->GetNDF()<1000){
              Double_t A    = hSignal[icomb][i]->GetFunction("fitfun")->GetParameter(0);
              Double_t mean = hSignal[icomb][i]->GetFunction("fitfun")->GetParameter(1);
              Double_t sigm = hSignal[icomb][i]->GetFunction("fitfun")->GetParameter(2);
//              Double_t a0   = hSignal[icomb][i]->GetFunction("fitfun")->GetParameter(3);
//              Double_t a1   = hSignal[icomb][i]->GetFunction("fitfun")->GetParameter(4);
//              Double_t a2   = 0;
//              if (polN == 2)
//                a2 = hSignal[icomb][i]->GetFunction("fitfun")->GetParameter(5);
              
              Double_t eA    = hSignal[icomb][i]->GetFunction("fitfun")->GetParError(0);
              Double_t emean = hSignal[icomb][i]->GetFunction("fitfun")->GetParError(1);
              Double_t esigm = hSignal[icomb][i]->GetFunction("fitfun")->GetParError(2);
//              Double_t ea0   = hSignal[icomb][i]->GetFunction("fitfun")->GetParError(3);
//              Double_t ea1   = hSignal[icomb][i]->GetFunction("fitfun")->GetParError(4);
//              Double_t ea2   = 0;
//              if (polN == 2)
//                ea2 = hSignal[icomb][i]->GetFunction("fitfun")->GetParError(5);
              pLegendIM[i] = new TLegend(0.55,0.65,0.95,0.93);
              pLegendIM[i]->SetTextSize(0.035);
              pLegendIM[i]->SetFillColor(10);
              pLegendIM[i]->SetBorderSize(1);
              pLegendIM[i]->SetHeader(Form("#pi %s",leg[icomb].Data()));
              pLegendIM[i]->AddEntry("",Form("A = %2.2f #pm %2.2f ",A,eA),"");
              pLegendIM[i]->AddEntry("",Form("#mu = %2.3f #pm %2.3f ",mean,emean),"");
              pLegendIM[i]->AddEntry("",Form("#sigma = %2.3f #pm %2.3f ",sigm,esigm),"");
              //pLegendIM[i]->AddEntry("",Form("p_{0} = %2.1f #pm %2.1f ",a0,ea0),"");
              //pLegendIM[i]->AddEntry("",Form("p_{1} = %2.1f #pm %2.1f ",a1,ea1),"");
              //if(polN==2)pLegendIM[i]->AddEntry("",Form("p_{2} = %2.1f#pm %2.1f ",a2,ea2),"");
              
              
              mesonPt[icomb] [i] = A*sigm / mStep * TMath::Sqrt(TMath::TwoPi());
              Double_t ePi0 =
              TMath::Power(eA/A,2) +
              TMath::Power(esigm/sigm,2);
              emesonPt[icomb][i]= TMath::Sqrt(ePi0)* mesonPt[icomb][i];
              
              emesonPt[icomb][i] = TMath::Min(emesonPt[icomb][i], TMath::Sqrt(mesonPt[icomb] [i]));
              
              mesonPt [icomb] [i]/=(nEvt*(exPt[i]*2));
              emesonPt[icomb] [i]/=(nEvt*(exPt[i]*2));
              
              //cout << "N(pi0) fit  = " <<  mesonPt[icomb][i] << "+-"<<  emesonPt[icomb][i] << " mstep "<<mStep<<endl;
              
              mesonMass [icomb][i] = mean *1000.;
              emesonMass[icomb][i] = emean*1000.;
              
              mesonWidth [icomb][i] = sigm *1000.;
              emesonWidth[icomb][i] = esigm*1000.;
              
            } //Good fit
            else{
              
              mesonPt[icomb][i] =-1;
              emesonPt[icomb][i] = -1;
              
              mesonMass[icomb][i] = -1;
              emesonMass[icomb][i] = -1;
              
              mesonWidth[icomb][i] = -1;
              emesonWidth[icomb][i] = -1;
            }
          }
          else{
            mesonPt[icomb][i] =-1;
            emesonPt[icomb][i] = -1;
            
            mesonMass[icomb][i] = -1;
            emesonMass[icomb][i] = -1;
            
            mesonWidth[icomb][i] = -1;
            emesonWidth[icomb][i] = -1;
          }
          
          Double_t mass  = mesonMass [icomb][i];
          Double_t width = mesonWidth[icomb][i];
          Double_t mMinBin =   mass - 2 * width;
          Double_t mMaxBin =   mass + 2 * width;
          if(mass <= 0 || width <= 0)
          {
            if(part=="Pi0")
            {
              mMinBin = 0.115; mMaxBin = 0.3 ;
            }
            else{
              mMinBin = 0.4;  mMaxBin = 0.9 ; 
            }
            
          }
          
          if(hSignal[icomb][i])
          {
            mesonPtBC[icomb] [i] = hSignal[icomb][i]->Integral(hSignal[icomb][i]->FindBin(mMinBin),  hSignal[icomb][i]->FindBin(mMaxBin)) ;
          
            emesonPtBC[icomb] [i] = TMath::Sqrt(mesonPtBC[icomb] [i]);
          
            mesonPtBC[icomb]  [i]/=(nEvt*(exPt[i]*2));
            emesonPtBC[icomb] [i]/=(nEvt*(exPt[i]*2));
             //printf("Bin %d, [%1.1f,%1.1f]\n",i,xPt[i],xPt[i+1]);
          }
          
//          if(filename.Contains("imu") && icomb==0)
//          {
//            mesonPtBCPu [i] = hIMPu[i]->Integral(hIMPu[i]->FindBin(mMinBin),  hIMPu[i]->FindBin(mMaxBin)) ;
//            
//            emesonPtBCPu[i] = TMath::Sqrt(mesonPtBCPu [i]);
//            
//            mesonPtBCPu [i]/=(nEvt*(exPt[i]*2));
//            emesonPtBCPu[i]/=(nEvt*(exPt[i]*2));
//            //printf("Bin %d, [%1.1f,%1.1f]\n",i,xPt[i],xPt[i+1]);
//
//          }
          
          printf("N(pi0) BC = %2.3e +- %2.4e\n",mesonPtBC[icomb][i],emesonPtBC[icomb][i]);
          if(filename.Contains("imu") && icomb==0)
          {
            printf("N(pi0) Pu = %2.3e +- %2.4e\n",mesonPtBCPu[i],emesonPtBCPu[i]);
            if(mesonPtBCPu[i]>0)printf("\t ratio BC / Pu %f \n",mesonPtBC[icomb][i]/mesonPtBCPu[i]);
          }
        }

        
        if(nMax > nPairCut)
        {
          hIM[icomb][i]->Draw("HE");

          hMix[icomb][i]->Draw("same");
          //if(hMixCorrected[icomb][i])hMixCorrected[icomb][i] ->Draw("same");

          if(hSignal[icomb][i]) hSignal[icomb][i] -> Draw("same");

          if(hSignal[icomb][i] && hSignal[icomb][i]->GetFunction("fitfun") && pLegendIM[i]) pLegendIM[i]->Draw();
        }

      }//mix on
      else if(nMax > nPairCut)
      {
        //printf("HERE 1\n");

        hIM[icomb][i]->Draw("HE");
        //printf("HERE 2\n");

        //if(hIM[icomb][i]->GetFunction("fitfun") && pLegendIM[i])pLegendIM[i]->Draw();
        //printf("HERE 3\n");

      }
    } //pT bins
    
    cIMModi->Print(Form("IMfigures/%s_%s_Mgg_%s_%s_%s.%s",prodname.Data(),calorimeter.Data(),part.Data(),filename.Data(), hname[icomb].Data(),fileFormat.Data()));
    
  }// combinations
    
  
  //xxxx Real / Mixxxxx
  if(mix)
  {
    for(Int_t icomb = 0; icomb< ncomb; icomb++)
    {
      printf("Do real/mix\n");

      TCanvas  * cRat = new TCanvas(Form("CombinationRatio_%d\n", icomb), Form("CombinationRatio_%d\n", icomb), 1200, 1200) ;
      cRat->Divide(col, col);
     
      for(Int_t i = 0; i < nPt; i++){
        cRat->cd(i+1) ;
        //gPad->SetGridy();
        //gPad->SetLog();
        if(!hRatio[icomb][i]) {
          //printf("No ratio in pt bin %d continue\n",i);
          continue;
        }
                
        Double_t nMax =0;
        if(part=="Pi0")
          nMax = hRatio[icomb][i]->Integral(hRatio[icomb][i]->FindBin(0.005),
                                            hRatio[icomb][i]->FindBin(0.20));
        else
          nMax = hRatio[icomb][i]->Integral(hRatio[icomb][i]->FindBin(0.4),
                                            hRatio[icomb][i]->FindBin(0.6));
        
        hRatio[icomb][i]->SetMaximum(nMax/4);
        hRatio[icomb][i]->SetMinimum(1e-6);
        hRatio[icomb][i]->SetAxisRange(mmin,mmax,"X");

        hRatio[icomb][i]->Draw();
        
        //if(line3[icomb][i]) line3[icomb][i]->Draw("same");
      }
      
      
     cRat->Print(Form("IMfigures/%s_%s_MggRatio_%s_%s_%s.%s",prodname.Data(),calorimeter.Data(),part.Data(),filename.Data(), hname[icomb].Data(),fileFormat.Data()));

    }
  }    
  
  //------------------------------
  // Fit parameters final plotting
  //------------------------------
  
  gStyle->SetOptTitle(0);
  
  // Put fit results in TGraphErrors
  //
  TGraphErrors* gPt   [ncombFix];
  TGraphErrors* gPtBC [ncombFix];
  TGraphErrors* gMass [ncombFix];
  TGraphErrors* gWidth[ncombFix];
  
  TString particleN = " #pi^{0}";
  if(part=="Eta") particleN = " #eta";
  
  for(Int_t icomb = 0; icomb<ncomb; icomb++)
  {
    //cout<<"icomb "<<icomb<<" " <<gPt[icomb]<<" "<<comment[icomb]<<endl;
    gPt[icomb] =  new TGraphErrors(nPt,xPt,mesonPt[icomb],exPt,emesonPt[icomb]);
    gPt[icomb] ->SetName(Form("gPt_%s",hname[icomb].Data()));
    gPt[icomb] ->GetHistogram()->SetTitle(Form("p_{T} of reconstructed %s, %s ",particleN.Data(), comment[icomb].Data()));
    gPt[icomb] ->GetHistogram()->SetXTitle("p_{T} (GeV/c)    ");
    gPt[icomb] ->GetHistogram()->SetYTitle(Form("dN_{%s}/dp_{T} (GeV/c)^{-1}  / N_{events}  ",particleN.Data()));
    gPt[icomb] ->GetHistogram()->SetAxisRange(0.,30);
    gPt[icomb] ->GetHistogram()->SetTitleOffset(1.5,"Y");
    gPt[icomb] ->SetMarkerStyle(20);
    gPt[icomb] ->SetMarkerSize(1);
    gPt[icomb] ->SetMarkerColor(1);
    //   gPt[icomb] ->GetHistogram()->SetMaximum(1e8);
    //   gPt[icomb] ->GetHistogram()->SetMinimum(1e-8);
    
    gPtBC[icomb] =  new TGraphErrors(nPt,xPt,mesonPtBC[icomb],exPt,emesonPtBC[icomb]);
    gPtBC[icomb] ->SetName(Form("gPtBC_%s",hname[icomb].Data()));
    gPtBC[icomb] ->GetHistogram()->SetTitle(Form("p_{T} of reconstructed %s, %s ",particleN.Data(), comment[icomb].Data()));
    gPtBC[icomb] ->GetHistogram()->SetXTitle("p_{T} (GeV/c)    ");
    gPtBC[icomb] ->GetHistogram()->SetYTitle(Form("dN_{%s}/dp_{T} (GeV/c)^{-1}  / N_{events}  ",particleN.Data()));
    gPtBC[icomb] ->GetHistogram()->SetAxisRange(0.,30);
    gPtBC[icomb] ->GetHistogram()->SetTitleOffset(1.5,"Y");
    gPtBC[icomb] ->SetMarkerStyle(20);
    gPtBC[icomb] ->SetMarkerSize(1);
    gPtBC[icomb] ->SetMarkerColor(1);
    
    
    //   gPtBC[icomb] ->GetHistogram()->SetMaximum(1e8);
    //   gPtBC[icomb] ->GetHistogram()->SetMinimum(1e-8);
    
    gMass[icomb] =  new TGraphErrors(nPt,xPt,mesonMass[icomb],exPt,emesonMass[icomb]);
    gMass[icomb] ->SetName(Form("gMass_%s",hname[icomb].Data()));
    gMass[icomb] ->GetHistogram()->SetTitle(Form("mass vs p_{T} of reconstructed %s, %s ",particleN.Data(), comment[icomb].Data()));
    gMass[icomb] ->GetHistogram()->SetXTitle("p_{T} (GeV/c)    ");
    gMass[icomb] ->GetHistogram()->SetYTitle(Form("mass_{%s} (MeV/c)^{2}    ",particleN.Data()));
    //gMass[icomb] ->GetHistogram()->SetAxisRange(0.,30);
    gMass[icomb] ->GetHistogram()->SetTitleOffset(1.5,"Y");
    gMass[icomb] ->SetMarkerStyle(20);
    gMass[icomb] ->SetMarkerSize(1);
    gMass[icomb] ->SetMarkerColor(1);
    
    
    gWidth[icomb] =  new TGraphErrors(nPt,xPt,mesonWidth[icomb],exPt,emesonWidth[icomb]);
    gWidth[icomb] ->SetName(Form("gWidth_%s",hname[icomb].Data()));
    gWidth[icomb] ->GetHistogram()->SetTitle(Form("Width vs p_{T} of reconstructed %s, %s ",particleN.Data(), comment[icomb].Data()));
    gWidth[icomb] ->GetHistogram()->SetXTitle("p_{T} (GeV/c)    ");
    gWidth[icomb] ->GetHistogram()->SetYTitle(Form("#sigma_{%s} (MeV/c)^{2}    ",particleN.Data()));
    //gWidth[icomb] ->GetHistogram()->SetAxisRange(0.,30);
    gWidth[icomb] ->GetHistogram()->SetTitleOffset(1.5,"Y");
    gWidth[icomb] ->SetMarkerStyle(20);
    gWidth[icomb] ->SetMarkerSize(1);
    gWidth[icomb] ->SetMarkerColor(1);
    
  }
  
  TGraphErrors *gPtBCPu = 0;
  if(filename.Contains("imu"))
  {
    gPtBCPu =  new TGraphErrors(nPt,xPt,mesonPtBCPu,exPt,emesonPtBCPu);
    gPtBCPu ->SetName(Form("gPtBC_%s_Pure",hname[0].Data()));
    gPtBCPu ->GetHistogram()->SetTitle(Form("p_{T} of reconstructed %s, %s ",particleN.Data(), comment[0].Data()));
    gPtBCPu ->GetHistogram()->SetXTitle("p_{T} (GeV/c)    ");
    gPtBCPu ->GetHistogram()->SetYTitle(Form("dN_{%s}/dp_{T} (GeV/c)^{-1}  / N_{events}  ",particleN.Data()));
    //gPtBCPu ->GetHistogram()->SetAxisRange(0.,30);
    gPtBCPu ->GetHistogram()->SetTitleOffset(1.6,"Y");
    gPtBCPu ->SetMarkerStyle(24);
    gPtBCPu ->SetMarkerSize(1);
    gPtBCPu ->SetMarkerColor(2);
  }
  
  
  //PLOTS
  //
  
  //SM
  TLegend pLegendFitM(0.15,0.53,0.3,0.93);
  pLegendFitM.SetTextSize(0.03);
  pLegendFitM.SetFillColor(10);
  pLegendFitM.SetBorderSize(1);
  
  TLegend pLegendFitW(0.4,0.53,0.6,0.93);
  pLegendFitW.SetTextSize(0.03);
  pLegendFitW.SetFillColor(10);
  pLegendFitW.SetBorderSize(1);
  
  TCanvas *cFitSM = new TCanvas("cFitSM","cFitSM",1200,600);
  cFitSM->Divide(2, 1);
  
  cFitSM->cd(1);
  //......................................................................
  gPad->SetGridx();
  gPad->SetGridy();
  
  if(part=="Pi0")
  {
    gMass[0]->SetMaximum(180);
    gMass[0]->SetMinimum(120);
  }
  else
  {
    gMass[0]->SetMaximum(660);
    gMass[0]->SetMinimum(420);
  }
  
  gMass[0]->SetMarkerStyle(modStyleIndex[0]);
  gMass[0]->SetMarkerColor(modColorIndex[0]);
  gMass[0]->SetLineColor(modColorIndex[0]);
  
  gMass[0]->Draw("AP");
  
  pLegendFitM.AddEntry(gMass[0],Form("%s",leg[0].Data()),"P");
  
  if(drawAllCombi)
  {
    for(Int_t icomb = 1; icomb <11; icomb++)
    {
      //printf("option %d\n",icomb);
      pLegendFitM.AddEntry(gMass[icomb],Form("%s",leg[icomb].Data()),"P");
      gMass[icomb]->SetMarkerStyle(modStyleIndex[icomb]);
      gMass[icomb]->SetMarkerColor(modColorIndex[icomb]);
      gMass[icomb]->SetLineColor(modColorIndex[icomb]);
      gMass[icomb]->Draw("P");
    }
  } 
  
  gMass[0]->Draw("P");
  if(ncomb > 1) pLegendFitM.Draw();
  
  cFitSM->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
 
  if(part=="Pi0")
  {
    gWidth[0]->SetMaximum(16);
    gWidth[0]->SetMinimum(6);
  }
  else
  {
    gWidth[0]->SetMaximum(60);
    gWidth[0]->SetMinimum(0);
  }
  
  gWidth[0]->SetMarkerStyle(modStyleIndex[0]);
  gWidth[0]->SetMarkerColor(modColorIndex[0]);
  gWidth[0]->SetLineColor  (modColorIndex[0]);  
  
  gWidth[0]->Draw("AP");   
  pLegendFitW.AddEntry(gWidth[0],Form("%s",leg[0].Data()),"P");
  
  if(drawAllCombi)
  {
    for(Int_t icomb = 1; icomb <11; icomb++)
    {
      //printf("option %d\n",icomb);
      pLegendFitW.AddEntry(gWidth[icomb],Form("%s",leg[icomb].Data()),"P");
      gWidth[icomb]->SetMarkerStyle(modStyleIndex[icomb]);
      gWidth[icomb]->SetMarkerColor(modColorIndex[icomb]);
      gWidth[icomb]->SetLineColor(modColorIndex[icomb]);
      gWidth[icomb]->Draw("P");
      
    }
  }
  
  gWidth[0]->Draw("P"); 
  if(ncomb > 1) pLegendFitW.Draw();
  
  cFitSM->Print(Form("IMfigures/%s_%s_SingleSM_MassWidth_%s_%s.%s",prodname.Data(),calorimeter.Data(),part.Data(),filename.Data(),fileFormat.Data()));
  
  
//  if(drawAllCombi)
//  {
//    //Pi0
//    //Sector
//    TLegend pLegendFitMS(0.75,0.15,0.95,0.45);
//    pLegendFitMS.SetTextSize(0.03);
//    pLegendFitMS.SetFillColor(10);
//    pLegendFitMS.SetBorderSize(1);
//    
//    
//    TLegend pLegendFitWS(0.7,0.6,0.9,0.9);
//    pLegendFitWS.SetTextSize(0.03);
//    pLegendFitWS.SetFillColor(10);
//    pLegendFitWS.SetBorderSize(1);
//    
//    TCanvas *cFitS = new TCanvas("cFitS","cFitS",1200,600);
//    cFitS->Divide(2, 1);
//    
//    cFitS->cd(1);
//    //  //......................................................................
//    //  TCanvas *c5 = new TCanvas("c5","c5",0,0,600,600);
//    gPad->SetGridx();
//    gPad->SetGridy();

//    gMass[0]->SetMaximum(0.15);
//    
//    //printf("Fitting Mass plot\n");
//    gMass[0]->Draw("AP");
//    pLegendFitMS.AddEntry(gMass[0],Form("%s",leg[0].Data()),"P");
//    
//    for(Int_t icomb = 11; icomb <16; icomb++){
//      //printf("option %d\n",icomb);
//      pLegendFitMS.AddEntry(gMass[icomb],Form("%s",leg[icomb].Data()),"P");
//      gMass[icomb]->SetMarkerStyle(modStyleIndex[icomb]);
//      gMass[icomb]->SetMarkerColor(modColorIndex[icomb]);
//      gMass[icomb]->SetLineColor(modColorIndex[icomb]);
//      gMass[icomb]->Draw("P");
//    }
//    pLegendFitMS.Draw();
//    
//    cFitS->cd(2);
//    gPad->SetGridx();
//    gPad->SetGridy();
//    gWidth[0]->Draw("AP");
//    //printf("Width plots\n");
//    pLegendFitWS.AddEntry(gWidth[0],Form("%s",leg[0].Data()),"P");
//    
//    for(Int_t icomb = 11; icomb <16; icomb++){
//      //printf("option %d\n",icomb);
//      pLegendFitWS.AddEntry(gWidth[icomb],Form("%s",leg[icomb].Data()),"P");
//      
//      gWidth[icomb]->SetMarkerStyle(modStyleIndex[icomb]);
//      gWidth[icomb]->SetMarkerColor(modColorIndex[icomb]);
//      gWidth[icomb]->SetLineColor(modColorIndex[icomb]);
//      gWidth[icomb]->Draw("P");
//      
//    }
//    pLegendFitWS.Draw();
//    
//    cFitS->Print(Form("IMfigures/%s_%s_Sector_MassWidth_%s_%s.%s",prodname.Data(),calorimeter.Data(),part.Data(),filename.Data(),fileFormat.Data()));
//
//    //Side
//    TLegend pLegendFitMSi(0.7,0.15,0.9,0.45);
//    pLegendFitMSi.SetTextSize(0.03);
//    pLegendFitMSi.SetFillColor(10);
//    pLegendFitMSi.SetBorderSize(1);
//    
//    
//    TLegend pLegendFitWSi(0.7,0.5,0.9,0.9);
//    pLegendFitWSi.SetTextSize(0.03);
//    pLegendFitWSi.SetFillColor(10);
//    pLegendFitWSi.SetBorderSize(1);
//    
//    TCanvas *cFitSi = new TCanvas("cFitSi","cFitSi",1200,600);
//    cFitSi->Divide(2, 1);
//    //printf("Side \n");
//    cFitSi->cd(1);
//    //  //......................................................................
//    //  TCanvas *c5 = new TCanvas("c5","c5",0,0,600,600);
//    gPad->SetGridx();
//    gPad->SetGridy();
//    //  c5->SetLeftMargin(0.12);
//    //  c5->SetRightMargin(0.08);
//    //printf("Fitting Mass plot\n");
//    gMass[0]->Draw("AP");
//    pLegendFitMSi.AddEntry(gMass[0],Form("%s",leg[0].Data()),"P");
//    
//    for(Int_t icomb = 16; icomb <24; icomb++){
//      //printf("option %d\n",icomb);
//      pLegendFitMSi.AddEntry(gMass[icomb],Form("%s",leg[icomb].Data()),"P");
//      gMass[icomb]->SetMarkerStyle(modStyleIndex[icomb]);
//      gMass[icomb]->SetMarkerColor(modColorIndex[icomb]);
//      gMass[icomb]->SetLineColor(modColorIndex[icomb]);
//      gMass[icomb]->Draw("P");
//    }
//    pLegendFitMSi.Draw();
//    
//    cFitSi->cd(2);
//    gPad->SetGridx();
//    gPad->SetGridy();
//    //  c6->SetLeftMargin(0.12);
//    //  c6->SetRightMargin(0.08);
//    gWidth[0]->Draw("AP");
//    //printf("Width plots\n");
//    pLegendFitWSi.AddEntry(gWidth[0],Form("%s",leg[0].Data()),"P");
//    
//    for(Int_t icomb = 16; icomb <24; icomb++){
//      //printf("option w %d\n",icomb);
//      pLegendFitWSi.AddEntry(gWidth[icomb],Form("%s",leg[icomb].Data()),"P");
//      
//      gWidth[icomb]->SetMarkerStyle(modStyleIndex[icomb]);
//      gWidth[icomb]->SetMarkerColor(modColorIndex[icomb]);
//      gWidth[icomb]->SetLineColor(modColorIndex[icomb]);
//      gWidth[icomb]->Draw("P");
//      
//    }
//    pLegendFitWSi.Draw();
//    
//    cFitSi->Print(Form("IMfigures/%s_%s_Side_MassWidth_%s_%s.%s",prodname.Data(),calorimeter.Data(),part.Data(),filename.Data(),fileFormat.Data()));
//
//  }    
  
  // ---- Pt ----
  
  TLegend pLegendFitPt(0.6,0.6,0.9,0.9);
  pLegendFitPt.SetTextSize(0.03);
  pLegendFitPt.SetFillColor(10);
  pLegendFitPt.SetBorderSize(1);  
  
  TCanvas *cFitPt = new TCanvas("cFitPt","cFitPt",600,600);
  cFitPt->Divide(1, 1);
  
  cFitPt->cd(1);
  
  //   gPad->SetGridx();
  //   gPad->SetGridy();
  gPad->SetLogy();
 
  gPt[0]->SetTitle(Form(" Number of %s vs p_{T}, %s ",particle.Data(),comment[0].Data()));
  gPt[0]->SetMarkerStyle(modStyleIndex[0]);
  gPt[0]->SetMarkerColor(modColorIndex[0]);
  gPt[0]->SetLineColor(modColorIndex[0]);

  gPt[0]->Draw("AP");
  
  //   for(Int_t icomb = 1; icomb <ncomb; icomb++){
  //     gPt[icomb]->SetMarkerStyle(modStyleIndex[icomb]);
  //     gPt[icomb]->SetMarkerColor(modColorIndex[icomb]);
  //     gPt[icomb]->SetLineColor(modColorIndex[icomb]);
  //     gPt[icomb]->Draw("AP");
  //     //pLegendFitPt.AddEntry(gPt[icomb],Form("%s",leg[icomb].Data()),"P");
  //   }
  
  //pLegendFitPt.Draw();
  cFitPt->Print(Form("IMfigures/%s_%s_Pt_%s_%s.%s",prodname.Data(),calorimeter.Data(),part.Data(),filename.Data(),fileFormat.Data()));

  
  //
  
  TGraphErrors *gEff  = 0;
  TGraphErrors *gPrim = 0;
//  if(filename.Contains("simu"))
//  {
//    TH2D* h2Den = (TH2D*) file ->Get("AnaPi0_TM1_PrimPi0Rapidity");
//    TH1D* hY    = h2Den ->ProjectionY("PtYproj",-1,-1);
//    Int_t binP07 = hY->FindBin(+0.65);
//    Int_t binM07 = hY->FindBin(-0.65);
//    TH1D* hDen   = (TH1D*)h2Den->ProjectionX("PrimaryPt" ,binM07   ,binP07   );
//    
//    for(Int_t ibin = 0; ibin < nPt; ibin++ )
//    {
//      Double_t ptMin = xPtLimits[ibin];
//      Double_t ptMax = xPtLimits[ibin+1];
//      primMeson [ibin] = hDen->Integral(hDen->FindBin(ptMin), hDen->FindBin(ptMax));
//      eprimMeson[ibin] = TMath::Sqrt(primMeson[ibin]);
//      //printf("Bin %d, [%1.1f,%1.1f], content num %f, content den %f times evt %e, bin size %f\n",ibin,xPt[ibin],xPt[ibin+1],mesonPtBC[0][ibin],denBin,nEvt,exPt[ibin]*2);
//      
//      primMeson [ibin]/=(nEvt*(exPt[ibin]*2));
//      eprimMeson[ibin]/=(nEvt*(exPt[ibin]*2));
//      if(filename.Contains("pp_7TeV_Pi0"))
//      {
//        primMeson[ibin]*=1.e6;
//        primMeson[ibin]*=360./100.;
//        
//        //eprimMeson[ibin]*=1.e6*1.e6;
//        eprimMeson[ibin]*=360./100.;
//        
//        mesonPtBC [0][ibin] /= (0.62 +xPt[ibin]*0.0017 );
//        emesonPtBC[0][ibin] /= (0.62 +xPt[ibin]*0.0017 );
//        
//        primMeson [ibin] /= (0.56 +xPt[ibin]*0.0096 );
//        eprimMeson[ibin] /= (0.56 +xPt[ibin]*0.0096 );
//      }
//      
//      if(primMeson[ibin] <=0) continue;
//      
//      mesonEffPt [0][ibin] = mesonPtBC [0][ibin] / primMeson[ibin] ;
//      
//      if(mesonEffPt[0][ibin] < 0)
//      {
//        mesonEffPt [0][ibin] = 0;
//        emesonEffPt[0][ibin] = 0;
//      }
//      else       emesonEffPt[0][ibin] = mesonEffPt[0][ibin] * TMath::Sqrt(1./(primMeson[ibin]*nEvt)+1./(mesonPtBC[0][ibin]*nEvt)*2);
//      
//      printf("Bin %d, [%1.1f,%1.1f], Num %2.3e, Den %2.3e, eDen %e, Eff %f, err %f\n",ibin,xPt[ibin],xPt[ibin+1],
//             mesonPtBC[0][ibin],primMeson[ibin],eprimMeson[ibin], mesonEffPt[0][ibin],emesonEffPt[0][ibin]);
//    }
//    
//    gEff =  new TGraphErrors(nPt,xPt,mesonEffPt[0],exPt,emesonEffPt[0]);
//    gEff->SetName("EffxAcc");
//    gEff->GetHistogram()->SetYTitle("#epsilon_{Reco #times PID #times Acc}");
//
//    
//    gPrim =  new TGraphErrors(nPt,xPt,primMeson,exPt,eprimMeson);
//    gPrim->SetName("Primary");
//    gPrim->GetHistogram()->SetYTitle("Primary");
//  }
  
  
  //
  if(mix)
  {
    TCanvas *cFitPtBC = new TCanvas("cFitPtBC","cFitPtBC",600,600);
    cFitPtBC->Divide(1, 1);
    
    cFitPtBC->cd(1);
    
    //   gPad->SetGridx();
    //   gPad->SetGridy();
    gPad->SetLogy();
    //   //  c3->SetLeftMargin(0.12);
    //   //  c3->SetRightMargin(0.08);
    gPtBC[0]->SetTitle(Form(" Number of %s vs p_{T}, bin counting, %s ",particle.Data(),comment[0].Data()));
    gPtBC[0]->SetMarkerStyle(modStyleIndex[0]);
    gPtBC[0]->SetMarkerColor(modColorIndex[0]);
    gPtBC[0]->SetLineColor(modColorIndex[0]);
    //  if(part=="Pi0"){
    //    gPtBC[0]->SetMaximum(1e-1);
    //    gPtBC[0]->SetMinimum(1e-6);
    //  }
    //  else{
    //    gPtBC[0]->SetMaximum(1e-1);
    //    gPtBC[0]->SetMinimum(1e-6);
    //  }

    gPtBC[0]->Draw("AP");
    
      gPt[0]->SetMarkerColor(2);
      gPt[0]->SetMarkerStyle(22);
      
      gPt[0]->Draw("P");
      if(filename.Contains("imu"))
      {
        gPtBCPu->SetMaximum(1e-1);
        gPtBCPu->SetMinimum(1e-6);
        gPtBCPu->Draw("P");
        
    //    gPrim->SetMarkerStyle(25);
    //    gPrim->SetMarkerColor(4);
    //    gPrim->Draw("P");
      }
      
      //   for(Int_t icomb = 1; icomb <ncomb; icomb++){
      //     gPtBC[icomb]->SetMarkerStyle(modStyleIndex[icomb]);
      //     gPtBC[icomb]->SetMarkerColor(modColorIndex[icomb]);
      //     gPtBC[icomb]->SetLineColor(modColorIndex[icomb]);
      //     gPtBC[icomb]->Draw("AP");
      //     //pLegendFitPtBC.AddEntry(gPtBC[icomb],Form("%s",leg[icomb].Data()),"P");
      //   }
      
      //pLegendFitPtBC.Draw();
      
      TLegend legBC(0.5,0.5,0.9,0.9);
      legBC.AddEntry(gPt[0],"Fit","P");
      legBC.AddEntry(gPtBC[0],"Integrated","P");
      legBC.AddEntry(gPtBCPu ,"Integrated pure","P");
      //  legBC.AddEntry(gPrim ,"Primary","P");
      legBC.Draw();
    
    cFitPtBC->Print(Form("IMfigures/%s_%s_PtBC_%s_%s.%s",prodname.Data(),calorimeter.Data(),part.Data(),filename.Data(),fileFormat.Data()));
  }
  
  //--------------------
  // Store plots in file
  //--------------------
  
  TFile * fout = new TFile(Form("IMfigures/%s_%s_MassWidthPtHistograms_%s_%s_%s.root",prodname.Data(),calorimeter.Data(),part.Data(),filename.Data(),hname[0].Data()), "recreate");

  for(Int_t icomb = 0; icomb <ncomb; icomb++)
  {
    for(Int_t ipt = 0; ipt < nPt; ipt++)
    {
      //if(hIM[icomb][ipt]) { hIM[icomb][ipt]->Scale(1./nEvt); hIM[icomb][ipt]->Write(); }
      if(hIM[icomb][ipt]) { hIM[icomb][ipt]->Write(); }
      if(mix)
      {
//        if(hMix         [icomb][ipt]) { hMix         [icomb][ipt]->Scale(1./nEvt); hMix         [icomb][ipt]->Write();}
//        if(hMixCorrected[icomb][ipt]) { hMixCorrected[icomb][ipt]->Scale(1./nEvt); hMixCorrected[icomb][ipt]->Write();}
//        if(hRatio       [icomb][ipt]) { hRatio       [icomb][ipt]->Scale(1./nEvt); hRatio       [icomb][ipt]->Write();}
//        if(hSignal      [icomb][ipt]) { hSignal      [icomb][ipt]->Scale(1./nEvt); hSignal      [icomb][ipt]->Write();}
        if(hMix         [icomb][ipt]) { hMix         [icomb][ipt]->Write();}
        if(hMixCorrected[icomb][ipt]) { hMixCorrected[icomb][ipt]->Write();}
        if(hRatio       [icomb][ipt]) { hRatio       [icomb][ipt]->Write();}
        if(hSignal      [icomb][ipt]) { hSignal      [icomb][ipt]->Write();}
      }
    }
    
    gPt   [icomb]->Write();
    gMass [icomb]->Write();
    gWidth[icomb]->Write();  
    if(mix) gPtBC [icomb]->Write();
  }
  
  if(gEff)    gEff->Write();
  if(gPtBCPu) gPtBCPu->Write();
  //  if(gPrim)   gPrim->Write();
  hRe[0]->Write();
  if(hMi[0])hMi[0]->Write();
  
  fout->Close();
  
}

/// Open the file and the list and the number of analyzed events
/// 
//-----------------------------------------------------------------------------
Bool_t GetFileAndEvents( TString prodname, TString filename, 
                         TString dirName , TString listName )
{
  fil = new TFile(Form("%s/%s.root",prodname.Data(),filename.Data()),"read");
 
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
/// Initialize the fitting function
///
//-----------------------------------------------------------------------------
void SetFitFun()
{
  //Fitting function
  //if(mix) polN = 0;
  
  if(part == "Pi0")
  {
    if ( polN == 0)
      fitfun = new TF1("fitfun",pi0massP0,0.100,0.250,4);
    else if      (polN == 1)
      fitfun = new TF1("fitfun",pi0massP1,0.100,0.250,5);
    else if (polN == 2)
      fitfun = new TF1("fitfun",pi0massP2,0.100,0.250,6);
    else if (polN == 3)
      fitfun = new TF1("fitfun",pi0massP3,0.100,0.250,7);
    else
    {
      printf("*** <<< Set Crystal Ball!!! >>> ***\n");
      fitfun = new TF1("fitfun",CrystalBall,0.100,0.250,6);
    }
    
    if(polN < 4)
    {
      fitfun->SetParLimits(0,  nPairCut/5,nPairCut*1.e4);
      fitfun->SetParLimits(1,  0.105,0.185);
      fitfun->SetParLimits(2,  0.001,0.040);
    }
    else
    {
      fitfun->SetParLimits(0,  nPairCut/5,nPairCut*1.e6);
      fitfun->SetParLimits(2,  0.105,0.185);
      fitfun->SetParLimits(1,  0.001,0.040);
      fitfun->SetParLimits(3,  0.001,0.040);
      fitfun->SetParLimits(4,  0,10);
      fitfun->SetParLimits(5,  1,1.e+6);
    }
    
  } // Pi0
  else  // Eta
  {
    if ( polN == 0)
      fitfun = new TF1("fitfun",pi0massP0,0.400,0.650,4);
    else if      (polN == 1)
      fitfun = new TF1("fitfun",pi0massP1,0.400,0.650,5);
    else if (polN == 2)
      fitfun = new TF1("fitfun",pi0massP2,0.400,0.650,6);
    else if (polN == 3)
      fitfun = new TF1("fitfun",pi0massP3,0.400,0.650,7);
    if(polN < 4){
      fitfun->SetParLimits(0,  nPairCut/10,1.e+6);
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

  if (polN > 1)
  {
    fitfun->SetParName(5,"a_{2}");
    if (polN > 2) fitfun->SetParName(6,"a_{3}"); 
  }  
  
}

//-----------------------------------------------------------------------------
Double_t pi0massP3(Double_t *x, Double_t *par)
{
  Double_t gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                      (2*par[2]*par[2]) );
  Double_t back = par[3] + par[4]*x[0] + par[5]*x[0]*x[0] + par[6]*x[0]*x[0]*x[0];
  return gaus+back;
}

//-----------------------------------------------------------------------------
Double_t pi0massP2(Double_t *x, Double_t *par)
{
  Double_t gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                      (2*par[2]*par[2]) );
  Double_t back = par[3] + par[4]*x[0] + par[5]*x[0]*x[0];
  return gaus+back;
}

//-----------------------------------------------------------------------------
Double_t pi0massP1(Double_t *x, Double_t *par)
{
  Double_t gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                      (2*par[2]*par[2]) );
  Double_t back = par[3] + par[4]*x[0];
  return gaus+back;
}

//-----------------------------------------------------------------------------
Double_t pi0massP0(Double_t *x, Double_t *par)
{
  Double_t gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) /
                                      (2*par[2]*par[2]) );
  Double_t back = 0;//par[3];
  return gaus+back;
}

//-----------------------------------------------------------------------------
Double_t truncatedPolEta(Double_t *x, Double_t *par)
{
  if ((x[0] > 0.425 && x[0] < 0.65) ) {
    TF1::RejectPoint();
    return 0;
  }
  
  return par[0] + par[1]*x[0];// + x[0]*x[0]*par[2]  + x[0]*x[0]*x[0]*par[3];//+ x[0]*x[0]*x[0]*x[0]*par[4];
}

//-----------------------------------------------------------------------------
Double_t truncatedPolPi0(Double_t *x, Double_t *par)
{
  if ((x[0] > 0.07 && x[0] < 0.2)) {
    TF1::RejectPoint();
    return 0;
  }
  
  return par[0] + par[1]*x[0];// + x[0]*x[0]*par[2] + x[0]*x[0]*x[0]*par[3] ;//+ x[0]*x[0]*x[0]*x[0]*par[4];
}

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
