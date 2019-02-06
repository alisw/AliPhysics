/// \file DrawPtHardBins.C
/// \ingroup CaloTrackCorrMacrosQAPtHard
/// \brief Plot pT hard dependent histograms of analysis EMCal PWG-GA QA wagon
///
/// Macro to plot few selected histograms
/// to QA MC productions done with pT hard bins (pythia jet-jet)./
/// Do plots with the calculated scaling factor per pT hard bin and without the scaling.
/// The plots are:
/// * EMCal or DCal reconstructed cluster energy spectra
/// * Hybrid track pT spectra, and their components (global or no SPD tracks)
/// * Cluster eta/phi hit/acceptance plots
/// * Track phi distribution
/// * Generated pi0 pT spectra, considering also those falling in the EMCal or DCal acceptance
/// * Generated photon pT spectra, considering also those falling in the EMCal or DCal acceptance
/// * Reconstructed invariant mass plots 2D and projections for EMCal/DCal
///
/// To execute:root -q -b -l DrawPtHardBins.C'(1,50,20,1,kFALSE)'
///
/// The input files Scaled.root and NotScaled.root are obtained executing the script:
/// * DownloadExtractScaleMergePtHardAnalysisFiles.sh
///
/// This script: 
///  1 downloads the analysis QA train output files per pT hard bin and run, 
///  1 merges the files per run in each pT hard bin,  
///  1 extracts the histograms of the EMCAL QA wagon in 2 separate files, in one of them it scales the histograms with the proper cross section
///  1 adds the different pT hard bins
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

//---------------------------------------------------------
// Set includes and declare methods for compilation

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMath.h"
#include "TGaxis.h"
#include "TLine.h"

#endif

///
/// Main method, produce the plots.
///
///
/// Input:
/// \param minE: minimum energy/pT on plots with E/pT spectra 
/// \param maxE: maximum energy/pT on plots with E/pT spectra 
/// \param nBin: total number of pT hard bins
/// \param firstBin: first pT hard bin (sometimes it is 0)
/// \param scaleHisto: bool, if true the do the scaling of the unscaled plots
//_______________________________________________________________________
void DrawPtHardBins
(  
 Int_t minE = 1,
 Int_t maxE = 70,
 const Int_t nBin = 20,
 Int_t firstBin = 1,
 Bool_t scaleHisto = kFALSE
 )
{
  TH1F * hNEvents;
  TH1F * hXsec;
  TH1F * hTrials;

  TH1F * hPtHard[nBin][2];
  TH1F * hPtHardSum[2];

  TH1F * hCent[nBin][2];
  TH1F * hCentSum[2];
  
  TH1F * hClusterE[nBin][2];
  TH1F * hClusterESum[2];
  
  TH1F * hClusterD[nBin][2];
  TH1F * hClusterDSum[2];

  TH1F * hTrackPt[nBin][2][3];
  TH1F * hTrackPtSum[2][3];

  TH1F * hPi0[nBin][2];
  TH1F * hPi0Sum[2];
  
  TH1F * hPi0E[nBin][2];
  TH1F * hPi0ESum[2];

  TH1F * hPi0D[nBin][2];
  TH1F * hPi0DSum[2];

  TH1F * hGam[nBin][2];
  TH1F * hGamSum[2];

  TH1F * hGamE[nBin][2];
  TH1F * hGamESum[2];
  
  TH1F * hGamD[nBin][2];
  TH1F * hGamDSum[2];
  
  TH2F * hEtaPhi     [nBin][2];
  TH2F * hCellEtaPhi [nBin][2];
  TH2F * hTrackEtaPhi[nBin][2];
  
  const Int_t nCent = 10;
  TString cenLegend[] = {"0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"};
  TH2F * hIMEMCal[nBin][2][nCent];
  TH2F * hIMDCal [nBin][2][nCent];
  
  TH2F * hEtaPhiSum     [2];
  TH2F * hCellEtaPhiSum [2];
  TH2F * hTrackEtaPhiSum[2];
  
  TH2F * hIMEMCalSum[2][nCent];
  TH2F * hIMDCalSum [2][nCent];
    
  TH2F* hTrackPhiGlobal[nBin][2];
  TH2F* hTrackPhiNoSPD [nBin][2];  
  
  TH2F* hTrackPhiGlobalSum[2];
  TH2F* hTrackPhiNoSPDSum [2];
  
  TFile * f[nBin][2];
  TFile * fTot [2];
  
  Int_t color[] = 
  { kRed   -3,    kRed, kRed   +3, kBlue  -3,   kBlue, kBlue  +3, kGreen -3, kGreen , kGreen +3,
    kViolet-3, kViolet, kViolet+3, kOrange-3, kOrange, kOrange+3, kYellow-3, kYellow, kYellow+3,
    kMagenta-3, kMagenta,kMagenta+3};
  
  Double_t scale[nBin];

  for(Int_t k = 0; k < 2; k++)
  {
    if ( k==1 ) fTot[k] = TFile::Open("Scaled.root"   ,"read");
    else        fTot[k] = TFile::Open("NotScaled.root","read");
    
    for(Int_t i = 0; i < nBin; i++)
    {
      if ( k==1 ) f[i][k] = TFile::Open(Form("%d/ScaledMerged.root"   ,i+firstBin),"read");
      else        f[i][k] = TFile::Open(Form("%d/NotScaledMerged.root",i+firstBin),"read");
    
      if(!f[i][k]) continue;
      
      //printf("i %d, k %d, f %p\n",i,k,f[i][k]);
      
      hPtHard[i][k] = (TH1F*) f[i][k]->Get("hPtHard");
      if(hPtHard[i][k])
      {
        hPtHard[i][k]->SetLineColor(color[i]);
        hPtHard[i][k]->SetLineWidth(2);
        //hPtHard[i][k]->SetAxisRange(minE, maxE,"X");
      }
 
      hCent[i][k] = (TH1F*) f[i][k]->Get("hCentrality");
      if(hCent[i][k])
      {
        hCent[i][k]->SetLineColor(color[i]);
        hCent[i][k]->SetLineWidth(2);
        //hCent[i][k]->SetAxisRange(minE, maxE,"X");
      }
      
      // Recover scaling parameters
      if ( k==0 )
      {
        hNEvents = (TH1F*) f[i][k]->Get("hNEvents");
        hXsec    = (TH1F*) f[i][k]->Get("hXsec");
        hTrials  = (TH1F*) f[i][k]->Get("hTrials");
        
        if(hXsec)
        {
          scale[i] = (hXsec->GetBinContent(1)/hXsec->GetEntries()) / 
          (hTrials->GetBinContent(1)/ hNEvents->GetBinContent(1)) / hNEvents->GetBinContent(1);
          //       hXsec->GetBinContent(1)/hTrials->GetBinContent(1);
          //        1. / hNEvents->GetBinContent(1);
          //        1./  hPtHard[i][k]->GetEntries();
          //        hXsec->GetBinContent(1)/hTrials->GetBinContent(1)/hPtHard[i][k]->GetEntries();
          
          printf("bin i %d, events %2.3e, pT hard entries %2.3e (fraction of pT hard %2.4f),"
                 "chunks %2.0f, xsec %2.3e, trails %2.0f, xsec/chunks %2.3e, trials/nevents %2.3e, scale %2.3e \n",
                 i, hNEvents->GetBinContent(1),hPtHard[i][k]->GetEntries(), hPtHard[i][k]->GetEntries()/hNEvents->GetBinContent(1),
                 hXsec->GetEntries(), hXsec->GetBinContent(1), hTrials->GetBinContent(1),
                 hXsec->GetBinContent(1)/hXsec->GetEntries(),hTrials->GetBinContent(1)/ hNEvents->GetBinContent(1), scale[i] );
        }
        else
        {
          scale[i] = 1./ hNEvents->GetBinContent(1);
        }
      }
      

      hClusterE[i][k] = (TH1F*) f[i][k]->Get("AnaPhoton_Calo0_hEPhoton");
      hClusterE[i][k]->SetLineColor(color[i]);
      hClusterE[i][k]->SetLineWidth(2);
      hClusterE[i][k]->SetAxisRange(minE, maxE,"X");
      
      hClusterD[i][k] = (TH1F*) f[i][k]->Get("AnaPhoton_Calo1_hEPhoton");
      if(hClusterD[i][k])
      {
        hClusterD[i][k]->SetLineColor(color[i]);
        hClusterD[i][k]->SetLineWidth(2);
        hClusterD[i][k]->SetAxisRange(minE, maxE,"X");
      }
      
      for(Int_t j=0; j<3; j++)
      {
        if(!hTrackPt[i][k][j]) continue;
        if(j==0) hTrackPt[i][k][j] = (TH1F*) f[i][k]->Get("AnaHadrons_hPt");
        if(j==1) hTrackPt[i][k][j] = (TH1F*) f[i][k]->Get("AnaHadrons_hPtSPDRefit");
        if(j==2) hTrackPt[i][k][j] = (TH1F*) f[i][k]->Get("AnaHadrons_hPtNoSPDRefit");
        hTrackPt[i][k][j]->SetLineColor(color[i]);
        hTrackPt[i][k][j]->SetLineWidth(2);
        hTrackPt[i][k][j]->SetLineStyle(j);
        hTrackPt[i][k][j]->SetAxisRange(minE, maxE,"X");
      }
      
      hPi0[i][k] = (TH1F*) f[i][k]->Get("AnaPi0_Calo0_hPrimPi0Pt");
      if(hPi0[i][k])
      {
        hPi0[i][k]->SetLineColor(color[i]);
        hPi0[i][k]->SetLineWidth(2);
        hPi0[i][k]->SetAxisRange(minE, maxE,"X");
      }
      
      hPi0E[i][k] = (TH1F*) f[i][k]->Get("AnaPi0_Calo0_hPrimPi0PtInCalo");
      if(hPi0E[i][k])
      {
        hPi0E[i][k]->SetLineColor(color[i]);
        hPi0E[i][k]->SetLineWidth(2);
        hPi0E[i][k]->SetAxisRange(minE, maxE,"X");
      }
      
      hPi0D[i][k] = (TH1F*) f[i][k]->Get("AnaPi0_Calo1_hPrimPi0PtInCalo");
      if(hPi0D[i][k])
      {
        hPi0D[i][k]->SetLineColor(color[i]);
        hPi0D[i][k]->SetLineWidth(2);
        hPi0D[i][k]->SetAxisRange(minE, maxE,"X");
      }
      
      hGam[i][k] = (TH1F*) f[i][k]->Get("AnaPhoton_Calo0_hPtPrim_MCPhoton");
      if(hGam[i][k])
      {
        hGam[i][k]->SetLineColor(color[i]);
        hGam[i][k]->SetLineWidth(2);
        hGam[i][k]->SetAxisRange(minE, maxE,"X");
      }
      
      hGamE[i][k] = (TH1F*) f[i][k]->Get("AnaPhoton_Calo0_hPtPrimAcc_MCPhoton");
      if(hGamE[i][k])
      {
        hGamE[i][k]->SetLineColor(color[i]);
        hGamE[i][k]->SetLineWidth(2);
        hGamE[i][k]->SetAxisRange(minE, maxE,"X");
      }
      
      hGamD[i][k] = (TH1F*) f[i][k]->Get("AnaPhoton_Calo1_hPtPrimAcc_MCPhoton");
      if(hGamD[i][k])
      {
        hGamD[i][k]->SetLineColor(color[i]);
        hGamD[i][k]->SetLineWidth(2);
        hGamD[i][k]->SetAxisRange(minE, maxE,"X");
      }
      
      for(Int_t icent = 0; icent < nCent; icent++)
      {
        hIMEMCal[i][k][icent] = (TH2F*) f[i][k]->Get(Form("AnaPi0_Calo0_hRe_cen%d_pidbit0_asy0_dist1",icent));
        hIMDCal [i][k][icent] = (TH2F*) f[i][k]->Get(Form("AnaPi0_Calo1_hRe_cen%d_pidbit0_asy0_dist1",icent));
      }
      
      hTrackPhiGlobal[i][k] = (TH2F*) f[i][k]->Get("AnaHadrons_hEtaPhiSPDRefitPt02");
      hTrackPhiNoSPD [i][k] = (TH2F*) f[i][k]->Get("AnaHadrons_hEtaPhiNoSPDRefitPt02");
      
      hEtaPhi     [i][k] = (TH2F*) f[i][k]->Get("hEMCALReaderEtaPhi");
      hCellEtaPhi [i][k] = (TH2F*) f[i][k]->Get("QA_Cell_hGridCells");
      hTrackEtaPhi[i][k] = (TH2F*) f[i][k]->Get("AnaHadrons_hEtaPhiNegative");
      hTrackEtaPhi[i][k]->Add((TH2F*) f[i][k]->Get("AnaHadrons_hEtaPhiPositive"));

      if(k==0)
      {
        if(hPtHard  [i][k]) hPtHard  [i][k]->Sumw2();
        if(hCent    [i][k]) hCent    [i][k]->Sumw2();
        hClusterE[i][k]->Sumw2();
        if(hClusterD[i][k]) hClusterD[i][k]->Sumw2();
        
        if(hPi0 [i][k]) hPi0 [i][k]->Sumw2();
        if(hPi0E[i][k]) hPi0E[i][k]->Sumw2();
        if(hPi0D[i][k]) hPi0D[i][k]->Sumw2();     
        if(hGam [i][k]) hGam [i][k]->Sumw2();
        if(hGamE[i][k]) hGamE[i][k]->Sumw2();
        if(hGamD[i][k]) hGamD[i][k]->Sumw2();
        
        for(Int_t j = 0; j < 3; j++) 
        {
          if(hTrackPt[i][k][j])hTrackPt[i][k][j]->Sumw2();
        }
        
        //hEtaPhi     [i][k]->Sumw2();
        //hCellEtaPhi [i][k]->Sumw2();
        //hTrackEtaPhi[i][k]->Sumw2();
        
        for(Int_t icent = 0; icent < nCent; icent++)
        {
          if(hIMEMCal[i][k][icent]) hIMEMCal[i][k][icent]->Sumw2();
          if(hIMDCal [i][k][icent]) hIMDCal [i][k][icent]->Sumw2();
        }
        
        hTrackPhiNoSPD [i][k]->Sumw2();
        hTrackPhiGlobal[i][k]->Sumw2();
      }
      
      // Recover the summed histograms, or scale sum 
      if ( k==1 || (k==0 && !scaleHisto))
      {
        hPtHardSum[k] = (TH1F*) fTot[k]->Get("hPtHard");
        if(hPtHardSum[k])
        {
          hPtHardSum[k]->SetLineColor(1);
          hPtHardSum[k]->SetLineWidth(2);
          //hPtHardSum[k]->SetAxisRange(minE, maxE,"X");
        }
 
        hCentSum[k] = (TH1F*) fTot[k]->Get("hCentrality");
        if(hCentSum[k])
        {
          hCentSum[k]->SetLineColor(1);
          hCentSum[k]->SetLineWidth(2);
          //hCentSum[k]->SetAxisRange(minE, maxE,"X");
        }
        
        hClusterESum[k] = (TH1F*) fTot[k]->Get("AnaPhoton_Calo0_hEPhoton");
        hClusterESum[k]->SetLineColor(1);
        hClusterESum[k]->SetLineWidth(2);
        hClusterESum[k]->SetAxisRange(minE, maxE,"X");
        
        hClusterDSum[k] = (TH1F*) fTot[k]->Get("AnaPhoton_Calo1_hEPhoton");
        if(hClusterDSum[k])
        {
          hClusterDSum[k]->SetLineColor(1);
          hClusterDSum[k]->SetLineWidth(2);
          hClusterDSum[k]->SetAxisRange(minE, maxE,"X");
        }
        
        for(Int_t j = 0; j < 3; j++)
        {
          if(j==0) hTrackPtSum[k][j] = (TH1F*) fTot[k]->Get("AnaHadrons_hPt");
          if(j==1) hTrackPtSum[k][j] = (TH1F*) fTot[k]->Get("AnaHadrons_hPtSPDRefit");
          if(j==2) hTrackPtSum[k][j] = (TH1F*) fTot[k]->Get("AnaHadrons_hPtNoSPDRefit");
          
          if(!hTrackPtSum[k][j]) continue;

          hTrackPtSum[k][j]->SetLineColor(1);
          hTrackPtSum[k][j]->SetLineWidth(2);
          hTrackPtSum[k][j]->SetLineStyle(j);
          hTrackPtSum[k][j]->SetAxisRange(minE, maxE,"X");
        }
        
        hPi0Sum[k] = (TH1F*) fTot[k]->Get("AnaPi0_Calo0_hPrimPi0Pt");
        if(hPi0Sum[k])
        {
          hPi0Sum[k]->SetLineColor(1);
          hPi0Sum[k]->SetLineWidth(2);
          hPi0Sum[k]->SetAxisRange(minE, maxE,"X");
        }
        
        hPi0ESum[k] = (TH1F*) fTot[k]->Get("AnaPi0_Calo0_hPrimPi0PtInCalo");
        if(hPi0ESum[k])
        {
          hPi0ESum[k]->SetLineColor(1);
          hPi0ESum[k]->SetLineWidth(2);
          hPi0ESum[k]->SetAxisRange(minE, maxE,"X");
        }
        
        hPi0DSum[k] = (TH1F*) fTot[k]->Get("AnaPi0_Calo1_hPrimPi0PtInCalo");
        if(hPi0DSum[k])
        {
          hPi0DSum[k]->SetLineColor(1);
          hPi0DSum[k]->SetLineWidth(2);
          hPi0DSum[k]->SetAxisRange(minE, maxE,"X");
        }

        hGamSum[k] = (TH1F*) fTot[k]->Get("AnaPhoton_Calo0_hPtPrim_MCPhoton");
        if(hGamSum[k])
        {
          hGamSum[k]->SetLineColor(1);
          hGamSum[k]->SetLineWidth(2);
          hGamSum[k]->SetAxisRange(minE, maxE,"X");
        }
        
        hGamESum[k] = (TH1F*) fTot[k]->Get("AnaPhoton_Calo0_hPtPrimAcc_MCPhoton");
        if(hGamESum[k])
        {
          hGamESum[k]->SetLineColor(1);
          hGamESum[k]->SetLineWidth(2);
          hGamESum[k]->SetAxisRange(minE, maxE,"X");
        }
        
        hGamDSum[k] = (TH1F*) fTot[k]->Get("AnaPhoton_Calo1_hPtPrimAcc_MCPhoton");
        if(hGamDSum[k])
        {
          hGamDSum[k]->SetLineColor(1);
          hGamDSum[k]->SetLineWidth(2);
          hGamDSum[k]->SetAxisRange(minE, maxE,"X");
        }
        
        for(Int_t icent = 0; icent < nCent; icent++)
        {
          hIMEMCalSum[k][icent] = (TH2F*) fTot[k]->Get(Form("AnaPi0_Calo0_hRe_cen%d_pidbit0_asy0_dist1",icent));
          hIMDCalSum [k][icent] = (TH2F*) fTot[k]->Get(Form("AnaPi0_Calo1_hRe_cen%d_pidbit0_asy0_dist1",icent));
        }
        
        hTrackPhiGlobalSum[k] = (TH2F*) fTot[k]->Get("AnaHadrons_hEtaPhiSPDRefitPt02");
        hTrackPhiNoSPDSum [k] = (TH2F*) fTot[k]->Get("AnaHadrons_hEtaPhiNoSPDRefitPt02");
        
        hEtaPhiSum     [k] = (TH2F*) fTot[k]->Get("hEMCALReaderEtaPhi");
        hCellEtaPhiSum [k] = (TH2F*) fTot[k]->Get("QA_Cell_hGridCells");
        hTrackEtaPhiSum[k] = (TH2F*) fTot[k]->Get("AnaHadrons_hEtaPhiNegative");
        hTrackEtaPhiSum[k]->Add((TH2F*) fTot[k]->Get("AnaHadrons_hEtaPhiPositive"));
        
//        if(k==0)
//        {
//          hPtHardSum  [k]->Sumw2();
//          hClusterESum[k]->Sumw2();
//          hClusterDSum[k]->Sumw2();
//          hPi0Sum     [k]->Sumw2();
//          hPi0ESum    [k]->Sumw2();
//          if(hPi0DSum[k])hPi0DSum[k]->Sumw2();     
//          hGamSum     [k]->Sumw2();
//          hGamESum    [k]->Sumw2();
//          if(hGamDSum[k])hGamDSum[k]->Sumw2();
//          for(Int_t j = 0; j < 3; j++) hTrackPtSum[k][j]->Sumw2();
//          //hEtaPhiSum     [k]->Sumw2();
//          //hCellEtaPhiSum [k]->Sumw2();
//          //hTrackEtaPhiSum[k]->Sumw2();
//          hIMEMCalSum [k]->Sumw2();
//          if(hIMDCalSum[k])hIMDCalSum[k]->Sumw2();
//          hTrackPhiNoSPDSum [k]->Sumw2();
//          hTrackPhiGlobalSum[k]->Sumw2();
//        }
        
      }
      // Scaler and merge
      else if ( scaleHisto && k == 0 )
      {
        if(hPtHard  [i][k])hPtHard[i][k]->Scale(scale[i]);
        if(hCent    [i][k])hCent  [i][k]->Scale(scale[i]);
        hClusterE[i][k]->Scale(scale[i]);
        if(hClusterD[i][k])hClusterD[i][k]->Scale(scale[i]);
        
        if(hPi0 [i][k])hPi0 [i][k]->Scale(scale[i]);
        if(hPi0E[i][k])hPi0E[i][k]->Scale(scale[i]);
        if(hPi0D[i][k])hPi0D[i][k]->Scale(scale[i]);     
        if(hGam [i][k])hGam [i][k]->Scale(scale[i]);
        if(hGamE[i][k])hGamE[i][k]->Scale(scale[i]);
        if(hGamD[i][k])hGamD[i][k]->Scale(scale[i]);
        
        for(Int_t j = 0; j < 3; j++) hTrackPt[i][k][j]->Scale(scale[i]);
        
        hEtaPhi     [i][k]->Scale(scale[i]);
        hCellEtaPhi [i][k]->Scale(scale[i]);
        hTrackEtaPhi[i][k]->Scale(scale[i]);
        
        for(Int_t icent = 0; icent < nCent; icent++)
        {
          if(hIMEMCal[i][k][icent]) hIMEMCal[i][k][icent]->Scale(scale[i]);
          if(hIMDCal [i][k][icent]) hIMDCal [i][k][icent]->Scale(scale[i]);
        }
        
        hTrackPhiNoSPD [i][k]->Scale(scale[i]);
        hTrackPhiGlobal[i][k]->Scale(scale[i]);

        if ( i == 0 ) 
        {
          if(hPtHardSum[k])
          {
            hPtHardSum[k] = (TH1F*) hPtHard[i][k]->Clone("hPtHardSum");
            hPtHardSum[k]->SetLineColor(1);
          }
 
          if(hCentSum[k])
          {
            hCentSum[k] = (TH1F*) hPtHard[i][k]->Clone("hCentSum");
            hCentSum[k]->SetLineColor(1);
          }
          
          hClusterESum[k] = (TH1F*) hClusterE[i][k]->Clone("hClusterESum");
          hClusterESum[k]->SetLineColor(1);

          hClusterDSum[k] = (TH1F*) hClusterD[i][k]->Clone("hClusterDSum");
          if(hClusterDSum[k])hClusterDSum[k]->SetLineColor(1);

          hPi0Sum [k] = (TH1F*) hPi0 [i][k]->Clone("hPi0Sum");
          if(hPi0Sum[k])hPi0Sum [k]->SetLineColor(1);
          
          hPi0ESum[k] = (TH1F*) hPi0E[i][k]->Clone("hPi0ESum");
          if(hPi0ESum[k])hPi0ESum[k]->SetLineColor(1);
          
          hPi0DSum[k] = (TH1F*) hPi0D[i][k]->Clone("hPi0DSum");
          if(hPi0DSum[k])hPi0DSum[k]->SetLineColor(1);

          hGamSum [k] = (TH1F*) hPi0 [i][k]->Clone("hGamSum");
          if(hGamSum[k])hGamSum [k]->SetLineColor(1);
          
          hGamESum[k] = (TH1F*) hPi0E[i][k]->Clone("hGamESum");
          if(hGamESum[k])hGamESum[k]->SetLineColor(1);
          
          hGamDSum[k] = (TH1F*) hPi0D[i][k]->Clone("hGamDSum");
          if(hGamDSum[k])hGamDSum[k]->SetLineColor(1);
          
          for(Int_t j = 0; j < 3; j++)
          {
            if(!hTrackPtSum[k][j]) continue;
            hTrackPtSum[k][j] = (TH1F*) hTrackPt[i][k][j]->Clone(Form("%sSum",hTrackPt[i][k][j]->GetName()));
            hTrackPtSum[k][j]->SetLineColor(1);
          }
          
          hEtaPhiSum     [k] = (TH2F*) hEtaPhi [i][k]->Clone("hEtaPhiSum");
          hCellEtaPhiSum [k] = (TH2F*) hEtaPhi [i][k]->Clone("hCellEtaPhiSum");
          hTrackEtaPhiSum[k] = (TH2F*) hTrackEtaPhi [i][k]->Clone("hTrackEtaPhiSum");
          
          for(Int_t icent = 0; icent < nCent; icent++)
          {
            if(hIMEMCal[i][k][icent]) hIMEMCalSum[k][icent] = (TH2F*) hIMEMCal[i][k][icent]->Clone(Form("hIMEMCalSum_Cent%d",icent));
            else hIMEMCalSum[k][icent] = 0;
            if(hIMDCal [i][k][icent]) hIMDCalSum [k][icent] = (TH2F*) hIMDCal [i][k][icent]->Clone(Form("hIMDCalSum_Cent%d" ,icent));
            else hIMDCalSum [k][icent] = 0;
          }
          
          hTrackPhiGlobalSum[k] = (TH2F*) hTrackPhiGlobal[i][k]->Clone("hTrackPhiGlobalSum");
          hTrackPhiNoSPDSum [k] = (TH2F*) hTrackPhiNoSPD [i][k]->Clone("hTrackPhiNoSPDSum");
        }
        else 
        {
          if(hPtHardSum[k])hPtHardSum[k]->Add(hPtHard[i][k]);
          if(hCentSum  [k])hCentSum  [k]->Add(hCent  [i][k]);

          hClusterESum[k]->Add(hClusterE[i][k]);
          if(hClusterD[i][k])hClusterDSum[k]->Add(hClusterD[i][k]);

          if(hPi0 [i][k])hPi0Sum [k]->Add(hPi0 [i][k]);
          if(hPi0E[i][k])hPi0ESum[k]->Add(hPi0E[i][k]);
          if(hPi0D[i][k])hPi0DSum[k]->Add(hPi0D[i][k]);

          if(hGam [i][k])hGamSum [k]->Add(hGam [i][k]);
          if(hGamD[i][k])hGamESum[k]->Add(hGamE[i][k]);
          if(hGamD[i][k])hGamDSum[k]->Add(hGamD[i][k]);
          
          for(Int_t j = 0; j < 3; j++) 
          {
            if(!hTrackPtSum[k][j]) continue;
            hTrackPtSum[k][j]->Add(hTrackPt[i][k][j]);
          }
          
          hEtaPhiSum     [k]->Add(hEtaPhi     [i][k]);
          hCellEtaPhiSum [k]->Add(hCellEtaPhi [i][k]);
          hTrackEtaPhiSum[k]->Add(hTrackEtaPhi[i][k]);
          
          for(Int_t icent = 0; icent < nCent; icent++)
          {
            if(hIMEMCal[i][k][icent]) hIMEMCalSum[k][icent]->Add(hIMEMCal[i][k][icent]);
            if(hIMDCal [i][k][icent]) hIMDCalSum [k][icent]->Add(hIMDCal [i][k][icent]);
          }
          
          hTrackPhiNoSPDSum [k]->Add(hTrackPhiNoSPD [i][k]);
          hTrackPhiGlobalSum[k]->Add(hTrackPhiGlobal[i][k]);
        }
      }
      
    } // pT hard loop
    
  } // loop scaled/non scaled
  
  
    
  TString scaleCase [] = {"NotScaled" ,"Scaled"};
  TString scaleTitle[] = {"Not scaled","Scaled"};
  TString trackType [] = {"All","Global","noSPD"};
  TString trackTitle[] = {"All hybrid","Global","No SPD"};

  for(Int_t k = 0; k < 2; k++)
  {
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadRightMargin(0.02);
    
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadLeftMargin(0.14);
    
    gStyle->SetOptTitle(1);
    gStyle->SetTitleOffset(1.6,"Y");
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(000000);
    
    TLegend l(0.8,0.3, 0.95, 0.95);
    l.SetBorderSize(0);
    l.SetFillColor(0);
    l.SetTextSize(0.04);
    l.SetHeader(scaleTitle[k]);
    
    //
    // CLUSTER spectrum
    //
    // EMCal
    TCanvas * cE = new TCanvas(Form("cClusterE%d",k),Form("cluster EMCal %s", scaleCase[k].Data()), 200,200);
    
    gPad->SetLogy();
    //gPad->SetLogx();
    hClusterESum[k]->SetTitle(Form("EMCal cluster, %s",scaleTitle[k].Data()));
    hClusterESum[k]->Draw("H");
    //hClusterESum[k]->SetMinimum(1);
    l.AddEntry(hClusterESum[k],"Sum","L");
    
    for(Int_t i = 0; i < nBin; i++)
    {
      if(!f[i][k]) continue;

      if(!hClusterE[i][k]) continue;
      hClusterE[i][k]->Draw("H same");
      l.AddEntry(hClusterE[i][k],Form("Bin %d",i+firstBin),"L");
    }
    
    hClusterESum[k]->Draw("H same");
    
    l.Draw();
    
    cE->Print(Form("Cluster_Energy_EMCal_%s.eps",scaleCase[k].Data()));
    
    // DCal
    if(hClusterDSum[k])
    {
      TCanvas * cD = new TCanvas(Form("cClusterD%d",k),Form("cluster DCal %s", scaleCase[k].Data()), 200,200);
      
      gPad->SetLogy();
      //gPad->SetLogx();
      
      hClusterDSum[k]->SetTitle(Form("DCal cluster, %s",scaleTitle[k].Data()));
      hClusterDSum[k]->Draw("H");
      //hClusterDSum[k]->SetMinimum(1);
      
      for(Int_t i = 0; i < nBin; i++)
      {
        if(!f[i][k]) continue;

        hClusterD[i][k]->Draw("H same");
      }
      
      hClusterDSum[k]->Draw("H same");

      l.Draw();
      
      cD->Print(Form("Cluster_Energy_DCal_%s.eps",scaleCase[k].Data()));
    }
    
    //
    // Parton PT hard spectrum
    //
    if(hPtHardSum[k])
    {
      TCanvas * cHard = new TCanvas(Form("cPtHard%d",k),Form("pT Hard %s", scaleCase[k].Data()), 200,200);
    
      gPad->SetLogy();
      //gPad->SetLogx();
    
      hPtHardSum[k]->SetTitle(Form("Generated parton hard-pT, %s",scaleTitle[k].Data()));
      hPtHardSum[k]->Draw("H");
      //hPtHardSum[k]->SetMinimum(1);
    
      for(Int_t i = 0; i < nBin; i++)
      {
        if(!f[i][k]) continue;

        if(!hPtHard[i][k]) continue;
        hPtHard[i][k]->Draw("H same");
      }
    
      hPtHardSum[k]->Draw("H same");

      l.Draw();
    
      cHard->Print(Form("PtHard_%s.eps",scaleCase[k].Data()));
    }
 
    //
    // Centrality
    //
    if(hCentSum[k])
    {
      TCanvas * cCent = new TCanvas(Form("cCent%d",k),Form("Centrality %s", scaleCase[k].Data()), 200,200);
      
      gPad->SetLogy();
      //gPad->SetLogx();
      
      hCentSum[k]->SetTitle(Form("Centrality, %s",scaleTitle[k].Data()));
      hCentSum[k]->Draw("H");
      //hCentSum[k]->SetMinimum(1);
      
      for(Int_t i = 0; i < nBin; i++)
      {
        if(!f[i][k]) continue;
        
        if(!hCent[i][k]) continue;
        hCent[i][k]->Draw("H same");
      }
      
      hCentSum[k]->Draw("H same");
      
      l.Draw();
      
      cCent->Print(Form("Centrality_%s.eps",scaleCase[k].Data()));
    }
    
    //
    // TRACK spectra
    //
    for(Int_t j=0; j<3; j++)
    {
      if(!hTrackPtSum[k][j]) continue;

      TCanvas * cTr = new TCanvas(Form("cTrackPt%d_Type%d",k,j),
                                  Form("Track Pt Type %s, %s",trackTitle[j].Data(),scaleTitle[k].Data()), 
                                  200,200);
      gPad->SetLogy();
      //gPad->SetLogx();
      
      hTrackPtSum[k][j]->SetTitle(Form("%s tracks, %s",trackTitle[j].Data(),scaleTitle[k].Data()));
      hTrackPtSum[k][j]->Draw("H");
      //hTrackPtSum[k][j]->SetMinimum(1);  
      for(Int_t i = 0; i < nBin; i++)
      {
        if(!f[i][k]) continue;

        if(!hTrackPt[i][k][j]) continue;
        hTrackPt[i][k][j]->Draw("H same");
      }
      
      l.Draw();
      
      cTr->Print(Form("TrackPt_%s_%s.eps",scaleCase[k].Data(),trackType[j].Data()));
    }
    
    //
    // Generated Pi0 spectrum
    //
    // No acceptance selection
    if(hPi0Sum[k])
    {
      TCanvas * cPi0 = new TCanvas(Form("cPi0%d",k),Form("Generated Pi0 %s", scaleCase[k].Data()), 200,200);
    
      gPad->SetLogy();
      //gPad->SetLogx();
    
      hPi0Sum[k]->SetTitle(Form("Generated #pi^{0}, %s",scaleTitle[k].Data()));
      hPi0Sum[k]->Draw("H");
      //hPi0Sum[k]->SetMinimum(1);
    
      for(Int_t i = 0; i < nBin; i++)
      {
        if(!f[i][k]) continue;

        if(!hPi0[i][k]) continue;
        hPi0[i][k]->Draw("H same");
      }
    
      l.Draw();
    
      cPi0->Print(Form("GeneratedPi0_Pt_%s.eps",scaleCase[k].Data()));
    }
    
    // EMCal
    if(hPi0ESum[k])
    {
      TCanvas * cPi0E = new TCanvas(Form("cPi0E%d",k),Form("Generated Pi0 in EMCal acceptance %s", scaleCase[k].Data()), 200,200);
    
      gPad->SetLogy();
      //gPad->SetLogx();
    
      hPi0ESum[k]->SetTitle(Form("Generated #pi^{0} in EMCal, %s",scaleTitle[k].Data()));
      hPi0ESum[k]->Draw("H");
      //hPi0ESum[k]->SetMinimum(1);
    
      for(Int_t i = 0; i < nBin; i++)
      {
        if(!f[i][k]) continue;

        if(!hPi0E[i][k]) continue;
        hPi0E[i][k]->Draw("H same");
      }
    
      l.Draw();
    
      cPi0E->Print(Form("GeneratedPi0_EMCal_Pt_%s.eps",scaleCase[k].Data()));
    }
    
    // DCal
    if(hPi0DSum[k])
    {
      TCanvas * cPi0D = new TCanvas(Form("cPi0D%d",k),Form("Generated Pi0 in DCal acceptance %s", scaleCase[k].Data()), 200,200);
      
      gPad->SetLogy();
      //gPad->SetLogx();
      
      hPi0DSum[k]->SetTitle(Form("Generated #pi^{0} in DCal, %s",scaleTitle[k].Data()));
      hPi0DSum[k]->Draw("H");
      //hPi0DSum[k]->SetMinimum(1);
      
      for(Int_t i = 0; i < nBin; i++)
      {
        if(!f[i][k]) continue;

        if(!hPi0D[i][k]) continue;
        hPi0D[i][k]->Draw("H same");
      }
      
      l.Draw();
      
      cPi0D->Print(Form("GeneratedPi0_DCal_Pt_%s.eps",scaleCase[k].Data()));
    }
    
    //
    // Generated Gamma spectrum
    //
    // No acceptance selection
    if(hGamSum[k])
    {
      TCanvas * cGam = new TCanvas(Form("cGamma%d",k),Form("Generated Gamma %s", scaleCase[k].Data()), 200,200);
    
      gPad->SetLogy();
      //gPad->SetLogx();
    
      hGamSum[k]->SetTitle(Form("Generated #gamma, %s",scaleTitle[k].Data()));
      hGamSum[k]->Draw("H");
      //hGamSum[k]->SetMinimum(1);
    
      for(Int_t i = 0; i < nBin; i++)
      {
        if(!f[i][k]) continue;

        if(!hGam[i][k]) continue;
        hGam[i][k]->Draw("H same");
      }
    
      l.Draw();
    
      cGam->Print(Form("GeneratedGam_Pt_%s.eps",scaleCase[k].Data()));
    }
  
    if(hGamESum[k])
    {
      // EMCal
      TCanvas * cGamE = new TCanvas(Form("cGammaE%d",k),Form("Generated Gamma in EMCal acceptance %s", scaleCase[k].Data()), 200,200);
    
      gPad->SetLogy();
      //gPad->SetLogx();
    
      hGamESum[k]->SetTitle(Form("Generated #gamma in EMCal acceptance, %s",scaleTitle[k].Data()));
      hGamESum[k]->Draw("H");
      //hGamESum[k]->SetMinimum(1);
    
      for(Int_t i = 0; i < nBin; i++)
      {
        if(!f[i][k]) continue;

        if(!hGamE[i][k]) continue;
        hGamE[i][k]->Draw("H same");
      }
    
      l.Draw();
    
      cGamE->Print(Form("GeneratedGamma_EMCal_Pt_%s.eps",scaleCase[k].Data()));
    }
  
    // DCal
    if(hGamDSum[k])
    {
      TCanvas * cGamD = new TCanvas(Form("cGammaD%d",k),Form("Generated Gamma in DCal acceptance %s", scaleCase[k].Data()), 200,200);
      
      gPad->SetLogy();
      //gPad->SetLogx();
      
      hGamDSum[k]->SetTitle(Form("Generated #gamma in DCal acceptance, %s",scaleTitle[k].Data()));
      hGamDSum[k]->Draw("H");
      //hGamDSum[k]->SetMinimum(1);
      
      for(Int_t i = 0; i < nBin; i++)
      {
        if(!f[i][k]) continue;

        if(!hGamD[i][k]) continue;
        hGamD[i][k]->Draw("H same");
      }
      
      l.Draw();
      
      cGamD->Print(Form("GeneratedGamma_DCal_Pt_%s.eps",scaleCase[k].Data()));
    }
    
    ////////////////////////////
    // Inv. Mass. Projections //
    ////////////////////////////
    TLine pi0Mass(0.135,0,0.135,1e6);
    pi0Mass.SetLineColor(2);
    pi0Mass.SetLineStyle(2);
    
    Bool_t allCent = kFALSE;
    for(Int_t icent = 1; icent < nCent; icent++) { if ( hIMEMCalSum[k][icent] ) allCent = kTRUE; }
    if ( !allCent ) cenLegend[0] = "";
    
    for(Int_t icent = 0; icent < nCent; icent++)
    {
      if ( !hIMEMCalSum[k][icent] ) continue;
      
      if ( icent > 0 && !allCent ) continue;
      
      Int_t binmin = hIMEMCalSum[k][icent]->GetXaxis()->FindBin(5);
      Int_t binmax = hIMEMCalSum[k][icent]->GetXaxis()->FindBin(10);
      //printf("Bins Min %d, Max %d\n",binmin,binmax);
    
      TH1F * hInvMassEMC = (TH1F*) hIMEMCalSum[k][icent]->ProjectionY(Form("hEMCInvMass%d_Cent%d",k,icent),binmin,binmax);
      hInvMassEMC->SetLineColor(1);
    
      TH1F * hInvMassDMC = 0;
      if(hIMDCalSum[k])
      {
        hInvMassDMC = (TH1F*) hIMDCalSum[k][icent] ->ProjectionY(Form("hDMCInvMass%d_Cent%d",k,icent),binmin,binmax);
        hInvMassDMC->SetLineColor(4);
      }
    
      TCanvas * cIMProj = new TCanvas(Form("cIMProj%d_Cent%d",k,icent),
                                      Form("DCal/EMCa; Inv. Mass; %s %s",scaleCase[k].Data(),cenLegend[icent].Data()), 
                                      200,200);
    
      //gPad->SetLogz();
      gPad->SetGridx();
    
      hInvMassEMC->SetAxisRange(0.05, 0.3);
      hInvMassEMC->SetMinimum(1);
    
      hInvMassEMC->SetTitle(Form("Cluster pair M in #pi^{0} region, %s %s",scaleTitle[k].Data(),cenLegend[icent].Data()));

      Double_t intE = hInvMassEMC->Integral();
      if(intE > 0) hInvMassEMC->Scale(1./intE);
      Double_t maxE = hInvMassEMC->GetMaximum();
      hInvMassEMC->SetYTitle("1/N dN / dM ");
      hInvMassEMC->Draw("H");
    
      TLegend lim(0.68,0.6,0.98,0.8);
      lim.SetHeader("5 < E < 10 GeV");
      lim.AddEntry(hInvMassEMC,"EMCal","L");
    
      if(hInvMassDMC)
      {
        Double_t intD = hInvMassDMC->Integral();
        if(intD > 0) hInvMassDMC->Scale(1./intD);
      
        Double_t maxD = hInvMassDMC->GetMaximum();
        if(maxD > maxE) hInvMassEMC->SetMaximum(hInvMassEMC->GetMaximum()*1.1);
      
        hInvMassDMC->Draw("H same");
        lim.AddEntry(hInvMassDMC,"DCal","L");
      }
      pi0Mass.Draw();
      lim.Draw();
    
      cIMProj->Print(Form("InvMassDCalEMCal_5_10GeV_%s_Cen%d.eps",scaleCase[k].Data(),icent));
    
      // Projection 2
      binmin = hIMEMCalSum[k][icent]->GetXaxis()->FindBin(2);
      binmax = hIMEMCalSum[k][icent]->GetXaxis()->FindBin(4);
      //printf("Bins Min %d, Max %d\n",binmin,binmax);
    
      hInvMassEMC = (TH1F*) hIMEMCalSum[k][icent]->ProjectionY(Form("hEMCInvMass2_%d_Cen%d",k,icent),binmin,binmax);
      hInvMassEMC->SetLineColor(1);
    
      hInvMassDMC = 0;
      if(hIMDCalSum[k])
      {
        hInvMassDMC = (TH1F*) hIMDCalSum[k][icent] ->ProjectionY(Form("hDMCInvMass2_%d_Cen%d",k,icent),binmin,binmax);
        hInvMassDMC->SetLineColor(4);
      }
    
      TCanvas * cIMProj2 = new TCanvas(Form("cIMProj2_%d",k),
                                       Form("DCal/EMCa; Inv. Mass; %s %s",scaleCase[k].Data(),cenLegend[icent].Data()), 
                                       200,200);
    
      //gPad->SetLogy();
      gPad->SetGridx();
    
      hInvMassEMC->SetAxisRange(0.05, 0.3);
      hInvMassEMC->SetMinimum(1);
    
      hInvMassEMC->SetTitle(Form("Cluster pair M in #pi^{0} region, %s %s",scaleTitle[k].Data(),cenLegend[icent].Data()));
    
      intE = hInvMassEMC->Integral();
      if(intE > 0) hInvMassEMC->Scale(1./intE);
      maxE = hInvMassEMC->GetMaximum();
      hInvMassEMC->SetYTitle("1/N dN / dM ");
      hInvMassEMC->Draw("H");
    
      TLegend lim2(0.68,0.6,0.98,0.8);
      lim2.SetHeader("2 < E < 4 GeV");
      lim2.AddEntry(hInvMassEMC,"EMCal","L");
    
      if(hInvMassDMC)
      {
        Double_t intD = hInvMassDMC->Integral();
        if(intD > 0) hInvMassDMC->Scale(1./intD);
      
        Double_t maxD = hInvMassDMC->GetMaximum();
        if(maxD > maxE) hInvMassEMC->SetMaximum(hInvMassEMC->GetMaximum()*1.1);

        Double_t minD = hInvMassDMC->GetMinimum();
        if(minD < hInvMassEMC->GetMinimum()) hInvMassEMC->SetMinimum(hInvMassDMC->GetMinimum());
       
        hInvMassDMC->Draw("H same");
        lim2.AddEntry(hInvMassDMC,"DCal","L");
      }
    
      pi0Mass.Draw();
      lim2.Draw();
    
      cIMProj2->Print(Form("InvMassDCalEMCal_2_4GeV_%s_Cen%d.eps",scaleCase[k].Data(),icent));
    }
    
    ////////////////////////////
    // Track Phi Projections //
    ////////////////////////////
    
    TCanvas * cPhiProj = new TCanvas(Form("cPhiProj%d",k),Form("Track phi; %s",scaleCase[k].Data()), 200,200);

    TH1F* hPhiSPD   = (TH1F*)hTrackPhiGlobalSum[k]->ProjectionY(Form("%s_hTrackPhiSPD"  ,scaleCase[k].Data()),0,1000);
    TH1F* hPhiNoSPD = (TH1F*)hTrackPhiNoSPDSum [k]->ProjectionY(Form("%s_hTrackPhiNoSPD",scaleCase[k].Data()),0,1000);
    TH1F* hPhi      = (TH1F*)hPhiSPD->Clone(Form("%s_hTrackPhi",scaleCase[k].Data()));
    hPhi->Add(hPhiNoSPD);
    
    Float_t normFactor = 1./hPhi->Integral();
    hPhi     ->Scale(normFactor);
    hPhiSPD  ->Scale(normFactor);
    hPhiNoSPD->Scale(normFactor);
    
    hPhi     ->SetTitle(Form("Hybrid track type #varphi, 0.2<#it{p}_{T}<2 GeV/#it{c}, %s",scaleTitle[k].Data()));
    hPhi     ->SetLineColor(1);
    hPhiSPD  ->SetLineColor(2);
    hPhiNoSPD->SetLineColor(4);
    
    hPhi     ->SetMinimum(0);
    hPhi     ->SetMaximum(hPhi->GetMaximum()*1.3);
    hPhi     ->SetTitleOffset(1.5,"Y");
    hPhi     ->SetYTitle("Entries");
    
    TGaxis::SetMaxDigits(3);
    
    hPhi     ->Draw("H");
    hPhiSPD  ->Draw("Hsame");
    hPhiNoSPD->Draw("Hsame");
    
    TLegend lphi(0.2,0.75,0.4,0.89);
    lphi.SetTextSize(0.04);
    lphi.AddEntry(hPhi,"Sum","L");
    lphi.AddEntry(hPhiSPD  ,"SPD+Refit","L");
    lphi.AddEntry(hPhiNoSPD,"No SPD+Refit","L");
    lphi.SetBorderSize(0);
    lphi.SetFillColor(0);
    lphi.Draw();
    
    cPhiProj->Print(Form("TrackPhi_%s.eps",scaleCase[k].Data()));

    //////////////
    // Acceptance
    //////////////
    
    gStyle->SetPadRightMargin(0.12);

    
    TCanvas * cEtaPhi = new TCanvas(Form("cEtaPhi%d",k),Form("cluster Eta/phi, %s",scaleCase[k].Data()), 200,200);
    
    gPad->SetLogz();
    //gPad->SetLogx();
    
    hEtaPhiSum[k]->SetAxisRange(-1,1,"X");
    hEtaPhiSum[k]->SetAxisRange( 1,6,"Y");
    
    hEtaPhiSum[k]->SetTitle(Form("EMCal/DCal Cluster acceptance, %s",scaleTitle[k].Data()));

    hEtaPhiSum[k]->Draw("colz");
    
    cEtaPhi->Print(Form("EtaPhi_Cluster_%s.eps",scaleCase[k].Data()));
    
    TCanvas * cCellEtaPhi = new TCanvas(Form("cCellEtaPhi%d",k),Form("cell Eta/phi, %s",scaleCase[k].Data()), 200,200);
    
    gPad->SetLogz();
    //gPad->SetLogx();
    
    hCellEtaPhiSum[k]->SetTitle(Form("EMCal/DCal cell acceptance, %s",scaleTitle[k].Data()));
    
    hCellEtaPhiSum[k]->Draw("colz");
    
    cCellEtaPhi->Print(Form("EtaPhi_Cell_%s.eps",scaleCase[k].Data()));
    
    TCanvas * cTrackEtaPhi = new TCanvas(Form("cTrackEtaPhi%d",k),Form("track Eta/phi, %s",scaleCase[k].Data()), 200,200);
    
    gPad->SetLogz();
    //gPad->SetLogx();
    
    hTrackEtaPhiSum[k]->SetTitle(Form("Hybrid track acceptance, %s",scaleTitle[k].Data()));
    
    hTrackEtaPhiSum[k]->Draw("colz");
    
    cTrackEtaPhi->Print(Form("EtaPhi_Track_%s.eps",scaleCase[k].Data()));

    for(Int_t icent = 0; icent < nCent; icent++)
    {
      if ( !hIMEMCalSum[k][icent] ) continue;
      
      if ( icent > 0 && !allCent ) continue;
      //////////////////////
      // EMCal Invariant mass
      //////////////////////
      TCanvas * cIM = new TCanvas(Form("cIM%d_Cen%d",k,icent),
                                  Form("EMCal Inv. Mass, %s %s",scaleCase[k].Data(),cenLegend[icent].Data()), 
                                  200,200);
    
      gPad->SetLogz();
      //gPad->SetLogx();
    
      hIMEMCalSum[k][icent]->SetAxisRange(4,20,"X");
    
      hIMEMCalSum[k][icent]->SetTitle(Form("EMCal cluster invariant mass, %s %s",
                                           scaleTitle[k].Data(),cenLegend[icent].Data()));
    
      hIMEMCalSum[k][icent]->Draw("colz");
    
      cIM->Print(Form("InvMassEMCal_%s_Cen%d.eps",scaleCase[k].Data(),icent));
    
      //////////////////////
      // DCal Invariant mass
      //////////////////////
      if(hIMDCalSum[k])
      {
        TCanvas * cIMD = new TCanvas(Form("cIMD%d_Cen%d",k,icent),
                                     Form("DCal Inv. Mass, %s %s",scaleCase[k].Data(),cenLegend[icent].Data()),
                                     200,200);
      
        gPad->SetLogz();
        //gPad->SetLogx();
      
        hIMDCalSum[k][icent]->SetAxisRange(4,20,"X");
      
        hIMDCalSum[k][icent]->SetTitle(Form("DCal cluster invariant mass, %s %s",
                                            scaleTitle[k].Data(),cenLegend[icent].Data()));
      
        hIMDCalSum[k][icent]->Draw("colz");
      
        cIMD->Print(Form("InvMassDCal_%s_Cen%d.eps",scaleCase[k].Data(),icent));
      }
   } // centrality
    
  } 
}


