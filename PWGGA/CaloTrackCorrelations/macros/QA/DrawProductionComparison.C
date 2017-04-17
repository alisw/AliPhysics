///
/// \file DrawProductionComparison.C
/// \ingroup CaloTrackCorrMacrosQA
/// \brief Plot N productions comparison analysis of QA histograms from EMCal PWG-GA wagon
///
///
/// Macro to plot comparison of different
/// distributions (spectra, correlations)
/// produced in QA trains but different data
/// Based on the plots provided by DrawAnaCaloTrackQA.C
///
/// To execute: root -q -b -l DrawProductionComparison.C'("Pi0IM_GammaTrackCorr_EMCAL_default","AnalysisResults.root")'
/// The input files must be placed in different directoried,
/// each one defined in the string array "prod"
///     TString prod[] = {"AOD142","AOD115","ESD"};
/// that has to be modified inside the macro.
/// The number of productions has to be specified
///     const Int_t nProd = 3;
/// There is no limitation to the amount of productions
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

// Some global variables
const Int_t nProd = 3; /// total number of productions

TString prod[] = {"AOD142","AOD115","ESD"}; /// productions directory name

TList *list[nProd];

TFile *file[nProd];

TString histoTag = "";

Int_t color[]={kBlack,kRed,kBlue,kOrange+1,kYellow+1,kGreen+2,kCyan+1,kViolet,kMagenta+2,kGray};

Float_t nEvents[nProd] = 0;

//_______________________________________________________________________
/// Main method, produce the plots with the comparisons
///
/// \param listName: Name of list with histograms in file (same for all productions)
/// \param fileName: File name (same for all productions)
//_______________________________________________________________________
void DrawProductionComparison(TString listName = "Pi0IM_GammaTrackCorr_EMCAL_default",
                              TString fileName = "AnalysisResults.root")
{
  printf("Open <%s>; Get List : <%s>\n",fileName.Data(),listName.Data());
  
  histoTag = listName;
  
  //Access the file and list of histograms, global variables
  GetFileAndList(fileName, listName);
  
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(000000);
  gStyle->SetPadRightMargin(0.15);
  //gStyle->SetPadTopMargin(0.02);
  //gStyle->SetPadLeftMargin(0.15);
  gStyle->SetTitleFontSize(0.06);
  
  // Declare the different histograms, arrays input is production
  TH1F* hRaw [nProd];
  TH1F* hCorr[nProd];
  TH1F* hTM  [nProd];
  TH1F* hShSh[nProd];
  
  TH1F* hRatRaw [nProd-1];
  TH1F* hRatCorr[nProd-1];
  TH1F* hRatTM  [nProd-1];
  TH1F* hRatShSh[nProd-1];

  TH1F* hCen   [nProd];
  TH1F* hRatCen[nProd-1];
  
  TH1F* hVertex[3][nProd];
  TH1F* hRatVertex[3][nProd-1];

  TH2F* h2TrackMatchResEtaNeg[nProd];
  TH2F* h2TrackMatchResEtaPos[nProd];
  TH2F* h2TrackMatchResPhiNeg[nProd];
  TH2F* h2TrackMatchResPhiPos[nProd];
  
  TH1F* hTrackMatchResEtaNeg[nProd];
  TH1F* hTrackMatchResEtaPos[nProd];
  TH1F* hTrackMatchResPhiNeg[nProd];
  TH1F* hTrackMatchResPhiPos[nProd];
  
  TH1F* hRatTrackMatchResEtaNeg[nProd-1];
  TH1F* hRatTrackMatchResEtaPos[nProd-1];
  TH1F* hRatTrackMatchResPhiNeg[nProd-1];
  TH1F* hRatTrackMatchResPhiPos[nProd-1];

  TH1F * hTrackPt[nProd] ;
  TH1F * hTrackPtSPD[nProd] ;
  TH1F * hTrackPtNoSPD[nProd] ;
  TH1F * hRatTrackPt[nProd-1] ;
  TH1F * hRatTrackPtSPD[nProd-1] ;
  TH1F * hRatTrackPtNoSPD[nProd-1] ;
  
  TH2F * hTrackEtaPhi[nProd] ;
  TH2F * hTrackEtaPhiSPD[nProd] ;
  TH2F * hTrackEtaPhiNoSPD[nProd] ;
  TH1F * hTrackPhi[nProd] ;
  TH1F * hTrackPhiSPD[nProd] ;
  TH1F * hTrackPhiNoSPD[nProd] ;
  TH1F * hRatTrackPhi[nProd-1] ;
  TH1F * hRatTrackPhiSPD[nProd-1] ;
  TH1F * hRatTrackPhiNoSPD[nProd-1] ;
  
  TH2F* h2XE[nProd];
  TH2F* h2XEUE[nProd];
  TH1F* hXE[nProd];
  TH1F* hXEUE[nProd];
  TH1F* hRatXE[nProd-1];
  TH1F* hRatXEUE[nProd-1];
  
  // Fill the histograms array for each of the productions, do the comparison ratios
  for(Int_t iprod = 0; iprod <  nProd; iprod++)
  {
    // Calorimeter Clusters
    {
      hRaw [iprod] = (TH1F*) GetHisto("AnaPhoton_hPt_Cut_0_Open"    ,iprod);
      hCorr[iprod] = (TH1F*) GetHisto("AnaPhoton_hPt_Cut_4_NCells"  ,iprod);
      hTM  [iprod] = (TH1F*) GetHisto("AnaPhoton_hPt_Cut_7_Matching",iprod);
      hShSh[iprod] = (TH1F*) GetHisto("AnaPhoton_hPt_Cut_9_PID"     ,iprod);

      hRaw [iprod]->Sumw2();
      hCorr[iprod]->Sumw2();
      hTM  [iprod]->Sumw2();
      hShSh[iprod]->Sumw2();
      
      hRaw [iprod]->Scale(1./nEvents[iprod]);
      hCorr[iprod]->Scale(1./nEvents[iprod]);
      hTM  [iprod]->Scale(1./nEvents[iprod]);
      hShSh[iprod]->Scale(1./nEvents[iprod]);
      
      hRaw[iprod]->SetMarkerColor(color[iprod]);
      hRaw[iprod]->SetMarkerStyle(24);
      
      hCorr[iprod]->SetTitle("Cluster spectra with/out cuts");
      hCorr[iprod]->SetYTitle("1/N_{events} dN/dp_{T}");
      hCorr[iprod]->SetTitleOffset(1.5,"Y");
      hCorr[iprod]->SetMarkerColor(color[iprod]);
      hCorr[iprod]->SetMarkerStyle(20);
      hCorr[iprod]->SetAxisRange(0.,30.,"X");
      //hCorr[iprod]->SetMaximum(1.1);
      //hCorr[iprod]->SetMinimum(0);
      
      hTM  [iprod]->SetMarkerColor(color[iprod]);
      hTM  [iprod]->SetMarkerStyle(21);
      
      hShSh[iprod]->SetMarkerColor(color[iprod]);
      hShSh[iprod]->SetMarkerStyle(22);
      
      hRaw [iprod]->SetTitleOffset(1.5,"Y");
      hTM  [iprod]->SetTitleOffset(1.5,"Y");
      hShSh[iprod]->SetTitleOffset(1.5,"Y");
      hCorr[iprod]->SetTitleOffset(1.5,"Y");
      
      if(iprod > 0)
      {
        hRatRaw [iprod-1] = (TH1F*)hRaw [iprod]->Clone(Form("hRatRaw%s" ,prod[iprod].Data()));
        hRatCorr[iprod-1] = (TH1F*)hCorr[iprod]->Clone(Form("hRatCorr%s",prod[iprod].Data()));
        hRatTM  [iprod-1] = (TH1F*)hTM  [iprod]->Clone(Form("hRatTM%s"  ,prod[iprod].Data()));
        hRatShSh[iprod-1] = (TH1F*)hShSh[iprod]->Clone(Form("hRatShSh%s",prod[iprod].Data()));
        
        hRatRaw [iprod-1]->Divide(hRatRaw [iprod-1],hRaw [0],1.000,1,"B");
        hRatCorr[iprod-1]->Divide(hRatCorr[iprod-1],hCorr[0],0.975,1,"B");
        hRatTM  [iprod-1]->Divide(hRatTM  [iprod-1],hTM  [0],0.950,1,"B");
        hRatShSh[iprod-1]->Divide(hRatShSh[iprod-1],hShSh[0],0.925,1,"B");
      }
    }
    
    // Cluster-Track Matching Residuals
    {
      h2TrackMatchResEtaNeg[iprod] = (TH2F*) GetHisto("QA_hTrackMatchedDEtaNeg",iprod);
      h2TrackMatchResEtaPos[iprod] = (TH2F*) GetHisto("QA_hTrackMatchedDEtaPos",iprod);
      h2TrackMatchResPhiNeg[iprod] = (TH2F*) GetHisto("QA_hTrackMatchedDPhiNeg",iprod);
      h2TrackMatchResPhiPos[iprod] = (TH2F*) GetHisto("QA_hTrackMatchedDPhiPos",iprod);

      Float_t binMin = hCorr[iprod]->FindBin(0.5);
      Float_t binMax = hCorr[iprod]->FindBin(2);
      hTrackMatchResEtaNeg[iprod] = (TH1F*) h2TrackMatchResEtaNeg[iprod]->ProjectionY(Form("TMProjEtaNeg%s",prod[iprod].Data()),binMin, binMax);
      hTrackMatchResEtaPos[iprod] = (TH1F*) h2TrackMatchResEtaPos[iprod]->ProjectionY(Form("TMProjEtaPos%s",prod[iprod].Data()),binMin, binMax);
      hTrackMatchResPhiNeg[iprod] = (TH1F*) h2TrackMatchResPhiNeg[iprod]->ProjectionY(Form("TMProjPhiNeg%s",prod[iprod].Data()),binMin, binMax);
      hTrackMatchResPhiPos[iprod] = (TH1F*) h2TrackMatchResPhiPos[iprod]->ProjectionY(Form("TMProjPhiPos%s",prod[iprod].Data()),binMin, binMax);
      
      hTrackMatchResEtaNeg[iprod]->SetXTitle("#Delta #eta");
      hTrackMatchResEtaNeg[iprod]->SetYTitle("entries / N events");
      hTrackMatchResEtaNeg[iprod]->SetTitle("Track-cluster #eta residuals, 0.5 < E < 2 GeV");
      hTrackMatchResEtaNeg[iprod]->SetAxisRange(-0.05,0.05,"X");
      hTrackMatchResEtaNeg[iprod]->Sumw2();
      hTrackMatchResEtaNeg[iprod]->SetMarkerStyle(24);
      hTrackMatchResEtaNeg[iprod]->SetMarkerColor(color[iprod]);
      
      hTrackMatchResEtaPos[iprod]->Sumw2();
      hTrackMatchResEtaPos[iprod]->SetAxisRange(-0.05,0.05,"X");
      hTrackMatchResEtaPos[iprod]->SetMarkerStyle(25);
      hTrackMatchResEtaPos[iprod]->SetMarkerColor(color[iprod]);
      
      hTrackMatchResPhiNeg[iprod]->SetXTitle("#Delta #phi");
      hTrackMatchResPhiNeg[iprod]->SetTitle("Track-cluster #phi residuals, 0.5 < E < 2 GeV");
      hTrackMatchResPhiNeg[iprod]->SetYTitle("entries / N events");
      hTrackMatchResPhiNeg[iprod]->SetAxisRange(-0.05,0.05,"X");
      hTrackMatchResPhiNeg[iprod]->Sumw2();
      hTrackMatchResPhiNeg[iprod]->SetMarkerStyle(24);
      hTrackMatchResPhiNeg[iprod]->SetMarkerColor(color[iprod]);
      
      hTrackMatchResPhiPos[iprod]->Sumw2();
      hTrackMatchResPhiPos[iprod]->SetAxisRange(-0.05,0.05,"X");
      hTrackMatchResPhiPos[iprod]->SetMarkerStyle(25);
      hTrackMatchResPhiPos[iprod]->SetMarkerColor(color[iprod]);
      
      hTrackMatchResEtaNeg[iprod]->Scale(1./nEvents[iprod]);
      hTrackMatchResEtaPos[iprod]->Scale(1./nEvents[iprod]);
      hTrackMatchResPhiNeg[iprod]->Scale(1./nEvents[iprod]);
      hTrackMatchResPhiPos[iprod]->Scale(1./nEvents[iprod]);

      hTrackMatchResEtaNeg[iprod]->SetTitleOffset(1.5,"Y");
      hTrackMatchResEtaPos[iprod]->SetTitleOffset(1.5,"Y");
      hTrackMatchResPhiNeg[iprod]->SetTitleOffset(1.5,"Y");
      hTrackMatchResPhiPos[iprod]->SetTitleOffset(1.5,"Y");
      
      if(iprod > 0)
      {
        hRatTrackMatchResPhiPos[iprod-1] = (TH1F*)hTrackMatchResPhiPos[iprod]->Clone(Form("hRatPhiPos%s",prod[iprod].Data()));
        hRatTrackMatchResPhiNeg[iprod-1] = (TH1F*)hTrackMatchResPhiNeg[iprod]->Clone(Form("hRatPhiNeg%s",prod[iprod].Data()));
        hRatTrackMatchResEtaPos[iprod-1] = (TH1F*)hTrackMatchResEtaPos[iprod]->Clone(Form("hRatEtaPos%s",prod[iprod].Data()));
        hRatTrackMatchResEtaNeg[iprod-1] = (TH1F*)hTrackMatchResEtaNeg[iprod]->Clone(Form("hRatEtaNeg%s",prod[iprod].Data()));
        
        hRatTrackMatchResPhiPos[iprod-1]->Divide(hRatTrackMatchResPhiPos[iprod-1],hTrackMatchResPhiPos[0],1.000,1,"B");
        hRatTrackMatchResPhiNeg[iprod-1]->Divide(hRatTrackMatchResPhiNeg[iprod-1],hTrackMatchResPhiNeg[0],1.000,1,"B");
        hRatTrackMatchResEtaPos[iprod-1]->Divide(hRatTrackMatchResEtaPos[iprod-1],hTrackMatchResEtaPos[0],1.000,1,"B");
        hRatTrackMatchResEtaNeg[iprod-1]->Divide(hRatTrackMatchResEtaNeg[iprod-1],hTrackMatchResEtaNeg[0],1.000,1,"B");
      }
    }
    
    // Hybrid Tracks
    {
      hTrackPt         [iprod] = (TH1F*) GetHisto("AnaHadrons_hPt"                  ,iprod);
      hTrackPtSPD      [iprod] = (TH1F*) GetHisto("AnaHadrons_hPtSPDRefit"          ,iprod);
      hTrackPtNoSPD    [iprod] = (TH1F*) GetHisto("AnaHadrons_hPtNoSPDRefit"        ,iprod);
      hTrackEtaPhiSPD  [iprod] = (TH2F*) GetHisto("AnaHadrons_hEtaPhiSPDRefitPt02"  ,iprod);
      hTrackEtaPhiNoSPD[iprod] = (TH2F*) GetHisto("AnaHadrons_hEtaPhiNoSPDRefitPt02",iprod);
      hTrackEtaPhi     [iprod] = (TH2F*) GetHisto("AnaHadrons_hEtaPhiPositive"      ,iprod);
      hTrackEtaPhi     [iprod]->Add((TH2F*) GetHisto("AnaHadrons_hEtaPhiNegative"   ,iprod));
      
      hTrackPhiSPD     [iprod] = (TH1F*)hTrackEtaPhiSPD  [iprod]->ProjectionY(Form("hTrackPhiSPD%s"  ,prod[iprod].Data()),0,1000);
      hTrackPhiNoSPD   [iprod] = (TH1F*)hTrackEtaPhiNoSPD[iprod]->ProjectionY(Form("hTrackPhiNoSPD%s",prod[iprod].Data()),0,1000);
      hTrackPhi        [iprod] = (TH1F*)hTrackEtaPhi     [iprod]->ProjectionY(Form("hTrackPhi%s"     ,prod[iprod].Data()),0,1000);
      
      hTrackPt     [iprod]->Sumw2();
      hTrackPtSPD  [iprod]->Sumw2();
      hTrackPtNoSPD[iprod]->Sumw2();
      
      hTrackPt     [iprod]->Scale(1./nEvents[iprod]);
      hTrackPtSPD  [iprod]->Scale(1./nEvents[iprod]);
      hTrackPtNoSPD[iprod]->Scale(1./nEvents[iprod]);
      
      hTrackPhi     [iprod]->Sumw2();
      hTrackPhiSPD  [iprod]->Sumw2();
      hTrackPhiNoSPD[iprod]->Sumw2();

      hTrackPhi     [iprod]->Scale(1./nEvents[iprod]);
      hTrackPhiSPD  [iprod]->Scale(1./nEvents[iprod]);
      hTrackPhiNoSPD[iprod]->Scale(1./nEvents[iprod]);
      
      hTrackPt[iprod]->SetTitle("Track spectra with/out SPD");
      hTrackPt[iprod]->SetYTitle("1/N_{events} dN/dp_{T}");
      hTrackPt[iprod]->SetTitleOffset(1.5,"Y");
      hTrackPt[iprod]->SetMarkerColor(color[iprod]);
      hTrackPt[iprod]->SetMarkerStyle(20);
      hTrackPt[iprod]->SetAxisRange(0.,30.,"X");
      //hTrackPt[iprod]->SetMaximum(1.1);
      //hTrackPt[iprod]->SetMinimum(0);
      
      hTrackPtSPD[iprod]->SetMarkerColor(color[iprod]);
      hTrackPtSPD[iprod]->SetMarkerStyle(26);
      
      hTrackPtNoSPD[iprod]->SetMarkerColor(color[iprod]);
      hTrackPtNoSPD[iprod]->SetMarkerStyle(25);
      
      hTrackPhi[iprod]->SetTitle("Track #phi with/out SPD");
      hTrackPhi[iprod]->SetYTitle("1/N_{events} dN/d#phi");
      hTrackPhi[iprod]->SetTitleOffset(1.5,"Y");
      hTrackPhi[iprod]->SetMarkerColor(color[iprod]);
      hTrackPhi[iprod]->SetMarkerStyle(20);
      hTrackPhi[iprod]->SetAxisRange(0.,30.,"X");
      //hTrackPhi[iprod]->SetMaximum(1.1);
      //hTrackPhi[iprod]->SetMinimum(0);
      
      hTrackPhiSPD[iprod]->SetMarkerColor(color[iprod]);
      hTrackPhiSPD[iprod]->SetMarkerStyle(26);
      
      hTrackPhiNoSPD[iprod]->SetMarkerColor(color[iprod]);
      hTrackPhiNoSPD[iprod]->SetMarkerStyle(25);
      
      if(iprod > 0)
      {
        hRatTrackPhi     [iprod-1] = (TH1F*)hTrackPhi     [iprod]->Clone(Form("hRatTrackPhi%s"     ,prod[iprod].Data()));
        hRatTrackPhiNoSPD[iprod-1] = (TH1F*)hTrackPhiNoSPD[iprod]->Clone(Form("hRatTrackPhiNoSPD%s",prod[iprod].Data()));
        hRatTrackPhiSPD  [iprod-1] = (TH1F*)hTrackPhiSPD  [iprod]->Clone(Form("hRatTrackPhiSPD%s"  ,prod[iprod].Data()));
        
        hRatTrackPhi     [iprod-1]->Divide(hRatTrackPhi     [iprod-1],hTrackPhi     [0],1.000,1,"B");
        hRatTrackPhiSPD  [iprod-1]->Divide(hRatTrackPhiSPD  [iprod-1],hTrackPhiSPD  [0],1.000,1,"B");
        hRatTrackPhiNoSPD[iprod-1]->Divide(hRatTrackPhiNoSPD[iprod-1],hTrackPhiNoSPD[0],1.000,1,"B");
        
        hRatTrackPt     [iprod-1] = (TH1F*)hTrackPt     [iprod]->Clone(Form("hRatTrackPt%s"     ,prod[iprod].Data()));
        hRatTrackPtNoSPD[iprod-1] = (TH1F*)hTrackPtNoSPD[iprod]->Clone(Form("hRatTrackPtNoSPD%s",prod[iprod].Data()));
        hRatTrackPtSPD  [iprod-1] = (TH1F*)hTrackPtSPD  [iprod]->Clone(Form("hRatTrackPtSPD%s"  ,prod[iprod].Data()));
        
        hRatTrackPt     [iprod-1]->Divide(hRatTrackPt     [iprod-1],hTrackPt     [0],1.000,1,"B");
        hRatTrackPtSPD  [iprod-1]->Divide(hRatTrackPtSPD  [iprod-1],hTrackPtSPD  [0],1.000,1,"B");
        hRatTrackPtNoSPD[iprod-1]->Divide(hRatTrackPtNoSPD[iprod-1],hTrackPtNoSPD[0],1.000,1,"B");
      }
    }
  
    // photon-track correlation
    {
      h2XE   [iprod]= (TH2F*) GetHisto("AnaPhotonHadronCorr_hXECharged",iprod);
      h2XEUE [iprod]= (TH2F*) GetHisto("AnaPhotonHadronCorr_hXEUeCharged",iprod);
      
      Float_t minClusterE = 8;
      TH1F * hTrigger = (TH1F*) GetHisto("AnaPhotonHadronCorr_hPtTrigger",iprod);
      Int_t minClusterEBin = hTrigger->FindBin(minClusterE);
      Float_t nTrig = hTrigger->Integral(minClusterE,100000);
      
      hXE  [iprod] = (TH1F*)h2XE  [iprod]->ProjectionY(Form("hXE%s"  ,prod[iprod].Data()),minClusterEBin,1000);
      hXEUE[iprod] = (TH1F*)h2XEUE[iprod]->ProjectionY(Form("hXEUE%s",prod[iprod].Data()),minClusterEBin,1000);
      
      hXE  [iprod]->Sumw2();
      hXEUE[iprod]->Sumw2();
      
      hXE  [iprod]->Scale(1./nTrig);
      hXEUE[iprod]->Scale(1./nTrig);
      
      hXE[iprod]->SetTitle(Form("#gamma-hadron x_{E}, p_{T,Trig}>%2.1f GeV/c",minClusterE));
      hXE[iprod]->SetYTitle("1/N_{trigger} dN/dx_{E}");
      hXE[iprod]->SetTitleOffset(1.5,"Y");
      hXE[iprod]->SetMarkerColor(color[iprod]);
      hXE[iprod]->SetMarkerStyle(20);
      hXE[iprod]->SetAxisRange(0.,1.,"X");
      //hXE[iprod]->SetMaximum(1.1);
      //hXE[iprod]->SetMinimum(0);
      
      hXEUE[iprod]->SetMarkerColor(color[iprod]);
      hXEUE[iprod]->SetMarkerStyle(25);
      
      if(iprod > 0)
      {
        hRatXE  [iprod-1] = (TH1F*)hXE  [iprod]->Clone(Form("hRatXE%s"  ,prod[iprod].Data()));
        hRatXEUE[iprod-1] = (TH1F*)hXEUE[iprod]->Clone(Form("hRatXEUE%s",prod[iprod].Data()));
        
        hRatXE  [iprod-1]->Divide(hRatXE  [iprod-1],hXE  [0],1.000,1,"B");
        hRatXEUE[iprod-1]->Divide(hRatXEUE[iprod-1],hXEUE[0],1.000,1,"B");
      }
    }

    //Centrality
    {
      hCen[iprod] = (TH1F*) GetHisto("hCentrality",iprod);
      
      hCen[iprod]->Sumw2();
      
      hCen[iprod]->Scale(1./nEvents[iprod]);
      
      hCen[iprod]->SetTitle("Centrality");
      hCen[iprod]->SetYTitle("1/N_{events} dN/d centrality");
      hCen[iprod]->SetTitleOffset(1.5,"Y");
      hCen[iprod]->SetMarkerColor(color[iprod]);
      hCen[iprod]->SetMarkerStyle(20);
      //hCen[iprod]->SetAxisRange(0.,30.,"X");
      //hCen[iprod]->SetMaximum(1.1);
      //hCen[iprod]->SetMinimum(0);
      
      if(iprod > 0)
      {
        hRatCen[iprod-1] = (TH1F*)hCen[iprod]->Clone(Form("hRatCen%s" ,prod[iprod].Data()));
        
        hRatCen[iprod-1]->Divide(hRatCen[iprod-1],hCen [0],1.000,1,"B");
       }
    }

    //Vertex
    {
      hVertex[0][iprod] = (TH1F*) GetHisto("hZVertex",iprod);
      hVertex[1][iprod] = (TH1F*) GetHisto("hYVertex",iprod);
      hVertex[2][iprod] = (TH1F*) GetHisto("hXVertex",iprod);
      
      hVertex[0][iprod]->Sumw2();
      hVertex[1][iprod]->Sumw2();
      hVertex[2][iprod]->Sumw2();
      
      for(Int_t ivertex = 0; ivertex < 3; ivertex++)
      {
        //hVertex[ivertex][iprod]->Sumw2();
        
        hVertex[ivertex][iprod]->Scale(1./nEvents[iprod]);
        
        //hVertex[ivertex][iprod]->SetTitle("Centrality");
        hVertex[ivertex][iprod]->SetYTitle("1/N_{events} dN/ d vertex");
        hVertex[ivertex][iprod]->SetTitleOffset(1.5,"Y");
        hVertex[ivertex][iprod]->SetMarkerColor(color[iprod]);
        hVertex[ivertex][iprod]->SetLineColor(color[iprod]);
        hVertex[ivertex][iprod]->SetMarkerStyle(20);
        if(ivertex==0)hVertex[ivertex][iprod]->SetAxisRange(-10,10.,"X");
        else          hVertex[ivertex][iprod]->SetAxisRange(-1.5,1.5,"X");
        //hVertex[ivertex][iprod]->SetMaximum(1.1);
        //hVertex[ivertex][iprod]->SetMinimum(0);
        
        if(iprod > 0)
        {
          hRatVertex[ivertex][iprod-1] = (TH1F*)hVertex[ivertex][iprod]->Clone(Form("hRatVertex%s_%d" ,prod[iprod].Data(),ivertex));
          
          hRatVertex[ivertex][iprod-1]->Divide(hRatVertex[ivertex][iprod-1],hVertex[ivertex][0],1.000,1,"B");
        }
      }
    }

  }
  
  /////////////////
  // Make the plots
  /////////////////
  
  //Legend for productions
  TLegend lprod(0.3,0.475,0.84,0.675);
  lprod.SetTextSize(0.04);
  lprod.SetBorderSize(0);
  lprod.SetFillColor(0);
  
  //Calimeter Cluster
  {
    TCanvas * ccalo = new TCanvas(Form("Cluster_%s",histoTag.Data()),"",1000,500);
    ccalo->Divide(2,1);
    
    ccalo->cd(1);
    gPad->SetLogy();
    
    hCorr[0]->Draw();
    for(Int_t iprod = 0; iprod <  nProd; iprod++)
    {
      hRaw [iprod]->Draw("same");
      hCorr[iprod]->Draw("same");
      hTM  [iprod]->Draw("same");
      hShSh[iprod]->Draw("same");
      
      lprod.AddEntry(hRaw[iprod],prod[iprod],"P");
    }
    
    lprod.Draw();
    
    TLegend lcl(0.35,0.7,0.84,0.89);
    lcl.SetTextSize(0.04);
    lcl.SetBorderSize(0);
    lcl.SetFillColor(0);
    lcl.AddEntry(hRaw [0],"Raw","P");
    lcl.AddEntry(hCorr[0],"No Exotics + non lin.","P");
    lcl.AddEntry(hTM  [0],  "+ Track matching","P");
    lcl.AddEntry(hShSh[0],"+ #lambda^{2}_{0} < 0.4","P");
    lcl.Draw();
    
    ccalo->cd(2);
    //gPad->SetLogy();
    
    hRatCorr[0]->SetTitle("Cluster spectra ratio");
    hRatCorr[0]->SetYTitle(Form("data X / %s",prod[0].Data()));
    hRatCorr[0]->SetMinimum(0.850);
    hRatCorr[0]->SetMaximum(1.025);
    hRatCorr[0]->Draw("");
    
    for(Int_t iprod = 0; iprod <  nProd-1; iprod++)
    {
      hRatRaw [iprod]->Draw("same");
      hRatCorr[iprod]->Draw("same");
      hRatTM  [iprod]->Draw("same");
      hRatShSh[iprod]->Draw("same");
    }
    
    TLine l1(0,1,30,1);
    TLine l2(0,0.975,30,0.975);
    TLine l3(0,0.95,30,0.95);
    TLine l4(0,0.925,30,0.925);
    
    l1.Draw("same");
    l2.Draw("same");
    l3.Draw("same");
    l4.Draw("same");
    
    ccalo->Print(Form("%s_ClusterSpectraComp.eps",histoTag.Data()));
  }
  
  //Cluster-Track Matching Residual
  {
    TLine l0(0,hTrackMatchResEtaNeg[0]->GetMinimum(),0,hTrackMatchResEtaNeg[0]->GetMaximum()*1.2);
    
    TLegend lres(0.6,0.75,0.84,0.89);
    lres.SetTextSize(0.04);
    //lres.SetBorderSize(0);
    lres.SetFillColor(0);
    lres.AddEntry(hTrackMatchResEtaNeg[0],"Negative","P");
    lres.AddEntry(hTrackMatchResEtaPos[0],"Positive","P");
    lres.Draw();
    
    TCanvas * ccalo2 = new TCanvas(Form("MatchingResiduals_%s",histoTag.Data()),"",500,500);
    ccalo2->Divide(2,2);
    
    ccalo2->cd(1);
    //gPad->SetLogy();
    
    hTrackMatchResEtaPos[0]->Draw("");
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      hTrackMatchResEtaNeg[iprod]->Draw("same");
      hTrackMatchResEtaPos[iprod]->Draw("same");
    }
    
    l0.Draw("same");
    lres.Draw();
    lprod.Draw();
    ccalo2->cd(2);
    
    hTrackMatchResPhiPos[0]->Draw("");
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      hTrackMatchResPhiNeg[iprod]->Draw("same");
      hTrackMatchResPhiPos[iprod]->Draw("same");
    }
    
    l0.Draw("same");
    
    ccalo2->cd(3);
    //gPad->SetLogy();
    
    hRatTrackMatchResEtaPos[0]->SetMaximum(1.25);
    hRatTrackMatchResEtaPos[0]->SetMinimum(0.98);
    hRatTrackMatchResEtaPos[0]->Draw("");
    hRatTrackMatchResEtaPos[0]->SetYTitle(Form("data X / %s",prod[0].Data()));
    for(Int_t iprod = 0; iprod < nProd-1; iprod++)
    {
      hRatTrackMatchResEtaNeg[iprod]->Draw("same");
      hRatTrackMatchResEtaPos[iprod]->Draw("same");
    }
    
    //l0.Draw("same");
    
    ccalo2->cd(4);
    
    hRatTrackMatchResPhiPos[0]->SetMaximum(1.20);
    hRatTrackMatchResPhiPos[0]->SetMinimum(0.95);
    hRatTrackMatchResPhiPos[0]->Draw("");
    hRatTrackMatchResPhiPos[0]->SetYTitle(Form("data X / %s",prod[0].Data()));
    for(Int_t iprod = 0; iprod < nProd-1; iprod++)
    {
      hRatTrackMatchResPhiNeg[iprod]->Draw("same");
      hRatTrackMatchResPhiPos[iprod]->Draw("same");
    }
    
    ccalo2->Print(Form("%s_MatchingResidualsComp.eps",histoTag.Data()));
  }
  
  // Hybrid tracks
  {
    TLegend ltrack(0.6,0.75,0.84,0.89);
    ltrack.SetTextSize(0.04);
    //ltrack.SetBorderSize(0);
    ltrack.SetFillColor(0);
    ltrack.AddEntry(hTrackPt     [0],"All","P");
    ltrack.AddEntry(hTrackPtSPD  [0],"SPD","P");
    ltrack.AddEntry(hTrackPtNoSPD[0],"No SPD","P");
    
    TCanvas * ctrack = new TCanvas(Form("TrackHisto_%s",histoTag.Data()),"",1500,1500);
    ctrack->Divide(2,2);
    ctrack->cd(1);
    gPad->SetLogy();
    hTrackPt[0]->Draw("");
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      hTrackPt     [iprod]->Draw("same");
      hTrackPtSPD  [iprod]->Draw("same");
      hTrackPtNoSPD[iprod]->Draw("same");
    }
    
    ltrack.Draw();
    lprod.Draw();
    
    ctrack->cd(2);
    
    hRatTrackPt[0]->SetMaximum(1.05);
    hRatTrackPt[0]->SetMinimum(0.95);
    hRatTrackPt[0]->Draw("");
    hRatTrackPt[0]->SetYTitle(Form("data X / %s",prod[0].Data()));
    for(Int_t iprod = 0; iprod < nProd-1; iprod++)
    {
      hRatTrackPt     [iprod]->Draw("same");
      hRatTrackPtSPD  [iprod]->Draw("same");
      hRatTrackPtNoSPD[iprod]->Draw("same");
    }
    
    ctrack->cd(3);
    hTrackPhi[0]->SetMaximum(3.);
    hTrackPhi[0]->SetMinimum(0.);
    hTrackPhi[0]->Draw("");
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      hTrackPhi     [iprod]->Draw("same");
      hTrackPhiSPD  [iprod]->Draw("same");
      hTrackPhiNoSPD[iprod]->Draw("same");
    }
    
    ctrack->cd(4);
    //gPad->SetLogy();
    
    hRatTrackPhi[0]->SetMaximum(1.05);
    hRatTrackPhi[0]->SetMinimum(0.95);
    hRatTrackPhi[0]->Draw("");
    hRatTrackPhi[0]->SetYTitle(Form("data X / %s",prod[0].Data()));
    for(Int_t iprod = 0; iprod < nProd-1; iprod++)
    {
      hRatTrackPhi     [iprod]->Draw("same");
      hRatTrackPhiSPD  [iprod]->Draw("same");
      hRatTrackPhiNoSPD[iprod]->Draw("same");
    }
    
    ctrack->Print(Form("%s_TrackComp.eps",histoTag.Data()));
  }
  
  // XE
  {
    TLegend lxe(0.6,0.75,0.84,0.89);
    lxe.SetTextSize(0.04);
    //lxe.SetBorderSize(0);
    lxe.SetFillColor(0);
    lxe.AddEntry(hXE  [0],"Signal+bkg","P");
    lxe.AddEntry(hXEUE[0],"Und. Event","P");
    
    TCanvas * cxe = new TCanvas(Form("XEHisto_%s",histoTag.Data()),"",1000,500);
    cxe->Divide(2,1);
    cxe->cd(1);
    gPad->SetLogy();
    hXE[0]->Draw("");
    for(Int_t iprod = 0; iprod < nProd; iprod++)
    {
      hXE  [iprod]->Draw("same");
      hXEUE[iprod]->Draw("same");
    }
    
    lxe.Draw();
    lprod.Draw();
    
    cxe->cd(2);
    
    hRatXE[0]->SetMaximum(1.05);
    hRatXE[0]->SetMinimum(0.95);
    hRatXE[0]->Draw("");
    hRatXE[0]->SetYTitle(Form("data X / %s",prod[0].Data()));
    for(Int_t iprod = 0; iprod < nProd-1; iprod++)
    {
      hRatXE  [iprod]->Draw("same");
      hRatXEUE[iprod]->Draw("same");
    }
    
    cxe->Print(Form("%s_XEComp.eps",histoTag.Data()));
  }
  
  // Centrality
  {
    TCanvas * ccen = new TCanvas(Form("Centrality_%s",histoTag.Data()),"",1000,500);
    ccen->Divide(2,1);
    
    ccen->cd(1);
    //gPad->SetLogy();
    
    hCen[0]->Draw();
    for(Int_t iprod = 0; iprod <  nProd; iprod++)
    {
      hCen[iprod]->Draw("same");
    }
    
    lprod.Draw();
    
    ccen->cd(2);
    //gPad->SetLogy();
    
    hRatCen[0]->SetTitle("Centrality");
    hRatCen[0]->SetYTitle(Form("data X / %s",prod[0].Data()));
    hRatCen[0]->SetMinimum(0.95);
    hRatCen[0]->SetMaximum(1.05);
    hRatCen[0]->Draw("");
    
    for(Int_t iprod = 0; iprod <  nProd-1; iprod++)
    {
      hRatCen [iprod]->Draw("same");
    }
    
    TLine l1(0,1,100,1);
    
    l1.Draw("same");

    ccen->Print(Form("%s_CentralityComp.eps",histoTag.Data()));
  }

  // Vertex
  {
    TCanvas * cvertex = new TCanvas(Form("Vertex_%s",histoTag.Data()),"",3*500,2*500);
    cvertex->Divide(3,2);
    Int_t npannel = 1;
    for(Int_t ivertex = 0; ivertex < 3; ivertex++)
    {
      cvertex->cd(npannel);
      //gPad->SetLogy();
      
      hVertex[ivertex][0]->Draw();
      for(Int_t iprod = 0; iprod <  nProd; iprod++)
      {
        hVertex[ivertex][iprod]->Draw("same");
      }
      
      lprod.Draw();
      
      cvertex->cd(npannel+3);
      gPad->SetGridy();
      
      //hRatVertex[ivertex][0]->SetTitle("");
      hRatVertex[ivertex][0]->SetYTitle(Form("data X / %s",prod[0].Data()));
      hRatVertex[ivertex][0]->SetMinimum(0.90);
      hRatVertex[ivertex][0]->SetMaximum(1.10);
      hRatVertex[ivertex][0]->Draw("");
      
      for(Int_t iprod = 0; iprod <  nProd-1; iprod++)
      {
        hRatVertex[ivertex][iprod]->Draw("same");
      }

      npannel++;
    }
    cvertex->Print(Form("%s_VertexComp.eps",histoTag.Data()));
  }
}

//_____________________________________________________
/// Access the file and list with histograms and number
/// of analyzed events per each production.
//_____________________________________________________
void GetFileAndList(TString fileName, TString listName)
{
  for(Int_t iprod = 0; iprod < nProd; iprod++)
  {
    file[iprod]  = new TFile(Form("%s/%s",prod[iprod].Data(),fileName.Data()),"read");
    if(file[iprod]->Get("hNEvents"))
    {
      nEvents[iprod] = ((TH1F*)file[iprod]->Get("hNEvents"))->GetEntries();
      printf("%s: nEvents %e\n",prod[iprod].Data(),nEvents[iprod]);
    }

    TDirectory * dir = (TDirectory*) file[iprod]->Get(listName);
    if(dir)
    {
      list[iprod] = (TList*) dir->Get(listName);
      nEvents[iprod] = ((TH1F*)list[iprod]->FindObject("hNEvents"))->GetEntries();
      printf("%s: nEvents %e\n",prod[iprod],nEvents[iprod]);
    }
    else list[iprod] = 0;
  }
}

//___________________________________
/// Check if the list is available,
/// if not get the histo directly from file
//___________________________________
TObject * GetHisto(TString histoName, Int_t iprod)
{
  if(list[iprod]) return list[iprod]->FindObject(histoName);
  else            return file[iprod]->Get       (histoName);
}
