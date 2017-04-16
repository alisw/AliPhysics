/// \file DrawAnaCaloTrackQA.C
/// \ingroup CaloTrackCorrMacrosQA
/// \brief Plot analysis QA histograms from EMCal PWG-GA wagon
///
/// Macro to plot few selected histograms
/// to QA data productions at 0th order
/// Analysis performed with the wagon
/// AddTaskPi0IMGammaCorrQA.C
/// It generates 5 eps plots, each containing 2 to 4 canvases
///
/// To execute: root -q -b -l DrawAnaCaloTrackQA.C'("Pi0IM_GammaTrackCorr_EMCAL_default","AnalysisResults.root")'
/// The input list name might change depending on the wagon / data type
/// In case output file is too large, possiblity to dump the list content in a sepate file:  export = kTRUE
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

// Some global variables
TList *list = 0;
TFile *file = 0;
TString histoTag = "";
Int_t color[]={kBlack,kRed,kOrange+1,kYellow+1,kGreen+2,kBlue,kCyan+1,kViolet,kMagenta+2,kGray};

///
/// Main method, produce the plots for the 5 different types of analysis:
/// * Calorimeter QA in CaloQA method
/// * Track QA in TrackQA method
/// * Invariant mass plots in Pi0QA method
/// * Cluster-track correlation plots in CorrelationQA method
/// * Dedicated generated particles QA in MCQA method
/// Input:
/// \param listName: Name of list with histograms in file
/// \param fileName: File name
/// \param export: export list with histograms to separate file, intereting in case of big output file.
//_______________________________________________________________________
void DrawAnaCaloTrackQA(TString listName = "Pi0IM_GammaTrackCorr_EMCAL_default",
                        TString fileName = "AnalysisResults.root",
                        Bool_t export = kFALSE)
{
  printf("Open <%s>; Get List : <%s>; Export list? <%d>\n",fileName.Data(),listName.Data(),export);
  
  histoTag = listName;
  
  //Access the file and list of histograms, global variables
  GetFileAndList(fileName, listName, export);
  
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(000000);
  gStyle->SetPadRightMargin(0.15);
  //gStyle->SetPadTopMargin(0.02);
  //gStyle->SetPadLeftMargin(0.15);
  gStyle->SetTitleFontSize(0.06);

  //Plot basic Calorimeter QA
  CaloQA();

  //Plot basic Track QA
  TrackQA();

  //Plot basic Pi0 QA
  Pi0QA();

  //Plot basic correlation QA
  CorrelationQA();
  
  // MC basic QA plots, cluster origins (only if it run on MC)
  MCQA();
}

///
/// Plot basic calorimeter QA histograms.
///
//______________________________________
void CaloQA()
{ 
  TCanvas * ccalo = new TCanvas(Form("CaloHisto_%s",histoTag.Data()),"",1000,1000);
  ccalo->Divide(2,2);
  
  ccalo->cd(1);
  gPad->SetLogy();
  
  TH1F* hCellAmplitude = (TH1F*) GetHisto("QA_hAmplitude");
  TH1F* hClusterEnergy = (TH1F*) GetHisto("QA_hE");
  
  hClusterEnergy->SetYTitle("entries");
  hClusterEnergy->SetTitle("Cluster-cell energy spectra");
  hClusterEnergy->Sumw2();
  hClusterEnergy->SetMarkerColor(1);
  hClusterEnergy->SetMarkerStyle(20);
  hClusterEnergy->SetAxisRange(0.,50.,"X");
  hClusterEnergy->Draw();
  
  hCellAmplitude->Sumw2();
  hCellAmplitude->SetMarkerColor(4);
  hCellAmplitude->SetMarkerStyle(25);
  hCellAmplitude->Draw("same");
  
  TLegend l(0.25,0.7,0.83,0.85);
  l.SetTextSize(0.04);
  l.AddEntry(hClusterEnergy,"Cluster (no exotic+non lin.)","P");
  l.AddEntry(hCellAmplitude,"Cell","P");
  l.SetBorderSize(0);
  l.SetFillColor(0);
  l.Draw();
  
  
  ccalo->cd(2);
  //gPad->SetLogy();
  
  TH1F* hRaw  = (TH1F*) GetHisto("AnaPhoton_hPt_Cut_0_Open");
  TH1F* hCorr = (TH1F*) GetHisto("AnaPhoton_hPt_Cut_4_NCells");
  TH1F* hTM   = (TH1F*) GetHisto("AnaPhoton_hPt_Cut_7_Matching");
  TH1F* hShSh = (TH1F*) GetHisto("AnaPhoton_hPt_Cut_9_PID");
  
  hRaw->Sumw2();
  
  hCorr->SetTitle("Ratio after cluster cuts application");
  hCorr->SetYTitle("Selected clusters / Raw clusters");
  hCorr->SetTitleOffset(1.5,"Y");
  hCorr->Sumw2();
  hCorr->SetMarkerColor(1);
  hCorr->SetMarkerStyle(20);
  hCorr->Divide(hRaw);
  hCorr->SetAxisRange(0.,30.,"X");
  hCorr->SetMaximum(1.1);
  hCorr->SetMinimum(0);
  hCorr->Draw();
  
  hTM  ->Sumw2();
  hTM  ->SetMarkerColor(2);
  hTM  ->SetMarkerStyle(21);
  hTM  ->Divide(hRaw);
  hTM  ->Draw("same");
  
  hShSh->Sumw2();
  hShSh->SetMarkerColor(4);
  hShSh->SetMarkerStyle(22);
  hShSh->Divide(hRaw);
  hShSh->Draw("same");
  
  TLegend l2(0.45,0.8,0.95,0.93);
  l2.SetTextSize(0.04);
  l2.AddEntry(hCorr,"No Exotics + non lin.","P");
  l2.AddEntry(hTM,  "+ Track matching","P");
  l2.AddEntry(hShSh,"+ #lambda^{2}_{0} < 0.4","P");
  l2.SetBorderSize(0);
  l2.SetFillColor(0);
  l2.Draw();
  
  
  // Plot track-matching residuals
  // first test did not have this histogram, add protection
  TH2F* hTrackMatchResEtaPhi = (TH2F*) GetHisto("QA_hTrackMatchedDEtaDPhiPos");
  if(hTrackMatchResEtaPhi)
  {
    hTrackMatchResEtaPhi ->Add( (TH2F*) GetHisto("QA_hTrackMatchedDEtaDPhiNeg") );
    
    ccalo->cd(3);
    gPad->SetLogz();
    
    hTrackMatchResEtaPhi->SetAxisRange(-0.025,0.025,"X");
    hTrackMatchResEtaPhi->SetAxisRange(-0.025,0.025,"Y");
    hTrackMatchResEtaPhi->SetTitleOffset(1.5,"Y");
    hTrackMatchResEtaPhi->SetTitle("Track-cluster residual #Delta #phi vs #Delta #eta, E > 0.5 GeV");
    hTrackMatchResEtaPhi->SetXTitle("#Delta #eta");
    hTrackMatchResEtaPhi->SetYTitle("#Delta #phi");
    hTrackMatchResEtaPhi->SetZTitle("entries");
    hTrackMatchResEtaPhi->Draw("colz");
    
    ccalo->cd(4);
    gPad->SetLogy();
    
    TH2F* h2TrackMatchResEtaNeg = (TH2F*) GetHisto("QA_hTrackMatchedDEtaNeg");
    TH2F* h2TrackMatchResEtaPos = (TH2F*) GetHisto("QA_hTrackMatchedDEtaPos");
    TH2F* h2TrackMatchResPhiNeg = (TH2F*) GetHisto("QA_hTrackMatchedDPhiNeg");
    TH2F* h2TrackMatchResPhiPos = (TH2F*) GetHisto("QA_hTrackMatchedDPhiPos");
    
    Float_t binMin = hCorr->FindBin(0.5);
    TH1F* hTrackMatchResEtaNeg = (TH1F*) h2TrackMatchResEtaNeg->ProjectionY("TMProjEtaNeg",binMin, 1000);
    TH1F* hTrackMatchResEtaPos = (TH1F*) h2TrackMatchResEtaPos->ProjectionY("TMProjEtaPos",binMin, 1000);
    TH1F* hTrackMatchResPhiNeg = (TH1F*) h2TrackMatchResPhiNeg->ProjectionY("TMProjPhiNeg",binMin, 1000);
    TH1F* hTrackMatchResPhiPos = (TH1F*) h2TrackMatchResPhiPos->ProjectionY("TMProjPhiPos",binMin, 1000);
    
    hTrackMatchResEtaNeg->SetXTitle("#Delta #eta, #Delta #phi");
    hTrackMatchResEtaNeg->SetYTitle("entries");
    hTrackMatchResEtaNeg->SetTitle("Track-cluster residuals, E > 1 GeV");
    hTrackMatchResEtaNeg->SetAxisRange(-0.05,0.05,"X");
    hTrackMatchResEtaNeg->Sumw2();
    hTrackMatchResEtaNeg->SetMarkerStyle(25);
    hTrackMatchResEtaNeg->SetMarkerColor(2);
    hTrackMatchResEtaNeg->Draw("");
    
    hTrackMatchResEtaPos->Sumw2();
    hTrackMatchResEtaPos->SetMarkerStyle(25);
    hTrackMatchResEtaPos->SetMarkerColor(4);
    hTrackMatchResEtaPos->Draw("same");
    
    hTrackMatchResPhiNeg->Sumw2();
    hTrackMatchResPhiNeg->SetMarkerStyle(24);
    hTrackMatchResPhiNeg->SetMarkerColor(2);
    hTrackMatchResPhiNeg->Draw("same");
    
    hTrackMatchResPhiPos->Sumw2();
    hTrackMatchResPhiPos->SetMarkerStyle(24);
    hTrackMatchResPhiPos->SetMarkerColor(4);
    hTrackMatchResPhiPos->Draw("same");
    
    TLine l0(0,hTrackMatchResEtaNeg->GetMinimum(),0,hTrackMatchResEtaNeg->GetMaximum()*1.2);
    l0.Draw("same");
    
    TLegend l3(0.55,0.7,0.83,0.85);
    l3.SetTextSize(0.04);
    l3.AddEntry(hTrackMatchResEtaNeg,"#Delta #eta, Negative","P");
    l3.AddEntry(hTrackMatchResEtaPos,"#Delta #eta, Positive","P");
    l3.AddEntry(hTrackMatchResPhiNeg,"#Delta #phi, Negative","P");
    l3.AddEntry(hTrackMatchResPhiPos,"#Delta #phi, Positive","P");
    l3.SetBorderSize(0);
    l3.SetFillColor(0);
    l3.Draw();
  }
  ccalo->Print(Form("%s_CaloHisto.eps",histoTag.Data()));
  
  
  TCanvas * ccalo2 = new TCanvas(Form("CaloHisto2_%s",histoTag.Data()),"",500,500);
  ccalo2->Divide(2,2);
  
  ccalo2->cd(3);
//  gPad->SetLogz();
//  TH2F* hCellAmpId   = (TH2F*) GetHisto("QA_hAmpId");
//  hCellAmpId->SetTitle("Cell Id vs energy");
//  hCellAmpId->SetYTitle("Cell Id");
//  //hCellAmpId->SetAxisRange(300.,900.,"Y");
//  hCellAmpId->SetAxisRange(0.,30.,"X");
//  hCellAmpId->SetTitleOffset(1.5,"Y");
//  hCellAmpId->Draw("colz");

  gPad->SetLogz();
  
  TH2F* hClusterTime   = (TH2F*) GetHisto("QA_hClusterTimeEnergy");
  hClusterTime->SetTitle("Cluster energy vs time");
  hClusterTime->SetYTitle("time (ns)");
  hClusterTime->SetAxisRange(300.,900.,"Y");
  hClusterTime->SetAxisRange(0.,30.,"X");
  hClusterTime->SetTitleOffset(1.5,"Y");
  hClusterTime->Draw("colz");
  
  ccalo2->cd(1);
  
  TH2F* hCellActivity  = (TH2F*) GetHisto("QA_hGridCells");
  hCellActivity->SetTitle("Hits per cell (E > 0.2 GeV)");
  hCellActivity->SetTitleOffset(1.5,"Y");
  hCellActivity->Draw("colz");
  
  ccalo2->cd(2);
  
  TH2F* hCellActivity  = (TH2F*) GetHisto("QA_hGridCells");
  TH2F* hCellActivityE = (TH2F*) GetHisto("QA_hGridCellsE");
  hCellActivityE->SetTitle("Mean energy per cell (E > 0.2 GeV)");
  hCellActivityE->Divide(hCellActivity);
  hCellActivityE->SetTitleOffset(1.5,"Y");
  hCellActivityE->Draw("colz");
  
  ccalo2->cd(4);
  //gPad->SetLogz();
  TH2F* hClusterActivity  = (TH2F*) GetHisto("AnaPhoton_hEtaPhi");
  hClusterActivity->SetTitle("Clusters activity (E > 0.5 GeV)");
  hClusterActivity->SetTitleOffset(1.5,"Y");
  hClusterActivity->Draw("colz");
  
  ccalo2->Print(Form("%s_CaloHisto2.eps",histoTag.Data()));
}

///
/// Plot basic hybrid tracks histograms.
///
//______________________________________
void TrackQA()
{
  TCanvas * ctrack = new TCanvas(Form("TrackHisto_%s",histoTag.Data()),"",1000,500);
  ctrack->Divide(2,1);
  
  ctrack->cd(1);
  //gPad->SetLogz();
  TH2F * hTrackEtaPhi = (TH2F*) GetHisto("AnaHadrons_hEtaPhiNegative");
  hTrackEtaPhi ->Add((TH2F*) GetHisto("AnaHadrons_hEtaPhiNegative"));
  hTrackEtaPhi ->SetAxisRange(-0.9,0.9,"X");
  hTrackEtaPhi ->SetTitle("Hybrid tracks #eta vs #phi (p_{T} > 0.2 GeV)");
  hTrackEtaPhi ->Draw("colz");
  
  ctrack->cd(2);
  //gPad->SetLogy();
  TH2F * hTrackEtaPhiSPD   = (TH2F*) GetHisto("AnaHadrons_hEtaPhiSPDRefitPt02");
  TH2F * hTrackEtaPhiNoSPD = (TH2F*) GetHisto("AnaHadrons_hEtaPhiNoSPDRefitPt02");
  
  TH1F* hPhiSPD   = (TH1F*)hTrackEtaPhiSPD  ->ProjectionY("hTrackPhiSPD"  ,0,1000);
  TH1F* hPhiNoSPD = (TH1F*)hTrackEtaPhiNoSPD->ProjectionY("hTrackPhiNoSPD",0,1000);
  //TH1F* hPhi      = (TH1F*)hTrackEtaPhi     ->ProjectionY("hTrackPhi"     ,0,1000);
  TH1F* hPhi      = hPhiSPD->Clone("hTrackPhi");
  hPhi->Add(hPhiNoSPD);
  hPhi     ->SetTitle("Hybrid track in #phi, composition, p_{T} > 0.2 GeV");
  hPhi     ->SetLineColor(1);
  hPhiSPD  ->SetLineColor(2);
  hPhiNoSPD->SetLineColor(4);
  
  hPhi     ->SetMinimum(1);
  hPhi     ->SetMaximum(hPhi->GetMaximum()*1.2);
  
  hPhi     ->Draw("H");
  hPhiSPD  ->Draw("Hsame");
  hPhiNoSPD->Draw("Hsame");
  
  TLegend l(0.2,0.75,0.4,0.89);
  l.SetTextSize(0.04);
  l.AddEntry(hPhi,"Sum","L");
  l.AddEntry(hPhiSPD  ,"SPD+Refit","L");
  l.AddEntry(hPhiNoSPD,"No SPD+Refit","L");
  l.SetBorderSize(0);
  l.SetFillColor(0);
  l.Draw();
//  ctrack->cd(3);
//  gPad->SetLogz();
//  
//  TH2F* hPtDCAxy = (TH2F*) GetHisto("AnaHadrons_hPtDCAxy");
//  hPtDCAxy->SetAxisRange(-1,1,"Y");
//  hPtDCAxy->SetAxisRange(0,30,"X");
//  hPtDCAxy->Draw("colz");
//  
//  ctrack->cd(4);
//  gPad->SetLogz();
//  
//  TH2F* hPtDCAz = (TH2F*) GetHisto("AnaHadrons_hPtDCAz");
//  hPtDCAz->SetAxisRange(-1,1,"Y");
//  hPtDCAz->SetAxisRange(0,30,"X");
//  hPtDCAz->Draw("colz");
  
  ctrack->Print(Form("%s_TrackHisto.eps",histoTag.Data()));
}

///
/// Plot basic invariant mass QA
///
//_____________________________
void Pi0QA()
{
  TCanvas * cpi0 = new TCanvas(Form("Pi0Histo_%s",histoTag.Data()),"",500,500);
  cpi0->Divide(2,2);
  
  TH2F* hMassE[10];
  TH2F* hMixMassE[10];
  for(Int_t icen = 0; icen < 10; icen++)
  {
    hMassE   [icen] = (TH2F*) GetHisto(Form("AnaPi0_hRe_cen%d_pidbit0_asy0_dist1",icen));
    hMixMassE[icen] = (TH2F*) GetHisto(Form("AnaPi0_hMi_cen%d_pidbit0_asy0_dist1",icen));
  }
  
  // 2D Invariant mass vs E, in PbPb from 60 to 100 %, all in pp
  cpi0->cd(1);
  gPad->SetLogz();
  TH2F* h2DMass;
  
  if(hMassE[1]) // Plot centrality from 60 to 100%
  {
    h2DMass = (TH2F*) hMassE[6]->Clone("h2DMass");
    for(Int_t icen = 7; icen < 10; icen++) h2DMass->Add(hMassE[icen]);
    h2DMass->SetTitle("Invariant mass vs pair E, Cen: 60-100%");
  }
  else
  {
    h2DMass = (TH2F*) hMassE[0]->Clone("hMassProj");
    h2DMass->SetTitle("Invariant mass vs cluster pair E");
  }
  
  h2DMass->SetTitleOffset(1.6,"Y");
  h2DMass->SetAxisRange(0.0,0.7,"Y");
  h2DMass->SetAxisRange(0,30,"X");
  h2DMass->Draw("colz");
  
  // Pi0 Invariant mass projection, in PbPb 6 centrality bins from 0 to 50%, all in pp
  cpi0->cd(2);
  TH1F* hMass   [10];
  TH1F* hMix    [10];
  TH1F* hMassEta[10];
  TH1F* hMassPi0[10];
  
  //Init to 0
  for(Int_t icen=0; icen<10; icen++ )
  {
    hMass   [icen] = 0;
    hMix    [icen] = 0;
    hMassEta[icen] = 0;
    hMassPi0[icen] = 0;
  }
  
  TH1F * hX = (TH1F*) hMassE[0]->ProjectionX("hEPairCen0",0,10000);
  Int_t binmin = hX->FindBin(2);  // Project histo from 2 GeV pairs
  Int_t binmax = hX->FindBin(10); // Project histo up to 10 GeV pairs
  Float_t maxPi0 = 0;
  Float_t maxEta = 0;
  for(Int_t icen = 0; icen < 6; icen++)
  {
    if(!hMassE[icen]) continue;

    hMass[icen] = (TH1F*) hMassE   [icen]->ProjectionY(Form("hMassCen%d",icen),binmin,binmax);
    hMix [icen] = (TH1F*) hMixMassE[icen]->ProjectionY(Form("hMixCen%d" ,icen),binmin,binmax);
    hMass[icen]->Sumw2();
    hMix [icen]->Sumw2();
    
    hMassPi0[icen] = (TH1F*) hMass[icen]->Clone(Form("hMassPi0Cen%d",icen));
    hMassEta[icen] = (TH1F*) hMass[icen]->Clone(Form("hMassEtaCen%d",icen));
    
    hMassPi0[icen]->Divide(hMix[icen]);
    hMassPi0[icen]->Fit("pol0","Q","",0.25,0.35);
    Float_t scale = 1;
    if(hMassPi0[icen]->GetFunction("pol0")) scale = hMassPi0[icen]->GetFunction("pol0")->GetParameter(0);
    //printf("Scale factor %f for cen %d\n",scale,icen);
    hMassPi0[icen]->Scale(1./scale);
    hMassPi0[icen]->SetMarkerStyle(24);
    hMassPi0[icen]->SetMarkerColor(color[icen]);
    hMassPi0[icen]->SetLineColor(color[icen]);
    hMassPi0[icen]->SetAxisRange(0.04,0.24);
    hMassPi0[icen]->SetMarkerSize(0.5);
    
    hMassEta[icen]->Rebin(4);
    hMix    [icen]->Rebin(4);
    hMassEta[icen]->Divide(hMix[icen]);
    hMassEta[icen]->SetMarkerStyle(25);
    hMassEta[icen]->SetMarkerColor(color[icen]);
    hMassEta[icen]->SetLineColor(color[icen]);
    hMassEta[icen]->SetAxisRange(0.4,0.9);
    hMassEta[icen]->SetMarkerSize(0.5);
    hMassEta[icen]->Scale(1./scale);
    
    if(maxEta < hMassEta[icen]->GetMaximum()) maxEta = hMassEta[icen]->GetMaximum();
    if(maxPi0 < hMassPi0[icen]->GetMaximum()) maxPi0 = hMassPi0[icen]->GetMaximum();
  }

  //gPad->SetLogy();
  //gPad->SetGridy();
  hMassPi0[0]->SetMinimum(0.8);
  hMassPi0[0]->SetTitleOffset(1.6,"Y");
  hMassPi0[0]->SetYTitle("Real / Mixed");
  hMassPi0[0]->SetTitle("#pi^{0} peak, 2 < E_{pair}< 10 GeV");
  hMassPi0[0]->Draw();
  
  if(hMass[1]) // PbPb
  {
    hMassPi0[0]->SetMaximum(maxPi0*1.2);
    hMassPi0[5]->Draw("Hsame");
    hMassPi0[4]->Draw("Hsame");
    hMassPi0[3]->Draw("Hsame");
    hMassPi0[2]->Draw("Hsame");
    hMassPi0[1]->Draw("Hsame");
    hMassPi0[0]->Draw("Hsame");
    //hMass[6]->Draw("Hsame");
    //hMass[7]->Draw("same");
    //hMass[8]->Draw("same");
    //hMass[9]->Draw("same");
    
    TLegend l(0.12,0.6,0.4,0.85);
    l.SetTextSize(0.04);
    l.AddEntry(hMassPi0[0],"0-10%","P");
    l.AddEntry(hMassPi0[1],"10-20%","P");
    l.AddEntry(hMassPi0[2],"20-30%","P");
    l.AddEntry(hMassPi0[3],"30-40%","P");
    l.AddEntry(hMassPi0[4],"40-70%","P");
    l.AddEntry(hMassPi0[5],"50-60%","P");
    l.SetBorderSize(0);
    l.SetFillColor(0);
    l.Draw();
  }

  TLine l1(0.04,1,0.24,1);
  l1.Draw("same");
  
  // Pi0 invariant mass per EMCal super module
  cpi0->cd(3);
  
  TH1F* hSM   [10];
  TH1F* hMixSM[10];
  binmin = hX->FindBin(4);  // Project histo from 3 GeV pairs
  binmax = hX->FindBin(20); // Project histo up to 20 GeV pairs
  Float_t maxSM = 0;

  for(Int_t ism = 0; ism < 10; ism++)
  {
    TH2F* hTmpSM = (TH2F*) GetHisto(Form("AnaPi0_hReMod_%d",ism));
    if(!hTmpSM) hTmpSM = (TH2F*) GetHisto(Form("QA_hIM_Mod%d",ism));
    
    hSM[ism] = (TH1F*) hTmpSM->ProjectionY(Form("hMassSM%d",ism),binmin,binmax);
    hSM[ism]->Sumw2();
    hSM[ism]->SetMarkerStyle(26);
    hSM[ism]->Rebin(2);
    //hSM[ism]->Scale(1./hSM[ism]->Integral(0,10000));
    hSM[ism]->SetMarkerColor(color[ism]);
    hSM[ism]->SetLineColor(color[ism]);
    hSM[ism]->SetMarkerSize(0.5);

    TH2F* hTmpMixSM = (TH2F*) GetHisto(Form("AnaPi0_hMiMod_%d",ism));
    if(hTmpMixSM)
    {
      hMixSM[ism] = (TH1F*) hTmpMixSM->ProjectionY(Form("hMassMixSM%d",ism),binmin,binmax);
      hMixSM[ism]->Sumw2();
      hMixSM[ism]->Rebin(2);
      hSM[ism]->Divide(hMixSM[ism]);
      hSM[ism]->Fit("pol0","Q","",0.25,0.35);
      Float_t scale = 1;
      if(hSM[ism]->GetFunction("pol0")) scale = hSM[ism]->GetFunction("pol0")->GetParameter(0);
      //printf("Scale factor %f for cen %d\n",scale,icen);
      hSM[ism]->Scale(1./scale);
      hSM[ism]->SetYTitle("Real / Mixed");
    }
    
    if(maxSM < hSM[ism]->GetMaximum()) maxSM = hSM[ism]->GetMaximum();
  }
  
  hSM[0]->SetTitle("#pi^{0} peak in Modules, 4 < E_{pair}< 10 GeV");
  hSM[0]->SetTitleOffset(1.6,"Y");
  hSM[0]->SetAxisRange(0.04,0.24);
  hSM[0]->SetMaximum(maxSM*1.2);
  hSM[0]->SetMinimum(0.8);

  hSM[0]->Draw("H");
  TLegend lsm(0.12,0.5,0.35,0.85);
  lsm.SetTextSize(0.04);
  lsm.AddEntry(hSM[0],Form("Mod %d",0),"P");
  
  for(Int_t ism = 1; ism < 10; ism++)
  {
    hSM[ism]->Draw("Hsame");
    lsm.AddEntry(hSM[ism],Form("Mod %d",ism),"P");
  }
  
  lsm.SetBorderSize(0);
  lsm.SetFillColor(0);
  lsm.Draw();
  
  l1.Draw("same");
  
  // Pi0 Invariant mass projection, in PbPb 6 centrality bins from 0 to 50%, all in pp
  cpi0->cd(4);
  
  //gPad->SetLogy();
  //gPad->SetGridy();
  hMassEta[0]->SetMinimum(0.8);
  hMassEta[0]->SetTitleOffset(1.6,"Y");
  hMassEta[0]->SetYTitle("Real / Mixed");
  hMassEta[0]->SetTitle("#eta peak, 2 < E_{pair}< 10 GeV");
  hMassEta[0]->Draw("H");
  
  if(hMass[1]) // PbPb
  {
    hMassEta[0]->SetMaximum(maxEta*1.2);
    hMassEta[5]->Draw("Hsame");
    hMassEta[4]->Draw("Hsame");
    hMassEta[3]->Draw("Hsame");
    hMassEta[2]->Draw("Hsame");
    hMassEta[1]->Draw("Hsame");
    hMassEta[0]->Draw("Hsame");
    
    
    TLegend l2(0.12,0.6,0.4,0.85);
    l2.SetTextSize(0.04);
    l2.AddEntry(hMassEta[0],"0-10%","P");
    l2.AddEntry(hMassEta[1],"10-20%","P");
    l2.AddEntry(hMassEta[2],"20-30%","P");
    l2.AddEntry(hMassEta[3],"30-40%","P");
    l2.AddEntry(hMassEta[4],"40-70%","P");
    l2.AddEntry(hMassEta[5],"50-60%","P");
    l2.SetBorderSize(0);
    l2.SetFillColor(0);
    l2.Draw();
  }

  cpi0->Print(Form("%s_Pi0Histo.eps",histoTag.Data()));
}

///
/// Plot basic cluster-track correlation histograms.
///
//__________________________________________________
void CorrelationQA()
{
  TCanvas * cCorrelation = new TCanvas(Form("CorrelationHisto_%s",histoTag.Data()),"",1000,500);
  cCorrelation->Divide(2,1);
  
  Float_t minClusterE = 8;
  Float_t assocBins[] = {0.5,2.,5.,10.,20.};
  Int_t nAssocBins = 4;
  
  TH1F * hTrigger = (TH1F*) GetHisto("AnaPhotonHadronCorr_hPtTrigger");
  Int_t minClusterEBin = hTrigger->FindBin(minClusterE);
  Float_t nTrig = hTrigger->Integral(minClusterE,100000);
  
  //Azimuthal correlation
  cCorrelation->cd(1);
  gPad->SetLogy();
  TH1F* hDeltaPhi[4];
  
  TLegend l(0.35,0.6,0.83,0.85);
  l.SetHeader(Form("p_{T,T} > %2.1f GeV/c",minClusterE));
  l.SetTextSize(0.04);
  l.SetBorderSize(0);
  l.SetFillColor(0);

  for(Int_t ibin = 0; ibin < nAssocBins; ibin++ )
  {
    TH2F* hDeltaPhiE = (TH2F*) GetHisto(Form("AnaPhotonHadronCorr_hDeltaPhiPtAssocPt%2.1f_%2.1f",assocBins[ibin],assocBins[ibin+1]));
    hDeltaPhi[ibin]  = (TH1F*) hDeltaPhiE->ProjectionY(Form("DeltaPhi%2.1f",assocBins[ibin]),minClusterEBin,10000);
    hDeltaPhi[ibin]->Sumw2();
    hDeltaPhi[ibin]->Rebin(2);
    hDeltaPhi[ibin]->Scale(1./nTrig);

    hDeltaPhi[ibin]->Fit("pol0","Q","",1,2);
    Float_t scale = 1;
    if(hDeltaPhi[ibin]->GetFunction("pol0"))
    {
      scale = hDeltaPhi[ibin]->GetFunction("pol0")->GetParameter(0);
      hDeltaPhi[ibin]->GetFunction("pol0")->SetRange(6,7); // move from plot
    }
    hDeltaPhi[ibin]->Scale(1./scale);
    //printf("ibin %d, scale %f\n",ibin,scale);

    hDeltaPhi[ibin]->SetAxisRange(-1.6,4.7);

    hDeltaPhi[ibin]->SetMarkerStyle(24);
    hDeltaPhi[ibin]->SetMarkerColor(color[ibin]);
    hDeltaPhi[ibin]->SetLineColor(color[ibin]);
    hDeltaPhi[ibin]->SetTitleOffset(1.5,"Y");
    hDeltaPhi[ibin]->SetYTitle("N_{pairs} / N_{trig} / ZYAM");
    hDeltaPhi[ibin]->SetTitle("#gamma (#lambda_{0}^{2} < 0.4, neutral cluster) trigger");
    
    l.AddEntry(hDeltaPhi[ibin],Form("%2.1f< p_{T,A}< %2.1f GeV/c",assocBins[ibin],assocBins[ibin+1]),"P");
  }
  
  hDeltaPhi[2]->SetMaximum(hDeltaPhi[2]->GetMaximum()*10);
  hDeltaPhi[2]->SetMinimum(0.8);
  
   hDeltaPhi[2]->Draw("H");
   hDeltaPhi[1]->Draw("Hsame");
   hDeltaPhi[3]->Draw("Hsame");
   hDeltaPhi[0]->Draw("Hsame");
  
  l.Draw("same");
  
  
  // xE correlation
  cCorrelation->cd(2);
  gPad->SetLogy();
  
  TLegend l2(0.35,0.6,0.83,0.85);
  l2.SetHeader(Form("p_{T,T} > %2.1f GeV/c",minClusterE));
  l2.SetTextSize(0.04);
  l2.SetBorderSize(0);
  l2.SetFillColor(0);
  
  TH2F* hEXE   = (TH2F*) GetHisto("AnaPhotonHadronCorr_hXECharged");
  TH2F* hEXEUE = (TH2F*) GetHisto("AnaPhotonHadronCorr_hXEUeCharged");
  
  TH1F* hXE  = (TH1F*) hEXE->ProjectionY(Form("XE%2.1fGeV",minClusterE),minClusterEBin,10000);
  hXE->Sumw2();
  hXE->Rebin(2);
  hXE->Scale(1./nTrig);
  hXE->SetAxisRange(0,1);
  hXE->SetMarkerStyle(24);
  hXE->SetMarkerColor(1);
  hXE->SetLineColor(1);
  hXE->SetTitleOffset(1.5,"Y");
  hXE->SetYTitle("N_{pairs} / N_{trig}");
  hXE->SetTitle("#gamma (#lambda_{0}^{2} < 0.4, neutral cluster) trigger");
  l2.AddEntry(hXE,"raw x_{E}","P");
  hXE->Draw();

  TH1F* hXEUE  = (TH1F*) hEXEUE->ProjectionY(Form("XEUE%2.1fGeV",minClusterE),minClusterEBin,10000);
  hXEUE->Sumw2();
  hXEUE->Rebin(2);
  hXEUE->Scale(1./nTrig);
  hXEUE->SetAxisRange(0,1);
  hXEUE->SetMarkerStyle(25);
  hXEUE->SetMarkerColor(2);
  hXEUE->SetLineColor(2);
  l2.AddEntry(hXEUE,"raw Und. Event x_{E}","P");
  hXEUE->Draw("same");
  
  l2.Draw("same");

  cCorrelation->Print(Form("%s_CorrelationHisto.eps",histoTag.Data()));
}

///
/// Plot basic generated particle distribution histograms.
///
//________________________________________________________
void MCQA()
{
  TH2F* h2ClusterPho = (TH2F*) GetHisto("QA_hRecoMCE_Photon_Match0");    // not track-matched
  TH2F* h2ClusterPi0 = (TH2F*) GetHisto("QA_hRecoMCE_Pi0_Match0");       // not track-matched
  TH2F* h2ClusterEta = (TH2F*) GetHisto("QA_hRecoMCE_Eta_Match0");       // not track-matched
  TH2F* h2ClusterEle = (TH2F*) GetHisto("QA_hRecoMCE_Electron_Match1");  // Track-matched
  
  if(!h2ClusterPho) return;
  
//  TH1F* hPrimPho = (TH1F*) GetHisto("QA_hGenMCAccE_Photon");
//  TH1F* hPrimPi0 = (TH1F*) GetHisto("QA_hGenMCAccE_Pi0");
//  TH1F* hPrimEta = (TH1F*) GetHisto("QA_hGenMCAccE_Eta");

  TH1F* hPrimPho = (TH1F*) GetHisto("AnaPhoton_hPtPrim_MCPhoton");
  TH1F* hPrimPi0 = (TH1F*) GetHisto("AnaPi0_hPrimPi0Pt");
  TH1F* hPrimEta = (TH1F*) GetHisto("AnaPi0_hPrimEtaPt");
  
  TCanvas * cmc = new TCanvas(Form("MCHisto_%s",histoTag.Data()),"",1000,1000);
  cmc->Divide(2,2);
  
  cmc->cd(1);
  gPad->SetLogy();
  
  TH1F* hClusterPho = (TH1F*) h2ClusterPho->ProjectionX("ClusterPho",0,1000);
  TH1F* hClusterPi0 = (TH1F*) h2ClusterPi0->ProjectionX("ClusterPi0",0,1000);
  TH1F* hClusterEta = (TH1F*) h2ClusterEta->ProjectionX("ClusterEta",0,1000);
  
  hClusterPho->SetTitle("Cluster origin spectra, primary spectra in Calo acceptance");
  hClusterPho->Sumw2();
  hClusterPho->SetMarkerColor(1);
  hClusterPho->SetMarkerStyle(20);
  hClusterPho->SetAxisRange(0.,50.,"X");
  //hClusterPho->SetXTitle("E_{rec,gen} (GeV)");
  hClusterPho->SetXTitle("E_{rec}, p_{T,gen} (GeV)");
  hClusterPho->Draw("");

  hClusterPi0->Sumw2();
  hClusterPi0->SetMarkerColor(4);
  hClusterPi0->SetMarkerStyle(21);
  hClusterPi0->Draw("same");

  hClusterEta->Sumw2();
  hClusterEta->SetMarkerColor(2);
  hClusterEta->SetMarkerStyle(22);
  hClusterEta->Draw("same");

  hPrimPho->Sumw2();
  hPrimPho->SetMarkerColor(1);
  hPrimPho->SetMarkerStyle(24);
  hPrimPho->Draw("same");
  
  hPrimPi0->Sumw2();
  hPrimPi0->SetMarkerColor(4);
  hPrimPi0->SetMarkerStyle(25);
  hPrimPi0->Draw("same");
  
  hPrimEta->Sumw2();
  hPrimEta->SetMarkerColor(2);
  hPrimEta->SetMarkerStyle(26);
  hPrimEta->Draw("same");
  
  TLegend l(0.45,0.6,0.83,0.89);
  l.SetTextSize(0.04);
  l.AddEntry(hClusterPho,"#gamma cluster","P");
  l.AddEntry(hClusterPi0,"#pi^{0} (merged) cluster","P");
  l.AddEntry(hClusterEta,"#eta (merged) cluster","P");
  l.AddEntry(hPrimPho,"#gamma generated","P");
  l.AddEntry(hPrimPi0,"#pi^{0} generated","P");
  l.AddEntry(hPrimEta,"#eta generated","P");
  l.SetBorderSize(0);
  l.SetFillColor(0);
  l.Draw();
  
  
  cmc->cd(2);
  gPad->SetLogy();
  TH1F* hRatPho = (TH1F*) hClusterPho->Clone("hGenRecoPho");
  TH1F* hRatPi0 = (TH1F*) hClusterPi0->Clone("hGenRecoPi0");
  TH1F* hRatEta = (TH1F*) hClusterEta->Clone("hGenRecoEta");
  
  hRatPho->Divide(hPrimPho);
  hRatPi0->Divide(hPrimPi0);
  hRatEta->Divide(hPrimEta);
  
  hRatPho->SetTitle("Reconstructed cluster / Generated particle in Calo acc.");
  hRatPho->SetYTitle("Ratio");
  hRatPho->SetXTitle("E(GeV)");
  hRatPho->SetMinimum(1e-3);
  hRatPho->SetMaximum(10);
  hRatPho->Draw("");
  hRatPi0->Draw("same");
  hRatEta->Draw("same");

  TLegend l2(0.15,0.7,0.3,0.85);
  l2.SetTextSize(0.04);
  l2.AddEntry(hRatPho,"#gamma","P");
  l2.AddEntry(hRatPi0,"#pi^{0}","P");
  l2.AddEntry(hRatEta,"#eta","P");
  l2.SetBorderSize(0);
  l2.SetFillColor(0);
  l2.Draw();

  cmc->cd(3);
  //gPad->SetLogy();

  TH2F* h2PrimPhoPhi = (TH2F*) GetHisto("AnaPhoton_hPhiPrim_MCPhoton");
  TH2F* h2PrimPi0Phi = (TH2F*) GetHisto("AnaPi0_hPrimPi0Phi");
  TH2F* h2PrimEtaPhi = (TH2F*) GetHisto("AnaPi0_hPrimEtaPhi");
  
  Int_t binMin = hPrimPho->FindBin(3);
  
  TH1F* hPrimPhoPhi = (TH1F*) h2PrimPhoPhi->ProjectionY("PrimPhoPhi",binMin,1000);
  TH1F* hPrimPi0Phi = (TH1F*) h2PrimPi0Phi->ProjectionY("PrimPi0Phi",binMin,1000);
  TH1F* hPrimEtaPhi = (TH1F*) h2PrimEtaPhi->ProjectionY("PrimEtaPhi",binMin,1000);

  hPrimPhoPhi->Sumw2();
  hPrimPi0Phi->Sumw2();
  hPrimEtaPhi->Sumw2();

  hPrimPhoPhi->Scale(1./hPrimPhoPhi->Integral(0,1000));
  hPrimPi0Phi->Scale(1./hPrimPi0Phi->Integral(0,1000));
  hPrimEtaPhi->Scale(1./hPrimEtaPhi->Integral(0,1000));

  Float_t maxPhi = hPrimPhoPhi->GetMaximum();
  if(maxPhi < hPrimPi0Phi->GetMaximum()) maxPhi =  hPrimPi0Phi->GetMaximum();
  if(maxPhi < hPrimEtaPhi->GetMaximum()) maxPhi =  hPrimEtaPhi->GetMaximum();

  Float_t minPhi = hPrimPhoPhi->GetMinimum();
  if(minPhi > hPrimPi0Phi->GetMinimum()) minPhi =  hPrimPi0Phi->GetMinimum();
  if(minPhi > hPrimEtaPhi->GetMinimum()) minPhi =  hPrimEtaPhi->GetMinimum();

  hPrimPi0Phi->SetMaximum(maxPhi*1.1);
  hPrimPi0Phi->SetMinimum(minPhi);
  TGaxis::SetMaxDigits(3);

  hPrimPi0Phi->SetYTitle("1/total entries dN/d#phi");
  hPrimPi0Phi->SetTitle("Generated particles #phi for E > 3 GeV");
  hPrimPi0Phi->SetTitleOffset(1.6,"Y");
  hPrimPi0Phi->SetMarkerColor(4);
  hPrimPi0Phi->SetMarkerStyle(21);
  hPrimPi0Phi->Draw("");
  
  hPrimPhoPhi->SetMarkerColor(1);
  hPrimPhoPhi->SetMarkerStyle(20);
  Float_t scale = TMath::RadToDeg();
  ScaleXaxis(hPrimPhoPhi, TMath::RadToDeg());
  hPrimPhoPhi->Draw("same");

  hPrimEtaPhi->SetMarkerColor(2);
  hPrimEtaPhi->SetMarkerStyle(22);
  hPrimEtaPhi->Draw("same");

  cmc->cd(4);
  //gPad->SetLogy();
  
  TH2F* h2PrimPhoEta = (TH2F*) GetHisto("AnaPhoton_hYPrim_MCPhoton");
  TH2F* h2PrimPi0Eta = (TH2F*) GetHisto("AnaPi0_hPrimPi0Rapidity");
  TH2F* h2PrimEtaEta = (TH2F*) GetHisto("AnaPi0_hPrimEtaRapidity");

  h2PrimPhoEta->Sumw2();
  h2PrimEtaEta->Sumw2();
  h2PrimPi0Eta->Sumw2();
  
  Int_t binMin = hPrimPho->FindBin(3);
  
  TH1F* hPrimPhoEta = (TH1F*) h2PrimPhoEta->ProjectionY("PrimPhoEta",binMin,1000);
  TH1F* hPrimPi0Eta = (TH1F*) h2PrimPi0Eta->ProjectionY("PrimPi0Eta",binMin,1000);
  TH1F* hPrimEtaEta = (TH1F*) h2PrimEtaEta->ProjectionY("PrimEtaEta",binMin,1000);
  
  hPrimPhoEta->Scale(1./hPrimPhoEta->Integral(0,1000));
  hPrimPi0Eta->Scale(1./hPrimPi0Eta->Integral(0,1000));
  hPrimEtaEta->Scale(1./hPrimEtaEta->Integral(0,1000));
  
  Float_t maxEta = hPrimPhoEta->GetMaximum();
  if(maxEta < hPrimPi0Eta->GetMaximum()) maxEta =  hPrimPi0Eta->GetMaximum();
  if(maxEta < hPrimEtaEta->GetMaximum()) maxEta =  hPrimEtaEta->GetMaximum();
  
  Float_t minEta = hPrimPhoEta->GetMinimum();
  if(minEta > hPrimPi0Eta->GetMinimum()) minEta =  hPrimPi0Eta->GetMinimum();
  if(minEta > hPrimEtaEta->GetMinimum()) minEta =  hPrimEtaEta->GetMinimum();
  
  hPrimPi0Eta->SetMaximum(maxEta*1.1);
  hPrimPi0Eta->SetMinimum(minEta);
  TGaxis::SetMaxDigits(3);
  
  hPrimPi0Eta->SetYTitle("1/total entries dN/d#eta");
  hPrimPi0Eta->SetTitle("Generated particles #eta for E > 3 GeV");
  hPrimPi0Eta->SetTitleOffset(1.6,"Y");
  hPrimPi0Eta->SetMarkerColor(4);
  hPrimPi0Eta->SetMarkerStyle(21);
  hPrimPi0Eta->Draw("");
  
  hPrimPhoEta->SetMarkerColor(1);
  hPrimPhoEta->SetMarkerStyle(20);
  Float_t scale = TMath::RadToDeg();
  hPrimPhoEta->Draw("same");
  
  hPrimEtaEta->SetMarkerColor(2);
  hPrimEtaEta->SetMarkerStyle(22);
  hPrimEtaEta->Draw("same");
  
  cmc->Print(Form("%s_MCHisto.eps",histoTag.Data()));
}

///
/// Open the file and list containing the histograms
///
//____________________________________________________________________
void GetFileAndList(TString fileName, TString listName, Bool_t export)
{
  file  = new TFile(fileName,"read");
  
  TDirectory * dir = (TDirectory*) file->Get(listName);
  if(dir)
  {
    list = (TList*) dir->Get(listName);
    if(export)
    {
      TFile * outputFile = new TFile("AnalysisResultsList.root","RECREATE");
      list->Write();
      outputFile->Close();
    }
  }
}

///
/// Check if the list is available,
/// if not get the histo directly from file
///
//___________________________________
TObject * GetHisto(TString histoName)
{
  if(list) return list->FindObject(histoName);
  else     return file->Get       (histoName);
}

///
/// Scale axis by a constant factor
/// used just to scale degrees to rad in a single histogram in the MC case
///
//___________________________________________________
void ScaleAxis(TAxis *a, Double_t scale)
{
  if (!a) return; // just a precaution
  if (a->GetXbins()->GetSize())
  {
    // an axis with variable bins
    // note: bins must remain in increasing order, hence the "Scale"
    // function must be strictly (monotonically) increasing
    TArrayD X(*(a->GetXbins()));
    for(Int_t i = 0; i < X.GetSize(); i++) X[i] = scale*X[i];
    a->Set((X.GetSize() - 1), X.GetArray()); // new Xbins
  }
  else
  {
    // an axis with fix bins
    // note: we modify Xmin and Xmax only, hence the "Scale" function
    // must be linear (and Xmax must remain greater than Xmin)
    a->Set(a->GetNbins(),
           -           scale*a->GetXmin(), // new Xmin
           -           scale*a->GetXmax()); // new Xmax
  }
  return;
}

///
/// Scale x axis by a constant factor
/// used just to scale degrees to rad in a single histogram in the MC case
///
//___________________________________________________
void ScaleXaxis(TH1 *h, Double_t scale)
{
  if (!h) return; // just a precaution
  ScaleAxis(h->GetXaxis(), scale);
  return;
}



