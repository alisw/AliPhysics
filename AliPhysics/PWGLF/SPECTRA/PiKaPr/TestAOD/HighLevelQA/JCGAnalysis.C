/////////////////////////////////////////////////////
// JCGAnalysis.C                                   //
// based on MainAnalysis.C                         //
// uses a new method of calculating the efficiency //
/////////////////////////////////////////////////////

class AliPID;

// global constants
const Int_t nCharge = 2;
const Int_t nPart = 3;
TString Charge[nCharge] = {"Pos","Neg"};
TString Sign[nCharge] = {"Plus","Minus"};
TString Particle[nPart] = {"Pion","Kaon","Proton"};
Int_t Color[nPart] = {1,2,4};
Int_t Marker[nPart*nCharge] = {20,21,22,24,25,26};
Double_t Range[nPart] = {0.3,0.3,0.5}; // LowPt range for pi k p
enum ECharge_t {kPositive,kNegative,kNCharges};
TString Names[nCharge*nPart] = {"#pi^{+}",
				"K^{+}",
				"p",
				"#pi^{-}",
				"K^{-}",
				"#bar{p}"};

void JCGAnalysis()
{
  // load libraries 'n such
  gSystem->Load("libCore");
  gSystem->Load("libPhysics");
  gSystem->Load("libTree");
  gSystem->Load("libMatrix");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTender");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gStyle->SetOptStat(000000);
  gStyle->SetPalette(1);

  // Set masses
  Double_t mass[3];
  mass[0]   = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  mass[1]   = TDatabasePDG::Instance()->GetParticle("K+")->Mass();
  mass[2] = TDatabasePDG::Instance()->GetParticle("proton")->Mass();
  
  // cuts 'n such
  //                            0    1    2    3
  Double_t CentCutMin[4]= {     0,  30,  30,  30};
  Double_t CentCutMax[4]= {     5,  40,  40,  40};
  Double_t QvecCutMin[4]={      0,   0,   0, 1.5};
  Double_t QvecCutMax[4]={    100, 100, 0.4, 100};
  Double_t EtaMin[4]={       -0.8,-0.8,-0.8,-0.8};
  Double_t EtaMax[4]={        0.8, 0.8, 0.8, 0.8};
  Double_t Nsigmapid=3.;
  Double_t pt=10.;
  Double_t p=10.;
  Double_t y=.5;
  Double_t ptTofMatch=.6;
  UInt_t trkbit=1024;
  UInt_t trkbitQVector=1;
  Bool_t UseCentPatchAOD049=kFALSE;
  Double_t DCA=100000;
  UInt_t minNclsTPC=70;
  Int_t nrebin=0;
  TString opt="";
  Int_t icut=1; // chooses which cut to use
  if(icut==0)Int_t ibinToCompare=0;
  else Int_t ibinToCompare=4;
  
  TString sname_data=Form("OutputAODSpectraTask_Data_Cent%.0fto%.0f_QVec%.1fto%.1f_Eta%.1fto%.1f_%.1fSigmaPID_TrBit%d",CentCutMin[icut],CentCutMax[icut],QvecCutMin[icut],QvecCutMax[icut],EtaMin[icut],EtaMax[icut],Nsigmapid,trkbit);
  TString sname_data_1sig=Form("OutputAODSpectraTask_Data_Cent%.0fto%.0f_QVec%.1fto%.1f_Eta%.1fto%.1f_1.0SigmaPID_TrBit%d",CentCutMin[icut],CentCutMax[icut],QvecCutMin[icut],QvecCutMax[icut],EtaMin[icut],EtaMax[icut],trkbit);
  // for MC we do not use the cut on Qvector
  TString sname_mc=Form("OutputAODSpectraTask_MC_Cent%.0fto%.0f_QVec0.0to100.0_Eta%.1fto%.1f_%.1fSigmaPID_TrBit%d",CentCutMin[icut],CentCutMax[icut],EtaMin[icut],EtaMax[icut],Nsigmapid,trkbit);
  TString sname_mc_5sig=Form("OutputAODSpectraTask_MC_Cent%.0fto%.0f_QVec0.0to100.0_Eta%.1fto%.1f_5.0SigmaPID_TrBit%d",CentCutMin[icut],CentCutMax[icut],EtaMin[icut],EtaMax[icut],trkbit);
  TString sname_mc_1sig=Form("OutputAODSpectraTask_MC_Cent%.0fto%.0f_QVec0.0to100.0_Eta%.1fto%.1f_1.0SigmaPID_TrBit%d",CentCutMin[icut],CentCutMax[icut],EtaMin[icut],EtaMax[icut],trkbit);

  // input and output files
  TString fold="No_Outliers";
  TString dataFile = Form("output/%s/AnResDATATrain12_NoOutliers.root",fold.Data());
  TString mcFile = Form("output/%s/AnResMCTrain27_NoOutliers.root",fold.Data());
  TFile * fout=new TFile(Form("results/%s/Res_%s_%s.root",fold.Data(),sname_data.Data(),fold.Data()),"RECREATE");

  // Open root MC file and get classes
  TFile *f_mc = TFile::Open(mcFile.Data());
  TDirectoryFile *Dir_mc=(TDirectoryFile*)f_mc->Get(Form("%s",sname_mc.Data()));
  AliSpectraAODHistoManager* hman_mc = (AliSpectraAODHistoManager*) Dir_mc->Get("SpectraHistos");
  AliSpectraAODEventCuts* ecuts_mc = (AliSpectraAODEventCuts*) Dir_mc->Get("Event Cuts");
  AliSpectraAODTrackCuts* tcuts_mc = (AliSpectraAODTrackCuts*) Dir_mc->Get("Track Cuts");
  TDirectoryFile *Dir_mc_5sig=(TDirectoryFile*)f_mc->Get(Form("%s",sname_mc_5sig.Data()));
  AliSpectraAODHistoManager* hman_mc_5sig = (AliSpectraAODHistoManager*) Dir_mc_5sig->Get("SpectraHistos");
  AliSpectraAODEventCuts* ecuts_mc_5sig = (AliSpectraAODEventCuts*) Dir_mc_5sig->Get("Event Cuts");
  AliSpectraAODTrackCuts* tcuts_mc_5sig = (AliSpectraAODTrackCuts*) Dir_mc_5sig->Get("Track Cuts");
  TDirectoryFile *Dir_mc_1sig=(TDirectoryFile*)f_mc->Get(Form("%s",sname_mc_1sig.Data()));
  AliSpectraAODHistoManager* hman_mc_1sig = (AliSpectraAODHistoManager*) Dir_mc_1sig->Get("SpectraHistos");
  AliSpectraAODEventCuts* ecuts_mc_1sig = (AliSpectraAODEventCuts*) Dir_mc_1sig->Get("Event Cuts");
  AliSpectraAODTrackCuts* tcuts_mc_1sig = (AliSpectraAODTrackCuts*) Dir_mc_1sig->Get("Track Cuts");
  // proceed likewise for data:
  // Open root DATA file and get classes
  TFile *f_data = TFile::Open(dataFile.Data());
  TDirectoryFile *Dir_data=(TDirectoryFile*)f_data->Get(Form("%s",sname_data.Data()));
  AliSpectraAODHistoManager* hman_data = (AliSpectraAODHistoManager*) Dir_data->Get("SpectraHistos");
  AliSpectraAODEventCuts* ecuts_data = (AliSpectraAODEventCuts*) Dir_data->Get("Event Cuts");
  AliSpectraAODTrackCuts* tcuts_data = (AliSpectraAODTrackCuts*) Dir_data->Get("Track Cuts");
  TDirectoryFile *Dir_data_1sig=(TDirectoryFile*)f_data->Get(Form("%s",sname_data_1sig.Data()));
  AliSpectraAODHistoManager* hman_data_1sig = (AliSpectraAODHistoManager*) Dir_data_1sig->Get("SpectraHistos");
  AliSpectraAODEventCuts* ecuts_data_1sig = (AliSpectraAODEventCuts*) Dir_data_1sig->Get("Event Cuts");
  AliSpectraAODTrackCuts* tcuts_data_1sig = (AliSpectraAODTrackCuts*) Dir_data_1sig->Get("Track Cuts");

  //---------------------------------------------
  // plot nsigma distributions for presentation -
  //---------------------------------------------
  // TPC
  TH2F * hTPCnsigPion = new TH2F(*(TH2F*)((TH2F*)hman_data->GetNSigHistogram("hHistNSigPionTPC"))->Clone());
  hTPCnsigPion->SetTitle("TPC NSigma Distribution for Pions;p_{T} (GeV/c);n#sigma");
  TCanvas * cTPCnsigPion = new TCanvas("cTPCnsigPion","cTPCnsigPion");
  gPad->SetLogz();
  hTPCnsigPion->GetXaxis()->SetRangeUser(0,5);
  hTPCnsigPion->GetYaxis()->SetRangeUser(-10,40);
  hTPCnsigPion->DrawCopy("COLZ");
  // TOF
  TH2F * hTOFnsigPion = new TH2F(*(TH2F*)((TH2F*)hman_data->GetNSigHistogram("hHistNSigPionTOF"))->Clone());
  hTOFnsigPion->SetTitle("TOF NSigma Distribution for Pions;p_{T} (GeV/c);n#sigma");
  TCanvas * cTOFnsigPion = new TCanvas("cTOFnsigPion","cTOFnsigPion");
  gPad->SetLogz();
  hTOFnsigPion->GetYaxis()->SetRangeUser(-10,40);
  hTOFnsigPion->GetXaxis()->SetRangeUser(0,5);
  hTOFnsigPion->DrawCopy("COLZ");
  // TPCTOF
  TH2F * hTPCTOFnsigPion = new TH2F(*(TH2F*)((TH2F*)hman_data->GetNSigHistogram("hHistNSigPionTPCTOF"))->Clone());
  hTPCTOFnsigPion->SetTitle("TPC+TOF NSigma Distribution for Pions;p_{T} (GeV/c);n#sigma");
  TCanvas * cTPCTOFnsigPion = new TCanvas("cTPCTOFnsigPion","cTPCTOFnsigPion");
  gPad->SetLogz();
  hTPCTOFnsigPion->GetYaxis()->SetRangeUser(-5,40);
  hTPCTOFnsigPion->GetXaxis()->SetRangeUser(0,5);
  hTPCTOFnsigPion->DrawCopy("COLZ");


  //------------------------------------------------------
  // plot the MC True Generated Primaries for comparison -
  //------------------------------------------------------
  TH1F * MCTruthSpectra[nCharge*nPart];
  Double_t events_mc = ecuts_mc->NumberOfEvents();
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  TString hname = Form("hHistPtGenTruePrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());
	  MCTruthSpectra[index] = (TH1F*)((TH1F*)hman_mc->GetPtHistogram1D(hname.Data(),1,1)->Clone());
	  MCTruthSpectra[index]->SetTitle("Generated Primaries Spectra (MC Truth);p_{T} (GeV/c);");
	  MCTruthSpectra[index]->SetMarkerStyle(Marker[index]);
	  MCTruthSpectra[index]->SetMarkerColor(Color[ipart]);
	  MCTruthSpectra[index]->SetLineColor(Color[ipart]);
	  MCTruthSpectra[index]->Scale(1./events_mc,"width");
	}
    }
  // neg/pos ratios
  TH1F * NegPosRatio[nPart];
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      NegPosRatio[ipart] = (TH1F*)MCTruthSpectra[ipart+nPart]->Clone();
      NegPosRatio[ipart]->Divide(MCTruthSpectra[ipart]);
      NegPosRatio[ipart]->SetTitle("Neg/Pos;p_{T} (GeV/c);");
    }
  // drawing
  TCanvas * cMCTruthSpectra = new TCanvas("cMCTruthSpectra","cMCTruthSpectra");
  cMCTruthSpectra->Divide(1,2);
  TLegend * lMCTruthSpectra = new TLegend(.69,.69,.89,.99);
  lMCTruthSpectra->SetFillColor(0);
  TLegend * lPosNegRatios = new TLegend (.69,.69,.99,.99);
  lPosNegRatios->SetFillColor(0);
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  cMCTruthSpectra->cd(1);
	  gPad->SetLogy();
	  if (index == 0) MCTruthSpectra[index]->DrawCopy("EP");
	  else MCTruthSpectra[index]->DrawCopy("EPsame");
	  lMCTruthSpectra->AddEntry(MCTruthSpectra[index],Names[index].Data(),"p");

	  cMCTruthSpectra->cd(2);
	  NegPosRatio[ipart]->GetYaxis()->SetRangeUser(0.4,1.1);
	  if (index == 0) NegPosRatio[ipart]->DrawCopy("hist][");
	  else NegPosRatio[ipart]->DrawCopy("hist][same");
	  if (icharge == 0) lPosNegRatios->AddEntry(NegPosRatio[ipart],Particle[ipart].Data(),"l");
	}
    }
  cMCTruthSpectra->cd(1);
  lMCTruthSpectra->DrawClone();
  cMCTruthSpectra->cd(2);
  lPosNegRatios->DrawClone();


  //-----------------------------------
  // Get the (normalized) raw spectra -
  //-----------------------------------
  TH1F * Spectra[nCharge*nPart];
  Double_t events_data =  ecuts_data->NumberOfEvents();
  Double_t events_data_1sig =  ecuts_data_1sig->NumberOfEvents();
  Double_t events_mc =  ecuts_mc->NumberOfEvents();
  Double_t events_mc_1sig = ecuts_mc_1sig->NumberOfEvents();
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  TString hname = Form("hHistPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data());
	  Spectra[index] = (TH1F*)((TH1F*)hman_data->GetPtHistogram1D(hname.Data(),-1,-1))->Clone();
	  Spectra[index]->SetTitle(Form("Spectra_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
          Spectra[index]->SetMarkerStyle(Marker[index]);
          Spectra[index]->SetMarkerColor(Color[ipart]); 
          Spectra[index]->SetLineColor(Color[ipart]);
	  Spectra[index]->Scale(1./events_data,"width"); // normalization
	}
    }
  // Set bins outisde the range to contain zeros
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  for (Int_t ibin=0; ibin<Spectra[index]->GetNbinsX(); ibin++)
	    {
	      if (Spectra[index]->GetBinCenter(ibin) < Range[ipart])
		{
		  Spectra[index]->SetBinContent(ibin, 0);
		  Spectra[index]->SetBinError(ibin, 0);
		}
	    }
	}
    }
  // save the raw spectra for use later
  TH1F * RawSpectra[nCharge*nPart];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  RawSpectra[index] = new TH1F(*(TH1F*)Spectra[index]->Clone());
	  RawSpectra[index]->SetTitle("Raw Spectra");
	}
    }

  // plot the raw spectra from DATA
  // also include allch
  TCanvas * cRawSpectra = new TCanvas("cRawSpectra","cRawSpectra");
  TLegend * lRawSpectra = new TLegend(.69,.69,.99,.99);
  lRawSpectra->SetFillColor(0);
  TH1F * hRawAllCh = (TH1F*)((TH1F*)hman_data->GetPtHistogram1D("hHistPtRec",-1,-1))->Clone();
  hRawAllCh->Scale(1./ecuts_data->NumberOfEvents(),"width");
  hRawAllCh->SetTitle("Raw Spectra, Data;p_{T} (GeV/c);");
  hRawAllCh->SetMarkerColor(kGreen+3);
  hRawAllCh->SetLineColor(kGreen+3);
  hRawAllCh->SetMarkerStyle(34);
  hRawAllCh->DrawCopy();
  lRawSpectra->AddEntry(hRawAllCh,"All Particles","p");
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  RawSpectra[index]->SetTitle("Raw Spectra, Data;p_{T} (GeV/c);");
	  RawSpectra[index]->SetMarkerColor(Color[ipart]);
	  RawSpectra[index]->SetLineColor(Color[ipart]);
	  RawSpectra[index]->SetMarkerStyle(Marker[index]);
	  RawSpectra[index]->DrawCopy("same");
	  lRawSpectra->AddEntry(RawSpectra[index],Names[index].Data(),"p");
	}
    }
  lRawSpectra->DrawClone();

  // also get MC raw for comparison purposes
  TH1F * MCSpectra[nCharge*nPart];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  TString hname = Form("hHistPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data());
	  MCSpectra[index] = (TH1F*)((TH1F*)hman_mc->GetPtHistogram1D(hname.Data(),-1,-1))->Clone();
	  MCSpectra[index]->SetTitle(Form("MCSpectra_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
	  MCSpectra[index]->SetMarkerStyle(Marker[index]);
	  MCSpectra[index]->SetMarkerColor(Color[ipart]);
	  MCSpectra[index]->SetLineColor(Color[ipart]);
	  MCSpectra[index]->Scale(1./events_mc,"width"); // normalization
	}
    }
  // Set bins outisde the range to contain zeros
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  for (Int_t ibin=0; ibin<MCSpectra[index]->GetNbinsX(); ibin++)
	    {
	      if (MCSpectra[index]->GetBinCenter(ibin) < Range[ipart])
		{
		  MCSpectra[index]->SetBinContent(ibin, 0);
		  MCSpectra[index]->SetBinError(ibin, 0);
		}
	    }
	}
    }
  // save the raw spectra for use later
  TH1F * MCRawSpectra[nCharge*nPart];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  MCRawSpectra[index] = new TH1F(*(TH1F*)Spectra[index]->Clone());
	  MCRawSpectra[index]->SetTitle("MC Raw Spectra");
	}
    }


  //----------------
  // MC Correction -
  //----------------
  // calculate the percentage of primaries reconstructed w/ 5 sigma PID
  Printf("\n\n-> Calculating MC Correction Factors");
  TH1F * CorrFact[nPart*nCharge];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  TString hname = Form("hHistPtRecTruePrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());
	  Printf("Getting %s",hname.Data());
	  CorrFact[index] = new TH1F(*(TH1F*)((TH1F*) hman_mc_5sig->GetPtHistogram1D(hname.Data(),-1,-1))->Clone());
	  CorrFact[index]->SetName(Form("CorrFact_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
	  CorrFact[index]->SetTitle("MC Correction Factor");
	  CorrFact[index]->SetMarkerStyle(Marker[index]);
	  CorrFact[index]->SetMarkerColor(Color[ipart]);
	  CorrFact[index]->SetLineColor(Color[ipart]);
	  hname=Form("hHistPtGenTruePrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());
	  Printf("... and divide it by %s",hname.Data());
	  CorrFact[index]->Divide(CorrFact[index],(TH1F*)((TH1F*)hman_mc_5sig->GetPtHistogram1D(hname.Data(),1,1))->Clone(),1,1,"B");//binomial error
	  if (index == 0)
	    {
	      TCanvas * cCorrFactMC = new TCanvas("cCorrFactMC","cCorrFactMC");
	      TLegend * lCorrFactMC = new TLegend(.69,.69,.99,.99);
	      lCorrFactMC->SetFillColor(0);
	      CorrFact[index]->DrawCopy("EP");
	    }
	  else CorrFact[index]->DrawCopy("EPsame");
	  lCorrFactMC->AddEntry(CorrFact[index],Names[index].Data(),"p");
	} // end loop over ipart
    } // end loop over icharge
  lCorrFactMC->DrawClone();
  // Correct spectra with MC correction factor
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  Spectra[index]->Divide(CorrFact[index]);
	}
    }
  // save this version of the spectra for use later
  TH1F * SpectraStep1[nCharge*nPart];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  SpectraStep1[index] = new TH1F(*(TH1F*)Spectra[index]->Clone());
	  SpectraStep1[index]->SetTitle("Spectra After MC Correction");
	}
    }
  // do the same with MCSpectra
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  MCSpectra[index]->Divide(CorrFact[index]);
	}
    }
  // save this version of the spectra for use later
  TH1F * MCSpectraStep1[nCharge*nPart];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  MCSpectraStep1[index] = new TH1F(*(TH1F*)MCSpectra[index]->Clone());
	  MCSpectraStep1[index]->SetTitle("MCSpectra After MC Correction");
	}
    }


  //--------------------------------------
  // Correct spectra for PID nsigma cuts -
  //--------------------------------------
  // for the TPC only, assume the nsigma distribution is purely gaussian, but account for the actual means and sigmas
  // the correction factor is the integral of a normal distribution from -3sqrt(2) to 3sqrt(2)
  TF1 * fGaussTPC[nPart];
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      fGaussTPC[ipart] = new TF1(Form("fGaussTPC%s",Particle[ipart].Data()),"gaus",-25,25);
      fGaussTPC[ipart]->SetParameter(0,1); // set a default height for now - will be normalized later
    }
  // pions
  fGaussTPC[0]->SetParameter(1,0.04);
  fGaussTPC[0]->SetParameter(2,1.02);
  // kaons
  fGaussTPC[1]->SetParameter(1,-0.06);
  fGaussTPC[1]->SetParameter(2,0.97);
  // protons
  fGaussTPC[2]->SetParameter(1,0.25);
  fGaussTPC[2]->SetParameter(2,0.72);
  // normalization
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      fGaussTPC[ipart]->SetParameter(0, 1./fGaussTPC[ipart]->Integral(-20,20));
    }
  // calculate correction factors
  Float_t PIDCorrFactTPC[nPart];
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      PIDCorrFactTPC[ipart] = fGaussTPC[ipart]->Integral(-3*TMath::Sqrt(2), 3*TMath::Sqrt(2));
    }
  cout << "\n\n\n--- TPC ONLY CORRECTION FACTORS: ---\n"
       << "Pions: " << PIDCorrFactTPC[0] << endl
       << "Kaons: " << PIDCorrFactTPC[1] << endl
       << "Protons: " << PIDCorrFactTPC[2] << endl << endl << endl;
  // for the TPC+TOF range, get the correction factor from integrating the TPC+TOF nsigma distribution
  TFile * fTPCTOFin = TFile::Open("TPCTOFnsigma.root");
  TH1F * hTPCTOFnsigma[nPart];
  Float_t PIDCorrFactTPCTOF[nPart];
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      hTPCTOFnsigma[ipart] = (TH1F*)((TH1F*)((TCanvas*)fTPCTOFin->Get("cTPCTOFnsigmaPredicted"))->FindObject(Form("hTPCTOFnsigma%s",Particle[ipart].Data())))->Clone();
      PIDCorrFactTPCTOF[ipart] = hTPCTOFnsigma[ipart]->Integral(hTPCTOFnsigma[ipart]->FindBin(0), hTPCTOFnsigma[ipart]->FindBin(3));
      PIDCorrFactTPCTOF[ipart] /= hTPCTOFnsigma[ipart]->Integral(hTPCTOFnsigma[ipart]->FindBin(0), hTPCTOFnsigma[ipart]->FindBin(20));
    }
  cout << "\n\n\n-- TPC+TOF CORRECTION FACTORS: ---\n"
       << "Pions: " << PIDCorrFactTPCTOF[0] << endl
       << "Kaons: " << PIDCorrFactTPCTOF[1] << endl
       << "Protons: " << PIDCorrFactTPCTOF[2] << endl << endl << endl;

  // calculate seperate correction factors for MC b/c the parameters will be different...
  TF1 * fGaussTPCMC[nPart];
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      fGaussTPCMC[ipart] = new TF1(Form("fGaussTPCMC%s",Particle[ipart].Data()),"gaus",-25,25);
      fGaussTPCMC[ipart]->SetParameter(0,1); // set a default height for now - will be normalized later
    }
  // pions
  fGaussTPCMC[0]->SetParameter(1,-0.03);
  fGaussTPCMC[0]->SetParameter(2,1.35);
  // kaons
  fGaussTPCMC[1]->SetParameter(1,0.14);
  fGaussTPCMC[1]->SetParameter(2,1.27);
  // protons
  fGaussTPCMC[2]->SetParameter(1,-0.17);
  fGaussTPCMC[2]->SetParameter(2,1.03);
  // normalization
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      fGaussTPCMC[ipart]->SetParameter(0, 1./fGaussTPCMC[ipart]->Integral(-20,20));
    }
  // calculate correction factors
  Float_t PIDCorrFactTPCMC[nPart];
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      PIDCorrFactTPCMC[ipart] = fGaussTPCMC[ipart]->Integral(-3*TMath::Sqrt(2), 3*TMath::Sqrt(2));
    }
  cout << "\n\n\n--- TPC ONLY CORRECTION FACTORS (MC): ---\n"
       << "Pions: " << PIDCorrFactTPCMC[0] << endl
       << "Kaons: " << PIDCorrFactTPCMC[1] << endl
       << "Protons: " << PIDCorrFactTPCMC[2] << endl << endl << endl;
  // for the TPC+TOF range, get the correction factor from integrating the TPC+TOF nsigma distribution
  TFile * fTPCTOFin = TFile::Open("TPCTOFnsigma.root");
  TH1F * hTPCTOFnsigmaMC[nPart];
  Float_t PIDCorrFactTPCTOFMC[nPart];
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      hTPCTOFnsigmaMC[ipart] = (TH1F*)((TH1F*)((TCanvas*)fTPCTOFin->Get("cTPCTOFnsigmaPredictedMC"))->FindObject(Form("hTPCTOFnsigmaMC%s",Particle[ipart].Data())))->Clone();
      PIDCorrFactTPCTOFMC[ipart] = hTPCTOFnsigmaMC[ipart]->Integral(hTPCTOFnsigmaMC[ipart]->FindBin(0), hTPCTOFnsigmaMC[ipart]->FindBin(3));
      PIDCorrFactTPCTOFMC[ipart] /= hTPCTOFnsigmaMC[ipart]->Integral(hTPCTOFnsigmaMC[ipart]->FindBin(0), hTPCTOFnsigmaMC[ipart]->FindBin(20));
    }
  cout << "\n\n\n-- TPC+TOF CORRECTION FACTORS (MC): ---\n"
       << "Pions: " << PIDCorrFactTPCTOFMC[0] << endl
       << "Kaons: " << PIDCorrFactTPCTOFMC[1] << endl
       << "Protons: " << PIDCorrFactTPCTOFMC[2] << endl << endl << endl;


  // correct spectra, taking into account the range where the TOF is valid
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  for (Int_t ibin=1; ibin<Spectra[index]->GetNbinsX(); ibin++)
	    {
	      if (ibin < Spectra[index]->FindBin(0.6))
		Spectra[index]->SetBinContent(ibin,Spectra[index]->GetBinContent(ibin) / PIDCorrFactTPC[ipart]);		  
	      else
		Spectra[index]->SetBinContent(ibin,Spectra[index]->GetBinContent(ibin) / PIDCorrFactTPCTOF[ipart]);
	    }  
	} // end for loop over ipart
    } // end for loop over icharge
  // save this version of the spectra for use later
  TH1F * SpectraStep2[nCharge*nPart];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  SpectraStep2[index] = new TH1F(*(TH1F*)Spectra[index]->Clone());
	  SpectraStep2[index]->SetTitle("Spectra After PID Correction");
	}
    }
  // also do this for MCSpectra
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  for (Int_t ibin=1; ibin<MCSpectra[index]->GetNbinsX(); ibin++)
	    {
	      if (ibin < MCSpectra[index]->FindBin(0.6))
		MCSpectra[index]->SetBinContent(ibin,MCSpectra[index]->GetBinContent(ibin) / PIDCorrFactTPCMC[ipart]);		  
	      else
		MCSpectra[index]->SetBinContent(ibin,MCSpectra[index]->GetBinContent(ibin) / PIDCorrFactTPCTOFMC[ipart]);
	    }  
	} // end for loop over ipart
    } // end for loop over icharge
  // save this version of the spectra for use later
  TH1F * MCSpectraStep2[nCharge*nPart];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  MCSpectraStep2[index] = new TH1F(*(TH1F*)MCSpectra[index]->Clone());
	  MCSpectraStep2[index]->SetTitle("MCSpectra After PID Correction");
	}
    }

  
  //-------------------------------------------------------------------------
  // Correct spectra for contamination from other Pions, Kaons, and Protons -
  //-------------------------------------------------------------------------
  // calculate the correction factors first
  TH1F * ContCorrFact[nPart*nCharge];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  TString hname = Form("hHistPtRecTrue%s%s",Particle[ipart].Data(),Sign[icharge].Data());
	  ContCorrFact[index] = new TH1F(*(TH1F*)((TH1F*)hman_mc->GetPtHistogram1D(hname.Data(),-1,-1))->Clone());
	  ContCorrFact[index]->SetTitle("Contamination Correction Factor");
	  ContCorrFact[index]->SetMarkerStyle(Marker[index]);
	  ContCorrFact[index]->SetMarkerColor(Color[ipart]);
	  ContCorrFact[index]->SetLineColor(Color[ipart]);
	  hname = Form("hHistPtRecSigma%s%s",Particle[ipart].Data(),Sign[icharge].Data());
	  ContCorrFact[index]->Divide(ContCorrFact[index],(TH1F*)((TH1F*)hman_mc->GetPtHistogram1D(hname.Data(),-1,-1))->Clone(),-1,-1,"B"); // binomial error
	  if (index == 0)
	    {
	      TCanvas * cContCorrFact = new TCanvas("cContCorrFact","cContCorrFact");
	      TLegend * lContCorrFact = new TLegend(.69,.69,.99,.99);
	      ContCorrFact[index]->DrawCopy("EP");
	    }
	  else ContCorrFact[index]->DrawCopy("EPsame");
	  lContCorrFact->AddEntry(ContCorrFact[index],Names[index].Data(),"p");
	} // end for loop over ipart
    } // end for loop over icharge
  lContCorrFact->DrawClone();
  // then use them to correct the spectra
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for(Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  Spectra[index]->Multiply(ContCorrFact[index]);
	}
    }
  // save this version of the spectra for use later
  TH1F * SpectraStep3[nCharge*nPart];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  SpectraStep3[index] = new TH1F(*(TH1F*)Spectra[index]->Clone());
	  SpectraStep3[index]->SetTitle("Spectra After Contamination Correction");
	}
    }
  // also do this for MCSpectra
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for(Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  MCSpectra[index]->Multiply(ContCorrFact[index]);
	}
    }
  // save this version of the spectra for use later
  TH1F * MCSpectraStep3[nCharge*nPart];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  MCSpectraStep3[index] = new TH1F(*(TH1F*)MCSpectra[index]->Clone());
	  MCSpectraStep3[index]->SetTitle("MCSpectra After Contamination Correction");
	}
    }


  //--------------------------------------------------
  // Correct spectra for the Geant-Fluka discrepancy -
  //--------------------------------------------------
  GFCorrection(Spectra,tcuts_data,fout);
  // save this version of the spectra for use later
  TH1F * SpectraStep4[nCharge*nPart];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  SpectraStep4[index] = new TH1F(*(TH1F*)Spectra[index]->Clone());
	  SpectraStep4[index]->SetTitle("Spectra After Geant-Fluka Correction");
	}
    }
  // NOTE! DO NOT USE GF FOR MCSPECTRA!


  //--------------------------------------------------
  // Correct spectra for the TOF matching efficiency -
  //--------------------------------------------------
  //Matching efficiency in data and Monte Carlo
  TCanvas *cMatchingEff=new TCanvas("cMatchingEff","cMatchingEff",700,500);
  TH1F *hMatcEffPos_data=(TH1F*)tcuts_data->GetHistoNMatchedPos();
  hMatcEffPos_data->Divide((TH1F*)tcuts_data->GetHistoNSelectedPos());
  hMatcEffPos_data->SetTitle("Matching Eff Pos - data");
  TH1F *hMatcEffNeg_data=(TH1F*)tcuts_data->GetHistoNMatchedNeg();
  hMatcEffNeg_data->Divide((TH1F*)tcuts_data->GetHistoNSelectedNeg());
  hMatcEffNeg_data->SetTitle("Matching Eff Neg - data");
  hMatcEffNeg_data->SetLineColor(2);
  /*  
  // for my (John Groh) report, plot the matching efficiency just for data
  TCanvas * cMatchEffDataOnly = new TCanvas("cMatchEffDataOnly","cMatchEffDataOnly");
  hMatcEffPos_data->SetTitle("TOF Matching Efficiency;p_{T} (GeV/c);");
  hMatcEffPos_data->DrawCopy("lhist");
  hMatcEffNeg_data->SetTitle("TOF Matching Efficiency;p_{T} (GeV/c);");
  hMatcEffNeg_data->DrawCopy("lhistsame");
  TLegend * lMatcEffDataOnly = new TLegend(.69,.69,.99,.99);
  lMatcEffDataOnly->SetFillColor(0);
  lMatcEffDataOnly->AddEntry(hMatcEffPos_data,"Positive Particles","l");
  lMatcEffDataOnly->AddEntry(hMatcEffNeg_data,"Negative Particles","l");
  lMatcEffDataOnly->DrawClone();
  */

  TH1F *hMatcEffPos_mc=(TH1F*)tcuts_mc->GetHistoNMatchedPos();
  hMatcEffPos_mc->Divide((TH1F*)tcuts_mc->GetHistoNSelectedPos());
  hMatcEffPos_mc->SetTitle("Matching Eff Pos - mc");
  hMatcEffPos_mc->SetLineStyle(2);
  TH1F *hMatcEffNeg_mc=(TH1F*)tcuts_mc->GetHistoNMatchedNeg();
  hMatcEffNeg_mc->Divide((TH1F*)tcuts_mc->GetHistoNSelectedNeg());
  hMatcEffNeg_mc->SetTitle("Matching Eff Neg - mc");
  hMatcEffNeg_mc->SetLineColor(2);
  hMatcEffNeg_mc->SetLineStyle(2);
  cMatchingEff->Divide(1,2);
  cMatchingEff->cd(1);
  gPad->SetGridy();
  gPad->SetGridx();
  hMatcEffPos_data->DrawClone("lhist");
  hMatcEffNeg_data->DrawClone("lhistsame");
  hMatcEffPos_mc->DrawClone("lhistsame");
  hMatcEffNeg_mc->DrawClone("lhistsame");
  gPad->BuildLegend()->SetFillColor(0);
  hMatcEffPos_data->Divide(hMatcEffPos_mc);
  hMatcEffNeg_data->Divide(hMatcEffNeg_mc);
  cMatchingEff->cd(2);
  gPad->SetGridy();
  gPad->SetGridx();
  hMatcEffPos_data->DrawClone("lhist");
  hMatcEffNeg_data->DrawClone("lhistsame");
  TF1 *pol0MatchPos_data=new TF1("pol0MatchPos_data","pol0",2.5,5);
  hMatcEffPos_data->Fit("pol0MatchPos_data","MNR");
  pol0MatchPos_data->DrawClone("same");
  TF1 *pol0MatchNeg_data=new TF1("pol0MatchNeg_data","pol0",2.5,5);
  hMatcEffNeg_data->Fit("pol0MatchNeg_data","MNR");
  pol0MatchNeg_data->SetLineColor(2);
  pol0MatchNeg_data->DrawClone("same");
  Float_t ScalingMatchingPos=pol0MatchPos_data->GetParameter(0);
  Float_t ScalingMatchingNeg=pol0MatchNeg_data->GetParameter(0);
  //Correction spectra for matching efficiency
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      for (Int_t ibin=1; ibin<=Spectra[ipart]->GetNbinsX(); ibin++)
	{
	  Float_t ptspectra=Spectra[ipart]->GetBinCenter(ibin);
	  if (ptspectra < tcuts_data->GetPtTOFMatching()) continue;
	  Spectra[ipart]->SetBinContent(ibin, (Spectra[ipart]->GetBinContent(ibin) / ScalingMatchingPos));
	  Spectra[ipart+3]->SetBinContent(ibin, (Spectra[ipart+3]->GetBinContent(ibin) / ScalingMatchingNeg));
	}
    }
  // save this version of the spectra for use later
  TH1F * SpectraStep5[nCharge*nPart];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  SpectraStep5[index] = new TH1F(*(TH1F*)Spectra[index]->Clone());
	  SpectraStep5[index]->SetTitle("Spectra After Matching Efficiency Correction");
	}
    }
  // NOTE! DO NOT USE THE MATCHING EFFICIENCY CORRECTION FOR MCSPECTRA!


  //-----------------------------------------------------
  // Secondaries correction - use MC, not DCA from data -
  //------------------------------------------------------
  // (correction factor for secondaries) = (true primaries) / (primaries)
  TH1F * SecondaryCorrFact[nCharge*nPart];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  cout << "index = " << index << endl;
	  TString hname = Form("hHistPtRecTruePrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data());
	  SecondaryCorrFact[index] = (TH1F*)((TH1F*)hman_mc->GetPtHistogram1D(hname.Data(),-1,-1))->Clone();
	  hname = Form("hHistPtRecTrue%s%s",Particle[ipart].Data(),Sign[icharge].Data());
	  SecondaryCorrFact[index]->Divide(SecondaryCorrFact[index], (TH1F*)((TH1F*)hman_mc->GetPtHistogram1D(hname.Data(),-1,-1))->Clone(), 1, 1, "B"); // binomial error!
	  SecondaryCorrFact[index]->SetName(Form("SecondaryCorrFact_%s%s",Particle[ipart].Data(),Sign[icharge].Data()));
	  SecondaryCorrFact[index]->SetTitle("Correction Factor for Secondaries");
	  SecondaryCorrFact[index]->SetStats(kFALSE);
	  SecondaryCorrFact[index]->SetMarkerStyle(Marker[index]);
	  SecondaryCorrFact[index]->SetMarkerColor(Color[ipart]);
	  SecondaryCorrFact[index]->SetLineColor(Color[ipart]);
	  if (index == 0)
	    {
	      TCanvas * cSecCorrFactMC = new TCanvas("cSecCorrFactMC","cSecCorrFactMC");
	      TLegend * lSecCorrFactMC = new TLegend(.69,.69,.99,.99);
	      lSecCorrFactMC->SetFillColor(0);
	      SecondaryCorrFact[index]->DrawCopy("EP");		  
	    }
	  else SecondaryCorrFact[index]->DrawCopy("EPsame");
	  lSecCorrFactMC->AddEntry(SecondaryCorrFact[index],Names[index],"ep");
	} // end for loop over ipart
    } // end for loop over icharge
  lSecCorrFactMC->DrawClone();

  // correct Spectra with this
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  Spectra[index]->Multiply(SecondaryCorrFact[index]);
	}
    }

  // also correct MCSpectra with this
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  MCSpectra[index]->Multiply(SecondaryCorrFact[index]);
	}
    }
  // save this version of the spectra for use later
  TH1F * SpectraStep6[nCharge*nPart];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  SpectraStep6[index] = new TH1F(*(TH1F*)Spectra[index]->Clone());
	  SpectraStep6[index]->SetTitle("Spectra After Secondaries Correction");
	}
    }
  // also save this version of MCSpectra
  TH1F * MCSpectraStep6[nCharge*nPart];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  MCSpectraStep6[index] = new TH1F(*(TH1F*)MCSpectra[index]->Clone());
	  MCSpectraStep6[index]->SetTitle("MCSpectra After Secondaries Correction");
	}
    }
  
  //-------------------------
  // draw the final spectra -
  //-------------------------
  TCanvas * cSpectra = new TCanvas("cSpectra","cSpectra",100,100,700,500);
  gPad->SetGridy();
  gPad->SetLogy();
  TLegend * lSpectra = new TLegend(.69,.69,.99,.99);
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  Spectra[index]->SetTitle("Final Spectra");
	  if (index == 0) Spectra[index]->DrawClone("EP");
	  else Spectra[index]->DrawClone("EPsame");
	  lSpectra->AddEntry(Spectra[index],Names[index].Data(),"p");
	}
    }
  lSpectra->DrawClone();

  //-----------------------------------------------
  // compare MCSpectra with hHistPtGenTruePrimary -
  //-----------------------------------------------
  TCanvas * cCompareMCRecWithMCTrue = new TCanvas("cCompareMCRecWithMCTrue","cCompareMCRecWithMCTrue",250,175,700,500);
  cCompareMCRecWithMCTrue->Divide(3,2);
  TH1F * hMCTruth[nPart*nCharge];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  hMCTruth[index] = new TH1F(*(TH1F*)((TH1F*)hman_mc->GetPtHistogram1D(Form("hHistPtGenTruePrimary%s%s",Particle[ipart].Data(),Sign[icharge].Data()),1,1))->Clone());
	  hMCTruth[index]->Scale(1./ecuts_mc->NumberOfEvents(),"width"); // normalize it first!
	  hMCTruth[index]->SetMarkerStyle(Marker[index]);
	  hMCTruth[index]->SetMarkerColor(Color[ipart]);
	  hMCTruth[index]->SetLineColor(Color[ipart]);
	}
    }
  // set bins outside range to zero
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  for (Int_t ibin=0; ibin<hMCTruth[index]->GetNbinsX(); ibin++)
	    {
	      if (hMCTruth[index]->GetBinCenter(ibin) < Range[ipart])
		{
		  hMCTruth[index]->SetBinContent(ibin, 0);
		  hMCTruth[index]->SetBinError(ibin, 0);
		}
	    }
	}
    }


  TH1F * hRatioMCRecTrue[nPart*nCharge];  
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      for (Int_t icharge=0; icharge<nCharge; icharge++)
	{
	  Int_t index = ipart + 3*icharge;	  
	  hRatioMCRecTrue[index] = (TH1F*)MCSpectra[index]->Clone();
	  hRatioMCRecTrue[index]->Divide(hMCTruth[index]);
	  hRatioMCRecTrue[index]->SetTitle("(Corrected MC Spectra) / (MC Truth Spectra);p_{T} (GeV/c);");
	  if (icharge == 1) hRatioMCRecTrue[index]->SetLineStyle(2);
	}
    }
  for (Int_t ipart=0; ipart<nPart; ipart++)
    {
      cCompareMCRecWithMCTrue->cd(ipart+1);
      gPad->SetLogy();
      MCSpectra[ipart]->SetMarkerColor(kGreen+3);
      MCSpectra[ipart+nPart]->SetMarkerColor(kGreen+2);
      MCSpectra[ipart]->SetMarkerSize(0.75);
      MCSpectra[ipart+nPart]->SetMarkerSize(0.75);
      MCSpectra[ipart]->SetTitle(Form("MC Corrected and Truth Spectra, %ss;p_{T} (GeV/c);",Particle[ipart].Data()));
      MCSpectra[ipart+nPart]->SetTitle(Form("MC Corrected and Truth Spectra, %ss;p_{T} (GeV/c);",Particle[ipart].Data()));
      hMCTruth[ipart]->SetTitle(Form("MC Corrected and Truth Spectra, %ss;p_{T} (GeV/c);",Particle[ipart].Data()));
      hMCTruth[ipart+nPart]->SetTitle(Form("MC Corrected and Truth Spectra, %ss;p_{T} (GeV/c);",Particle[ipart].Data()));
      hMCTruth[ipart]->SetMarkerSize(0.75);
      hMCTruth[ipart]->SetMarkerSize(0.75);
      hMCTruth[ipart]->DrawCopy("EP");
      hMCTruth[ipart+nPart]->DrawCopy("EPsame");
      MCSpectra[ipart]->DrawCopy("EPsame");
      MCSpectra[ipart+nPart]->DrawCopy("EPsame");
      TLegend * lCompare = new TLegend(.69,.59,.99,.91);
      lCompare->SetFillColor(0);
      lCompare->AddEntry(hMCTruth[ipart],Form("%s, Truth",Names[ipart].Data()),"p");
      lCompare->AddEntry(hMCTruth[ipart+nPart],Form("%s, Truth",Names[ipart+nPart].Data()),"p");
      lCompare->AddEntry(MCSpectra[ipart],Form("%s, Corrected",Names[ipart].Data()),"p");
      lCompare->AddEntry(MCSpectra[ipart+nPart],Form("%s, Corrected",Names[ipart+nPart].Data()),"p");
      lCompare->DrawClone();

      cCompareMCRecWithMCTrue->cd(ipart+1+nPart);
      hRatioMCRecTrue[ipart]->GetYaxis()->SetRangeUser(.97,1.03);
      hRatioMCRecTrue[ipart]->DrawCopy("hist][");
      hRatioMCRecTrue[ipart+nPart]->DrawCopy("hist][same");
      TLegend * lCompareRatio = new TLegend(.79,.79,.99,.99);
      lCompareRatio->SetFillColor(0);
      lCompareRatio->AddEntry(hRatioMCRecTrue[ipart],Names[ipart].Data(),"l");
      lCompareRatio->AddEntry(hRatioMCRecTrue[ipart+nPart],Names[ipart+nPart].Data(),"l");
      lCompareRatio->DrawClone();
    }


  //---------------------------------------------------------
  // plot the evolution of the spectra afte each correction -
  //---------------------------------------------------------
  TCanvas * cEvolutionOfSpectra = new TCanvas("cEvolutionOfSpectra","cEvolutionOfSpectra",300,200,700,500);
  cEvolutionOfSpectra->Divide(4,2);
  TLegend * lEvolutionOfSpectra = new TLegend(.69,.69,.99,.99);
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  lEvolutionOfSpectra->AddEntry(Spectra[index],Names[index].Data(),"p");
	}
    }
  cEvolutionOfSpectra->cd(1);
  gPad->SetLogy();
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  RawSpectra[index]->SetMarkerSize(.5);
	  if (index == 0) RawSpectra[index]->DrawClone("EP");
	  else RawSpectra[index]->DrawClone("EPsame");
	}
    }
  lEvolutionOfSpectra->DrawClone();
  cEvolutionOfSpectra->cd(2);
  gPad->SetLogy();
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  SpectraStep1[index]->SetMarkerSize(.5);
	  if (index == 0) SpectraStep1[index]->DrawClone("EP");
	  else SpectraStep1[index]->DrawClone("EPsame");
	}
    }
  lEvolutionOfSpectra->DrawClone();
  cEvolutionOfSpectra->cd(3);
  gPad->SetLogy();
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
      	{
	  Int_t index = ipart + 3*icharge;
	  SpectraStep2[index]->SetMarkerSize(.5);
	  if (index == 0) SpectraStep2[index]->DrawClone("EP");
	  else SpectraStep2[index]->DrawClone("EPsame");
	}
    }
  lEvolutionOfSpectra->DrawClone();
  cEvolutionOfSpectra->cd(4);
  gPad->SetLogy();
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  SpectraStep3[index]->SetMarkerSize(.5);
	  if (index == 0) SpectraStep3[index]->DrawClone("EP");
	  else SpectraStep3[index]->DrawClone("EPsame");
	}
    }
  lEvolutionOfSpectra->DrawClone();
  cEvolutionOfSpectra->cd(5);
  gPad->SetLogy();
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  SpectraStep4[index]->SetMarkerSize(.5);
	  if (index == 0) SpectraStep4[index]->DrawClone("EP");
	  else SpectraStep4[index]->DrawClone("EPsame");
	}
    }
  lEvolutionOfSpectra->DrawClone();
  cEvolutionOfSpectra->cd(6);
  gPad->SetLogy();
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  SpectraStep5[index]->SetMarkerSize(.5);
	  if (index == 0) SpectraStep5[index]->DrawClone("EP");
	  else SpectraStep5[index]->DrawClone("EPsame");
	}
    }
  lEvolutionOfSpectra->DrawClone();
  cEvolutionOfSpectra->cd(7);
  gPad->SetLogy();
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  SpectraStep6[index]->SetMarkerSize(.5);
	  if (index == 0) SpectraStep6[index]->DrawClone("EP");
	  else SpectraStep6[index]->DrawClone("EPsame");
	}
    }
  lEvolutionOfSpectra->DrawClone();
  cEvolutionOfSpectra->cd(8);
  gPad->SetLogy();
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  Spectra[index]->SetMarkerSize(.5);
	  if (index == 0) Spectra[index]->DrawClone("EP");
	  else Spectra[index]->DrawClone("EPsame");
	}
    }
  lEvolutionOfSpectra->DrawClone();
  

  //-----------------------------------------
  // Compare spectra with published spectra -
  //-----------------------------------------
  // also add another canvas that just plots the ratios b/c I can't use the spectra in my presentation (hasn't been published yet)
  TCanvas * cRatioCombined = new TCanvas("cRatioCombined","cRatioCombined");
  cRatioCombined->Divide(3,1);
  TCanvas *CratioComb=new TCanvas("CratioComb","CratioComb",150,125,700,500);
  CratioComb->Divide(3,2);
  TString nameComb[6] = {Form("cent%d_pion_plus",ibinToCompare),Form("cent%d_kaon_plus",ibinToCompare),Form("cent%d_proton_plus",ibinToCompare),
			 Form("cent%d_pion_minus",ibinToCompare),Form("cent%d_kaon_minus",ibinToCompare),Form("cent%d_proton_minus",ibinToCompare)};
  TFile *fComb=new TFile("Combined05/SPECTRA_COMB_20120730.root");
  TLegend * lComp[nPart];
  TLegend * lRatio[nPart];
  for(Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for(Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  TH1F * htmp = (TH1F*)((TH1F*)Spectra[index])->Clone("");
	  htmp->SetMarkerSize(1);
	  htmp->GetXaxis()->SetRangeUser(0,5);
	  TH1F * hcomb = (TH1F*)fComb->Get(nameComb[index].Data())->Clone();
	  CratioComb->cd(ipart+1);
	  gPad->SetGridy();
	  gPad->SetGridx();
	  htmp->SetTitle("Comparison with Combined Analysis;p_{T} (GeV/c);");
	  if(icharge == 0) htmp->DrawCopy();
	  else htmp->DrawCopy("same");
	  hcomb->SetMarkerStyle(Marker[index]);
	  hcomb->DrawCopy("same");
	  if (icharge == 0) {lComp[ipart] = new TLegend(.49,.49,.89,.89); lComp[ipart]->SetFillColor(0);}
	  lComp[ipart]->AddEntry(htmp,Form("%s, This Analysis",Names[index].Data()),"p");
	  lComp[ipart]->AddEntry(hcomb,Form("%s, Combined Analysis",Names[index].Data()),"lep");
	  if (icharge == 1) lComp[ipart]->DrawClone();
 
	  CratioComb->cd(ipart+nPart+1);
	  gPad->SetGridy();
	  gPad->SetGridx();
	  htmp->Divide(hcomb);
	  htmp->SetTitle("(This Analysis) / (Combined Analysis);p_{T} (GeV/c);");
	  htmp->GetYaxis()->SetRangeUser(0.8, 1.2);
	  if (icharge == 0)
	    {
	      CratioComb->cd(ipart+nPart+1);
	      htmp->DrawClone("hist][");
	      cRatioCombined->cd(ipart+1);
	      htmp->DrawClone("hist][");
	      lRatio[ipart] = new TLegend(.69,.69,.89,.89);
	      lRatio[ipart]->SetFillColor(0);
	    }
	  else 
	    {
	      htmp->SetLineStyle(2);
	      CratioComb->cd(ipart+nPart+1);
	      htmp->DrawClone("histsame][");
	      cRatioCombined->cd(ipart+1);
	      htmp->DrawClone("histsame][");
	    }
	  lRatio[ipart]->AddEntry(htmp,Names[index].Data(),"l");
	  if (icharge == 1)
	    {
	      CratioComb->cd(ipart+nPart+1);
	      lRatio[ipart]->DrawClone();
	      cRatioCombined->cd(ipart+1);
	      lRatio[ipart]->DrawClone();    
	    }
	}
    }
  //comparison with charged hadron
  Printf("\n\n-> ChargedHadron comparison");
  TH1F *hChHad_data=(TH1F*)((TH1F*)hman_data->GetPtHistogram1D("hHistPtRec",-1,-1))->Clone();
  TH1F* hMatchCorrectionAllCh=(TH1F*)hMatcEffPos_data->Clone("hMatchCorrectionAllCh");
  hMatchCorrectionAllCh->Add(hMatcEffNeg_data); //correction for Matching efficiency!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  hMatchCorrectionAllCh->Scale(0.5);
  for(Int_t ibin=1;ibin<hChHad_data->GetNbinsX();ibin++){
    Float_t ptch=hChHad_data->GetBinCenter(ibin);
    if(ptch<tcuts_data->GetPtTOFMatching())continue;
    hChHad_data->SetBinContent(ibin,2*(hChHad_data->GetBinContent(ibin)/(ScalingMatchingPos+ScalingMatchingNeg)));
  }
  //fraction of sec in MC
  TH1F *hPrimRec_mc=(TH1F*)((TH1F*)hman_mc->GetPtHistogram1D("hHistPtRecPrimary",-1,-1))->Clone();
  TH1F *hAllRec_mc=(TH1F*)((TH1F*)hman_mc->GetPtHistogram1D("hHistPtRec",-1,-1))->Clone();
  for(Int_t ibin=1;ibin<=hChHad_data->GetNbinsX();ibin++){
    Double_t en_data=hChHad_data->GetBinContent(ibin);
    Double_t en_mc=hAllRec_mc->GetBinContent(ibin);
    Double_t prim_mc=hPrimRec_mc->GetBinContent(ibin);
    if(en_mc!=0)hChHad_data->SetBinContent(ibin,en_data-(en_data*(en_mc-prim_mc)*1.2/en_mc));
    //Printf("Before: %.0f After: %.0f  fraction: %.1f",en_data,hChHad_data->GetBinContent(ibin),hChHad_data->GetBinContent(ibin)/en_data);
  }
  hPrimRec_mc->Divide(hAllRec_mc);
  //efficiency for primaries
  TH1F *hEff_mc=(TH1F*)((TH1F*)hman_mc->GetPtHistogram1D("hHistPtRecPrimary",-1,-1))->Clone();
  hEff_mc->Divide((TH1F*)((TH1F*)hman_mc->GetPtHistogram1D("hHistPtGen",1,1))->Clone());
  TCanvas *cAllChFactors=new TCanvas("cAllChFactors","cAllChFactors",700,500);
  cAllChFactors->Divide(1,2);
  cAllChFactors->cd(1);
  gPad->SetGridy();
  gPad->SetGridx();
  hPrimRec_mc->SetTitle("Prim/All, charged hadron pure MC");
  hPrimRec_mc->DrawClone("lhist");
  gPad->BuildLegend()->SetFillColor(0);
  cAllChFactors->cd(2);
  gPad->SetGridy();
  gPad->SetGridx();
  hEff_mc->SetTitle("Efficiency for Primaries, charged hadron pure MC");
  hEff_mc->DrawClone("lhist");
  gPad->BuildLegend()->SetFillColor(0);
  //Printf("--------%f ",((TH1F*)hman_mc->GetPtHistogram1D("hHistPtGen",1,1))->GetEntries()/1.6/ecuts_mc->NumberOfEvents());
  hChHad_data->Scale(1./events_data,"width");//NORMALIZATION
  hChHad_data->Divide(hEff_mc);//Efficiency
  hChHad_data->Scale(1./(TMath::Abs(tcuts_data->GetEtaMin())+TMath::Abs(tcuts_data->GetEtaMax())));
  hChHad_data->SetTitle("All Ch from AOD");
  hChHad_data->SetName("AllCh");
  TCanvas *cAllCh=new TCanvas("cAllCh","cAllCh",700,500);
  cAllCh->Divide(1,2);
  cAllCh->cd(1);
  gPad->SetGridy();
  gPad->SetGridx();
  hChHad_data->GetYaxis()->SetRangeUser(0,2500);
  hChHad_data->DrawClone();
  TFile *fCh=new TFile("ChargedHadron/SPECTRA_UNID_110906.root");
  TH1F *hCh=fCh->Get("hpt_c0_5");
  hCh->SetTitle("All Ch from Jacek");
  hCh->SetMarkerColor(2);
  hCh->SetLineColor(2);
  //invariant yield
  for(Int_t ibin=0;ibin<hCh->GetNbinsX();ibin++){
    hCh->SetBinContent(ibin,hCh->GetBinContent(ibin)*(2*TMath::Pi()*hCh->GetBinCenter(ibin)));
    hCh->SetBinError(ibin,hCh->GetBinError(ibin)*(2*TMath::Pi()*hCh->GetBinCenter(ibin)));
  }
  hCh->DrawClone("same");
  gPad->BuildLegend()->SetFillColor(0);
  TH1F *gRatio=AliPWGHistoTools::DivideHistosDifferentBins(hChHad_data,hCh);
  gRatio->SetMaximum(1.3);
  gRatio->SetMinimum(.7);
  cAllCh->cd(2);
  gPad->SetGridy();
  gPad->SetGridx();
  gRatio->Print("all");
  gRatio->SetBinContent(48,1.02);
  gRatio->SetBinContent(49,1.022);
  gRatio->SetBinContent(50,1.021);
  gRatio->SetBinContent(51,1.026);
  gRatio->SetBinContent(52,1.027);
  gRatio->SetBinError(48,0.056);
  gRatio->SetBinError(49,0.056);
  gRatio->SetBinError(50,0.056);
  gRatio->SetBinError(51,0.056);
  gRatio->SetBinError(52,0.056);
  gRatio->GetYaxis()->SetRangeUser(-1.5,1.5);
  gRatio->DrawClone("");
  
  //fitting
  TCanvas *cFitChargHad=new TCanvas("cFitChargHad","cFitChargHad",700,500);
  gPad->SetGridy();
  gPad->SetGridx();
  hChHad_data->DrawClone();
  //Fitting the sum of all particles
  AliPWGFunc * fm = new AliPWGFunc;
  fm->SetVarType(AliPWGFunc::kdNdpt);
  Float_t fitmin = 0.3;
  Float_t fitmax = 3;
  TF1 * func = 0;
  Int_t normPar = 0;
  func = fm->GetBGBW(0.13,0.6,0.3, 1, 1e5);
  func->SetParLimits(1, 0.1, 0.99);
  func->SetParLimits(2, 0.01, 1);
  func->SetParLimits(3, 0.01, 2);
  TH1F * hToFit = hChHad_data;
  hToFit->Fit(func,"N","VMRN",fitmin,fitmax);
  // // if(!AliPWGHistoTools::Fit(hToFit,func,fitmin,fitmax)) {
  // //   cout << " FIT ERROR " << endl;
  // //   return;      
  // // }
  Double_t yieldTools, yieldETools;
  Double_t partialYields[3],partialYieldsErrors[3]; 
  AliPWGHistoTools::GetYield(hToFit, func, yieldTools, yieldETools,1, 5, partialYields,partialYieldsErrors);
  func->DrawClone("same");   
  Printf("TOTAL YIELD (AOD Charged Hadron) : %f +- %f",yieldTools,yieldETools);
  //Fit All Charged
  hToFit = hCh;
  hToFit->Fit(func,"N","VMRN",fitmin,fitmax);
  AliPWGHistoTools::GetYield(hToFit, func, yieldTools, yieldETools,1, 5, partialYields,partialYieldsErrors);
  func->SetLineColor(2);
  hCh->DrawClone("same");
  func->DrawClone("same");   
  gPad->BuildLegend()->SetFillColor(0);
  Printf("TOTAL YIELD (JACEK): %f +- %f",yieldTools,yieldETools);
  //sumID vs AllCh
  //Convert spectra to dNdeta and sum
  TH1F * hsum = (TH1F*) Spectra[0]->Clone();
  hsum->Reset("all");
  hsum->SetMarkerSize(1);
  Double_t epsilon = 0.0001;
  for(Int_t icharge = 0; icharge < 2; icharge++){
    for(Int_t ipart = 0; ipart < 3; ipart++){
      Int_t index=ipart+3*icharge;
      TH1F *htmp =(TH1F*)Spectra[index]->Clone("htmp");
      Int_t nbin = htmp->GetNbinsX();
      for(Int_t ibin = 1; ibin <= nbin; ibin++){
	Double_t pt = htmp->GetBinCenter(ibin);
	Double_t eta=0.8;//////////////////eta range///////////////////////////////////////
	Double_t jacobian =eta2y(pt,mass[ipart],eta)/eta;
	//Printf("jacobian: %f pt:%f   BinContent:%f  mass:%f",jacobian,pt,htmp->GetBinContent(ibin),mass[ipart]);
	htmp->SetBinContent(ibin,htmp->GetBinContent(ibin)*jacobian);
	htmp->SetBinError(ibin,htmp->GetBinError(ibin)*jacobian);
	Int_t ibinSum = hsum->FindBin(pt);
	if ( htmp->GetBinContent(ibin) > 0 && 
	     (TMath::Abs(htmp->GetBinLowEdge(ibin)   - hsum->GetBinLowEdge(ibinSum)) > epsilon || 
	      TMath::Abs(htmp->GetBinLowEdge(ibin+1) - hsum->GetBinLowEdge(ibinSum+1)) )
	     ) {
	  cout << "DISCREPANCY IN BIN RANGES" << endl;
	  cout << pt << " " << ibinSum << " " << ibin  << "; " << h->GetBinContent(ibin) << endl
	       << h->GetBinLowEdge(ibin) << "-"  << h->GetBinLowEdge(ibin+1) << endl
	       << hsum->GetBinLowEdge(ibinSum) << "-"  << hsum->GetBinLowEdge(ibinSum+1) << endl;
	  cout << "" << endl;	    
	}
      }
      hsum->Add(htmp);
    }
  }
  hsum->SetTitle("Sum ID from AOD");
  TCanvas *cChargHadComp=new TCanvas("cChargHadComp","cChargHadComp",200,150,700,500);
  cChargHadComp->Divide(1,2);
  cChargHadComp->cd(1);
  gPad->SetGridy();
  gPad->SetGridx();
  hsum->SetTitle("Charged Hadron Comparison;p_{T} (GeV/c);");
  hsum->DrawClone();
  hToFit = hsum;
  hToFit->Fit(func,"N","VMRN",fitmin,fitmax);
  AliPWGHistoTools::GetYield(hToFit, func, yieldTools, yieldETools,1, 5, partialYields,partialYieldsErrors);
  func->SetLineColor(2);
  Printf("TOTAL YIELD (Pi+K+p): %f +- %f",yieldTools,yieldETools);
  hChHad_data->SetMarkerColor(2);
  hChHad_data->SetLineColor(2);
  hChHad_data->DrawClone("same");
  TLegend * lChargHadComp = new TLegend(.69,.69,.89,.89);
  lChargHadComp->AddEntry(hsum,"Sum ID from AOD","p");
  lChargHadComp->AddEntry(hChHad_data,"All Ch from AOD","p");
  lChargHadComp->DrawClone();
  cChargHadComp->cd(2);
  gPad->SetGridy();
  gPad->SetGridx();
  TH1F *hRatio=AliPWGHistoTools::DivideHistosDifferentBins(hsum,hChHad_data);
  hRatio->SetTitle("(Sum ID from AOD) / (All Ch from AOD);p_{T} (GeV/c);");
  hRatio->SetMaximum(1.2);
  hRatio->SetMinimum(.8);
  hRatio->DrawClone("");


  //----------------------------------------------------------
  // Compare each correction step with the published spectra -
  //----------------------------------------------------------
  TCanvas * cCompareToPublishedAtEachStep = new TCanvas("cCompareToPublishedAtEachStep","cCompareToPublishedAtEachStep",350,225,700,500);
  cCompareToPublishedAtEachStep->Divide(4,2);
  TLegend * lCompareToPublishedAtEachStep = new TLegend(.69,.69,.99,.99);
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  lCompareToPublishedAtEachStep->AddEntry(Spectra[index],Names[index].Data(),"p");
	}
    }
  TH1F * hComb[nPart*nCharge];
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  fComb->cd();
	  hComb[index] = (TH1F*)fComb->Get(nameComb[index].Data())->Clone();
	}
    }
  cCompareToPublishedAtEachStep->cd(1);
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  RawSpectra[index]->Divide(hComb[index]);
	  if (icharge == 1) RawSpectra[index]->SetLineStyle(2);
	  if (index == 0)
	    {
	      RawSpectra[index]->GetYaxis()->SetRangeUser(0,1);
	      RawSpectra[index]->DrawCopy("hist][");
	    }
	  else RawSpectra[index]->DrawCopy("hist][same");
	}
    }
  lCompareToPublishedAtEachStep->DrawClone();
  cCompareToPublishedAtEachStep->cd(2);
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  SpectraStep1[index]->Divide(hComb[index]);
	  if (icharge == 1) SpectraStep1[index]->SetLineStyle(2);
	  if (index == 0)
	    {
	      SpectraStep1[index]->GetYaxis()->SetRangeUser(.5,1.5);
	      SpectraStep1[index]->DrawCopy("hist][");
	    }
	  else SpectraStep1[index]->DrawCopy("hist][same");
	}
    }
  lCompareToPublishedAtEachStep->DrawClone();
  cCompareToPublishedAtEachStep->cd(3);
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  SpectraStep2[index]->Divide(hComb[index]);
	  if (icharge == 1) SpectraStep2[index]->SetLineStyle(2);
	  if (index == 0)
	    {
	      SpectraStep2[index]->GetYaxis()->SetRangeUser(.5,1.5);
	      SpectraStep2[index]->DrawCopy("hist][");
	    }
	  else SpectraStep2[index]->DrawCopy("hist][same");
	}
    }
  lCompareToPublishedAtEachStep->DrawClone();
  cCompareToPublishedAtEachStep->cd(4);
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  SpectraStep3[index]->Divide(hComb[index]);
	  if (icharge == 1) SpectraStep3[index]->SetLineStyle(2);
	  if (index == 0)
	    {
	      SpectraStep3[index]->GetYaxis()->SetRangeUser(.5,1.5);
	      SpectraStep3[index]->DrawCopy("hist][");
	    }
	  else SpectraStep3[index]->DrawCopy("hist][same");
	}
    }
  lCompareToPublishedAtEachStep->DrawClone();
  cCompareToPublishedAtEachStep->cd(5);
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  SpectraStep4[index]->Divide(hComb[index]);
	  if (icharge == 1) SpectraStep4[index]->SetLineStyle(2);
	  if (index == 0)
	    {
	      SpectraStep4[index]->GetYaxis()->SetRangeUser(.5,1.5);
	      SpectraStep4[index]->DrawCopy("hist][");
	    }
	  else SpectraStep4[index]->DrawCopy("hist][same");
	}
    }
  lCompareToPublishedAtEachStep->DrawClone();
  cCompareToPublishedAtEachStep->cd(6);
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  SpectraStep5[index]->Divide(hComb[index]);
	  if (icharge == 1) SpectraStep5[index]->SetLineStyle(2);
	  if (index == 0)
	    {
	      SpectraStep5[index]->GetYaxis()->SetRangeUser(.5,1.5);
	      SpectraStep5[index]->DrawCopy("hist][");
	    }
	  else SpectraStep5[index]->DrawCopy("hist][same");
	}
    }
  lCompareToPublishedAtEachStep->DrawClone();
  cCompareToPublishedAtEachStep->cd(7);
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  SpectraStep6[index]->Divide(hComb[index]);
	  if (icharge == 1) SpectraStep6[index]->SetLineStyle(2);
	  if (index == 0)
	    {
	      SpectraStep6[index]->GetYaxis()->SetRangeUser(.5,1.5);
	      SpectraStep6[index]->DrawCopy("hist][");
	    }
	  else SpectraStep6[index]->DrawCopy("hist][same");
	}
    }
  lCompareToPublishedAtEachStep->DrawClone();
  cCompareToPublishedAtEachStep->cd(8);
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  Spectra[index]->Divide(hComb[index]);
	  if (icharge == 1) Spectra[index]->SetLineStyle(2);
	  if (index == 0)
	    {
	      Spectra[index]->GetYaxis()->SetRangeUser(.5,1.5);
	      Spectra[index]->DrawCopy("hist][");
	    }
	  else Spectra[index]->DrawCopy("hist][same");
	}
    }
  lCompareToPublishedAtEachStep->DrawClone();



  //---------------------------------------------------------
  // For MC, compare each correction step with the MC truth -
  //---------------------------------------------------------
  TCanvas * cCompareMCtoMCTruthAtEachStep = new TCanvas("cCompareMCtoMCTruthAtEachStep","cCompareMCtoMCTruthAtEachStep",350,225,700,500);
  cCompareMCtoMCTruthAtEachStep->Divide(4,2);
  TLegend * lCompareMCtoMCTruthAtEachStep = new TLegend(.69,.69,.99,.99);
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  lCompareMCtoMCTruthAtEachStep->AddEntry(MCSpectra[index],Names[index].Data(),"p");
	}
    }
  cCompareMCtoMCTruthAtEachStep->cd(1);
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  MCRawSpectra[index]->Divide(hMCTruth[index]);
	  if (icharge == 1) RawSpectra[index]->SetLineStyle(2);
	  if (index == 0)
	    {
	      MCRawSpectra[index]->GetYaxis()->SetRangeUser(0,1);
	      MCRawSpectra[index]->DrawCopy("hist][");
	    }
	  else MCRawSpectra[index]->DrawCopy("hist][same");
	}
    }
  lCompareMCtoMCTruthAtEachStep->DrawClone();
  cCompareMCtoMCTruthAtEachStep->cd(2);
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  MCSpectraStep1[index]->Divide(hMCTruth[index]);
	  if (icharge == 1) MCSpectraStep1[index]->SetLineStyle(2);
	  if (index == 0)
	    {
	      MCSpectraStep1[index]->GetYaxis()->SetRangeUser(.5,1.5);
	      MCSpectraStep1[index]->DrawCopy("hist][");
	    }
	  else MCSpectraStep1[index]->DrawCopy("hist][same");
	}
    }
  lCompareMCtoMCTruthAtEachStep->DrawClone();
  cCompareMCtoMCTruthAtEachStep->cd(3);
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  MCSpectraStep2[index]->Divide(hMCTruth[index]);
	  if (icharge == 1) MCSpectraStep2[index]->SetLineStyle(2);
	  if (index == 0)
	    {
	      MCSpectraStep2[index]->GetYaxis()->SetRangeUser(.5,1.5);
	      MCSpectraStep2[index]->DrawCopy("hist][");
	    }
	  else MCSpectraStep2[index]->DrawCopy("hist][same");
	}
    }
  lCompareMCtoMCTruthAtEachStep->DrawClone();
  cCompareMCtoMCTruthAtEachStep->cd(4);
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  MCSpectraStep3[index]->Divide(hMCTruth[index]);
	  if (icharge == 1) MCSpectraStep3[index]->SetLineStyle(2);
	  if (index == 0)
	    {
	      MCSpectraStep3[index]->GetYaxis()->SetRangeUser(.5,1.5);
	      MCSpectraStep3[index]->DrawCopy("hist][");
	    }
	  else MCSpectraStep3[index]->DrawCopy("hist][same");
	}
    }
  lCompareMCtoMCTruthAtEachStep->DrawClone();
  /*  cCompareMCtoMCTruthAtEachStep->cd(5);
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  MCSpectraStep4[index]->Divide(hMCTruth[index]);
	  if (icharge == 1) MCSpectraStep4[index]->SetLineStyle(2);
	  if (index == 0)
	    {
	      MCSpectraStep4[index]->GetYaxis()->SetRangeUser(.5,1.5);
	      MCSpectraStep4[index]->DrawCopy("hist][");
	    }
	  else MCSpectraStep4[index]->DrawCopy("hist][same");
	}
    }
    lCompareMCtoMCTruthAtEachStep->DrawClone();
  cCompareMCtoMCTruthAtEachStep->cd(6);
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  MCSpectraStep5[index]->Divide(hMCTruth[index]);
	  if (icharge == 1) MCSpectraStep5[index]->SetLineStyle(2);
	  if (index == 0)
	    {
	      MCSpectraStep5[index]->GetYaxis()->SetRangeUser(.5,1.5);
	      MCSpectraStep5[index]->DrawCopy("hist][");
	    }
	  else MCSpectraStep5[index]->DrawCopy("hist][same");
	}
    }
    lCompareMCtoMCTruthAtEachStep->DrawClone();*/
  cCompareMCtoMCTruthAtEachStep->cd(7);
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  MCSpectraStep6[index]->Divide(hMCTruth[index]);
	  if (icharge == 1) MCSpectraStep6[index]->SetLineStyle(2);
	  if (index == 0)
	    {
	      MCSpectraStep6[index]->GetYaxis()->SetRangeUser(.5,1.5);
	      MCSpectraStep6[index]->DrawCopy("hist][");
	    }
	  else MCSpectraStep6[index]->DrawCopy("hist][same");
	}
    }
  lCompareMCtoMCTruthAtEachStep->DrawClone();
  cCompareMCtoMCTruthAtEachStep->cd(8);
  for (Int_t icharge=0; icharge<nCharge; icharge++)
    {
      for (Int_t ipart=0; ipart<nPart; ipart++)
	{
	  Int_t index = ipart + 3*icharge;
	  MCSpectra[index]->Divide(hMCTruth[index]);
	  if (icharge == 1) MCSpectra[index]->SetLineStyle(2);
	  if (index == 0)
	    {
	      MCSpectra[index]->GetYaxis()->SetRangeUser(.5,1.5);
	      MCSpectra[index]->DrawCopy("hist][");
	    }
	  else MCSpectra[index]->DrawCopy("hist][same");
	}
    }
  lCompareMCtoMCTruthAtEachStep->DrawClone();



  // close unwanted canvases
  //  cGFCorrection->Close();
  //cMatchingEff->Close();
  //cContCorrFact->Close();
  //cCorrFactMC->Close();
  //cSecCorrFactMC->Close();
  cAllCh->Close();
  cFitChargHad->Close();
  cAllChFactors->Close();

  // save wanted canvases
  fout->cd();
  cSpectra->Write();
  CratioComb->Write();
  cMatchingEff->Write();
  cEvolutionOfSpectra->Write();  
  cCompareToPublishedAtEachStep->Write();
}



///////////////////////
Double_t eta2y(Double_t pt, Double_t mass, Double_t eta){
  Double_t mt = TMath::Sqrt(mass * mass + pt * pt);
  return TMath::ASinH(pt / mt * TMath::SinH(eta));
}

///////////////////////
void GFCorrection(TH1F **Spectra,AliSpectraAODTrackCuts* tcuts_data, TFile * fout){
  //Geant/Fluka Correction
  Printf("\nGF correction for Kaons");
  //Getting GF For Kaons in TPC
  TGraph *gGFCorrectionKaonPlus=new TGraph();
  gGFCorrectionKaonPlus->SetName("gGFCorrectionKaonPlus");
  gGFCorrectionKaonPlus->SetTitle("gGFCorrectionKaonPlus");
  TGraph *gGFCorrectionKaonMinus=new TGraph();
  gGFCorrectionKaonMinus->SetName("gGFCorrectionKaonMinus");
  gGFCorrectionKaonMinus->SetTitle("gGFCorrectionKaonMinus");
  TString fnameGeanFlukaK="GFCorrection/correctionForCrossSection.321.root";
  TFile *fGeanFlukaK= new TFile(fnameGeanFlukaK.Data());
  TH1F *hGeantFlukaKPos=(TH1F*)fGeanFlukaK->Get("gHistCorrectionForCrossSectionParticles");
  TH1F *hGeantFlukaKNeg=(TH1F*)fGeanFlukaK->Get("gHistCorrectionForCrossSectionAntiParticles");
  //getting GF func for Kaons with TOF
  TF1 *fGFKPosTracking;
  fGFKPosTracking = TrackingEff_geantflukaCorrection(3,kPositive);
  TF1 *fGFKNegTracking;
  fGFKNegTracking = TrackingEff_geantflukaCorrection(3,kNegative);
  TF1 *fGFKPosMatching;
  fGFKPosMatching = TOFmatchMC_geantflukaCorrection(3,kPositive);
  TF1 *fGFKNegMatching;
  fGFKNegMatching = TOFmatchMC_geantflukaCorrection(3,kNegative);
  for(Int_t binK=0;binK<=Spectra[1]->GetNbinsX();binK++){
    if(Spectra[1]->GetBinCenter(binK)<tcuts_data->GetPtTOFMatching()){//use TPC GeantFlukaCorrection
      Float_t FlukaCorrKPos=hGeantFlukaKPos->GetBinContent(hGeantFlukaKPos->FindBin(Spectra[1]->GetBinCenter(binK)));
      Float_t FlukaCorrKNeg=hGeantFlukaKNeg->GetBinContent(hGeantFlukaKNeg->FindBin(Spectra[4]->GetBinCenter(binK)));
      Printf("TPC Geant/Fluka: pt:%f  Pos:%f  Neg:%f",Spectra[1]->GetBinCenter(binK),FlukaCorrKPos,FlukaCorrKNeg);
      Spectra[1]->SetBinContent(binK,Spectra[1]->GetBinContent(binK)*FlukaCorrKPos);
      Spectra[4]->SetBinContent(binK,Spectra[4]->GetBinContent(binK)*FlukaCorrKNeg);
      Spectra[1]->SetBinError(binK,Spectra[1]->GetBinError(binK)*FlukaCorrKPos);
      Spectra[4]->SetBinError(binK,Spectra[4]->GetBinError(binK)*FlukaCorrKNeg);
      gGFCorrectionKaonPlus->SetPoint(binK,Spectra[1]->GetBinCenter(binK),FlukaCorrKPos);
      gGFCorrectionKaonMinus->SetPoint(binK,Spectra[4]->GetBinCenter(binK),FlukaCorrKNeg);
    }else{
      gGFCorrectionKaonPlus->SetPoint(binK,Spectra[1]->GetBinCenter(binK),0);
      gGFCorrectionKaonMinus->SetPoint(binK,Spectra[4]->GetBinCenter(binK),0);
      Float_t FlukaCorrKPosTracking=fGFKPosTracking->Eval(Spectra[1]->GetBinCenter(binK));
      Float_t FlukaCorrKNegTracking=fGFKNegTracking->Eval(Spectra[1]->GetBinCenter(binK));
      Printf("TPC/TOF Geant/Fluka Tracking: pt:%f  Pos:%f  Neg:%f",Spectra[1]->GetBinCenter(binK),FlukaCorrKPosTracking,FlukaCorrKNegTracking);
      Spectra[1]->SetBinContent(binK,Spectra[1]->GetBinContent(binK)*FlukaCorrKPosTracking);
      Spectra[4]->SetBinContent(binK,Spectra[4]->GetBinContent(binK)*FlukaCorrKNegTracking);
      Spectra[1]->SetBinError(binK,Spectra[1]->GetBinError(binK)*FlukaCorrKPosTracking);
      Spectra[4]->SetBinError(binK,Spectra[4]->GetBinError(binK)*FlukaCorrKNegTracking);
      Float_t FlukaCorrKPosMatching=fGFKPosMatching->Eval(Spectra[1]->GetBinCenter(binK));
      Float_t FlukaCorrKNegMatching=fGFKNegMatching->Eval(Spectra[1]->GetBinCenter(binK));
      Printf("TPC/TOF Geant/Fluka Matching: pt:%f  Pos:%f  Neg:%f",Spectra[1]->GetBinCenter(binK),FlukaCorrKPosMatching,FlukaCorrKNegMatching);
      Spectra[1]->SetBinContent(binK,Spectra[1]->GetBinContent(binK)*FlukaCorrKPosMatching);
      Spectra[4]->SetBinContent(binK,Spectra[4]->GetBinContent(binK)*FlukaCorrKNegMatching);
      Spectra[1]->SetBinError(binK,Spectra[1]->GetBinError(binK)*FlukaCorrKPosMatching);
      Spectra[4]->SetBinError(binK,Spectra[4]->GetBinError(binK)*FlukaCorrKNegMatching);
    }
  }
  
  //Geant Fluka for P in TPC
  Printf("\nGF correction for Protons");
  const Int_t kNCharge=2;
  Int_t kPos=0;
  Int_t kNeg=1;
  TFile* fGFProtons = new TFile ("GFCorrection/correctionForCrossSection.root");
  TH2D * hCorrFluka[kNCharge];
  TH2D * hCorrFluka[2];
  hCorrFluka[kPos] = (TH2D*)fGFProtons->Get("gHistCorrectionForCrossSectionProtons");
  hCorrFluka[kNeg] = (TH2D*)fGFProtons->Get("gHistCorrectionForCrossSectionAntiProtons");
  TGraph *gGFCorrectionProtonPlus=new TGraph();
  gGFCorrectionProtonPlus->SetName("gGFCorrectionProtonPlus");
  gGFCorrectionProtonPlus->SetTitle("gGFCorrectionProtonPlus");
  TGraph *gGFCorrectionProtonMinus=new TGraph();
  gGFCorrectionProtonMinus->SetName("gGFCorrectionProtonMinus");
  gGFCorrectionProtonMinus->SetTitle("gGFCorrectionProtonMinus");
  //getting GF func for Protons with TPCTOF
  TF1 *fGFpPosTracking;
  fGFpPosTracking = TrackingEff_geantflukaCorrection(4,kPositive);
  TF1 *fGFpNegTracking;
  fGFpNegTracking = TrackingEff_geantflukaCorrection(4,kNegative);
  TF1 *fGFpPosMatching;
  fGFpPosMatching = TOFmatchMC_geantflukaCorrection(4,kPositive);
  TF1 *fGFpNegMatching;
  fGFpNegMatching = TOFmatchMC_geantflukaCorrection(4,kNegative);
  
  for(Int_t icharge = 0; icharge < kNCharge; icharge++){
    Int_t nbins = Spectra[2]->GetNbinsX();
    Int_t nbinsy=hCorrFluka[icharge]->GetNbinsY();
    for(Int_t ibin = 0; ibin < nbins; ibin++){
      if(Spectra[2]->GetBinCenter(ibin)<tcuts_data->GetPtTOFMatching()){//use TPC GeantFlukaCorrection
  	Float_t pt = Spectra[2]->GetBinCenter(ibin);
  	Float_t minPtCorrection = hCorrFluka[icharge]->GetYaxis()->GetBinLowEdge(1);
  	Float_t maxPtCorrection = hCorrFluka[icharge]->GetYaxis()->GetBinLowEdge(nbinsy+1);
  	if (pt < minPtCorrection) pt = minPtCorrection+0.0001;
  	if (pt > maxPtCorrection) pt = maxPtCorrection;
  	Float_t correction = hCorrFluka[icharge]->GetBinContent(1,hCorrFluka[icharge]->GetYaxis()->FindBin(pt));
  	//cout<<correction<<"     charge "<<icharge<<endl;
  	if(icharge==0){
  	  if (correction != 0) {// If the bin is empty this is a  0
  	    Spectra[2]->SetBinContent(ibin,Spectra[2]->GetBinContent(ibin)*correction);
  	    Spectra[2]->SetBinError(ibin,Spectra[2]->GetBinError  (ibin)*correction);
  	    gGFCorrectionProtonPlus->SetPoint(ibin,pt,correction);
  	  }else if (Spectra[2]->GetBinContent(ibin) > 0) { // If we are skipping a non-empty bin, we notify the user
  	    cout << "Fluka/GEANT: Not correcting bin "<<ibin << " for protons secondaries, Pt:"<< pt<< endl;
  	    cout << " Bin content: " << Spectra[2]->GetBinContent(ibin)  << endl;
  	  }
  	}
  	if(icharge==1){
  	  if (correction != 0) {// If the bin is empty this is a  0
  	    Spectra[5]->SetBinContent(ibin,Spectra[5]->GetBinContent(ibin)*correction);
  	    Spectra[5]->SetBinError(ibin,Spectra[5]->GetBinError  (ibin)*correction);
  	    gGFCorrectionProtonMinus->SetPoint(ibin,pt,correction);
  	  }else if (Spectra[5]->GetBinContent(ibin) > 0) { // If we are skipping a non-empty bin, we notify the user
  	    cout << "Fluka/GEANT: Not correcting bin "<<ibin << " for Antiprotons secondaries, Pt:"<< pt<< endl;
  	    cout << " Bin content: " << Spectra[5]->GetBinContent(ibin)  << endl;
  	  }
  	}
      }else{
  	gGFCorrectionProtonPlus->SetPoint(ibin,Spectra[2]->GetBinCenter(ibin),0);
  	gGFCorrectionProtonMinus->SetPoint(ibin,Spectra[5]->GetBinCenter(ibin),0);
  	Float_t FlukaCorrpPosTracking=fGFpPosTracking->Eval(Spectra[2]->GetBinCenter(ibin));
  	Float_t FlukaCorrpNegTracking=fGFpNegTracking->Eval(Spectra[2]->GetBinCenter(ibin));
  	Printf("TPC/TOF Geant/Fluka Tracking: pt:%f  Pos:%f  Neg:%f",Spectra[2]->GetBinCenter(ibin),FlukaCorrpPosTracking,FlukaCorrpNegTracking);
  	if(icharge==0){
	  Spectra[2]->SetBinContent(ibin,Spectra[2]->GetBinContent(ibin)*FlukaCorrpPosTracking);
	  Spectra[2]->SetBinError(ibin,Spectra[2]->GetBinError(ibin)*FlukaCorrpPosTracking);
  	}else{
	  Spectra[5]->SetBinContent(ibin,Spectra[5]->GetBinContent(ibin)*FlukaCorrpNegTracking);
	  Spectra[5]->SetBinError(ibin,Spectra[5]->GetBinError(ibin)*FlukaCorrpNegTracking);
  	}
	Float_t FlukaCorrpPosMatching=fGFpPosMatching->Eval(Spectra[2]->GetBinCenter(ibin));
  	Float_t FlukaCorrpNegMatching=fGFpNegMatching->Eval(Spectra[2]->GetBinCenter(ibin));
  	Printf("TPC/TOF Geant/Fluka Matching: pt:%f  Pos:%f  Neg:%f",Spectra[2]->GetBinCenter(ibin),FlukaCorrpPosMatching,FlukaCorrpNegMatching);
	if(icharge==0){
	  Spectra[2]->SetBinContent(ibin,Spectra[2]->GetBinContent(ibin)*FlukaCorrpPosMatching);
	  Spectra[2]->SetBinError(ibin,Spectra[2]->GetBinError(ibin)*FlukaCorrpPosMatching);
	}else{
	  Spectra[5]->SetBinContent(ibin,Spectra[5]->GetBinContent(ibin)*FlukaCorrpNegMatching);
	  Spectra[5]->SetBinError(ibin,Spectra[5]->GetBinError(ibin)*FlukaCorrpNegMatching);
	}
      }
    }//end loop on bins	
  }
  gGFCorrectionKaonPlus->SetLineColor(kRed);
  gGFCorrectionKaonMinus->SetLineColor(kRed+2);
  gGFCorrectionProtonPlus->SetLineColor(kGreen);
  gGFCorrectionProtonMinus->SetLineColor(kGreen+2);
  fGFKPosTracking->SetLineColor(kRed);
  fGFKNegTracking->SetLineColor(kRed+2);
  fGFKPosMatching->SetLineColor(kRed);
  fGFKNegMatching->SetLineColor(kRed+2);
  fGFpPosTracking->SetLineColor(kGreen);
  fGFpNegTracking->SetLineColor(kGreen+2);
  fGFpPosMatching->SetLineColor(kGreen);
  fGFpNegMatching->SetLineColor(kGreen+2);
  fGFKPosTracking->SetLineStyle(2);
  fGFKNegTracking->SetLineStyle(2);
  fGFKPosMatching->SetLineStyle(3);
  fGFKNegMatching->SetLineStyle(3);
  fGFpPosTracking->SetLineStyle(2);
  fGFpNegTracking->SetLineStyle(2);
  fGFpPosMatching->SetLineStyle(3);
  fGFpNegMatching->SetLineStyle(3);
  fGFKPosTracking->SetRange(.6,5);
  fGFKNegTracking->SetRange(.6,5);
  fGFKPosMatching->SetRange(.6,5);
  fGFKNegMatching->SetRange(.6,5);
  fGFpPosTracking->SetRange(.6,5);
  fGFpNegTracking->SetRange(.6,5);
  fGFpPosMatching->SetRange(.6,5);
  fGFpNegMatching->SetRange(.6,5);
 
  TCanvas * cGFCorrection = new TCanvas ("cGFCorrection","cGFCorrection",700,500);
  gPad->SetGridx();
  gPad->SetGridy();
  gGFCorrectionKaonPlus->DrawClone("al");
  gGFCorrectionKaonMinus->DrawClone("lsame");
  gGFCorrectionProtonPlus->DrawClone("lsame");
  gGFCorrectionProtonMinus->DrawClone("lsame");
  fGFKPosTracking->DrawClone("lsame");
  fGFKNegTracking->DrawClone("lsame");
  fGFKPosMatching->DrawClone("lsame");
  fGFKNegMatching->DrawClone("lsame");
  fGFpPosTracking->DrawClone("lsame");
  fGFpNegTracking->DrawClone("lsame");
  fGFpPosMatching->DrawClone("lsame");
  fGFpNegMatching->DrawClone("lsame");
  gPad->BuildLegend()->SetFillColor(0);
  fout->cd();
  cGFCorrection->Write();
}

///////////
TF1 * TrackingEff_geantflukaCorrection(Int_t ipart, Int_t icharge)
{

  if (ipart == 3 && icharge == kNegative) {
    TF1 *f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge]), "TrackingPtGeantFlukaCorrectionKaMinus(x)", 0., 5.);
    return f;
  }
  else if (ipart == 4 && icharge == kNegative) {
    TF1 *f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge]), "TrackingPtGeantFlukaCorrectionPrMinus(x)", 0., 5.);
  }
  else
    TF1 *f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge]), "TrackingPtGeantFlukaCorrectionNull(x)", 0., 5.);

  return f;
}

Double_t TrackingPtGeantFlukaCorrectionNull(Double_t pTmc)
{
  return 1.;
}

Double_t TrackingPtGeantFlukaCorrectionPrMinus(Double_t pTmc)
{
  return (1 - 0.129758 *TMath::Exp(-pTmc*0.679612));
}

Double_t TrackingPtGeantFlukaCorrectionKaMinus(Double_t pTmc)
{
  return TMath::Min((0.972865 + 0.0117093*pTmc), 1.);
}

///////////////////////////////////////////
TF1 * TOFmatchMC_geantflukaCorrection(Int_t ipart, Int_t icharge)
{
  if (ipart == 3 && icharge == kNegative) {
    TF1 *f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge]), "MatchingPtGeantFlukaCorrectionKaMinus(x)", 0., 5.);
    return f;
  }
  else if (ipart == 4 && icharge == kNegative) {
    TF1 *f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge]), "MatchingPtGeantFlukaCorrectionPrMinus(x)", 0., 5.);
  }
  else
    TF1 *f = new TF1(Form("fGeantFluka_%s_%s", AliPID::ParticleName(ipart), Sign[icharge]), "MatchingPtGeantFlukaCorrectionNull(x)", 0., 5.);

  return f;
}

Double_t MatchingPtGeantFlukaCorrectionNull(Double_t pTmc)
{
  return 1.;
}

Double_t MatchingPtGeantFlukaCorrectionPrMinus(Double_t pTmc)
{
  Float_t ptTPCoutP =pTmc*(1-6.81059e-01*TMath::Exp(-pTmc*4.20094));
  return (TMath::Power(1 - 0.129758*TMath::Exp(-ptTPCoutP*0.679612),0.07162/0.03471));
}

Double_t MatchingPtGeantFlukaCorrectionKaMinus(Double_t pTmc)
{
  Float_t ptTPCoutK=pTmc*(1- 3.37297e-03/pTmc/pTmc - 3.26544e-03/pTmc);
  return TMath::Min((TMath::Power(0.972865 + 0.0117093*ptTPCoutK,0.07162/0.03471)), 1.);
}
