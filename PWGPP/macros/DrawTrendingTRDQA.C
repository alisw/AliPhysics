/*
 author: Ionut-Cristian Arsene
 email: i.c.arsene@cern.ch
 
 Macro to plot trends from an input root file containing a tree.
 Usage:
 DrawTrendingTRDQA("trendingFile.root")
 
 How to add a new trending plot:
 1. If a new trending variables needs to be added, its name must be added
    in the gTrendNames array and also a corresponding entry in the Trends enumeration.
    The name must match the name of its corresponding branch in the trending tree.
    If no new trending variable is added, proceed to step 2
 2. Define the needed histogram in DefineHistograms()
   2.1 If just a plot of the trending variable as a function of run number is needed, then
       a) define the histogram, 
       b) set the uniqueID for the TH1D (with the corresponding Trends index needed to be filled) 
       c) set the uniqueID for the y axis which will represent the variable used as error bar (optional)
       d) add it to the output list
       The histogram will then be filled and drawn automatically.
   2.2 If the plot is more complicated then defined at 2.1 (e.g. use quantities derived from the existing trends, 
                                                           correlations between various trends, etc.):
       a) define the histogram
       b) add it to the output list
       c) fill the histogram in the FillHistograms() function in the section for exceptions.
       d) drawing will be done automatically
 */


// Trending variable indices
enum Trends {
  kTPCTRDmatchEffPosAll=0, kTPCTRDmatchEffPosAllErr,
  kTPCTRDmatchEffNegAll, kTPCTRDmatchEffNegAllErr,
  kTRDTOFmatchEffPosAll, kTRDTOFmatchEffPosAllErr,
  kTRDTOFmatchEffNegAll, kTRDTOFmatchEffNegAllErr,
  kAvTRDtrkltsPerTrackAll, kAvTRDtrkltsPerTrackAllErr,
  kAvNclsPerTrackAll, kAvNclsPerTrackAllErr,
  kPHplateauHeight, kPHplateauHeightErr,
  kPHplateauSlope, kPHplateauSlopeErr,
  kQtotLandauMPV1GeVAll, kQtotLandauWidth1GeVAll,
  kPHplateauHeightAbsolute, kPHplateauHeightErrAbsolute,
  kPHplateauSlopeAbsolute, kPHplateauSlopeErrAbsolute,
  kQtotLandauMPV1GeVAllAbsolute, kQtotLandauWidth1GeVAllAbsolute,
  kTRDcheckDET_NTracksEvent, kTRDcheckDET_NTracksEventRMS,
  kTRDcheckDET_NTracksSector, 
  kTRDcheckDET_NClustersTrack, kTRDcheckDET_NClustersTrackRMS, 
  kTRDcheckDET_NClustersTracklet, kTRDcheckDET_NClustersTrackletRMS, 
  kTRDcheckDET_NTrackletsTrack, kTRDcheckDET_NTrackletsTrackRMS, 
  kTRDcheckDET_ChargeCluster, kTRDcheckDET_ChargeClusterRMS, 
  kTRDcheckDET_ChargeTracklet, kTRDcheckDET_ChargeTrackletRMS, 
  kTRDcheckDET_PHplateau, kTRDcheckDET_PHslope, kTRDcheckDET_PHamplificationPeak,
  kMeanExB, kRmsExB,
  kMeanGainFactor, kRmsGainFactor,
  kMeanT0, kRmsT0,
  kMeanVdrift, kRmsVdrift,
  kBeamIntensityA, kBeamIntensityC,
  kNtrends,
  kRun
};

// Trending variable names
// These names must match the names of the branches in the trending tree
// since branches are detected automatically based on the strings in this array
TString gTrendNames[kNtrends] = {
  "TPCTRDmatchEffPosAll", "TPCTRDmatchEffPosAllErr",
  "TPCTRDmatchEffNegAll", "TPCTRDmatchEffNegAllErr",
  "TRDTOFmatchEffPosAll", "TRDTOFmatchEffPosAllErr",
  "TRDTOFmatchEffNegAll", "TRDTOFmatchEffNegAllErr",
  "AvTRDtrkltsPerTrackAll", "AvTRDtrkltsPerTrackAllErr",
  "AvNclsPerTrackAll", "AvNclsPerTrackAllErr",
  "PHplateauHeight", "PHplateauHeightErr",
  "PHplateauSlope", "PHplateauSlopeErr",
  "QtotLandauMPV1GeVAll", "QtotLandauWidth1GeVAll",
  "PHplateauHeightAbsolute", "PHplateauHeightErrAbsolute",
  "PHplateauSlopeAbsolute", "PHplateauSlopeErrAbsolute",
  "QtotLandauMPV1GeVAllAbsolute", "QtotLandauWidth1GeVAllAbsolute",
  "TRDcheckDET_NTracksEvent", "TRDcheckDET_NTracksEventRMS",
  "TRDcheckDET_NTracksSector", 
  "TRDcheckDET_NClustersTrack", "TRDcheckDET_NClustersTrackRMS", 
  "TRDcheckDET_NClustersTracklet", "TRDcheckDET_NClustersTrackletRMS", 
  "TRDcheckDET_NTrackletsTrack", "TRDcheckDET_NTrackletsTrackRMS", 
  "TRDcheckDET_ChargeCluster", "TRDcheckDET_ChargeClusterRMS", 
  "TRDcheckDET_ChargeTracklet", "TRDcheckDET_ChargeTrackletRMS", 
  "TRDcheckDET_PHplateau", "TRDcheckDET_PHslope", "TRDcheckDET_PHamplificationPeak",
  "meanExB", "rmsExB",
  "meanGainFactor", "rmsGainFactor",
  "meanT0", "rmsT0",
  "meanVdrift", "rmsVdrift",
  "beamIntensityA", "beamIntensityC"
};


// prototypes
void DefineHistograms(TList* outList, Int_t nRuns);
void FillHistograms(Int_t irun, TList* hList, Double_t* values, Bool_t* branchFound);
void DrawAllHistograms(TList* hList);
void SetDrawStyle(TObject* obj, TString drawOption="E1",
	          Int_t markerStyle=24, Int_t markerColor=4, Double_t markerSize=2.0,
	          Int_t lineStyle=1, Int_t lineColor=4, Double_t lineWidth=2.0);

//________________________________________________________________________
void DrawTrendingTRDQA(TString trendingFilename="trending.root") {
  //
  // Draw the TRD QA trending
  //
  TFile* trendingFile = TFile::Open(trendingFilename.Data(), "READ");
  TTree* tree = (TTree*)trendingFile->Get("trending");
  if(!tree){
    cout << "E-DrawTrendingTRDQA.C: Cannot get the trending tree!" << endl;
    return;
  }
    
  Int_t run=0;
  tree->SetBranchAddress("run",&run);
  
  // Array which will hold the information from the tree
  // Note that its size is kNtrends+1. The extra element is reserved for the run number and 
  // is always located at the back of the array (element at index kNtrends)
  Double_t trends[kNtrends+1]={0.0};
  
  // Detect branches in the tree and assign an address
  Bool_t branchFound[kNtrends]; for(Int_t i=0;i<kNtrends;++i) branchFound[i]=kFALSE;
  for(Int_t i=0;i<kNtrends;++i) {
    TBranch* branch = tree->FindBranch(gTrendNames[i].Data());
    if(!branch) continue;
    branch->SetAddress(&trends[i]);
    branchFound[i]=kTRUE;
  }
  
  // Define histograms and add them to a list
  Int_t nRuns = tree->GetEntries();
  TList histList;
  histList.SetOwner(kTRUE);
  DefineHistograms(&histList, nRuns);
  
  // loop over all entries in the tree (one entry per run)
  for(Int_t i=0; i<nRuns; ++i) {
    tree->GetEntry(i);
    trends[kNtrends] = run;
    FillHistograms(i, &histList, trends, branchFound);
  }  // end loop over runs
  
  // Draw trending histograms
  // Here the canvases with one drawn histogram each are saved as PNG files.
  // The names of the PNG files use the name of the histogram
  DrawAllHistograms(&histList);
  
  // save histograms into a root file
  TFile* save=new TFile("PeriodTRDQAtrends.hist.root", "RECREATE");
  histList.Write();
  save->Close();
}


//________________________________________________________________________
void DrawAllHistograms(TList* l) {
  //
  // Draw all histograms in the list
  //
  if(!l) return;
  TCanvas* c=0x0;
  for(Int_t i=0; i<l->GetEntries(); ++i) {
    TObject* o=l->At(i);
    c=new TCanvas(Form("canvas_%s",o->GetName()), Form("%s", o->GetName()), 800.,600.);
    c->SetLeftMargin(0.15); c->SetBottomMargin(0.15); c->SetTopMargin(0.08); c->SetRightMargin(0.03);
    ((TH1*)o)->SetStats(kFALSE);
    
    // some exceptions
    if(TString(o->GetName()).Contains("BeamIntensity")) c->SetLogy();
    
    o->Draw();
    c->Print(Form("%s.png", TString(o->GetName()).Remove(0,1).Data()));
    delete c;
  }
}


//________________________________________________________________________
void FillHistograms(Int_t irun, TList* l, Double_t* v, Bool_t* branchFound) {
  //
  // Fill all histograms in the list
  //
  if(!l) return;
  for(Int_t i=0; i<l->GetEntries(); ++i) {
    TObject* o=l->At(i);
    if(TString(o->IsA()->GetName()).Contains("TH1D")) {
      TH1D* h1=(TH1D*)o;
      h1->GetXaxis()->SetBinLabel(irun+1, Form("%.0f", v[kNtrends]));
      
      // Deal with the exceptions (no assigned space in the Trends enum)
      if(TString(h1->GetName()).Contains("hTPCTRDmatchChargeAsymm")) {
	if(branchFound[kTPCTRDmatchEffPosAll] && branchFound[kTPCTRDmatchEffNegAll]) {
	  Double_t p=v[kTPCTRDmatchEffPosAll], m=v[kTPCTRDmatchEffNegAll], 
	       dp=v[kTPCTRDmatchEffPosAllErr], dm=v[kTPCTRDmatchEffNegAllErr];
	  if(TMath::Abs(p+m)>1.0e-6) {
	    h1->SetBinContent(irun+1, 100.0*2.0*(p-m)/(p+m));
	    h1->SetBinError(irun+1, 100.0*4.0*TMath::Sqrt(m*m*dp*dp+p*p*dm*dm)/(p+m)/(p+m));
	  }
	}
	continue;
      }
      if(TString(h1->GetName()).Contains("hTRDTOFmatchChargeAsymm")) {
	if(branchFound[kTRDTOFmatchEffPosAll] && branchFound[kTRDTOFmatchEffNegAll]) {
	  Double_t p=v[kTRDTOFmatchEffPosAll], m=v[kTRDTOFmatchEffNegAll], 
	       dp=v[kTRDTOFmatchEffPosAllErr], dm=v[kTRDTOFmatchEffNegAllErr];
	  if(TMath::Abs(p+m)>1.0e-6) {
	    h1->SetBinContent(irun+1, 100.0*2.0*(p-m)/(p+m));
	    h1->SetBinError(irun+1, 100.0*4.0*TMath::Sqrt(m*m*dp*dp+p*p*dm*dm)/(p+m)/(p+m));
	  }
	}
	continue;
      }
      
      // Fill the rest of the histograms
      Int_t var = h1->GetUniqueID()-1;
      Int_t errY = h1->GetYaxis()->GetUniqueID()-1;
      if(var>=0 && branchFound[var]) 
	h1->SetBinContent(irun+1, v[var]);
      if(errY>=0 && branchFound[errY]) h1->SetBinError(irun+1, v[errY]);
      
    }  // end if(TH1D)
  }  // end loop over histograms
}


//________________________________________________________________________
void DefineHistograms(TList* outList, Int_t nRuns) {
  //
  // Define trending histograms
  //
  TH1D* hTPCTRDmatchPos=new TH1D("hTPCTRDmatchPos","TPC-TRD matching efficiency at pt=1GeV/c, positive tracks;run;matching efficiency",
                                 nRuns, 0., nRuns);
  hTPCTRDmatchPos->SetUniqueID(kTPCTRDmatchEffPosAll+1); hTPCTRDmatchPos->GetYaxis()->SetUniqueID(kTPCTRDmatchEffPosAllErr+1);
  SetDrawStyle(hTPCTRDmatchPos); outList->Add(hTPCTRDmatchPos);
  
  TH1D* hTPCTRDmatchNeg=new TH1D("hTPCTRDmatchNeg","TPC-TRD matching efficiency at pt=1GeV/c, negative tracks;run;matching efficiency",
                                 nRuns, 0., nRuns);
  hTPCTRDmatchNeg->SetUniqueID(kTPCTRDmatchEffNegAll+1); hTPCTRDmatchNeg->GetYaxis()->SetUniqueID(kTPCTRDmatchEffNegAllErr+1);
  SetDrawStyle(hTPCTRDmatchNeg); outList->Add(hTPCTRDmatchNeg);
  
  TH1D* hTPCTRDmatchChargeAsymm=new TH1D("hTPCTRDmatchChargeAsymm","TPC-TRD matching efficiency at pt=1GeV/c, charge asymmetry;run;percents",
                                 nRuns, 0., nRuns);
  SetDrawStyle(hTPCTRDmatchChargeAsymm); outList->Add(hTPCTRDmatchChargeAsymm);
  
  TH1D* hTRDTOFmatchPos=new TH1D("hTRDTOFmatchPos","TRD-TOF matching efficiency at pt=1GeV/c, positive tracks;run;matching efficiency",
                                 nRuns, 0., nRuns);
  hTRDTOFmatchPos->SetUniqueID(kTRDTOFmatchEffPosAll+1); hTRDTOFmatchPos->GetYaxis()->SetUniqueID(kTRDTOFmatchEffPosAllErr+1);
  SetDrawStyle(hTRDTOFmatchPos); outList->Add(hTRDTOFmatchPos);
  
  TH1D* hTRDTOFmatchNeg=new TH1D("hTRDTOFmatchNeg","TRD-TOF matching efficiency at pt=1GeV/c, negative tracks;run;matching efficiency",
                                 nRuns, 0., nRuns);
  hTRDTOFmatchNeg->SetUniqueID(kTRDTOFmatchEffNegAll+1); hTRDTOFmatchNeg->GetYaxis()->SetUniqueID(kTRDTOFmatchEffNegAllErr+1);
  SetDrawStyle(hTRDTOFmatchNeg); outList->Add(hTRDTOFmatchNeg);
  
  TH1D* hTRDTOFmatchChargeAsymm=new TH1D("hTRDTOFmatchChargeAsymm","TRD-TOF matching efficiency at pt=1GeV/c, charge asymmetry;run;percents",
                                 nRuns, 0., nRuns);
  SetDrawStyle(hTRDTOFmatchChargeAsymm); outList->Add(hTRDTOFmatchChargeAsymm);

  TH1D* hTrkltsPerTrack=new TH1D("hTrkltsPerTrack","Average no. tracklets per TRD track at pt=1GeV/c;run;#tracklets",
                                 nRuns, 0., nRuns);
  hTrkltsPerTrack->SetUniqueID(kAvTRDtrkltsPerTrackAll+1); hTrkltsPerTrack->GetYaxis()->SetUniqueID(kAvTRDtrkltsPerTrackAllErr+1);
  SetDrawStyle(hTrkltsPerTrack); outList->Add(hTrkltsPerTrack);
  
  TH1D* hClsPerTrack=new TH1D("hClsPerTrack","Average no. clusters per TRD track at pt=1GeV/c;run;#tracklets",
                              nRuns, 0., nRuns);
  SetDrawStyle(hClsPerTrack); outList->Add(hClsPerTrack);
  hClsPerTrack->SetUniqueID(kAvNclsPerTrackAll+1); hClsPerTrack->GetYaxis()->SetUniqueID(kAvNclsPerTrackAllErr+1);
  
  TH1D* hPHplateauHeight=new TH1D("hPHplateauHeight","PH plateau height from slices;run;PH",
                                  nRuns, 0., nRuns);
  hPHplateauHeight->SetUniqueID(kPHplateauHeightAbsolute+1); hPHplateauHeight->GetYaxis()->SetUniqueID(kPHplateauHeightErrAbsolute+1);
  SetDrawStyle(hPHplateauHeight); outList->Add(hPHplateauHeight);
  
  TH1D* hPHplateauSlope=new TH1D("hPHplateauSlope","PH plateau slope from slices;run;slope",
                                  nRuns, 0., nRuns);
  hPHplateauSlope->SetUniqueID(kPHplateauSlopeAbsolute+1); hPHplateauSlope->GetYaxis()->SetUniqueID(kPHplateauSlopeErrAbsolute+1);
  SetDrawStyle(hPHplateauSlope); outList->Add(hPHplateauSlope);
  
  TH1D* hQtotLandauMPV=new TH1D("hQtotLandauMPV","Landau MPV for the tracklet charge distribution (Q_{tot}) ;run;Q_{tot}",
                                  nRuns, 0., nRuns);
  hQtotLandauMPV->SetUniqueID(kQtotLandauMPV1GeVAllAbsolute+1);
  SetDrawStyle(hQtotLandauMPV); outList->Add(hQtotLandauMPV);
  
  TH1D* hQtotLandauWidth=new TH1D("hQtotLandauWidth","Landau width for the tracklet charge distribution (Q_{tot});run;Q_{tot}",
                                  nRuns, 0., nRuns);
  hQtotLandauWidth->SetUniqueID(kQtotLandauWidth1GeVAllAbsolute+1);
  SetDrawStyle(hQtotLandauWidth); outList->Add(hQtotLandauWidth);
  
  TH1D* hNTracksPerEvent_DET=new TH1D("hNTracksPerEvent_DET", "Number of TRD tracks per event, (TRDcheckDET);run;#tracks",
                                  nRuns, 0., nRuns);
  hNTracksPerEvent_DET->SetUniqueID(kTRDcheckDET_NTracksEvent+1); hNTracksPerEvent_DET->GetYaxis()->SetUniqueID(kTRDcheckDET_NTracksEventRMS+1);
  SetDrawStyle(hNTracksPerEvent_DET); outList->Add(hNTracksPerEvent_DET);
  
  TH1D* hNTracksPerSector_DET=new TH1D("hNTracksPerSector_DET", "Number of TRD tracks per sector (TRDcheckDET);run;#tracks",
                                  nRuns, 0., nRuns);
  hNTracksPerSector_DET->SetUniqueID(kTRDcheckDET_NTracksSector+1);
  SetDrawStyle(hNTracksPerSector_DET); outList->Add(hNTracksPerSector_DET);
  
  TH1D* hNClustersPerTrack_DET=new TH1D("hNClustersPerTrack_DET", "Number of clusters per track (TRDcheckDET);run;#clusters",
                                  nRuns, 0., nRuns);
  hNClustersPerTrack_DET->SetUniqueID(kTRDcheckDET_NClustersTrack+1); hNClustersPerTrack_DET->GetYaxis()->SetUniqueID(kTRDcheckDET_NClustersTrackRMS+1);
  SetDrawStyle(hNClustersPerTrack_DET); outList->Add(hNClustersPerTrack_DET);
  
  TH1D* hNClustersPerTracklet_DET=new TH1D("hNClustersPerTracklet_DET", "Number of clusters per tracklet (TRDcheckDET);run;#clusters",
                                  nRuns, 0., nRuns);
  hNClustersPerTracklet_DET->SetUniqueID(kTRDcheckDET_NClustersTracklet+1); hNClustersPerTracklet_DET->GetYaxis()->SetUniqueID(kTRDcheckDET_NClustersTrackletRMS+1);
  SetDrawStyle(hNClustersPerTracklet_DET); outList->Add(hNClustersPerTracklet_DET);

  TH1D* hNTrackletsPerTrack_DET=new TH1D("hNTrackletsPerTrack_DET", "Number of tracklets per track (TRDcheckDET);run;#tracklets",
                                  nRuns, 0., nRuns);
  hNTrackletsPerTrack_DET->SetUniqueID(kTRDcheckDET_NTrackletsTrack+1); hNTrackletsPerTrack_DET->GetYaxis()->SetUniqueID(kTRDcheckDET_NTrackletsTrackRMS+1);
  SetDrawStyle(hNTrackletsPerTrack_DET); outList->Add(hNTrackletsPerTrack_DET);
  
  TH1D* hClusterCharge_DET=new TH1D("hClusterCharge_DET", "Cluster charge (TRDcheckDET);run;#charge",
                                  nRuns, 0., nRuns);
  hClusterCharge_DET->SetUniqueID(kTRDcheckDET_ChargeCluster+1); hClusterCharge_DET->GetYaxis()->SetUniqueID(kTRDcheckDET_ChargeClusterRMS+1);
  SetDrawStyle(hClusterCharge_DET); outList->Add(hClusterCharge_DET);
  
  TH1D* hTrackletCharge_DET=new TH1D("hTrackletCharge_DET", "Tracklet charge (TRDcheckDET);run;#charge",
                                  nRuns, 0., nRuns);
  hTrackletCharge_DET->SetUniqueID(kTRDcheckDET_ChargeTracklet+1); hTrackletCharge_DET->GetYaxis()->SetUniqueID(kTRDcheckDET_ChargeTrackletRMS+1);
  SetDrawStyle(hTrackletCharge_DET); outList->Add(hTrackletCharge_DET);
  
  TH1D* hPHplateau_DET=new TH1D("hPHplateau_DET", "PH plateau (TRDcheckDET);run;#charge",
                                  nRuns, 0., nRuns);
  hPHplateau_DET->SetUniqueID(kTRDcheckDET_PHplateau+1);
  SetDrawStyle(hPHplateau_DET); outList->Add(hPHplateau_DET);
  
  TH1D* hPHslope_DET=new TH1D("hPHslope_DET", "PH slope (TRDcheckDET);run;#charge/timebin",
                                  nRuns, 0., nRuns);
  hPHslope_DET->SetUniqueID(kTRDcheckDET_PHslope+1);
  SetDrawStyle(hPHslope_DET); outList->Add(hPHslope_DET);
  
  TH1D* hPHamplifPeak_DET=new TH1D("hPHamplifPeak_DET", "PH amplification peak position (TRDcheckDET);run;#timebin",
                                  nRuns, 0., nRuns);
  hPHamplifPeak_DET->SetUniqueID(kTRDcheckDET_PHamplificationPeak+1);
  SetDrawStyle(hPHamplifPeak_DET); outList->Add(hPHamplifPeak_DET);
  
  TH1D* hMeanExB=new TH1D("hMeanExB", "Mean ExB over all chambers from OCDB;run;ExB correction",
                                  nRuns, 0., nRuns);
  hMeanExB->SetUniqueID(kMeanExB+1); hMeanExB->GetYaxis()->SetUniqueID(kRmsExB+1);
  SetDrawStyle(hMeanExB); outList->Add(hMeanExB);
  
  TH1D* hMeanGainFactor=new TH1D("hMeanGainFactor", "Mean gain factor over all chambers from OCDB;run;gain factor",
                                  nRuns, 0., nRuns);
  hMeanGainFactor->SetUniqueID(kMeanGainFactor+1); hMeanGainFactor->GetYaxis()->SetUniqueID(kRmsGainFactor+1);
  SetDrawStyle(hMeanGainFactor); outList->Add(hMeanGainFactor);
  
  TH1D* hMeanT0=new TH1D("hMeanT0", "Mean T0 over all chambers from OCDB;run;T0",
                                  nRuns, 0., nRuns);
  hMeanT0->SetUniqueID(kMeanT0+1); hMeanT0->GetYaxis()->SetUniqueID(kRmsT0+1);
  SetDrawStyle(hMeanT0); outList->Add(hMeanT0);
  
  TH1D* hMeanVdrift=new TH1D("hMeanVdrift", "Mean drift velocity over all chambers from OCDB;run;drift velocity",
                                  nRuns, 0., nRuns);
  hMeanVdrift->SetUniqueID(kMeanVdrift+1); hMeanVdrift->GetYaxis()->SetUniqueID(kRmsVdrift+1);
  SetDrawStyle(hMeanVdrift); outList->Add(hMeanVdrift);
  
  TH1D* hBeamIntensityA=new TH1D("hBeamIntensityA", "Beam A intensity from OCDB;run;beam intensity",
                                  nRuns, 0., nRuns);
  hBeamIntensityA->SetUniqueID(kBeamIntensityA+1);
  SetDrawStyle(hBeamIntensityA); outList->Add(hBeamIntensityA);
  
  TH1D* hBeamIntensityC=new TH1D("hBeamIntensityC", "Beam C intensity from OCDB;run;beam intensity",
                                  nRuns, 0., nRuns);
  hBeamIntensityC->SetUniqueID(kBeamIntensityC+1);
  SetDrawStyle(hBeamIntensityC); outList->Add(hBeamIntensityC);
}


//________________________________________________________________________
void SetDrawStyle(TObject* obj, 
		  TString drawOption /*="E1"*/,
	          Int_t markerStyle /*=24*/, Int_t markerColor /*=4*/, Double_t markerSize /*=2.0*/,
	          Int_t lineStyle /*=1*/, Int_t lineColor /*=4*/, Double_t lineWidth /*=2.0*/) {
  //
  // Set draw options
  //
  if(TString(obj->IsA()->GetName()).Contains("TH1D")) {
    TH1D* o=(TH1D*)obj;
    o->SetDrawOption(drawOption.Data());
    o->SetMarkerStyle(markerStyle);
    o->SetMarkerColor(markerColor);
    o->SetMarkerSize(markerSize);
    o->SetLineStyle(lineStyle);
    o->SetLineColor(lineColor);
    o->SetLineWidth(lineWidth);
    o->GetXaxis()->SetTitleOffset(1.6);
    o->GetYaxis()->SetNdivisions(507);
    o->GetXaxis()->SetLabelFont(42);
    o->GetXaxis()->SetTitleFont(42);
    o->GetYaxis()->SetLabelFont(42);
    o->GetYaxis()->SetTitleFont(42);
    o->SetTitleFont(42);
  }
}
