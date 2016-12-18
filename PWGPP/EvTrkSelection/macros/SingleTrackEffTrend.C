

// Connecting macro to calcualate QA eff
// Jitendra, Zaida

Int_t SingleTrackEffTrend(char *infile="", Int_t runN=117222) {

  gStyle->SetOptStat(0);
  if(!infile) return -1;
  TStopwatch timer;
  timer.Start();
    
  char input_file[300];
  printf("Running TrackEffTrend for file: %s and Run: %i \n \n",infile,runN);
  gROOT->LoadMacro("./CalcSingleTrackEffQA.C");

  //TGrid::Connect("alien://"); //May be for central need
  Int_t run = runN;
  
    AliCFContainer *LoadSlice1 = CalcSingleTrackEffQAC(infile,"NchFbit0");
    if(!LoadSlice1) {
        Printf("ERROR: (SingleTrackEffTrend.C) Directory or container not found. Skipping this file!");
        return -1;
    }
    CalcSingleTrackEffQAV(LoadSlice1, "pt", 0.5,  24.0,   2,  0,  "Nch",  kFALSE, runN);
    Float_t AvgpTeff_1_2_Nch_fbit0_type2   = CalcSingleTrackEffQAV(LoadSlice1,  "pt", 1.0,  2.0,   2,  0,  "Nch",  kTRUE, runN); //runN is not used here
    Float_t AvgpTeff_2_6_Nch_fbit0_type2   = CalcSingleTrackEffQAV(LoadSlice1,  "pt", 2.0,  6.0,   2,  0,  "Nch",  kTRUE, runN); //runN is not used here
    Float_t AvgpTeff_6_8_Nch_fbit0_type2   = CalcSingleTrackEffQAV(LoadSlice1,  "pt", 6.0,  8.0,   2,  0,  "Nch",  kTRUE, runN); //runN is not used here
    Float_t AvgpTeff_8_10_Nch_fbit0_type2  = CalcSingleTrackEffQAV(LoadSlice1,  "pt", 8.0,  10.0,  2,  0,  "Nch",  kTRUE, runN); //runN is not used here
    CalcSingleTrackEffQAV(LoadSlice1, "eta", -0.8,  0.8,   2,  0,  "Nch",  kFALSE, runN); //range not matters
    Float_t Avgetaeff_dot8_Nch_fbit0_type2    = CalcSingleTrackEffQAV(LoadSlice1, "eta", -0.8,  0.8,   2,  0,  "Nch",  kTRUE, runN); //runN is not used here
    Float_t Avgetaeff_dot5_Nch_fbit0_type2    = CalcSingleTrackEffQAV(LoadSlice1, "eta", -0.5,  0.5,   2,  0,  "Nch",  kTRUE, runN); //runN is not used here

    Printf("\n \n");
    AliCFContainer *LoadSlice2 = CalcSingleTrackEffQAC(infile,"PionFbit0");
    CalcSingleTrackEffQAV(LoadSlice2, "pt", 0.5,  24.0,   2,  0,  "Pion",  kFALSE, runN); // Range doen not matter
    Float_t AvgpTeff_1_2_Pion_fbit0_type2   = CalcSingleTrackEffQAV(LoadSlice2,  "pt", 1.0,  2.0,   2,  0,  "Pion",  kTRUE, runN); //runN is not used here
    Float_t AvgpTeff_2_6_Pion_fbit0_type2   = CalcSingleTrackEffQAV(LoadSlice2,  "pt", 2.0,  6.0,   2,  0,  "Pion",  kTRUE, runN); //runN is not used here
    Float_t AvgpTeff_6_8_Pion_fbit0_type2    = CalcSingleTrackEffQAV(LoadSlice2,  "pt", 6.0,  8.0,  2,  0,  "Pion",  kTRUE, runN); //runN is not used here
    Float_t AvgpTeff_8_10_Pion_fbit0_type2  = CalcSingleTrackEffQAV(LoadSlice2,  "pt", 8.0,  10.0,  2,  0,  "Pion",  kTRUE, runN); //runN is not used here
    CalcSingleTrackEffQAV(LoadSlice2, "eta", -0.8,  0.8,   2,  0,  "Pion",  kFALSE, runN); //range not matters
    Float_t Avgetaeff_dot8_Pion_fbit0_type2    = CalcSingleTrackEffQAV(LoadSlice2, "eta", -0.8,  0.8,   2,  0,  "Pion",  kTRUE, runN); //runN is not used here
    Float_t Avgetaeff_dot5_Pion_fbit0_type2    = CalcSingleTrackEffQAV(LoadSlice2, "eta", -0.5,  0.5,   2,  0,  "Pion",  kTRUE, runN); //runN is not used here

    Printf("\n \n");
    AliCFContainer *LoadSlice3 = CalcSingleTrackEffQAC(infile,"KaonFbit0");
    CalcSingleTrackEffQAV(LoadSlice3, "pt", 0.5,  24.0,   2,  0,  "Kaon",  kFALSE, runN); // Range doen not matter
    Float_t AvgpTeff_1_2_Kaon_fbit0_type2   = CalcSingleTrackEffQAV(LoadSlice3,  "pt", 1.0,  2.0,   2,  0,  "Kaon",  kTRUE, runN); //runN is not used here
    Float_t AvgpTeff_2_6_Kaon_fbit0_type2   = CalcSingleTrackEffQAV(LoadSlice3,  "pt", 2.0,  6.0,   2,  0,  "Kaon",  kTRUE, runN); //runN is not used here
    Float_t AvgpTeff_6_8_Kaon_fbit0_type2   = CalcSingleTrackEffQAV(LoadSlice3,  "pt", 6.0,  8.0,   2,  0,  "Kaon",  kTRUE, runN); //runN is not used here
    Float_t AvgpTeff_8_10_Kaon_fbit0_type2  = CalcSingleTrackEffQAV(LoadSlice3,  "pt", 8.0,  10.0,  2,  0,  "Kaon",  kTRUE, runN); //runN is not used here
    CalcSingleTrackEffQAV(LoadSlice3, "eta", -0.8,  0.8,   2,  0,  "Kaon",  kFALSE, runN); //range not matters
    Float_t Avgetaeff_dot8_Kaon_fbit0_type2    = CalcSingleTrackEffQAV(LoadSlice3, "eta", -0.8,  0.8,   2,  0,  "Kaon",  kTRUE, runN); //runN is not used here
    Float_t Avgetaeff_dot5_Kaon_fbit0_type2    = CalcSingleTrackEffQAV(LoadSlice3, "eta", -0.5,  0.5,   2,  0,  "Kaon",  kTRUE, runN); //runN is not used here

    Printf("\n \n");
    AliCFContainer *LoadSlice4 = CalcSingleTrackEffQAC(infile,"ElectronFbit0");
    CalcSingleTrackEffQAV(LoadSlice4, "pt", 0.5,  24.0,   2,  0,  "Electron",  kFALSE, runN); // Range doen not matter
    Float_t AvgpTeff_1_2_Elec_fbit0_type2   = CalcSingleTrackEffQAV(LoadSlice4,  "pt", 1.0,  2.0,   2,  0,  "Electron",  kTRUE, runN); //runN is not used here
    Float_t AvgpTeff_2_6_Elec_fbit0_type2   = CalcSingleTrackEffQAV(LoadSlice4,  "pt", 2.0,  6.0,   2,  0,  "Electron",  kTRUE, runN); //runN is not used here
    Float_t AvgpTeff_6_8_Elec_fbit0_type2    = CalcSingleTrackEffQAV(LoadSlice4,  "pt", 6.0,  8.0,  2,  0,  "Electron",  kTRUE, runN); //runN is not used here
    Float_t AvgpTeff_8_10_Elec_fbit0_type2  = CalcSingleTrackEffQAV(LoadSlice4,  "pt", 8.0,  10.0,  2,  0,  "Electron",  kTRUE, runN); //runN is not used here
    CalcSingleTrackEffQAV(LoadSlice4, "eta", -0.8,  0.8,   2,  0,  "Electron",  kFALSE, runN); //range not matters
    Float_t Avgetaeff_dot8_Elec_fbit0_type2    = CalcSingleTrackEffQAV(LoadSlice4, "eta", -0.8,  0.8,   2,  0,  "Electron",  kTRUE, runN); //runN is not used here
    Float_t Avgetaeff_dot5_Elec_fbit0_type2    = CalcSingleTrackEffQAV(LoadSlice4, "eta", -0.5,  0.5,   2,  0,  "Electron",  kTRUE, runN); //runN is not used here

    
    Printf("\n \n");
    AliCFContainer *LoadSlice5 = CalcSingleTrackEffQAC(infile,"ProtonFbit0");
    CalcSingleTrackEffQAV(LoadSlice5, "pt", 0.5,  24.0,   2,  0,  "Proton",  kFALSE, runN); // Range doen not matter
    Float_t AvgpTeff_1_2_Prot_fbit0_type2   = CalcSingleTrackEffQAV(LoadSlice5,  "pt", 1.0,  2.0,   2,  0,  "Proton",  kTRUE, runN); //runN is not used here
    Float_t AvgpTeff_2_6_Prot_fbit0_type2   = CalcSingleTrackEffQAV(LoadSlice5,  "pt", 2.0,  6.0,   2,  0,  "Proton",  kTRUE, runN); //runN is not used here
    Float_t AvgpTeff_6_8_Prot_fbit0_type2    = CalcSingleTrackEffQAV(LoadSlice5,  "pt", 6.0,  8.0,  2,  0,  "Proton",  kTRUE, runN); //runN is not used here
    Float_t AvgpTeff_8_10_Prot_fbit0_type2  = CalcSingleTrackEffQAV(LoadSlice5,  "pt", 8.0,  10.0,  2,  0,  "Proton",  kTRUE, runN); //runN is not used here
    CalcSingleTrackEffQAV(LoadSlice5, "eta", -0.8,  0.8,   2,  0,  "Proton",  kFALSE, runN); //range not matters
    Float_t Avgetaeff_dot8_Prot_fbit0_type2    = CalcSingleTrackEffQAV(LoadSlice5, "eta", -0.8,  0.8,   2,  0,  "Proton",  kTRUE, runN); //runN is not used here
    Float_t Avgetaeff_dot5_Prot_fbit0_type2    = CalcSingleTrackEffQAV(LoadSlice5, "eta", -0.5,  0.5,   2,  0,  "Proton",  kTRUE, runN); //runN is not used here
    
    
    Printf("\n \n");
    AliCFContainer *LoadSlice6 = CalcSingleTrackEffQAC(infile,"NchFbit4");
    CalcSingleTrackEffQAV(LoadSlice6, "pt", 0.5,  24.0,   2,  4,  "Nch",  kFALSE, runN);
    Float_t AvgpTeff_1_2_Nch_fbit4_type2    = CalcSingleTrackEffQAV(LoadSlice6,  "pt", 1.0,  2.0,   2,  4,  "NchFB4",  kTRUE, runN); //runN is not used here
    Float_t AvgpTeff_2_6_Nch_fbit4_type2    = CalcSingleTrackEffQAV(LoadSlice6,  "pt", 2.0,  6.0,   2,  4,  "NchFB4",  kTRUE, runN); //runN is not used here
    Float_t AvgpTeff_6_8_Nch_fbit4_type2    = CalcSingleTrackEffQAV(LoadSlice6,  "pt", 6.0,  8.0,   2,  4,  "NchFB4",  kTRUE, runN); //runN is not used here
    Float_t AvgpTeff_8_10_Nch_fbit4_type2   = CalcSingleTrackEffQAV(LoadSlice6,  "pt", 8.0,  10.0,  2,  4,  "NchFB4",  kTRUE, runN); //runN is not used here
    CalcSingleTrackEffQAV(LoadSlice6, "eta", -0.8,  0.8,   2,  0,  "Nch",  kFALSE, runN); //range not matters
    Float_t Avgetaeff_dot8_Nch_fbit4_type2    = CalcSingleTrackEffQAV(LoadSlice6, "eta", -0.8,  0.8,   2,  0,  "NchFB4",  kTRUE, runN); //runN is not used here
    Float_t Avgetaeff_dot5_Nch_fbit4_type2    = CalcSingleTrackEffQAV(LoadSlice6, "eta", -0.5,  0.5,   2,  0,  "NchFB4",  kTRUE, runN); //runN is not used here

    
    //Preparaing treding tree
    TFile* fouttrend =  TFile::Open("trending.root","recreate");
    TTree* teff = new TTree("trending","variables efficiency tree");
    teff->Branch("run",&run,"run/I");
    
    teff->Branch("AvgpTeff_1_2_Nch_fbit0_type2"   ,&AvgpTeff_1_2_Nch_fbit0_type2,   "AvgpTeff_1_2_Nch_fbit0_type2/F");
    teff->Branch("AvgpTeff_2_6_Nch_fbit0_type2"   ,&AvgpTeff_2_6_Nch_fbit0_type2,   "AvgpTeff_2_6_Nch_fbit0_type2/F");
    teff->Branch("AvgpTeff_6_8_Nch_fbit0_type2"   ,&AvgpTeff_6_8_Nch_fbit0_type2,   "AvgpTeff_6_8_Nch_fbit0_type2/F");
    teff->Branch("AvgpTeff_8_10_Nch_fbit0_type2"  ,&AvgpTeff_8_10_Nch_fbit0_type2,  "AvgpTeff_8_10_Nch_fbit0_type2/F");
    teff->Branch("Avgetaeff_dot5_Nch_fbit0_type2"   ,&Avgetaeff_dot5_Nch_fbit0_type2,   "Avgetaeff_dot5_Nch_fbit0_type2/F");
    teff->Branch("Avgetaeff_dot8_Nch_fbit0_type2"   ,&Avgetaeff_dot8_Nch_fbit0_type2,   "Avgetaeff_dot8_Nch_fbit0_type2/F");

    
    teff->Branch("AvgpTeff_1_2_Pion_fbit0_type2"   ,&AvgpTeff_1_2_Pion_fbit0_type2,   "AvgpTeff_1_2_Pion_fbit0_type2/F");
    teff->Branch("AvgpTeff_2_6_Pion_fbit0_type2"   ,&AvgpTeff_2_6_Pion_fbit0_type2,   "AvgpTeff_2_6_Pion_fbit0_type2/F");
    teff->Branch("AvgpTeff_6_8_Pion_fbit0_type2"   ,&AvgpTeff_6_8_Pion_fbit0_type2,   "AvgpTeff_6_8_Pion_fbit0_type2/F");
    teff->Branch("AvgpTeff_8_10_Pion_fbit0_type2"  ,&AvgpTeff_8_10_Pion_fbit0_type2,  "AvgpTeff_8_10_Pion_fbit0_type2/F");
    teff->Branch("Avgetaeff_dot5_Pion_fbit0_type2"   ,&Avgetaeff_dot5_Pion_fbit0_type2,   "Avgetaeff_dot5_Pion_fbit0_type2/F");
    teff->Branch("Avgetaeff_dot8_Pion_fbit0_type2"   ,&Avgetaeff_dot8_Pion_fbit0_type2,   "Avgetaeff_dot8_Pion_fbit0_type2/F");

    teff->Branch("AvgpTeff_1_2_Kaon_fbit0_type2"   ,&AvgpTeff_1_2_Kaon_fbit0_type2,   "AvgpTeff_1_2_Kaon_fbit0_type2/F");
    teff->Branch("AvgpTeff_2_6_Kaon_fbit0_type2"   ,&AvgpTeff_2_6_Kaon_fbit0_type2,   "AvgpTeff_2_6_Kaon_fbit0_type2/F");
    teff->Branch("AvgpTeff_6_8_Kaon_fbit0_type2"   ,&AvgpTeff_6_8_Kaon_fbit0_type2,   "AvgpTeff_6_8_Kaon_fbit0_type2/F");
    teff->Branch("AvgpTeff_8_10_Kaon_fbit0_type2"  ,&AvgpTeff_8_10_Kaon_fbit0_type2,  "AvgpTeff_8_10_Kaon_fbit0_type2/F");
    teff->Branch("Avgetaeff_dot5_Kaon_fbit0_type2"   ,&Avgetaeff_dot5_Kaon_fbit0_type2,   "Avgetaeff_dot5_Kaon_fbit0_type2/F");
    teff->Branch("Avgetaeff_dot8_Kaon_fbit0_type2"   ,&Avgetaeff_dot8_Kaon_fbit0_type2,   "Avgetaeff_dot8_Kaon_fbit0_type2/F");
    
    teff->Branch("AvgpTeff_1_2_Elec_fbit0_type2"   ,&AvgpTeff_1_2_Elec_fbit0_type2,   "AvgpTeff_1_2_Elec_fbit0_type2/F");
    teff->Branch("AvgpTeff_2_6_Elec_fbit0_type2"   ,&AvgpTeff_2_6_Elec_fbit0_type2,   "AvgpTeff_2_6_Elec_fbit0_type2/F");
    teff->Branch("AvgpTeff_6_8_Elec_fbit0_type2"   ,&AvgpTeff_6_8_Elec_fbit0_type2,   "AvgpTeff_6_8_Elec_fbit0_type2/F");
    teff->Branch("AvgpTeff_8_10_Elec_fbit0_type2"  ,&AvgpTeff_8_10_Elec_fbit0_type2,  "AvgpTeff_8_10_Elec_fbit0_type2/F");
    teff->Branch("Avgetaeff_dot5_Elec_fbit0_type2"   ,&Avgetaeff_dot5_Elec_fbit0_type2,   "Avgetaeff_dot5_Elec_fbit0_type2/F");
    teff->Branch("Avgetaeff_dot8_Elec_fbit0_type2"   ,&Avgetaeff_dot8_Elec_fbit0_type2,   "Avgetaeff_dot8_Elec_fbit0_type2/F");

    teff->Branch("AvgpTeff_1_2_Prot_fbit0_type2"   ,&AvgpTeff_1_2_Prot_fbit0_type2,   "AvgpTeff_1_2_Prot_fbit0_type2/F");
    teff->Branch("AvgpTeff_2_6_Prot_fbit0_type2"   ,&AvgpTeff_2_6_Prot_fbit0_type2,   "AvgpTeff_2_6_Prot_fbit0_type2/F");
    teff->Branch("AvgpTeff_6_8_Prot_fbit0_type2"   ,&AvgpTeff_6_8_Prot_fbit0_type2,   "AvgpTeff_6_8_Prot_fbit0_type2/F");
    teff->Branch("AvgpTeff_8_10_Prot_fbit0_type2"  ,&AvgpTeff_8_10_Prot_fbit0_type2,  "AvgpTeff_8_10_Prot_fbit0_type2/F");
    teff->Branch("Avgetaeff_dot5_Prot_fbit0_type2"   ,&Avgetaeff_dot5_Prot_fbit0_type2,   "Avgetaeff_dot5_Prot_fbit0_type2/F");
    teff->Branch("Avgetaeff_dot8_Prot_fbit0_type2"   ,&Avgetaeff_dot8_Prot_fbit0_type2,   "Avgetaeff_dot8_Prot_fbit0_type2/F");
    
    
    teff->Branch("AvgpTeff_1_2_Nch_fbit4_type2"   ,&AvgpTeff_1_2_Nch_fbit4_type2,   "AvgpTeff_1_2_Nch_fbit4_type2/F");
    teff->Branch("AvgpTeff_2_6_Nch_fbit4_type2"   ,&AvgpTeff_2_6_Nch_fbit4_type2,   "AvgpTeff_2_6_Nch_fbit4_type2/F");
    teff->Branch("AvgpTeff_6_8_Nch_fbit4_type2"   ,&AvgpTeff_6_8_Nch_fbit4_type2,   "AvgpTeff_6_8_Nch_fbit4_type2/F");
    teff->Branch("AvgpTeff_8_10_Nch_fbit4_type2"  ,&AvgpTeff_8_10_Nch_fbit4_type2,  "AvgpTeff_8_10_Nch_fbit4_type2/F");
    teff->Branch("Avgetaeff_dot5_Nch_fbit4_type2"   ,&Avgetaeff_dot5_Nch_fbit4_type2,   "Avgetaeff_dot5_Nch_fbit4_type2/F");
    teff->Branch("Avgetaeff_dot8_Nch_fbit4_type2"   ,&Avgetaeff_dot8_Nch_fbit4_type2,   "Avgetaeff_dot8_Nch_fbit4_type2/F");
    
    
    teff->Fill();
    teff->Write();
    
    timer.Stop();
    timer.Print();
    
    fouttrend->Close();
    
    return 0;

}