

//Period level eff macro
//Jitendra Kumar, Zaida

void periodLevelQAEff(TString inputFile ="trending.root"){
    
    gStyle->SetOptStat(0);
    gStyle->SetPadGridX(0);
    gStyle->SetPadGridY(1);
    
    TFile* fTrend = TFile::Open(inputFile.Data());
    TTree* tTrend = (TTree*) fTrend->Get("trending");
    Int_t nRuns = tTrend->GetEntries();
    
    Int_t run = 0;
    tTrend->SetBranchAddress("run",&run);
    
    Float_t AvgpTeff_1_2_Nch_fbit0_type2 = 0.,AvgpTeff_2_6_Nch_fbit0_type2 = 0.,AvgpTeff_6_8_Nch_fbit0_type2 = 0., AvgpTeff_8_10_Nch_fbit0_type2=0.;
    Float_t Avgetaeff_dot5_Nch_fbit0_type2 = 0.,Avgetaeff_dot8_Nch_fbit0_type2 = 0;
    tTrend->SetBranchAddress("AvgpTeff_1_2_Nch_fbit0_type2",  &AvgpTeff_1_2_Nch_fbit0_type2);
    tTrend->SetBranchAddress("AvgpTeff_2_6_Nch_fbit0_type2",  &AvgpTeff_2_6_Nch_fbit0_type2);
    tTrend->SetBranchAddress("AvgpTeff_6_8_Nch_fbit0_type2",  &AvgpTeff_6_8_Nch_fbit0_type2);
    tTrend->SetBranchAddress("AvgpTeff_8_10_Nch_fbit0_type2", &AvgpTeff_8_10_Nch_fbit0_type2);
    tTrend->SetBranchAddress("Avgetaeff_dot5_Nch_fbit0_type2",  &Avgetaeff_dot5_Nch_fbit0_type2);
    tTrend->SetBranchAddress("Avgetaeff_dot8_Nch_fbit0_type2", &Avgetaeff_dot8_Nch_fbit0_type2);
    TH1F* hTrendAvgPteff_1_2_Nch_fbit0_type2  = new TH1F("hTrendAvgPteff_1_2_Nch_fbit0_type2",  "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgPteff_2_6_Nch_fbit0_type2  = new TH1F("hTrendAvgPteff_2_6_Nch_fbit0_type2",  "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgPteff_6_8_Nch_fbit0_type2  = new TH1F("hTrendAvgPteff_6_8_Nch_fbit0_type2",  "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgPteff_8_10_Nch_fbit0_type2 = new TH1F("hTrendAvgPteff_8_10_Nch_fbit0_type2", "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgetaeff_dot5_Nch_fbit0_type2 = new TH1F("hTrendAvgetaeff_dot5_Nch_fbit0_type2",  "Eta efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgetaeff_dot8_Nch_fbit0_type2 = new TH1F("hTrendAvgetaeff_dot8_Nch_fbit0_type2", "Eta efficiency Vs Runs", nRuns, 0, nRuns);
    
    Float_t AvgpTeff_1_2_Pion_fbit0_type2 = 0.,AvgpTeff_2_6_Pion_fbit0_type2 = 0.,AvgpTeff_6_8_Pion_fbit0_type2 = 0., AvgpTeff_8_10_Pion_fbit0_type2=0.;
    Float_t Avgetaeff_dot5_Pion_fbit0_type2 = 0.,Avgetaeff_dot8_Pion_fbit0_type2 = 0;
    tTrend->SetBranchAddress("AvgpTeff_1_2_Pion_fbit0_type2",  &AvgpTeff_1_2_Pion_fbit0_type2);
    tTrend->SetBranchAddress("AvgpTeff_2_6_Pion_fbit0_type2",  &AvgpTeff_2_6_Pion_fbit0_type2);
    tTrend->SetBranchAddress("AvgpTeff_6_8_Pion_fbit0_type2",  &AvgpTeff_6_8_Pion_fbit0_type2);
    tTrend->SetBranchAddress("AvgpTeff_8_10_Pion_fbit0_type2", &AvgpTeff_8_10_Pion_fbit0_type2);
    tTrend->SetBranchAddress("Avgetaeff_dot5_Pion_fbit0_type2",  &Avgetaeff_dot5_Pion_fbit0_type2);
    tTrend->SetBranchAddress("Avgetaeff_dot8_Pion_fbit0_type2", &Avgetaeff_dot8_Pion_fbit0_type2);
    TH1F* hTrendAvgPteff_1_2_Pion_fbit0_type2  = new TH1F("hTrendAvgPteff_1_2_Pion_fbit0_type2",  "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgPteff_2_6_Pion_fbit0_type2  = new TH1F("hTrendAvgPteff_2_6_Pion_fbit0_type2",  "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgPteff_6_8_Pion_fbit0_type2  = new TH1F("hTrendAvgPteff_6_8_Pion_fbit0_type2",  "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgPteff_8_10_Pion_fbit0_type2 = new TH1F("hTrendAvgPteff_8_10_Pion_fbit0_type2", "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgetaeff_dot5_Pion_fbit0_type2 = new TH1F("hTrendAvgetaeff_dot5_Pion_fbit0_type2",  "Eta efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgetaeff_dot8_Pion_fbit0_type2 = new TH1F("hTrendAvgetaeff_dot8_Pion_fbit0_type2", "Eta efficiency Vs Runs", nRuns, 0, nRuns);
    
    Float_t AvgpTeff_1_2_Kaon_fbit0_type2 = 0.,AvgpTeff_2_6_Kaon_fbit0_type2 = 0.,AvgpTeff_6_8_Kaon_fbit0_type2 = 0., AvgpTeff_8_10_Kaon_fbit0_type2=0.;
    Float_t Avgetaeff_dot5_Kaon_fbit0_type2 = 0.,Avgetaeff_dot8_Kaon_fbit0_type2 = 0;
    tTrend->SetBranchAddress("AvgpTeff_1_2_Kaon_fbit0_type2",  &AvgpTeff_1_2_Kaon_fbit0_type2);
    tTrend->SetBranchAddress("AvgpTeff_2_6_Kaon_fbit0_type2",  &AvgpTeff_2_6_Kaon_fbit0_type2);
    tTrend->SetBranchAddress("AvgpTeff_6_8_Kaon_fbit0_type2",  &AvgpTeff_6_8_Kaon_fbit0_type2);
    tTrend->SetBranchAddress("AvgpTeff_8_10_Kaon_fbit0_type2", &AvgpTeff_8_10_Kaon_fbit0_type2);
    tTrend->SetBranchAddress("Avgetaeff_dot5_Kaon_fbit0_type2",  &Avgetaeff_dot5_Kaon_fbit0_type2);
    tTrend->SetBranchAddress("Avgetaeff_dot8_Kaon_fbit0_type2", &Avgetaeff_dot8_Kaon_fbit0_type2);
    TH1F* hTrendAvgPteff_1_2_Kaon_fbit0_type2  = new TH1F("hTrendAvgPteff_1_2_Kaon_fbit0_type2",  "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgPteff_2_6_Kaon_fbit0_type2  = new TH1F("hTrendAvgPteff_2_6_Kaon_fbit0_type2",  "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgPteff_6_8_Kaon_fbit0_type2  = new TH1F("hTrendAvgPteff_6_8_Kaon_fbit0_type2",  "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgPteff_8_10_Kaon_fbit0_type2 = new TH1F("hTrendAvgPteff_8_10_Kaon_fbit0_type2", "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgetaeff_dot5_Kaon_fbit0_type2 = new TH1F("hTrendAvgetaeff_dot5_Kaon_fbit0_type2",  "Eta efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgetaeff_dot8_Kaon_fbit0_type2 = new TH1F("hTrendAvgetaeff_dot8_Kaon_fbit0_type2", "Eta efficiency Vs Runs", nRuns, 0, nRuns);
    
    Float_t AvgpTeff_1_2_Elec_fbit0_type2 = 0.,AvgpTeff_2_6_Elec_fbit0_type2 = 0.,AvgpTeff_6_8_Elec_fbit0_type2 = 0., AvgpTeff_8_10_Elec_fbit0_type2=0.;
    Float_t Avgetaeff_dot5_Elec_fbit0_type2 = 0.,Avgetaeff_dot8_Elec_fbit0_type2 = 0;
    tTrend->SetBranchAddress("AvgpTeff_1_2_Elec_fbit0_type2",  &AvgpTeff_1_2_Elec_fbit0_type2);
    tTrend->SetBranchAddress("AvgpTeff_2_6_Elec_fbit0_type2",  &AvgpTeff_2_6_Elec_fbit0_type2);
    tTrend->SetBranchAddress("AvgpTeff_6_8_Elec_fbit0_type2",  &AvgpTeff_6_8_Elec_fbit0_type2);
    tTrend->SetBranchAddress("AvgpTeff_8_10_Elec_fbit0_type2", &AvgpTeff_8_10_Elec_fbit0_type2);
    tTrend->SetBranchAddress("Avgetaeff_dot5_Elec_fbit0_type2",  &Avgetaeff_dot5_Elec_fbit0_type2);
    tTrend->SetBranchAddress("Avgetaeff_dot8_Elec_fbit0_type2", &Avgetaeff_dot8_Elec_fbit0_type2);
    TH1F* hTrendAvgPteff_1_2_Elec_fbit0_type2  = new TH1F("hTrendAvgPteff_1_2_Elec_fbit0_type2",  "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgPteff_2_6_Elec_fbit0_type2  = new TH1F("hTrendAvgPteff_2_6_Elec_fbit0_type2",  "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgPteff_6_8_Elec_fbit0_type2  = new TH1F("hTrendAvgPteff_6_8_Elec_fbit0_type2",  "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgPteff_8_10_Elec_fbit0_type2 = new TH1F("hTrendAvgPteff_8_10_Elec_fbit0_type2", "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgetaeff_dot5_Elec_fbit0_type2 = new TH1F("hTrendAvgetaeff_dot5_Elec_fbit0_type2",  "Eta efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgetaeff_dot8_Elec_fbit0_type2 = new TH1F("hTrendAvgetaeff_dot8_Elec_fbit0_type2", "Eta efficiency Vs Runs", nRuns, 0, nRuns);
    
    
    Float_t AvgpTeff_1_2_Prot_fbit0_type2 = 0.,AvgpTeff_2_6_Prot_fbit0_type2 = 0.,AvgpTeff_6_8_Prot_fbit0_type2 = 0., AvgpTeff_8_10_Prot_fbit0_type2=0.;
    Float_t Avgetaeff_dot5_Prot_fbit0_type2 = 0.,Avgetaeff_dot8_Prot_fbit0_type2 = 0;
    tTrend->SetBranchAddress("AvgpTeff_1_2_Prot_fbit0_type2",  &AvgpTeff_1_2_Prot_fbit0_type2);
    tTrend->SetBranchAddress("AvgpTeff_2_6_Prot_fbit0_type2",  &AvgpTeff_2_6_Prot_fbit0_type2);
    tTrend->SetBranchAddress("AvgpTeff_6_8_Prot_fbit0_type2",  &AvgpTeff_6_8_Prot_fbit0_type2);
    tTrend->SetBranchAddress("AvgpTeff_8_10_Prot_fbit0_type2", &AvgpTeff_8_10_Prot_fbit0_type2);
    tTrend->SetBranchAddress("Avgetaeff_dot5_Prot_fbit0_type2",  &Avgetaeff_dot5_Prot_fbit0_type2);
    tTrend->SetBranchAddress("Avgetaeff_dot8_Prot_fbit0_type2", &Avgetaeff_dot8_Prot_fbit0_type2);
    TH1F* hTrendAvgPteff_1_2_Prot_fbit0_type2  = new TH1F("hTrendAvgPteff_1_2_Prot_fbit0_type2",  "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgPteff_2_6_Prot_fbit0_type2  = new TH1F("hTrendAvgPteff_2_6_Prot_fbit0_type2",  "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgPteff_6_8_Prot_fbit0_type2  = new TH1F("hTrendAvgPteff_6_8_Prot_fbit0_type2",  "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgPteff_8_10_Prot_fbit0_type2 = new TH1F("hTrendAvgPteff_8_10_Prot_fbit0_type2", "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgetaeff_dot5_Prot_fbit0_type2 = new TH1F("hTrendAvgetaeff_dot5_Prot_fbit0_type2",  "Eta efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgetaeff_dot8_Prot_fbit0_type2 = new TH1F("hTrendAvgetaeff_dot8_Prot_fbit0_type2", "Eta efficiency Vs Runs", nRuns, 0, nRuns);
    
    
    Float_t AvgpTeff_1_2_Nch_fbit4_type2 = 0.,AvgpTeff_2_6_Nch_fbit4_type2 = 0.,AvgpTeff_6_8_Nch_fbit4_type2 = 0., AvgpTeff_8_10_Nch_fbit4_type2=0.;
    Float_t Avgetaeff_dot5_Nch_fbit4_type2 = 0.,Avgetaeff_dot8_Nch_fbit4_type2 = 0;
    tTrend->SetBranchAddress("AvgpTeff_1_2_Nch_fbit4_type2",  &AvgpTeff_1_2_Nch_fbit4_type2);
    tTrend->SetBranchAddress("AvgpTeff_2_6_Nch_fbit4_type2",  &AvgpTeff_2_6_Nch_fbit4_type2);
    tTrend->SetBranchAddress("AvgpTeff_6_8_Nch_fbit4_type2",  &AvgpTeff_6_8_Nch_fbit4_type2);
    tTrend->SetBranchAddress("AvgpTeff_8_10_Nch_fbit4_type2", &AvgpTeff_8_10_Nch_fbit4_type2);
    tTrend->SetBranchAddress("Avgetaeff_dot5_Nch_fbit4_type2",  &Avgetaeff_dot5_Nch_fbit4_type2);
    tTrend->SetBranchAddress("Avgetaeff_dot8_Nch_fbit4_type2", &Avgetaeff_dot8_Nch_fbit4_type2);
    TH1F* hTrendAvgPteff_1_2_Nch_fbit4_type2  = new TH1F("hTrendAvgPteff_1_2_Nch_fbit4_type2",  "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgPteff_2_6_Nch_fbit4_type2  = new TH1F("hTrendAvgPteff_2_6_Nch_fbit4_type2",  "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgPteff_6_8_Nch_fbit4_type2  = new TH1F("hTrendAvgPteff_6_8_Nch_fbit4_type2",  "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgPteff_8_10_Nch_fbit4_type2 = new TH1F("hTrendAvgPteff_8_10_Nch_fbit4_type2", "p_{T} efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgetaeff_dot5_Nch_fbit4_type2 = new TH1F("hTrendAvgetaeff_dot5_Nch_fbit4_type2",  "Eta efficiency Vs Runs", nRuns, 0, nRuns);
    TH1F* hTrendAvgetaeff_dot8_Nch_fbit4_type2 = new TH1F("hTrendAvgetaeff_dot8_Nch_fbit4_type2", "Eta efficiency Vs Runs", nRuns, 0, nRuns);
    

    
    for (Int_t r=0;r<nRuns;r++){
        
        tTrend->GetEntry(r);
        printf("... fetching run = %i \n",run);
        
        //1. Pt eff hist Nch + Filtebit0
        hTrendAvgPteff_1_2_Nch_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_1_2_Nch_fbit0_type2->SetBinContent(r+1, AvgpTeff_1_2_Nch_fbit0_type2);
        hTrendAvgPteff_2_6_Nch_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_2_6_Nch_fbit0_type2->SetBinContent(r+1, AvgpTeff_2_6_Nch_fbit0_type2);
        hTrendAvgPteff_6_8_Nch_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_6_8_Nch_fbit0_type2->SetBinContent(r+1, AvgpTeff_6_8_Nch_fbit0_type2);
        hTrendAvgPteff_8_10_Nch_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_8_10_Nch_fbit0_type2->SetBinContent(r+1, AvgpTeff_8_10_Nch_fbit0_type2);
        hTrendAvgetaeff_dot5_Nch_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgetaeff_dot5_Nch_fbit0_type2->SetBinContent(r+1, Avgetaeff_dot5_Nch_fbit0_type2);
        hTrendAvgetaeff_dot8_Nch_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgetaeff_dot8_Nch_fbit0_type2->SetBinContent(r+1, Avgetaeff_dot8_Nch_fbit0_type2);
        
        
        
        //2. Pt eff hist Pions + Filtebit0
        hTrendAvgPteff_1_2_Pion_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_1_2_Pion_fbit0_type2->SetBinContent(r+1, AvgpTeff_1_2_Pion_fbit0_type2);
        hTrendAvgPteff_2_6_Pion_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_2_6_Pion_fbit0_type2->SetBinContent(r+1, AvgpTeff_2_6_Pion_fbit0_type2);
        hTrendAvgPteff_6_8_Pion_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_6_8_Pion_fbit0_type2->SetBinContent(r+1, AvgpTeff_6_8_Pion_fbit0_type2);
        hTrendAvgPteff_8_10_Pion_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_8_10_Pion_fbit0_type2->SetBinContent(r+1, AvgpTeff_8_10_Pion_fbit0_type2);
        hTrendAvgetaeff_dot5_Pion_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgetaeff_dot5_Pion_fbit0_type2->SetBinContent(r+1, Avgetaeff_dot5_Pion_fbit0_type2);
        hTrendAvgetaeff_dot8_Pion_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgetaeff_dot8_Pion_fbit0_type2->SetBinContent(r+1, Avgetaeff_dot8_Pion_fbit0_type2);
        
        
        //3. Pt eff hist Kaons + Filtebit0
        hTrendAvgPteff_1_2_Kaon_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_1_2_Kaon_fbit0_type2->SetBinContent(r+1, AvgpTeff_1_2_Kaon_fbit0_type2);
        hTrendAvgPteff_2_6_Kaon_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_2_6_Kaon_fbit0_type2->SetBinContent(r+1, AvgpTeff_2_6_Kaon_fbit0_type2);
        hTrendAvgPteff_6_8_Kaon_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_6_8_Kaon_fbit0_type2->SetBinContent(r+1, AvgpTeff_6_8_Kaon_fbit0_type2);
        hTrendAvgPteff_8_10_Kaon_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_8_10_Kaon_fbit0_type2->SetBinContent(r+1, AvgpTeff_8_10_Kaon_fbit0_type2);
        hTrendAvgetaeff_dot5_Kaon_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgetaeff_dot5_Kaon_fbit0_type2->SetBinContent(r+1, Avgetaeff_dot5_Kaon_fbit0_type2);
        hTrendAvgetaeff_dot8_Kaon_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgetaeff_dot8_Kaon_fbit0_type2->SetBinContent(r+1, Avgetaeff_dot8_Kaon_fbit0_type2);
        
        
        
        //4. Pt eff hist Electros + Filtebit0
        hTrendAvgPteff_1_2_Elec_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_1_2_Elec_fbit0_type2->SetBinContent(r+1, AvgpTeff_1_2_Elec_fbit0_type2);
        hTrendAvgPteff_2_6_Elec_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_2_6_Elec_fbit0_type2->SetBinContent(r+1, AvgpTeff_2_6_Elec_fbit0_type2);
        hTrendAvgPteff_6_8_Elec_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_6_8_Elec_fbit0_type2->SetBinContent(r+1, AvgpTeff_6_8_Elec_fbit0_type2);
        hTrendAvgPteff_8_10_Elec_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_8_10_Elec_fbit0_type2->SetBinContent(r+1, AvgpTeff_8_10_Elec_fbit0_type2);
        hTrendAvgetaeff_dot5_Elec_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgetaeff_dot5_Elec_fbit0_type2->SetBinContent(r+1, Avgetaeff_dot5_Elec_fbit0_type2);
        hTrendAvgetaeff_dot8_Elec_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgetaeff_dot8_Elec_fbit0_type2->SetBinContent(r+1, Avgetaeff_dot8_Elec_fbit0_type2);
        
        
        
        //5. Pt eff hist Proton + Filtebit0
        hTrendAvgPteff_1_2_Prot_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_1_2_Prot_fbit0_type2->SetBinContent(r+1, AvgpTeff_1_2_Prot_fbit0_type2);
        hTrendAvgPteff_2_6_Prot_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_2_6_Prot_fbit0_type2->SetBinContent(r+1, AvgpTeff_2_6_Prot_fbit0_type2);
        hTrendAvgPteff_6_8_Prot_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_6_8_Prot_fbit0_type2->SetBinContent(r+1, AvgpTeff_6_8_Prot_fbit0_type2);
        hTrendAvgPteff_8_10_Prot_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_8_10_Prot_fbit0_type2->SetBinContent(r+1, AvgpTeff_8_10_Prot_fbit0_type2);
        hTrendAvgetaeff_dot5_Prot_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgetaeff_dot5_Prot_fbit0_type2->SetBinContent(r+1, Avgetaeff_dot5_Prot_fbit0_type2);
        hTrendAvgetaeff_dot8_Prot_fbit0_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgetaeff_dot8_Prot_fbit0_type2->SetBinContent(r+1, Avgetaeff_dot8_Prot_fbit0_type2);
        
        
        //6. Pt eff hist Nch + Filtebit4
        hTrendAvgPteff_1_2_Nch_fbit4_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_1_2_Nch_fbit4_type2->SetBinContent(r+1, AvgpTeff_1_2_Nch_fbit4_type2);
        hTrendAvgPteff_2_6_Nch_fbit4_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_2_6_Nch_fbit4_type2->SetBinContent(r+1, AvgpTeff_2_6_Nch_fbit4_type2);
        hTrendAvgPteff_6_8_Nch_fbit4_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_6_8_Nch_fbit4_type2->SetBinContent(r+1, AvgpTeff_6_8_Nch_fbit4_type2);
        hTrendAvgPteff_8_10_Nch_fbit4_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgPteff_8_10_Nch_fbit4_type2->SetBinContent(r+1, AvgpTeff_8_10_Nch_fbit4_type2);
        hTrendAvgetaeff_dot5_Nch_fbit4_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgetaeff_dot5_Nch_fbit4_type2->SetBinContent(r+1, Avgetaeff_dot5_Nch_fbit4_type2);
        hTrendAvgetaeff_dot8_Nch_fbit4_type2->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
        hTrendAvgetaeff_dot8_Nch_fbit4_type2->SetBinContent(r+1, Avgetaeff_dot8_Nch_fbit4_type2);
        
        
    }
    
    printf(".....DONE !!");
    gSystem->Exec("mkdir -p Trending/PDFs");
    gSystem->Exec("mkdir -p Trending/root");
    
    //Pt Nch+Fitbit0 eff hist
    TFile* fPtEffTredingVsRun_Nch_fbit0_type2 = new TFile("Trending/root/PlotTredingPtVsRun_Nch_fbit0_type2.root","recreate");
    TCanvas* cTrendingNchPt = new TCanvas("cTrendingNch","p_{T} Eff Vs Runs <Nch_fbit0>",900,600);
    gPad->SetMargin(0.08,0.02,0.08,0.08);
    
    hTrendAvgPteff_1_2_Nch_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_1_2_Nch_fbit0_type2);
    StyleLatex("hTrendAvgPteff_1_2_Nch_fbit0_type2",1.0, 2.0, 0, "Nch", 0.6);
    hTrendAvgPteff_1_2_Nch_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Nch_fbit0_type2.pdf(");
    
    hTrendAvgPteff_2_6_Nch_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_2_6_Nch_fbit0_type2);
    StyleLatex("hTrendAvgPteff_2_6_Nch_fbit0_type2",2.0, 6.0, 0, "Nch", 0.6);
    hTrendAvgPteff_2_6_Nch_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Nch_fbit0_type2.pdf");
    
    hTrendAvgPteff_6_8_Nch_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_6_8_Nch_fbit0_type2);
    StyleLatex("hTrendAvgPteff_6_8_Nch_fbit0_type2",6.0, 8.0, 0, "Nch", 0.6);
    hTrendAvgPteff_6_8_Nch_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Nch_fbit0_type2.pdf");
    
    hTrendAvgPteff_8_10_Nch_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_8_10_Nch_fbit0_type2);
    StyleLatex("hTrendAvgPteff_8_10_Nch_fbit0_type2",8.0, 10.0, 0, "Nch", 0.6);
    hTrendAvgPteff_8_10_Nch_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Nch_fbit0_type2.pdf)");
    fPtEffTredingVsRun_Nch_fbit0_type2->Close();
    
    TFile* fEtaEffTredingVsRun_Nch_fbit0_type2 = new TFile("Trending/root/PlotTredingEtaVsRun_Nch_fbit0_type2.root","recreate");
    TCanvas* cTrendingNchEta = new TCanvas("cTrendingNchEta","Eta Eff Vs Runs <Nch_fbit0>",900,600);
    gPad->SetMargin(0.08,0.02,0.08,0.08);
    
    hTrendAvgetaeff_dot5_Nch_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgetaeff_dot5_Nch_fbit0_type2);
    StyleLatex("hTrendAvgetaeff_dot5_Nch_fbit0_type2",-0.5, 0.5, 0, "Nch", 0.6);
    hTrendAvgetaeff_dot5_Nch_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingEtaVsRun_Nch_fbit0_type2.pdf(");
    
    hTrendAvgetaeff_dot8_Nch_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgetaeff_dot8_Nch_fbit0_type2);
    StyleLatex("hTrendAvgetaeff_dot8_Nch_fbit0_type2",-0.8, 0.8, 0, "Nch", 0.6);
    hTrendAvgetaeff_dot8_Nch_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingEtaVsRun_Nch_fbit0_type2.pdf)");
    fEtaEffTredingVsRun_Nch_fbit0_type2->Close();
    
    
    
    //Pt  Pion+Fitbit0 eff hist
    TFile* fPtEffTredingVsRun_Pion_fbit0_type2 = new TFile("Trending/root/PlotTredingPtVsRun_Pion_fbit0_type2.root","recreate");
    TCanvas* cTrendingPionPt = new TCanvas("cTrendingPionPt","p_{T} Eff Vs Runs <Pion_fbit0>",900,600);
    gPad->SetMargin(0.08,0.02,0.08,0.08);
    
    hTrendAvgPteff_1_2_Pion_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_1_2_Pion_fbit0_type2);
    StyleLatex("hTrendAvgPteff_1_2_Pion_fbit0_type2",1.0, 2.0, 0, "Pion", 0.6);
    hTrendAvgPteff_1_2_Pion_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Pion_fbit0_type2.pdf(");
    
    hTrendAvgPteff_2_6_Pion_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_2_6_Pion_fbit0_type2);
    StyleLatex("hTrendAvgPteff_2_6_Pion_fbit0_type2",2.0, 6.0, 0, "Pion", 0.6);
    hTrendAvgPteff_2_6_Pion_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Pion_fbit0_type2.pdf");
    
    hTrendAvgPteff_6_8_Pion_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_6_8_Pion_fbit0_type2);
    StyleLatex("hTrendAvgPteff_6_8_Pion_fbit0_type2",6.0, 8.0, 0, "Pion", 0.6);
    hTrendAvgPteff_6_8_Pion_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Pion_fbit0_type2.pdf");
    
    hTrendAvgPteff_8_10_Pion_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_8_10_Pion_fbit0_type2);
    StyleLatex("hTrendAvgPteff_8_10_Pion_fbit0_type2",8.0, 10.0, 0, "Pion", 0.6);
    hTrendAvgPteff_8_10_Pion_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Pion_fbit0_type2.pdf)");
    fPtEffTredingVsRun_Pion_fbit0_type2->Close();
    
    
    TFile* fEtaEffTredingVsRun_Pion_fbit0_type2 = new TFile("Trending/root/PlotTredingEtaVsRun_Pion_fbit0_type2.root","recreate");
    TCanvas* cTrendingPionEta = new TCanvas("cTrendingPionEta","Eta Eff Vs Runs <Pion_fbit0>",900,600);
    gPad->SetMargin(0.08,0.02,0.08,0.08);
    
    hTrendAvgetaeff_dot5_Pion_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgetaeff_dot5_Pion_fbit0_type2);
    StyleLatex("hTrendAvgetaeff_dot5_Pion_fbit0_type2",-0.5, 0.5, 0, "Pion", 0.6);
    hTrendAvgetaeff_dot5_Pion_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingEtaVsRun_Pion_fbit0_type2.pdf(");
    
    hTrendAvgetaeff_dot8_Pion_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgetaeff_dot8_Pion_fbit0_type2);
    StyleLatex("hTrendAvgetaeff_dot8_Pion_fbit0_type2",-0.8, 0.8, 0, "Pion", 0.6);
    hTrendAvgetaeff_dot8_Pion_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingEtaVsRun_Pion_fbit0_type2.pdf)");
    fEtaEffTredingVsRun_Pion_fbit0_type2->Close();
    
    
    
    //Pt<>Eta  Kaon+Fitbit0 eff hist
    TFile* fPtEffTredingVsRun_Kaon_fbit0_type2 = new TFile("Trending/root/PlotTredingPtVsRun_Kaon_fbit0_type2.root","recreate");
    TCanvas* cTrendingKaonPt = new TCanvas("cTrendingKaonPt","p_{T} Eff Vs Runs <Kaon_fbit0>",900,600);
    gPad->SetMargin(0.08,0.02,0.08,0.08);
    
    hTrendAvgPteff_1_2_Kaon_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_1_2_Kaon_fbit0_type2);
    StyleLatex("hTrendAvgPteff_1_2_Kaon_fbit0_type2",1.0, 2.0, 0, "Kaon", 0.6);
    hTrendAvgPteff_1_2_Kaon_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Kaon_fbit0_type2.pdf(");
    
    hTrendAvgPteff_2_6_Kaon_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_2_6_Kaon_fbit0_type2);
    StyleLatex("hTrendAvgPteff_2_6_Kaon_fbit0_type2",2.0, 6.0, 0, "Kaon", 0.6);
    hTrendAvgPteff_2_6_Kaon_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Kaon_fbit0_type2.pdf");
    
    hTrendAvgPteff_6_8_Kaon_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_6_8_Kaon_fbit0_type2);
    StyleLatex("hTrendAvgPteff_6_8_Kaon_fbit0_type2",6.0, 8.0, 0, "Kaon", 0.6);
    hTrendAvgPteff_6_8_Kaon_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Kaon_fbit0_type2.pdf");
    
    hTrendAvgPteff_8_10_Kaon_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_8_10_Kaon_fbit0_type2);
    StyleLatex("hTrendAvgPteff_8_10_Kaon_fbit0_type2",8.0, 10.0, 0, "Kaon", 0.6);
    hTrendAvgPteff_8_10_Kaon_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Kaon_fbit0_type2.pdf)");
    fPtEffTredingVsRun_Kaon_fbit0_type2->Close();
    
    TFile* fEtaEffTredingVsRun_Kaon_fbit0_type2 = new TFile("Trending/root/PlotTredingEtaVsRun_Kaon_fbit0_type2.root","recreate");
    TCanvas* cTrendingKaonEta = new TCanvas("cTrendingKaonEta","Eta Eff Vs Runs <Kaon_fbit0>",900,600);
    gPad->SetMargin(0.08,0.02,0.08,0.08);
    
    hTrendAvgetaeff_dot5_Kaon_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgetaeff_dot5_Kaon_fbit0_type2);
    StyleLatex("hTrendAvgetaeff_dot5_Kaon_fbit0_type2",-0.5, 0.5, 0, "Kaon", 0.6);
    hTrendAvgetaeff_dot5_Kaon_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingEtaVsRun_Kaon_fbit0_type2.pdf(");
    
    hTrendAvgetaeff_dot8_Kaon_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgetaeff_dot8_Kaon_fbit0_type2);
    StyleLatex("hTrendAvgetaeff_dot8_Kaon_fbit0_type2",-0.8, 0.8, 0, "Kaon", 0.6);
    hTrendAvgetaeff_dot8_Kaon_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingEtaVsRun_Kaon_fbit0_type2.pdf)");
    fEtaEffTredingVsRun_Kaon_fbit0_type2->Close();
    
    //Pt<>Eta  Elec+Fitbit0 eff hist
    TFile* fPtEffTredingVsRun_Elec_fbit0_type2 = new TFile("Trending/root/PlotTredingPtVsRun_Elec_fbit0_type2.root","recreate");
    TCanvas* cTrendingElecPt = new TCanvas("cTrendingElecPt","p_{T} Eff Vs Runs <Elec_fbit0>",900,600);
    gPad->SetMargin(0.08,0.02,0.08,0.08);
    
    hTrendAvgPteff_1_2_Elec_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_1_2_Elec_fbit0_type2);
    StyleLatex("hTrendAvgPteff_1_2_Elec_fbit0_type2",1.0, 2.0, 0, "Elec", 0.6);
    hTrendAvgPteff_1_2_Elec_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Elec_fbit0_type2.pdf(");
    
    hTrendAvgPteff_2_6_Elec_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_2_6_Elec_fbit0_type2);
    StyleLatex("hTrendAvgPteff_2_6_Elec_fbit0_type2",2.0, 6.0, 0, "Elec", 0.6);
    hTrendAvgPteff_2_6_Elec_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Elec_fbit0_type2.pdf");
    
    hTrendAvgPteff_6_8_Elec_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_6_8_Elec_fbit0_type2);
    StyleLatex("hTrendAvgPteff_6_8_Elec_fbit0_type2",6.0, 8.0, 0, "Elec", 0.6);
    hTrendAvgPteff_6_8_Elec_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Elec_fbit0_type2.pdf");
    
    hTrendAvgPteff_8_10_Elec_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_8_10_Elec_fbit0_type2);
    StyleLatex("hTrendAvgPteff_8_10_Elec_fbit0_type2",8.0, 10.0, 0, "Elec", 0.6);
    hTrendAvgPteff_8_10_Elec_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Elec_fbit0_type2.pdf)");
    fPtEffTredingVsRun_Elec_fbit0_type2->Close();
    
    
    TFile* fEtaEffTredingVsRun_Elec_fbit0_type2 = new TFile("Trending/root/PlotTredingEtaVsRun_Elec_fbit0_type2.root","recreate");
    TCanvas* cTrendingElecEta = new TCanvas("cTrendingElecEta","Eta Eff Vs Runs <Elec_fbit0>",900,600);
    gPad->SetMargin(0.08,0.02,0.08,0.08);
    
    hTrendAvgetaeff_dot5_Elec_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgetaeff_dot5_Elec_fbit0_type2);
    StyleLatex("hTrendAvgetaeff_dot5_Elec_fbit0_type2",-0.5, 0.5, 0, "Elec", 0.6);
    hTrendAvgetaeff_dot5_Elec_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingEtaVsRun_Elec_fbit0_type2.pdf(");
    
    hTrendAvgetaeff_dot8_Elec_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgetaeff_dot8_Elec_fbit0_type2);
    StyleLatex("hTrendAvgetaeff_dot8_Elec_fbit0_type2",-0.8, 0.8, 0, "Elec", 0.6);
    hTrendAvgetaeff_dot8_Elec_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingEtaVsRun_Elec_fbit0_type2.pdf)");
    fEtaEffTredingVsRun_Elec_fbit0_type2->Close();
    
    
    
    //Pt<>Eta  Proton+Fitbit0 eff hist
    TFile* fPtEffTredingVsRun_Prot_fbit0_type2 = new TFile("Trending/root/PlotTredingPtVsRun_Prot_fbit0_type2.root","recreate");
    TCanvas* cTrendingProtPt = new TCanvas("cTrendingProtPt","p_{T} Eff Vs Runs <Elec_fbit0>",900,600);
    gPad->SetMargin(0.08,0.02,0.08,0.08);
    
    hTrendAvgPteff_1_2_Prot_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_1_2_Prot_fbit0_type2);
    StyleLatex("hTrendAvgPteff_1_2_Prot_fbit0_type2",1.0, 2.0, 0, "Proton", 0.6);
    hTrendAvgPteff_1_2_Prot_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Prot_fbit0_type2.pdf(");
    
    hTrendAvgPteff_2_6_Prot_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_2_6_Prot_fbit0_type2);
    StyleLatex("hTrendAvgPteff_2_6_Prot_fbit0_type2",2.0, 6.0, 0, "Proton", 0.6);
    hTrendAvgPteff_2_6_Prot_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Prot_fbit0_type2.pdf");
    
    hTrendAvgPteff_6_8_Prot_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_6_8_Prot_fbit0_type2);
    StyleLatex("hTrendAvgPteff_6_8_Prot_fbit0_type2",6.0, 8.0, 0, "Proton", 0.6);
    hTrendAvgPteff_6_8_Prot_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Prot_fbit0_type2.pdf");
    
    hTrendAvgPteff_8_10_Prot_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_8_10_Prot_fbit0_type2);
    StyleLatex("hTrendAvgPteff_8_10_Prot_fbit0_type2",8.0, 10.0, 0, "Proton", 0.6);
    hTrendAvgPteff_8_10_Prot_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Prot_fbit0_type2.pdf)");
    fPtEffTredingVsRun_Prot_fbit0_type2->Close();
    
    TFile* fEtaEffTredingVsRun_Prot_fbit0_type2 = new TFile("Trending/root/PlotTredingEtaVsRun_Prot_fbit0_type2.root","recreate");
    TCanvas* cTrendingProtcEta = new TCanvas("cTrendingProtEta","Eta Eff Vs Runs <Elec_fbit0>",900,600);
    gPad->SetMargin(0.08,0.02,0.08,0.08);
    
    hTrendAvgetaeff_dot5_Prot_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgetaeff_dot5_Prot_fbit0_type2);
    StyleLatex("hTrendAvgetaeff_dot5_Prot_fbit0_type2",-0.5, 0.5, 0, "Proton", 0.6);
    hTrendAvgetaeff_dot5_Prot_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingEtaVsRun_Prot_fbit0_type2.pdf(");
    
    hTrendAvgetaeff_dot8_Prot_fbit0_type2->Draw("PH");
    StyleHist(hTrendAvgetaeff_dot8_Prot_fbit0_type2);
    StyleLatex("hTrendAvgetaeff_dot8_Prot_fbit0_type2",-0.8, 0.8, 0, "Proton", 0.6);
    hTrendAvgetaeff_dot8_Prot_fbit0_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingEtaVsRun_Prot_fbit0_type2.pdf)");
    fEtaEffTredingVsRun_Prot_fbit0_type2->Close();

    
    
    //Pt<>Eta  Nch+Fitbit4 eff hist
    TFile* fPtEffTredingVsRun_Nch_fbit4_type2 = new TFile("Trending/root/PlotTredingPtVsRun_Nch_fbit4_type2.root","recreate");
    TCanvas* cTrendingNchPtf4 = new TCanvas("cTrendingNchPtF4","p_{T} Eff Vs Runs <Nch_fbit4>",900,600);
    gPad->SetMargin(0.08,0.02,0.08,0.08);
    
    hTrendAvgPteff_1_2_Nch_fbit4_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_1_2_Nch_fbit4_type2);
    StyleLatex("hTrendAvgPteff_1_2_Nch_fbit4_type2",1.0, 2.0, 4, "Nch", 0.6);
    hTrendAvgPteff_1_2_Nch_fbit4_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Nch_fbit4_type2.pdf(");
    
    hTrendAvgPteff_2_6_Nch_fbit4_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_2_6_Nch_fbit4_type2);
    StyleLatex("hTrendAvgPteff_2_6_Nch_fbit4_type2",2.0, 6.0, 4, "Nch", 0.6);
    hTrendAvgPteff_2_6_Nch_fbit4_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Nch_fbit4_type2.pdf");
    
    hTrendAvgPteff_6_8_Nch_fbit4_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_6_8_Nch_fbit4_type2);
    StyleLatex("hTrendAvgPteff_6_8_Nch_fbit4_type2",6.0, 8.0, 4, "Nch", 0.6);
    hTrendAvgPteff_6_8_Nch_fbit4_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Nch_fbit4_type2.pdf");
    
    hTrendAvgPteff_8_10_Nch_fbit4_type2->Draw("PH");
    StyleHist(hTrendAvgPteff_8_10_Nch_fbit4_type2);
    StyleLatex("hTrendAvgPteff_8_10_Nch_fbit4_type2",8.0, 10.0, 4, "Nch", 0.6);
    hTrendAvgPteff_8_10_Nch_fbit4_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingPtVsRun_Nch_fbit4_type2.pdf)");
    fPtEffTredingVsRun_Nch_fbit4_type2->Close();
    
    
    TFile* fEtaEffTredingVsRun_Nch_fbit4_type2 = new TFile("Trending/root/PlotTredingEtaVsRun_Nch_fbit4_type2.root","recreate");
    TCanvas* cTrendingNchEtaf4 = new TCanvas("cTrendingNchEta4","Eta Eff Vs Runs <Nch_fbit4>",900,600);
    gPad->SetMargin(0.08,0.02,0.08,0.08);
    
    hTrendAvgetaeff_dot5_Nch_fbit4_type2->Draw("PH");
    StyleHist(hTrendAvgetaeff_dot5_Nch_fbit4_type2);
    StyleLatex("hTrendAvgetaeff_dot5_Nch_fbit4_type2",-0.5, 0.5, 4, "Nch", 0.6);
    hTrendAvgetaeff_dot5_Nch_fbit4_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingEtaVsRun_Nch_fbit4_type2.pdf(");
    
    hTrendAvgetaeff_dot8_Nch_fbit4_type2->Draw("PH");
    StyleHist(hTrendAvgetaeff_dot8_Nch_fbit4_type2);
    StyleLatex("hTrendAvgetaeff_dot8_Nch_fbit4_type2",-0.8, 0.8, 4, "Nch", 0.6);
    hTrendAvgetaeff_dot8_Nch_fbit4_type2->Write();
    gPad->Print("Trending/PDFs/PlotTredingEtaVsRun_Nch_fbit4_type2.pdf)");
    fEtaEffTredingVsRun_Nch_fbit4_type2->Close();
    
}


//________________| Setting Name for histogram
TH1F *StyleHist(TH1F *h)
{
    h->SetMinimum(0.20);
    h->SetMaximum(1.10);
    h->SetMarkerSize(1.5);
    h->SetMarkerColor(4);
    h->SetMarkerStyle(20);
    h->LabelsDeflate("Y");
    h->GetYaxis()->SetTitle("Efficiency");
    h->GetXaxis()->SetTitle("Run Number");
    h->LabelsOption("u", "X");
    //h->SetTitle("");
    
    return h;
}


//________________| TLatex
void StyleLatex(TString str="title", Double_t fmin, Double_t fmax, const Int_t fbit=0, const char *fPartype="", const Double_t xvalue =5){
    
    TLatex Tl(0.5,0.5,"u");
    Tl.SetTextFont(42);
    Tl.SetTextAlign(12);
    Tl.SetTextSize(0.035);
    Tl.SetTextColor(46);
    Tl.DrawLatex(xvalue,0.64, Form("Configurations"));
    Tl.DrawLatex(xvalue,0.62, "___________________");
    Tl.DrawLatex(xvalue,0.58, Form("1. Filterbit: %d", fbit));
    Tl.DrawLatex(xvalue,0.54, Form("2. Particle Type: %s", fPartype));
    
    if(str.Contains("PT") || str.Contains("PT") || str.Contains("Pt")){
        Tl.DrawLatex(xvalue,0.50, Form("3. %.1f #leq #eta #leq %.1f ", -0.9, 0.9));
        Tl.DrawLatex(xvalue,0.46, Form("4. %.1f #leq Zvtx #leq %.1f (cm)", -10., 10.));
        Tl.DrawLatex(xvalue,0.42, Form("5. %.1f #leq p_{T} #leq %.1f (GeV/c)", fmin, fmax));
    }
    else if(str.Contains("eta") || str.Contains("ETA") || str.Contains("Eta")){
        Tl.DrawLatex(xvalue,0.50, Form("3. %.1f #leq p_{T} #leq %.1f  (GeV/c)", 0.5, 24.));
        Tl.DrawLatex(xvalue,0.46, Form("4. %.1f #leq Zvtx #leq %.1f (cm)", -10., 10.));
        Tl.DrawLatex(xvalue,0.42, Form("5. %.1f #leq Eta #leq %.1f ", fmin, fmax));
    }
    else if(str.Contains("zvtx") || str.Contains("ZVTX") || str.Contains("Zvtx")){
        Tl.DrawLatex(xvalue,0.50, Form("3. %.1f #leq #eta #leq %.1f ", -0.9, 0.9));
        Tl.DrawLatex(xvalue,0.46, Form("4. %.1f #leq p_{T} #leq %.1f  (GeV/c)", 0., 24.));
        Tl.DrawLatex(xvalue,0.42, Form("5. %.1f #leq AvgZv #leq %.1f (cm)", fmin, fmax));
    }
    else {
        Tl.DrawLatex(xvalue,0.50, Form("3. %.1f #leq Zvtx #leq %.1f (cm)",  -10.,  10.));
        Tl.DrawLatex(xvalue,0.46, Form("4. %.1f #leq #eta #leq %.1f ", -0.9,  0.9));
        Tl.DrawLatex(xvalue,0.42, Form("5. %.1f #leq p_{T} #leq %.1f (GeV/c)",  0.0,  24.0));
    }    
}

