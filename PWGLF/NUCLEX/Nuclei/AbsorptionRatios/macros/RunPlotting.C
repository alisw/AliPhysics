#ifdef __CLING__
#include "Functions.C"
#include "RawParticleSpectra.C"
#include "Secondary.C"
#include "Summary.C"
#endif

void SaveRootFile(TString Name, TH1F *Histo){
    TFile *myfile = TFile::Open(Name.Data(),"RECREATE");
    Histo->Write();
}

void RunPlotting(){
    
    Bool_t doEffCorrection = kFALSE;
    Bool_t doPurCorrection = kTRUE;
    Bool_t doSecCorrection = kTRUE;
    
    TString FolderName[4] = {"ProtonTrackCuts","AntiProtonTrackCuts","DeuteronTrackCuts","AntiDeuteronTrackCuts"};
    TString FolderNameMC[4] = {"ProtonTrkCutsMC","AntiProtonTrkCutsMC","DeuteronTrkCutsMC","AntiDeuteronTrkCutsMC"};
    TString Particle[4] = {"Proton","Antiproton","Deuteron","Antideuteron"};
    TString FileNameStandard;
    TString FileNameEnriched;
    TString FileNameData;
    FileNameStandard = "/Users/lucascordova/Documents/Dokumente/Studium/TUM/werkstudent/AntiParticleAnalysis/190326_LightN/MC_f2b_cent_q/LightNAnalysis/AnalysisResults.root";
    FileNameEnriched = "/Users/lucascordova/Documents/Dokumente/Studium/TUM/werkstudent/AntiParticleAnalysis/190114_LightN/MC_d10_cent/LightNAnalysis/AnalysisResults.root";
    FileNameData = "/Users/lucascordova/Documents/Dokumente/Studium/TUM/werkstudent/AntiParticleAnalysis/190326_LightN/Data_qwSDD/LightNAnalysis/AnalysisResults.root";
    TFile* FileMCStandard=TFile::Open(FileNameStandard.Data(),"READ"); //number 1
    TFile* FileMCEnriched=TFile::Open(FileNameEnriched.Data(),"READ"); //number 2
    TFile* FileData=TFile::Open(FileNameData.Data(),"READ"); // number 3
    
    //Keys to get Data from directories
    TDirectoryFile *ParticleDirectory[4];
    TList *trackCutsData[4];
    TList *tmpListKeys[4];
    TString ListNameResults[4];
    TList *ListResults[4];
    
    //Keys to get MC Data from direcories
    TDirectoryFile *ParticleDirectoryMC[4];
    TList *trackCutsMC[4];
    TList *tmpListKeysMC[4];
    TString ListNameResultsMC[4];
    TList *ListResultsMC[4];
    TList *ListMCDCA[4];
    
    //momentum spectra for different Steps
    TH1F *MomentumData[4];
    TH1F *RawSpectrum[4];
    TH1F *PrimarySpectrum[4];
    
    //Histos for Signal extraction
    TH2F *Mass2sqData[4];
    TH1F *QAMassFit[4];
    
    //Histos for DCAFits
    TH2F *DCAData[4];
    TH2F *PrimaryMC[4];
    TH2F *SecondaryMC[4];
    TH2F *MaterialMC[4];
    TH1F *fractions_primary[4];
    
    
    TH1F* Ratio_primary[2];
    float binsx_p[28] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0};
    Ratio_primary[0] = new TH1F(Form("Ratios_primary_%s",Particle[0].Data()), Form("Ratios_primary_%s",Particle[0].Data()),45, 0, 4.5);
    Ratio_primary[1] = new TH1F(Form("Ratios_primary_%s",Particle[2].Data()), Form("Ratios_primary_%s",Particle[2].Data()), 27 , binsx_p);
    
    for(Int_t i=0; i<4; i++){
        
        //Getting Data
        ParticleDirectory[i] = (TDirectoryFile*)(FileData->FindObjectAny(FolderName[i].Data()));
        if(!ParticleDirectory[i]){
            std::cout << "Data Directory Missing" << std::endl;
        }
        tmpListKeys[i] = ParticleDirectory[i]->GetListOfKeys();
        ListNameResults[i] = tmpListKeys[i]->At(0)->GetName();
        ParticleDirectory[i]->GetObject(ListNameResults[i], ListResults[i]);
        trackCutsData[i] = (TList*)(ListResults[i]->FindObject("after"));
        
        //Getting MC (Satandad and Enriched)
        if(i<2){
            ParticleDirectoryMC[i] = (TDirectoryFile*)(FileMCStandard->FindObjectAny(FolderNameMC[i].Data()));
            if(!ParticleDirectoryMC[i]){
                std::cout << "MCStandard Directory Missing" << std::endl;
            }
        }else{
            ParticleDirectoryMC[i] = (TDirectoryFile*)(FileMCEnriched->FindObjectAny(FolderNameMC[i].Data()));
            if(!ParticleDirectoryMC[i]){
                std::cout << "MCEnriched Directory Missing" << std::endl;
            }
        }
        tmpListKeysMC[i] = ParticleDirectoryMC[i]->GetListOfKeys();
        ListNameResultsMC[i] = tmpListKeysMC[i]->At(0)->GetName();
        ParticleDirectoryMC[i]->GetObject(ListNameResultsMC[i], ListResultsMC[i]);
        
        //Spectra for different steps
        MomentumData[i] = (TH1F*)trackCutsData[i]->FindObject("pDist_after");
        if(i<2)RawSpectrum[i] = new TH1F(Form("hRawSpectrum %s",Particle[i].Data()), Form("hRawSpectrum %s",Particle[i].Data()),45,0,4.5);
        if(i>1)RawSpectrum[i] = new TH1F(Form("hRawSpectrum %s",Particle[i].Data()), Form("hRawSpectrum %s",Particle[i].Data()), 27 , binsx_p);
        if(i<2)PrimarySpectrum[i] = new TH1F(Form("PrimarySpectrum %s",Particle[i].Data()), Form("PrimarySpectrum %s",Particle[i].Data()),45,0,4.5);
        if(i>1)PrimarySpectrum[i] = new TH1F(Form("PrimarySpectrum %s",Particle[i].Data()), Form("PrimarySpectrum %s",Particle[i].Data()), 27 , binsx_p);
        
        
        //MassSq vs. momentum For signal extraction
        Mass2sqData[i] = (TH2F*)trackCutsData[i]->FindObject("Mass2sq_after");
        if(i<2)QAMassFit[i]  = new TH1F(Form("QAMassFit %s",Particle[i].Data()), Form("QAMassFit %s",Particle[i].Data()), 45,0,4.5);
        if(i>1)QAMassFit[i]  = new TH1F(Form("QAMassFit %s",Particle[i].Data()), Form("QAMassFit %s",Particle[i].Data()), 27 , binsx_p);
        
        //Data and MC for DCA Fits
        DCAData[i] = (TH2F*)(ListResults[i]->FindObject("DCAXYPBinningTot"));
        ListMCDCA[i] = (TList*)ListResultsMC[i]->FindObject("DCAPtBinning");
        PrimaryMC[i] = (TH2F*)ListMCDCA[i]->FindObject("DCAPtBinningPri");
        SecondaryMC[i] = (TH2F*)ListMCDCA[i]->FindObject("DCAPtBinningSec");
        MaterialMC[i] = (TH2F*)ListMCDCA[i]->FindObject("DCAPtBinningMat");
        if(i<2)fractions_primary[i] = new TH1F(Form("fractions primary %s",Particle[i].Data()), Form("fractions primary %s",Particle[i].Data()),45, 0, 4.5);
        if(i>1)fractions_primary[i] = new TH1F(Form("fractions primary %s",Particle[i].Data()), Form("fractions primary %s",Particle[i].Data()), 27 , binsx_p);
        
        gSystem->Exec(Form("mkdir %s",Particle[i].Data()));
    }
    
    //Extract the raw particle spectra from TPC singal and TOF m2 disctribution (RawParticleSpectra.C)
    RawParticleSpectraProton("Proton",MomentumData[0],RawSpectrum[0],Mass2sqData[0],QAMassFit[0]);
    RawParticleSpectraProton("Antiproton",MomentumData[1],RawSpectrum[1],Mass2sqData[1],QAMassFit[1]);
    RawParticleSpectraDeuteron("Deuteron",MomentumData[2],RawSpectrum[2],Mass2sqData[2],QAMassFit[2]);
    RawParticleSpectraDeuteron("Antideuteron",MomentumData[3],RawSpectrum[3],Mass2sqData[3],QAMassFit[3]);
    
    //Extract the primary fraction from DCAxy disctribution (Signal.C)
    MakePrimFracProton("Proton",fractions_primary[0],DCAData[0],PrimaryMC[0],SecondaryMC[0],MaterialMC[0]);
    MakePrimFracProton("Antiproton",fractions_primary[1],DCAData[1],PrimaryMC[1],SecondaryMC[1],MaterialMC[1]);
    MakePrimFracDeuteron("Deuteron",fractions_primary[2],DCAData[2]);
    
    //Get the primary spectra (Summary.C)
    PrimarySpectrum[0] = GetPrimarySpectra("Proton",RawSpectrum[0],fractions_primary[0]);
    PrimarySpectrum[1] = GetPrimarySpectra("Antiproton",RawSpectrum[1],fractions_primary[1]);
    PrimarySpectrum[2] = GetPrimarySpectra("Deuteron",RawSpectrum[2],fractions_primary[2]);
    PrimarySpectrum[3] = GetPrimarySpectra("Antideuteron",RawSpectrum[3],fractions_primary[3]);
    
    //Get the primary particle ratios (Summary.C)
    Ratio_primary[0] = GetPrimaryRatios("Antiproton/Proton",PrimarySpectrum[0],PrimarySpectrum[1]);
    Ratio_primary[1] = GetPrimaryRatios("Antideuteron/Deuteron",PrimarySpectrum[2],PrimarySpectrum[3]);
    TCanvas *cRatio = new TCanvas("cRatio","cRatio",0,0,1200,600);
    cRatio->Divide(2,1);
    for(Int_t i = 0; i < 2; i++){
        cRatio->cd(i+1);
        Ratio_primary[i]->Draw();
    }
    
    //Plot some Histograms to see if everything worked or not (Summary.C)
    MakeQAPlot("QAMassFit", "#it{p} (GeV/c)", "x100%",QAMassFit[0],QAMassFit[1],QAMassFit[2],QAMassFit[3]);
    MakeQAPlot("RawMomentum", "#it{p} (GeV/c)", "Counts",RawSpectrum[0],RawSpectrum[1],RawSpectrum[2],RawSpectrum[3]);
    MakeQAPlot("Primary Fraction", "#it{p} (GeV/c)", "x100%",fractions_primary[0],fractions_primary[1],fractions_primary[2],fractions_primary[3]);
    MakeQAPlot("RawPrimarySpectra", "#it{p} (GeV/c)", "Counts",PrimarySpectrum[0],PrimarySpectrum[1],PrimarySpectrum[2],PrimarySpectrum[3]);
    
    //Save the important histograms in .root files
    TString saveIn = "std";
    gSystem->Exec(Form("mkdir %s",saveIn.Data()));
    for (Int_t i = 0;i<4 ; i++){
        SaveRootFile(Form("%s/MassFitResult_%s.root",saveIn.Data(),Particle[i].Data()),QAMassFit[i]);
        SaveRootFile(Form("%s/Prim_Frac%s.root",saveIn.Data(),Particle[i].Data()),fractions_primary[i]);
        SaveRootFile(Form("%s/Prim_Spectra_%s.root",saveIn.Data(),Particle[i].Data()),PrimarySpectrum[i]);
    }
    SaveRootFile(Form("%s/Ratio_antip_p.root",saveIn.Data()),Ratio_primary[0]);
    SaveRootFile(Form("%s/Ratio_antid_d.root",saveIn.Data()),Ratio_primary[1]);
}
