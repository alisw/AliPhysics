AliAnalysisTaskSEDvsEventShapes *AddTaskDvsEventShapes(Int_t system=0,
                                                       Bool_t readMC=kFALSE,
                                                       Int_t MCOption=0,
                                                       Int_t pdgMeson=411,
                                                       TString finDirname="Loose",
                                                       TString filename="",
                                                       TString finAnObjname="AnalysisCuts",
                                                       Bool_t CalculateSphericity=kFALSE, // set true to calculate sphericity
                                                       Int_t SoSparseChecks=0,/*0 mult, 1 multUncorr, 2 NoPid, 3 All*/
                                                       Double_t ptMin=0.15,
                                                       Double_t ptMax=10.,
                                                       Double_t etaMin=-0.8,
                                                       Double_t etaMax=0.8,
                                                       Int_t minMult=3,
                                                       Double_t phiStepSizeDeg=0.1,
                                                       Int_t filtbit1=256,
                                                       Int_t filtbit2=512,
                                                       TString estimatorFilename="",
                                                       Double_t refMult=9.26,
                                                       Bool_t subtractDau=kFALSE,
                                                       Bool_t subtractDauFromSphero=kFALSE,
                                                       Int_t NchWeight=0,
                                                       Int_t recoEstimator = AliAnalysisTaskSEDvsEventShapes::kNtrk10,
                                                       Int_t MCEstimator = AliAnalysisTaskSEDvsEventShapes::kEta10,
                                                       Bool_t isPPbData=kFALSE)
{
    //
    // Macro for the AliAnalysisTaskSE for D candidates vs Multiplicity as a function of Event shape variables
    // Invariant mass histogram in Sphero(i)city, pt and multiplicity bins in a 3D histogram
    //   different estimators implemented
    //==============================================================================
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskDvsEventShapes", "No analysis manager to connect to.");
    }
    
    Bool_t stdcuts=kFALSE;
    TFile* filecuts;
    if( filename.EqualTo("") ) {
        stdcuts=kTRUE;
    } else {
        filecuts=TFile::Open(filename.Data());
        if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
            AliFatal("Input file not found : check your cut object");
        }
    }
    
    
    //Analysis Task
    AliRDHFCuts *analysiscuts=0x0;
    
    TString Name="";
    if(pdgMeson==411){
        if(stdcuts) {
            analysiscuts = new AliRDHFCutsDplustoKpipi();
            if (system == 0) analysiscuts->SetStandardCutsPP2010();
            else analysiscuts->SetStandardCutsPbPb2011();
        }
        else analysiscuts = (AliRDHFCutsDplustoKpipi*)filecuts->Get(finAnObjname);
        Name="Dplus";
    }else if(pdgMeson==421){
        if(stdcuts) {
            analysiscuts = new AliRDHFCutsD0toKpi();
            if (system == 0) analysiscuts->SetStandardCutsPP2010();
            else analysiscuts->SetStandardCutsPbPb2011();
        }
        else analysiscuts = (AliRDHFCutsD0toKpi*)filecuts->Get(finAnObjname);
        Name="D0";
    }else if(pdgMeson==413){
        if(stdcuts) {
            analysiscuts = new AliRDHFCutsDStartoKpipi();
            if (system == 0) analysiscuts->SetStandardCutsPP2010();
            else analysiscuts->SetStandardCutsPbPb2011();
        }
        else analysiscuts = (AliRDHFCutsDStartoKpipi*)filecuts->Get(finAnObjname);
        Name="DStar";
    }
    
    AliAnalysisTaskSEDvsEventShapes *dMultTask = new AliAnalysisTaskSEDvsEventShapes("dEvtShapeAnalysis",pdgMeson,analysiscuts,isPPbData);
    dMultTask->SetReadMC(readMC);
    dMultTask->SetDebugLevel(0);
    dMultTask->SetUseBit(kTRUE);
    dMultTask->SetDoImpactParameterHistos(kFALSE);
    dMultTask->SetFillSoSparseForMultUncorrNoPid(SoSparseChecks); //Set Fill THnSparse for Spherocity
    dMultTask->SetEventShapeParameters(ptMin, ptMax, etaMin, etaMax, minMult, phiStepSizeDeg, filtbit1, filtbit2); //parameters to calculate Sphero(i)city
    dMultTask->SetCalculationsForSphericity(CalculateSphericity);
    dMultTask->SetSubtractTrackletsFromDaughters(subtractDau);
    dMultTask->SetRecomputeSpherocityWithoutDau(subtractDauFromSphero);
    dMultTask->SetMultiplicityEstimator(recoEstimator);
    dMultTask->SetMCPrimariesEstimator(MCEstimator);
    dMultTask->SetMCOption(MCOption);
    if(isPPbData) dMultTask->SetIsPPbData();
    
    if(NchWeight){
        TH1F *hNchPrimaries = NULL;
        TH1F *hMeasNchPrimaries = NULL;
        if(NchWeight==1){
            if(isPPbData) {
                hNchPrimaries = (TH1F*)filecuts->Get("hNtrUnCorrEvWithDWeight"); // MC distribution
            }
            else hNchPrimaries = (TH1F*)filecuts->Get("hGenPrimaryParticlesInelGt0");
            if(hNchPrimaries) {
                dMultTask->UseMCNchWeight(NchWeight);
                dMultTask->SetHistoNchWeight(hNchPrimaries);
            } else {
                AliFatal("Histogram for Nch multiplicity weights not found");
                return 0x0;
            }
            hMeasNchPrimaries = (TH1F*)filecuts->Get("hMeasNtrUnCorrEvWithD"); // data distribution
            if(hMeasNchPrimaries) {
                dMultTask->SetMeasuredNchHisto(hMeasNchPrimaries);
            }
        }
        else if(NchWeight==2){
            hNchPrimaries = (TH1F*)filecuts->Get("hNtrUnCorrEvWithDWeight"); // MC distribution
            hMeasNchPrimaries = (TH1F*)filecuts->Get("hMeasNtrUnCorrEvWithD"); // data distribution
            if(hNchPrimaries && hMeasNchPrimaries) {
                dMultTask->UseMCNchWeight(NchWeight);
                dMultTask->SetHistoNchWeight(hNchPrimaries);
                dMultTask->SetMeasuredNchHisto(hMeasNchPrimaries);
            } else {
                AliFatal("Histogram for Ntrk multiplicity weights not found");
                return 0x0;
            }
        }
    }
    
    
    if(pdgMeson==421) {
        dMultTask->SetMassLimits(1.5648,2.1648);
        dMultTask->SetNMassBins(200);
    }else if(pdgMeson==411)dMultTask->SetMassLimits(pdgMeson,0.2);
    
    if(estimatorFilename.EqualTo("") ) {
        printf("Estimator file not provided, multiplcity corrected histograms will not be filled\n");
    } else{
        
        TFile* fileEstimator=TFile::Open(estimatorFilename.Data());
        if(!fileEstimator)  {
            AliFatal("File with multiplicity estimator not found\n");
            return;
        }
        
        dMultTask->SetReferenceMultiplcity(refMult);
        
        const Char_t* profilebasename="SPDmult10";
        if(recoEstimator==AliAnalysisTaskSEDvsEventShapes::kVZEROA || recoEstimator==AliAnalysisTaskSEDvsEventShapes::kVZEROAEq) profilebasename="VZEROAmult";
        else if(recoEstimator==AliAnalysisTaskSEDvsEventShapes::kVZERO || recoEstimator==AliAnalysisTaskSEDvsEventShapes::kVZEROEq) profilebasename="VZEROMmult";
        cout<<endl<<endl<<" profilebasename="<<profilebasename<<endl<<endl;
        
        if (isPPbData) {    //Only use two profiles if pPb
            const Char_t* periodNames[2] = {"LHC13b", "LHC13c"};
            TProfile* multEstimatorAvg[2];
            for(Int_t ip=0; ip<2; ip++) {
                cout<< " Trying to get "<<Form("%s_%s",profilebasename,periodNames[ip])<<endl;
                multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
                if (!multEstimatorAvg[ip]) {
                    AliFatal(Form("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]));
                    return;
                }
            }
            dMultTask->SetMultiplVsZProfileLHC13b(multEstimatorAvg[0]);
            dMultTask->SetMultiplVsZProfileLHC13c(multEstimatorAvg[1]);
        }
        else {
            const Char_t* periodNames[4] = {"LHC10b", "LHC10c", "LHC10d", "LHC10e"};
            TProfile* multEstimatorAvg[4];
            for(Int_t ip=0; ip<4; ip++) {
                multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
                if (!multEstimatorAvg[ip]) {
                    AliFatal(Form("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]));
                    return;
                }
            }
            dMultTask->SetMultiplVsZProfileLHC10b(multEstimatorAvg[0]);
            dMultTask->SetMultiplVsZProfileLHC10c(multEstimatorAvg[1]);
            dMultTask->SetMultiplVsZProfileLHC10d(multEstimatorAvg[2]);
            dMultTask->SetMultiplVsZProfileLHC10e(multEstimatorAvg[3]);
        }
    }
    mgr->AddTask(dMultTask);
    
    // Create containers for input/output
    
    TString inname = "cinput";
    TString outname = "coutput";
    TString cutsname = "coutputCuts";
    TString normname = "coutputNorm";
    TString profname = "coutputProf";
    TString effname = "coutputEffCorr";
    
    inname += Name.Data();
    outname += Name.Data();
    cutsname += Name.Data();
    normname += Name.Data();
    profname += Name.Data();
    effname += Name.Data();
    inname += finDirname.Data();
    outname += finDirname.Data();
    cutsname += finDirname.Data();
    normname += finDirname.Data();
    profname += finDirname.Data();
    effname += finDirname.Data();
    
    AliAnalysisDataContainer *cinput = mgr->CreateContainer(inname,TChain::Class(),AliAnalysisManager::kInputContainer);
    
    TString outputfile = AliAnalysisManager::GetCommonFileName();
    outputfile += ":PWG3_D2H_DEvtShape_";
    outputfile += Name.Data(); 
    outputfile += finDirname.Data(); 
    
    AliAnalysisDataContainer *coutputCuts = mgr->CreateContainer(cutsname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutput = mgr->CreateContainer(outname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutputNorm = mgr->CreateContainer(normname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutputProf = mgr->CreateContainer(profname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    if(readMC) AliAnalysisDataContainer *coutputEffCorr = mgr->CreateContainer(effname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());

    
    mgr->ConnectInput(dMultTask,0,mgr->GetCommonInputContainer());
    
    mgr->ConnectOutput(dMultTask,1,coutput);
    mgr->ConnectOutput(dMultTask,2,coutputCuts);
    mgr->ConnectOutput(dMultTask,3,coutputNorm);
    mgr->ConnectOutput(dMultTask,4,coutputProf);
    if(readMC) mgr->ConnectOutput(dMultTask,5,coutputEffCorr);
    
    return dMultTask;
}
