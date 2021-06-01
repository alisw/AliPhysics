
/*_____________________________________________________________
 
 Run analysis macro for HF Correlations OnFlySim task
 Jitendra Kumar (jitendra.kumar@cern.ch)
 Andrea Rossi   (andrea.rossi@cern.ch)
 _____________________________________________________________*/


void SourceEnv_Libs(){
    
    gSystem->Setenv("GEN_TOTAL_EVENTS" , "50");
    gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_PHYSICS/macros -I$ALICE_PHYSICS/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF/ -I$ALICE_PHYSICS/PWGHF/base -I$ALICE_PHYSICS/PWG/FLOW/Base -I$ALICE_PHYSICS/PWG/FLOW/Tasks -I$ALICE_PHYSICS/PWGHF/vertexingHF -I$ALICE_PHYSICS/PWGHF/correlationHF -I$ALICE_PHYSICS/PWGHF/correlationHF/macros  -g");
    
    gSystem->Load("libOADB.so");
    gSystem->Load("libCORRFW.so");
    gSystem->Load("libPWGflowBase.so");
    gSystem->Load("libPWGflowTasks.so");
    gSystem->Load("libPWGHFbase.so");
    gSystem->Load("libPWGHFvertexingHF.so");
    gSystem->Load("libPWGHFcorrelationHF.so");
    
    // MC generator libraries
    gSystem->Load("liblhapdf");
    gSystem->Load("libpythia6_4_25");
    gSystem->Load("libqpythia");
    gSystem->Load("libEGPythia6");
    gSystem->Load("libAliPythia6");
    gSystem->Load("libAliPythia8");
}


//_____________________| Run Analysis
void RunHFCorrOnFlySim(Long64_t nEvents = 700000){
    
    SourceEnv_Libs();
    
    AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager","My Manager");
    mgr->SetDebugLevel(11);
    mgr->SetFileInfoLog("fileinfo.log");
    
    AliMCGenHandler* mcInputHandler = new AliMCGenHandler();
    
    const Int_t seed = atoi(gSystem->Getenv("CONFIG_SEED"));
    AliGenerator *genPythia = NULL;
    //genPythia = AddPythia6GenSettings("Perugia0");
    genPythia = AddPythia8GenSettings();
    mcInputHandler->SetGenerator(genPythia);
    mcInputHandler->SetSeedMode(1);
    mcInputHandler->SetSeed(seed);
    
    mgr->SetInputEventHandler(new AliDummyHandler());
    mgr->SetMCtruthEventHandler(mcInputHandler);
    
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGHF/correlationHF/macros/AddTaskHFCorrOnFlySim.C");
    AliAnalysisHFCorrOnFlySim *task = AddTaskHFCorrOnFlySim("Pythia8wBoost");
    task->SetDoCCbarOpeningAngleStudies(kTRUE);
    task->SetOpeningAngleEdges(1.,2.);    
    
    if(!mgr->InitAnalysis())return;
    mgr->PrintStatus();
    mgr->EventLoop(nEvents);
    return;
    
    
}



//__________________________________________| PYTHIA8
AliGenerator* AddPythia8GenSettings(){
    

    Float_t e_cms = 13000.0;
    Bool_t cr = kTRUE;
    Int_t kF=1, kProcess=0;
    Double_t ptHardMin = 0.0, ptHardMax = 1.0;
    Bool_t withBoost = kFALSE;

    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));

    
    AliGenPythiaPlus* gener = new AliGenPythiaPlus(AliPythia8::Instance());
    if(kProcess==0) gener->SetProcess(kPyMbDefault);
    else if(kProcess==1) {
        gener->SetProcess(kPyJets);
        if(ptHardMin>0.)AliPythia8::Instance()->SetPtHardRange(ptHardMin,ptHardMax);
    }

    gener->SetEnergyCMS(e_cms); // in GeV
    //gener->SetEventListRange(-1, -1);
    
    if(withBoost){
        gener->SetProjectile("p",1,1);
        gener->SetTarget("n",208,82);
        //genP->SetUseLorentzBoost(kTRUE);
    }
    
    AliPythia8::Instance()->ReadString("Tune:pp = 5");//CR
    AliPythia8::Instance()->ReadString("Random:setSeed = on");
    AliPythia8::Instance()->ReadString("Random:seed = 0");
    if(cr)(AliPythia8::Instance())->ReadString("ColourReconnection:reconnect = on");
    else (AliPythia8::Instance())->ReadString("ColourReconnection:reconnect = off");
    
    AliPythia8::Instance()->ReadString(Form("MultipartonInteractions:kFactor = %i", kF));
    
    gener->Init();
    return gener;
    
}




//__________________________________________| PYTHIA6
AliGenerator* AddPythia6GenSettings(TString PythiaTune = "Perugia2011"){
    
    Float_t e_cms = 13000.0;
    Bool_t cr = kTRUE;
    Bool_t withBoost = kFALSE;

    AliGenPythia* genP = new AliGenPythia(-1);
    //genP->SetMomentumRange(0, 999999.);
    //genP->SetThetaRange(0., 180.);
    //genP->SetYRange(-20.,20.);
    //genP->SetPtRange(0,1000.);
    //genP->SetVertexSmear(kPerEvent);
    genP->SetProcess(kPyMbDefault);
    genP->SetEnergyCMS(e_cms);
    //genP->SetEventListRange(-1,-1);
    if(withBoost){
        genP->SetProjectile("p",1,1);
        genP->SetTarget("n",208,82);
        genP->SetUseLorentzBoost(kTRUE);
        //genP->SetNuclearPDF(19);
        //genP->SetUseNuclearPDF(kTRUE);
    }

    if(PythiaTune == "Perugia0"){
        genP->SetTune(320);
        if(!cr) genP->SetTune(324);
        cout << "---Perugia0 tunes--" << endl;
    
    }else if(PythiaTune == "Perugia2010"){
        genP->SetTune(327);
        if(!cr) genP->SetTune(324);
        cout << "---Perugia2010 tunes--" << endl;
        
    }else if(PythiaTune == "Perugia2011"){
        genP->SetTune(350);
        if(!cr) genP->SetTune(354);
        cout << "---Perugia2011 tunes--" << endl;
        
    }else if(PythiaTune == "Perugia2012"){
        genP->SetTune(370);
        if(!cr)genP->SetTune(375);
        cout << "---Perugia2012 tunes--" << endl;

    }else cout << "Select Proper tunes" << endl;
    
    //genP->UseNewMultipleInteractionsScenario();
        
    //genP->SetCrossingAngle(0,0.000280); //! WARNING for crossing angle seeting (Phi distributions weird)

    //Float_t sigmaz = 6.245; // [cm]
    //Float_t sigmax = 0.0025;
    //Float_t sigmay = 0.0029;
    //genP->SetOrigin(0., 0., 0.); // Taken from OCDB
    //genP->SetSigma(sigmax, sigmay, sigmaz);

    genP->Init();
    //genP->Print();
    return genP;
}


