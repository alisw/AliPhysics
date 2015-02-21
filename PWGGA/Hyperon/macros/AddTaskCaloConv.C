AliAnalysisTaskCaloConv * AddTaskCaloConv(){
    
  
    
    //Macro to add class CaloConv (conversion+calorimeters pi0 analysis) to train
    //Argument is the path to the PHOS recalibration parameters (file with OCDB entry)
    //Default path to the file with unit recalibration == no recalibnration
    //If file does not exist, no recalibration histograms will be filled

  Int_t isHeavyIon = 0;
    
    // Get the pointer to the existing analysis manager via the static access method.
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskCaloConv", "No analysis manager to connect to.");
        return NULL;
    }
    
    // Check the analysis type using the event handlers connected to the analysis manager.

    // ================== GetInputEventHandler =============================
    AliVEventHandler *inputHandler=mgr->GetInputEventHandler();
    
    if (!mgr->GetInputEventHandler()) {
        ::Error("AddTaskCaloConv", "This task requires an input event handler");
        return NULL;
    }
    
    //=========  Set Cutnumber for V0Reader ================================
    /*    TString cutnumber = "0000000002084000002200000000";
	  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer(); */
        //=========  Set Cutnumber for V0Reader ================================
        TString cutnumberPhoton = "002084000002200000000";
        TString cutnumberEvent = "0000000";
        AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();


    
    //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
    if( !(AliV0ReaderV1*)mgr->GetTask("V0ReaderV1") ){
        AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1("V0ReaderV1");
        
        fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
        //fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
        //fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);
        
        if (!mgr) {
            Error("AddTaskCaloConv", "No analysis manager found.");
            return;
        }
        
        // Set AnalysisCut Number
	/*        AliConversionCuts *fCuts=NULL;
        if(cutnumber!=""){
            fCuts= new AliConversionCuts(cutnumber.Data(),cutnumber.Data());
            fCuts->SetPreSelectionCutFlag(kTRUE);
            if(fCuts->InitializeCutsFromCutString(cutnumber.Data())){
                fV0ReaderV1->SetConversionCuts(fCuts);
                fCuts->SetFillCutHistograms("",kTRUE);
            }
	    } */
 
                AliConvEventCuts *fEventCuts=NULL;
                if(cutnumberEvent!=""){
                        fEventCuts= new AliConvEventCuts(cutnumberEvent.Data(),cutnumberEvent.Data());
                        fEventCuts->SetPreSelectionCutFlag(kTRUE);
                        if(fEventCuts->InitializeCutsFromCutString(cutnumberEvent.Data())){
                                fV0ReaderV1->SetEventCuts(fEventCuts);
                                fEventCuts->SetFillCutHistograms("",kTRUE);
                        }
                }

                
                // Set AnalysisCut Number
                AliConversionPhotonCuts *fCuts=NULL;
                if(cutnumberPhoton!=""){
                        fCuts= new AliConversionPhotonCuts(cutnumberPhoton.Data(),cutnumberPhoton.Data());
                        fCuts->SetPreSelectionCutFlag(kTRUE);
                        fCuts->SetIsHeavyIon(isHeavyIon);
                        if(fCuts->InitializeCutsFromCutString(cutnumberPhoton.Data())){
                                fV0ReaderV1->SetConversionCuts(fCuts);
                                fCuts->SetFillCutHistograms("",kTRUE);
                        }
                }



       fV0ReaderV1->Init();
        AliLog::SetGlobalLogLevel(AliLog::kInfo);
        
        
        //connect input V0Reader
        mgr->AddTask(fV0ReaderV1);
        mgr->ConnectInput(fV0ReaderV1,0,cinput);
        
    }
 
    
    AliAnalysisTaskCaloConv *task = new AliAnalysisTaskCaloConvCorr("CaloConv");
    
    //mgr->AddTask(task);
 /*
    TDirectory* saveDir = gDirectory;
    TGrid a;
    
    // TFile *fBadMap = TFile::Open("./EMCALBadChannels.root") ;
    
    //  TFile *fBadMap = TFile::Open("./BadMap_LHC10b.root");
    if(a.IsConnected()){
        // if ( 1>0 ) {
        
        TFile *fBadMap = TFile::Open("alien:///alice/cern.ch/user/p/prsnko/BadMaps/BadMap_LHC10b.root") ;
        //    TFile *fBadMap = TFile::Open("./BadMap_LHC10b.root");
        
        //    TFile *fBadMap = TFile::Open("../BadMap_LHC11a_pp2760.root");
        // jisongTFile *fBadMap = TFile::Open("../BadMap_LHC10de_Majority.root");
        
        if(fBadMap->IsOpen()){
            if (saveDir) saveDir->cd(); else gROOT->cd();
            
            printf("Adding PHOS and EMCAL bad maps \n") ;
            char key[55] ;
            for(Int_t mod=1;mod<4; mod++){
                sprintf(key,"PHOS_BadMap_mod%d",mod) ;
                TH2I * h = (TH2I*)fBadMap->Get(key) ;
                if(h)
                    task->SetPHOSBadMap(mod,h) ;
            }
            for(Int_t sm=0; sm<5; sm++){
                sprintf(key,"EMCAL_BadMap_mod%d",sm) ;
                TH2I * h = (TH2I*)fBadMap->Get(key) ;
                if(h)
                    task->SetEMCALBadMap(mod,h) ;
            }
            
            fBadMap->Close() ;
        }
    }
    else{
        printf("Can not open Bad Map file \n") ;
    }
    
   */
    // Create containers for input/output
     AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    //mgr->ConnectInput(task, 0, cinput);
    TString outputfile = AliAnalysisManager::GetCommonFileName();
    
    
    AliAnalysisDataContainer *coutput = mgr->CreateContainer("CaloConv", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:PWGGA_CaloConv",outputfile.Data()));
    
  
    
    mgr->AddTask(task);
    
    mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutput);
   // mgr->ConnectOutput(task, 2, coutput2);
    
    return task ;
    
}
