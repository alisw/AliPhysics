void CreateAOD()
{
  // Main
  
  // LoadLibraries
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libXMLIO.so");
  gSystem->Load("libMatrix.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWG3base.so");
  gSystem->Load("libPWG3muon.so");
  gSystem->Load("libTENDER.so");
//gSystem->Load("libTENDERSupplies.so"); 

  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_ROOT/PHOS");
  gROOT->LoadMacro("AliPHOSDigitDecalibrate.cxx++") ;

  //Data
  TChain *chain       = new TChain("esdTree") ;
  chain->Add("AliESDs.root");
  
  if(chain){
    AliLog::SetGlobalLogLevel(AliLog::kError);//Minimum prints on screen
    
    //--------------------------------------
    // Make the analysis manager
    //-------------------------------------
    AliAnalysisManager *mgr  = new AliAnalysisManager("Manager", "Manager");
    
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mcHandler->SetReadTR(kFALSE);//Do not search TrackRef file
    mgr->SetMCtruthEventHandler(mcHandler);
    
    AliAODHandler* aodoutHandler   = new AliAODHandler();
    aodoutHandler->SetOutputFileName("AliAOD.root");
    ////aodoutHandler->SetCreateNonStandardAOD();
    mgr->SetOutputEventHandler(aodoutHandler);
    
    AliESDInputHandler *esdHandler = new AliESDInputHandler();
    esdHandler->SetReadFriends(kFALSE);
    mgr->SetInputEventHandler(esdHandler);
    
    //mgr->SetDebugLevel(-1); // For debugging, do not uncomment if you want no messages.
    
    //-------------------------------------------------------------------------
    //Define task, put here any other task that you want to use.
    //-------------------------------------------------------------------------
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1;
    
    coutput1 = mgr->GetCommonOutputContainer();



    //========= Add tender to the ANALYSIS manager and set default storage =====

    AliTender *tender=new AliTender("AnalysisTender");
    tender->SetDefaultCDBStorage("local://$ALICE_ROOT/OCDB"); //comment this you acess grid
    mgr->AddTask(tender);

    //========= Attach PHOS supply ======
    AliPHOSDigitDecalibrate *PHOSSupply=new AliPHOSDigitDecalibrate("PHOStender");
    tender->AddSupply(PHOSSupply);
    TFile * fDecalib = TFile::Open("Decalibraiton.root") ;
    char key[55] ;
    for(Int_t mod=0; mod<5; mod++){
      sprintf(key,"DecalibrMod%d",mod) ;
      TH2F * h = (TH2F*)fDecalib->Get(key) ;
      PHOSSupply->SetDecalibration(mod,h) ;
    }
    fDecalib->Close() ;

    //================================================
    //              data containers
    //================================================

    //            define output containers, please use 'username'_'somename'
    AliAnalysisDataContainer *coutput1 =
      mgr->CreateContainer("tender_event", AliESDEvent::Class(),
                           AliAnalysisManager::kExchangeContainer,"default_tender");
 
    //           connect containers
    mgr->ConnectInput  (tender,  0, mgr->GetCommonInputContainer() );
    mgr->ConnectOutput (tender,  1, coutput1);






    
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskESDFilter.C");
    AliAnalysisTaskESDfilter *taskesdfilter = AddTaskESDFilter(kTRUE);
    
    //-----------------------
    // Run the analysis
    //-----------------------    
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->StartAnalysis("local",chain);
    
    cout <<" Analysis ended sucessfully "<< endl ;
    
  }
  else cout << "Chain was not produced ! "<<endl;
  
}
