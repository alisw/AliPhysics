//Create by Christine Nattrass
//University of Tennessee at Knoxville

//by default this runs locally
//With the argument true this submits jobs to the grid
//As written this requires an xml script tag.xml in the ~/et directory on the grid to submit jobs
void runHadEt(bool submit = false, bool data = false, Int_t dataset = 20100, Int_t test = 1, Int_t material = 0, Bool_t altV0Scale = kFALSE, bool runCompiledVersion = kTRUE, int simflag = 0,Bool_t finecentbins = kFALSE) {
    TStopwatch timer;
    timer.Start();
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libXMLIO.so");

    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libAOD.so");

    gSystem->Load("libANALYSIS");
    gSystem->Load("libOADB.so");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libANALYSISaliceBase");
    gSystem->Load("libCORRFW");
    gSystem->Load("libPWGUDbase.so");//libANALYSISaliceBase, libCORRFW, libOADB

//     gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$CMAKE_INSTALL_PREFIX/include");
//     gSystem->AddIncludePath("-I$ALICE_ROOT/PWGUD/base");
    if(runCompiledVersion){
    gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
    gSystem->AddIncludePath("-I$ALICE_PHYSICS/PWGUD ${AliPhysics_SOURCE_DIR}/OADB/COMMON/MULTIPLICITY");
      cout<<"Using compiled version"<<endl;
      gSystem->Load("libPWGLFtotEt.so");
    }
    else{
      cout<<"Not using compiled version"<<endl;
      gROOT->LoadMacro("macros/loadLocalHadLibraries.C");
      loadLocalHadLibraries();
//   gROOT->ProcessLine(".L AliAnalysisEtCutsLocal.cxx+g");
//   gROOT->ProcessLine(".L AliAnalysisHadEtCorrectionsLocal.cxx+g");
//    gROOT->ProcessLine(".L AliAnalysisEtCommonLocal.cxx+g");
//   gROOT->ProcessLine(".L AliAnalysisHadEtLocal.cxx+g");
//      gROOT->ProcessLine(".L AliAnalysisHadEtMonteCarloLocal.cxx+g");
//      gROOT->ProcessLine(".L AliAnalysisHadEtReconstructedLocal.cxx+g");
//        gROOT->ProcessLine(".L AliAnalysisTaskTransverseEnergyLocal.cxx+g");
//        gROOT->ProcessLine(".L AliAnalysisTaskHadEtLocal.cxx+g");
    }

  char *kTreeName = "esdTree" ;
  TChain * chain   = new TChain(kTreeName,"myESDTree") ;
  if(submit){      
    gSystem->Load("libNetx.so") ; 
    gSystem->Load("libgapiUI.so");
    gSystem->Load("libRAliEn.so"); 
    TGrid::Connect("alien://") ;
  }
  bool PbPb = false;
  if(dataset ==20100 || dataset==2011 || dataset==2015){
    bool PbPb = true;
    if(data){
      cout<<"Testing new data"<<endl;
//        chain->Add("/data/LHC15o/000244827/15000244827038.113/AliESDs.root");
       chain->Add("/data/LHC10h/pass2/10000137366041.860/AliESDs.root");//Data Pb+Pb
       chain->Add("/data/LHC10h/pass2/10000137366041.870/AliESDs.root");//Data Pb+Pb
       chain->Add("/data/LHC10h/pass2/10000137366041.880/AliESDs.root");//Data Pb+Pb
//        chain->Add("/data/LHC10h/pass2/10000137366041.890/AliESDs.root");//Data Pb+Pb
//        chain->Add("/data/LHC10h/pass2/10000137366041.900/AliESDs.root");//Data Pb+Pb
    }
    else{
      //chain->Add("/data/LHC12d3/168464/201/AliESDs.root");//HIJING with embedded signals - works, full acceptance
      //chain->Add("/data/LHC10h2/137161/001/AliESDs.root");//DPMJET
      //chain->Add("/data/LHC11a9a/137366/001/AliESDs.root");//AMPT
//       chain->Add("/data/LHC10h8/137161/999/AliESDs.root");//Hijing Pb+Pb
//       chain->Add("/data/LHC10h8/137161/111/AliESDs.root");//Hijing Pb+Pb
//       chain->Add("/data/LHC10h8/137161/222/AliESDs.root");//Hijing Pb+Pb
      //chain->Add("/data/LHC11a4_bis/137161/999/AliESDs.root");//Hijing Pb+Pb
      //chain->Add("/data/LHC10h12/999/AliESDs.root");//Hijing Pb+Pb
      cout<<"I am here "<<endl;
      //chain->Add("/data/LHC15k1/225717/999/AliESDs.root");//HIJING 2015
 
      //       chain->Add("/data/LHC15k1_plus/244918/9999/AliESDs.root");
//       chain->Add("/data/LHC15k1/244918/9999/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/001/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/002/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/003/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/004/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/005/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/006/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/007/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/008/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/009/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/010/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/012/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/013/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/014/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/015/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/016/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/017/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/018/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/019/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/020/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/021/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/022/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/023/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/024/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/025/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/026/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/027/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/029/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/030/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/031/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/032/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/034/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/035/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/036/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/037/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/038/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/039/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/040/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/041/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/042/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/043/AliESDs.root");
// chain->Add("/data/LHC15k1/225717/243984/044/AliESDs.root");

      chain->Add("/data/LHC11a10a_bis/139465/003/AliESDs.root");
 chain->Add("/data/LHC11a10a_bis/139465/004/AliESDs.root");
 chain->Add("/data/LHC11a10a_bis/139465/006/AliESDs.root");
//  chain->Add("/data/LHC11a10a_bis/139465/007/AliESDs.root");
//  chain->Add("/data/LHC11a10a_bis/139465/008/AliESDs.root");
//  chain->Add("/data/LHC11a10a_bis/139465/009/AliESDs.root");
//  chain->Add("/data/LHC11a10a_bis/139465/010/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/011/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/012/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/013/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/014/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/015/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/016/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/017/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/018/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/019/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/020/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/021/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/022/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/023/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/024/AliESDs.root");

// chain->Add("/data/LHC11a10a_bis/139465/025/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/026/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/027/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/028/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/029/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/030/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/031/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/032/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/033/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/034/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/035/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/036/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/037/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/038/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/039/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/040/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/041/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/042/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/043/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/044/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/045/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/046/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/047/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/048/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/049/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/050/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/051/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/052/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/053/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/054/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/055/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/056/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/057/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/058/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/059/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/060/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/061/AliESDs.root");
// chain->Add("/data/LHC11a10a_bis/139465/062/AliESDs.root");

    }
  } 
  else{
     if(dataset==2013){//pPb 5 TeV
       if(data){
	 chain->Add("/data/LHC13b/13000195483082.95/AliESDs.root");
       }
       else{
	 //cout<<"Yes I am analyzing the correct file"<<endl;
	 //chain->Add("/data/LHC12c1b/111/AliESDs.root");
	 chain->Add("/data/LHC13b3/9999/AliESDs.root");
       }
     }
     else{
       if(data){
	 //chain->Add("/data/LHC10dpass2/10000126403050.70/AliESDs.root");//data
	 //chain->Add("/data/LHC11a/11000146856042.90/AliESDs.root");//data pp 2.76 TeV w/SDD pass 2
	 //chain->Add("/data/LHC11a/11000146860043.90/AliESDs.root");//data pp 2.76 TeV w/SDD pass 3
	 chain->Add("/data/LHC10c/10000120683048.90/AliESDs.root");//Data Pb+Pb
       }
       else{
	 cerr<<"Hello I am here 75"<<endl;
	 //chain->Add("/data/LHC10d15/1821/AliESDs.root");//simulation p+p
	 chain->Add("/data/LHC11b10a/001/AliESDs.root");
       }
     }
  }

  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TotEtManager");
  if(submit){

    gROOT->LoadMacro("CreateAlienHandlerHadEt.C");
    AliAnalysisGrid *alienHandler = CreateAlienHandlerHadEt(dataset,data,test,material,altV0Scale,runCompiledVersion,simflag,finecentbins);//integer dataset, boolean isData, bool submit-in-test-mode, bool use alternatve V0 scaling
      if (!alienHandler) return;
      mgr->SetGridHandler(alienHandler);
  }

  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  AliMCEventHandler* handler = new AliMCEventHandler;
  if(!data){
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }


  gROOT->LoadMacro(" $ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(!data);
  if(!physSelTask) { Printf("no physSelTask"); return; }
  AliPhysicsSelection *physSel = physSelTask->GetPhysicsSelection();
  //physSel->AddCollisionTriggerClass("+CINT1B-ABCE-NOPF-ALL");// #3119 #769");
  //physSelTask->AddCollisionTriggerClass(AliVEvent::kINT7);// #3119 #769");


  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask * multtask;// = AddTaskMultSelection(kFALSE); // user mode:

  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");

  AliCentralitySelectionTask *centTask;
  
  if(PbPb){
    //if(2015){
      multtask = AddTaskMultSelection(kFALSE); // user mode:
      multtask->SetSelectedTriggerClass(AliVEvent::kINT7);
      //if(!data){
	//cout<<"Setting up for MC LHC15k1"<<endl;
	//multtask -> SetAlternateOADBforEstimators ( "LHC15k1" );
      //}
//     }
//     else{
//       AliCentralitySelectionTask *centTask = AddTaskCentrality();
//       if(!data){
// 	cout<<"Setting up centrality for MC"<<endl;
// 	centTask->SetMCInput();
//       }
//       else{
// 	cout<<"Setting up centrality for data"<<endl;
//       }
//     }
  }
  Int_t pass = 2;
   if(dataset==20100){//PbPb 2.76 TeV
     if(simflag==0){
       if(finecentbins){
	 if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC11a10a_bis.FineCentralityBins.PbPb.139465.ForData.root","corrections.root",kTRUE);}
	 else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC11a10a_bis.FineCentralityBins.PbPb.139465.ForSimulations.root","corrections.root",kTRUE);}
       }
       else{
	 if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC11a10a_bis.PbPb.139465.ForData.root","corrections.root",kTRUE);}
	 else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC11a10a_bis.PbPb.139465.ForSimulations.root","corrections.root",kTRUE);}
	 if(material==11){//AMPT only available for this run
	   cout<<"Running AMPT"<<endl;
	   gSystem->CopyFile("rootFiles/corrections/corrections.LHC11a10a_bis.PbPb.137366.ForData.root","corrections.root",kTRUE);
	 }
       }
     }
     if(simflag==1){
       if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC11a10a_bis.PbPb.Run138442.ForData.root","corrections.root",kTRUE);}
       else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC11a10a_bis.PbPb.Run138442.ForSimulations.root","corrections.root",kTRUE);}
     }
     if(simflag==2){
       if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC11a10a_bis.PbPb.Run138364.ForData.root","corrections.root",kTRUE);}
       else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC11a10a_bis.PbPb.Run138364.ForSimulations.root","corrections.root",kTRUE);}
     }
     if(simflag==3){
       if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC11a10a_bis.PbPb.Run138396.ForData.root","corrections.root",kTRUE);}
       else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC11a10a_bis.PbPb.Run138396.ForSimulations.root","corrections.root",kTRUE);}
     }
     if(simflag==4){
       if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC11a10a_bis.PbPb.Run137722.ForData.root","corrections.root",kTRUE);}
       else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC11a10a_bis.PbPb.Run137722.ForSimulations.root","corrections.root",kTRUE);}
     }
     if(simflag==5){
       if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC11a10a_bis.PbPb.Run137366.ForData.root","corrections.root",kTRUE);}
       else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC11a10a_bis.PbPb.Run137366.ForSimulations.root","corrections.root",kTRUE);}
     }
     if(simflag==6){
       if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC11a10a_bis.PbPb.Run137161.ForData.root","corrections.root",kTRUE);}
       else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC11a10a_bis.PbPb.Run137161.ForSimulations.root","corrections.root",kTRUE);}
     }
     gSystem->CopyFile("ConfigHadEtMonteCarloPbPb.C","ConfigHadEtMonteCarlo.C",kTRUE);
     gSystem->CopyFile("ConfigHadEtReconstructedPbPb.C","ConfigHadEtReconstructed.C",kTRUE);
     if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC11a10a_bis.PbPb.ForData.root","corrections.root",kTRUE);}
     else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC11a10a_bis.PbPb.ForSimulations.root","corrections.root",kTRUE);}
     pass = 2;
     //centTask->SetPass(1);
   }
   else{
     if(dataset==2011){//PbPb 2.76 TeV 2011
       gSystem->CopyFile("ConfigHadEtMonteCarloPbPb2011.C","ConfigHadEtMonteCarlo.C",kTRUE);
       gSystem->CopyFile("ConfigHadEtReconstructedPbPb2011.C","ConfigHadEtReconstructed.C",kTRUE);
       if(finecentbins){
	 if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC13e1abc.FineCentralityBins.PbPb.168464.ForData.root","corrections.root",kTRUE);}
	 else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC13e1abc.FineCentralityBins.PbPb.168464.ForSimulations.root","corrections.root",kTRUE);}
       }
       else{
	 if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC13e1abc.PbPb.168464.ForData.root","corrections.root",kTRUE);}
	 else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC13e1abc.PbPb.168464.ForSimulations.root","corrections.root",kTRUE);}
       }
     }
     if(dataset==2015){//PbPb 5 TeV 2015
       cout<<"I am here!  I should "<<endl;
       pass = 1;
       gSystem->CopyFile("ConfigHadEtMonteCarloPbPb2015.C","ConfigHadEtMonteCarlo.C",kTRUE);
       gSystem->CopyFile("ConfigHadEtReconstructedPbPb2015.C","ConfigHadEtReconstructed.C",kTRUE);
//        if(finecentbins){
// 	 if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC13e1abc.FineCentralityBins.PbPb.168464.ForData.root","corrections.root",kTRUE);}
// 	 else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC13e1abc.FineCentralityBins.PbPb.168464.ForSimulations.root","corrections.root",kTRUE);}
//        }
//        else{
// 	 if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1_plus.PbPb.244918.ForData.2015.root","corrections.root",kTRUE);}
// 	 else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1_plus.PbPb.244918.ForSimulations.2015.root","corrections.root",kTRUE);}
       
       if(simflag<6){
 	 if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1.PbPb.244918.ForData.2015.root","corrections.root",kTRUE);}
 	 else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1.PbPb.244918.ForSimulations.2015.root","corrections.root",kTRUE);}
       }
       if(simflag==6){
 	 if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1_plus21.PbPb.244975.ForSimulations.2015.root","corrections.root",kTRUE);}
 	 else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1_plus21.PbPb.244975.ForData.2015.root","corrections.root",kTRUE);}
       }
       if(simflag==7){
 	 if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1_plus21.PbPb.244980.ForSimulations.2015.root","corrections.root",kTRUE);}
 	 else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1_plus21.PbPb.244980.ForData.2015.root","corrections.root",kTRUE);}
       }
       if(simflag==8){
 	 if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1_plus21.PbPb.244982.ForSimulations.2015.root","corrections.root",kTRUE);}
 	 else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1_plus21.PbPb.244982.ForData.2015.root","corrections.root",kTRUE);}
       }
       if(simflag==9){
 	 if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1_plus21.PbPb.244983.ForSimulations.2015.root","corrections.root",kTRUE);}
 	 else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1_plus21.PbPb.244983.ForData.2015.root","corrections.root",kTRUE);}
       }
       if(simflag==10){
 	 if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1_plus21.PbPb.245064.ForSimulations.2015.root","corrections.root",kTRUE);}
 	 else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1_plus21.PbPb.245064.ForData.2015.root","corrections.root",kTRUE);}
       }
       if(simflag==11){
 	 if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1_plus21.PbPb.245066.ForSimulations.2015.root","corrections.root",kTRUE);}
 	 else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1_plus21.PbPb.245066.ForData.2015.root","corrections.root",kTRUE);}
       }
       if(simflag==12){
 	 if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1_plus21.PbPb.245068.ForSimulations.2015.root","corrections.root",kTRUE);}
 	 else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1_plus21.PbPb.245068.ForData.2015.root","corrections.root",kTRUE);}
       }
       if(simflag==13){
 	 if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1_plus21.PbPb.246390.ForSimulations.2015.root","corrections.root",kTRUE);}
 	 else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1_plus21.PbPb.246390.ForData.2015.root","corrections.root",kTRUE);}
       }
       if(simflag==14){
 	 if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1_plus21.PbPb.246391.ForSimulations.2015.root","corrections.root",kTRUE);}
 	 else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1_plus21.PbPb.246391.ForData.2015.root","corrections.root",kTRUE);}
       }
// 	 if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1.PbPb.225717.ForData.2015.root","corrections.root",kTRUE);}
// 	 else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC15k1.PbPb.225717.ForData.2015.root","corrections.root",kTRUE);}
//        }
     }
     if(dataset==2009){//pp 900 GeV
       gSystem->CopyFile("ConfigHadEtMonteCarlopp900GeV.C","ConfigHadEtMonteCarlo.C",kTRUE);
       gSystem->CopyFile("ConfigHadEtReconstructedpp900GeV.C","ConfigHadEtReconstructed.C",kTRUE);
       if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC11b1a.pp.ForData.root","corrections.root",kTRUE);}
       else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC11b1a.pp.ForSimulations.root","corrections.root",kTRUE);}
       pass = 3;
     }
     if(dataset==20111){//pp 2.76 TeV
       if(altV0Scale)gSystem->CopyFile("ConfigHadEtMonteCarlopp276TeVAlt.C","ConfigHadEtMonteCarlo.C",kTRUE);
       else{gSystem->CopyFile("ConfigHadEtMonteCarlopp276TeV.C","ConfigHadEtMonteCarlo.C",kTRUE);}
       gSystem->CopyFile("ConfigHadEtReconstructedpp276TeV.C","ConfigHadEtReconstructed.C",kTRUE);
       if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC11b10a.pp.ForData.root","corrections.root",kTRUE);}
       else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC11b10a.pp.ForSimulations.root","corrections.root",kTRUE);}
       pass = 4;
     }
     if(dataset==2010){//pp 7 TeV
       gSystem->CopyFile("ConfigHadEtMonteCarlopp7TeV.C","ConfigHadEtMonteCarlo.C",kTRUE);
       gSystem->CopyFile("ConfigHadEtReconstructedpp7TeV.C","ConfigHadEtReconstructed.C",kTRUE);
       if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC10e20.pp.ForData.root","corrections.root",kTRUE);}
       else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC10e20.pp.ForSimulations.root","corrections.root",kTRUE);}
       pass = 2;
     }
     if(dataset==2012){//pp 8 TeV
       gSystem->CopyFile("ConfigHadEtMonteCarlopp8TeV.C","ConfigHadEtMonteCarlo.C",kTRUE);
       gSystem->CopyFile("ConfigHadEtReconstructedpp8TeV.C","ConfigHadEtReconstructed.C",kTRUE);
       if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC12c1b.pp.ForData.root","corrections.root",kTRUE);}
       else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC12c1b.pp.ForSimulations.root","corrections.root",kTRUE);}
       pass = 1;
     }
     if(dataset==2013){//pPb 5 TeV
       gSystem->CopyFile("ConfigHadEtMonteCarlopPb5TeV.C","ConfigHadEtMonteCarlo.C",kTRUE);
       gSystem->CopyFile("ConfigHadEtReconstructedpPb5TeV.C","ConfigHadEtReconstructed.C",kTRUE);
       if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC13b3.pPb.ForData.root","corrections.root",kTRUE);}
       else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC13b3.pPb.ForSimulations.root","corrections.root",kTRUE);}
       pass = 3;
     }
   }
   
   
   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  cerr<<"I am adding PID response task 169"<<endl;
  //AliAnalysisTask *AddTaskPIDResponse(Bool_t isMC=kFALSE, Bool_t autoMCesd=kTRUE,
//                                     Bool_t tuneOnData=kFALSE, Int_t recoPass=2,
//                                     Bool_t cachePID=kFALSE, TString detResponse="",
//                                     Bool_t useTPCEtaCorrection = kFALSE);
  AliAnalysisTaskPIDResponse *taskPID=AddTaskPIDResponse(!data,kTRUE,kTRUE,pass);
  //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
  //AddTaskPIDqa();
  if(runCompiledVersion){
    AliAnalysisTaskHadEt *task2 = new AliAnalysisTaskHadEt("TaskHadEtCompiled",!data);
    if(!data) task2->SetMcData();
    //Add thing here to select collision type!!
    if(dataset==2013 || dataset==2015){//pPb 5 TeV
      if(dataset==2015){
	task2->SelectCollisionCandidates(AliVEvent::kINT7 ) ;
	//task2->SelectCollisionCandidates(AliVEvent::kAny ) ;
      }
      else{
	task2->SelectCollisionCandidates(AliVEvent::kINT7 ) ;
      }
    }
    else{   if(dataset!=20100){task2->SelectCollisionCandidates(AliVEvent::kMB ) ;}
    }
    //if(dataset!=20100){task2->SelectCollisionCandidates(AliVEvent::kMB ) ;}
    mgr->AddTask(task2);
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("out2", TList::Class(), AliAnalysisManager::kOutputContainer,"Et.ESD.new.sim.root");
    mgr->ConnectInput(task2,0,cinput1);
    mgr->ConnectOutput(task2,1,coutput2);
    mgr->ConnectInput(task2,0,cinput1);
    mgr->ConnectOutput(task2,1,coutput2);
    
  }
  else{
    AliAnalysisTaskHadEtLocal *task2 = new AliAnalysisTaskHadEtLocal("TaskHadEtLocal",!data);
    if(!data) task2->SetMcData();
    //Add thing here to select collision type!!
    if(dataset==2013 ||dataset==2015){//pPb 5 TeV
      //task2->SelectCollisionCandidates(AliVEvent::kAnyINT ) ;
      task2->SelectCollisionCandidates(AliVEvent::kINT7 ) ;
    }
    else{   if(dataset!=20100){task2->SelectCollisionCandidates(AliVEvent::kMB ) ;}}
    //if(dataset!=20100){task2->SelectCollisionCandidates(AliVEvent::kMB ) ;}
    mgr->AddTask(task2);
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("out2", TList::Class(), AliAnalysisManager::kOutputContainer,"Et.ESD.new.sim.root");
    mgr->ConnectInput(task2,0,cinput1);
    mgr->ConnectOutput(task2,1,coutput2);
    mgr->ConnectInput(task2,0,cinput1);
    mgr->ConnectOutput(task2,1,coutput2);
  }
   
  
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  if(submit){
    mgr->StartAnalysis("grid");
  }
  else{
    mgr->SetDebugLevel(0);
    mgr->StartAnalysis("local",chain);
  }

  timer.Stop();
  timer.Print();
}
