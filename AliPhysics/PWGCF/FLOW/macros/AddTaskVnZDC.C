#include<TList.h>

void AddTaskVnZDC(Double_t dcentrMin, Double_t dcentrMax, Double_t dptMin = 0.2, Double_t dptMax = 10.0,
		  Double_t detaMin = -0.8, Double_t detaMax = 0.8, Int_t iCharge = 1,
		  TString sanalysisTypeUser = "AOD", Int_t iAODfilterBit = 768, TString sDataSet = "2010",
                  TString sAnalysisDef = "recenter1", TString sanalysisType = "AUTOMATIC", TString sEvTrigger = "MB",
                  Bool_t bEventCutsQA = kFALSE, Bool_t bTrackCutsQA = kFALSE, Bool_t bUseVZERO = kTRUE,
                  Bool_t bPileUp = kFALSE, Bool_t bPileUpTight = kFALSE, Bool_t useParFiles = kFALSE,
                  Bool_t recent = kFALSE,Bool_t SetFBEffi = kFALSE,
		  TString ZDCRecenterFile1="alien:///alice/cern.ch/user/m/mhaque/calib_files/recenter1_zdc_ver1.root",
		  TString EfficiencyFB768="alien:///alice/cern.ch/user/m/mhaque/calib_files/FB_hijing_eff2010.root",
                  TString sCentrEstimator = "V0", Double_t dVertexRange = 10., Double_t dDCAxy = 2.4, Double_t dDCAz = 3.2,
                  Double_t dMinClusTPC = 70, const char *suffix = "")
{
	gSystem->Load("libGeom");
	gSystem->Load("libVMC");
	gSystem->Load("libXMLIO");
	gSystem->Load("libPhysics");
	gSystem->Load("libCore.so");
	gSystem->Load("libTree.so");
	gSystem->Load("libSTEERBase.so");
	gSystem->Load("libESD.so");
	gSystem->Load("libAOD.so");
	gSystem->Load("libANALYSIS.so");
	gSystem->Load("libANALYSISalice.so");
	gSystem->Load("libOADB.so");
	gSystem->Load("libPWGflowBase.so");

	if(!useParFiles)
        gSystem->Load("libPWGflowTasks.so");

	//rihan:is this my local AliPhysics or from Cern Grid?
	gSystem->AddIncludePath("-I. -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include -g ");
	gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/OCDB -I$ALICE_ROOT/STEER/macros -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/TRD -I$ALICE_ROOT/ZDC -I$ALICE_ROOT/macros -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/OADB $ALICE_PHYSICS/OADB/macros -I$ALICE_PHYSICS/PWGGA -I$ALICE_PHYSICS/PWGCF -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/TENDER -I$ALICE_PHYSICS/TENDER/Tender -I$ALICE_PHYSICS/TENDER/TenderSupplies -I$ALICE_PHYSICS/PARfiles -I$ALICE_PHYSICS/PWGCF/FLOW/macros I$ALICE_PHYSICS/PWGPP/ZDC -g ");

	//the manager is static, so get the existing manager via the static method
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

	if(!mgr){
	  printf("\n!! No Manager !! \n check runGrid.C macro --> Rihan \n");
	  return;
	}
	if(!mgr->GetInputEventHandler()){
	   printf("\n !!No input event handler!!\n check runGrid.C \n");
	   return;
	}


	//gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"); // Needed for LHC2015o
	//AliMultSelectionTask *MultSelection = AddTaskMultSelection(kFALSE);            // kFALSE == User mode, kTRUE == Calibration mode
        //MultSelection->SetSelectedTriggerClass(AliVEvent::kINT7);
        //MultSelection->SetSelectedTriggerClass(AliVEvent::kMB);
        //mgr->AddTask(MultSelection);


	TString taskFEname = "FlowEventTaskCh";                // charged particle (not proton)
	Bool_t bCutsQA = (Bool_t)(bEventCutsQA || bTrackCutsQA); // if any QA is 1, bCutsQA = 1

	// create instance of the class: because possible qa plots are added in a second output slot,
	AliAnalysisTaskFlowEvent *taskFE_charge = new AliAnalysisTaskFlowEvent(taskFEname,"",bCutsQA);

	taskFE_charge->SetQAOn(bCutsQA);
	taskFE_charge->SetAnalysisType(sanalysisType); //sanalysisType = AUTOMATIC see the initializers!!

	// add the task to the manager. Added later.
	// mgr->AddTask(taskFE_charge);



	// set the trigger selection
	if(sEvTrigger == "Cen"){
		taskFE_charge->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
	}
	else if(sEvTrigger == "SemiCen"){
		taskFE_charge->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kSemiCentral);
	}
	else if(sEvTrigger == "MB"){
	       if(sDataSet == "2010"){
		taskFE_charge->SelectCollisionCandidates(AliVEvent::kMB);
		}
	       else if(sDataSet == "2011"){
		taskFE_charge->SelectCollisionCandidates(AliVEvent::kMB);
		}
               else{
		taskFE_charge->SelectCollisionCandidates(AliVEvent::kINT7);
		}
	}
	else{
	       printf("Rihan, you did not set any Event trigger.!! \n check AddTaskVnZDC.C\n");
	       return;
          }





	//--------------- define the event cuts object ---------------------
	 AliFlowEventCuts *cutsEvent = new AliFlowEventCuts("EventCuts");
         cutsEvent->SetCheckPileup(kFALSE);
	 cutsEvent->SetPrimaryVertexZrange(-dVertexRange, dVertexRange);      // vertex-z cut
         cutsEvent->SetQA(bEventCutsQA);                                      // enable the qa plots
         cutsEvent->SetCutTPCmultiplicityOutliersAOD(kTRUE); 	              // multiplicity outlier cut

	if(sDataSet=="2015")
	{
	 cutsEvent->SetCentralityPercentileRange(dcentrMin, dcentrMax, kTRUE);
	}
        else
	{
	 cutsEvent->SetCentralityPercentileRange(dcentrMin, dcentrMax);
	 cutsEvent->SetCheckPileup(kTRUE);
	}

      //method used for centrality determination
       if(sCentrEstimator=="V0"||sCentrEstimator=="V0M")
          cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kV0);

	if(sCentrEstimator=="TPC")
	  cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kTPConly);

	if(sDataSet == "2011"){
	  cutsEvent->SetLHC11h(kTRUE);
	 }
	  else if(sDataSet == "2010"){
          cutsEvent->SetLHC10h(kTRUE);
	 }


      
  AliFlowTrackCuts* RefMultCuts = new AliFlowTrackCuts("RefMultCuts");
  RefMultCuts->SetParamType(AliFlowTrackCuts::kAODFilterBit);
  RefMultCuts->SetAODfilterBit(iAODfilterBit);
  RefMultCuts->SetMinimalTPCdedx(-999999999);
  RefMultCuts->SetMaxDCAToVertexXY(dDCAxy);
  RefMultCuts->SetMaxDCAToVertexZ(dDCAz);
  RefMultCuts->SetMinNClustersTPC(dMinClusTPC);
  RefMultCuts->SetMinChi2PerClusterTPC(0.1);
  RefMultCuts->SetMaxChi2PerClusterTPC(4.);
  RefMultCuts->SetPtRange(dptMin,dptMax);
  RefMultCuts->SetEtaRange(detaMin,detaMax);
  RefMultCuts->SetAcceptKinkDaughters(kFALSE);


        cutsEvent->SetRefMultCuts(RefMultCuts);
        cutsEvent->SetRefMultMethod(AliFlowEventCuts::kTPConly);




	taskFE_charge->SetCutsEvent(cutsEvent);	//pass these cuts to your flow event task
        //-------------------------- --------------------------------







	//Define Track Cuts for RP and POI
	AliFlowTrackCuts *cutsRP  =  new AliFlowTrackCuts("RP cuts");
	cutsRP->SetParamType(AliFlowTrackCuts::kAODFilterBit);//sets how we want to select the tracks.
	cutsRP->SetAODfilterBit(iAODfilterBit);
	cutsRP->SetEtaGap(-0.0,0.0);                          //rihan:check this gap **************
	cutsRP->SetPhiMin(0.);
	cutsRP->SetPhiMax(TMath::TwoPi());
      //cutsRP->SetVZEROgainEqualizationPerRing(kTRUE);       //options for the reweighting VZERO
      //cutsRP->SetApplyRecentering(kTRUE);                   //check flow manual for Qcum
	cutsRP->SetDivSigma(kTRUE);
	cutsRP->SetQA(bTrackCutsQA);
        cutsRP->SetMinimalTPCdedx(-99999);
	cutsRP->SetMinNClustersTPC(dMinClusTPC);
	cutsRP->SetPtRange(dptMin,dptMax);
	cutsRP->SetEtaRange(detaMin,detaMax);
      //cutsRP->SetCharge(iCharge);
      //cutsRP->SetMaxDCAToVertexXY(dDCAxy);
      //cutsRP->SetMaxDCAToVertexZ(dDCAz);
	cutsRP->SetRequireCharge(kTRUE);
        cutsRP->SetMinChi2PerClusterTPC(0.1);
	cutsRP->SetMaxChi2PerClusterTPC(4.0);
	cutsRP->SetAcceptKinkDaughters(kFALSE);
      //cutsRP->SetMaxDCAToVertexXY(dDCAxy);
      //cutsRP->SetMaxDCAToVertexZ(dDCAz);

	AliFlowTrackCuts *cutsPOI_charge   = new AliFlowTrackCuts("POI cuts Prot");
	cutsPOI_charge->SetParamType(AliFlowTrackCuts::kAODFilterBit);
	cutsPOI_charge->SetAODfilterBit(iAODfilterBit);
	cutsPOI_charge->SetMinimalTPCdedx(-99999);
	cutsPOI_charge->SetMinNClustersTPC(dMinClusTPC);
	cutsPOI_charge->SetPtRange(dptMin,dptMax);
	cutsPOI_charge->SetEtaRange(detaMin, detaMax);
      //cutsPOI_charge->SetCharge(iCharge);
      //cutsPOI_charge->SetMaxDCAToVertexXY(dDCAxy);
      //cutsPOI_charge->SetMaxDCAToVertexZ(dDCAz);
	cutsPOI_charge->SetRequireCharge(kTRUE);
	cutsPOI_charge->SetQA(bTrackCutsQA);
        cutsPOI_charge->SetMinChi2PerClusterTPC(0.1);
	cutsPOI_charge->SetMaxChi2PerClusterTPC(4.0);
	cutsPOI_charge->SetAcceptKinkDaughters(kFALSE);
      //cutsPOI_charge->SetMaxDCAToVertexXY(dDCAxy);
      //cutsPOI_charge->SetMaxDCAToVertexZ(dDCAz);
      //cutsPOI_charge->SetRequireStrictTOFTPCagreement(kTRUE);


	taskFE_charge->SetCutsRP(cutsRP);
	taskFE_charge->SetCutsPOI(cutsPOI_charge);
	taskFE_charge->SetSubeventEtaRange(-10.0,-2.0, 2.0,10.0);


	// Finally add the task to the manager
	mgr->AddTask(taskFE_charge);



	//get the default name of the output file ("AnalysisResults.root")
	TString file = AliAnalysisManager::GetCommonFileName();

	// Create a data container for the output of the 'FlowEvent' task. The output of the task (AliFlowEventSimple object)
	// will be passed to the 'AliFlowAnalysis' tasks. Note that we use a "kExchangeContainer here", which exchanges data
	// between classes of the analysis chain, but is not written to the output file.

	AliAnalysisDataContainer *coutputFE = mgr->CreateContainer("FECont",AliFlowEventSimple::Class(),AliAnalysisManager::kExchangeContainer);
	AliAnalysisDataContainer    *cinput = mgr->GetCommonInputContainer();  //AOD event

	mgr->ConnectInput( taskFE_charge, 0, cinput); 	        //connect the input data (AOD ?) to the flow event task
	mgr->ConnectOutput(taskFE_charge, 1, coutputFE); 	//get the output of taskFE_charge to a exchange container.

	// create an additional container for the QA output of the flow event task
	// the QA histograms will be stored in a sub-folder of the output file called 'QA'

	TString taskFEQA = file;      // file is the common outfile filename
	taskFEQA   += ":QAcharge";

	AliAnalysisDataContainer *coutputFEQA = mgr->CreateContainer("FEContQA",TList::Class(),AliAnalysisManager::kOutputContainer,taskFEQA.Data());
	mgr->ConnectOutput(taskFE_charge, 2, coutputFEQA);          // kOutputContainer: written to the output file


       







	//================= Here is the actual TaskMacro object ===============

	AliAnalysisTaskVnZDC *taskQC_prot = new AliAnalysisTaskVnZDC("TaskQC_charge", kFALSE);

	if(sEvTrigger == "Cen"){
	          taskQC_prot->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
		}
	else if(sEvTrigger == "SemiCen"){
		  taskQC_prot->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kSemiCentral);
		}
	 else if(sEvTrigger == "MB")
	{
	       if(sDataSet == "2010")
		 {
		  taskQC_prot->SelectCollisionCandidates(AliVEvent::kMB);
		 }
	       if(sDataSet == "2011")
		 {
		  taskQC_prot->SelectCollisionCandidates(AliVEvent::kMB);
		 }
	          else if(sDataSet == "2015")
		 {
		  taskQC_prot->SelectCollisionCandidates(AliVEvent::kINT7);
		 }
	  }
          else if(sEvTrigger == "Any")
	  {
	          taskQC_prot->SelectCollisionCandidates(AliVEvent::kAny);
          }



	    taskQC_prot->SetHarmonic(2);
	  //taskQC_prot->SetRunFlag(77);              //this is set automatically
	    taskQC_prot->SetUsePhiWeights(kFALSE);    //Rihan
	    taskQC_prot->SetDataSet(sDataSet);        //Rihan
	    taskQC_prot->SetAnalysisSet(sAnalysisDef); //Rihan
	    taskQC_prot->SetRejectPileUp(bPileUp);     //Rihan
	    taskQC_prot->SetRejectPileUpTight(bPileUpTight); //kTRUE:700,kFALSE:15000
            //(if 'SetRejectPileUp' is 'kFALSE' then 'SetRejectPileUpTight' does not apply)

	    //read and pass the efficiency histograms:



      if(SetFBEffi){
          TFile* FBefficiency = TFile::Open(EfficiencyFB768,"READ");
	  //std::cout<<"\n Info: FB efficiency requested... \n"<<std::endl;
          if(!FBefficiency) {
	   std::cout<<"\n ERROR: FB efficiency file not found! \n"<<std::endl;
           exit(1);
         } 
	 /* TList* FBEffiList = new TList();
          FBEffiList->SetOwner(kTRUE);
	  TH1D *h[10];
	  for(int i=0;i<10;i++){
	    h[i] = (TH1D *) FBefficiency->Get(Form("eff_unbiased_%d",i));
            FBEffiList->Add(h[i]);
	  }*/

       TList* FBEffiList = dynamic_cast<TList*>(FBefficiency->FindObjectAny("fMcEffiHij"));
       TList* FBEffiListUse = FBEffiList->Clone();

       //FBefficiency->Close();
  
       if(FBEffiListUse) {
       //std::cout<<"\n\n Efficiency Histograms found\n:"<<FBEffiListUse->ls()<<"\n\n"<<std::endl;
         taskQC_prot->SetFBEffiList(FBEffiListUse);
        }
         else{
	  std::cout<<"\n\n !!!!**** ERROR: FB Efficiency Histograms not found !!\n check name in AddTaskMacro ?? ****!!!\n"<<std::endl;
          exit(1);
        }
      }



      //sanity:
      if(sAnalysisDef=="recenter2"){
	recent=kTRUE;
        SetFBEffi=kTRUE;
      }
      if(sAnalysisDef=="recenter1"){
	recent=kFALSE;
        SetFBEffi=kFALSE;
      }

     if(recent) // use ZDC recentering
      {
       TFile* ZDCESEFile1 = TFile::Open(ZDCRecenterFile1,"READ");
       //std::cout<<"\n Info: ZDC recenter histograms requested... \n"<<std::endl;
         if(!ZDCESEFile1) {
	  std::cout<<"\n ERROR: ZDC recenter file1 not found! \n"<<std::endl;
          exit(1);
         }

       TList* ZDCESEList  = dynamic_cast<TList*>(ZDCESEFile1->FindObjectAny("recenterZDC"));
 
       const TList* ZDCESEListUse = new TList();
       ZDCESEListUse = (TList* ) ZDCESEList->Clone();
       ZDCESEFile1->Close();

       if(ZDCESEList) {
	 //std::cout<<"\n@@@@@@@@@\n Recenter 1 Histograms found ("<<ZDCESEList->ls()<<")\n\n"<<std::endl;
         //taskQC_prot->SetZDCESEList(ZDCESEList);
         taskQC_prot->SetZDCESEList(ZDCESEListUse);
         //std::cout<<"\n@@@@@@@@@\n Recenter 1 Histograms found ("<<ZDCESEListUse->ls()<<")\n\n"<<std::endl;
       }
        else{
	 std::cout<<"\n\n !!!!**** ERROR: Root file is there but TList not found !!\n check name in AddTaskMacro ?? ****!!!\n"<<std::endl;
         exit(1);
        }

       /*
        TFile* ZDCESEFile2 = TFile::Open(ZDCRecenterFile2,"READ");
         if(!ZDCESEFile2) {
	  std::cout<<"\n ERROR: ZDC recenter file2 not found! \n"<<std::endl;
          exit(1);
         }
       TList* ZDCESEList2 = dynamic_cast<TList*>(ZDCESEFile2->FindObjectAny("recenterZDC"));
       TList* ZDCESEListUse2 = ZDCESEList2->Clone();

       ZDCESEFile2->Close();

       if(ZDCESEListUse2) {
         taskQC_prot->SetZDCESEList2(ZDCESEListUse2);
	 //std::cout<<"\n@@@@@@@@@\n Recenter 2 Histograms found ("<<ZDCESEListUse->ls()<<")\n\n"<<std::endl;
       }
        else{
	 std::cout<<"\n\n !!!!**** ERROR: Root file is there but TList not found !!\n check name in AddTaskMacro ?? ****!!!\n"<<std::endl;
         exit(1);
        }*/

      }





      mgr->AddTask(taskQC_prot);                      // connect the task to the analysis manager
      mgr->ConnectInput(taskQC_prot, 0, cinput);      // AOD event.!!
      mgr->ConnectInput(taskQC_prot, 1, coutputFE);   // connect the output of the flow event task to the flow analysis task


      TString outputSP = file;      // file is the common outfile filename
              outputSP += ":ZDCgains";

      AliAnalysisDataContainer *coutputSP = mgr->CreateContainer("QAotherZDC",TList::Class(),AliAnalysisManager::kOutputContainer,outputSP.Data());
      mgr->ConnectOutput(taskQC_prot, 1, coutputSP);

      AliAnalysisDataContainer *coutputSP = mgr->CreateContainer("recenterZDC",TList::Class(),AliAnalysisManager::kOutputContainer,outputSP.Data());
      mgr->ConnectOutput(taskQC_prot, 2, coutputSP);

	if(!mgr->InitAnalysis())  // check if we can initialize the manager
        {
	  printf("Rihan,your analysis manager is not set correctly.\n Check AddTaskMacro.C");
	  return;
	}

      printf("\n ... AddTaskVnZDC called succesfully ....  \n");

	return;

}//----------- main ends ---------------
