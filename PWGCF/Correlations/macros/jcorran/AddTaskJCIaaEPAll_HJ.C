#if defined(__CLING__)
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include <PWGCF/Correlations/macros/jcorran/AddTaskJCatalyst.C>
#include <PWGCF/Correlations/macros/jcorran/AddTaskJFlowBaseTask.C>
#endif

//_____________________________________________________________________
AliAnalysisTask *AddTaskJCIaaEPAll_HJ(TString taskName="JCIAA", TString CatalystName="JCatalystTaskEP",TString FlowBaseTaskName="JFlowBaseTask", double startEP=0., const int NEPbins=6, Bool_t isMC=kFALSE,  UInt_t flags = 0, Int_t FilterBit = 768 , double eta_min = -0.8, double eta_max = 0.8, double pt_min = 0.0, double pt_max = 100.0, int EPdetID=0, TString cardName="card.input", TString jtrigg="hadron", TString jassoc="hadron", TString cardSetting="", TString inclusFileName=""){
	// Load Custom Configuration and parameters
	// override values with parameters

	cout<<"### DEBUG Input is "<< cardName <<"\t"<<jtrigg<<"\t"<<jassoc<<"\t"<<inclusFileName<<"\t"<<"#########"<<endl;
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

        gROOT->LoadMacro("$ALICE_PHYSICS/PWGCF/Correlations/macros/jcorran/AddTaskJCatalyst.C");
        gROOT->LoadMacro("$ALICE_PHYSICS/PWGCF/Correlations/macros/jcorran/AddTaskJFlowBaseTask_HJ.C");
        AddTaskJCatalyst(CatalystName.Data(),  flags, FilterBit , eta_min, eta_max, pt_min, pt_max, 0);
        AddTaskJFlowBaseTask_HJ(FlowBaseTaskName.Data(), CatalystName.Data(), isMC);

	// 15 degree...double EPbins[NEPbins+1] = {0,15,30,45,60,75,90};
	double EPbins[NEPbins+1];
        for(int i=0;i<=NEPbins;i++) EPbins[i] = startEP + 15.*i;

	//==== Set up di-hadron correlation jT task ====
	AliJCIaaEPTask *myTask[NEPbins];
	for(int i=0;i<NEPbins;i++) {
		TString mytaskName = Form("%s_EP%.0f_%.0f",taskName.Data(),EPbins[i],EPbins[i+1]);
		myTask[i] = new AliJCIaaEPTask(mytaskName.Data(),"JOD");
		myTask[i]->SetDebugLevel(5);
		myTask[i]->SetJFlowBaseTaskName(FlowBaseTaskName.Data());  // AliJCatalystTask has this name hard coded
		myTask[i]->SetEPDector( EPdetID );
		myTask[i]->SetEPmin(EPbins[i]);
		myTask[i]->SetEPmax(EPbins[i+1]);
		cout << myTask[i]->GetName() << endl;
	}

	// === Set up JCard ====
	AliJCard *card = new AliJCard(cardName.Data());
	card->PrintOut();
	card->ReadLine( cardSetting.Data() );
	card->ReCompile();
	card->PrintOut();

	// === Create analysis object ===

	AliJIaaAna *fAna[NEPbins];
	for(int i=0;i<NEPbins;i++) {
		fAna[i] = new AliJIaaAna( kFALSE );
		fAna[i]->SetCard( card );
		fAna[i]->SetTrigger( jtrigg.Data() );
		fAna[i]->SetAssoc( jassoc.Data() );
		fAna[i]->SetEnableEP( true );
		if( inclusFileName ) fAna[i]->SetInclusiveFile(inclusFileName.Data());
	}

	for(int i=0;i<NEPbins;i++) {
		myTask[i]->SetAnalysis( fAna[i] );
		mgr->AddTask((AliAnalysisTask*) myTask[i]);
	}


	// Create containers for input/output
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();


	// Connect input/output
	for(int i=0;i<NEPbins;i++) {
		mgr->ConnectInput(myTask[i], 0, cinput);
		AliAnalysisDataContainer *jHist = mgr->CreateContainer(Form("%scontainer",myTask[i]->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",AliAnalysisManager::GetCommonFileName(), myTask[i]->GetName()));
		mgr->ConnectOutput(myTask[i], 1, jHist );
	}

	return myTask[0];
}

