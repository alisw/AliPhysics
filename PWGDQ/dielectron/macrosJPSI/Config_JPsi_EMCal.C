///*******************************************************
///Config Description
/// August 13th 2018 - Cristiane Jahnke
///*******************************************************

AliAnalysisTask_JPsi_EMCal* Config_JPsi_EMCal(
											
Bool_t isMC=kFALSE, 
Bool_t isAOD = kTRUE,
char * period  = "16l",
Int_t trigger_index=0,
Int_t config=0,
Bool_t isTender,
Bool_t is_ESparse,
Bool_t is_MSparse
                                              
)

{
////____________
///Task config
   // printf("Config loaded properly\n");
	AliAnalysisTask_JPsi_EMCal *task = new AliAnalysisTask_JPsi_EMCal();
	printf("task loaded ------------------------ %p\n ", task);
    
	
	task->SetAODanalysis(isAOD);

	if(period == "11d")task->SetPeriod2011();
    
    if(isTender) task->SetUseTender();
    
    if(is_ESparse)task->Set_Fill_ESparse();
    if(is_MSparse)task->Set_Fill_MSparse();
    
    
    //event cuts
    task->SetVertexCut(10.0);
	
	
    //to separate trigger threshold
	if(trigger_index==3) task->SetEMCalTriggerEG1();
    if(trigger_index==6) task->SetEMCalTriggerEG1();
	if(trigger_index==4) task->SetEMCalTriggerEG2();
	
	if(trigger_index==7) task->SetEMCalTriggerDG1();
	if(trigger_index==8) task->SetEMCalTriggerDG2();

   //track cuts
	task->SetEtaCut(-0.9,0.9);
    task->SetPtCutMainEle(1.0);
    task->SetPtCutPartner(1.0);
    task->SetRejectKinkMother(kTRUE);
    task->SetTPCandITSrefit(kTRUE);
    task->SetITSncls(2);
    task->SetITSpixel(1); //1 kAny, 2 kBoth, 3 kFirst
    task->SetTPCncls(85);
    task->SetTPCnclsPID(85);
    task->SetTPCchi2(4.0);
    task->SetDCACut(1.0,3.0); //xy, z
    
    
    //PID cuts
    task->SetTPCnsigmaCut(-2.25,3.0);
	task->SetEoverPCut(0.8,1.3);
	task->SetEnergyCut(1);
    
    if(trigger_index==3)task->SetEnergyCut(7);//eg1 16l
    if(trigger_index==4 || trigger_index==8)task->SetEnergyCut(5);//eg2
	if(trigger_index==6 || trigger_index==7)task->SetEnergyCut(10);//eg1
    
   // if(period=="17h2h")task->SetEnergyCut(7);//eg1 16l to test
    
    if(trigger_index==30)task->SetEnergyCut(7);//eg1 16l to test
    if(trigger_index==60)task->SetEnergyCut(11);//eg1 16k to test
    if(trigger_index==40)task->SetEnergyCut(5);//eg1 16l to test
    
    
    //Jpsi analysis
    task->SetMassCut(2.92, 3.16);
    
    
//______________________________________
///Created by the user
    
	
//______________________________________
///Particle identification
	//AliHFEpid *pid = task->GetPID();

//______________________________________
//In the case of a simulation
	if(isMC)
	{
	 // pid->SetHasMCData(kTRUE);
	  task->SetMCanalysis();
	}
//______________________________________

//______________________________________________________
//Configure PID
	//_________________________
	//TPC PID
	//pid->AddDetector("TPC", 1);				//Add TPC PID
	
	//_________________________
	//Configure TPC cut
	//Defaul = -1 to 3 sigmas
	//Note that it is also possible to define a model instead of a constant
	//--------->For this change the "cut model"
	/*
	Double_t params[4];
	char *cutmodel;
	cutmodel = "pol0";
	
    params[0] = -2.25;
	
    Double_t max=3.0;
	
	pid->ConfigureTPCdefaultCut(cutmodel,params,max);

	*/
    
///_______________________________________________________________________________________________________________

	//printf("*************************************\n");
	//printf("Configuring standard Task:\n");
	//pid->PrintStatus();
	//printf("*************************************\n");

	return task;
}
