///*******************************************************
/// Config Description
/// January 22, 2021 - Cristiane Jahnke
/// cristiane.jahnke@cern.ch
/// TPC calibrations for 2017 and 2018 data
///*******************************************************

//isMC,isAOD, period,trigger_index, config, isTender, is_ESparse, is_ESparseTPC, is_EventsEG1, is_EventsEG2, isMultiAnalysis, is_MSparse, is_TPCcalibration

AliAnalysisTask_JPsi_EMCal* Config_JPsi_EMCal(
											
Bool_t isMC=kFALSE, 
Bool_t isAOD = kTRUE,
char * period  = "16l",
Int_t trigger_index=0,
Int_t config=0,
Bool_t isTender,
Bool_t is_ESparse,
Bool_t is_ESparseTPC,
Bool_t is_EventsEG1,
Bool_t is_EventsEG2,
Bool_t isMultiAnalysis,
Bool_t is_MSparse,
Bool_t is_TPCcalibration
                                            
)

{
////____________
///Task config
   // printf("Config loaded properly\n");
	AliAnalysisTask_JPsi_EMCal *task = new AliAnalysisTask_JPsi_EMCal();
	printf("task loaded ------------------------ %p\n ", task);
    
	
	task->SetAODanalysis(isAOD);

	
    
    if(isTender) task->SetUseTender();
    if(isMultiAnalysis) task->SetMultiAnalysis();
    
    if(is_ESparse)task->Set_Fill_ESparse();
    if(is_ESparseTPC)task->Set_Fill_ESparseTPC();
    if(is_MSparse)task->Set_Fill_MSparse();
    
    if(is_TPCcalibration)task->Set_TPCCalibration();
    
    if(is_EventsEG1)task->Set_Select_trigger_events1();
    if(is_EventsEG2)task->Set_Select_trigger_events2();
    
    //event cuts
    task->SetVertexCut(10.0);
	
	
    //to separate trigger threshold
    if(!isMC){
        if(trigger_index==3) task->SetEMCalTriggerEG1();
        if(trigger_index==6) task->SetEMCalTriggerEG1();
        if(trigger_index==4) task->SetEMCalTriggerEG2();
	
        if(trigger_index==7) task->SetEMCalTriggerDG1();
        if(trigger_index==8) task->SetEMCalTriggerDG2();
    
        if(trigger_index==10) task->SetEMCalTriggerEG1DG1();
        if(trigger_index==11) task->SetEMCalTriggerEG2DG2();
    }
//========================================================================================
   //track cuts
    task->SetPtCutMainEle(1.0);
    task->SetRejectKinkMother(kTRUE);
    task->SetTPCandITSrefit(kTRUE);
    task->SetTPCnclsPID(85);//not used inside the task
    
    
    task->SetTPCchi2(4.0);
    task->SetITSchi2(36.0);
    
    
    if(config==1)task->SetEtaCut(-0.8,0.8);
	else task->SetEtaCut(-0.9,0.9);
    
    if(config==2)task->SetPtCutPartner(2.0);
    else task->SetPtCutPartner(1.0);
   
    if(config==3)task->SetITSncls(3);
    else if(config==4)task->SetITSncls(4);
    else task->SetITSncls(2);
    
    if(config==5)task->SetITSpixel(2); //1 kAny, 2 kBoth, 3 kFirst
    else if(config==6)task->SetITSpixel(3); //1 kAny, 2 kBoth, 3 kFirst
    else task->SetITSpixel(1); //1 kAny, 2 kBoth, 3 kFirst
    
   /* if(config==7)task->SetTPCncls(70);
    else if(config==8)task->SetTPCncls(80);
    else if(config==9)task->SetTPCncls(90);
    else if(config==10)task->SetTPCncls(110);*/
    
    task->SetTPCncls(85);//not used anymore
    //now we are using TPCncrossed rows
    task->SetTPCnCrossedRows(70);
    
    if(config==11)task->SetDCACut(0.2,0.4); //xy, z
    else if(config==12)task->SetDCACut(0.5,3.0); //xy, z
    else if(config==13)task->SetDCACut(1.0,4.0); //xy, z
    else if(config==14)task->SetDCACut(1.0,2.0); //xy, z
    else task->SetDCACut(1.0,3.0); //xy, z
//=========================================================================================
    
    //PID cuts
    if(config==15)task->SetTPCnsigmaCut(-3.0,3.0);
    else if(config==16)task->SetTPCnsigmaCut(-2.5,3.0);
    else if(config==17)task->SetTPCnsigmaCut(-1.5,3.0);
    else if(config==18)task->SetTPCnsigmaCut(-1.0,3.0);
    else if(config==19)task->SetTPCnsigmaCut(0,3.0);
    else if(config==20)task->SetTPCnsigmaCut(1.0,3.0);
    
    else if(config==21)task->SetTPCnsigmaCut(-1.5,2.5);
    else if(config==22)task->SetTPCnsigmaCut(-1.5,4.0);
    else task->SetTPCnsigmaCut(-2.25,3.0);
    
	if(config==23)task->SetEoverPCut(0.75,1.3);
    else if(config==24)task->SetEoverPCut(0.85,1.3);
    else if(config==25)task->SetEoverPCut(0.9,1.3);
    else if(config==26)task->SetEoverPCut(0.8,1.4);
    else task->SetEoverPCut(0.8,1.3);
    
    
	task->SetEnergyCut(1);
    
    if(trigger_index==3)task->SetEnergyCut(7);//eg1 16l
    
    if(trigger_index==4 || trigger_index==8 || trigger_index==11){
        
        if(config==27)task->SetEnergyCut(4.5);//eg2
        else if(config==28)task->SetEnergyCut(5.5);//eg2
        else task->SetEnergyCut(5);//eg2
    }
    if(trigger_index==6 || trigger_index==7 || trigger_index==10){
        if(config==29)task->SetEnergyCut(9.5);//eg1
        else if(config==30)task->SetEnergyCut(10.5);//eg1
        else task->SetEnergyCut(10);//eg1
    }
    
   // if(period=="17h2h")task->SetEnergyCut(7);//eg1 16l to test
    
    if(trigger_index==30)task->SetEnergyCut(7);//eg1 16l to test
    if(trigger_index==60)task->SetEnergyCut(10);//eg1 16k to test
    if(trigger_index==40)task->SetEnergyCut(5);//eg1 16l to test
    
    
    //Jpsi analysis (cut online only for efficiencies in MC)
    //For data, cut is applied offline
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
