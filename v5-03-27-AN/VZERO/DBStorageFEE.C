
void DBStorageFEE(){
	
	AliCDBManager *man = AliCDBManager::Instance();
	
	AliCDBStorage *storLoc;
	man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
	
	// Setting gain and pedestal values :
	
	AliVZEROTriggerData *triggerData = new AliVZEROTriggerData();
	
	const int N = AliVZEROTriggerData::kNCIUBoards;
	
	UShort_t Clk1Win1[N] = {1,1,1,1,3,3,3,3};
	UShort_t Clk2Win1[N] = {3,3,3,3,14,14,14,14};
	triggerData->SetClk1Win1(Clk1Win1);
	triggerData->SetClk2Win1(Clk2Win1);
	
	UShort_t Clk1Win2[N] = {7,7,7,7,7,7,7,7};
	UShort_t Clk2Win2[N] = {14,14,14,14,14,14,14,14};
	triggerData->SetClk1Win2(Clk1Win1);
	triggerData->SetClk2Win2(Clk2Win1);
		
	UShort_t delayClk1Win1[N] = {0,0,0,0,300,300,300,300}; // 1 unit correspond to 10 ps delay
	UShort_t delayClk2Win1[N] = {0,0,0,0,200,200,200,200}; // 1 unit correspond to 10 ps delay
	triggerData->SetDelayClk1Win1(delayClk1Win1);
	triggerData->SetDelayClk2Win1(delayClk2Win1);
	
	UShort_t delayClk1Win2[N] = {100,100,100,100,100,100,100,100}; // 1 unit correspond to 10 ps delay
	UShort_t delayClk2Win2[N] = {200,200,200,200,200,200,200,200}; // 1 unit correspond to 10 ps delay
	triggerData->SetDelayClk1Win2(delayClk1Win2);
	triggerData->SetDelayClk2Win2(delayClk2Win2);
	
	UShort_t LatchWin1[N] = {16,16,16,16,16,16,16,16};
	triggerData->SetLatchWin1(LatchWin1);
	
	UShort_t LatchWin2[N] = {16,16,16,16,16,16,16,16};
	triggerData->SetLatchWin2(LatchWin2);
	
	UShort_t ResetWin1[N] = {16,16,16,16,16,16,16,16};
	triggerData->SetResetWin1(ResetWin1);
	
	UShort_t ResetWin2[N] = {16,16,16,16,16,16,16,16};
	triggerData->SetResetWin2(ResetWin2);

	Bool_t PedestalSubtraction[N] = {1,1,1,1,1,1,1,1};
	triggerData->SetPedestalSubtraction(PedestalSubtraction);
	
	triggerData->SetBBAThreshold(1);
	triggerData->SetBBCThreshold(1);

	triggerData->SetBGAThreshold(1);
	triggerData->SetBGCThreshold(1);
	
	triggerData->SetBBAForBGThreshold(1);
	triggerData->SetBBCForBGThreshold(1);
	
	triggerData->SetCentralityV0AThrLow(100);
	triggerData->SetCentralityV0AThrHigh(500);
	
	triggerData->SetCentralityV0CThrLow(100);
	triggerData->SetCentralityV0CThrHigh(500);
	
	triggerData->SetMultV0AThrLow(2);
	triggerData->SetMultV0AThrHigh(10);
	
	triggerData->SetMultV0CThrLow(2);
	triggerData->SetMultV0CThrHigh(10);
	
	for(int ibrd =0;ibrd<8;ibrd++){
		for(int ich =0;ich<8;ich++){
			triggerData->SetEnableTiming(kTRUE,ibrd,ich);
			triggerData->SetEnableCharge(kTRUE,ibrd,ich);
			triggerData->SetDelayHit(0,ibrd,ich);
			for(int iint=0;iint<2;iint++){
				triggerData->SetPedestal(16,iint,ibrd,ich);
				triggerData->SetPedestalCut(17,iint,ibrd,ich);
			}
		}
	}
	
	for(int i =0;i<5;i++) triggerData->SetTriggerSelected(i, i);
	
	// Creation of the object VZERO Trigger Configuration as a MetaData
	
	TObjString str("VZERO Trigger Configuration");      // object that will be stored
	
	AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
	
	AliCDBId id("VZERO/Trigger/Data",0,AliCDBRunRange::Infinity());
	
	//md->SetObjectClassName("VZERO calibration parameters"); automatically 
	//set to AliVZEROCalibData by the CDB classes during storage 
	md->SetResponsible("Brigitte Cheynis");
	md->SetBeamPeriod(0);
	md->SetAliRootVersion("v4-18-Release");
	md->SetComment("Prototype");
	md->PrintMetaData();
	
	storLoc = man->GetDefaultStorage();
	storLoc->Put(triggerData, id, md);
	
	storLoc->Delete();
	delete md;
	
}
