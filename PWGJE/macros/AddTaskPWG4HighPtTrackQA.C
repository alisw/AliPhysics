void AddTaskPWG4HighPtTrackQA(TString year = "2010", TString prodType = "LHC10h",Bool_t isPbPb = kTRUE, Int_t iAODanalysis = 0, Bool_t bReduced = kTRUE, Int_t filterBit = 272)
{

  if(iAODanalysis==0) { //run on ESDs
    if(year.Contains("2010")) {
      if(bReduced)
	AddTaskPWG4HighPtTrackQAAllReduced(prodType.Data(),isPbPb,iAODanalysis);
      else
	AddTaskPWG4HighPtTrackQAAll(prodType.Data(),isPbPb,iAODanalysis);
    }
    else if(year.Contains("2011")) {
      if(bReduced)
	AddTaskPWG4HighPtTrackQAAllReduced2011(prodType.Data(),isPbPb,iAODanalysis);
      else
	AddTaskPWG4HighPtTrackQAAll2011(prodType.Data(),isPbPb,iAODanalysis);
    }
    else if(year.Contains("2013")) {
      AddTaskPWG4HighPtTrackQApPb();
    }

  }
  else if(iAODanalysis==1) {
    AddTaskPWG4HighPtTrackQAAOD(prodType.Data(),isPbPb,iAODanalysis,filterBit); 
  }

}

void AddTaskPWG4HighPtTrackQApPb(char *prodType = "LHC13b") {

  AliPWG4HighPtTrackQA *taskTrackQA05cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,kFALSE,0,10,0,5,AliVEvent::kINT7);
  AliPWG4HighPtTrackQA *taskTrackQA75cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,kFALSE,0,10,7,5,AliVEvent::kINT7);
  
  if(!strcmp(prodType,"LHC13d") || !strcmp(prodType,"LHC13e") || !strcmp(prodType,"LHC13f")) {
    AliPWG4HighPtTrackQA *taskTrackQA05cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,kFALSE,0,10,0,5,AliVEvent::kEMCEJE);
    AliPWG4HighPtTrackQA *taskTrackQA75cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,kFALSE,0,10,7,5,AliVEvent::kEMCEJE);
  }

}


void AddTaskPWG4HighPtTrackQAAll(char *prodType = "LHC10h",Bool_t isPbPb=kTRUE, Int_t iAODanalysis = 0) 
{    

  int cent = 10;
  
  AliPWG4HighPtTrackQA *taskTrackQA00cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0);
  AliPWG4HighPtTrackQA *taskTrackQA01cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,1);
  //  AliPWG4HighPtTrackQA *taskTrackQA02cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,2);
  // AliPWG4HighPtTrackQA *taskTrackQA10cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,1,0);
  // AliPWG4HighPtTrackQA *taskTrackQA11cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,1,1);
  //  AliPWG4HighPtTrackQA *taskTrackQA20cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,2,0);
  //  AliPWG4HighPtTrackQA *taskTrackQA21cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,2,1);
  //  AliPWG4HighPtTrackQA *taskTrackQA40cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,4,0);
  //  AliPWG4HighPtTrackQA *taskTrackQA41cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,4,1);
  //  AliPWG4HighPtTrackQA *taskTrackQA50cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,5,0);
  //  AliPWG4HighPtTrackQA *taskTrackQA60cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,6,0);
  AliPWG4HighPtTrackQA *taskTrackQA70cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,0);
  AliPWG4HighPtTrackQA *taskTrackQA71cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,1);
  AliPWG4HighPtTrackQA *taskTrackQA72cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,2);

  if(isPbPb) {
    for(cent=0; cent<4; cent++) {
      AliPWG4HighPtTrackQA *taskTrackQA00 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0);
      AliPWG4HighPtTrackQA *taskTrackQA01 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,1);
      //    AliPWG4HighPtTrackQA *taskTrackQA02 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,2);
      // AliPWG4HighPtTrackQA *taskTrackQA10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,1,0);
      // AliPWG4HighPtTrackQA *taskTrackQA11 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,1,1);
      //      AliPWG4HighPtTrackQA *taskTrackQA20 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,2,0);
      //      AliPWG4HighPtTrackQA *taskTrackQA21 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,2,1);
      //      AliPWG4HighPtTrackQA *taskTrackQA40 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,4,0);
      //      AliPWG4HighPtTrackQA *taskTrackQA41 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,4,1);
      //      AliPWG4HighPtTrackQA *taskTrackQA50 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,5,0);
      //      AliPWG4HighPtTrackQA *taskTrackQA60 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,6,0);
      AliPWG4HighPtTrackQA *taskTrackQA70 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,0);
      AliPWG4HighPtTrackQA *taskTrackQA71 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,1);
      AliPWG4HighPtTrackQA *taskTrackQA72 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,2);
    }
  }

}

void AddTaskPWG4HighPtTrackQAAll2011(char *prodType = "LHC11h",Bool_t isPbPb=kTRUE, Int_t iAODanalysis = 0) 
{    

  Int_t cent = 10;
  UInt_t iPhysicsSelectionFlag = AliVEvent::kMB;
  UInt_t iPhysicsSelectionFlagCentral = AliVEvent::kCentral;
  UInt_t iPhysicsSelectionFlagSemiCentral = AliVEvent::kSemiCentral;

  
  AliPWG4HighPtTrackQA *taskTrackQA00cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0,iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA01cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,1,iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA70cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,0,iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA71cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,1,iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA72cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,2,iPhysicsSelectionFlag);
  
  AliPWG4HighPtTrackQA *taskTrackQA05cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,5,iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA74cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,4,iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA75cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,5,iPhysicsSelectionFlag);

  if(isPbPb) {
    for(cent=0; cent<4; cent++) {
      AliPWG4HighPtTrackQA *taskTrackQA00 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0,iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA01 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,1,iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA70 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,0,iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA71 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,1,iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA72 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,2,iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA05 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,5,iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA74 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,4,iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA75 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,5,iPhysicsSelectionFlag);

      if(cent==0) {
	AliPWG4HighPtTrackQA *taskTrackQA00 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0,iPhysicsSelectionFlagCentral);
	AliPWG4HighPtTrackQA *taskTrackQA01 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,1,iPhysicsSelectionFlagCentral);
	AliPWG4HighPtTrackQA *taskTrackQA70 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,0,iPhysicsSelectionFlagCentral);
	AliPWG4HighPtTrackQA *taskTrackQA71 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,1,iPhysicsSelectionFlagCentral);
	AliPWG4HighPtTrackQA *taskTrackQA72 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,2,iPhysicsSelectionFlagCentral);
	AliPWG4HighPtTrackQA *taskTrackQA05 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,5,iPhysicsSelectionFlagCentral);
	AliPWG4HighPtTrackQA *taskTrackQA74 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,4,iPhysicsSelectionFlagCentral);
	AliPWG4HighPtTrackQA *taskTrackQA75 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,5,iPhysicsSelectionFlagCentral);
      }
      else {
	AliPWG4HighPtTrackQA *taskTrackQA00 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0,iPhysicsSelectionFlagSemiCentral);
	AliPWG4HighPtTrackQA *taskTrackQA01 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,1,iPhysicsSelectionFlagSemiCentral);
	AliPWG4HighPtTrackQA *taskTrackQA70 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,0,iPhysicsSelectionFlagSemiCentral);
	AliPWG4HighPtTrackQA *taskTrackQA71 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,1,iPhysicsSelectionFlagSemiCentral);
	AliPWG4HighPtTrackQA *taskTrackQA72 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,2,iPhysicsSelectionFlagSemiCentral);
	AliPWG4HighPtTrackQA *taskTrackQA05 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,5,iPhysicsSelectionFlagSemiCentral);
	AliPWG4HighPtTrackQA *taskTrackQA74 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,4,iPhysicsSelectionFlagSemiCentral);
	AliPWG4HighPtTrackQA *taskTrackQA75 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,5,iPhysicsSelectionFlagSemiCentral);
      }

    }
  }

}

void AddTaskPWG4HighPtTrackQAAllReduced(char *prodType = "LHC11h",Bool_t isPbPb=kTRUE, Int_t iAODanalysis = 0) 
{    

  int cent = 10;
  
  if(isPbPb) {
    for(cent=0; cent<4; cent++) {
      AliPWG4HighPtTrackQA *taskTrackQA00 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0);
      AliPWG4HighPtTrackQA *taskTrackQA01 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,1);
      AliPWG4HighPtTrackQA *taskTrackQA70 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,0);
      AliPWG4HighPtTrackQA *taskTrackQA71 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,1);
      AliPWG4HighPtTrackQA *taskTrackQA72 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,2);
    }
  }

}


void AddTaskPWG4HighPtTrackQALHC11hLTS(char *prodType = "LHC10h",Bool_t isPbPb=kTRUE, Int_t iAODanalysis = 0) 
{    

  Int_t cent = 10;
  UInt_t iPhysicsSelectionFlag = AliVEvent::kMB|AliVEvent::kCentral|AliVEvent::kSemiCentral;
  UInt_t iPhysicsSelectionFlagEMCEJE = AliVEvent::kEMCEJE;
  
  AliPWG4HighPtTrackQA *taskTrackQA00cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0,iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA01cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,1,iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA70cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,0,iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA71cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,1,iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA72cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,2,iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA05cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,5,iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA74cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,4,iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA75cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,5,iPhysicsSelectionFlag);
  AliPWG4HighPtTrackQA *taskTrackQA40cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,4,0,iPhysicsSelectionFlag);

  AliPWG4HighPtTrackQA *taskTrackQAEMCJE00cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0,iPhysicsSelectionFlagEMCEJE);
  AliPWG4HighPtTrackQA *taskTrackQAEMCJE01cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,1,iPhysicsSelectionFlagEMCEJE);
  AliPWG4HighPtTrackQA *taskTrackQAEMCJE70cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,0,iPhysicsSelectionFlagEMCEJE);
  AliPWG4HighPtTrackQA *taskTrackQAEMCJE71cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,1,iPhysicsSelectionFlagEMCEJE);
  AliPWG4HighPtTrackQA *taskTrackQAEMCJE72cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,2,iPhysicsSelectionFlagEMCEJE);
  AliPWG4HighPtTrackQA *taskTrackQAEMCJE05cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,5,iPhysicsSelectionFlagEMCEJE);
  AliPWG4HighPtTrackQA *taskTrackQAEMCJE74cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,4,iPhysicsSelectionFlagEMCEJE);
  AliPWG4HighPtTrackQA *taskTrackQAEMCJE75cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,5,iPhysicsSelectionFlagEMCEJE);
  AliPWG4HighPtTrackQA *taskTrackQAEMCJE40cent10 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,4,0,iPhysicsSelectionFlagEMCEJE);

  if(isPbPb) {
    for(cent=0; cent<4; cent++) {

      AliPWG4HighPtTrackQA *taskTrackQA00 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0,iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA01 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,1,iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA70 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,0,iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA71 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,1,iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA72 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,2,iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA05 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,5,iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA74 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,4,iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA75 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,5,iPhysicsSelectionFlag);
      AliPWG4HighPtTrackQA *taskTrackQA40 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,4,0,iPhysicsSelectionFlag);

      AliPWG4HighPtTrackQA *taskTrackQAEMCJE00 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0,iPhysicsSelectionFlagEMCEJE);
      AliPWG4HighPtTrackQA *taskTrackQAEMCJE01 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,1,iPhysicsSelectionFlagEMCEJE);
      AliPWG4HighPtTrackQA *taskTrackQAEMCJE70 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,0,iPhysicsSelectionFlagEMCEJE);
      AliPWG4HighPtTrackQA *taskTrackQAEMCJE71 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,1,iPhysicsSelectionFlagEMCEJE);
      AliPWG4HighPtTrackQA *taskTrackQAEMCJE72 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,2,iPhysicsSelectionFlagEMCEJE);
      AliPWG4HighPtTrackQA *taskTrackQAEMCJE05 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,5,iPhysicsSelectionFlagEMCEJE);
      AliPWG4HighPtTrackQA *taskTrackQAEMCJE74 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,4,iPhysicsSelectionFlagEMCEJE);
      AliPWG4HighPtTrackQA *taskTrackQAEMCJE75 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,5,iPhysicsSelectionFlagEMCEJE);
      AliPWG4HighPtTrackQA *taskTrackQAEMCJE40 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,4,0,iPhysicsSelectionFlagEMCEJE);

    }
  }

}




void AddTaskPWG4HighPtTrackQAAllReduced2011(char *prodType = "LHC10h",Bool_t isPbPb=kTRUE, Int_t iAODanalysis = 0) 
{    

  int cent = 10;

  UInt_t iPhysicsSelectionFlagCentral = AliVEvent::kCentral;
  UInt_t iPhysicsSelectionFlagSemiCentral = AliVEvent::kSemiCentral;

  AliPWG4HighPtTrackQA *taskTrackQA00C = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0,iPhysicsSelectionFlagCentral);
  AliPWG4HighPtTrackQA *taskTrackQA01C = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,1,iPhysicsSelectionFlagCentral);
  AliPWG4HighPtTrackQA *taskTrackQA21C = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,2,1,iPhysicsSelectionFlagCentral);
  AliPWG4HighPtTrackQA *taskTrackQA70C = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,0,iPhysicsSelectionFlagCentral);
  AliPWG4HighPtTrackQA *taskTrackQA71C = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,1,iPhysicsSelectionFlagCentral);
  AliPWG4HighPtTrackQA *taskTrackQA72C = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,2,iPhysicsSelectionFlagCentral);
  AliPWG4HighPtTrackQA *taskTrackQA05C = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,5,iPhysicsSelectionFlagCentral); 
  AliPWG4HighPtTrackQA *taskTrackQA74C = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,4,iPhysicsSelectionFlagCentral);
  AliPWG4HighPtTrackQA *taskTrackQA75C = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,5,iPhysicsSelectionFlagCentral);

  AliPWG4HighPtTrackQA *taskTrackQA00SC = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0,iPhysicsSelectionFlagSemiCentral);
  AliPWG4HighPtTrackQA *taskTrackQA01SC = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,1,iPhysicsSelectionFlagSemiCentral);
  AliPWG4HighPtTrackQA *taskTrackQA21SC = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,2,1,iPhysicsSelectionFlagSemiCentral);
  AliPWG4HighPtTrackQA *taskTrackQA70SC = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,0,iPhysicsSelectionFlagSemiCentral);
  AliPWG4HighPtTrackQA *taskTrackQA71SC = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,1,iPhysicsSelectionFlagSemiCentral);
  AliPWG4HighPtTrackQA *taskTrackQA72SC = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,2,iPhysicsSelectionFlagSemiCentral);
  AliPWG4HighPtTrackQA *taskTrackQA05SC = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,5,iPhysicsSelectionFlagSemiCentral);
  AliPWG4HighPtTrackQA *taskTrackQA74SC = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,4,iPhysicsSelectionFlagSemiCentral);
  AliPWG4HighPtTrackQA *taskTrackQA75SC = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,5,iPhysicsSelectionFlagSemiCentral);

}

void AddTaskPWG4HighPtTrackQAAOD(char *prodType = "LHC11h",Bool_t isPbPb=kTRUE, Int_t iAODanalysis = 1, Int_t filterBit = 768) 
{  
  UInt_t iPhysicsSelectionFlagMB = AliVEvent::kMB; 
  UInt_t iPhysicsSelectionFlagCentral = AliVEvent::kCentral;
  UInt_t iPhysicsSelectionFlagSemiCentral = AliVEvent::kSemiCentral;
  UInt_t iPhysicsSelectionFlagINT7 = AliVEvent::kINT7; 
  UInt_t iPhysicsSelectionFlagEMCEJE = AliVEvent::kEMCEJE; 

  Int_t cent = 10;

  Int_t filterBit1 = 256; //standard global tracks
  Int_t filterBit2 = 512; //complementary tracks

  Bool_t bIncludeNoITS = kFALSE;

  TString strRunPeriod = TString(prodType);
  strRunPeriod.ToLower();

  if (strRunPeriod == "lhc10h" ||strRunPeriod == "lhc11h" || strRunPeriod == "lhc13b" || strRunPeriod == "lhc13c" || strRunPeriod == "lhc13d" || strRunPeriod == "lhc13e" || strRunPeriod == "lhc13f" || strRunPeriod == "lhc13g" || strRunPeriod == "lhc12g" || strRunPeriod == "lhc12a15e" || strRunPeriod == "lhc13b4" || strRunPeriod == "lhc13b4_fix" || strRunPeriod == "lhc13b4_plus" || strRunPeriod == "lhc12a15f") {
    filterBit  = 768;
    filterBit1 = 256;
    filterBit2 = 512;
    bIncludeNoITS = kFALSE;
  }
  else if (strRunPeriod == "lhc11a" || strRunPeriod == "lhc10hold" || runPeriod.Contains("lhc12a15a") || runPeriod.Contains("lhc11a2")) {
    filterBit  = 272; 
    filterBit1 = 16;
    filterBit2 = 256;
    bIncludeNoITS = kTRUE;
  }

  if(isPbPb) {
    for(cent=0; cent<4; cent++) {
      AliPWG4HighPtTrackQA *taskTrackQAMB = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0,iPhysicsSelectionFlagMB);
      taskTrackQAMB->SetFilterMask(filterBit);
      taskTrackQAMB->SetIncludeNoITS(bIncludeNoITS);

      AliPWG4HighPtTrackQA *taskTrackQAMB1 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,5,iPhysicsSelectionFlagMB);
      taskTrackQAMB1->SetFilterMask(filterBit1);
      taskTrackQAMB1->SetIncludeNoITS(bIncludeNoITS);
      
      AliPWG4HighPtTrackQA *taskTrackQAMB2 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,5,iPhysicsSelectionFlagMB);
      taskTrackQAMB2->SetFilterMask(filterBit2);
      taskTrackQAMB2->SetIncludeNoITS(bIncludeNoITS);
    }

    cent = 10;

    if(strRunPeriod.Contains("lhc11h")) {
      AliPWG4HighPtTrackQA *taskTrackQAC = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0,iPhysicsSelectionFlagCentral);
      taskTrackQAC->SetFilterMask(filterBit);
      taskTrackQAC->SetIncludeNoITS(bIncludeNoITS);

      AliPWG4HighPtTrackQA *taskTrackQAC1 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,5,iPhysicsSelectionFlagCentral);
      taskTrackQAC1->SetFilterMask(filterBit1);
      taskTrackQAC1->SetIncludeNoITS(bIncludeNoITS);

      AliPWG4HighPtTrackQA *taskTrackQAC2 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,5,iPhysicsSelectionFlagCentral);
      taskTrackQAC2->SetFilterMask(filterBit2);
      taskTrackQAC2->SetIncludeNoITS(bIncludeNoITS);

      AliPWG4HighPtTrackQA *taskTrackQASC = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0,iPhysicsSelectionFlagSemiCentral);
      taskTrackQASC->SetFilterMask(filterBit);
      taskTrackQASC->SetIncludeNoITS(bIncludeNoITS);

      AliPWG4HighPtTrackQA *taskTrackQASC1 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,5,iPhysicsSelectionFlagSemiCentral);
      taskTrackQASC1->SetFilterMask(filterBit1);
      taskTrackQASC1->SetIncludeNoITS(bIncludeNoITS);

      AliPWG4HighPtTrackQA *taskTrackQASC2 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,5,iPhysicsSelectionFlagSemiCentral);
      taskTrackQASC2->SetFilterMask(filterBit2);
      taskTrackQASC2->SetIncludeNoITS(bIncludeNoITS);
    }
  }
  else {
    cent = 10;

    if(strRunPeriod.Contains("lhc13")) {
      AliPWG4HighPtTrackQA *taskTrackQAMB = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0,iPhysicsSelectionFlagINT7);
      taskTrackQAMB->SetFilterMask(filterBit);
      taskTrackQAMB->SetIncludeNoITS(bIncludeNoITS);

      AliPWG4HighPtTrackQA *taskTrackQAMB1 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,5,iPhysicsSelectionFlagINT7);
      taskTrackQAMB1->SetFilterMask(filterBit1);
      taskTrackQAMB1->SetIncludeNoITS(bIncludeNoITS);

      AliPWG4HighPtTrackQA *taskTrackQAMB2 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,5,iPhysicsSelectionFlagINT7);
      taskTrackQAMB2->SetFilterMask(filterBit2);
      taskTrackQAMB2->SetIncludeNoITS(bIncludeNoITS);

      if(strRunPeriod.EqualTo("LHC13d") || strRunPeriod.EqualTo("LHC13e") || strRunPeriod.EqualTo("LHC13f") || strRunPeriod.EqualTo("LHC13g")) {
	AliPWG4HighPtTrackQA *taskTrackQAEMCEJE = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0,iPhysicsSelectionFlagEMCEJE);
	taskTrackQAEMCEJE->SetFilterMask(filterBit);
	taskTrackQAEMCEJE->SetIncludeNoITS(bIncludeNoITS);

	AliPWG4HighPtTrackQA *taskTrackQAEMCEJE1 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,5,iPhysicsSelectionFlagEMCEJE);
	taskTrackQAEMCEJE1->SetFilterMask(filterBit1);
	taskTrackQAEMCEJE1->SetIncludeNoITS(bIncludeNoITS);

	AliPWG4HighPtTrackQA *taskTrackQAEMCEJE2 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,5,iPhysicsSelectionFlagEMCEJE);
	taskTrackQAEMCEJE2->SetFilterMask(filterBit2);
	taskTrackQAEMCEJE2->SetIncludeNoITS(bIncludeNoITS);
      }
    }
    else {
      AliPWG4HighPtTrackQA *taskTrackQAMB = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,0,iPhysicsSelectionFlagMB);
      taskTrackQAMB->SetFilterMask(filterBit);
      taskTrackQAMB->SetIncludeNoITS(bIncludeNoITS);

      AliPWG4HighPtTrackQA *taskTrackQAMB1 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,0,5,iPhysicsSelectionFlagMB);
      taskTrackQAMB1->SetFilterMask(filterBit1);
      taskTrackQAMB1->SetIncludeNoITS(bIncludeNoITS);
      

      AliPWG4HighPtTrackQA *taskTrackQAMB2 = ConfigureTaskPWG4HighPtTrackQA(prodType,isPbPb,iAODanalysis,cent,7,5,iPhysicsSelectionFlagMB);
      taskTrackQAMB2->SetFilterMask(filterBit2);
      taskTrackQAMB2->SetIncludeNoITS(bIncludeNoITS);
    }
  }
}

AliPWG4HighPtTrackQA* ConfigureTaskPWG4HighPtTrackQA(char *prodType = "LHC10e14",Bool_t isPbPb=kTRUE,Int_t iAODanalysis = 0, Int_t centClass = 0, Int_t trackType = 0, Int_t cuts = 0, UInt_t iPhysicsSelectionFlag = AliVEvent::kMB)
{
  /*
    trackType: 0 = global
               1 = TPC stand alone
               2 = TPC stand alone constrained to SPD vertex
               4 = TPC stand alone constrained to SPD vertex with QA track selection on global tracks
	       5 = Hybrid tracks: constrained TPConly for which no tight ITS is available
               6 = Hybrid tracks: constrained loose global for which no tight ITS is available
    cuts:      0 (global) = standard ITSTPC2010 a la RAA analysis
               1 (global) = ITSrefit, no SPD requirements -> standard for jet analysis
               2 (global) = ITSrefit + no hits in SPD
	       3 (global) = standard ITS tight cuts with nCrossed rows cut for hybrid tracks
               0 (TPC)    = standard TPC + NClusters>70
               1 (TPC)    = standard TPC + NClusters>0 --> to study new TPC QA recommendations
               0 (hybrid 5) = constrained TPConly for which no tight ITS is available
               0 (hybrid 6) = constrained loose global for which no tight ITS is available
   */

  //Load common track cut class
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/macros/CreateTrackCutsPWGJE.C");

  // Creates HighPtTrackQA analysis task and adds it to the analysis manager.
  
  // A. Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskPWG4HighPtQMC", "No analysis manager to connect to.");
    return NULL;
  }  

  // B. Check the analysis type using the event handlers connected to the analysis
  //    manager. The availability of MC handler can also be checked here.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddPWG4TaskHighPtTrackQA", "This task requires an input event handler");
    return NULL;
  }  

  // C. Create the task, add it to manager.
  //===========================================================================
 
  //CREATE THE  CUTS -----------------------------------------------
  //Use AliESDtrackCuts
  AliESDtrackCuts *trackCuts = new AliESDtrackCuts("AliESDtrackCuts","Standard Cuts");
  AliESDtrackCuts *trackCutsReject = 0x0;
  AliESDtrackCuts *trackCutsTPConly = new AliESDtrackCuts("AliESDtrackCutsTPConly","TPC only Cuts");

  //Standard Cuts
  //Set track cuts for global tracks
  if(trackType==0 && cuts==0) {
    // tight global tracks - RAA analysis
    trackCuts = CreateTrackCutsPWGJE(1000);
  }
  if(trackType==0 && cuts==1) {
    //Cuts global tracks with ITSrefit requirement and SPDrequirement for jet analysis
    trackCuts = CreateTrackCutsPWGJE(10001006);
  }
  if(trackType==0 && cuts==5) {
    //Cuts global tracks with ITSrefit requirement and SPDrequirement for jet analysis + NCrossedRowsCut>70 recommended in 2011
    trackCuts = CreateTrackCutsPWGJE(10001008);
  }
  if(trackType==0 && cuts==6) {
    //Cuts global tracks with ITSrefit requirement and no SPDrequirement for jet analysis + NCrossedRowsCut>70 recommended in 2011
    trackCuts = CreateTrackCutsPWGJE(10051008);
  }
  
  if(trackType==0 && cuts==2) {
    //Cuts global tracks with ITSrefit requirement but without SPD
    trackCuts = CreateTrackCutsPWGJE(10011006);
  }
  if(trackType==7 && cuts==0) {
    // tight global tracks
    trackCuts = CreateTrackCutsPWGJE(10041006);
    trackCutsReject = CreateTrackCutsPWGJE(1006);
    trackCutsReject->SetEtaRange(-0.9,0.9);
    trackCutsReject->SetPtRange(0.15, 1e10);
  }
  if(trackType==7 && cuts==4) {
    // tight global tracks +  NCrossedRowsCut>120 recommended in 2011
    trackCuts = CreateTrackCutsPWGJE(10041008);
    trackCutsReject = CreateTrackCutsPWGJE(1008);
    trackCutsReject->SetEtaRange(-0.9,0.9);
    trackCutsReject->SetPtRange(0.15, 1e10);
  }
  if(trackType==7 && cuts==1) {
    // tight global tracks
    trackCuts = CreateTrackCutsPWGJE(10011006);
  }
  if(trackType==7 && cuts==5) {
    // tight global tracks  + NCrossedRowsCut>70 recommended in 2011
    trackCuts = CreateTrackCutsPWGJE(10011008);
  }
  if(trackType==7 && cuts==2) {
    // no requirements on SPD and ITSrefit failed
    trackCuts = CreateTrackCutsPWGJE(10041006);       //no ITSrefit requirement filter 256
    trackCutsReject = CreateTrackCutsPWGJE(10001006); //ITSrefit requirement filter 16
    trackCutsReject->SetEtaRange(-0.9,0.9);
    trackCutsReject->SetPtRange(0.15, 1e10);
  }
  if(trackType==7 && cuts==6) {
    // no requirements on SPD and ITSrefit failed
    trackCuts = CreateTrackCutsPWGJE(10041008);       //no ITSrefit requirement filter 256
    trackCutsReject = CreateTrackCutsPWGJE(10001008); //ITSrefit requirement filter 16
    trackCutsReject->SetEtaRange(-0.9,0.9);
    trackCutsReject->SetPtRange(0.15, 1e10);
  }

  if(trackType==1 && cuts==0) {
    //Set track cuts for TPConly tracks
    trackCuts = CreateTrackCutsPWGJE(2001);
  }
  if(trackType==1 && cuts==1) {
    //Set track cuts for TPConly tracks
    trackCuts = CreateTrackCutsPWGJE(10032001);
  }

  if(trackType==2 && cuts==0) {
    //Set track cuts for TPConly constrained tracks
    trackCuts = CreateTrackCutsPWGJE(2001);
  }
  if(trackType==2 && cuts==1) {
    //Set track cuts for TPConly constrained tracks w/o cut on NClusters or NCrossedRows
    trackCuts = CreateTrackCutsPWGJE(10032001);
  }

  if(trackType==4 && cuts==0) {
     //	      Set track cuts for TPConly constrained tracks
    trackCuts = CreateTrackCutsPWGJE(2001);
  }
  if(trackType==4 && cuts==1) {
     //	      Set track cuts for TPConly constrained tracks
    trackCuts = CreateTrackCutsPWGJE(10032001);
  }
  if(trackType==5 || trackType==6) {
    // tight global tracks
    trackCuts = CreateTrackCutsPWGJE(1003);

    trackCutsReject = CreateTrackCutsPWGJE(10021003); 
    
    trackCutsTPConly = CreateTrackCutsPWGJE(2002);

    trackCutsReject->SetEtaRange(-0.9,0.9);
    trackCutsReject->SetPtRange(0.15, 1e10);
    
    trackCutsTPConly->SetEtaRange(-0.9,0.9);
    trackCutsTPConly->SetPtRange(0.15, 1e10);

  }

  trackCuts->SetEtaRange(-0.9,0.9);
  trackCuts->SetPtRange(0.15, 1e10);
  
  TString trigName = "";
  if (iPhysicsSelectionFlag == AliVEvent::kAnyINT)
    trigName += "kAnyINT";
  else if (iPhysicsSelectionFlag == AliVEvent::kAny)
    trigName += "kAny";
  else if(iPhysicsSelectionFlag == AliVEvent::kINT7)
    trigName += "kINT7";
  else if(iPhysicsSelectionFlag == AliVEvent::kMB)
    trigName += "kMB";
  else if(iPhysicsSelectionFlag == AliVEvent::kCentral)
    trigName += "kCentral";
  else if(iPhysicsSelectionFlag == AliVEvent::kSemiCentral)
    trigName += "kSemiCentral";
  else if(iPhysicsSelectionFlag == AliVEvent::kEMC7)
    trigName += "kEMC7";
  else if(iPhysicsSelectionFlag == AliVEvent::kEMCEJE)
    trigName += "kEMCEJE";
  else if(iPhysicsSelectionFlag == AliVEvent::kEMCEGA)
    trigName += "kEMCEGA";


  //Create the task
  AliPWG4HighPtTrackQA *taskPWG4TrackQA = new AliPWG4HighPtTrackQA(Form("AliPWG4HighPtTrackQACent%dTrack%dCuts%d%s",centClass,trackType,cuts,trigName.Data()));
  taskPWG4TrackQA->SetTrackType(trackType);
  taskPWG4TrackQA->SetCuts(trackCuts);
  taskPWG4TrackQA->SetCutsITSLoose(trackCutsReject);
  taskPWG4TrackQA->SetCutsTPConly(trackCutsTPConly);
  
  taskPWG4TrackQA->SetPtMax(100.);
 
  if(iAODanalysis)
    taskPWG4TrackQA->SetDataType(AliPWG4HighPtTrackQA::kAOD);
  else
    taskPWG4TrackQA->SetDataType(AliPWG4HighPtTrackQA::kESD);

  if(isPbPb) {
    taskPWG4TrackQA->SetIsPbPb(kTRUE);
    taskPWG4TrackQA->SetCentralityClass(centClass);
  }
  //  taskPWG4TrackQA->SetSigmaConstrainedMax(5.);

  cout << "iPhysicsSelectionFlag: " << iPhysicsSelectionFlag << endl;
  taskPWG4TrackQA->SelectCollisionCandidates(iPhysicsSelectionFlag);




  // E. Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += Form(":PWG4_HighPtTrackQACent%dTrackType%dCuts%d%s",centClass,trackType,cuts,trigName.Data());

  AliAnalysisDataContainer *cout_histQAtrack = 0x0;
  TString contName = Form("qa_histsQAtrackCent%dType%dcuts%d%s",centClass,trackType,cuts,trigName.Data());
  cout_histQAtrack = mgr->CreateContainer(contName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer,outputfile);

  mgr->AddTask(taskPWG4TrackQA);
  mgr->ConnectInput(taskPWG4TrackQA,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskPWG4TrackQA,1,cout_histQAtrack);

  // Return task pointer at the end
  return taskPWG4TrackQA;
}
