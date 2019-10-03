AliAnalysisTaskJetFFMoments* AddTaskJetFFMoments(char* bGen = "AODMC2b", char* bRec1 = "AOD", UInt_t filterMask = 128, UInt_t iPhysicsSelectionFlag = AliVEvent::kAny, Char_t *jf = "ANTIKT", Float_t radius = 0.4, Int_t nSkip = 0, Bool_t kReadJetBranch = kTRUE, char* suffix ="", Float_t ptTrackCut = 0.15, Float_t etaTrackWindow = 0.9, Float_t ptJetCut = 5.,  char* anaJetType = "leading", Float_t vertexWindow = 10., Double_t ffmPower = 2, Int_t bType = -1, Double_t bcut1 = 0.4, Double_t bcut2 = TMath::Pi(), double mu = 25, Int_t nUsedJets = 8, char* jfTask = "clusters", Bool_t kRandom = kFALSE);

AliAnalysisTaskJetFFMoments* AddTaskJetFFMoments(Float_t radius = 0.4, Float_t ptTrackCut = 0.15, Char_t *jf = "ANTIKT", Float_t ptJetCut = 5, char* anaJetType = "leading", Double_t ffmPower = 2, char* bGen = "", Int_t bType = -1, Double_t bcut1 = 0.4, Double_t bcut2 = TMath::Pi(), double mu = 25, Int_t nUsedJets = 8, char* jfTask = "clusters", Int_t nSkip = 0, Bool_t kRandom = kFALSE, char* suffix);

AliAnalysisTaskJetFFMoments *AddTaskJetFFMoments(Float_t radius = 0.4, char* bGen = "KINE2B", Double_t ffmPower = 2, char* bRec1 = "KINEDET", char* suffix ="", Int_t bType = -1, Bool_t kReadJetBranch = kFALSE, char* anaJetType="leading", Float_t ptJetCut = 5., Double_t bcut1 = 0.4, Double_t bcut2 = TMath::Pi(), double mu = 25, Int_t nUsedJets = 8, char* jfTask = "clusters", Char_t *jf = "ANTIKT", Int_t nSkip = 0, Bool_t kRandom = kFALSE);

AliAnalysisTaskJetFFMoments *AddTaskJetFFMoments(Float_t radius = 0.4, char* bGen = "KINE2B", Double_t ffmPower = 2, char* bRec1 = "KINEDET", Int_t bType = -1, Bool_t kReadJetBranch = kFALSE, char* anaJetType="leading", Float_t ptJetCut = 5., Double_t bcut1 = 0.4, Double_t bcut2 = TMath::Pi(), double mu = 25, Int_t nUsedJets = 8, char* jfTask = "clusters", Char_t *jf = "ANTIKT", Int_t nSkip = 0, Bool_t kRandom = kFALSE, char* suffix);

AliAnalysisTaskJetFFMoments *AddTaskJetFFMoments(Float_t radius, Float_t ptTrackCut, Char_t *jf, Float_t ptJetCut, char* anaJetType,
                                                 Double_t ffmPower, char* bGen, Int_t bType, Double_t bcut1, Double_t bcut2, double mu,
                                                 Int_t nUsedJets, char* jfTask, Int_t nSkip, Bool_t kRandom, char* suffix)
{

AddTaskJetFFMoments(bGen, "AOD", AliAnalysisManager::GetGlobalInt("kHighPtFilterMask",gDebug), AliAnalysisManager::GetGlobalInt("kPhysicsSelectionFlag",gDebug), jf, radius, nSkip, kTRUE, suffix , ptTrackCut,  AliAnalysisManager::GetGlobalDbl("kTrackEtaWindow",gDebug), ptJetCut, anaJetType, AliAnalysisManager::GetGlobalDbl("kVertexWindow",gDebug), ffmPower, bType, bcut1, bcut2, mu, nUsedJets, jfTask, kRandom);

}


AliAnalysisTaskJetFFMoments *AddTaskJetFFMoments(Float_t radius, char* bGen, Double_t ffmPower,
                                                 char* bRec1, char* suffix ,Int_t bType, Bool_t kReadJetBranch, char* anaJetType, Float_t ptJetCut,
                                                 Double_t bcut1, Double_t bcut2, double mu,
                                                 Int_t nUsedJets, char* jfTask,  Char_t *jf, Int_t nSkip, Bool_t kRandom)
{

AddTaskJetFFMoments(bGen, bRec1, 0 , 0 , jf, radius, nSkip,kReadJetBranch, suffix , 0.15, 0.9 , ptJetCut, anaJetType, 10 , ffmPower, bType, bcut1, bcut2, mu, nUsedJets, jfTask, kRandom);

}

AliAnalysisTaskJetFFMoments *AddTaskJetFFMoments(Float_t radius, char* bGen, Double_t ffmPower,
                                                 char* bRec1, Int_t bType, Bool_t kReadJetBranch, char* anaJetType, Float_t ptJetCut,
                                                 Double_t bcut1, Double_t bcut2, double mu,
                                                 Int_t nUsedJets, char* jfTask,  Char_t *jf, Int_t nSkip, Bool_t kRandom, char* suffix)
{

AddTaskJetFFMoments(bGen, bRec1, 0 , 0 , jf, radius, nSkip,kReadJetBranch, suffix , 0.15, 0.9 , ptJetCut, anaJetType, 10 , ffmPower, bType, bcut1, bcut2, mu, nUsedJets, jfTask, kRandom);

}

AliAnalysisTaskJetFFMoments *AddTaskJetFFMoments(char* bGen, char* bRec1, UInt_t filterMask, UInt_t iPhysicsSelectionFlag, 
						 Char_t *jf, Float_t radius, Int_t nSkip, Bool_t kReadJetBranch, char* suffix, 
                                                 Float_t ptTrackCut, Float_t etaTrackWindow, Float_t ptJetCut, char* anaJetType,  
                                                 Float_t vertexWindow, Double_t ffmPower, Int_t bType, Double_t bcut1, Double_t bcut2, double mu,
                                                 Int_t nUsedJets, char* jfTask, Bool_t kRandom)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
 AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
     ::Error("AddTaskJetFFMoments", "No analysis manager to connect to.");
     return NULL;
  }  

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskJetFFMoments", "This task requires an input event handler");
     return NULL;
  }

//Jet we need to analysis
 TString typeGen(bGen);
 TString typeRec1(bRec1);
 TString jetFindingTask(jfTask);

  if(typeGen.Length() != 0) {
    if(!typeGen.Contains("AODextra") && !typeGen.Contains("AODMCextra") && !typeGen.Contains("AODMC2b")) {
      typeGen.ToUpper();
    }
  }

  if(typeRec1.Length() != 0) {
    if(!typeRec1.Contains("AODextra") && !typeRec1.Contains("AODMCextra") && !typeRec1.Contains("AODMC2b")) {
      typeRec1.ToUpper();
    }
  }

  // Create the task and configure it.
  //===========================================================================
  // Define jet branch name 
  Int_t kBackgroundMode = 0;
  TString cAdd = "UndefinedJetType";
  cAdd += Form("_%s", jf);
  cAdd += Form("%02d", (int)((radius+0.01)*10.));
  cAdd += Form("_B%d", (int)kBackgroundMode);
  cAdd += Form("_Filter%05d", filterMask);
  cAdd += Form("_Cut%05d", (int)(1000.*ptTrackCut));
  if(nSkip>=0) cAdd += Form("_Skip%02d", nSkip);
  // random cone bkg is not supported yet!!!!!!!!!!
  if(kRandom) {
    cAdd += Form("_random");
    return NULL;
  }

  // Define analysis name
  TString kAnaName = Form("JetFFMoments_%s", cAdd.Data());
  if(typeGen.Length() != 0) {//MC
    if( typeRec1.Length() != 0)		kAnaName.ReplaceAll("UndefinedJetType","MC");
    else                                kAnaName.ReplaceAll("UndefinedJetType","Kine");
  } else {//Data
    if( typeRec1.Length() != 0)		kAnaName.ReplaceAll("UndefinedJetType","Data");
    else {
     printf("No jets in analysis!!\n");
     return NULL;
    }
  }

  TString sSuffix(suffix);
  if(bType!= -1) sSuffix += Form("_bkg%d", bType);
  int ffmpowerInt; if (std::floor(ffmPower) == ffmPower) ffmpowerInt = (int) ffmPower; else ffmpowerInt = (int) (ffmPower*10);
  if(ffmPower!= 2) sSuffix += Form("_Power%d",ffmpowerInt);
  if(sSuffix.Length() != 0) kAnaName += Form("_%s", sSuffix.Data());

  // Define task parameters
  AliAnalysisTaskJetFFMoments* ffm = new  AliAnalysisTaskJetFFMoments(kAnaName.Data());
  ffm->SetFilterMask(filterMask); 
  ffm->SetVtxCuts(vertexWindow, 1);
  ffm->SetTrackPtMin(ptTrackCut);
  ffm->SetTrackEtaWindow(etaTrackWindow);
  ffm->SetJetPtMin(ptJetCut);
  ffm->SetJetEtaWindow(ffm->GetTrackEtaMin()+radius, ffm->GetTrackEtaMax()-radius);
  ffm->SetDoJetReco(!kReadJetBranch); 
  ffm->SetNUsedJets(nUsedJets);
  ffm->SetRparam(radius);
  ffm->SetGhostArea(0.005);
  ffm->SetGhostEtaMax(etaTrackWindow);
  ffm->SetFFMAxis(750);
  ffm->SetTracksInJetMethod(0);
  ffm->SetAxisForTracks(1000,0,100);
  ffm->SetJetStructureAxis();
  if(typeGen.Contains("AODMC2b"))      ffm->SetTrackTypeGen(AliAnalysisTaskJetFFMoments::kTrackAODMCChargedAcceptance);
  else if(typeGen.Contains("AODMC2"))  ffm->SetTrackTypeGen(AliAnalysisTaskJetFFMoments::kTrackAODMCCharged);
  else if(typeGen.Contains("AODMC"))   ffm->SetTrackTypeGen(AliAnalysisTaskJetFFMoments::kTrackAODMCAll);
  else if(typeGen.Contains("KINE2B"))  ffm->SetTrackTypeGen(AliAnalysisTaskJetFFMoments::kTrackKineChargedAcceptance);
  else if(typeGen.Contains("KINE2"))   ffm->SetTrackTypeGen(AliAnalysisTaskJetFFMoments::kTrackKineCharged);
  else if(typeGen.Contains("KINE"))    ffm->SetTrackTypeGen(AliAnalysisTaskJetFFMoments::kTrackKineAll);
  else {if(typeGen.Length()) Printf("trackType Gen %s not found", typeGen.Data());}
  if(typeRec1.Contains("AOD"))         ffm->SetTrackTypeRec(AliAnalysisTaskJetFFMoments::kTrackAOD);
  else if(typeRec1.Contains("KINEDET")) {ffm->SetTrackTypeRec(AliAnalysisTaskJetFFMoments::kTrackKineChargedAcceptanceDet);}
  else {if(typeRec1.Length()) Printf("trackType Rec %s not found", typeRec1.Data());} 
 ffm->SetFFMScalePower(ffmPower);
 ffm->SetAnaJetType(anaJetType);
 if(bType==0)       {ffm->SetFFBckgMode(1); // FF: Perp; FFM: Calculate ffm bckg from a doughnut
                    ffm->SetFFMBckgTypeAndBounds(bType,bcut1*radius,bcut2*radius,mu);}
 else if(bType==1)  {ffm->SetFFBckgMode(1);  // FF: Perp; FFM: Calculate bckg from a rectangle around jet axis (radius,pi)
                    ffm->SetFFMBckgTypeAndBounds(bType,bcut1,bcut2,mu); }
 else if(bType==2)  {ffm->SetFFBckgMode(1);  // FF: Perp; FFM:  Calculate bckg in an eta window
                    ffm->SetFFMBckgTypeAndBounds(bType,bcut1,bcut2,mu); }
 else if(bType==3) {ffm->SetFFBckgMode(1); // FF: perp ; FFM: perpFFM
                    ffm->SetFFMBckgTypeAndBounds(bType,bcut1,bcut2,mu); }
 else if(bType==4) {ffm->SetFFBckgMode(2); // FF: perp2; FFM: perp2FFM
                    ffm->SetFFMBckgTypeAndBounds(bType,bcut1,bcut2,mu); }
 else if(bType==5) {ffm->SetFFBckgMode(1); // FF: perp; FFM: RapPhiRange
                    ffm->SetFFMBckgTypeAndBounds(bType,bcut1,bcut2,mu); }
 else {ffm->SetFFBckgMode(0); ffm->SetFFMBckgTypeAndBounds(bType);} // switch off bckg calc for FF and FFM

  // Define reading branch
  if(kReadJetBranch) { 
   //define reading branch name
   cAdd = Form("%s%s", jetFindingTask.Data(), cAdd.Data());
   TString kBranchName =cAdd;
    if(typeGen.Length() != 0) {//MC
     kBranchName.ReplaceAll("UndefinedJetType", typeGen.Data());
     ffm->SetJetBranches(0, kBranchName.Data());
     printf("Gen BranchName: %s\n", kBranchName.Data());
    }
    kBranchName =cAdd;
    if(typeRec1.Length() != 0) {
     kBranchName.ReplaceAll("UndefinedJetType", typeRec1.Data());
     ffm->SetJetBranches(1, kBranchName.Data());
     printf("Rec1 BranchName 1: %s\n", kBranchName.Data());
    }
   cAdd.ReplaceAll("clusters","");
  } else {// rec jets
   //define rec algorithm for jet rec
    switch (jf) {
     case "ANTIKT":
      ffm->SetAlgorithm(2); // antikt from fastjet/JetDefinition.hh
      break;
     case "CA":
      ffm->SetAlgorithm(1); // CA from fastjet/JetDefinition.hh
      break;
     case "KT":
      ffm->SetAlgorithm(0); // kt from fastjet/JetDefinition.hh
      break;
     default:
      ::Error("AddTaskJetFFMoments", "Wrong jet finder selected\n");
      return 0;
    }
  }


  //define analysis name
  if(typeGen.Length() != 0) {//MC
    if( typeRec1.Length() != 0)         cAdd.ReplaceAll("UndefinedJetType","MC");
    if( typeRec1.Contains("DET"))       cAdd.ReplaceAll("MC","KINE");
  } else {//Data
    ffm->SetHistosLevel(1);// no MC, can not do the correction from MC
    if( typeRec1.Length() != 0)         cAdd.ReplaceAll("UndefinedJetType","Data");
    else {
     printf("No jets in analysis!!\n");
     return NULL;
    }
  }


 TString type = mgr->GetInputEventHandler()->GetDataType();
 if(type == "AOD") {
   // Assume all jet are produced already
   ffm->SetAODTrackInput(kTRUE);
   ffm->SetAODMCInput(kTRUE);
 }

 TString Ana_Name = "FFMoments";
 if( cAdd.Data() != "" ) Ana_Name += Form("_%s", cAdd.Data());
 if( sSuffix.Length() != 0) Ana_Name += Form("_%s", sSuffix.Data());

 // Write new jet branch in 
 if(!kReadJetBranch){ 
   ffm->SetJetOutputBranch(Ana_Name.Data());
   ffm->SetFilterPt(5.); // store only jets / clusters above a certain threshold
 }

  if(typeGen.Length() == 0 && iPhysicsSelectionFlag) {
  // printf("iPhysicsSelectionFlag %d\n", iPhysicsSelectionFlag);
   ffm->SelectCollisionCandidates(iPhysicsSelectionFlag);
  }

 mgr->AddTask(ffm);

 // Create ONLY the output containers for the data produced by the task.
 // Get and connect other common input/output containers via the manager as below
 //==============================================================================
 AliAnalysisDataContainer *coutput1_ffm = mgr->CreateContainer(Ana_Name.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:PWGJE_%s", AliAnalysisManager::GetCommonFileName(), Ana_Name.Data()));

 mgr->ConnectInput  (ffm, 0, mgr->GetCommonInputContainer());
 if(mgr->GetCommonOutputContainer()) mgr->ConnectOutput (ffm, 0, mgr->GetCommonOutputContainer());
 mgr->ConnectOutput (ffm, 1, coutput1_ffm );
 
 return ffm;
}
