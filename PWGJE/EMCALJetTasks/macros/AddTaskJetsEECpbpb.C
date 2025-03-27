AliAnalysisTaskJetsEECpbpb* AddTaskJetsEECpbpb(
                             const char * njetsBase,
                             const char * njetsUS,
                             const char * njetsTrue,
                             const char * njetsPartLevel,
                             const Double_t R,
                             const char * nrhoBase,
                             const char * ntracks,
                             const char * ntracksUS,
                             const char * ntracksPartLevel,
                             const char * nclusters,
                             const char * ntracksTrue,
                             const char *type,
                             const char *CentEst,
                             Int_t       pSel,
                             TString     trigClass      = "",
                             TString     kEmcalTriggers = "",
                             TString     tag            = "",
                             AliAnalysisTaskJetsEECpbpb::JetShapeType jetShapeType = AliAnalysisTaskJetsEECpbpb::kMCTrue,
                             AliAnalysisTaskJetsEECpbpb::JetShapeSub jetShapeSub = AliAnalysisTaskJetsEECpbpb::kNoSub,
                             AliAnalysisTaskJetsEECpbpb::JetSelectionType jetSelection = AliAnalysisTaskJetsEECpbpb::kInclusive,
                             Float_t minpTHTrigger =0.,  Float_t maxpTHTrigger =0., Float_t acut =0.6, AliAnalysisTaskJetsEECpbpb::DerivSubtrOrder derivSubtrOrder = AliAnalysisTaskJetsEECpbpb::kSecondOrder){
 
cout<<"In the task MANAGER"<<endl;
// Get the pointer to the existing analysis manager via the static access method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskJetsEECpbpb","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      Error("AddTaskJetsEECpbpb", "This task requires an input event handler");
      return NULL;
    }

cout<<"Setting Wagon Names"<<endl;

TString wagonName1 = Form("JetsEECpbpb_%s_TC%s%s",njetsBase,trigClass.Data(),tag.Data());
TString wagonName2 = Form("JetsEECpbpb_%s_TC%s%sTree",njetsBase,trigClass.Data(),tag.Data());


TString configBasePath;

if(jetShapeType==AliAnalysisTaskJetsEECpbpb::kDetEmbPartPythia){

cout<<"---embedding:Getting config file"<<endl;
Bool_t getFromAlien=kTRUE;
TString configFileName = "Config_anraiJetsEECpbpb.C";
Int_t settingType = 0;
Int_t year = 2018;
TString periodName = "18q"; 
Int_t passIndex = 3;
const char* suffix = "";
Int_t containerNameMode=0;
TString combinedName;
combinedName.Form("anrai_%s", suffix);
TString configFilePath(configBasePath+configFileName);

if (getFromAlien)
  {
    std::cout<< "Info::anrai: Settings for grid testing"<<endl;
    gSystem->Exec(Form("alien_cp /alice/cern.ch/user/a/anrai/%s file:%s",configFileName.Data(),configFileName.Data()));
    configBasePath=Form("%s/",gSystem->pwd());
  } else {
    std::cout << " Info::anrai: Settings for local testing " << std::endl;
    configBasePath = "";
    suffix = "test";
    containerNameMode=0;
  }
  AliAnalysisTaskJetsEECpbpb *task = new AliAnalysisTaskJetsEECpbpb(wagonName1.Data());
  
  #ifdef __CLING__
      std::stringstream triggermakeradd;
      triggermakeradd << ".x " << configFilePath.Data() << "(";
      triggermakeradd << (getFromAlien ? "kTRUE" : "kFALSE") << ", ";
      triggermakeradd << settingType << ", ";
      triggermakeradd << year << ", ";
      triggermakeradd << "\"" << periodName.Data() << "\"" << ", ";
      triggermakeradd << passIndex << ", ";
      triggermakeradd << "\"" << combinedName.Data() << "\"";
      triggermakeradd << ")";
      std::string triggermakeraddstring = triggermakeradd.str();
      std::cout << triggermakeraddstring << std::endl;
      task = (AliAnalysisTaskJetsEECpbpb*)gROOT->ProcessLine(triggermakeraddstring.c_str());
  #else
  gROOT->LoadMacro(configFilePath.Data());
  task = Config_anraiJetsEECpbpb(getFromAlien,settingType,year,periodName,passIndex,combinedName);
  #endif
 
  std::cout << " Info::anrai: Configpath:  " << configFilePath << " year = " << year << " --- period name = " << periodName << " --- pass = " << passIndex << " --- settingType = " << settingType << std::endl;
  task->SetCaloTriggerPatchInfoName(kEmcalTriggers.Data());
  task->SetCentralityEstimator(CentEst);
  task->SelectCollisionCandidates(pSel);
  task->SetUseAliAnaUtils(kFALSE);
  task->SetJetShapeType(jetShapeType);

  mgr->AddTask(task);
  
  //Connnect input
  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer() );

  //Connect output
  TString contName1(wagonName1);
  TString contName2(wagonName2);

  if (jetShapeType == AliAnalysisTaskJetsEECpbpb::kMCTrue) contName1 += "_MCTrue";
  if (jetShapeType == AliAnalysisTaskJetsEECpbpb::kData) contName1 += "_Data";
  if (jetShapeType == AliAnalysisTaskJetsEECpbpb::kGenOnTheFly) contName1 += "_GenOnTheFly";
  if (jetShapeType == AliAnalysisTaskJetsEECpbpb::kPythiaDef) contName1 +="_PythiaDef";
  if (jetShapeSub == AliAnalysisTaskJetsEECpbpb::kNoSub) contName1 += "_NoSub";
  if (jetShapeSub == AliAnalysisTaskJetsEECpbpb::kConstSub) contName1 += "_ConstSub";
  if (jetShapeSub == AliAnalysisTaskJetsEECpbpb::kEventSub) contName1 += "_EventSub";
  if (jetShapeSub == AliAnalysisTaskJetsEECpbpb::kDerivSub) contName1 += "_DerivSub";
  
  if (jetSelection == AliAnalysisTaskJetsEECpbpb::kInclusive) contName1 += "_Incl";
 
  if (jetShapeType == AliAnalysisTaskJetsEECpbpb::kMCTrue) contName2 += "_MCTrue";
  if (jetShapeType == AliAnalysisTaskJetsEECpbpb::kData) contName2 += "_Data";
  if (jetShapeType == AliAnalysisTaskJetsEECpbpb::kGenOnTheFly) contName2 += "_GenOnTheFly";
  if (jetShapeType == AliAnalysisTaskJetsEECpbpb::kPythiaDef) contName2 +="_PythiaDef";
  if (jetShapeSub == AliAnalysisTaskJetsEECpbpb::kNoSub) contName2 += "_NoSub";
  if (jetShapeSub == AliAnalysisTaskJetsEECpbpb::kConstSub) contName2 += "_ConstSub";
  if (jetShapeSub == AliAnalysisTaskJetsEECpbpb::kEventSub) contName2 += "_EventSub";
  if (jetShapeSub == AliAnalysisTaskJetsEECpbpb::kDerivSub) contName2 += "_DerivSub";
  
  if (jetSelection == AliAnalysisTaskJetsEECpbpb::kInclusive) contName2 += "_Incl";

  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,1,coutput1);

  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contName2.Data(), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,2,coutput2);

  return task;


}
else{
  //For running over data
  //Configure jet tagger task
  AliAnalysisTaskJetsEECpbpb *task = new AliAnalysisTaskJetsEECpbpb(wagonName1.Data());
  task->SetJetShapeType(jetShapeType);
  task->SetJetShapeSub(jetShapeSub);
  task->SetJetSelection(jetSelection);
  task->SetDerivativeSubtractionOrder(derivSubtrOrder);
  //task->SetNCentBins(4);
  task->SetGeneratorLevelName(ntracksTrue);//for truth mc info
  task->SetDetectorLevelName(ntracks);//for det mc info

  TString thename(njetsBase);
  //if(thename.Contains("Sub")) task->SetIsConstSub(kTRUE);
  //task->SetVzRange(-10.,10.);

  AliParticleContainer *trackCont;// = task->AddTrackContainer(ntracks);
 
  cout<<"Finding Track container"<<endl;
  if ((jetShapeSub==AliAnalysisTaskJetsEECpbpb::kConstSub || jetShapeSub==AliAnalysisTaskJetsEECpbpb::kEventSub ) && ((jetShapeType==AliAnalysisTaskJetsEECpbpb::kData) || (jetShapeType==AliAnalysisTaskJetsEECpbpb::kDetEmbPartPythia) || (jetShapeType==AliAnalysisTaskJetsEECpbpb::kPythiaDef))){
    trackCont = task->AddParticleContainer(ntracks);}
  else trackCont = task->AddTrackContainer(ntracks);

  //Printf("tracks() = %s, trackCont =%p", ntracks, trackCont);
  AliParticleContainer *trackContUS  = task->AddTrackContainer(ntracksUS);


  //Printf("tracksUS() = %s", ntracksUS);
  AliParticleContainer *trackContTrue = task->AddMCParticleContainer(ntracksTrue);


  //Printf("ntracksTrue() = %s, trackContTrue=%p ", ntracksTrue, trackContTrue);
   if (jetShapeType==AliAnalysisTaskJetsEECpbpb::kDetEmbPartPythia) trackContTrue->SetIsEmbedding(true);


cout<<" Adding track containers for part level "<<endl;

  AliParticleContainer *trackContPartLevel=0;
  if ((jetShapeSub==AliAnalysisTaskJetsEECpbpb::kConstSub) && ((jetShapeType==AliAnalysisTaskJetsEECpbpb::kMCTrue) || (jetShapeType==AliAnalysisTaskJetsEECpbpb::kPythiaDef))){
    trackContPartLevel = task->AddParticleContainer(ntracksPartLevel);
  }
  else {trackContPartLevel = task->AddMCParticleContainer(ntracksPartLevel);}

cout<<" Track container for partLevel added "<<endl;
  if (jetShapeType==AliAnalysisTaskJetsEECpbpb::kDetEmbPartPythia){ 
    trackContPartLevel->SetIsEmbedding(true);
    }

   //Printf("ntracksPartLevel() = %s, trackContPartLevel=%p ", ntracksPartLevel, trackContPartLevel);
  
  cout<<"Other Track containers are also set"<<endl;


  AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters);
  
  AliJetContainer *jetContBase=0x0;
  AliJetContainer *jetContReco=0x0;
  AliJetContainer *jetContUS=0x0;
  AliJetContainer *jetContTrue=0x0;
  AliJetContainer *jetContPart=0x0;
  TString strType(type);

  if ((jetShapeType==AliAnalysisTaskJetsEECpbpb::kMCTrue || (jetShapeType==AliAnalysisTaskJetsEECpbpb::kGenOnTheFly))) {
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackContPartLevel);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
    }
  }
  
  if (jetShapeType==AliAnalysisTaskJetsEECpbpb::kData){
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
      if(jetShapeSub==AliAnalysisTaskJetsEECpbpb::kConstSub) jetContBase->SetAreaEmcCut(-2);
    }
  }
  

  if (jetShapeType==AliAnalysisTaskJetsEECpbpb::kDetEmbPartPythia){

    cout<<"you chose detembPartPythia"<<endl;

    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetName("EmbJets");
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
     
      if(jetShapeSub==AliAnalysisTaskJetsEECpbpb::kConstSub) jetContBase->SetAreaEmcCut(-2);
    }

    jetContTrue = task->AddJetContainer(njetsTrue,strType,R);
    if(jetContTrue) {
      jetContTrue->SetName("GenJets");
      jetContTrue->SetRhoName(nrhoBase);
      jetContTrue->ConnectParticleContainer(trackContTrue);
      jetContTrue->SetPercAreaCut(acut);
    
    }

    if(jetShapeSub==AliAnalysisTaskJetsEECpbpb::kConstSub || jetShapeSub==AliAnalysisTaskJetsEECpbpb::kEventSub){
      jetContUS=task->AddJetContainer(njetsUS,strType,R);
      if(jetContUS) {
        jetContUS->SetRhoName(nrhoBase);
        jetContUS->ConnectParticleContainer(trackContUS);
        jetContUS->SetPercAreaCut(acut);
       
      }
    }
 
     jetContPart = task->AddJetContainer(njetsPartLevel,strType,R);
      if(jetContPart) {

        jetContPart->SetRhoName(nrhoBase);
        jetContPart->ConnectParticleContainer(trackContPartLevel);
        jetContPart->SetPercAreaCut(acut);
        
      }
  }
  
  if (jetShapeType==AliAnalysisTaskJetsEECpbpb::kPythiaDef){
    
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
    }
    
    jetContTrue = task->AddJetContainer(njetsTrue,strType,R);
    if(jetContTrue) {
      jetContTrue->SetRhoName(nrhoBase);
      jetContTrue->ConnectParticleContainer(trackContTrue);
      jetContTrue->SetPercAreaCut(acut);
      
    }
    
    if(jetShapeSub==AliAnalysisTaskJetsEECpbpb::kConstSub){
      jetContUS=task->AddJetContainer(njetsUS,strType,R);
      if(jetContUS) {
        jetContUS->SetRhoName(nrhoBase);
        jetContUS->ConnectParticleContainer(trackContUS);
        jetContUS->SetPercAreaCut(acut);
        
      }
    }
    
    jetContPart = task->AddJetContainer(njetsPartLevel,strType,R);
    if(jetContPart) {
      jetContPart->SetRhoName(nrhoBase);
      jetContPart->ConnectParticleContainer(trackContPartLevel);
      jetContPart->SetPercAreaCut(acut);
    }
    
  }

  task->SetCaloTriggerPatchInfoName(kEmcalTriggers.Data());
  task->SetCentralityEstimator(CentEst);
  task->SelectCollisionCandidates(pSel);
  task->SetUseAliAnaUtils(kFALSE);

  mgr->AddTask(task);
  
  //Connnect input
  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer() );

  //Connect output
  TString contName1(wagonName1);
  TString contName2(wagonName2);

  if (jetShapeType == AliAnalysisTaskJetsEECpbpb::kMCTrue) contName1 += "_MCTrue";
  if (jetShapeType == AliAnalysisTaskJetsEECpbpb::kData) contName1 += "_Data";
  if (jetShapeType == AliAnalysisTaskJetsEECpbpb::kGenOnTheFly) contName1 += "_GenOnTheFly";
  if (jetShapeType == AliAnalysisTaskJetsEECpbpb::kPythiaDef) contName1 +="_PythiaDef";
  if (jetShapeSub == AliAnalysisTaskJetsEECpbpb::kNoSub) contName1 += "_NoSub";
  if (jetShapeSub == AliAnalysisTaskJetsEECpbpb::kConstSub) contName1 += "_ConstSub";
  if (jetShapeSub == AliAnalysisTaskJetsEECpbpb::kEventSub) contName1 += "_EventSub";
  if (jetShapeSub == AliAnalysisTaskJetsEECpbpb::kDerivSub) contName1 += "_DerivSub";
  
  if (jetSelection == AliAnalysisTaskJetsEECpbpb::kInclusive) contName1 += "_Incl";
 
  if (jetShapeType == AliAnalysisTaskJetsEECpbpb::kMCTrue) contName2 += "_MCTrue";
  if (jetShapeType == AliAnalysisTaskJetsEECpbpb::kData) contName2 += "_Data";
  if (jetShapeType == AliAnalysisTaskJetsEECpbpb::kGenOnTheFly) contName2 += "_GenOnTheFly";
  if (jetShapeType == AliAnalysisTaskJetsEECpbpb::kPythiaDef) contName2 +="_PythiaDef";
  if (jetShapeSub == AliAnalysisTaskJetsEECpbpb::kNoSub) contName2 += "_NoSub";
  if (jetShapeSub == AliAnalysisTaskJetsEECpbpb::kConstSub) contName2 += "_ConstSub";
  if (jetShapeSub == AliAnalysisTaskJetsEECpbpb::kEventSub) contName2 += "_EventSub";
  if (jetShapeSub == AliAnalysisTaskJetsEECpbpb::kDerivSub) contName2 += "_DerivSub";
  
  if (jetSelection == AliAnalysisTaskJetsEECpbpb::kInclusive) contName2 += "_Incl";

  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
   mgr->ConnectOutput(task,1,coutput1);


   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contName2.Data(), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,2,coutput2);

  return task;
}
 

}

