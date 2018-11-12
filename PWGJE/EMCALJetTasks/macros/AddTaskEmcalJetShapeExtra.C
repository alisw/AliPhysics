AliAnalysisTaskEmcalJetShapeExtra* AddTaskEmcalJetShapeExtra(const char * njetsBase,
                                             const char * njetsUS,
            
						     const char * njetsTrue,
                                             const char * njetsPartLevel,
        
						     const Double_t R,
                                             const char * nrhoBase,

						     const char * ntracks, 
                                             const char * ntracksUS,
                                             const char *ntracksPartLevel,
             
						     const char * nclusters,
						     const char * ntracksTrue,
						     const char *type,				      
						     const char *CentEst,
						     Int_t       pSel,
						     TString     trigClass      = "",
						     TString     kEmcalTriggers = "",
						     TString     tag            = "",
							     AliAnalysisTaskEmcalJetShapeExtra::JetShapeType jetShapeType=AliAnalysisTaskEmcalJetShapeExtra::kPythiaDef,
							     AliAnalysisTaskEmcalJetShapeExtra::JetShapeSub jetShapeSub=AliAnalysisTaskEmcalJetShapeExtra::kNoSub
                                             
)
{
 

  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskEmcalJetShapeExtra","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kTRUE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

 
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskEmcalJetShapeExtra", "This task requires an input event handler");
      return NULL;
    }

  TString wagonName1 = Form("JetQGTaggings_%s_TC%s%s",njetsBase,trigClass.Data(),tag.Data());
  TString wagonName2 = Form("JetQGTaggings_%s_TC%s%sTree",njetsBase,trigClass.Data(),tag.Data());
  //Configure jet tagger task
  AliAnalysisTaskEmcalJetShapeExtra *task = new AliAnalysisTaskEmcalJetShapeExtra(wagonName1.Data());

  //task->SetNCentBins(4);
   TString thename(njetsBase);
  //if(thename.Contains("Sub")) task->SetIsConstSub(kTRUE);
 task->SetVzRange(-10.,10.);

  AliParticleContainer *trackCont;// = task->AddTrackContainer(ntracks);
 
  if ((jetShapeType==AliAnalysisTaskEmcalJetShapeExtra::kData)  || (jetShapeType==AliAnalysisTaskEmcalJetShapeExtra::kPythiaDef)){
    trackCont = task->AddParticleContainer(ntracks);}
  else trackCont = task->AddTrackContainer(ntracks);

  
  Printf("tracks() = %s, trackCont =%p", ntracks, trackCont);
  AliParticleContainer *trackContUS  = task->AddTrackContainer(ntracksUS);
  Printf("tracksUS() = %s", ntracksUS);
  AliParticleContainer *trackContTrue = task->AddMCParticleContainer(ntracksTrue);
  Printf("ntracksTrue() = %s, trackContTrue=%p ", ntracksTrue, trackContTrue);
  
  AliParticleContainer *trackContPartLevel=0;
  
    
    if ((jetShapeSub==AliAnalysisTaskEmcalJetShapeExtra::kConstSub) &&  (jetShapeType==AliAnalysisTaskEmcalJetShapeExtra::kPythiaDef)){
    trackContPartLevel = task->AddParticleContainer(ntracksPartLevel);
  }
  else trackContPartLevel = task->AddMCParticleContainer(ntracksPartLevel);
  
   Printf("ntracksPartLevel() = %s, trackContPartLevel=%p ", ntracksPartLevel, trackContPartLevel);
  

  AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters);
  
  AliJetContainer *jetContBase=0x0;
  AliJetContainer *jetContUS=0x0;
  AliJetContainer *jetContTrue=0x0;
 AliJetContainer *jetContPart=0x0;
  TString strType(type);


  
  if (jetShapeType==AliAnalysisTaskEmcalJetShapeExtra::kData){
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
    //  jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(0.6);
      if(jetShapeSub==AliAnalysisTaskEmcalJetShapeExtra::kConstSub) jetContBase->SetAreaEmcCut(-2);
    }    
  }
  

  if (jetShapeType==AliAnalysisTaskEmcalJetShapeExtra::kPythiaDef){
    
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(0.6);
    }
    
    jetContTrue = task->AddJetContainer(njetsTrue,strType,R);
    if(jetContTrue) {
      jetContTrue->SetRhoName(nrhoBase);
      jetContTrue->ConnectParticleContainer(trackContTrue);
      jetContTrue->SetPercAreaCut(0.6);
      
      }
    
    if(jetShapeSub==AliAnalysisTaskEmcalJetShapeExtra::kConstSub){
      jetContUS=task->AddJetContainer(njetsUS,strType,R);
      if(jetContUS) {
        jetContUS->SetRhoName(nrhoBase);
        jetContUS->ConnectParticleContainer(trackContUS);
        jetContUS->SetPercAreaCut(0.6);
        
      }
        }
    
  jetContPart = task->AddJetContainer(njetsPartLevel,strType,R);
    if(jetContPart) {
      jetContPart->SetRhoName(nrhoBase);
      jetContPart->ConnectParticleContainer(trackContPartLevel);
      jetContPart->SetPercAreaCut(0.6);
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

 
  if (jetShapeType == AliAnalysisTaskEmcalJetShapeExtra::kData) contName1 += "_Data";
 
  if (jetShapeType == AliAnalysisTaskEmcalJetShapeExtra::kPythiaDef) contName1 +="_PythiaDef";
  if (jetShapeSub == AliAnalysisTaskEmcalJetShapeExtra::kNoSub) contName1 += "_NoSub";
  if (jetShapeSub == AliAnalysisTaskEmcalJetShapeExtra::kConstSub) contName1 += "_ConstSub";
 
  if (jetShapeType == AliAnalysisTaskEmcalJetShapeExtra::kData) contName2 += "_Data";
 
  if (jetShapeType == AliAnalysisTaskEmcalJetShapeExtra::kPythiaDef) contName2 +="_PythiaDef";
  if (jetShapeSub == AliAnalysisTaskEmcalJetShapeExtra::kNoSub) contName2 += "_NoSub";
  if (jetShapeSub == AliAnalysisTaskEmcalJetShapeExtra::kConstSub) contName2 += "_ConstSub";


  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
   mgr->ConnectOutput(task,1,coutput1);


   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contName2.Data(), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,2,coutput2);
   
  return task;  

}

