AliAnalysisTaskJetHBOM *AddTaskJetHBOM(char* bRec = "AOD",char* bGen = "",UInt_t filterMask = 272, UInt_t iPhysicsSelectionFlag = AliVEvent::kAny,Char_t *jf = "KT", Float_t radius = 0.4,Int_t kWriteAOD = kFALSE,char* deltaFile = "",Float_t ptTrackCut = 0.15, Float_t etaTrackWindow = 0.9,Float_t vertexWindow = 10.,TString effLoc = "$ALICE_PHYSICS/OADB/PWGJE/HBOM/fastMCInput_LHC10h_110719a.root",Int_t fNHBOM = 0, Int_t constCone = kFALSE, Float_t constConePhi = 0, Float_t constConeEta = 0);

Int_t kBackgroundModeCl = 0;
Float_t kPtTrackCutCl = 0.15;
Float_t kTrackEtaWindowCl = 0.8;
Float_t kVertexWindowCl = 10;


AliAnalysisTaskJetHBOM *AddTaskJetHBOM(char* bRec,char* bGen ,UInt_t filterMask,UInt_t iPhysicsSelectionFlag,Char_t *jf,Float_t radius,Int_t kWriteAOD,char *deltaFile,Float_t ptTrackCut,Float_t etaTrackWindow,Float_t vertexWindow,TString effLoc,Int_t fNHBOM, Int_t constCone, Float_t constConePhi,Float_t constConeEta)
 {
   //if constCone is true, the random Cone positon is set to constConePhi and constConeEta. Else the cone is random set and the parameters constConePhi and constConeEta are irrelevant

 // Creates a jet fider task, configures it and adds it to the analysis manager.
   kPtTrackCutCl = ptTrackCut;
   kTrackEtaWindowCl = etaTrackWindow;
   kVertexWindowCl = vertexWindow;

   TString outputFile(deltaFile);
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
       ::Error("AddTaskJetHBOM", "No analysis manager to connect to.");
       return NULL;
    }  

    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskJetHBOM", "This task requires an input event handler");
       return NULL;
    }

    TString type = mgr->GetInputEventHandler()->GetDataType();
    TString typeRec(bRec);
    TString typeGen(bGen);
    if(!typeRec.Contains("AODextra")) {
      typeGen.ToUpper();
      typeRec.ToUpper();
    }
    cout << "typeRec: " << typeRec << endl;
    // Create the task and configure it.
    //===========================================================================


    TString cAdd = "";
    cAdd += Form("%02d_",(int)((radius+0.01)*10.));
    cAdd += Form("B%d",(int)kBackgroundModeCl);
    cAdd += Form("_Filter%05d",filterMask);
    cAdd += Form("_Cut%05d",(int)(1000.*kPtTrackCutCl));
    cAdd += Form("_hbom%02d",fNHBOM);
    if(constCone){
      cAdd += Form("_constConePhi%02dEta%02d",constConePhi,constConeEta);
    }
    Printf("%s %E",cAdd.Data(),kPtTrackCutCl);
    AliAnalysisTaskJetHBOM* hbom = new  AliAnalysisTaskJetHBOM(Form("JetHBOM%s_%s%s",bRec,jf,cAdd.Data()));
      
   hbom->SetFilterMask(filterMask); 
   //   hbom->SetUseGlobalSelection(kTRUE); 
   hbom->SetVtxCuts(kVertexWindowCl,1);//sets fVtxZCut and fVtxR2Cut
   if(type == "AOD"){
     // Assume all jet are produced already
     hbom->SetAODTrackInput(kTRUE);
     hbom->SetAODMCInput(kTRUE);
   }

   if(typeRec.Contains("AODMC2b")){// work down from the top AODMC2b -> AODMC2 -> AODMC -> AOD
     hbom->SetTrackTypeRec(AliAnalysisTaskJetHBOM::kTrackAODMCChargedAcceptance);
     hbom->SetTrackPtCut(kPtTrackCutCl);
     hbom->SetTrackEtaWindow(kTrackEtaWindowCl);
   }
   else if (typeRec.Contains("AODMC2")){
     hbom->SetTrackTypeRec(AliAnalysisTaskJetHBOM::kTrackAODMCCharged);
     hbom->SetTrackPtCut(kPtTrackCutCl);
     hbom->SetTrackEtaWindow(5);
   }
   else if (typeRec.Contains("AODMC")){
     hbom->SetTrackTypeRec(AliAnalysisTaskJetHBOM::kTrackAODMCAll);
     hbom->SetTrackPtCut(kPtTrackCutCl);
     hbom->SetTrackEtaWindow(5);
   }
   else if (typeRec.Contains("AODextraonly")) {
     hbom->SetTrackTypeRec(AliAnalysisTaskJetHBOM::kTrackAODextraonly);
     hbom->SetTrackPtCut(kPtTrackCutCl);
     hbom->SetTrackEtaWindow(kTrackEtaWindowCl);
   }
   else if (typeRec.Contains("AODextra")) {
     cout << "AliAnalysisTaskJetHBOM::kTrackAODextra: " << AliAnalysisTaskJetHBOM::kTrackAODextra << endl;
     hbom->SetTrackTypeRec(AliAnalysisTaskJetHBOM::kTrackAODextra);
     hbom->SetTrackPtCut(kPtTrackCutCl);
     hbom->SetTrackEtaWindow(kTrackEtaWindowCl);
   }
   else if (typeRec.Contains("AOD")) {
     hbom->SetTrackTypeRec(AliAnalysisTaskJetHBOM::kTrackAOD);
     hbom->SetTrackPtCut(kPtTrackCutCl);
     hbom->SetTrackEtaWindow(kTrackEtaWindowCl);
   }

   hbom->SetRparam(radius);
   hbom->SetGhostArea(0.005);
   hbom->SetGhostEtamax(kTrackEtaWindowCl);

   switch (jf) {
   case "ANTIKT":
     hbom->SetAlgorithm(2); // antikt from fastjet/JetDefinition.hh
     break;
   case "CA":
     hbom->SetAlgorithm(1); // CA from fastjet/JetDefinition.hh
     break;
   case "KT":
     hbom->SetAlgorithm(0); // kt from fastjet/JetDefinition.hh
     break;
   default:
     ::Error("AddTaskJetHBOM", "Wrong jet finder selected\n");
     return 0;
   }

   //Constant Cone analysis
   if(constCone){
     hbom->SetRandConePos(constConeEta,constConePhi);
   }

   
   if(kWriteAOD){
     if(outputFile.Length())hbom->SetJetOutputFile(outputFile);
     hbom->SetJetOutputBranch(Form("hbom%s_%s%s",bRec,jf,cAdd.Data()));
     hbom->SetJetOutputMinPt(0); // store only jets above a certain threshold
   }

   //sets number of detector hits
   hbom->SetfNHBOM(fNHBOM);

   //physics Selection
   if(iPhysicsSelectionFlag)hbom->SelectCollisionCandidates(iPhysicsSelectionFlag);

   mgr->AddTask(hbom);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   AliAnalysisDataContainer *coutput1_Spec = mgr->CreateContainer(Form("hbom_%s_%s_%s%s",bRec,bGen,jf,cAdd.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:PWGJE_hbom_%s_%s_%s%s",AliAnalysisManager::GetCommonFileName(),bRec,bGen,jf,cAdd.Data()));

   mgr->ConnectInput  (hbom, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (hbom, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (hbom,  1, coutput1_Spec );
   
   //loads efficiencies
   hbom->SetEfficiencyPt(GetEfficiencyPt(effLoc));
   hbom->SetEfficiencyPhi(GetEfficiencyPhi(effLoc));
   
   return hbom;
}

//loads single track pT efficiency from root file
TH1F *GetEfficiencyPt(TString effLoc){
  TFile *fIn = 0;
  TH1F *hEffPt = 0; 

  if(!fIn)fIn = TFile::Open(effLoc.Data());
  if(!fIn)Printf("%s%d no input data",(char*)__FILE__,__LINE__);
  if(!hEffPt)hEffPt = (TH1F*)fIn->Get("hSingleTrackEffPt");
  if(!hEffPt)Printf("%s%d no single track efficiency spectrum available",(char*)__FILE__,__LINE__);
  gROOT->cd();

  TH1F *hEffPtClone = (TH1F*)hEffPt->Clone(hEffPt->GetName()); 
  fIn->Close();
  return hEffPtClone;


}

//loads phi-pT efficiency from root file
TH2D *GetEfficiencyPhi(TString effLoc){
  TFile *fIn = 0;
  TH2D *hPhiPt = 0; 

  if(!fIn)fIn = TFile::Open(effLoc.Data());
  if(!fIn)Printf("%s%d no input data",(char*)__FILE__,__LINE__);
  if(!hPhiPt)hPhiPt = (TH2D*)fIn->Get("h2TrackPtPhiNorm");
  if(!hPhiPt) cout<<"Could not load h2TrackPtPhiNorm"<<endl; 
  if(!hPhiPt)Printf("%s%d no phi-pt efficiency spectrum available",(char*)__FILE__,__LINE__);

  gROOT->cd();
  TH2D *hPhiPtClone = (TH2D*)hPhiPt->Clone(hPhiPt->GetName()); 
  fIn->Close();
  return hPhiPtClone;

}
