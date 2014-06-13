AliAnalysisTaskHJetSpectra* AddTaskHJetSpectra(
  Double_t            jetRadius               = 0.4,  //radius of analyzed jets
  Double_t            jetRadiusBg             = 0.3,  //radius of jets to be removed when estimating cell median bg
  Int_t               trigger                 = AliVEvent::kINT7,  //trigger
  Int_t               isMC                    = 0,  //MC flag
  Double_t            randomConeR             = 0.4,  //random cone for deltaPt + perp cone bg
  const char*         containerSuffix         = "",   //tag to the name of container
  const char*         usedTracks              = "PicoTracks",  //tracks
  const char*         centralityType          = "V0A",   //centrality
  Double_t            trackEtaWindow          = 0.9,   //pseudorapidity range for tracks
  Bool_t              useVertexCut            = kTRUE,  // vertex cut
  Bool_t              usePileUpCut            = kTRUE, // discard pile up event
  Int_t               numberOfCentralityBins  = 1,     // the number of centrality bins
  Double_t            ttLow                   = 8.0,      // trigger hardron low pT
  Double_t            ttHigh                  = 50.0,     // trigger hadron high pT
  Int_t               ttType                  = 0,        // 0= single inclusive hadron trigger, else inclusive hadron trigger
  Double_t            dphi                    = 0.6, // |Delta phi_jet, trigger|< pi-0.6
  Bool_t              binning                 = 0    //binning of jet histograms 0=2GeV width  1=1GeV width
){


   // #### Detect the demanded trigger with its readable name
   TString triggerName(Form("Trigger_%i", trigger));
   if(trigger == AliVEvent::kAnyINT)
      triggerName = "kAnyINT";
   else if(trigger == AliVEvent::kAny)
      triggerName = "kAny";
   else if(trigger == AliVEvent::kINT7)
      triggerName = "kINT7";
   else if(trigger == AliVEvent::kMB)
      triggerName = "kMB";
   else if(trigger == AliVEvent::kEMC7)
      triggerName = "kEMC7";
   else if(trigger == AliVEvent::kEMCEJE)
      triggerName = "kEMCEJE";
   else if(trigger == AliVEvent::kEMCEGA)
      triggerName = "kEMCEGA";
 
   // #### DEFINE MANAGER AND DATA CONTAINER NAMES
   AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
   if (!manager) {
     ::Error("AddTaskHJetSpectra.C", "No analysis manager to connect to.");
     return NULL;
   }
 
   TString containerNameSuffix("");
   if(strcmp(containerSuffix,""))  containerNameSuffix = Form("_%s", containerSuffix);
 
   TString myContName("");
   if(isMC){
      myContName = Form("AnalysisR0%2.0f_%s_MC%s_Dphi%02d_T%d_Ptt%d_%d", 
         jetRadius*100, triggerName.Data(), containerNameSuffix.Data(),
         TMath::Nint(10*dphi), ttType, TMath::Nint(ttLow), TMath::Nint(ttHigh));
   }else{
      myContName = Form("AnalysisR0%2.0f_%s%s_Dphi%02d_T%d_Ptt%d_%d",    
         jetRadius*100, triggerName.Data(), containerNameSuffix.Data(),
         TMath::Nint(10*dphi), ttType, TMath::Nint(ttLow), TMath::Nint(ttHigh));
   }

   // #### ADD NECESSARY JET FINDER TASKS
   enum AlgoType {kKT, kANTIKT};
   enum JetType  {kFULLJETS, kCHARGEDJETS, kNEUTRALJETS};

   if(jetRadius < 0.1 || jetRadiusBg < 0.1) return NULL;
 
   gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
   AliEmcalJetTask* jetFinderTask   = AddTaskEmcalJet(usedTracks,"",kANTIKT,jetRadius,  kCHARGEDJETS,0.150,0.300); //FK//
   AliEmcalJetTask* jetFinderTaskBg = AddTaskEmcalJet(usedTracks,"",kANTIKT,jetRadiusBg,kCHARGEDJETS,0.150,0.300); //FK//excl from bg

   // #### DEFINE EXTERN CMS RHO TASK
   TString myRhoName("ExternalRhoTask");
   AliEmcalJetTask* jetFinderRho   = AddTaskEmcalJet(usedTracks,"", kANTIKT, 0.4, kCHARGEDJETS,0.150,0.300); // anti-kt
   AliEmcalJetTask* jetFinderRhoKT = AddTaskEmcalJet(usedTracks,"", kKT,     0.4, kCHARGEDJETS,0.150,0.300); // kt
   jetFinderRhoKT->SetMinJetPt(0);

   gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskRhoSparse.C");
   AliAnalysisTaskRhoSparse* rhotask = AddTaskRhoSparse(jetFinderRhoKT->GetName(), 
                                                        jetFinderRho->GetName(), 
                                                        usedTracks,   //pico trakcs
                                                                "",   //calo clusters
                                                        myRhoName.Data(), 
                                                               0.4,  //jet radius
                                                             "TPC",  //cut type
                                                                0.,  //jet area cut
                                                               15.,  //jet pt cut ????????
                                                                 0,  //enareacut 
                                                                 0,  //sfunc
                                                                 1,  //excl Jets  //FK// ????????
                                                            kFALSE,   //no histo
                                                   myRhoName.Data(),  //task name
                                                             kTRUE); //claculate rho CMS

   // #### DEFINE ANALYSIS TASK
   AliAnalysisTaskHJetSpectra *task = new AliAnalysisTaskHJetSpectra(
                                           Form("HJetSpectra_%s_%s_TT", jetFinderTask->GetName(), triggerName.Data()), 
                                           usedTracks, 
                                           jetFinderTask->GetName(), 
                                           jetFinderTaskBg->GetName());
 
   // #### Task preferences
   task->SetAcceptanceWindows(trackEtaWindow, jetRadius, jetRadiusBg);
   task->SetUsePileUpCut(usePileUpCut);
   task->SetUseDefaultVertexCut(useVertexCut);
   task->SetSignalJetMinArea(0.6*jetRadius*jetRadius*TMath::Pi()); //To same puziva Mata
   task->SetRandConeRadius(randomConeR);
   task->SelectCollisionCandidates(trigger);
   task->SetCentralityType(centralityType); 
   task->SetNumberOfCentralityBins(numberOfCentralityBins);
   task->SetExternalRhoTaskName(myRhoName.Data());
 
   task->SetTT(ttLow, ttHigh);
   task->SetTTType(ttType); 
   task->SetDphi(dphi); // |Delta phi_jet, trigger|< pi-0.6
   task->SetDoubleBinPrecision(binning); 
   task->SetMinPtOfJetsToBeRemovedInBg(15.0); 
   //task->SetMC(isMC);
   task->SetNofRandomCones(1);

   if(isMC) task->SetAnalyzeMC(isMC);

   // output container
   contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:ChJetSpectra%s", AliAnalysisManager::GetCommonFileName(), myContName.Data()));
 
 
   // #### ADD ANALYSIS TASK
   manager->AddTask(task);
   manager->ConnectInput(task, 0, manager->GetCommonInputContainer());
   manager->ConnectOutput(task, 1, contHistos);
 
   return task;
}
