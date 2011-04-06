//
// *** Configuration script for phi->KK analysis with 2010 runs ***
// 
// A configuration script for RSN package needs to define the followings:
//
// (1) decay tree of each resonance to be studied, which is needed to select
//     true pairs and to assign the right mass to all candidate daughters
// (2) cuts at all levels: single daughters, tracks, events
// (3) output objects: histograms or trees
//
Bool_t RsnConfigMonitorTPC
(
   AliRsnAnalysisTask *task,
   Bool_t              isMC,
   Bool_t              useCentrality,
   AliRsnCutSet       *eventCuts
)
{
   void myError  (const char *msg) {::Error  ("RsnConfigPhi", msg);}
   void myWarning(const char *msg) {::Warning("RsnConfigPhi", msg);}
   void myInfo   (const char *msg) {::Info   ("RsnConfigPhi", msg);}

   if (!task) myError("NULL task");
   
   const char *suffix = "tpcstd";
      
   // ==================================================================================================================
   // == DEFINITIONS ===================================================================================================
   // ==================================================================================================================
   
   // daughter definitions
   AliRsnDaughterDef *tracks = new AliRsnDaughterDef(AliRsnDaughter::kTrack, 0);

   // computation objects:
   // these are the objects which drive the computations, whatever it is (histos or tree filling)
   // and all tracks passed to them will be given the mass according to the reference pair definition
   // we create two unlike-sign pair computators, one for all pairs and another for true pairs (useful in MC)
   AliRsnLoopDaughter *loop[3];
   loop[0] = new AliRsnLoopDaughter(Form("%s_quality", suffix), 0, tracks);
   loop[1] = new AliRsnLoopDaughter(Form("%s_pion"   , suffix), 0, tracks);
   loop[2] = new AliRsnLoopDaughter(Form("%s_kaon"   , suffix), 0, tracks);

   // assign the ID of the entry lists to be used by each pair to get selected daughters
   // in our case, the AliRsnInputHandler contains only one list for selecting kaons
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   AliMultiInputEventHandler *multi = dynamic_cast<AliMultiInputEventHandler*>(mgr->GetInputEventHandler());
   if (!multi) {
      myError("Needed a multi input handler!");
      return kFALSE;
   }
   TObjArray *array = multi->InputEventHandlers();
   AliRsnInputHandler *rsn = (AliRsnInputHandler*)array->FindObject("rsnInputHandler");
   if (!rsn) {
      myError("Needed an RSN event handler");
      return kFALSE;
   }
   AliRsnDaughterSelector *sel = rsn->GetSelector();
   Int_t id[3];
   id[0] = sel->GetID("qualityTPC", kTRUE);
   id[1] = sel->GetID("pionTPC", kTRUE);
   id[2] = sel->GetID("kaonTPC", kTRUE);
   
   // ==================================================================================================================
   // == COMPUTED VALUES & OUTPUTS =====================================================================================
   // ==================================================================================================================
   
   // axes
   AliRsnValueStd *axisMomTPC = new AliRsnValueStd("pTPC", AliRsnValueStd::kTrackPtpc     , 0.0,   5.0, 0.01 );
   AliRsnValueStd *axisSigTPC = new AliRsnValueStd("sTPC", AliRsnValueStd::kTrackTPCsignal, 0.0, 500.0, 2.0  );
   AliRsnValueStd *axisMomTOF = new AliRsnValueStd("pTOF", AliRsnValueStd::kTrackP        , 0.0,   5.0, 0.01 );
   AliRsnValuePID *axisSigTOF = new AliRsnValuePID("sTOF", AliRsnValuePID::kTOFnsigma     , AliPID::kKaon, -200.0, 200.0, 10.0);
   
   // create outputs
   AliRsnListOutput *outTPC = new AliRsnListOutput("outTPC", AliRsnListOutput::kHistoDefault);
   AliRsnListOutput *outTOF = new AliRsnListOutput("outTOF", AliRsnListOutput::kHistoDefault);
   
   // add values to outputs
   outTPC->AddValue(axisMomTPC);
   outTPC->AddValue(axisSigTPC);
   outTOF->AddValue(axisMomTOF);
   outTOF->AddValue(axisSigTOF);
   
   // ==================================================================================================================
   // == CONCLUSION ====================================================================================================
   // ==================================================================================================================
   
   for (Int_t i = 0; i < 3; i++) {
      // list ID
      if (id[i] < 0) {
         myError(Form("Required entry list for loop #%d is not added in the selector", i));
         return kFALSE;
      }
      myInfo(Form("Required entry list for loop #%d [%s] is %d", i, loop[i]->GetName(), id[i]));
      loop[i]->SetListID(id[i]);
      
      // event cuts
      loop[i]->SetEventCuts(eventCuts);
      
      // outputs
      loop[i]->AddOutput(outTPC);
      loop[i]->AddOutput(outTOF);
      
      // add to task
      task->Add(loop[i]);
   }
   
   return kTRUE;
}
