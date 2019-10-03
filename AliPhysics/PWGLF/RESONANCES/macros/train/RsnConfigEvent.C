//
// Configuration script for monitor task with 2010 runs
//
// It contains a class definition where the cuts for each object
// are defined separately, functions are initialized and so on.
// This is used in the main function (named after the file name),
// which is called by the 'AddTask' function.
//

Bool_t RsnConfigEvent
(
   AliRsnAnalysisTask *task,
   Bool_t              isMC,
   Bool_t              useCentrality,
   AliRsnCutSet       *eventCuts
)
{
   if (!task) myError("NULL task");
   
   // ----------------------------------------------------------------------------------------------
   // -- LOOP DEFINITION ---------------------------------------------------------------------------
   // ----------------------------------------------------------------------------------------------
   
   // loop on events
   AliRsnLoopEvent *loop = new AliRsnLoopEvent("evtLoop");
   
   // add cuts to loop
   loop->SetEventCuts(eventCuts);
   
   // add loop to task
   task->Add(loop);
   
   // ----------------------------------------------------------------------------------------------
   // -- OUTPUTS DEFINITION ------------------------------------------------------------------------
   // ----------------------------------------------------------------------------------------------
   
   Double_t mult[] = { 0.,  1.,  2.,  3.,  4.,  5.,   6.,   7.,   8.,   9.,  10., 11., 12., 13., 
                      14., 15., 16., 17., 18., 19.,  20.,  21.,  22.,  23.,  24.,  25., 30., 35., 
                      40., 50., 60., 70., 80., 90., 100., 120., 140., 160., 180., 200., 500.};
   Int_t    nmult  = sizeof(mult) / sizeof(mult[0]);

   // axes
   AliRsnValueStd *axisMultESD  = new AliRsnValueStd("MESD", AliRsnValueStd::kEventMultESDCuts , nmult, mult);
   AliRsnValueStd *axisMultSPD  = new AliRsnValueStd("MSPD", AliRsnValueStd::kEventMultSPD     , nmult, mult);
   AliRsnValueStd *axisMultMC   = new AliRsnValueStd("MMC" , AliRsnValueStd::kEventMultMC      , nmult, mult);
   AliRsnValueStd *axisMultTRK  = new AliRsnValueStd("MTRK", AliRsnValueStd::kEventMult        , nmult, mult);
   AliRsnValueStd *axisCentV0   = new AliRsnValueStd("CNT" , AliRsnValueStd::kEventCentralityV0, 0.0  , 100.0, 1.0);
   AliRsnValueStd *axisMultTest = new AliRsnValueStd("TEST", AliRsnValueStd::kEventMult        , 1, 1.0 , 1E10);
   
   // create output
   AliRsnListOutput *out = new AliRsnListOutput("evtOut", AliRsnListOutput::kHistoDefault);
   
   // always add test axis, used to know how many events had at least one track
   out->AddValue(axisMultTest);
   
   // add values to output:
   // centrality only if this is requested, otherwise multiplicity
   if (useCentrality) {
      ::Info("RsnConfigEvent", "Adding centrality");
      out->AddValue(axisCentV0);
   } else {
      ::Info("RsnConfigEvent", "Adding multiplicity");
      //out->AddValue(axisMultESD);
      //out->AddValue(axisMultSPD);
      out->AddValue(axisMultTRK);
      if (isMC) {
         out->AddValue(axisMultMC);
      }
   }
   
   // add output to loop
   loop->AddOutput(out);
   
   return kTRUE;
}
