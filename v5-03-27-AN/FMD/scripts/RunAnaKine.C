void
RunAnaKine(Bool_t primary_only=false, Bool_t segmented=true)
{
  // Below for v2=0.06, hijing para=10000
  //  no secondaries/segs | false  | true     
  //  --------------------+--------+-------
  //  false               | low    | low
  //  true                | ok     | ok
  //  
  // Below for v2=0.05, hijing para=10000
  //  no secondaries/segs | false  | true     
  //  --------------------+--------+-------
  //  false               | low    | low
  //  true                | ok     | ok
  //  
  // gROOT->LoadMacro("Compile.C");
  // gSystem->Unload("AliFMDAnaFlowKine_C.so");
  // gSystem->Unload("AliFMDAnaFlowRing_h.so");
  if (!Compile("AliFMDAnaFlowRing.h")) return;
  if (!Compile("AliFMDAnaFlowKine.C")) return;
  AliFMDAnaFlowKine ak(1, segmented, primary_only); 
  ak.Run(); 
  TBrowser b;
  b.Add(&ak);
}
