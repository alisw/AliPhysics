// Simulate events and merge with other events (simulated or real)
// at the level of sdigits.
// Input: sdigits of background file

// Full exercise can be done automatically with other alice systems
// in the examples in $ALICE_ROOT/test/merge

// Author: GCB (from example in $ALICE_ROOT/test/merge)

void simmerge() 
{
  
  AliSimulation sim;
  sim.SetRunSimulation(kTRUE);
  sim.SetMakeSDigits("EMCAL");
  sim.SetMakeDigits("EMCAL");
  sim.MergeWith("bkg/galice.root",1);

  sim.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  sim.SetSpecificStorage("GRP/GRP/Data",
			       Form("local://%s",gSystem->pwd()));
  sim.SetRunQA(":") ; 
  AliQA::SetQARefStorage("local://$ALICE_ROOT/OCDB") ;
  
  for (Int_t det = 0 ; det < AliQA::kNDET ; det++) {
    sim.SetQACycles(det, 1) ;
  }

  sim.Run(10);
}
