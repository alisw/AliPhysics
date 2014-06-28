void AddTaskAODPWG4HighPtSpectra(char *prodType = "LHC12a15e",
				 Bool_t isPbPb = kFALSE, 
				 UInt_t triggerMask = AliVEvent::kAny,
				 Bool_t bSelectHijingParticles = kFALSE,
				 Bool_t usePythiaxsec = kTRUE)
{
  gROOT->LoadMacro(gSystem->ExpandPathName("$ALICE_ROOT/PWGJE/macros/AddTaskPWG4HighPtSpectra.C"));
  AddTaskPWG4HighPtSpectraQA_AOD(prodType, isPbPb, triggerMask, bSelectHijingParticles, usePythiaxsec);
}
