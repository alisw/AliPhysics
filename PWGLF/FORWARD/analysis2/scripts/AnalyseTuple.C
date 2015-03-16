void
AnalyseTuple(Bool_t proof=true, Long64_t maxEvents=-1, const char* opt="")
{
  const char* fwd = "$(ALICE_PHYSICS)/PWGLF/FORWARD/analysis2";
  gSystem->AddIncludePath("-I${ANA_SRC} -I${ALICE_PHYSICS}/include");
  gROOT->Macro(Form("%s/scripts/LoadLibs.C",fwd));
  gROOT->LoadMacro(Form("%s/scripts/TupleSelector.C+%s",fwd,opt));

  if (proof) TupleSelector::Proof(maxEvents, opt);
  else       TupleSelector::Run(maxEvents);
}
