void rootlogon() {
   printf("\nWELCOME to ALICE\n\n");
   TInterpreter * mycint = gROOT->GetInterpreter();
   mycint->AddIncludePath("$ALICE_ROOT/include");
   mycint->AddIncludePath("$ALICE_ROOT/PYTHIA6");
   mycint->AddIncludePath("$ROOTSYS/include");
   mycint->AddIncludePath("$ALICE_ROOT");
   mycint->AddIncludePath("$ALICE_ROOT/HLT/BASE");
   mycint->AddIncludePath("$ALICE_ROOT/HLT/BASE/util");
   mycint->AddIncludePath("$ALICE_ROOT/HLT/JET");
   mycint->AddIncludePath("$ALICE_ROOT/HLT/JET/cone");

   gStyle->SetPalette(1);	
   gROOT->SetStyle("Plain");
}
