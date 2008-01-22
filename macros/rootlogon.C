void rootlogon() {
   printf("\nWELCOME to ALICE\n\n");
   TInterpreter * mycint = gROOT->GetInterpreter();
   mycint->AddIncludePath("$ALICE_ROOT/include");
}

