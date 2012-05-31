void T0Digit () 
{
// Dynamically link some shared libs

   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   }

  AliRunDigitizer * manager = new AliRunDigitizer(1,1);
  manager->SetInputStream(0,"galice.root");
  //  manager->SetOutputFile("digits.root");
  AliT0Digitizer *T0 = new AliT0Digitizer(manager);
  manager->Exec("");
}














