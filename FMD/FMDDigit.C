void FMDDigit () 
{
// Dynamically link some shared libs

   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   }

  AliRunDigitizer * manager = new AliRunDigitizer(1,1);
  manager->SetInputStream(0,"galice.root");
  AliFMDDigitizer *FMD = new AliFMDDigitizer(manager);
  manager->Exec("");
}














