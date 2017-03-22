// -*- C++ -*-

// this macro tests all STARLIGHT processes defined in $ALIDPG_ROOT/MC/GeneratorConfig.C
// for this the $ALIDPG_ROOT environment variable needs to be set
void testsl()
{
  gSystem->Exec("cat $ALIDPG_ROOT/MC/GeneratorConfig.C "
		"| awk '/STARLIGHT_PROCESSES_BEGIN/ {p=1; getline}"
		"       /STARLIGHT_PROCESSES_END/   {p=0} "
		"       // {if (p) {print $2}} '"
		" > sl_proc.txt");

  TTree t;
  t.ReadFile("sl_proc.txt", "process_name/C");
  Char_t process_name[1024];
  t.SetBranchAddress("process_name", &process_name);

  for (Int_t i=0, n=t.GetEntries(); i<n; ++i) {
    t.GetEntry(i);
    Printf("%s", process_name);
    // if (!TString(process_name).Contains("RhoPrime"))
    //   continue;
    gSystem->Exec(Form("aliroot -b -q .x simStarlight.C'(\"%s\")' | tee log.txt 2>&1", process_name));
  }

  t.ResetBranchAddresses();
}
