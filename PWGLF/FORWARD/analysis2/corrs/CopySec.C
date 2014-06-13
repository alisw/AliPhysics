void
CopySec(Bool_t fmd=true, UShort_t sys=2, UShort_t sNN=2760, Short_t fld=-5, ULong_t run=0)
{
  const char* fwd = "${ALICE_ROOT}/PWGLF/FORWARD/analysis2";
  if (!gROOT->GetClass("AliOADBForward"))
    gROOT->Macro(Form("%s/scripts/LoadLibs.C", fwd));

  TString   table   = "secondary";
  ULong_t   oldRun  = 118506;
  UShort_t  oldSys  = 1; // pp
  UShort_t  oldSNN  = 900;
  Short_t   oldFld  = +5;
  const char* det   = fmd ? "fmd" : "spd";
  AliOADBForward* db = new AliOADBForward;
  db->Open(Form("%s_corrections.root", det), table, true);

  if (!db->CopyEntry(table, 
		     oldRun, oldSys, oldSNN, oldFld, 
		     run,    sys,    sNN,    fld,    
		     false, false)) 
    Warning("CopySec", "Failed to copy %s %lu/%hu/%hu/%hd -> %lu/%hu/%hu/%hd",
	    det,
	    oldRun, oldSys, oldSNN, oldFld, 
	    run,    sys,    sNN,    fld);
}
//
// EOF
//







