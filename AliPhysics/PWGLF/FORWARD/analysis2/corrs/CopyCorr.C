Bool_t
CopyCorr(Bool_t      fmd, 
	 const char* table,
	 ULong_t     tgtRun, 
	 UShort_t    tgtSys,
	 UShort_t    tgtSNN, 
	 Short_t     tgtFld, 
	 ULong_t     srcRun, 
	 UShort_t    srcSys,
	 UShort_t    srcSNN, 
	 Short_t     srcFld,
	 Bool_t      mc)
{
  const char* fwd = "${ALICE_PHYSICS}/PWGLF/FORWARD/analysis2";
  if (!gROOT->GetClass("AliOADBForward"))
    gROOT->Macro(Form("%s/scripts/LoadLibs.C", fwd));

  TString tab(table);
  const char*  det        = fmd ? "fmd" : "spd";
  const char*  possible[] = { "elossfits", 
			     "secondary", 
			     "noisegain", 
			     "acceptance", 
			     "merging", 
			     "vertexbias", 
			     "doublehit",
			     0 };
  const char** pTest      = possible;
  while (*pTest) { 
    if (tab.EqualTo(*pTest)) break;
    pTest++;
  }
  if (!(*pTest)) { 
    Warning("CopyCorr", "Unknown table: %s", table);
    return false;
  }
	 
  if (tgtSys < 1)    tgtSys = srcSys;
  if (tgtSNN < 1)    tgtSNN = srcSNN;
  if (tgtFld >= 999) tgtFld = srcFld;

  AliOADBForward* db = new AliOADBForward;
  db->Open(Form("%s_corrections.root", det), tab, true);

  if (!db->CopyEntry(tab, 
		     srcRun, srcSys, srcSNN, srcFld, 
		     tgtRun, tgtSys, tgtSNN, tgtFld,    
		     mc, false)) {
    Error("CopySec", 
	  "Failed to copy %s %lu/%hu/%hu/%hd/%s -> %lu/%hu/%hu/%hd/%s",
	  det,
	  srcRun, srcSys, srcSNN, srcFld, (mc ? "MC" : "Real"), 
	  tgtRun, tgtSys, tgtSNN, tgtFld, (mc ? "MC" : "Real"));
    return false;
  }
  return true;
}
