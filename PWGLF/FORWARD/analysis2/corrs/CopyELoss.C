void
CopyELoss(ULong_t     tgtRun, 
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
  gROOT->LoadMacro(Form("%s/corrs/CopyCorr.C", fwd));

  CopyCorr(true, "elossfits", 
	   tgtRun, tgtSys, tgtSNN, tgtFld, 
	   srcRun, srcSys, srcSNN, srcFld, mc);
}
//
// EOF
//







