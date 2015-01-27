void
CopySec(Bool_t      fmd=true, 
	ULong_t     tgtRun=0, 
	UShort_t    tgtSys=2,
	UShort_t    tgtSNN=2760, 
	Short_t     tgtFld=-5, 
	ULong_t     srcRun=118506, 
	UShort_t    srcSys=1,
	UShort_t    srcSNN=900, 
	Short_t     srcFld=+5)
{

  const char* fwd = "${ALICE_PHYSICS}/PWGLF/FORWARD/analysis2";
  gROOT->LoadMacro(Form("%s/corrs/CopyCorr.C", fwd));

  CopyCorr(fmd, "secondary", 
	   tgtRun, tgtSys, tgtSNN, tgtFld, 
	   srcRun, srcSys, srcSNN, srcFld, false);
}
//
// EOF
//







