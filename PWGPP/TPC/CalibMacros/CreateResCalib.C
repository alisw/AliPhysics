void CreateResCalib(int run=245231
		    ,Long64_t tmin=0
		    ,Long64_t tmax= 9999999999
		    ,const char * resList="lst.txt"
		    ,const char * cdb="raw://"
		    ,int maxTracks = 5000000
		    ,Bool_t useTOFBC = kFALSE
	  )
{
  TString useTOFBCStr = gSystem->Getenv("useTOFBC");    
  if (!useTOFBCStr.IsNull()) {
    useTOFBCStr.ToLower();
    if      (useTOFBCStr.Contains("true")) useTOFBC = kTRUE;
    else if (useTOFBCStr.Contains("false")) useTOFBC = kFALSE;
    else {
      ::Info("AliTPCcalibAlignInterpolation::FitDrift",
	     "Cannot interpret Env.var {useTOFBC}=%s as true or false, keep default setting %d",gSystem->Getenv("useTOFBC"),useTOFBC);
    }
  }
  //
  AliTPCDcalibRes* clb = new AliTPCDcalibRes(run, tmin, tmax, resList);
  clb->SetOCDBPath(cdb);
  if (maxTracks>0) clb->SetMaxTracks(maxTracks);
  clb->SetUseTOFBC(useTOFBC);
  clb->ProcessFromDeltaTrees();
  clb->Save();
}

void postProcResCalib(int run=245231
	      ,Long64_t tmin=0
	      ,Long64_t tmax= 9999999999
	      )
{

  AliTPCDcalibRes* clb = new AliTPCDcalibRes(run, tmin, tmax);
  clb->SetMaxTracks(5000000);
  clb->SetOCDBPath("local:///cvmfs/alice.cern.ch/calibration/data/2015/OCDB");
  clb->ProcessFromLocalBinnedTrees();
  clb->Save();
}
