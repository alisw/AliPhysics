void CreateResCalib(int run=245231
	  ,Long64_t tmin=0
	  ,Long64_t tmax= 9999999999
	  ,const char * resList="lst.txt"	 
	  )
{

  AliTPCDcalibRes* clb = new AliTPCDcalibRes(run, tmin, tmax, resList);
  clb->SetOCDBPath("local:///cvmfs/alice.cern.ch/calibration/data/2015/OCDB");
  clb->ProcessFromDeltaTrees();
  clb->Save();
}

void postProcResCalib(int run=245231
	      ,Long64_t tmin=0
	      ,Long64_t tmax= 9999999999
	      )
{

  AliTPCDcalibRes* clb = new AliTPCDcalibRes(run, tmin, tmax);
  clb->SetOCDBPath("local:///cvmfs/alice.cern.ch/calibration/data/2015/OCDB");
  clb->ProcessFromLocalBinnedTrees();
  clb->Save();
}
