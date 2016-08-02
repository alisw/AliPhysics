// -*- C++ -*-
// $Id$


void MakeLumiRegion_4269() {
  gROOT->LoadMacro("AliLuminousRegionFit.cxx+");
  gROOT->LoadMacro("Util.C");

  const ScanData sd[] = {  
    { 4269,
      11,
      "root/4269/AnalysisResults_234040_allparams.root",
      "txt/4269/SepVStime.out",
      "Scan1",
      "Scan1X", 0, 1440534260, 1440535435, 0,
      "Scan1Y", 1, 1440535779, 1440536959, 0
    },
    
    { 4269,
      11,
      "root/4269/AnalysisResults_234040_allparams.root",
      "txt/4269/SepVStime.out",
      "Scan2",
      "Scan2X", 0, 1440537438, 1440538613, 0,
      "Scan2Y", 1, 1440538946, 1440540071, 0
    },
    
    { 4269,
      11,
      "root/4269/AnalysisResults_234045_allparams.root",
      "txt/4269/SepVStime.out",
      "ScanOffset",
      "ScanOffsetX", 0, 1440543600, 1440544900, 336,
      "ScanOffsetY", 1, 1440545550, 1440546900, 372
    }
  };

  const Int_t n = sizeof(sd)/sizeof(ScanData);

  // fill 4269
  const Int_t bcs[9] = {
    -1,   //all
    1465, //BCM5  1465H1L2098H
    1625, //BCM6  1625H1L1938H
    1785, //BCM7  1785H1L1778H
    1945, //BCM8  1945H1L1618H
    2065, //BCM9  2065H1L1498H
    3067, //BCM10 3067H1L496H
    3267, //BCM11 3267H1L296H
    3467  //BCM12 3467H1L96H
  };
  

  for (Int_t i=0; i<n; ++i) {
    CheckCopyFile(sd[i].vtxFileName);
    for (Int_t j=0; j<1; ++j) {
      Int_t bcid = bcs[j];
      AliLuminousRegionFit f(sd[i].fillNumber,
			     sd[i].minNumberOfTracks,
			     sd[i].vtxFileName,
			     sd[i].sepFileName);
      f.DoFit(sd[i].scanName1, sd[i].t1, sd[i].t2, sd[i].scanType1, sd[i].offset1, bcid);
      f.DoFit(sd[i].scanName2, sd[i].t3, sd[i].t4, sd[i].scanType2, sd[i].offset2, bcid);
    }
  }
}
