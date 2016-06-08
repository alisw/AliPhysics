// -*- C++ -*-
// $Id$

struct ScanData {
  Int_t fillNumber;
  Int_t minNumberOfTracks;
  const char* vtxFileName;
  const char* sepFileName;
  const char* scanName;

  const char* scanName1;
  Int_t scanType1;
  Double_t t1,t2;
  Double_t offset1;

  const char* scanName2;
  Int_t scanType2;
  Double_t t3,t4;
  Double_t offset2;
} ;

const ScanData sd[] = {  
  { 4634,
    11,
    "root/4634/AnalysisResults_VdM_244369.root",
    "txt/4634/SepVStime.out",
    "Scan1",
    "Scan1X", 0, 1447981655, 1447982680, 0,
    "Scan1Y", 1, 1447982863, 1447983886, 0                             
  },

  { 4634,
    11,
    "root/4634/AnalysisResults_VdM_244369.root",
    "txt/4634/SepVStime.out",
    "Scan2",
    "Scan2X", 0, 1447984119, 1447985146, 0,
    "Scan2Y", 1, 1447985327, 1447986351, 0
  },
  
  { 4634,
    11,
    "root/4634/AnalysisResults_VdM_244375.root",
    "txt/4634/SepVStime.out",
    "ScanOffset",
    "ScanOffsetX", 0, 1447988540, 1447989542, -457.5,
    "ScanOffsetY", 1, 1447989790, 1447990716, -457.5
  }
};

const Int_t n = sizeof(sd)/sizeof(ScanData);

void MakeLumiRegion_4634() {
  gROOT->LoadMacro("AliLuminousRegionFit.cxx+");

  for (Int_t i=0; i<n; ++i) {
    AliLuminousRegionFit f(sd[i].fillNumber,
			   sd[i].minNumberOfTracks,
			   sd[i].vtxFileName,
			   sd[i].sepFileName);
    f.DoFit(sd[i].scanName1, sd[i].t1, sd[i].t2, sd[i].scanType1, sd[i].offset1, -1);
    f.DoFit(sd[i].scanName2, sd[i].t3, sd[i].t4, sd[i].scanType2, sd[i].offset2, -1);
  }
}
