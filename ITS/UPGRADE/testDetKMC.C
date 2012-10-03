#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TString.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include "KMCDetector.h"
#endif

TObjArray* CreateTrackConditions(const char* nm);
TObjArray  summary;


TObject* gDet=0;

// ptBins to generate
//double pts[] = {0.15,0.175,0.2,0.225,0.25,0.3,0.4,0.5,0.8,1.1,1.6,2,4,8,13,20};
double pts[] = {0.2};
//double pts[] = {0.2, 1., 10.};

void testDetKMC(int nev=1000,  // n events to generate per bin
		int setVer=0,   // optional version to pass for detector building
		double mass=0.14,  // particle mass to generate
		double eta=0.45,   // eta to test
		Bool_t aliceNew=kFALSE,
		Bool_t tpc=kTRUE, Double_t eff=0.95,  // detector settings
		Double_t vtxConstr=-1., // use vertex constraint in the fit
		const char* sumOut="trsumDef.root" // output file name
		)
{
  if (!gROOT->GetClass("KMCDetector")) gROOT->LoadMacro("KMCDetector.cxx+");
  //
  // Detector, tracking initialization >>>
  //
  if (vtxConstr>0) KMCDetector::SetVtxConstraint(vtxConstr,vtxConstr);
  KMCLayer::SetDefEff(eff);
  //
  KMCDetector *itsdet = new KMCDetector((char*)"ALICE",(char*)"ITS");
  gDet = itsdet;
  KMCDetector &its = *itsdet;
  its.SetAvgRapidity(eta);
  if (aliceNew)     its.MakeAliceAllNew(tpc,1, setVer); // its sa
  else              its.MakeAliceCurrent(tpc,0); // its sa  
  //its.MakeAliceCurrent(0,1); // with tpc
  its.SetMaxSeedToPropagate(300);
  its.SetUseBackground();
  //
  // min hits per track to validate it
  its.RequireMinITSHits( its.GetNActiveITSLayers() - 3 );//(int)TMath::Max(4.0,(3./4.)*its.GetNActiveITSLayers()));
  printf("Min number of hits requested: %d\n",its.GetMinITSHits());
  //
  // max cluster-track chi2 
  its.RequireMaxChi2Cl(16.);
  // max chi2/NDF of the Kalman fit
  its.RequireMaxNormChi2NDF(5.);
  // penalty to chi2 from missing hit
  KMCProbe::SetMissingHitPenalty(6);
  //
  // Detector, tracking initialization <<<

  its.InitMCWatchHistos();
  summary.Clear();
  //
  int npt = sizeof(pts)/sizeof(double);
  printf("NPt bins = %d\n",npt);
  double ptMin = pts[0];
  double ptMax = pts[npt-1];
  //
  for (int ip=npt;ip--;) {
    //
    int nevR = nev;
    double pt = pts[ip];//1./(ptminI+dpt*ip);  // ptmin+ip*(ptmax-ptmin)/npt;

    double chi2Cl = 16.;
    if (npt>1) chi2Cl -= (16.-9.)*(1./pt-1./ptMax)/(1./ptMin-1./ptMax);
    //    if (pt<0.3) nevR = int(0.25*nevR);
    //    else if (pt<0.5) nevR = int(0.7*nevR);
    printf("Doing pt(%d)=%.3f for mass %.3f, %d ev | chi2Cl=%.2f\n",ip,pt,mass,nevR,chi2Cl);
    TString smn = Form("pt%d",ip);
    TObjArray* trSum = CreateTrackConditions(smn.Data());
    if (trSum) summary.AddLast(trSum);
    //
    its.SolveSingleTrack(mass, pt, eta, trSum, nevR);
  }
  //
  TString smout = sumOut;
  if (!smout.IsNull()) {
    gSystem->ExpandPathName(smout);
    TFile* fl = TFile::Open(smout.Data(),"recreate");
    if (!fl) {printf("Failed to open %s file\n",smout.Data()); return;}
    fl->WriteObject(&summary,"trSum","kSingleKey");
    fl->Close();
    delete fl;
  }
}

//_____________________________________________________________
TObjArray* CreateTrackConditions(const char* nm)
{
  // create set of track criteria to extract the summary
  KMCTrackSummary *sm = 0;
  TObjArray* arr = new TObjArray();
  //
  int nlr = KMCProbe::GetNITSLayers();
  //-------------------- put here user code -------------------
  // summary for ideal correct tracks
  sm = new KMCTrackSummary("corrAll","corrAll",nlr);
  sm->SetMinMaxClITS(nlr);
  sm->SetMinMaxClITSFake(0,0);
  sm->SetNamePrefix(nm);
  sm->AddPatternITSCorr( sm->Bits(1,0));    // the innermost 2 layers must have correct cluster
  sm->AddPatternITSCorr( sm->Bits(0,1));
  arr->AddLast(sm);
  //
  sm = new KMCTrackSummary("any","any",nlr);
  sm->SetMinMaxClITS(0);
  sm->SetMinMaxClITSFake(0,999);
  sm->SetNamePrefix(nm);
  arr->AddLast(sm);
  //
  // summary for good correct tracks
  sm = new KMCTrackSummary("corr3hSPD","corr3hSPD",nlr);
  sm->SetMinMaxClITS(3);
  sm->SetMinMaxClITSFake(0,0);
  sm->SetNamePrefix(nm);
  arr->AddLast(sm);
  //
  // summary for good correct tracks
  sm = new KMCTrackSummary("corr4h","corr4h",nlr);
  sm->SetMinMaxClITS(4);
  sm->SetMinMaxClITSFake(0,0);
  sm->SetNamePrefix(nm);
  arr->AddLast(sm);
  //
  // summary for good correct tracks with enough points at inner single layers
  sm = new KMCTrackSummary("corr4hs2in","corr4hs2in",nlr);
  sm->SetMinMaxClITS(4);
  sm->SetMinMaxClITSFake(0,0);
  sm->SetNamePrefix(nm);
  sm->AddPatternITSCorr( sm->Bits(1,0));    // the innermost 2 layers must have correct cluster
  sm->AddPatternITSCorr( sm->Bits(0,1));   
  arr->AddLast(sm);
  //
  // summary for good correct tracks with enough points at inner double layers
  sm = new KMCTrackSummary("corr4hd2in","corr5hd2in",nlr);
  sm->SetMinMaxClITS(4);
  sm->SetMinMaxClITSFake(0,0);
  sm->SetNamePrefix(nm);
  sm->AddPatternITSCorr( sm->Bits(1));    // the innermost 2 double layers must have correct cluster
  sm->AddPatternITSCorr( sm->Bits(0,1,1));   
  arr->AddLast(sm);
  //
  //
  // summary for golden tracks with AllNew RS setup: good pattern for offset and pt measurement
  sm = new KMCTrackSummary("corr4hdGoodPatt","corr4hdGoodPatt",nlr);
  sm->SetMinMaxClITS(4);
  sm->SetMinMaxClITSFake(0,0);
  sm->SetNamePrefix(nm);
  sm->AddPatternITSCorr( sm->Bits(1,1));    // the innermost 2 double layers must have correct cluster
  sm->AddPatternITSCorr( sm->Bits(0,0,1)); 
  //
  sm->AddPatternITSCorr( sm->Bits(0,0,0,1,1)); // measurement in the middle (would be prefereable 2 points!)
  //
  sm->AddPatternITSCorr( sm->Bits(0,0,0, 0,0, 1,1)); // measurement at large radius
  //
  arr->AddLast(sm);
  //
  // summary for golden tracks with curr setup: good pattern for offset and pt measurement
  sm = new KMCTrackSummary("corr4hdGoodPatt","corr4hdGoodPatt",nlr);
  sm->SetMinMaxClITS(4);
  sm->SetMinMaxClITSFake(0,0);
  sm->SetNamePrefix(nm);
  sm->AddPatternITSCorr( sm->Bits(1,0));    // the innermost 2 double layers must have correct cluster
  sm->AddPatternITSCorr( sm->Bits(0,1)); 
  //
  sm->AddPatternITSCorr( sm->Bits(0,0,1,1)); // measurement in the middle (would be prefereable 2 points!)
  //
  sm->AddPatternITSCorr( sm->Bits(0,0, 0,0, 1,1)); // measurement at large radius
  //
  arr->AddLast(sm);
  //
  // summary for tracks with at least 1 fake cluster
  sm = new KMCTrackSummary("fakeAny","fakeAny",nlr);
  sm->SetMinMaxClITS(4);
  sm->SetMinMaxClITSFake(1);
  sm->SetNamePrefix(nm);
  arr->AddLast(sm);
  //  
  // // summary for tracks with at least 1 fake cluster (but still acceptable?)
  // sm = new KMCTrackSummary("corrF","corrF",nlr);
  // sm->SetMinMaxClITS(4);
  // sm->SetMinMaxClITSCorr(3);              // at least 3 correct clusters
  // sm->SetMinMaxClITSFake(1);              // at least 1 fake clusters
  // //  /*
  // sm->AddPatternITSCorr( sm->Bits(1));    // the innermost 2 layers must have correct cluster
  // sm->AddPatternITSCorr( sm->Bits(0,1,1));   
  // //  */
  // //  sm->AddPatternITSCorr( sm->Bits(1,1));    // 
  // // sm->AddPatternITSCorr( sm->Bits(0,0,1,1,1,1));   

  // //  sm->AddPatternITSFakeExcl( sm->Bits(1,1,1)); // no fakes in innermost 3 layers
  // sm->SetNamePrefix(nm);
  // arr->AddLast(sm);
  //  
  //-----------------------------------------------------------
  return arr->GetEntries()>0 ? arr : 0;
}
