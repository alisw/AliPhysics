#include "TSpline.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TFile.h"
#include "TError.h"
#include "TROOT.h"
#include "TPRegexp.h"
#include "TClass.h"
#include "TKey.h"

#include "AliNDLocalRegression.h"
#include "AliOADBContainer.h"
#include "AliTPCPIDResponse.h"
#include "AliPID.h"

AliOADBContainer fContainer("TPCSplines");
Int_t fFailures=0;

Bool_t AddOADBObjectFromSplineFile(const TString fileName,
                                   const Int_t firstRun, const Int_t lastRun,
                                   const TString pass,
                                   const TString dEdxType="",
                                   const TString multCorr="",
                                   const TString resolution="",
                                   const TString pileupDefinition="",
                                   const AliTPCPIDResponse::EMultiplicityEstimator multEstimator = AliTPCPIDResponse::kNumberOfESDTracks);

Bool_t CheckMultiplicityCorrection(const TString& corrections);
TObjArray* SetupSplineArrayFromFile(const TString fileName);

AliOADBContainer* GetOADBContainer(const TString fileName);
TObject* GetObjectFromContainer(AliOADBContainer* c, const TString objName);

//______________________________________________________________________________
void MakeTPCPIDResponseOADB(TString outfile="$ALICE_PHYSICS/OADB/COMMON/PID/data/TPCPIDResponseOADB.root")
{

  //
  // ===| 2010 |================================================================
  //

  // ---| pass2 |---------------------------------------------------------------
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/10b.pass2/splines_10b.pass2.root", 114650, 117630, "2");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/10c.pass2/splines_10c.pass2.root", 117631, 121526, "2");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/10d.pass2/splines_10d.pass2.root", 121527, 126460, "2");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/10e.pass2/splines_10e.pass2.root", 126461, 130930, "2");

  // -- PbPb
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/LHC10h.pass2/TPCpidResponseFunctions.root",136782,139846, "2");

  // ---| pass3 |---------------------------------------------------------------
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/LHC10b.pass3/splines_10b.pass3.root",114650,117630, "3");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/LHC10c.pass3/splines_10c.pass3.root",117631,121526, "3");

  // ---| pass4 |---------------------------------------------------------------
  // first version of splines for the 2010 pass4 rereconstruction. No correction for the pt problem, yet
  // using the standard dEdx
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/10b.pass4/first/splines_10b.pass4.root", 114650, 117630, "4");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/10c.pass4/first/splines_10c.pass4.root", 117631, 121526, "4");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/10d.pass4/first/splines_10d.pass4.root", 121527, 126460, "4");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/10e.pass4/first/splines_10e.pass4.root", 126461, 130930, "4");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/10f.pass4/first/splines_10f.pass4.root", 130962, 135393, "4");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/10g.pass4/first/splines_10g.pass4.root", 135394, 136781, "4");


  //
  // ===| 2011 |================================================================
  //
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/11a_7TeV.pass1/splines_11a_7TeV.pass1.root"         , 139847, 146974, "1");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/11a_wo_SDD.pass2/splines_11a_without_SDD.pass2.root", 139847, 146974, "2");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/11a_wo_SDD.pass4/splines_11a_without_SDD.pass4.root", 139847, 146974, "4");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/LHC11a.pass3/splines_11aPass3.root"                          , 139847, 146974, "3");

  // use 11d for 11b-11f
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/LHC11d.pass1/splines_11d.root"                               , 146975, 165771, "2");

  // PbPb
  // AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/LHC11h.pass2/splines_11h.2.root"                             , 165772, 170718, "2");
  // splines as above, adding a triton parametrisation from Annalisa Mastroserio and Caio Lagana
  // this is done using the macro /u/wiechula/svn/train/PID/splines/annalisa_triton/myspline.C
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/annalisa_triton/splines_11h.2.root"                             , 165772, 170718, "2");

  // not needed any longer?
  //   AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/LHC11h.pass2.scenarios/splines_11h.allhigh.root",""));
  //   AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/LHC11h.pass2.scenarios/splines_11h.oroc.root",""));


  //
  // ===| 2012 |================================================================
  //

  // ---| pass1 |---------------------------------------------------------------
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/12a.pass1/splines_12a.pass1.root", 170719, 177311, "1");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/12b.pass1/splines_12b.pass1.root", 177312, 179356, "1");
  // pPb pilot
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/LHC12g.pass1/splines_12g.1.root", 188167, 188418, "1");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/LHC12g.pass2/splines_12g.2.root", 188167, 188418, "2");

  // ---| pass2 |---------------------------------------------------------------
  // second iteration of pass2 LHC12X splines
  // originally from /hera/alice/nschmidt/splines/Iteration2/*
  //
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/12a.pass2/splines_12a.pass2.root", 170719, 177311, "2");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/12b.pass2/splines_12b.pass2.root", 177312, 179356, "2");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/12c.pass2/splines_12c.pass2.root", 179357, 183173, "2");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/12d.pass2/splines_12d.pass2.root", 183174, 186345, "2");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/12e.pass2/splines_12e.pass2.root", 186346, 186635, "2");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/12f.pass2/splines_12f.pass2.root", 186636, 188166, "2");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/12g.pass2/splines_12g.pass2.root", 188419, 188719, "2"); // except pPb pilot runs
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/12h.pass2/splines_12h.pass2.root", 188720, 192738, "2");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/12i.pass2/splines_12i.pass2.root", 192739, 194479, "2"); //12i+j


  //
  // ===| 2013 |================================================================
  //

  // pPb - w/ multiplicity correction
  // used for 13a-d
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/13b.pass1.mult/splines_13b.pass1.root", 194480, 195874, "1", "",
                              "-5.906e-06,-5.064e-04,-3.521e-02,2.469e-02,0 ; -5.32e-06, 1.177e-05, -0.5 ; 0.,0.,0.,0.");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/13b.pass2.mult/splines_13b.pass2.root", 194480, 195874, "2", "",
                              "-5.906e-06,-5.064e-04,-3.521e-02,2.469e-02,0 ; -5.32e-06, 1.177e-05, -0.5 ; 0.,0.,0.,0.");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/13b.pass3.mult/splines_13b.pass3.root", 194480, 195874, "3", "",
                              "-5.906e-06,-5.064e-04,-3.521e-02,2.469e-02,0 ; -5.32e-06, 1.177e-05, -0.5 ; 0.,0.,0.,0.");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/nicolas/data/LHC13b.pass4/splines_13c.pass4.root", 194480, 195874, "4","",
                              "-5.906e-06,-5.064e-04,-3.521e-02,2.469e-02,0 ; -5.32e-06, 1.177e-05, -0.5 ; 0.,0.,0.,0.");

  // pPb high lumi used for 13e and 13f
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/13f.pass1/splines_13f.pass1.root", 195875, 197411, "1");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/13f.pass2/splines_13f.pass2.root", 195875, 197411, "2");

  //
  // ===| 2015 |================================================================
  //

  // first iteration of 15f
  // used for all 2015 so far
  // ---| pass1 |---------------------------------------------------------------
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/ben/data/15f.pass1/splines_15f.pass1.root",     208505, 228930, "1"); //a-f (+15g low field)
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/LHC15g.pass1_megarun/splines_LHC15g.pass1_megarun.root",     228931, 229245, "1"); //15g low field
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/nicolas/data/15h.pass1/v1/splines_15h.pass1.root", 229246, 235169, "1"); //includes 15g high field
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/nicolas/data/15i.pass1/v1/splines_15i.pass1.root", 235170, 236866, "1");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/nicolas/data/15j.pass1/v1/splines_15j.pass1.root", 236867, 239154, "1"); //j-k
//AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC15j.pass1/Iteration3/splines.root", 236867, 238653, "pass1"); // not needed, done by mistake. The ones by nicolas will stay in use
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/martin/data/LHC15l.pass1/splines_15l.pass1.root", 239155, 244299, "1"); //l-m
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/nicolas/data/15n.pass1/v1/splines_15n.pass1.root", 244300, 244639, "1");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/nicolas/data/15o.pass1_highIR/v1/splines_15o.pass1.root", 244640, 247173, "1"); // 15o high IR
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/nicolas/data/15o.pass1_highIR_pidfix/v1/splines_15o.pass1.root", 245145, 245554, "pass1_pidfix"); // 15o high IR 245145 - 245554

  // ---| pass2 |---------------------------------------------------------------
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/nicolas/data/15f.pass2/splines_15f.pass2.root", 208505, 229245, "2"); //a-f
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/martin/data/LHC15h.pass2/splines_15h.pass2.root", 229246, 235169, "2"); //includes 15g high field
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC15i.pass2/splines_15i_pass2.root",      235170, 236866, "2");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/nicolas/data/15n.pass2/v1/splines_15n.pass2.root", 244300, 244639, "2");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/nicolas/data/15o.pass2/v1/splines_15o.pass2.root", 244640, 247173, "2", "",
                              "-1.077e-06,-4.999e-05,-9.812e-03,2.492e-02,0 ; -5.562e-07, 1.644e-06, -0.5 ; 0.,0.,0.,0.");

  // ---| pass3 |---------------------------------------------------------------
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/nicolas/data/15n.pass3/v1/splines_15n.pass3.root", 244300, 244639, "3");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/martin/data/LHC15o.pass3_lowIR_pidfix/splines_15o.pass3.root", 244917, 246392, "pass3_lowIR_pidfix", "",
                              "-1.548980e-06,-3.174573e-05,-2.764576e-03,4.965499e-02,0 ; 5.362701e-07,-6.012158e-07, -0.5 ; 0.,0.,0.,0.");

  // ---| pass4 |---------------------------------------------------------------
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/martin/data/LHC15n.pass4/splines_15n.pass4.root", 244300, 244639, "4");


  //
  // ===| 2016 |================================================================
  //

  // for the moment use
  //   - 16d splines for 16a to 16d
  //   - 16e splines for 16e
  //   - 16f splines for the low field runs 253610 - 253883
  //   - 16g,h,i,j splines
  //   - 16k splines for 16k
  //   - 16l splines for 16l to 16n
  //   - 16o,p splines
  // ---| pass1 |---------------------------------------------------------------
  // low field data
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/martin/data/LHC16d.pass1/splines_16d.pass1.root",  247174, 252536, "1");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/martin/data/LHC16e.pass1/splines_16e.pass1.root",  252537, 253609, "1");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/martin/data/LHC16f.pass1/splines_16f.pass1.root",  253610, 253883, "1"); 
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/martin/data/LHC16g.pass1/splines_16g.pass1.root",  253884, 254338, "1"); 
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/martin/data/LHC16h.pass1/splines_16h.pass1.root",  254339, 255476, "1"); 
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/martin/data/LHC16i.pass1/splines_16i.pass1.root",  255477, 255798, "1"); 
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/martin/data/LHC16j.pass1/splines_16j.pass1.root",  255799, 256489, "1"); 
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/nicolas/data/LHC16k.pass1/v1/splines_16k.pass1.root", 256490, 258860, "1");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/nicolas/data/LHC16l.pass1/v1/splines_16l.pass1.root", 258861, 261928, "1");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/martin/data/LHC16o.pass1/splines_16o.pass1.root",  261929, 264611, "1"); 
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/martin/data/LHC16p.pass1/splines_16p.pass1.root",  264612, 264895, "1"); 

  // ---| pPb periods |---------------------------------------------------------
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/martin/data/LHC16q.pass1/splines_16q.pass1.root", 264896, 265533, "1", "",
                             "-1.609459e-06,-6.765851e-04,9.610860e-03,2.864834e-02,0 ; -2.121118e-07,-1.181542e-06, -0.5 ; 0.,0.,0.,0.");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/martin/data/LHC16r.pass1/splines_16r.pass1.root", 265534, 266329, "1", "",
                              "1.033846e-05,-1.292888e-03,3.029898e-02,1.936201e-02,0 ; -1.011865e-06,5.973020e-07, -0.5 ; 0.,0.,0.,0.");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/martin/data/LHC16s.pass1/splines_16s.pass1.root", 266330, 267138, "1", "",
                              "4.450472e-06,-1.932623e-03,6.953790e-02,1.940033e-02,0 ; -2.234788e-06,9.764794e-07, -0.5 ; 0.,0.,0.,0.");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/martin/data/LHC16q.pass1/splines_16q.pass1.root", 267139, 267388, "1", "",
                              "-1.609459e-06,-6.765851e-04,9.610860e-03,2.864834e-02,0 ; -2.121118e-07,-1.181542e-06, -0.5 ; 0.,0.,0.,0."); // same configuration as 16q, so use these splines

  // ---| pass2 |---------------------------------------------------------------
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC16k.pass2/Iteration3/splines.root", 256490, 258860, "2");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC16l.pass2/Iteration3/splines.root", 258861, 260199, "2");


  //
  // ===| 2017 |================================================================
  //
  //
  // For the moment use 17c for 17a-e
  //                    17f for 17f
  //                    17g for 17g (low field)
  //                    17h for 17h
  //                    17j for 17j
  //                    17k for 17k
  //                    17m for 17m
  //                    17n for 17n 
  //                    17o for 17o
  //                    17p for 17p
  //                    17q for 17q
  //                    17r for 17r
  //                    
  //                    Systems: 17c-m pp   13   TeV
  //                             17n   XeXe  5.4 TeV
  //                             17o   pp   13   TeV
  //                             17p   pp    5   TeV low IR (50kHz)
  //                             17q   pp    5   TeV high IR (200-600kHz)
  //                             17r   pp   13   TeV
  // ---| pass1 |---------------------------------------------------------------

  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC17c.pass1/Iteration2/splines.root", 268198, 270853, "1"); //a-e
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC17f.pass1/Iteration2/splines.root", 270854, 270870, "1");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC17g.pass1/splines_17G_pass1.root",  270871, 271785, "1");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC17h.pass1/splines_17H_PASS1.root", 271786, 273471, "1");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC17i.pass1/splines_17I_PASS1.root", 273472, 274518, "1");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC17j.pass1/splines_17j_pass1.root", 274519, 274677, "1");
  //AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC17k.pass1/splines_17K.pass1.root", 274678, 276533, "1");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC17k.pass1_2nd/Iteration3/splines.root", 274678, 276533, "1");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC17l.pass1/splines_17L_PASS1.root", 276534, 278765, "1");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC17m.pass1/splines_17M_PASS1.root", 278766, 280154, "1");

  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/martin/data/LHC17n.pass1/splines_17n.pass1.root",280155, 280246, "1", "",
                              "9.135699e-07,-8.855631e-05,-3.877750e-03,1.729602e-02,0 ; 1.029664e-06,-4.357166e-07, -0.5 ; 0.,0.,0.,0."); //Xe-Xe 5.4TeV

  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC17o.pass1/splines_17O_PASS1.root", 280247, 281969, "1"); // temporary for 17o onwards

  //AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC17p.pass1_CENT_wSDD/splines_17P_PASS1.root", 281970, 282350, "1");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC17p.pass1_CENT_wSDD_2nd/Iteration2/splines.root", 281970, 282350, "1");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC17q.pass1_CENT_wSDD/splines_17Q_PASS1.root", 282351, 282444, "1");

  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC17r.pass1/splines_17R_PASS1.root", 282445, 282955, "1"); // temporary for 17r onwards

  //
  // ===| 2018 |================================================================
  //
  // For the moment use 16p splines (pp 13 TeV with Neon) for 18a-b and 18e-
  //                    18b for 18a-b
  //                    18c for 18c
  //                    18d for 18d-h
  //                    18i for 18i -
  //
  // ---| pass1 |---------------------------------------------------------------
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC18b.pass1/Iteration2/splines.root",             282956, 285451, "1"); // 18a-b
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC18c.pass1_CENT_woSDD/Iteration2/splines.root",  285452, 285965, "1"); 
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC18d.pass1/Iteration3/splines.root",             285966, 286358, "1");
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC18e.pass1/Iteration3/splines.root",             286359, 288842, "1");// last run in 18e: 286969, "1"); // LHC18e-h

  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC18i.pass1/Iteration3/splines.root",             288843, 290110, "1"); // 18i-l last run in 18i: 288920, "1");

  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC18m.pass1/Iteration3/splines.root",             290111, 295242, "1"); // 18m-p last run in 18m: 293253, "1"); // LHC18m
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/jiyoung/LHC18r.pass1/Iteration4/splines.root",             295243, 999999, "1", "",
      "-1.945069e-07,-5.163672e-04;-3.168292e-11,2.773070e-08;3.529986e-06,4.000030e-04,8.761510e-02,1.453940e-02"                             ); // 18q-r 18r range: 296631, 999999, "1"); // LHC18r



  // ---| pass2 |---------------------------------------------------------------
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/mciupek/LHC18q.pass3/TPCPIDResponseOADB_2020_08_13_18q_pass3_It3_withDeuteron.root", 295243, 296630, "3", "",
      "", "", "", AliTPCPIDResponse::kNTPCTrackBeforeClean                             ); // 18q
  AddOADBObjectFromSplineFile("/u/wiechula/svn/train/PID/splines/mciupek/LHC18r.pass3/TPCPIDResponseOADB_2020_07_22_18r_pass3_It8_withDeuteron.root", 296631, 999999, "3", "",
      "", "", "", AliTPCPIDResponse::kNTPCTrackBeforeClean                             ); // 18r

/*
  // ---| local test |----------------------------------------------------------
  AddOADBObjectFromSplineFile("/data/Work/data/PID/testTransferFunction/10e/splines_10d.pass4.root", 194480, 195874, "1", "",
                              "-5.906e-06,-5.064e-04,-3.521e-02,2.469e-02,0 ; -5.32e-06, 1.177e-05, -0.5 ; 0.,0.,0.,0.");
*/

  if (fFailures) {
    Error("MakeSplineOADB","Process ended with %d fFailures. Not writing the container", fFailures);
    return;
  }

  fContainer.WriteToFile(outfile.Data());
}

//______________________________________________________________________________
TObjArray* SetupSplineArrayFromFile(const TString fileName)
{
  TFile f(fileName);
  if (!f.IsOpen() || f.IsZombie()) {
    Error("AddOADBObject","Could not open file '%s'",fileName.Data());
    return 0x0;
  }
  gROOT->cd();

  // ===| Spline Array |========================================================
  // A sorted array with the spline objects following the particle nameing convention
  //   in AliPID::EParticleType up to AliPID::kSPECIESC
  //   usually only e, pi, K, p will be available, for the muons the pion spline will be used
  //   for the light nuclei the proton spline will be used. This needs to be set up
  //   in the loading of the splines on order to save space in the file
  // At the position AliPID::kSPECIESC the spline without low momentum corrections is placed

  TObjArray *arr=new TObjArray(AliPID::kSPECIESC+1);
  arr->SetName("Splines");
  arr->SetOwner();

  TIter nextKey(f.GetListOfKeys());

  // ---| loop over all splines in the file |-----------------------------------
  TObject *o=0x0;
  while ( (o=nextKey()) ) {
    TSpline3 *sp=dynamic_cast<TSpline3*>(f.Get(o->GetName()));

    if (!sp) continue;
    const TString splineName(sp->GetName());

    // ---| check for default spline |------------------------------------------
    if (splineName.Contains("ALL")) {
      arr->AddAt(sp, AliPID::kSPECIESC);
      continue;
    }

    // ---| check for particles |------------------------------------------------
    TString partName;
    for (Int_t ispec=0; ispec<Int_t(AliPID::kSPECIESC); ++ispec) {
      partName=AliPID::ParticleName(ispec);
      partName.ToUpper();

      if (splineName.Contains(partName)) {
        arr->AddAt(sp, ispec);
        break;
      }
    }
  }

  return arr;
}

//______________________________________________________________________________
Bool_t AddOADBObjectFromSplineFile(const TString fileName,
                                   const Int_t firstRun, const Int_t lastRun,
                                   const TString pass,
                                   const TString dEdxType,
                                   const TString multCorr,
                                   const TString resolution,
                                   const TString pileupDefinition,
                                   const AliTPCPIDResponse::EMultiplicityEstimator multEstimator)
{

  // ---| Check if input file is an OADB file already, get the OADB container |-
  AliOADBContainer* contFromFile = GetOADBContainer(fileName);

  // ---| Master array for TPC PID response |-----------------------------------
  TObjArray *arrTPCPIDResponse = new TObjArray;
  arrTPCPIDResponse->SetName("TPCPIDResponse");
  arrTPCPIDResponse->SetOwner();

  // ---| Period and pass name deduced from the file name |---------------------
  TPRegexp regPeriod(".*/((LHC|)[0-9]{2}[a-z].pass[0-9_a-zA-Z]*)/.*");
  TObjArray* arrPeriod = regPeriod.MatchS(fileName);
  if (arrPeriod && arrPeriod->GetEntriesFast() == 3) {
    TString periodName = arrPeriod->At(1)->GetName();
    if (! periodName.BeginsWith("LHC") ) periodName.Prepend("LHC");
    printf("period name: %s\n", periodName.Data());
    arrTPCPIDResponse->SetName(periodName);
  }
  else {
    printf("Cannot extract period and pass from file name %s\n", fileName.Data());
  }

  delete arrPeriod;

  // ---| Splines |-------------------------------------------------------------
  if (contFromFile) {
    arrTPCPIDResponse->Add(GetObjectFromContainer(contFromFile, "Splines"));
  }
  else {
    TObjArray *arrSplines = SetupSplineArrayFromFile(fileName);
    if (!arrSplines) {
      ++fFailures;
      return kFALSE;
    }
    arrTPCPIDResponse->Add(arrSplines);
  }

  // ---| PID config |----------------------------------------------------------
  if (!dEdxType.IsNull()) {
    TNamed *pidConfig = new TNamed("dEdxType",dEdxType.Data());
    arrTPCPIDResponse->Add(pidConfig);
  }

  // ---| multiplicity correction |---------------------------------------------
  if (contFromFile && multCorr.IsNull()) {
    arrTPCPIDResponse->Add(GetObjectFromContainer(contFromFile, "MultiplicityCorrection"));
  }
  else {
    TObjArray *arrMultiplicityCorrection=0x0;
    if (!multCorr.IsNull()) {
      if ( (arrMultiplicityCorrection=AliTPCPIDResponse::GetMultiplicityCorrectionArrayFromString(multCorr)) ) {
        TNamed *multCorrConfig = new TNamed("MultiplicityCorrection",multCorr.Data());
        arrTPCPIDResponse->Add(multCorrConfig);
        delete arrMultiplicityCorrection;
      } else {
        ++fFailures;
      }
    }
  }

  // ---| resolution parametrisation |------------------------------------------
  if (!resolution.IsNull() && !resolution.EndsWith(".root")) {
    TNamed *resolutionParam = new TNamed("dEdxResolution", resolution.Data());
    arrTPCPIDResponse->Add(resolutionParam);
  }

  // ---| Pileup correction |---------------------------------------------------
  if (!pileupDefinition.IsNull()) {
    if (pileupDefinition.BeginsWith("fileName:")) {
      // if the string begins with 'fileName:' add the fileName to the array
      TString fileName(pileupDefinition);
      fileName.ReplaceAll("fileName:", "");
      TNamed *pileupCorrectionFile = new TNamed("PileupCorrectionFile", fileName.Data());
      arrTPCPIDResponse->Add(pileupCorrectionFile);
    } else {
      // else directly load the object from file and add it
      AliNDLocalRegression* pileupCorrection = AliTPCPIDResponse::GetPileupCorrectionFromFile(pileupDefinition);
      arrTPCPIDResponse->Add(pileupCorrection);
    }
  }
  fContainer.AppendObject(arrTPCPIDResponse, firstRun, lastRun, pass);

  // ---| multiplicity estimator |---
  if (multEstimator != AliTPCPIDResponse::kNumberOfESDTracks) {
    TNamed *multEstimatorDef = new TNamed("MultiplicityEstimator", Form("%d", (Int_t)multEstimator));
    arrTPCPIDResponse->Add(multEstimatorDef);
  }
  else {
    if (contFromFile) {
      TObject* multEst = GetObjectFromContainer(contFromFile, "MultiplicityEstimator");
      if (multEst) {
        arrTPCPIDResponse->Add(multEst);
      }
    }
  }

  // ---| TF1 sigma parametrization |-------------------------------------------
  if (resolution.EndsWith(".root")) {
    TFile* f = TFile::Open(resolution);
    if (!f->IsOpen() || f->IsZombie()) {
      Error("AddOADBObjectFromSplineFile", "Could not open '%s' to extract the TF1 sigma parametrization", resolution.Data());
    } else {
      TObject* tf1Sigma = f->Get("SigmaParametrization");
      TObject* tf1SigmaParams = f->Get("SigmaParametrizationParams");
      TObject* multEstimator = f->Get("MultiplicityNormalization");

      if (!tf1Sigma) {
        Fatal("AddOADBObjectFromSplineFile", "Could not get TF1 function with name 'SigmaParametrization' from file '%s'", resolution.Data());
      }

      if (!tf1SigmaParams) {
        Fatal("AddOADBObjectFromSplineFile", "Could not get TF1 function with name 'SigmaParametrization' from file '%s'", resolution.Data());
      }

      if (!multEstimator) {
        Fatal("AddOADBObjectFromSplineFile", "Could not get 'MultiplicityNormalization' from file '%s'", resolution.Data());
      }

      TObjArray* arrSigmaParam = new TObjArray;
      arrSigmaParam->SetName("SigmaParametrization");
      arrSigmaParam->Add(tf1Sigma);
      arrSigmaParam->Add(tf1SigmaParams);
      arrSigmaParam->Add(multEstimator);

      arrTPCPIDResponse->Add(arrSigmaParam);
    }
  } else {
    if (contFromFile) {
      TObject* tf1Sigma = GetObjectFromContainer(contFromFile, "SigmaParametrization");
      if (tf1Sigma) {
        arrTPCPIDResponse->Add(tf1Sigma);

        TObject* multEstimator = GetObjectFromContainer(contFromFile, "MultiplicityNormalization");
        if (!multEstimator) {
          Fatal("AddOADBObjectFromSplineFile", "'MultiplicityNormalization' missing in OADB container, always must be provided together with 'SigmaParametrization'");
        } else {
          arrTPCPIDResponse->Add(multEstimator);
        }
      }
    }
  }

  return kTRUE;
}

//______________________________________________________________________________
Bool_t CheckMultiplicityCorrection(const TString& corrections)
{
  // the corrections string is supposed to be in the format
  // parMC_0,parMC_1,parMC_2,parMC_3, parMC_4; parMCTT_0,parMCTT_1,parMCTT_2; parMSC_0,parMSC_1,parMSC_2, parMSC_3
  // where parMC are the parameters for AliTPCPIDResponse::SetParameterMultiplicityCorrection
  // parMCTT are the parameters for AliTPCPIDResponse::SetParameterMultiplicityCorrectionTanTheta
  // parMSC are the parameters for AliTPCPIDResponse::SetParameterMultiplicitySigmaCorrection

  const Int_t nSets=3;
  const Int_t nPars[nSets]={5,3,4};

  TObjArray *arrCorrectionSets = corrections.Tokenize(";");
  if (arrCorrectionSets->GetEntriesFast() != nSets){
    Error("CheckMultiplicityCorrection","Number of parameter sets not equal to %d. Please read documentation of this function", nSets);
    ++fFailures;
    delete arrCorrectionSets;
    return kFALSE;
  }

  for (Int_t iset=0 ;iset<nSets; ++iset) {
    TObjArray *arrTmp=((TObjString*)arrCorrectionSets->At(iset))->String().Tokenize(",");
    if (arrTmp->GetEntriesFast() != nPars[iset]){
      Error("CheckMultiplicityCorrection","Number of parameters in set %d not equal to %d. Please read documentation of this function", iset, nPars[iset]);
      ++fFailures;
      delete arrCorrectionSets;
      delete arrTmp;
      return kFALSE;
    }

  }

  delete arrCorrectionSets;
  return kTRUE;
}

TObject* GetObjectFromContainer(AliOADBContainer* c, const TString objName)
{
  TObjArray* arr = (TObjArray*)c->GetObjectByIndex(0);
  return arr->FindObject(objName);
}

AliOADBContainer* GetOADBContainer(const TString fileName)
{
  TFile f(fileName);
  if (!f.IsOpen() || f.IsZombie()) {
    Error("AddOADBObject","Could not open file '%s'",fileName.Data());
    return 0x0;
  }

  TList* keys = f.GetListOfKeys();
  if (keys->GetEntries()>0) {
    TKey* key = (TKey*)keys->At(0);
    if (TString(key->GetClassName()) == "AliOADBContainer") {
      AliOADBContainer* c = (AliOADBContainer*)f.Get(key->GetName());
      return c;
    }
  }

  return 0x0;
}
