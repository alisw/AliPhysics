#include "AliMultVariable.h"
#include "AliMultEstimator.h"
#include "AliMultSelection.h"
#include "AliMultSelectionCalibrator.h"
#include "AliMultSelectionCalibratorMC.h"

void CalibratePeriod_LHC18r(TString lPeriodName = "LHC18r"){
  
  //Load ALICE stuff
  TString gLibs[] =   {"STEER",
    "ANALYSIS", "ANALYSISalice", "ANALYSIScalib"};
  TString thislib = "lib";
  for(Int_t ilib = 0; ilib<4; ilib++){
    thislib="lib";
    thislib.Append(gLibs[ilib].Data());
    cout<<"Will load "<<thislib.Data()<<endl;
    gSystem->Load(thislib.Data());
  }
  gSystem->SetIncludePath("-I$ROOTSYS/include  -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  
  cout<<"Alive! "<<endl;
  
  //All fine, let's try the calibrator
  AliMultSelectionCalibrator *lCalib = new AliMultSelectionCalibrator("lCalib");
  AliMultSelection *lMultSelDefault = new AliMultSelection();
  //lCalib->SetMultSelection(lMultSelDefault);
  
  //============================================================
  // --- Definition of Boundaries ---
  //============================================================
  
  //Set Adaptive Percentile Boundaries, adjust if finer selection desired
  Double_t lDesiredBoundaries[1001];
  Long_t   lNDesiredBoundaries=0;
  lDesiredBoundaries[0] = 100;
  /*
   //From Low To High Multiplicity
   for( Int_t ib = 1; ib < 91; ib++) {
   lNDesiredBoundaries++;
   lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] - 1.0;
   }
   for( Int_t ib = 1; ib < 91; ib++) {
   lNDesiredBoundaries++;
   lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] - 0.1;
   }
   for( Int_t ib = 1; ib < 91; ib++) {
   lNDesiredBoundaries++;
   lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] - 0.01;
   }
   for( Int_t ib = 1; ib < 101; ib++) {
   lNDesiredBoundaries++;
   lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] - 0.001;
   }
   lNDesiredBoundaries++;
   */
  
  //Very simple 1%-wide bins all the way
  for( Int_t ib = 1; ib < 1000; ib++) {
    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] - 0.1;
  }
  lNDesiredBoundaries++;
  lDesiredBoundaries[lNDesiredBoundaries] = 0;
  //cout<<"set i = "<<lNDesiredBoundaries<<" at zero"<<endl;
  
  //cout<< "Dump for debug: "<<endl;
  //for( Int_t ib=0;ib<1001; ib++)
  //    cout<<"Boundary #"<<ib<<" at "<<lDesiredBoundaries[ib]<<endl;
  
  lCalib->SetBoundaries( lNDesiredBoundaries, lDesiredBoundaries );
  cout<<"Boundaries set. Will attempt calibration now... "<<endl;
  //return;
  //// ANCHOR POINT /////
  Double_t lDefaultV0MAnchor     = 0;
  Double_t lDefaultV0MPercentile = 0;
  
  Double_t lDefaultCL0Anchor     = 0;
  Double_t lDefaultCL0Percentile = 0;
  
  Double_t lDefaultCL1Anchor     = 0;
  Double_t lDefaultCL1Percentile = 0;
  
  if ( lPeriodName.Contains("LHC15o") ){
    lDefaultV0MAnchor     = 133.5;
    lDefaultV0MPercentile = 90.007;
    lDefaultCL0Anchor     = 35.5;
    lDefaultCL0Percentile = 90.04;
    lDefaultCL1Anchor     = 30.5;
    lDefaultCL1Percentile = 90.09;
  }
  if ( lPeriodName.Contains("LHC15m") ){
    lDefaultV0MAnchor     = 115.0;
    lDefaultV0MPercentile = 87.5;
    lDefaultCL0Anchor     = 39.5;
    lDefaultCL0Percentile = 88.9;
    lDefaultCL1Anchor     = 40.5;
    lDefaultCL1Percentile = 88.1;
  }
  if ( lPeriodName.Contains("LHC10h") ){
    lDefaultV0MAnchor     = 81.0;
    lDefaultV0MPercentile = 90.00;
    lDefaultCL0Anchor     = 23.5;
    lDefaultCL0Percentile = 90.20;
    lDefaultCL1Anchor     = 21.5;
    lDefaultCL1Percentile = 90.19;
  }
  
  if ( lPeriodName.Contains("LHC18r") ){
    lDefaultV0MAnchor     = 165;
    lDefaultV0MPercentile = 90.1;
    lDefaultCL0Anchor     = 32;
    lDefaultCL0Percentile = 90.11;
    lDefaultCL1Anchor     = 32;
    lDefaultCL1Percentile = 90.07;
  }
  
  ///// EVENT SELECTION /////
  
  lCalib->GetEventCuts()->SetVzCut(10.0);
  lCalib->GetEventCuts()->SetTriggerCut                (kTRUE );
  lCalib->GetEventCuts()->SetINELgtZEROCut             (kFALSE);
  lCalib->GetEventCuts()->SetTrackletsVsClustersCut    (kFALSE );
  lCalib->GetEventCuts()->SetRejectPileupInMultBinsCut (kFALSE);
  lCalib->GetEventCuts()->SetVertexConsistencyCut      (kFALSE);
  lCalib->GetEventCuts()->SetNonZeroNContribs          (kTRUE );
  
  //============================================================
  // --- Definition of Input Variables ---
  //============================================================
  
  lCalib->SetupStandardInput();
  
  //============================================================
  // --- Definition of Estimators ---
  //============================================================
  
    AliMultEstimator *fEstV0M = new AliMultEstimator("V0M", "", "((fAmplitude_V0A)+(fAmplitude_V0C))/((((9.19807e+03)+(1.13465e+01)*(fEvSel_VtxZ)+(-2.66087e-01)*TMath::Power((fEvSel_VtxZ),2)+(2.15205e-02)*TMath::Power((fEvSel_VtxZ),3))/(9.19807e+03)))");
  
  fEstV0M -> SetUseAnchor        ( kTRUE   ) ;
  fEstV0M -> SetAnchorPoint      ( lDefaultV0MAnchor     ) ;
  fEstV0M -> SetAnchorPercentile ( lDefaultV0MPercentile ) ;
  
  
  //Only do this in run 2, AD didn't exist in Run 1
  //Will also save space in the OADB for old datasets!
  //if( lPeriodName.Contains("LHC15") ){
  //lCalib->GetMultSelection() -> AddEstimator( fEstOnlineV0M );
  //lCalib->GetMultSelection() -> AddEstimator( fEstOnlineV0A );
  //lCalib->GetMultSelection() -> AddEstimator( fEstOnlineV0C );
  //lCalib->GetMultSelection() -> AddEstimator( fEstADM );
  //lCalib->GetMultSelection() -> AddEstimator( fEstADA );
  //lCalib->GetMultSelection() -> AddEstimator( fEstADC );
  //}
  
  //Universal: Tracking, etc
  //lCalib->GetMultSelection() -> AddEstimator( fEstnSPDClusters  );
  //lCalib->GetMultSelection() -> AddEstimator( fEstnSPDTracklets );
  //lCalib->GetMultSelection() -> AddEstimator( fEstRefMultEta5 );
  //lCalib->GetMultSelection() -> AddEstimator( fEstRefMultEta8 );
  
  
  //Univeral: V0
  lMultSelDefault -> AddEstimator( fEstV0M );
  
  //Needed as basis of MC calibrator -- keep this here, please!
  AliMultEstimator *fEstnSPDTracklets = new AliMultEstimator("SPDTracklets", "", "(fnTracklets)");
  fEstnSPDTracklets->SetIsInteger(kTRUE);
  lMultSelDefault -> AddEstimator( fEstnSPDTracklets );
  
  //Config central CL0
  AliMultEstimator *fEstCL0 = new AliMultEstimator("CL0", "", "(fnSPDClusters0)/(1 + ((fEvSel_VtxZ)-1.2)*(0.0092 + ((fEvSel_VtxZ)-1.2)*(-0.0018)))");
  fEstCL0 -> SetUseAnchor        ( kTRUE  ) ;
  fEstCL0 -> SetAnchorPoint      ( lDefaultCL0Anchor   ) ;
  fEstCL0 -> SetAnchorPercentile ( lDefaultCL0Percentile  ) ;
  lMultSelDefault -> AddEstimator( fEstCL0 );
  
  //Config central CL1
  AliMultEstimator *fEstCL1 = new AliMultEstimator("CL1", "", "(fnSPDClusters1)/(1 + ((fEvSel_VtxZ)-1.2)*(0.0057 - ((fEvSel_VtxZ)-1.2)*(-0.0013)))");
  fEstCL1 -> SetUseAnchor        ( kTRUE   ) ;
  fEstCL1 -> SetAnchorPoint      ( lDefaultCL1Anchor    ) ;
  fEstCL1 -> SetAnchorPercentile ( lDefaultCL1Percentile  ) ;
  lMultSelDefault -> AddEstimator( fEstCL1 );
  
  
  AliMultSelection *lMultSelArray[500]; //We need one AliMultSelection for each run above!
  
  TF1 *f = new TF1("f", "[0]+[1]*x", 0,50);
  f->SetParameter(0, 5.46681e+01/5.46681e+01);
  f->SetParameter(1, -2.94245e-02/5.46681e+01);
  
  Double_t lV0MAnchors[100];
  for(Int_t i=0; i<100; i++) lV0MAnchors[i] = 165*f->Eval(i);
  
  //Sweep array
  for(Int_t i=0; i<100; i++){
    lMultSelArray[i] = new AliMultSelection(lMultSelDefault);
    lMultSelArray[i]->GetEstimator("V0M")->SetAnchorPoint( lV0MAnchors[i] );
  }
  
  Int_t lCounter = 0;
  /*
   lCalib->AddRunRange( 295581, 295589, lMultSelArray[lCounter++] );
   lCalib->AddRunRange( 295610, 295615, lMultSelArray[lCounter++] );
   lCalib->AddRunRange( 295665, 295677, lMultSelArray[lCounter++] );
   lCalib->AddRunRange( 295711, 295726, lMultSelArray[lCounter++] );
   lCalib->AddRunRange( 295753, 295763, lMultSelArray[lCounter++] );
   
   lCalib->AddRunRange( 295786, 295791, lMultSelArray[lCounter++] );
   lCalib->AddRunRange( 295816, 295831, lMultSelArray[lCounter++] );
   lCalib->AddRunRange( 295853, 295881, lMultSelArray[lCounter++] );
   lCalib->AddRunRange( 295906, 295916, lMultSelArray[lCounter++] );
   lCalib->AddRunRange( 295936, 295947, lMultSelArray[lCounter++] );
   
   lCalib->AddRunRange( 296006, 296016, lMultSelArray[lCounter++] );
   lCalib->AddRunRange( 296060, 296074, lMultSelArray[lCounter++] );
   lCalib->AddRunRange( 296123, 296143, lMultSelArray[lCounter++] );
   lCalib->AddRunRange( 296191, 296198, lMultSelArray[lCounter++] );
   lCalib->AddRunRange( 296240, 296247, lMultSelArray[lCounter++] );
   
   lCalib->AddRunRange( 296270, 296280, lMultSelArray[lCounter++] );
   lCalib->AddRunRange( 296303, 296312, lMultSelArray[lCounter++] );
   lCalib->AddRunRange( 296375, 296383, lMultSelArray[lCounter++] );
   lCalib->AddRunRange( 296414, 296472, lMultSelArray[lCounter++] );
   lCalib->AddRunRange( 296508, 296517, lMultSelArray[lCounter++] );
   
   lCalib->AddRunRange( 296547, 296555, lMultSelArray[lCounter++] );
   lCalib->AddRunRange( 296594, 296594, lMultSelArray[lCounter++] );
   lCalib->AddRunRange( 296614, 296623, lMultSelArray[lCounter++] );
   */
  lCounter = 23; //start in the middle of the period (important for anchor point aging)
  //cout<<"lCounter = "<<lCounter<<endl;
  lCalib->AddRunRange( 296690, 296698, lMultSelArray[lCounter++] );
  lCalib->AddRunRange( 296749, 296752, lMultSelArray[lCounter++] );
  
  lCalib->AddRunRange( 296780, 296799, lMultSelArray[lCounter++] );
  lCalib->AddRunRange( 296835, 296852, lMultSelArray[lCounter++] );
  lCalib->AddRunRange( 296890, 296903, lMultSelArray[lCounter++] );
  lCalib->AddRunRange( 296930, 296941, lMultSelArray[lCounter++] );
  lCalib->AddRunRange( 296966, 296979, lMultSelArray[lCounter++] );
  lCalib->AddRunRange( 297028, 297036, lMultSelArray[lCounter++] ); //7474
  lCalib->AddRunRange( 297084, 297085, lMultSelArray[lCounter++] ); //7475
  lCalib->AddRunRange( 297116, 297133, lMultSelArray[lCounter++] ); //7477
  lCalib->AddRunRange( 297193, 297196, lMultSelArray[lCounter++] ); //7480
  lCalib->AddRunRange( 297218, 297222, lMultSelArray[lCounter++] ); //7481
  lCalib->AddRunRange( 297277, 297278, lMultSelArray[lCounter++] ); //7482
  
  lCalib->AddRunRange( 297310, 297336, lMultSelArray[lCounter++] ); //7483
  lCalib->AddRunRange( 297363, 297380, lMultSelArray[lCounter++] ); //7485
  
  lCalib->AddRunRange( 297403, 297415, lMultSelArray[lCounter++] ); //7486
  lCalib->AddRunRange( 297441, 297452, lMultSelArray[lCounter++] ); //7487
  lCalib->AddRunRange( 297479, 297485, lMultSelArray[lCounter++] ); //7488
  lCalib->AddRunRange( 297512, 297513, lMultSelArray[lCounter++] ); //7489
  
  lCalib->AddRunRange( 297537, 297558, lMultSelArray[lCounter++] ); //7490
  lCalib->AddRunRange( 297588, 297595, lMultSelArray[lCounter++] ); //7491
  lCalib->AddRunRange( 297623, 297624, lMultSelArray[lCounter++] ); //7492
  
  
  //return;
  cout<<"Inspect default object:"<<endl;
  lMultSelDefault->PrintInfo();
  
  lCalib->SetRunToUseAsDefault(297031);
  
  //============================================================
  // --- Definition of Input/Output ---
  //============================================================
  //lCalib -> SetInputFile  ( "~/Dropbox/MultSelCalib/LHC15o/MergedLHC15o.root");
  //lCalib -> SetInputFile  ( "/hera/alice/alberica/centrality/Trees/Singles/LHC15o/muon_calo_pass1_1301/files/AnalysisResults_244918.root");
  lCalib -> SetInputFile  ( "LHC18r.root");
  //lCalib -> SetInputFile  ( "~/Dropbox/MultSelCalib/LHC15o/files/AnalysisResults_245064.root");
  
  lCalib -> SetBufferFile ( "buffer-LHC18r.root" );
  lCalib -> SetOutputFile ( "OADB-LHC18r.root" );
  lCalib -> Calibrate     ();
  
}
