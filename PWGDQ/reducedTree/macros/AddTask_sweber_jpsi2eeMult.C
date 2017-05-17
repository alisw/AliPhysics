//void Setup(AliReducedAnalysisJpsi2eeMult* processor, TString prod="LHC10h");
//AliHistogramManager* SetupHistogramManager(TString prod="LHC10h");
//void DefineHistograms(AliHistogramManager* man, TString prod="LHC10h");

// #define VZEROBINS 
#define SYSTEMATICS 0




//__________________________________________________________________________________________
AliAnalysisTask* AddTask_sweber_jpsi2eeMult(Bool_t isAliRoot=kTRUE, Int_t runMode=1, TString prod="LHC10h"){    
   //
   // isAliRoot=kTRUE for ESD/AOD analysis in AliROOT, kFALSE for root analysis on reduced trees
   // runMode=1 (AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents)
   //               =2 (AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree)
   //
   //get the current analysis manager

  printf("INFO on AddTask_sweber_jpsi2eeMult(): (isAliRoot, runMode) :: (%d,%d) \n", isAliRoot, runMode);

  AliReducedAnalysisJpsi2eeMult* jpsi2eeAnalysis = new AliReducedAnalysisJpsi2eeMult("Jpsi2eeAnalysis","Jpsi->ee analysis");
  jpsi2eeAnalysis->Init();
  Setup(jpsi2eeAnalysis, prod);
  // initialize an AliAnalysisTask which will wrapp the AliReducedAnalysisJpsi2eeMult such that it can be run in an aliroot analysis train (e.g. LEGO, local analysis etc)
  AliAnalysisTaskReducedEventProcessor* task = new AliAnalysisTaskReducedEventProcessor("ReducedEventAnalysisManager", runMode);
  task->AddTask(jpsi2eeAnalysis);
  
  if(isAliRoot){
     AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
     if (!mgr) {
        Error("AddTask_iarsene_dst", "No analysis manager found.");
        return 0;
     }
     
     AliAnalysisDataContainer* cReducedEvent = NULL;
     if(runMode==AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents) {
       printf("INFO on AddTask_sweber_jpsi2eeMult(): use on the fly events ");
       cReducedEvent = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("ReducedEventDQ");
       if(!cReducedEvent) {
         printf("ERROR: In AddTask_sweber_jpsi2eeMult(), couldn't find exchange container with ReducedEvent! ");
         printf("             You need to use AddTask_iarsene_dst() in the on-the-fly reduced events mode!");
         return 0x0;
       }
     }
            
     mgr->AddTask(task);
      
     if(runMode==AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree) 
        mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
      
     if(runMode==AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents) 
       mgr->ConnectInput(task, 0, cReducedEvent);
  
     AliAnalysisDataContainer *cOutputHist = mgr->CreateContainer("jpsi2eeHistos", THashList::Class(),
                                                                  AliAnalysisManager::kOutputContainer, "dstAnalysisHistograms.root");
     mgr->ConnectOutput(task, 1, cOutputHist );
  }
  else {
    // nothing at the moment   
  }
  
  return task;
}


//_________________________________________________________________
void Setup(AliReducedAnalysisJpsi2eeMult* processor, TString prod /*="LHC10h"*/) {
  //
  // Configure the analysis task
  // Setup histograms, handlers, cuts, etc.
  //
  
  Int_t prefilterType = 4;
  
   processor->SetRunPairing(1);
   processor->SetRunEventMixing(1);
  
  TFile* trackletsFile = TFile::Open("/home/steffen/ALICE/trees/16L/testTask/multVtxMB_noPileupCut.root");
  AliReducedVarManager::SetMultiplicityProfile( (TH1*) trackletsFile->Get("multVtxMB") , 0 );
  trackletsFile->Close();
  
  
  TFile* trackletsFileOuter = TFile::Open("/home/steffen/ALICE/trees/16L/testTask/multVtxMB_outer.root");
  AliReducedVarManager::SetMultiplicityProfile( (TH1*) trackletsFileOuter->Get("multVtxMB"), 1 );
  trackletsFileOuter->Close();
  
  TFile* v0File = TFile::Open("/home/steffen/ALICE/trees/16L/testTask/v0Vtx.root");
  AliReducedVarManager::SetMultiplicityProfile( (TH1*) v0File->Get("slice_py_of_VZEROmult_ZVertexprof"), 2 );
  v0File->Close();
  
  
  
  TString runNumbers("256504;256506;256510;256512;256514;256552;256554;256556;256557;256560;256561;256562;256564;256565;256567;256589;256591;256592;256619;256620;256658;256676;256677;256681;256684;256691;256692;256694;256695;256697;256941;256942;256944;257011;257012;257021;257026;257028;257071;257077;257080;257082;257083;257084;257086;257092;257095;257100;257136;257137;257138;257139;257140;257141;257142;257144;257145;257204;257206;257209;257224;257260;257318;257320;257322;257330;257358;257364;257433;257457;257468;257474;257487;257488;257490;257491;257492;257530;257531;257537;257539;257540;257541;257560;257561;257562;257566;257587;257588;257590;257592;257594;257595;257601;257604;257605;257606;257630;257632;257635;257636;257642;257644;257682;257684;257685;257687;257688;257689;257691;257692;257694;257697;257724;257725;257727;257733;257734;257735;257737;257754;257757;257765;257773;257797;257798;257799;257800;257803;257804;257850;257851;257853;257855;257892;257893;257936;257937;257939;257957;257958;257960;257963;257979;257986;257989;257992;258003;258008;258012;258014;258017;258019;258039;258041;258042;258045;258048;258049;258053;258059;258060;258062;258063;258107;258108;258109;258113;258114;258117;258178;258197;258198;258202;258203;258204;258256;258257;258258;258270;258271;258273;258274;258278;258299;258301;258302;258303;258306;258307;258332;258336;258359;258387;258391;258393;258426;258452;258454;258456;258477;258499;258537;258962;258964;259088;259090;259091;259096;259099;259117;259118;259162;259164;259204;259257;259261;259263;259264;259269;259270;259271;259272;259273;259274;259302;259303;259305;259307;259334;259336;259339;259340;259341;259342;259378;259382;259388;259389;259394;259395;259396;259473;259477;259649;259650;259668;259697;259700;259703;259704;259705;259711;259713;259747;259748;259750;259751;259752;259756;259781;259788;259789;259822;259841;259842;259860;259866;259867;259868;259888;260010;260011;260014");
  
  AliReducedVarManager::SetRunNumbers( runNumbers );
  
   
  // Set event cuts
  AliReducedEventCut* evCut1 = new AliReducedEventCut("EventCuts","Event cuts");
  evCut1->AddCut(AliReducedVarManager::kVtxZ, -10.0, 10.0);
  evCut1->AddCut(AliReducedVarManager::kIsPhysicsSelection, 0.1, 2.);   // request physics selection
//   evCut1->AddCut(AliReducedVarManager::kIsSPDPileup5, .9, 1.1, kTRUE, AliReducedVarManager::kINT7Triggered, -.1 , .1 );   // in HM events: pileup cut with 5 contributors
//   evCut1->AddCut(AliReducedVarManager::kIsSPDPileup, .9, 1.1, kTRUE);  
   //   evCut1->AddCut(AliReducedVarManager::kINT7Triggered, 0.1, 2.);   // only MB events
//   evCut1->AddCut(AliReducedVarManager::kINT7Triggered, 0.1, 2., kTRUE);   // only HM events
  
  
  
//   evCut1->AddCut(AliReducedVarManager::kMultEstimatorOnlineV0M, 0., new TF1("V0cut", "-100. + 7.*x") ,   kTRUE, AliReducedVarManager::kVZEROTotalMult, 0., 10000. );   // VZERO mult vs. SPD n tracklets
  
   evCut1->AddCut(AliReducedVarManager::kIsSPDPileupMultBins, .9, 1.1, kTRUE);   // SPD pileup rejection in multiplicity bins
   evCut1->AddCut(AliReducedVarManager::kNVtxContributors, 0, 1., kTRUE);   // require 1 Vtx contributor
//   evCut1->AddCut(AliReducedVarManager::kMultEstimatorOnlineV0M, 0, 5500);
//    evCut1->AddCut( AliReducedVarManager::kNTracksPerTrackingStatus+ AliReducedVarManager::kTPCout,  0,new TF1("f","800 + 0.05*x") ,  kFALSE, AliReducedVarManager::kMultEstimatorOnlineV0M, 0., 5500.);
//   evCut1->AddCut(AliReducedVarManager::kVZERO, .9, 1.1, kTRUE);   // SPD pileup rejection in multiplicity bins
  
  
  
// evCut1->AddCut(AliReducedVarManager::kTZEROpileup, .9, 1.1, kTRUE  /*, AliReducedVarManager::kINT7Triggered, -.1 , .1 */);   // TZERO pileup rejection
//    evCut1->AddCut(AliReducedVarManager::kTZEROsatellite, .9, 1.1, kTRUE /*, AliReducedVarManager::kINT7Triggered, -.1 , .1 */);    // TZERO sattelite rejection
   
  processor->AddEventCut(evCut1);
  // Set track cuts

if(SYSTEMATICS){
    
  processor->AddTrackCut( SetTrackCuts("default") );
  processor->AddTrackCut( SetTrackCuts("strictPion1", 
  .9, 1., 3.,70.,3., 4., -2., 3., 3.75) ) ;
  processor->AddTrackCut( SetTrackCuts("strictPion2", 
  .9, 1., 3.,70.,3., 4., -2., 3., 4.) ) ;
  
  processor->AddTrackCut( SetTrackCuts("loosePion1", 
  .9, 1., 3.,70.,3., 4., -2., 3., 3.25) ) ;
  processor->AddTrackCut( SetTrackCuts("loosePion2", 
  .9, 1., 3.,70.,3., 4., -2., 3., 3.) ) ;
  
  processor->AddTrackCut( SetTrackCuts("strictProton1", 
  .9, 1., 3.,70.,3., 4., -2., 3., 3.5, 3.75) ) ;
  processor->AddTrackCut( SetTrackCuts("strictProton2", 
  .9, 1., 3.,70.,3., 4., -2., 3., 3.5, 4.) ) ;
  
  processor->AddTrackCut( SetTrackCuts("looseProton1", 
  .9, 1., 3.,70.,3., 4., -2., 3., 3.5, 3.25) ) ;
  processor->AddTrackCut( SetTrackCuts("looseProton2", 
  .9, 1., 3.,70.,3., 4., -2., 3., 3.5, 3.) ) ;
  
  processor->AddTrackCut( SetTrackCuts("strictElectronUp1", 
  .9, 1., 3.,70., 3., 4., -2., 2.75) ) ;
  processor->AddTrackCut( SetTrackCuts("strictElectronUp2", 
  .9, 1., 3.,70., 3., 4., -2., 2.5) ) ;
  
  processor->AddTrackCut( SetTrackCuts("looseElectronUp1", 
  .9, 1., 3.,70., 3., 4., -2., 3.25) ) ;
  processor->AddTrackCut( SetTrackCuts("looseElectronUp2", 
  .9, 1., 3.,70., 3., 4., -2., 3.5) ) ;
  
  processor->AddTrackCut( SetTrackCuts("strictElectronDown1", 
  .9, 1., 3.,70., 3., 4., -1.75) ) ;
  processor->AddTrackCut( SetTrackCuts("strictElectronDown2", 
  .9, 1., 3.,70., 3., 4., -1.5) ) ;
  
  processor->AddTrackCut( SetTrackCuts("looseElectronDown1", 
  .9, 1., 3.,70., 3., 4., -2.25) ) ;
  processor->AddTrackCut( SetTrackCuts("looseElectronDown2", 
  .9, 1., 3.,70., 3., 4., -2.5) ) ;
  
  
  
  processor->AddTrackCut( SetTrackCuts("strictChi21", 
  .9, 1., 3.,70., 3., 3.) ) ;
  processor->AddTrackCut( SetTrackCuts("strictChi22", 
  .9, 1., 3.,70., 3., 2.5) ) ;
  
  processor->AddTrackCut( SetTrackCuts("looseChi21", 
  .9, 1., 3.,70., 3., 5.) ) ;
  processor->AddTrackCut( SetTrackCuts("looseChi22", 
  .9, 1., 3.,70., 3., 10.) ) ;

  processor->AddTrackCut( SetTrackCuts("strictnclsTPC1", 
  .9, 1., 3.,80.) ) ;
  processor->AddTrackCut( SetTrackCuts("strictnclsTPC2", 
  .9, 1., 3.,90.) ) ;
  
  processor->AddTrackCut( SetTrackCuts("loosenclsTPC1", 
  .9, 1., 3.,60.) ) ;
  processor->AddTrackCut( SetTrackCuts("loosenclsTPC2", 
  .9, 1., 3.,50.) ) ;
  processor->AddTrackCut( SetTrackCuts("strictnclsITS1", 
  .9, 1., 3.,70., 4.) ) ;
  processor->AddTrackCut( SetTrackCuts("strictnclsITS2", 
  .9, 1., 3.,70., 5.) ) ;
  processor->AddTrackCut( SetTrackCuts("loosenclsITS1", 
  .9, 1., 3.,70., 2.) ) ;
  processor->AddTrackCut( SetTrackCuts("loosenclsITS2", 
  .9, 1., 3.,70., 1.) ) ;
  
  

  
  processor->AddTrackCut( SetTrackCuts("strictDCAxy1", 
  .9, .5) ) ;
  processor->AddTrackCut( SetTrackCuts("strictDCAxy2", 
  .9, .2) ) ;
  processor->AddTrackCut( SetTrackCuts("looseDCAxy1", 
  .9, 2.) ) ;
  processor->AddTrackCut( SetTrackCuts("looseDCAxy2", 
  .9, 3.) ) ;
  
  processor->AddTrackCut( SetTrackCuts("strictDCAz1", 
  .9, 1., 2.) ) ;
  processor->AddTrackCut( SetTrackCuts("strictDCAz2", 
  .9, 1., 1.) ) ;
  processor->AddTrackCut( SetTrackCuts("looseDCAz1", 
  .9, 1., 4.) ) ;
  processor->AddTrackCut( SetTrackCuts("looseDCAz2", 
  .9, 1., 5.) ) ;
  
  
  
  processor->AddTrackCut( SetTrackCuts("noKinkRejection", 
  .9, 1., 3.,70., 3., 4., -2., 3., 3.5, 3.5, kFALSE) ) ;
  
  
}
else{
  
  AliReducedTrackCut* trackCut1 = new AliReducedTrackCut("default","SPDany");
  trackCut1->AddCut(AliReducedVarManager::kPt, 1.,100.0);
  trackCut1->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  trackCut1->AddCut(AliReducedVarManager::kDcaXY, -1.,1.);
  trackCut1->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  trackCut1->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  trackCut1->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
  trackCut1->AddCut(AliReducedVarManager::kITSncls, 3., 9.);
  trackCut1->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.5, 1.0e+30);
  trackCut1->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.5, 1.0e+30);
  trackCut1->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2., 3.);
  trackCut1->SetRejectKinks();
  trackCut1->SetRequestITSrefit();
  trackCut1->SetRequestTPCrefit();
  trackCut1->SetRequestSPDany();
     processor->AddTrackCut(trackCut1);
     
     /*
  
  AliReducedTrackCut* trackCut2 = new AliReducedTrackCut("eta08","SPDany");
  trackCut2->AddCut(AliReducedVarManager::kPt, 1.,100.0);
  trackCut2->AddCut(AliReducedVarManager::kEta, -0.8, 0.8);
  trackCut2->AddCut(AliReducedVarManager::kDcaXY, -1.,1.);
  trackCut2->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  trackCut2->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  trackCut2->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
  trackCut2->AddCut(AliReducedVarManager::kITSncls, 3., 9.);
  trackCut2->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.5, 1.0e+30);
  trackCut2->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.5, 1.0e+30);
  trackCut2->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2., 3.);
  trackCut2->SetRejectKinks();
  trackCut2->SetRequestITSrefit();
  trackCut2->SetRequestTPCrefit();
  trackCut2->SetRequestSPDany();
     processor->AddTrackCut(trackCut2);
  
  AliReducedTrackCut* trackCut3 = new AliReducedTrackCut("eta07","SPDany");
  trackCut3->AddCut(AliReducedVarManager::kPt, 1.,100.0);
  trackCut3->AddCut(AliReducedVarManager::kEta, -0.7, 0.7);
  trackCut3->AddCut(AliReducedVarManager::kDcaXY, -1.,1.);
  trackCut3->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  trackCut3->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  trackCut3->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
  trackCut3->AddCut(AliReducedVarManager::kITSncls, 3., 9.);
  trackCut3->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.5, 1.0e+30);
  trackCut3->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.5, 1.0e+30);
  trackCut3->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2., 3.);
  trackCut3->SetRejectKinks();
  trackCut3->SetRequestITSrefit();
  trackCut3->SetRequestTPCrefit();
  trackCut3->SetRequestSPDany();
     processor->AddTrackCut(trackCut3);
  
  AliReducedTrackCut* trackCut4 = new AliReducedTrackCut("eta06","SPDany");
  trackCut4->AddCut(AliReducedVarManager::kPt, 1.,100.0);
  trackCut4->AddCut(AliReducedVarManager::kEta, -0.6, 0.6);
  trackCut4->AddCut(AliReducedVarManager::kDcaXY, -1.,1.);
  trackCut4->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  trackCut4->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  trackCut4->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
  trackCut4->AddCut(AliReducedVarManager::kITSncls, 3., 9.);
  trackCut4->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.5, 1.0e+30);
  trackCut4->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.5, 1.0e+30);
  trackCut4->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2., 3.);
  trackCut4->SetRejectKinks();
  trackCut4->SetRequestITSrefit();
  trackCut4->SetRequestTPCrefit();
  trackCut4->SetRequestSPDany();
     processor->AddTrackCut(trackCut4);
     
     
  AliReducedTrackCut* trackCut5 = new AliReducedTrackCut("eta01excl","SPDany");
  trackCut5->AddCut(AliReducedVarManager::kPt, 1.,100.0);
  trackCut5->AddCut(AliReducedVarManager::kEta, -0.9, 0.9 );
  trackCut5->AddCut(AliReducedVarManager::kEta, -0.1, 0.1, kTRUE);
  trackCut5->AddCut(AliReducedVarManager::kDcaXY, -1.,1.);
  trackCut5->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  trackCut5->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  trackCut5->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
  trackCut5->AddCut(AliReducedVarManager::kITSncls, 3., 9.);
  trackCut5->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.5, 1.0e+30);
  trackCut5->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.5, 1.0e+30);
  trackCut5->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2., 3.);
  trackCut5->SetRejectKinks();
  trackCut5->SetRequestITSrefit();
  trackCut5->SetRequestTPCrefit();
  trackCut5->SetRequestSPDany();
     processor->AddTrackCut(trackCut5);
     
     
  AliReducedTrackCut* trackCut6 = new AliReducedTrackCut("eta02excl","SPDany");
  trackCut6->AddCut(AliReducedVarManager::kPt, 1.,100.0);
  trackCut6->AddCut(AliReducedVarManager::kEta, -0.9, 0.9 );
  trackCut6->AddCut(AliReducedVarManager::kEta, -0.2, 0.2, kTRUE);
  trackCut6->AddCut(AliReducedVarManager::kDcaXY, -1.,1.);
  trackCut6->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  trackCut6->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  trackCut6->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
  trackCut6->AddCut(AliReducedVarManager::kITSncls, 3., 9.);
  trackCut6->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.5, 1.0e+30);
  trackCut6->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.5, 1.0e+30);
  trackCut6->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2., 3.);
  trackCut6->SetRejectKinks();
  trackCut6->SetRequestITSrefit();
  trackCut6->SetRequestTPCrefit();
  trackCut6->SetRequestSPDany();
     processor->AddTrackCut(trackCut6);
     
  AliReducedTrackCut* trackCut7 = new AliReducedTrackCut("eta10","SPDany");
  trackCut7->AddCut(AliReducedVarManager::kPt, 1.,100.0);
  trackCut7->AddCut(AliReducedVarManager::kEta, -1., 1. );
  trackCut7->AddCut(AliReducedVarManager::kDcaXY, -1.,1.);
  trackCut7->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  trackCut7->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  trackCut7->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
  trackCut7->AddCut(AliReducedVarManager::kITSncls, 3., 9.);
  trackCut7->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.5, 1.0e+30);
  trackCut7->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.5, 1.0e+30);
  trackCut7->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2., 3.);
  trackCut7->SetRejectKinks();
  trackCut7->SetRequestITSrefit();
  trackCut7->SetRequestTPCrefit();
  trackCut7->SetRequestSPDany();
     processor->AddTrackCut(trackCut7);
  */
     
     /*
  AliReducedTrackCut* trackCut2 = new AliReducedTrackCut("first","SPDfirst");
  trackCut2->AddCut(AliReducedVarManager::kPt, 1.,100.0);
  trackCut2->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  trackCut2->AddCut(AliReducedVarManager::kDcaXY, -1.,1.);
  trackCut2->AddCut(AliReducedVarManager::kDcaZ, -3.0,3.0);
  trackCut2->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  trackCut2->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
  trackCut2->AddCut(AliReducedVarManager::kITSncls, 3., 9.);
  trackCut2->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.5, 1.0e+30);
  trackCut2->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.5, 1.0e+30);
  trackCut2->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2., 3.);
  trackCut2->SetRejectKinks();
  trackCut2->SetRequestITSrefit();
  trackCut2->SetRequestTPCrefit();
  trackCut2->SetRequestSPDfirst();
     processor->AddTrackCut(trackCut2);*/
     
     
     
     /*
  
     
  
  AliReducedTrackCut* trackCut2 = new AliReducedTrackCut("lowPt","low pT");
  trackCut2->AddCut(AliReducedVarManager::kPt, 1.,2.0);
  trackCut2->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  trackCut2->AddCut(AliReducedVarManager::kDcaXY, -1.,1.);
  trackCut2->AddCut(AliReducedVarManager::kDcaZ, -3.,3.);
  trackCut2->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  trackCut2->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
  trackCut2->AddCut(AliReducedVarManager::kITSncls, 3., 9.);
  trackCut2->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.5, 1.0e+30);
  trackCut2->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.5, 1.0e+30);
  trackCut2->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2., 3.);
  trackCut2->SetRejectKinks();
  trackCut2->SetRequestITSrefit();
  trackCut2->SetRequestTPCrefit();
  trackCut2->SetRequestSPDany();
     processor->AddTrackCut(trackCut2);

     
     
  
  AliReducedTrackCut* trackCut3 = new AliReducedTrackCut("highPt","high pT");
  trackCut3->AddCut(AliReducedVarManager::kPt, 2.,100.0);
  trackCut3->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  trackCut3->AddCut(AliReducedVarManager::kDcaXY, -1.,1.);
  trackCut3->AddCut(AliReducedVarManager::kDcaZ, -3.,3.);
  trackCut3->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  trackCut3->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
  trackCut3->AddCut(AliReducedVarManager::kITSncls, 3., 9.);
  trackCut3->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.5, 1.0e+30);
  trackCut3->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.5, 1.0e+30);
  trackCut3->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2., 3.);
  trackCut3->SetRejectKinks();
  trackCut3->SetRequestITSrefit();
  trackCut3->SetRequestTPCrefit();
  trackCut3->SetRequestSPDany();
     processor->AddTrackCut(trackCut3);
     
     
     
     
     
  
//   AliReducedTrackCut* trackCut4 = new AliReducedTrackCut("strictChi2","tighter chi^2 cuts");
//   trackCut4->AddCut(AliReducedVarManager::kPt, 1.,100.0);
//   trackCut4->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
//   trackCut4->AddCut(AliReducedVarManager::kDcaXY, -1.,1.);
//   trackCut4->AddCut(AliReducedVarManager::kDcaZ, -3.,3.);
//   trackCut4->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
//   trackCut4->AddCut(AliReducedVarManager::kTPCchi2, 0., 1.4);
//   trackCut4->AddCut(AliReducedVarManager::kITSncls, 3., 9.);
//   trackCut4->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.5, 1.0e+30);
//   trackCut4->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.5, 1.0e+30);
//   trackCut4->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2., 3.);
//   trackCut4->SetRejectKinks();
//   trackCut4->SetRequestITSrefit();
//   trackCut4->SetRequestTPCrefit();
//   trackCut4->SetRequestSPDany();
//      processor->AddTrackCut(trackCut4);
     
     
    
  
  AliReducedTrackCut* trackCut5 = new AliReducedTrackCut("strictEta","|et| < 0.8");
  trackCut5->AddCut(AliReducedVarManager::kPt, 1.,100.0);
  trackCut5->AddCut(AliReducedVarManager::kEta, -0.8, 0.8);
  trackCut5->AddCut(AliReducedVarManager::kDcaXY, -1.,1.);
  trackCut5->AddCut(AliReducedVarManager::kDcaZ, -3.,3.);
  trackCut5->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  trackCut5->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
  trackCut5->AddCut(AliReducedVarManager::kITSncls, 3., 9.);
  trackCut5->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.5, 1.0e+30);
  trackCut5->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.5, 1.0e+30);
  trackCut5->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2., 3.);
  trackCut5->SetRejectKinks();
  trackCut5->SetRequestITSrefit();
  trackCut5->SetRequestTPCrefit();
  trackCut5->SetRequestSPDany();
     processor->AddTrackCut(trackCut5);
     
     
  
  AliReducedTrackCut* trackCut6 = new AliReducedTrackCut("TOF","TOF required");
  trackCut6->AddCut(AliReducedVarManager::kPt, .8,100.0);
  trackCut6->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  trackCut6->AddCut(AliReducedVarManager::kDcaXY, -1.,1.);
  trackCut6->AddCut(AliReducedVarManager::kDcaZ, -3.,3.);
  trackCut6->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
  trackCut6->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
  trackCut6->AddCut(AliReducedVarManager::kITSncls, 3., 9.);
  trackCut6->AddCut(AliReducedVarManager::kTOFbeta, 0.9, 1.1);
  trackCut6->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.5, 1.0e+30);
  trackCut6->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.5, 1.0e+30);
  trackCut6->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2., 3.);
  trackCut6->SetRejectKinks();
  trackCut6->SetRequestITSrefit();
  trackCut6->SetRequestTPCrefit();
  trackCut6->SetRequestSPDany();
     processor->AddTrackCut(trackCut6);
     
     
  
  AliReducedTrackCut* trackCut7 = new AliReducedTrackCut("strictCls","120 clusters required");
  trackCut7->AddCut(AliReducedVarManager::kPt, .8,100.0);
  trackCut7->AddCut(AliReducedVarManager::kEta, -0.9, 0.9);
  trackCut7->AddCut(AliReducedVarManager::kDcaXY, -1.,1.);
  trackCut7->AddCut(AliReducedVarManager::kDcaZ, -3.,3.);
  trackCut7->AddCut(AliReducedVarManager::kTPCncls, 120.,160.0);
  trackCut7->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
  trackCut7->AddCut(AliReducedVarManager::kITSncls, 3., 9.);
  trackCut7->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, 3.5, 1.0e+30);
  trackCut7->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.5, 1.0e+30);
  trackCut7->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -2., 3.);
  trackCut7->SetRejectKinks();
  trackCut7->SetRequestITSrefit();
  trackCut7->SetRequestTPCrefit();
  trackCut7->SetRequestSPDany();
     processor->AddTrackCut(trackCut7);*/
     
}
   
   
  if(prefilterType){
  
    AliReducedTrackCut* prefTrackCut2 = new AliReducedTrackCut("prefCutPt09","prefilter Pt selection");
    prefTrackCut2->AddCut(AliReducedVarManager::kPt, .1,100.0);
    prefTrackCut2->AddCut(AliReducedVarManager::kEta, -.9, .9);
    prefTrackCut2->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, -3., 3.0);
    prefTrackCut2->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, 3.5, 1e30);
    prefTrackCut2->SetRequestITSrefit();
    prefTrackCut2->SetRequestTPCrefit();
    // //    prefTrackCut2->SetRequestSPDany();
    // //   
        prefTrackCut2->AddCut(AliReducedVarManager::kTPCncls, 70.,160.0);
        prefTrackCut2->AddCut(AliReducedVarManager::kTPCchi2, 0., 4.);
    //   
      processor->AddPrefilterTrackCut(prefTrackCut2);
    
    if(prefilterType == 4){
      
      
      AliReducedTrackCut* prefPairCutPhiVDca = new AliReducedTrackCut("prefPairCutPhiVDca","prefilter pair cut");
      prefPairCutPhiVDca->AddCut(AliReducedVarManager::kMass, 0., .2 , kTRUE);
      processor->AddPrefilterPairCut(prefPairCutPhiVDca);
     
    }
      
      else{
      
    // //  // only apply on opposite-sign pairs
     AliReducedTrackCut* prefCutPairType = new AliReducedTrackCut("prefCutPairType","prefilter pair cut2");
     prefCutPairType->AddCut(AliReducedVarManager::kPairType, .9, 1.1, kTRUE );
     processor->AddPrefilterPairCut(prefCutPairType);
     
// // pair DCA cut -> if too large than random combination
      AliReducedTrackCut* prefPairCutDca = new AliReducedTrackCut("prefPairCutDca","prefilter pair cut");
      prefPairCutDca->AddCut(AliReducedVarManager::kPairDcaXYSqrt, 0.0, 1., kTRUE );
      processor->AddPrefilterPairCut(prefPairCutDca);
      
// wide mass cut: Dalitz electrons
//     AliReducedTrackCut* prefPairCutMassWide = new AliReducedTrackCut("prefPairCutMassWide","prefilter pair cut");
//     prefPairCutMassWide->AddCut(AliReducedVarManager::kMass, 0., 0.1  ,kTRUE);
//     processor->AddPrefilterPairCut( prefPairCutMassWide );
      }
      

    if(prefilterType == 1){
      
      
      AliReducedTrackCut* prefPairCutMassCorr = new AliReducedTrackCut("prefPairCutMassCorr","prefilter pair cut");
      prefPairCutMassCorr->AddCut( AliReducedVarManager::kMassDcaPtCorr, -0.01, .02 ,kTRUE  );
      processor->AddPrefilterPairCut(prefPairCutMassCorr);  
      
      
      AliReducedTrackCut* prefPairCutPhiV = new AliReducedTrackCut("prefPairCutPhiV","prefilter pair cut");
      prefPairCutPhiV->AddCut(   AliReducedVarManager::kPairPhiV, 0.0, 1. , kTRUE);
      processor->AddPrefilterPairCut(prefPairCutPhiV);
     
    }
    
    else if(prefilterType == 2){
      
      // // // wide cut on mass -> Dalitz electrons
//       

// // pair DCA cut -> if too large than random combination
      AliReducedTrackCut* prefPairCutDca = new AliReducedTrackCut("prefPairCutDca","prefilter pair cut");
      prefPairCutDca->AddCut(AliReducedVarManager::kPairDcaXYSqrt, 0.0, .5, kTRUE );
      processor->AddPrefilterPairCut(prefPairCutDca);

      
    }
    
    else if(prefilterType == 3){
     //    
// // // cut on phiV DCA dependent
      AliReducedTrackCut* prefPairCutPhiVDca = new AliReducedTrackCut("prefPairCutPhiVDca","prefilter pair cut");
      prefPairCutPhiVDca->AddCut( AliReducedVarManager::kMassDcaPtCorr, -0.02, .02 ,kTRUE/* AliReducedVarManager::kPairPhiV, 0., 1. */);
      prefPairCutPhiVDca->AddCut(AliReducedVarManager::kMass, 0., 0.1 , kTRUE ,AliReducedVarManager::kPairDcaXYSqrt, 0.0, .4  );
//       prefPairCutPhiVDca->AddCut(   AliReducedVarManager::kPairPhiV, 0., 1. , kTRUE, AliReducedVarManager::kPairDcaXYSqrt, .5, 1.);
//         prefPairCutPhiVDca->AddCut( AliReducedVarManager::kPairDcaXYSqrt, 0.0, .5, kTRUE, AliReducedVarManager::kPairPhiV, 1., 4. );
      processor->AddPrefilterPairCut(prefPairCutPhiVDca);
      
      
     
//     AliReducedTrackCut* prefPairCutMassWide = new AliReducedTrackCut("prefPairCutMassWide","prefilter pair cut");
//     prefPairCutMassWide->AddCut(AliReducedVarManager::kMass, 0., 0.1 , kTRUE    );
//     processor->AddPrefilterPairCut( prefPairCutMassWide );

    }
     
  }
     
  // Set pair cuts
   AliReducedTrackCut* pairCut1 = new AliReducedTrackCut("Ptpair","Pt pair selection");
   pairCut1->AddCut(AliReducedVarManager::kRap, -0.9,0.9);
//    pairCut1->AddCut(AliReducedVarManager::kEta, -0.9,0.9);
//    pairCut1->AddCut(AliReducedVarManager::kMassDcaPtCorr, -0.03,0.03);
   processor->AddPairCut(pairCut1);



  SetupHistogramManager(processor, prod, runNumbers);
  SetupMixingHandler(processor);
}


AliReducedTrackCut *SetTrackCuts(
  TString name = "default",
  Double_t eta = .9,
  Double_t dcaxy = 1.,
  Double_t dcaz = 3.,
  Double_t nclsTPC = 70.,
  Double_t nclsITS = 3.,
  Double_t chi2 = 4.,
  Double_t nsigmaElectronLow = -2.,
  Double_t nsigmaElectronUp = 3.,
  Double_t nSigmaPion = 3.5,
  Double_t nSigmaProton = 3.5,
  Bool_t rejectKinks = kTRUE,
  Int_t  spd = 2,
  Double_t pt =1.
){
  AliReducedTrackCut *trackCut1 = new AliReducedTrackCut(name.Data(), name.Data());
  trackCut1->AddCut(AliReducedVarManager::kPt, pt,100.0);
  trackCut1->AddCut(AliReducedVarManager::kEta, -1.*eta, eta);
  trackCut1->AddCut(AliReducedVarManager::kDcaXY, -1. * dcaxy, dcaxy);
  trackCut1->AddCut(AliReducedVarManager::kDcaZ, -1. * dcaz, dcaz);
  trackCut1->AddCut(AliReducedVarManager::kTPCncls, nclsTPC,160.0);
  trackCut1->AddCut(AliReducedVarManager::kITSncls, nclsITS,9.0);
  trackCut1->AddCut(AliReducedVarManager::kTPCchi2, 0., chi2);
  trackCut1->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kProton, nSigmaProton, 1.0e+30);
  trackCut1->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kPion, nSigmaPion, 1.0e+30);
  trackCut1->AddCut(AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron, nsigmaElectronLow, nsigmaElectronUp);
  if (rejectKinks ) trackCut1->SetRejectKinks();
  trackCut1->SetRequestITSrefit();
  trackCut1->SetRequestTPCrefit();
  if(spd ==1) trackCut1->SetRequestSPDfirst();
  if(spd ==2) trackCut1->SetRequestSPDany();
  if(spd ==3) trackCut1->SetRequestSPDboth();
  
  return trackCut1;
}




//_________________________________________________________________
void SetupMixingHandler(AliReducedAnalysisJpsi2eeMult* task) {
   //
   // setup the mixing handler
   //
   AliMixingHandler* handler = task->GetMixingHandler();
   handler->SetPoolDepth(200);
   handler->SetMixingThreshold(1.0);
   handler->SetDownscaleEvents(1);
   handler->SetDownscaleTracks(1);
//    handler->SetEventVariables(AliReducedVarManager::kSPDntracklets,AliReducedVarManager::kVtxZ,AliReducedVarManager::kINT7Triggered);
//     handler->SetEventVariables(AliReducedVarManager::kSPDntrackletsCorr,AliReducedVarManager::kVtxZ,AliReducedVarManager::kINT7Triggered);

    
    
   
   
#ifdef VZEROBINS
    const Int_t nTrkBinLimits = 10;
    Float_t trkBinLimits[nTrkBinLimits] ={ 0.,80., 160., 240.,  320., 400., 480., 560.,   700., 2000.};
    handler->SetEventVariables(AliReducedVarManager::kVZEROTotalMultCorr,AliReducedVarManager::kVZEROTotalMult,AliReducedVarManager::kINT7Triggered);
#else
    const Int_t nTrkBinLimits = 12;
    Float_t trkBinLimits[nTrkBinLimits] ={0., 1., 8., 16., 25., 35., 45.,  55., 65., 80., 110.  ,200.  };
//    handler->SetEventVariables(AliReducedVarManager::kSPDntrackletsCorr,AliReducedVarManager::kVZEROTotalMult,AliReducedVarManager::kINT7Triggered);
   
     handler->SetEventVariables(AliReducedVarManager::kSPDntrackletsCorr,AliReducedVarManager::kVtxZ,AliReducedVarManager::kINT7Triggered);
//    handler->SetEventVariables(AliReducedVarManager::kSPDntrackletsCorr,AliReducedVarManager::kTrackletsOverV0,AliReducedVarManager::kINT7Triggered);
#endif    
    
    
   const Int_t nZbinLimits = 2;//1;
   Float_t zLims[nZbinLimits] = {-10., 10. };
//     Double_t piOver8 = TMath::PiOver4() / 2;
//    const Int_t nZbinLimits = 9;
//    Float_t zLims[nZbinLimits] = {0., 1. * piOver8, 2. * piOver8, 3. * piOver8, 4. * piOver8, 5. * piOver8, 6. * piOver8, 7. * piOver8, 8. * piOver8 };
//    const Int_t nZbinLimits = 5;
//    Float_t zLims[nZbinLimits] = { 0., 410., 450., 800., 2000.  };
   
/*       
   const Int_t nTrkBinLimits =   8 ;
   Float_t trkBinLimits[nTrkBinLimits] ={ 0.,200,400.,600.,800.,1000.,1200.,2000. };*/
    
   
   
   
   const Int_t nTriggerbinLimits = 3;
   Float_t triggerLims[nTriggerbinLimits] = {0.  , 1., 2. };
   
   handler->SetCentralityLimits(nTrkBinLimits, trkBinLimits);
   handler->SetEventVertexLimits(nZbinLimits, zLims);
   handler->SetEventPlaneLimits(nTriggerbinLimits,  triggerLims);
}


//_________________________________________________________________
void SetupHistogramManager(AliReducedAnalysisJpsi2eeMult* task, TString prod , TString runNumbers ) {
  //
  // setup the histograms manager
  //
  AliReducedVarManager::SetDefaultVarNames();
  
  DefineHistograms(task, prod, runNumbers);
  
  AliReducedVarManager::SetUseVars(task->GetHistogramManager()->GetUsedVars());
}


//_________________________________________________________________
void DefineHistograms(AliReducedAnalysisJpsi2eeMult* task, TString prod, TString runNumbers /*="LHC10h"*/) {
  //
  // define histograms
  // NOTE: The DefineHistograms need to be called after the track cuts are defined because the name of the cuts
  //           are used in the histogram lists
   // NOTE: The name of the pair histogram lists for event mixing need to contain "PairMEPP", "PairMEPM", "PairMEMM", in their name, see below.
   //  TODO: make needed changes such that this becomes less prone to mistakes
   
  AliHistogramManager* man = task->GetHistogramManager(); 
   
    TString histClasses = "";
//   histClasses += "Event_BeforeCuts;";
   histClasses += "EventTag_BeforeCuts;";
//   histClasses += "EventTriggers_BeforeCuts;";
  
//   histClasses += "TrackStatusFlags_BeforeCuts;"; 
//   histClasses += "TrackITSclusterMap_BeforeCuts;";   
//   histClasses += "TrackTPCclusterMap_BeforeCuts;";   
  
  
  histClasses += "Event_AfterCuts;";   
   histClasses += "EventTag_AfterCuts;";
//   histClasses += "EventTriggers_AfterCuts;";
//   histClasses += "Track_BeforeCuts;";   
  for(Int_t i=0; i<task->GetNTrackCuts(); ++i) {
    TString cutName = task->GetTrackCutName(i);
//     histClasses += Form("Track_%s;", cutName.Data());
    histClasses += Form("Track+_%s;", cutName.Data());
    histClasses += Form("Track-_%s;", cutName.Data());
    
//   histClasses += Form("TrackStatusFlags_%s;", cutName.Data());
//   histClasses += Form("TrackITSclusterMap_%s;", cutName.Data());
//   histClasses += Form("TrackTPCclusterMap_%s;", cutName.Data());
    
    
    //histClasses += Form("Track_NoPrefilter_%s;", cutName.Data());
    //histClasses += Form("PairPrefilterSEPP_%s;", cutName.Data());
    //histClasses += Form("PairPrefilterSEPM_%s;", cutName.Data());
    //histClasses += Form("PairPrefilterSEMM_%s;", cutName.Data());
    //histClasses += Form("PairSEPP_%s;", cutName.Data());
    //histClasses += Form("PairSEPM_%s;", cutName.Data());
    //histClasses += Form("PairSEMM_%s;", cutName.Data());
    histClasses += Form("PairSEPP_%s;PairSEPM_%s;PairSEMM_%s;", cutName.Data(), cutName.Data(), cutName.Data());
    //histClasses += Form("PairMEPP_%s;", cutName.Data());
    //histClasses += Form("PairMEPM_%s;", cutName.Data());
    //histClasses += Form("PairMEMM_%s;", cutName.Data());
//     histClasses += Form("PairMEPP_%s;PairMEPM_%s;PairMEMM_%s;", cutName.Data(), cutName.Data(), cutName.Data());
    histClasses += Form("PairMEPM_%s;PairMEPP_%s;PairMEMM_%s;", cutName.Data(), cutName.Data(), cutName.Data());
  }
 
  TString classesStr(histClasses);
  TObjArray* arr=classesStr.Tokenize(";");
  
  cout << "Histogram classes included in the Histogram Manager" << endl;
  for(Int_t iclass=0; iclass<arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();
    man->AddHistClass(classStr.Data());
    // Event wise histograms
    if(classStr.Contains("Event_")) {
      
      cout << classStr.Data() << endl;
      
      if(!SYSTEMATICS){
      
//       man->AddHistogram(classStr.Data(),"Trk/V0", "ratio tracklets/ratio", kFALSE, 200, 0., 200., AliReducedVarManager::kTrackletsOverV0);
        
        
//       man->AddHistogram(classStr.Data(),"PtMax", "max. p_{T}", kFALSE, 200, 0., 20., AliReducedVarManager::kPtMax);
//       man->AddHistogram(classStr.Data(),"PhiOfPtMax", "phi of max. p_{T}", kFALSE, 800, -1., 7., AliReducedVarManager::kPhiOfptMax);
        
      
      man->AddHistogram(classStr.Data(),"RunID","Run IDs",kFALSE,274, 0, 274, AliReducedVarManager::kRunID,0,0,0,-1, 0,0,0,-1,runNumbers.Data());
      man->AddHistogram(classStr.Data(),"RunID_MB","Run IDs vs trigger",kFALSE, 
                     274, 0, 274, AliReducedVarManager::kRunID,2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered, 0,0,0,-1,runNumbers.Data());
      man->AddHistogram(classStr.Data(),"RunNo_VZEROmult","",kFALSE, 274, 0, 274, AliReducedVarManager::kRunID,
                                                                100, 0.0, 2000., AliReducedVarManager::kVZEROTotalMult, 0,0,0,-1,runNumbers.Data());

      
      
      
      man->AddHistogram(classStr.Data(),"T0_MB","",kTRUE,   2, 0.0, 2., AliReducedVarManager::kINT7Triggered, 2, 0, 2, AliReducedVarManager::kTZEROpileup );
      
      
      
      man->AddHistogram(classStr.Data(),"T0S_MB","",kTRUE,  2, 0.0, 2., AliReducedVarManager::kINT7Triggered , 2, 0, 2, AliReducedVarManager::kTZEROsatellite );
      
            
            
      
      
      man->AddHistogram(classStr.Data(),"RunNo_VZEROmultProfile","",kTRUE, 274, 0, 274, AliReducedVarManager::kRunID,
							2,0,2,AliReducedVarManager::kINT7Triggered,
                                                        100, 0.0, 1000., AliReducedVarManager::kVZEROTotalMult,runNumbers.Data());
      
      man->AddHistogram(classStr.Data(),"RunNo_NTracksTotalt","",0, 274, 0, 274, AliReducedVarManager::kRunID,
                                                                      100, 0.0, 1000., AliReducedVarManager::kNtracksTotal, 0,0,0,-1,runNumbers.Data() );
      
      man->AddHistogram(classStr.Data(),"RunNo_NTracksSelected","",0, 274, 0, 274, AliReducedVarManager::kRunID,
                                                                      20, 0.0, 20., AliReducedVarManager::kNtracksSelected, 0,0,0,-1,runNumbers.Data());
      
      man->AddHistogram(classStr.Data(),"RunNo_SPDtracklets","",0, 274, 0, 274, AliReducedVarManager::kRunID,
                                                                      200, 0.0, 200., AliReducedVarManager::kSPDntracklets, 0,0,0,-1,runNumbers.Data());
      man->AddHistogram(classStr.Data(),"RunNo_SPDtrackletsProfile","",kTRUE, 274, 0, 274, AliReducedVarManager::kRunID, 
							2,0,2,AliReducedVarManager::kINT7Triggered,
                                                                      200, 0.0, 200., AliReducedVarManager::kSPDntracklets, runNumbers.Data()    );
      
      man->AddHistogram(classStr.Data(),"RunNo_SPDtrackletsCorrected","",0, 274, 0, 274, AliReducedVarManager::kRunID,
                                                                      200, 0.0, 200., AliReducedVarManager::kSPDntrackletsCorr, 0,0,0,-1,runNumbers.Data());
      
      man->AddHistogram(classStr.Data(),"RunNo_SPDtrackletsCorrectedProfile","",kTRUE, 274, 0, 274, AliReducedVarManager::kRunID,
							2,0,2,AliReducedVarManager::kINT7Triggered,
                                                                      200, 0.0, 200., AliReducedVarManager::kSPDntrackletsCorr, runNumbers.Data());
      
      
      
      
      man->AddHistogram(classStr.Data(),"VtxX","Vtx X",kFALSE,300,-0.4,0.4,AliReducedVarManager::kVtxX);
      man->AddHistogram(classStr.Data(),"VtxY","Vtx Y",kFALSE,300,-0.4,0.4,AliReducedVarManager::kVtxY);
      man->AddHistogram(classStr.Data(),"VtxZ","Vtx Z",kFALSE,300,-15.,15.,AliReducedVarManager::kVtxZ);
      man->AddHistogram(classStr.Data(),"NTracksTotal","Number of total tracks per event",kFALSE,1000,0.,2000.,AliReducedVarManager::kNtracksTotal);
      man->AddHistogram(classStr.Data(),"NTracksSelected","Number of selected tracks per event",kFALSE,20,0.,20.,AliReducedVarManager::kNtracksSelected);


      
      
      
      man->AddHistogram(classStr.Data(),"SPDntracklets", "SPD #tracklets in |#eta|<1.0", kFALSE, 201, -1., 200., AliReducedVarManager::kSPDntracklets);
      man->AddHistogram(classStr.Data(),"SPDntrackletsCorrected", "corrected SPD #tracklets in |#eta|<1.0", kFALSE, 201, -1., 200., AliReducedVarManager::kSPDntrackletsCorr );
      man->AddHistogram(classStr.Data(),"SPDntrackletsCorrectedSmeared", "corrected SPD #tracklets in |#eta|<1.0, smeared", kFALSE, 201, -1., 200., AliReducedVarManager::kSPDntrackletsCorrSmear);
      
      
      
      man->AddHistogram(classStr.Data(),"SPDntrackletsOuter", "SPD #tracklets in outer eta region", kFALSE, 201, -1., 200., AliReducedVarManager::kSPDntrackletsOuter);
      man->AddHistogram(classStr.Data(),"SPDntrackletsOuterCorrected", "corrected SPD #tracklets in outer eta region", kFALSE, 201, -1., 200., AliReducedVarManager::kSPDntrackletsOuterCorr );
      man->AddHistogram(classStr.Data(),"SPDntrackletsOuterCorrectedSmeared", "corrected SPD #tracklets in outer eta region, smeared", kFALSE, 201, -1., 200., AliReducedVarManager::kSPDntrackletsOuterCorrSmear);

      
//       for(Int_t il=0; il<2; ++il)
//         man->AddHistogram(classStr.Data(), Form("SPDfiredChips_layer%d",il+1), Form("SPD fired chips in layer %d",il+1), 
// 			  kFALSE, 200, 0., 1000., AliReducedVarManager::kSPDFiredChips+il);
//       for(Int_t il=0; il<6; ++il)
//         man->AddHistogram(classStr.Data(), Form("ITSclusters_layer%d",il+1), Form("ITS clusters in layer %d",il+1), 
// 			  kFALSE, 100, 0., 1000., AliReducedVarManager::kITSnClusters+il);
//       man->AddHistogram(classStr.Data(), "SPDnSingleClusters", "SPD single clusters", 
// 			kFALSE, 500, 0., 500., AliReducedVarManager::kSPDnSingleClusters);	
//       man->AddHistogram(classStr.Data(),"VZEROmult", "", kFALSE,4000, 0.0, 4000., AliReducedVarManager::kVZEROTotalMult);
//       man->AddHistogram(classStr.Data(),"VZEROmult_v0m", "", kFALSE,
//                         1000, 0.0, 2000., AliReducedVarManager::kVZEROTotalMult,                        
//                         1000, 0.0, 10000., AliReducedVarManager::kMultEstimatorOnlineV0M
//       );
//       man->AddHistogram(classStr.Data(),"v0m", "", kFALSE,
//                         1000, 0.0, 10000., AliReducedVarManager::kMultEstimatorOnlineV0M
//       );
//      
//             man->AddHistogram(classStr.Data(),"VZEROmult_AC", "", kFALSE,
//                         1000, 0.0, 1000., AliReducedVarManager::kVZEROATotalMult,                        
//                         1000, 0.0, 1000., AliReducedVarManager::kVZEROCTotalMult
//       );
//      
//             man->AddHistogram(classStr.Data(),"v0m_AC", "", kFALSE,
//                         1000, 0.0, 10000., AliReducedVarManager::kMultEstimatorOnlineV0A,                        
//                         1000, 0.0, 10000., AliReducedVarManager::kMultEstimatorOnlineV0C
//       );
//       
      
      /*
      man->AddHistogram(classStr.Data(),"VZEROmult_VtxContributors", "", kFALSE, 
          200, 0.0, 2000., AliReducedVarManager::kVZEROTotalMult, 
          500, 0.0, 500., AliReducedVarManager::kNVtxContributors);
      
      
      
      man->AddHistogram(classStr.Data(),"VZEROmult_Ntracks", "", kFALSE, 
          100, 0.0, 1000., AliReducedVarManager::kVZEROTotalMult, 
          400, 0.0, 2000., AliReducedVarManager::kNtracksTotal);
      
      
      
      man->AddHistogram(classStr.Data(),"VZEROmult_NtracksTPCout", "", kFALSE, 
          100, 0.0, 2000., AliReducedVarManager::kVZEROTotalMult, 
          200, 0.0, 2000., AliReducedVarManager::kNTracksPerTrackingStatus+ AliReducedVarManager::kTPCout );
      
      man->AddHistogram(classStr.Data(),"SPDtracklets_NtracksTPCout", "", kFALSE, 
          200, 0.0, 200., AliReducedVarManager::kSPDntracklets, 
          2000, 0.0, 2000., AliReducedVarManager::kNTracksPerTrackingStatus+ AliReducedVarManager::kTPCout );
      
      man->AddHistogram(classStr.Data(),"VZEROmult_SPDtracklets", "", kFALSE, 
          200, 0.0, 2000., AliReducedVarManager::kVZEROTotalMult, 
          200, 0.0, 200., AliReducedVarManager::kSPDntracklets);
      
      */
      
      
      
      man->AddHistogram(classStr.Data(),"VZEROmult_SPDtracklets_MB", "", kFALSE, 
          200, 0.0, 2000., AliReducedVarManager::kVZEROTotalMult, 
          200, 0.0, 200., AliReducedVarManager::kSPDntracklets, 
            2,-.5 , 1.5, AliReducedVarManager::kINT7Triggered);

      man->AddHistogram(classStr.Data(),"VZEROmult_tracklets", "", kFALSE, 
          100, 0.0, 2000., AliReducedVarManager::kVZEROTotalMult, 
          200, 0.0, 200., AliReducedVarManager::kMultEstimatorSPDTracklets);
      
      
      man->AddHistogram(classStr.Data(),"VZEROmult_SPDclusters", "", kFALSE, 
          100, 0.0, 2000., AliReducedVarManager::kVZEROTotalMult, 
          200, 0.0, 800., AliReducedVarManager::kMultEstimatorSPDClusters);
      
      
   man->AddHistogram(classStr.Data(),"clusters_tracklets", "", kFALSE,
                     200, 0.0, 200.,  AliReducedVarManager::kMultEstimatorSPDTracklets,
                     200, 0.0, 1000.,  AliReducedVarManager::kMultEstimatorSPDClusters);
   
   man->AddHistogram(classStr.Data(),"clusters_v0m", "", kFALSE,
                     200, 0.0, 15000.,  AliReducedVarManager::kMultEstimatorOnlineV0M,
                     400, 0.0, 800.,  AliReducedVarManager::kMultEstimatorSPDClusters);
   
   man->AddHistogram(classStr.Data(),"tracklets_v0m", "", kFALSE,
                     200, 0.0, 15000.,  AliReducedVarManager::kMultEstimatorOnlineV0M,
                     200, 0.0, 200.,  AliReducedVarManager::kMultEstimatorSPDTracklets);
   
   man->AddHistogram(classStr.Data(),"SPDtracklets_v0m", "", kFALSE,
                     200, 0.0, 15000.,  AliReducedVarManager::kMultEstimatorOnlineV0M,
                     200, 0.0, 200.,  AliReducedVarManager::kSPDntracklets);
   
   man->AddHistogram(classStr.Data(),"tpcOut_v0m", "", kFALSE,
                     200, 0.0, 15000.,  AliReducedVarManager::kMultEstimatorOnlineV0M,
                     200, 0.0, 2000.,  AliReducedVarManager::kNTracksPerTrackingStatus+ AliReducedVarManager::kTPCout );
   
   man->AddHistogram(classStr.Data(),"tpcOut_tracklets", "", kFALSE,
                     200, 0.0, 200.,  AliReducedVarManager::kSPDntracklets,
                     400, 0.0, 4000.,  AliReducedVarManager::kNTracksPerTrackingStatus+ AliReducedVarManager::kTPCout );
   
   
      
      man->AddHistogram(classStr.Data(),"SPDtracklets_tracklets", "", kFALSE, 
          200, 0.0, 200., AliReducedVarManager::kSPDntracklets16, 
          200, 0.0, 200., AliReducedVarManager::kMultEstimatorSPDTracklets );

   
   
      
      
      
      man->AddHistogram(classStr.Data(),"Ntracks_SPDtracklets", "", kFALSE, 
          400, 0.0, 2000., AliReducedVarManager::kNtracksTotal, 
          200, 0.0, 200., AliReducedVarManager::kSPDntracklets);
      
      man->AddHistogram(classStr.Data(),"IsPhysicsSelection","Physics selection flag",kFALSE,
                        2,-0.5,1.5,AliReducedVarManager::kIsPhysicsSelection, 0,0.0,0.0,AliReducedVarManager::kNothing, 0,0.0,0.0,AliReducedVarManager::kNothing, "off;on");

            man->AddHistogram(classStr.Data(),"TZEROpileup","SPD pileup flag",kFALSE,
                        2,-0.5,1.5,AliReducedVarManager::kTZEROpileup, 0,0.0,0.0,AliReducedVarManager::kNothing, 0,0.0,0.0,AliReducedVarManager::kNothing, "off;on");
            
            
            man->AddHistogram(classStr.Data(),"TZEROsatellite","SPD pileup flag",kFALSE,
                        2,-0.5,1.5,AliReducedVarManager::kTZEROsatellite, 0,0.0,0.0,AliReducedVarManager::kNothing, 0,0.0,0.0,AliReducedVarManager::kNothing, "off;on");
            
            
            man->AddHistogram(classStr.Data(),"IsSPDpileup","SPD pileup flag",kFALSE,
                        2,-0.5,1.5,AliReducedVarManager::kIsSPDPileup, 0,0.0,0.0,AliReducedVarManager::kNothing, 0,0.0,0.0,AliReducedVarManager::kNothing, "off;on");           
            
            
            man->AddHistogram(classStr.Data(),"IsSPDpileupMultBins","SPD pileup flag (mult bins)",kFALSE,
                        2,-0.5,1.5,AliReducedVarManager::kIsSPDPileupMultBins, 0,0.0,0.0,AliReducedVarManager::kNothing, 0,0.0,0.0,AliReducedVarManager::kNothing, "off;on");
            
            
            man->AddHistogram(classStr.Data(),"IsSPDpileup5","SPD pileup flag",kFALSE,
                        2,-0.5,1.5,AliReducedVarManager::kIsSPDPileup5, 0,0.0,0.0,AliReducedVarManager::kNothing, 0,0.0,0.0,AliReducedVarManager::kNothing, "off;on");
            
            
            
            man->AddHistogram(classStr.Data(),"IsSPDpileup5VsTrigger","SPD pileup flag (5contr)",kFALSE,
                        2,-0.5,1.5,AliReducedVarManager::kIsSPDPileup5, 
                        2,-.5 , 1.5, AliReducedVarManager::kINT7Triggered );
            
            man->AddHistogram(classStr.Data(),"IsSPDpileupVsTrigger","SPD pileup flag",kFALSE,
                        2,-0.5,1.5,AliReducedVarManager::kIsSPDPileup, 
                        2,-.5 , 1.5, AliReducedVarManager::kINT7Triggered );
            
            
            
            
                  man->AddHistogram(classStr.Data(),"IsSPDpileup3:5","SPD pileup flag",kFALSE,
                        2,-0.5,1.5,AliReducedVarManager::kIsSPDPileup, 2,-0.5,1.5,AliReducedVarManager::kIsSPDPileup5);
            
      
      
      man->AddHistogram(classStr.Data(),"IsSPDpileup_SPDtracklets","",kFALSE,
                        2,-0.5,1.5,AliReducedVarManager::kIsSPDPileup,
                        200, 0.0, 200., AliReducedVarManager::kSPDntracklets); 
      man->AddHistogram(classStr.Data(),"IsSPDpileup_VZEROmult","",kFALSE,
                        2,-0.5,1.5,AliReducedVarManager::kIsSPDPileup,
                        100, 0.0, 1000., AliReducedVarManager::kVZEROTotalMult); 
      
      
      man->AddHistogram(classStr.Data(),"Zvertex_TPCout","",kFALSE,
                        20,-10.,10.,AliReducedVarManager::kVtxZ,
                        400, 0.0, 1200., AliReducedVarManager::AliReducedVarManager::kNTracksPerTrackingStatus+ AliReducedVarManager::kTPCout);
      man->AddHistogram(classStr.Data(),"Zvertex_TPCout_prof","",kTRUE,
                        20,-10.,10.,AliReducedVarManager::kVtxZ,
                        400, 0.0, 1200., AliReducedVarManager::AliReducedVarManager::kNTracksPerTrackingStatus+ AliReducedVarManager::kTPCout);
      
      
      
      
      man->AddHistogram(classStr.Data(),"Zvertex_SPDtracklets","",kFALSE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                        200,-10.,10.,AliReducedVarManager::kVtxZ,
                        200, 0.0, 200., AliReducedVarManager::kSPDntracklets);
      
      man->AddHistogram(classStr.Data(),"Zvertex_SPDtrackletsProf","",kTRUE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                        200,-10.,10.,AliReducedVarManager::kVtxZ,
                        200, 0.0, 200., AliReducedVarManager::kSPDntracklets);
      
      
      man->AddHistogram(classStr.Data(),"Zvertex_SPDtrackletsCorrected","",kFALSE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                        200,-10.,10.,AliReducedVarManager::kVtxZ,
                        200, 0.0, 200., AliReducedVarManager::kSPDntrackletsCorr);
      
      
      man->AddHistogram(classStr.Data(),"Zvertex_SPDtrackletsCorrectedSmeared","",kFALSE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                        200,-10.,10.,AliReducedVarManager::kVtxZ,
                        200, 0.0, 200., AliReducedVarManager::kSPDntrackletsCorrSmear);
      
      
      man->AddHistogram(classStr.Data(),"Zvertex_SPDtrackletsCorrectedProf","",kTRUE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                        200,-10.,10.,AliReducedVarManager::kVtxZ,
                        200, 0.0, 200., AliReducedVarManager::kSPDntrackletsCorr);
      
      
      man->AddHistogram(classStr.Data(),"Zvertex_SPDtrackletsCorrectedSmearedProf","",kTRUE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                        200,-10.,10.,AliReducedVarManager::kVtxZ,
                        200, 0.0, 200., AliReducedVarManager::kSPDntrackletsCorrSmear);
      
      
      /*
      man->AddHistogram(classStr.Data(),"Zvertex_SPDtrackletsOuter","",kFALSE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                        200,-10.,10.,AliReducedVarManager::kVtxZ,
                        200, 0.0, 200., AliReducedVarManager::kSPDntrackletsOuter);
      
      man->AddHistogram(classStr.Data(),"Zvertex_SPDtrackletsOuterProf","",kTRUE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                        200,-10.,10.,AliReducedVarManager::kVtxZ,
                        200, 0.0, 200., AliReducedVarManager::kSPDntrackletsOuter);*/
      
      
      /*
      
      man->AddHistogram(classStr.Data(),"Zvertex_SPDtrackletsOuterCorr","",kFALSE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                        200,-10.,10.,AliReducedVarManager::kVtxZ,
                        200, 0.0, 200., AliReducedVarManager::kSPDntrackletsOuterCorr);
      
      man->AddHistogram(classStr.Data(),"Zvertex_SPDtrackletsOuterCorrProf","",kTRUE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                        200,-10.,10.,AliReducedVarManager::kVtxZ,
                        200, 0.0, 200., AliReducedVarManager::kSPDntrackletsOuterCorr);
      
      
      
      
      man->AddHistogram(classStr.Data(),"Zvertex_SPDtrackletsOuterCorrSmeared","",kFALSE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                        200,-10.,10.,AliReducedVarManager::kVtxZ,
                        200, 0.0, 200., AliReducedVarManager::kSPDntrackletsOuterCorrSmear);
      
      man->AddHistogram(classStr.Data(),"Zvertex_SPDtrackletsOuterCorrSmearedProf","",kTRUE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                        200,-10.,10.,AliReducedVarManager::kVtxZ,
                        200, 0.0, 200., AliReducedVarManager::kSPDntrackletsOuterCorrSmear);
      
      
      
      */
      /*
      
      man->AddHistogram(classStr.Data(),"Zvertex_tracklets16","",kFALSE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                        200,-10.,10.,AliReducedVarManager::kVtxZ,
                        200, 0.0, 200., AliReducedVarManager::kSPDntracklets16);
      
      man->AddHistogram(classStr.Data(),"Zvertex_tracklets16Prof","",kTRUE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                        200,-10.,10.,AliReducedVarManager::kVtxZ,
                        200, 0.0, 200., AliReducedVarManager::kSPDntracklets16);*/
      
      
      /*
            
   man->AddHistogram(classStr.Data(),"v0m_ZVertexprof", "", kTRUE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                     200,-10.,10.,AliReducedVarManager::kVtxZ,
                     100, 0., 10000., AliReducedVarManager::kMultEstimatorOnlineV0M);
      
      
   
            
   man->AddHistogram(classStr.Data(),"VZEROmult_ZVertexprof", "", kTRUE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                     200,-10.,10.,AliReducedVarManager::kVtxZ,
                     100, 0., 10000., AliReducedVarManager::kVZEROTotalMult);
            
   man->AddHistogram(classStr.Data(),"VZEROmultCorr_ZVertexprof", "", kTRUE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                     200,-10.,10.,AliReducedVarManager::kVtxZ,
                     100, 0., 10000., AliReducedVarManager::kVZEROTotalMultCorr);
            
   man->AddHistogram(classStr.Data(),"VZEROmultCorrSmear_ZVertexprof", "", kTRUE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                     200,-10.,10.,AliReducedVarManager::kVtxZ,
                     100, 0., 10000., AliReducedVarManager::kVZEROTotalMultCorrSmear);*/
   
   
      
      
      
/*      
   man->AddHistogram(classStr.Data(),"MultEstV0M", "", kFALSE,1000, 0.0, 10000., AliReducedVarManager::kMultEstimatorOnlineV0M);
//   man->AddHistogram(classStr.Data(),"MultEstV0A", "", kFALSE,1000, 0.0, 6000.,  AliReducedVarManager::kMultEstimatorOnlineV0A);
//   man->AddHistogram(classStr.Data(),"MultEstV0C", "", kFALSE,1000, 0.0, 6000.,  AliReducedVarManager::kMultEstimatorOnlineV0C);
//   man->AddHistogram(classStr.Data(),"MultEstADM", "", kFALSE,1000, 0.0, 30000.,  AliReducedVarManager::kMultEstimatorADM);
//   man->AddHistogram(classStr.Data(),"MultEstADA", "", kFALSE,1000, 0.0, 10000.,  AliReducedVarManager::kMultEstimatorADA);
//   man->AddHistogram(classStr.Data(),"MultEstADC", "", kFALSE,1000, 0.0, 20000.,  AliReducedVarManager::kMultEstimatorADC);
   man->AddHistogram(classStr.Data(),"MultEstSPDclusters", "", kFALSE,1000, 0.0, 1000.,  AliReducedVarManager::kMultEstimatorSPDClusters);
   man->AddHistogram(classStr.Data(),"MultEstSPDtracklets", "", kFALSE,500, 0.0, 500.,  AliReducedVarManager::kMultEstimatorSPDTracklets);
   man->AddHistogram(classStr.Data(),"MultEstRefMult05", "", kFALSE,200, 0.0, 200.,  AliReducedVarManager::kMultEstimatorRefMult05);
   man->AddHistogram(classStr.Data(),"MultEstRefMult08", "", kFALSE,200, 0.0, 200.,  AliReducedVarManager::kMultEstimatorRefMult08);
      
   man->AddHistogram(classStr.Data(),"MultEstPercentileV0M", "", kFALSE,102, -1.0, 101., AliReducedVarManager::kMultEstimatorPercentileOnlineV0M);
//   man->AddHistogram(classStr.Data(),"MultEstPercentileV0A", "", kFALSE,102, -1.0, 101.,  AliReducedVarManager::kMultEstimatorPercentileOnlineV0A);
//   man->AddHistogram(classStr.Data(),"MultEstPercentileV0C", "", kFALSE,102, -1.0, 101.,  AliReducedVarManager::kMultEstimatorPercentileOnlineV0C);
//   man->AddHistogram(classStr.Data(),"MultEstPercentileADM", "", kFALSE,102, -1.0, 101.,  AliReducedVarManager::kMultEstimatorPercentileADM);
//   man->AddHistogram(classStr.Data(),"MultEstPercentileADA", "", kFALSE,102, -1.0, 101.,  AliReducedVarManager::kMultEstimatorPercentileADA);
//   man->AddHistogram(classStr.Data(),"MultEstPercentileADC", "", kFALSE,102, -1.0, 101.,  AliReducedVarManager::kMultEstimatorPercentileADC);
   man->AddHistogram(classStr.Data(),"MultEstPercentileSPDclusters", "", kFALSE,102, -1.0, 101.,  AliReducedVarManager::kMultEstimatorPercentileSPDClusters);
   man->AddHistogram(classStr.Data(),"MultEstPercentileSPDtracklets", "", kFALSE,102, -1.0, 101.,  AliReducedVarManager::kMultEstimatorPercentileSPDTracklets);
   man->AddHistogram(classStr.Data(),"MultEstPercentileRefMult05", "", kFALSE,102, -1.0, 101.,  AliReducedVarManager::kMultEstimatorPercentileRefMult05);
   man->AddHistogram(classStr.Data(),"MultEstPercentileRefMult08", "", kFALSE,102, -1.0, 101.,  AliReducedVarManager::kMultEstimatorPercentileRefMult08);
*/
  /* 
   man->AddHistogram(classStr.Data(),"MultEstV0A_ZVertex", "", kFALSE,
                     200,-10.,10.,AliReducedVarManager::kVtxZ,
                     100, 0., 10000.,  AliReducedVarManager::kMultEstimatorOnlineV0A);
  
   man->AddHistogram(classStr.Data(),"MultEstV0C_ZVertex", "", kFALSE,
                     200,-10.,10.,AliReducedVarManager::kVtxZ,
                     100, 0., 10000.,  AliReducedVarManager::kMultEstimatorOnlineV0C);

  */
  /* 
   man->AddHistogram(classStr.Data(),"MultEstADM_ZVertex", "", kFALSE,
                     200,-10.,10.,AliReducedVarManager::kVtxZ,
                     100, 0., 10000.,  AliReducedVarManager::kMultEstimatorADM);
  
   man->AddHistogram(classStr.Data(),"MultEstADA_ZVertex", "", kFALSE,
                     200,-10.,10.,AliReducedVarManager::kVtxZ,
                     100, 0., 10000.,  AliReducedVarManager::kMultEstimatorADA);
  
   man->AddHistogram(classStr.Data(),"MultEstADC_ZVertex", "", kFALSE,
                     200,-10.,10.,AliReducedVarManager::kVtxZ,
                     100, 0., 10000.,  AliReducedVarManager::kMultEstimatorADC);
   */
  /*
  
   man->AddHistogram(classStr.Data(),"MultEstSPDclusters_ZVertex", "", kFALSE,
                     200,-10.,10.,AliReducedVarManager::kVtxZ,
                     100, 0., 1000.,  AliReducedVarManager::kMultEstimatorSPDClusters);
  
   man->AddHistogram(classStr.Data(),"MultEstSPDtracklets_ZVertex", "", kFALSE,
                     200,-10.,10.,AliReducedVarManager::kVtxZ,
                     100, 0., 1000.,  AliReducedVarManager::kMultEstimatorSPDTracklets);
   
   man->AddHistogram(classStr.Data(),"MultEstRefMult05_ZVertex", "",kFALSE,
                     200,-10.,10.,AliReducedVarManager::kVtxZ, 
                     100, 0., 1000.,  AliReducedVarManager::kMultEstimatorRefMult05);
  
   man->AddHistogram(classStr.Data(),"MultEstRefMult08_ZVertex", "",kFALSE,
                     200,-10.,10.,AliReducedVarManager::kVtxZ, 
                     100, 0., 1000.,  AliReducedVarManager::kMultEstimatorRefMult08);
      

   
   man->AddHistogram(classStr.Data(),"MultEstRefMult08_Percentile", "", kFALSE,
                     200, 0.0, 200.,  AliReducedVarManager::kMultEstimatorRefMult08,
                     100, 0.0, 100.,  AliReducedVarManager::kMultEstimatorPercentileRefMult08);
   
   
   */
   
   
   
   
   
   
   
   man->AddHistogram(classStr.Data(),"SPDnTrackletsCorr_MB_Vertex", "", kFALSE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                     200, 0.0, 200.,  AliReducedVarManager::kSPDntrackletsCorr,
                     20, -10., 10., AliReducedVarManager::kVtxZ
                    );
   
   

      man->AddHistogram(classStr.Data(),"SPDnTrackletsCorrSmear_MB_Vertex", "", kFALSE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                     200, 0.0, 200.,  AliReducedVarManager::kSPDntrackletsCorrSmear,
                     20, -10., 10., AliReducedVarManager::kVtxZ
                    );
   
   
   man->AddHistogram(classStr.Data(),"SPDnTrackletsCorrSmear_MB_VZEROmult", "", kFALSE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                     200, 0.0, 200.,  AliReducedVarManager::kSPDntrackletsCorrSmear,
                     200, 0, 2000, AliReducedVarManager::kVZEROTotalMult
                    );
   
   
   
   
   
      
   
      
   man->AddHistogram(classStr.Data(),"VZEROmultCorr_MB", "", kFALSE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                     200, 0.0, 2000.,  AliReducedVarManager::kVZEROTotalMultCorr
                    );
   
   
   
   
      }
   
   man->AddHistogram(classStr.Data(),"SPDnTrackletsCorr_MB_VZEROmult", "", kFALSE,
                     2, 0.0, 2.,  AliReducedVarManager::kINT7Triggered,
                     200, 0.0, 200.,  AliReducedVarManager::kSPDntrackletsCorr,
                     200, 0, 2000, AliReducedVarManager::kVZEROTotalMult
                    );

   
    
    }  // end if className contains "Event"    
    
    
    
    
    else  if(classStr.Contains("EventTag_")) {
      
      man->AddHistogram(classStr.Data(), "EventTag", "Event Tag", kFALSE, 14,-1.0, 13.0, AliReducedVarManager::kEventTag);
      
    }
    // Track histograms
    else if(classStr.Contains("Track")  && !SYSTEMATICS  ) {
      cout << classStr.Data() << endl;
      man->AddHistogram(classStr.Data(), "Pt", "p_{T} distribution", kFALSE, 500, 0.0, 20.0, AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(), "Eta", "", kFALSE, 500, -1., 1., AliReducedVarManager::kEta);
      man->AddHistogram(classStr.Data(), "Phi", "", kFALSE, 100, 0.0, 6.3, AliReducedVarManager::kPhi);
      
        man->AddHistogram(classStr.Data(),"Phi_Eta","TPC N_{#sigma} electron vs. #eta",kFALSE,
                          100,-1.,1.,AliReducedVarManager::kEta,64,0.,6.4,AliReducedVarManager::kPhi);
      
      
      man->AddHistogram(classStr.Data(), "DCAxy", "DCAxy", kFALSE, 600, -3.0, 3.0, AliReducedVarManager::kDcaXY);
      
      
      man->AddHistogram(classStr.Data(), "DCAz", "DCAz", kFALSE, 500, -10.0, 10.0, AliReducedVarManager::kDcaZ);
      
      
      man->AddHistogram(classStr.Data(), "DCAz_DCAxy", "DCAz vs DCAxy", kFALSE, 
                        60, -3.0, 3.0, AliReducedVarManager::kDcaZ,
                        20, -1.0, 1.0, AliReducedVarManager::kDcaXY);
      
      
      man->AddHistogram(classStr.Data(), "DCAxy_Pt", "DCAxy_Pt", kFALSE, 90, 1.0, 10.0, AliReducedVarManager::kPt, 100, -1.0, 1.0, AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(), "DCAz_Pt", "DCAz_Pt", kFALSE, 90, 1.0, 10.0, AliReducedVarManager::kPt, 100, -3.0, 3.0, AliReducedVarManager::kDcaZ);
/*      
      man->AddHistogram(classStr.Data(), "DCAxy_P", "DCAxy_P", kFALSE, 90, 1.0, 10.0, AliReducedVarManager::kP, 200, -1.0, 1.0, AliReducedVarManager::kDcaXY);
      man->AddHistogram(classStr.Data(), "DCAz_P", "DCAz_P", kFALSE, 90, 1.0, 10.0, AliReducedVarManager::kP, 300, -3.0, 3.0, AliReducedVarManager::kDcaZ);*/
      
      

      
        man->AddHistogram(classStr.Data(),"ITSncls_Pt", "ITS nclusters vs. pT", kFALSE, 7,-0.5,6.5,AliReducedVarManager::kITSncls,
          90, 1.0, 10.0, AliReducedVarManager::kPt
          
        );


      
        man->AddHistogram(classStr.Data(),"ITSncls", "ITS nclusters", kFALSE, 7,-0.5,6.5,AliReducedVarManager::kITSncls);
//         man->AddHistogram(classStr.Data(),"Eta_Phi_ITSncls_prof","ITS <nclusters> vs (#eta,#phi)",kTRUE,
//                      36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 7, -0.5, 6.5, AliReducedVarManager::kITSncls);
	man->AddHistogram(classStr.Data(),"ITSchi2", "ITS #chi^{2}", kFALSE, 200,0.0,20.0, AliReducedVarManager::kITSchi2);
// 	man->AddHistogram(classStr.Data(),"Eta_Phi_ITSchi2_prof","ITS <#chi^{2}> vs (#eta,#phi)",kTRUE,
//                      36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 100, 0.0, 1000, AliReducedVarManager::kITSchi2);
      
        man->AddHistogram(classStr.Data(),"TPCncls","",kFALSE,160,-0.5,159.5,AliReducedVarManager::kTPCncls);
	man->AddHistogram(classStr.Data(),"TPCcrossedRows","", kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCcrossedRows);
	man->AddHistogram(classStr.Data(),"TPCnclsShared","",kFALSE, 160,-0.5,159.5,AliReducedVarManager::kTPCnclsShared);
	man->AddHistogram(classStr.Data(),"TPCchi2","",kFALSE, 510,-0.1,5.,AliReducedVarManager::kTPCchi2);
        man->AddHistogram(classStr.Data(),"Eta_Phi_TPCncls_prof","TPC <nclusters> vs (#eta,#phi)",kTRUE,
                     40, -1., 1., AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCncls);
        
        
        man->AddHistogram(classStr.Data(),"Eta_TPCncls","TPC <nclusters> vs #eta", kFALSE,
                     100, -1., 1., AliReducedVarManager::kEta, 161, -1., 160, AliReducedVarManager::kTPCncls);
// 	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCcrossedRows_prof","TPC <n crossed rows> vs (#eta,#phi)",kTRUE,
//                      36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCcrossedRows);
// 	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCnclsF_prof","TPC <nclusters findable> vs (#eta,#phi)",kTRUE,
//                      36, -0.9, 0.9, AliReducedVarManager::kEta, 90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 160, -0.5, 159.5, AliReducedVarManager::kTPCnclsF);
	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCchi2_prof","TPC <#chi^{2}> vs (#eta,#phi)",kTRUE,
                     100, -1., 1., AliReducedVarManager::kEta, 
		     90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 
		     160, -0.5, 159.5, AliReducedVarManager::kTPCchi2);
//         man->AddHistogram(classStr.Data(),"TPCsignal_Pin","TPC dE/dx vs. inner param P",kFALSE,
//                      100,0.0,10.0,AliReducedVarManager::kPin,200,-0.5,199.5,AliReducedVarManager::kTPCsignal);
// 	man->AddHistogram(classStr.Data(),"TPCsignalN","",kFALSE,160,-0.5,159.5,AliReducedVarManager::kTPCsignalN);
  
  
	man->AddHistogram(classStr.Data(),"Eta_Phi_TPCsignal_prof","TPC <dEdx> vs (#eta,#phi)",kTRUE,
                     100, -1., 1., AliReducedVarManager::kEta, 
		     90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 
                     160, -0.5, 159.5, AliReducedVarManager::kTPCsignal);   
  
  
  	man->AddHistogram(classStr.Data(),"Eta_Phi_nSigma_prof","TPC <nSigma> vs (#eta,#phi)",kTRUE,
                     100, -1., 1., AliReducedVarManager::kEta, 
		     90, 0.0, 2.0*TMath::Pi(), AliReducedVarManager::kPhi, 60, -3., 3.,
                      AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
  
  
        man->AddHistogram(classStr.Data(),"ratioClFind","ratio clusters/findable",kFALSE,
                     50,0.5,1.0,AliReducedVarManager::kTPCnclsRatio);
        /*
        man->AddHistogram(classStr.Data(),"ratioClFind_mult","ratio clusters/findable",kFALSE,
                     120,0.0,120.0,AliReducedVarManager::kSPDntrackletsCorrSmear,
                     50,0.5,1.0,AliReducedVarManager::kTPCnclsRatio);*/
  
        man->AddHistogram(classStr.Data(),"TPCnsigElectron_Pin","TPC N_{#sigma} electron vs. inner param P",kFALSE,
                     790,0.1,8.0,AliReducedVarManager::kPin,120,-3.0,3.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        
        man->AddHistogram(classStr.Data(),"TPCnsigElectron_Eta","TPC N_{#sigma} electron vs. #eta",kFALSE,
                          100,-1.,1.,AliReducedVarManager::kEta,120,-3.0,3.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        
        
        man->AddHistogram(classStr.Data(),"dEdx_Eta","dEdx vs. #eta",kFALSE,
                          100,-1.,1.,AliReducedVarManager::kEta,80,60.,100.,AliReducedVarManager::kTPCsignal);
        
        
        
        
        
        
                man->AddHistogram(classStr.Data(),"TPCnsigElectron_SPDnTracklets","",kFALSE,
                     60,0.0,120.0,AliReducedVarManager::kSPDntrackletsCorrSmear,60,-3.0,3.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);
        /*
                man->AddHistogram(classStr.Data(),"TPCnsigElectron_chi2","",kFALSE,
                     120,0.0,4.0,AliReducedVarManager::kTPCchi2,120,-3.0,3.0,AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron);*/
                
                
                man->AddHistogram(classStr.Data(),"dedx_SPDnTracklets","",kFALSE,
                     60,0.0,120.0,AliReducedVarManager::kSPDntrackletsCorrSmear,60,0.,100.,AliReducedVarManager::kTPCsignal);
                
//                 man->AddHistogram(classStr.Data(),"TPCncls_SPDnTracklets","",kFALSE,
//                      60,0.,120.,AliReducedVarManager::kSPDntrackletsCorrSmear,
//                      100,60.,160.,AliReducedVarManager::kTPCncls);
                
                /*
                man->AddHistogram(classStr.Data(),"TPCchi2_SPDnTracklets","",kFALSE,
                     60,0.0,120.0,AliReducedVarManager::kSPDntrackletsCorrSmear,100,0.,5.,AliReducedVarManager::kTPCchi2);
                
                
                man->AddHistogram(classStr.Data(),"eta_SPDnTracklets","",kFALSE,
                     60,0.0,120.0,AliReducedVarManager::kSPDntrackletsCorrSmear,180,-.9,.9,AliReducedVarManager::kEta);
                
                
                     man->AddHistogram(classStr.Data(),"sharedITS_SPDnTracklets","",kFALSE,
                     50,0.0,100.0,AliReducedVarManager::kSPDntrackletsCorrSmear,10,-1.,9., AliReducedVarManager::kITSnclsShared);
                                     
                     man->AddHistogram(classStr.Data(),"sharedRatioITS_SPDnTracklets","",kFALSE,
                     50,0.0,100.0,AliReducedVarManager::kSPDntrackletsCorrSmear,100,0.,1., AliReducedVarManager::kNclsSFracITS);
                
                     man->AddHistogram(classStr.Data(),"sharedRatio_SPDnTracklets","",kTRUE,
                     50,0.0,100.0,AliReducedVarManager::kSPDntrackletsCorrSmear,100,0.,1., AliReducedVarManager::kTPCnclsSharedRatio);*/

                     

                /*
                     man->AddHistogram(classStr.Data(),"sharedRatio_Vtx","",kTRUE,
                     20,-10.0,10.0,AliReducedVarManager::kVtxZ,100,0.,1., AliReducedVarManager::kTPCnclsSharedRatio);*/
                     
                
                man->AddHistogram(classStr.Data(),"chi2_eta","",kFALSE,
                      100,-1.,1.,AliReducedVarManager::kEta,
                     40,0.0,4.0,AliReducedVarManager::kTPCchi2);
                
                
                
                
                
                
      
        man->AddHistogram(classStr.Data(),"TOFbeta_P","TOF #beta vs P",kFALSE,
                     200,0.0,10.0,AliReducedVarManager::kP, 110,0.0,1.1,AliReducedVarManager::kTOFbeta);
      
        
        
        
// QA profiles vs. run



        
//       man->AddHistogram(classStr.Data(),"RunNo_ITSchi2Prof","",1, 274, 0, 274, AliReducedVarManager::kRunID,
//                                                                       100, 0.0, 100., AliReducedVarManager::kITSchi2,
//                        0,0,0,-1,runNumbers);
//       man->AddHistogram(classStr.Data(),"RunNo_TPCchi2Prof","",1, 274, 0, 274, AliReducedVarManager::kRunID,
//                                                                       100, 0.0, 10., AliReducedVarManager::kTPCchi2,
//                        0,0,0,-1,runNumbers);
      /*
      man->AddHistogram(classStr.Data(),"RunNo_goldenChi2Prof","",1, 274, 0, 274, AliReducedVarManager::kRunID,
                                                                      100, 0.0, 100., AliReducedVarManager::kChi2TPCConstrainedVsGlobal,
                       0,0,0,-1,runNumbers);*/
      
      
//       man->AddHistogram(classStr.Data(),"RunNo_TPCchi2","",0, 274, 0, 274, AliReducedVarManager::kRunID,
//                                                                       800, 0.0, 8., AliReducedVarManager::kTPCchi2,
//                        0,0,0,-1,runNumbers);
//        man->AddHistogram(classStr.Data(),"RunNo_TPCnclsProf","",1,274, 0, 274, AliReducedVarManager::kRunID,
//                                                                       160, 0.0, 160., AliReducedVarManager::kTPCncls,
//                        0,0,0,-1,runNumbers);
//        man->AddHistogram(classStr.Data(),"RunNo_EtaProf","",1, 274, 0, 274, AliReducedVarManager::kRunID,
//                                                                       180, -.9, .9, AliReducedVarManager::kEta,
//                        0,0,0,-1,runNumbers);
//                                                                       
//        man->AddHistogram(classStr.Data(),"RunNo_PhiProf","",1, 274, 0, 274, AliReducedVarManager::kRunID,
//                                                                       630, 0.0, 6.3, AliReducedVarManager::kPhi,
//                        0,0,0,-1,runNumbers);
//                                                                       
//        man->AddHistogram(classStr.Data(),"RunNo_PtProf","",1, 274, 0, 274, AliReducedVarManager::kRunID,
//                                                                       100, 0.0, 10., AliReducedVarManager::kPt,
//                        0,0,0,-1,runNumbers);
//        
//        man->AddHistogram(classStr.Data(),"RunNo_kTPCsignalProf","",1, 274, 0, 274, AliReducedVarManager::kRunID,
//                                                                       100, 0.0, 200., AliReducedVarManager::kTPCsignal,
//                        0,0,0,-1,runNumbers);
//        
//        man->AddHistogram(classStr.Data(),"RunNo_kTPCsignalNProf","",1, 274, 0, 274, AliReducedVarManager::kRunID,
//                                                                       160, 0.0, 160., AliReducedVarManager::kTPCsignalN,
//                        0,0,0,-1,runNumbers);
//        
//        man->AddHistogram(classStr.Data(),"RunNo_NsigmaElectronProf","",1, 274, 0, 274, AliReducedVarManager::kRunID,
//                                                                       100, -4, 4., AliReducedVarManager::kTPCnSig+AliReducedVarManager::kElectron,
//                        0,0,0,-1,runNumbers);
        
        
    }  // end if "TrackQA"
        
    // Histograms for pairs
    else if(classStr.Contains("Pair")) {
      
      
      
      
    const Int_t kNMassBins =  151;
    Double_t massBins[kNMassBins];
    for(Int_t i=0; i<kNMassBins; ++i) massBins[i] = 0.0 + i*0.04; 
    
    const Int_t kNPtBins = 10;
    Double_t ptBins[kNPtBins] = { 0., 1., 2., 3., 4., 5., 7., 10., 20., 100. };
    
    
//     const Int_t kNVtxBins = 7;
//     Double_t vtxBins[kNVtxBins] = {-10. ,-8., -5., 0., 5.,8., 10. };
    
    
    

   Double_t piOver16 = TMath::PiOver4()/4;
   const Int_t kNangleBins = 17;
   Double_t angleBins[kNangleBins] = { 0., 1. * piOver16, 2. * piOver16, 3. * piOver16, 4. * piOver16, 5. * piOver16, 6. * piOver16, 7. * piOver16, 8. * piOver16 , 9. * piOver16 , 10. * piOver16 , 11. * piOver16 , 12. * piOver16 , 13. * piOver16 , 14. * piOver16 , 15. * piOver16 , 16. * piOver16 };

//    const Int_t kNVtxBins = 19;
//    Double_t vtxBins[kNVtxBins] = {  -.9,-.8, -.7, -.6,-.5, -.4, -.3,-.2, -.1,0., .1,.2, .3,.4, .5,.6, .7,.8, .9  };
   const Int_t kNVtxBins = 2;
   Double_t vtxBins[kNVtxBins] = {  -10.,  10.  };
   
   const Int_t kNVtxBins2 = 6;
   Double_t vtxBins2[kNVtxBins2] = { 0., 5.,10.,15.,20.,100. };
   
//    const Int_t kNVtxBins2 = 3;
//    Double_t vtxBins2[kNVtxBins2] = {  -1.8,0., 1.8  };
    
    
        
    const Int_t kNVzeroBins = 3;
    Double_t vzeroBins[kNVzeroBins] = {0.,380.,1000.};
    
    
#ifdef VZEROBINS
   const Int_t nTrkBinLimits =   10;
   Double_t trkBinLimits[nTrkBinLimits] ={  0.,80., 160., 240.,  320., 400., 480., 560.,   700., 2000. };
#else
   const Int_t nTrkBinLimits =   12 ;
   Double_t trkBinLimits[nTrkBinLimits] ={ 0., 1., 8., 16., 25., 35., 45.,  55., 65., 80., 110.  ,200. };
#endif
//    const Int_t nTrkBinLimits =   8 ;
//    Double_t trkBinLimits[nTrkBinLimits] ={ 0.,200,400.,600.,800.,1000.,1200.,2000. };
    
    
    
    const Int_t kNMBBins = 3;
    Double_t mbBins[kNMBBins] = {0.,1.,2.};
    
    
    
    const Int_t nVars = 5;
    

#ifdef VZEROBINS
      Int_t vars[nVars] = {AliReducedVarManager::kMass, AliReducedVarManager::kPt, 
          AliReducedVarManager::kVZEROTotalMult, AliReducedVarManager::kVZEROTotalMultCorr, AliReducedVarManager::kINT7Triggered };
    
    #else
//       Int_t vars[nVars] = {AliReducedVarManager::kMass, AliReducedVarManager::kPt, 
//           AliReducedVarManager::kVtxZ, AliReducedVarManager::kSPDntrackletsCorr, AliReducedVarManager::kINT7Triggered };
//           
       Int_t vars[nVars] = {AliReducedVarManager::kMass, AliReducedVarManager::kPt, 
           AliReducedVarManager::kVtxZ, AliReducedVarManager::kSPDntrackletsCorr, AliReducedVarManager::kINT7Triggered };
//       Int_t vars[nVars] = {AliReducedVarManager::kMass, AliReducedVarManager::kPt, 
//           AliReducedVarManager::kTrackletsOverV0, AliReducedVarManager::kSPDntrackletsCorr, AliReducedVarManager::kINT7Triggered, AliReducedVarManager::kPairOpeningAngle};
#endif
//      Int_t vars[nVars] = {AliReducedVarManager::kMass, AliReducedVarManager::kPt, 
//          AliReducedVarManager::kVtxZ, AliReducedVarManager::kSPDntracklets, AliReducedVarManager::kINT7Triggered };
       

     TArrayD pairHistBinLimits[nVars];
     pairHistBinLimits[0] = TArrayD(kNMassBins,massBins);
     pairHistBinLimits[1] = TArrayD(kNPtBins,ptBins);
     pairHistBinLimits[2] = TArrayD(kNVtxBins,vtxBins);
     pairHistBinLimits[3] = TArrayD(nTrkBinLimits,trkBinLimits);
     pairHistBinLimits[4] = TArrayD(kNMBBins,mbBins);
//       pairHistBinLimits[5] = TArrayD(kNangleBins,angleBins);
    
      cout << classStr.Data() << endl;
    
      man->AddHistogram(classStr.Data(), "PairInvMass", "Differential pair inv. mass", nVars, vars, pairHistBinLimits);
      
      
      if(  !SYSTEMATICS && (1|| classStr.Contains("PairSE") || classStr.Contains("PairPrefilterSE") ) ) {
        
        man->AddHistogram(classStr.Data(), "PairType", "Pair type", kFALSE, 4, -0.5, 3.5, AliReducedVarManager::kPairType);
        man->AddHistogram(classStr.Data(), "Mass", "Invariant mass", kFALSE, 1100, 0., 11.0, AliReducedVarManager::kMass);
//         man->AddHistogram(classStr.Data(), "CorrectedMass", "Invariant mass, DCA corrected", kFALSE, 500, 0., 5.0, AliReducedVarManager::kMassDcaPtCorr);
        man->AddHistogram(classStr.Data(), "Pt", "", kFALSE, 200, 0.0, 20.0, AliReducedVarManager::kPt);
        man->AddHistogram(classStr.Data(), "Rapidity", "Rapidity", kFALSE, 300, -1.5, 1.5, AliReducedVarManager::kRap);
        man->AddHistogram(classStr.Data(), "Eta", "Eta", kFALSE, 240, -1.2, 1.2, AliReducedVarManager::kEta);
        man->AddHistogram(classStr.Data(), "Phi", "Azimuthal distribution", kFALSE, 315, 0.0, 6.3, AliReducedVarManager::kPhi);
        man->AddHistogram(classStr.Data(), "Phiv", "Phiv distribution", kFALSE, 630, 0., 3.15, AliReducedVarManager::kPairPhiV);
        
        
                
        man->AddHistogram(classStr.Data(), "PairDcaSqrt", "", kFALSE,
                          
                          400, 0.,2., AliReducedVarManager::kPairDcaXYSqrt);
        
        man->AddHistogram(classStr.Data(), "PairDca", "", kFALSE,
                          400, 0.,4., AliReducedVarManager::kPairDcaXY);
        
        
//         man->AddHistogram(classStr.Data(), "OpeningAngle", "opening angle distribution", kFALSE, 315, 0., 3.15, AliReducedVarManager::kPairOpeningAngle);
        
//
//         man->AddHistogram(classStr.Data(), "Phiv_Mass", "", kFALSE,
//                           316, -0.01, 3.15, AliReducedVarManager::kPairPhiV,
//                           300, 0.,6., AliReducedVarManager::kMass);
        
        
        man->AddHistogram(classStr.Data(), "Rap_Mass", "", kFALSE,
                          200, -1, 1., AliReducedVarManager::kRap,
                          300, 0.,6., AliReducedVarManager::kMass);
        
        man->AddHistogram(classStr.Data(), "Eta_Mass", "", kFALSE,
                          200, -4, 4., AliReducedVarManager::kEta,
                          300, 0.,6., AliReducedVarManager::kMass);
        
        man->AddHistogram(classStr.Data(), "Eta_Pt", "", kFALSE,
                          200, -2, 2., AliReducedVarManager::kEta,
                          200, 0.,20., AliReducedVarManager::kPt);
        
        
        
        
//         man->AddHistogram(classStr.Data(), "Phiv_OpeningAngle", "", kFALSE,
//                           316, -0.01, 3.15, AliReducedVarManager::kPairPhiV,
//                           200, 0.,.2, AliReducedVarManager::kPairOpeningAngle);
        
//                 man->AddHistogram(classStr.Data(), "Mass_CorrectedMass", "", kFALSE,
//                           200, 2.3, 3.3, AliReducedVarManager::kMass,
//                           200, 2.3, 3.3, AliReducedVarManager::kMassDcaPtCorr);
        
        
        
//         man->AddHistogram(classStr.Data(), "Phiv_OpeningAngleCorr", "", kFALSE,
//                           315, 0., 3.15, AliReducedVarManager::kPairPhiV,
//                           400, -.1,.3, AliReducedVarManager::kOpAngDcaPtCorr);
//         
//         
//         
//         
//         man->AddHistogram(classStr.Data(), "Phiv_PairDcaSqrt", "", kFALSE,
//                           315, 0., 3.15, AliReducedVarManager::kPairPhiV,
//                           400, 0.,1., AliReducedVarManager::kPairDcaXYSqrt);
//         
//         man->AddHistogram(classStr.Data(), "Phiv_PairDca", "", kFALSE,
//                           315, 0., 3.15, AliReducedVarManager::kPairPhiV,
//                           400, 0.,1., AliReducedVarManager::kPairDcaXY);
        
        
        
        
        
        man->AddHistogram(classStr.Data(), "OpeningAngleMass", "", kFALSE,
                          160, 0, 3.2, AliReducedVarManager::kPairOpeningAngle,
                          300, 0., 6., AliReducedVarManager::kMass);
        
        man->AddHistogram(classStr.Data(), "PtMass", "", kFALSE,
                          200, 0., 10., AliReducedVarManager::kPt,
                          300, 0., 6., AliReducedVarManager::kMass);


        man->AddHistogram(classStr.Data(), "OpeningAnglePt", "", kFALSE,
                          160, 0, 3.2, AliReducedVarManager::kPairOpeningAngle,
                          200, 0., 10., AliReducedVarManager::kPt);
        
        /*
        
        
                man->AddHistogram(classStr.Data(), "Phiv_CorrectedMass", "", kFALSE,
                          316, -0.01, 3.15, AliReducedVarManager::kPairPhiV,
                          200, -.1,.3, AliReducedVarManager::kMassDcaPtCorr);*/
  
  /*
      man->AddHistogram(classStr.Data(),"RunNo_Pt","",1, 274,0,274, AliReducedVarManager::kRunID,
                                                                      100, 0.0, 10., AliReducedVarManager::kPt);
      man->AddHistogram(classStr.Data(),"RunNo_InvMassProfile","",1, 274,0,274, AliReducedVarManager::kRunID,
                                                                      150, 0.0, 6., AliReducedVarManager::kMass);*/
  
/*
    man->AddHistogram(classStr.Data(),"RunNo_InvMass","",0, 274,0,274, AliReducedVarManager::kRunID,
                                                                      125, 0.0, 5., AliReducedVarManager::kMass);
      man->AddHistogram(classStr.Data(),"RunNo_Rapidity","",1, 274,0,274, AliReducedVarManager::kRunID,
                                                                      120, -1.2, 1.2, AliReducedVarManager::kRap);
      man->AddHistogram(classStr.Data(),"RunNo_Phi","",1, 274,0,274, AliReducedVarManager::kRunID,
                                                                     315, 0.0, 6.3, AliReducedVarManager::kPhi);
      
      
      
      man->AddHistogram(classStr.Data(),"Zvertex_InvMass","",1, 20, -10, 10, AliReducedVarManager::kVtxZ,
                                                                      150, 0.0, 6., AliReducedVarManager::kMass);
      
      man->AddHistogram(classStr.Data(),"SPDtracklets_InvMass","",1, 24, 0., 120, AliReducedVarManager::kSPDntracklets,
                                                                      150, 0.0, 6., AliReducedVarManager::kMass);
      
      
      man->AddHistogram(classStr.Data(),"Rapidity_InvMass","",1, 20, -1., 1., AliReducedVarManager::kRap,
                                                                      150, 0.0, 6., AliReducedVarManager::kMass);
      man->AddHistogram(classStr.Data(),"Eta_InvMassProf","",1, 80, -4., 4., AliReducedVarManager::kEta,
                                                                      150, 0.0, 6., AliReducedVarManager::kMass);
      man->AddHistogram(classStr.Data(),"Eta_InvMass","",0, 80, -4., 4., AliReducedVarManager::kEta,
                                                                      150, 0.0, 6., AliReducedVarManager::kMass);
      
      man->AddHistogram(classStr.Data(),"Eta_SPDtracklets","",0, 80, -4., 4., AliReducedVarManager::kEta,
                                                                      24, 0.0, 120., AliReducedVarManager::kSPDntracklets);
      
      
      man->AddHistogram(classStr.Data(),"Zvertex_SPDtracklets","",1, 20, -10, 10, AliReducedVarManager::kVtxZ,
                                                                      23, 0.0, 120., AliReducedVarManager::kSPDntracklets);
      
      man->AddHistogram(classStr.Data(),"Phi_InvMass","",1, 63, 0., 6.3, AliReducedVarManager::kPhi,
                                                                      150, 0.0, 6., AliReducedVarManager::kMass);
      
*/      
      /*
        man->AddHistogram(classStr.Data(), "CorrOpAng_pairDCAxySqrt_OneOverSqrtPt", "", 1,
                          200, 0., 1., AliReducedVarManager::kPairDcaXYSqrt,
                          50, 0.1, 0.7, AliReducedVarManager::kOneOverSqrtPt,
                          200, -1., .3, AliReducedVarManager::kOpAngDcaPtCorr);
        
        man->AddHistogram(classStr.Data(), "CorrOpAng_OneOverSqrtPt", "", 0,
                          200, -.1, .3, AliReducedVarManager::kOpAngDcaPtCorr,
                          50, 0.1, 0.7, AliReducedVarManager::kOneOverSqrtPt);
        
        man->AddHistogram(classStr.Data(), "CorrOpAng_pairDCAxySqrt", "", 0,
                          200, -.1, .3, AliReducedVarManager::kOpAngDcaPtCorr,
                          400, 0., 1., AliReducedVarManager::kPairDcaXYSqrt);
        
        
        
        man->AddHistogram(classStr.Data(), "OpAng_pairDCAxySqrt_OneOverSqrtPt", "", 1,
                          400, 0., 2., AliReducedVarManager::kPairDcaXYSqrt,
                          50, 0.1, 0.7, AliReducedVarManager::kOneOverSqrtPt,
                          200, -1., .3, AliReducedVarManager::kPairOpeningAngle);
        
        man->AddHistogram(classStr.Data(), "OpAng_OneOverSqrtPt", "", 0,
                          200, -.1, .3, AliReducedVarManager::kPairOpeningAngle,
                          50, 0.1, 0.7, AliReducedVarManager::kOneOverSqrtPt);
        
        man->AddHistogram(classStr.Data(), "OpAng_pairDCAxySqrt", "", 0,
                          200, -.1, .3, AliReducedVarManager::kPairOpeningAngle,
                          400, 0., 1., AliReducedVarManager::kPairDcaXYSqrt);
        
        
        
        
               
        man->AddHistogram(classStr.Data(), "CorrMass_pairDCAxySqrt_Pt", "", 1,
                          400, 0.,2., AliReducedVarManager::kPairDcaXYSqrt,
                           50, 0.6, 5.6, AliReducedVarManager::kPt,
                          300, -.1, .2, AliReducedVarManager::kMassDcaPtCorr);
        
               
        man->AddHistogram(classStr.Data(), "CorrMass_Pt", "", 0,
                          300, -.1, .2, AliReducedVarManager::kMassDcaPtCorr,
                           50, 0.6, 5.6, AliReducedVarManager::kPt);
        
               
        man->AddHistogram(classStr.Data(), "CorrMass_pairDCAxySqrt", "", 0,
                          300, -.1, .2, AliReducedVarManager::kMassDcaPtCorr,
                          400, 0., 1., AliReducedVarManager::kPairDcaXYSqrt);
               
        
        
        
        man->AddHistogram(classStr.Data(), "Mass_pairDCAxySqrt_Pt", "", 1,
                          400, 0., 2., AliReducedVarManager::kPairDcaXYSqrt,
                           50, 0.6, 5.6, AliReducedVarManager::kPt,
                          300, -.1, .2, AliReducedVarManager::kMass);
        
               
        man->AddHistogram(classStr.Data(), "Mass_Pt", "", 0,
                          300, -.1, .2, AliReducedVarManager::kMass,
                           50, 0.6, 5.6, AliReducedVarManager::kPt);
        
               
        man->AddHistogram(classStr.Data(), "Mass_pairDCAxySqrt", "", 0,
                          300, -.1, .2, AliReducedVarManager::kMass,
                          400, 0., 1., AliReducedVarManager::kPairDcaXYSqrt);
        
        */
        
        
        
  
      }   // end if "QA"
      
    }   // end if for Pair classes of histograms
  }  // end loop over histogram classes
}
