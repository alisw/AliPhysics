// event plane calibration task
// author: Jaap Onderwaater, jacobus.onderwaater@cern.ch
//         2014/Dec/06


//////////////////////////////////////////////////////////////
//
//    Event plane correction options:
//
//
//    Calibration method: 0: Calibration of Q/sqrt(M)
//                        1: Calibration of Q/M
//                        2: Calibration of Q/|Q|
//                        3: Calibration of Q
//
//    Calibration step  : 0: Raw
//                        1: Equalization
//                        2: Recentering
//                        3: Twist
//                        4: Scaling
//
//    Equalization method : 0: M/<M>
//                          1: 1+(M-<M>)/sigma(M)
//
//    Channel list      : Array of channel numbers that are included in Q-vector calculation
//
//    Twist and Scaling method: to be implemented
//                              0: Double harmonic track wise (advisable for TPC)
//                              1: Double harmonic Q wise
//                              2: Correlations
//
//    Event plane detector name : Name to give your event plane (has to be unique)
//
//    Correlation detector names: The detectors used to perform the twist and scaling with correlations
//
//
///////////////////////////////////////////////////////////////

#ifdef __ECLIPSE_IDE

#include "AliAnalysisTaskQnVectorAnalysis.h"

#endif // ifdef __ECLIPSE_IDE declaration and includes for the ECLIPSE IDE

using std::cout;
using std::endl;

#define VAR AliAnalysisTaskQnVectorAnalysis

void DefineHistogramsQnAnalysis(AliQnCorrectionsHistos* histos, TString histClass);


AliAnalysisTask* AddTaskQnVectorAnalysis(Bool_t bUseMultiplicity, Bool_t b2015DataSet) {

  AliAnalysisTaskQnVectorAnalysis* taskQn = new AliAnalysisTaskQnVectorAnalysis("QnAnalysis");

  /* let's establish the event cuts for event selection */
  AliQnCorrectionsCutsSet *eventCuts = new AliQnCorrectionsCutsSet();
  eventCuts->Add(new AliQnCorrectionsCutWithin(VAR::kVtxZ,zvertexMin,zvertexMax));
  if (bUseMultiplicity) {
    varForEventMultiplicity = VAR::kVZEROMultPercentile;
  }
  else {
    varForEventMultiplicity = VAR::kCentVZERO;
  }
  eventCuts->Add(new AliQnCorrectionsCutWithin(varForEventMultiplicity,centralityMin,centralityMax));
  taskQn->SetEventCuts(eventCuts);
  taskQn->SetCentralityVariable(varForEventMultiplicity);

  if (!b2015DataSet) {
    taskQn->SelectCollisionCandidates(AliVEvent::kMB);  // Events passing trigger and physics selection for analysis
  }
  else
    taskQn->SelectCollisionCandidates(AliVEvent::kMB|AliVEvent::kINT7);  // Events passing trigger and physics selection for analysis


  AliQnCorrectionsHistos* hists = taskQn->GetHistograms();
  TString histClass = "";
  histClass += "Event_NoCuts;";
  histClass += "Event_Analysis;";
  DefineHistogramsQnAnalysis(hists, histClass);

  return taskQn;
}




//__________________________________________________________________
void DefineHistogramsQnAnalysis(AliQnCorrectionsHistos* histos, TString histClass) {
  //
  // define the histograms
  //
  //#define VAR AliQnCorrectionsVarManager

  const Char_t* histClasses = histClass.Data();

  cout << "Defining histograms ..." << endl;
  cout << "histogram classes: " << histClass<< endl;

  //fHistosFile=new TFile(output,"RECREATE");

  TString classesStr(histClasses);
  TObjArray* arr=classesStr.Tokenize(";");

  const Int_t kNRunBins = 3000;
  Double_t runHistRange[2] = {137000.,140000.};

  for(Int_t iclass=0; iclass<arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();
    cout << "hist class: " << classStr.Data() << endl;

    // Event wise histograms
    if(classStr.Contains("Event")) {
      histos->AddHistClass(classStr.Data());
      histos->AddHistogram(classStr.Data(),"RunNo","Run numbers;Run", kFALSE, kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo);
      histos->AddHistogram(classStr.Data(),"BC","Bunch crossing;BC", kFALSE,3000,0.,3000.,VAR::kBC);
      histos->AddHistogram(classStr.Data(),"IsPhysicsSelection","Physics selection flag;;", kFALSE,
          2,-0.5,1.5,VAR::kIsPhysicsSelection, 0,0.0,0.0,VAR::kNothing, 0,0.0,0.0,VAR::kNothing, "off;on");

      histos->AddHistogram(classStr.Data(),"VtxZ","Vtx Z;vtx Z (cm)", kFALSE,300,0.0,0.0,VAR::kVtxZ);
      //histos->AddHistogram(classStr.Data(),"VtxZ","Vtx Z;vtx Z (cm)", kFALSE,300,-15.,15.,VAR::kVtxZ);
      histos->AddHistogram(classStr.Data(),"VtxX","Vtx X;vtx X (cm)", kFALSE,300,-1.,1.,VAR::kVtxX);
      histos->AddHistogram(classStr.Data(),"VtxY","Vtx Y;vtx Y (cm)", kFALSE,300,-1.,1.,VAR::kVtxY);


      histos->AddHistogram(classStr.Data(),"CentVZEROvsMultPVZERO","Multiplicity percentile (VZERO);multiplicity VZERO (percents);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kVZEROMultPercentile, 100, 0.0, 100.0, VAR::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentVZEROvsCentSPD","Centrality(VZERO);centrality VZERO (percents);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentVZERO, 100, 0.0, 100.0, VAR::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPCvsCentSPD","Centrality(TPC);centrality TPC (percents);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentTPC, 100, 0.0, 100.0, VAR::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPCvsCentVZERO","Centrality(TPC);centrality TPC (percents);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentTPC, 100, 0.0, 100.0, VAR::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentTPCvsCentZDC","Centrality(TPC);centrality TPC (percents);centrality ZDC (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentTPC, 100, 0.0, 100.0, VAR::kCentZDC);
      histos->AddHistogram(classStr.Data(),"CentZDCvsCentVZERO","Centrality(ZDC);centrality ZDC (percents);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentZDC, 100, 0.0, 100.0, VAR::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentZDCvsCentSPD","Centrality(ZDC);centrality ZDC (percents);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentZDC, 100, 0.0, 100.0, VAR::kCentSPD);

      histos->AddHistogram(classStr.Data(),"MultVZEROvsCentVZERO","Multiplicity;multiplicity VZERO;VZERO centrality", kFALSE,
          100, 0.0, 32000.0, VAR::kVZEROTotalMult, 100,0.,100., VAR::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"MultSPDvsCentPSD","Multiplicity;SPD tracklets;SPD centrality", kFALSE,
          100, 0.0, 3000.0, VAR::kSPDntracklets, 100,0.,100., VAR::kCentSPD);
      histos->AddHistogram(classStr.Data(),"MultTPCvsCentSPD","Multiplicity;TPC selected tracks;TPC centrality", kFALSE,
          100, 0.0, 3500.0, VAR::kNtracksSelected, 100,0.,100., VAR::kCentTPC);
      histos->AddHistogram(classStr.Data(),"MultZDCvsCentZDC","Multiplicity;multiplicity ZDC;ZDC centrality", kFALSE,
          100, 0.0, 300000.0, VAR::kZDCTotalEnergy, 100,0.,100., VAR::kCentZDC);


      histos->AddHistogram(classStr.Data(),"MultTPCvsMultVZERO","Multiplicity;tracks TPC;multiplicity VZERO", kFALSE,
          100, 0.0, 3500.0, VAR::kNtracksSelected, 100, 0.0, 32000.0, VAR::kVZEROTotalMult);
      histos->AddHistogram(classStr.Data(),"MultTPCvsMultSPD","Multiplicity;tracklets SPD;tracks TPC", kFALSE,
          100, 0.0, 3500.0, VAR::kNtracksSelected, 100, 0.0, 3000.0, VAR::kSPDntracklets);
      histos->AddHistogram(classStr.Data(),"MultSPDvsMultVZERO","Multiplicity;tracklets SPD;multiplicity VZERO", kFALSE,
          100, 0.0, 32000.0, VAR::kVZEROTotalMult, 100, 0.0, 3000.0, VAR::kSPDntracklets);
      histos->AddHistogram(classStr.Data(),"MultTPCvsMultZDC","Multiplicity;tracks TPC;energy ZDC", kFALSE,
          100, 0.0, 3500.0, VAR::kNtracksSelected, 100, 0.0, 300000.0, VAR::kZDCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultVZEROvsMultZDC","Multiplicity;multiplicity VZERO;energy ZDC", kFALSE,
          100, 0.0, 32000.0, VAR::kVZEROTotalMult, 100, 0.0, 300000.0, VAR::kZDCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultSPDvsMultZDC","Multiplicity;tracklets SPD;energy ZDC", kFALSE,
          100, 0.0, 3000.0, VAR::kSPDntracklets, 100, 0.0, 300000.0, VAR::kZDCTotalEnergy);



      histos->AddHistogram(classStr.Data(),"MultVZERO","Multiplicity;multiplicity VZERO", kFALSE,
          320, 0.0, 25000.0, VAR::kVZEROTotalMult);
      histos->AddHistogram(classStr.Data(),"MultVZEROA","Multiplicity;multiplicity VZEROA", kFALSE,
          250, 0.0, 0.0, VAR::kVZEROATotalMult);//10000.0
      histos->AddHistogram(classStr.Data(),"MultVZEROC","Multiplicity;multiplicity VZEROC", kFALSE,
          250, 0.0, 0.0, VAR::kVZEROCTotalMult);//15000.0
      histos->AddHistogram(classStr.Data(),"MultZDC","Multiplicity;multiplicity ZDC", kFALSE,
          200, 0.0, 300000.0, VAR::kZDCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultZDCA","Multiplicity;multiplicity ZDCA", kFALSE,
          200, 0.0, 150000.0, VAR::kZDCATotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultZDCC","Multiplicity;multiplicity ZDCC", kFALSE,
          200, 0.0, 150000.0, VAR::kZDCCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultFMD1","Multiplicity;multiplicity FMD1", kFALSE,
          300, 0.0, 10000.0, VAR::kFMD1TotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD2I","Multiplicity;multiplicity FMD2I", kFALSE,
          300, 0.0, 10000.0, VAR::kFMD2ITotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD2O","Multiplicity;multiplicity FMD2O", kFALSE,
          300, 0.0, 10000.0, VAR::kFMD2OTotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD3I","Multiplicity;multiplicity FMD3I", kFALSE,
          300, 0.0, 10000.0, VAR::kFMD3ITotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD3O","Multiplicity;multiplicity FMD3O", kFALSE,
          300, 0.0, 10000.0, VAR::kFMD3OTotalMult);
      histos->AddHistogram(classStr.Data(),"MultTZEROA","Multiplicity;multiplicity TZEROA", kFALSE,
          300, 0.0, 3000.0, VAR::kTZEROATotalMult);
      histos->AddHistogram(classStr.Data(),"MultTZEROC","Multiplicity;multiplicity TZEROC", kFALSE,
          300, 0.0, 3000.0, VAR::kTZEROCTotalMult);




      histos->AddHistogram(classStr.Data(),"MultPercentVZERO","Multiplicity percentile (VZERO);multiplicity VZERO (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kVZEROMultPercentile);
      histos->AddHistogram(classStr.Data(),"CentVZERO","Centrality(VZERO);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentSPD","Centrality(SPD);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPC","Centrality(TPC);centrality TPC (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentTPC);
      histos->AddHistogram(classStr.Data(),"CentZDC","Centrality(ZDC);centrality ZDC (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentZDC);

      histos->AddHistogram(classStr.Data(),"CentQuality","Centrality quality;centrality quality", kFALSE,
          100, -50.5, 49.5, VAR::kCentQuality);
      histos->AddHistogram(classStr.Data(),"CentVZERO_Run_prof","<Centrality(VZERO)> vs run;Run; centrality VZERO (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0.0, 100.0, VAR::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentSPD_Run_prof","<Centrality(SPD)> vs run;Run; centrality SPD (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0.0, 100.0, VAR::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPC_Run_prof","<Centrality(TPC)> vs run;Run; centrality TPC (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0.0, 100.0, VAR::kCentTPC);
      histos->AddHistogram(classStr.Data(),"CentZDC_Run_prof","<Centrality(ZDC)> vs run;Run; centrality ZDC (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0.0, 100.0, VAR::kCentZDC);


      histos->AddHistogram(classStr.Data(),"NV0sTotal","Number of V0 candidates per event;# pairs", kFALSE,
          1000,0.,30000.,VAR::kNV0total);
      histos->AddHistogram(classStr.Data(),"NV0sSelected","Number of selected V0 candidates per event;# pairs", kFALSE,
          1000,0.,10000.,VAR::kNV0selected);
      histos->AddHistogram(classStr.Data(),"NPairs","Number of candidates per event;# pairs", kFALSE,
          5000,0.,5000.,VAR::kNdielectrons);
      histos->AddHistogram(classStr.Data(),"NPairsSelected", "Number of selected pairs per event; #pairs", kFALSE,
          5000,0.,5000.,VAR::kNpairsSelected);
      histos->AddHistogram(classStr.Data(),"NTracksTotal","Number of total tracks per event;# tracks", kFALSE,
          1000,0.,30000.,VAR::kNtracksTotal);
      histos->AddHistogram(classStr.Data(),"NTracksSelected","Number of selected tracks per event;# tracks", kFALSE,
          1000,0.,30000.,VAR::kNtracksSelected);
      histos->AddHistogram(classStr.Data(),"SPDntracklets", "SPD #tracklets; tracklets", kFALSE,
          3000, -0.5, 2999.5, VAR::kSPDntracklets);
      histos->AddHistogram(classStr.Data(),"SPDnSingleClusters", "SPD #single clusters; tracklets", kFALSE,
          3000, -0.5, 2999.5, VAR::kSPDnSingleClusters);

      histos->AddHistogram(classStr.Data(),"NV0total_Run_prof", "<Number of total V0s> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNV0total);
      histos->AddHistogram(classStr.Data(),"NV0selected_Run_prof", "<Number of selected V0s> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNV0selected);
      histos->AddHistogram(classStr.Data(),"Ndielectrons_Run_prof", "<Number of dielectrons> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNdielectrons);
      histos->AddHistogram(classStr.Data(),"NpairsSelected_Run_prof", "<Number of selected pairs> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNpairsSelected);
      histos->AddHistogram(classStr.Data(),"NTracksTotal_Run_prof", "<Number of tracks> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNtracksTotal);
      histos->AddHistogram(classStr.Data(),"NTracksSelected_Run_prof", "<Number of selected tracks> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNtracksSelected);
      histos->AddHistogram(classStr.Data(),"SPDntracklets_Run_prof", "<SPD ntracklets> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kSPDntracklets);

      histos->AddHistogram(classStr.Data(),"VtxZ_CentVZERO","Centrality(VZERO) vs vtx. Z;vtx Z (cm); centrality VZERO (%)", kFALSE,
          300,-15.,15.,VAR::kVtxZ, 100, 0.0, 100.0, VAR::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"VtxZ_CentSPD","Centrality(SPD) vs vtx. Z;vtx Z (cm); centrality SPD (%)", kFALSE,
          300,-15.,15.,VAR::kVtxZ, 100, 0.0, 100.0, VAR::kCentSPD);
      histos->AddHistogram(classStr.Data(),"VtxZ_CentTPC","Centrality(TPC) vs vtx. Z;vtx Z (cm); centrality TPC (%)", kFALSE,
          300,-15.,15.,VAR::kVtxZ, 100, 0.0, 100.0, VAR::kCentTPC);
      continue;
    }  // end if className contains "Event"


    // Offline trigger histograms
    if(classStr.Contains("OfflineTriggers")) {
      histos->AddHistClass(classStr.Data());

      TString triggerNames = "";
      for(Int_t i=0; i<64; ++i) {triggerNames += Form("%s",VAR::fOfflineTriggerNames[i]); triggerNames+=";";}

      histos->AddHistogram(classStr.Data(), "Triggers", "Offline triggers fired; ; ;", kFALSE,
          64, -0.5, 63.5, VAR::kOfflineTrigger, 2, -0.5, 1.5, VAR::kOfflineTriggerFired, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data(), "off;on");
      histos->AddHistogram(classStr.Data(), "Triggers2", "Offline triggers fired; ; ;", kFALSE,
          64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 0, 0.0, 0.0, VAR::kNothing, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentVZERO_Triggers2", "Offline triggers fired vs centrality VZERO; ; centrality VZERO;", kFALSE,
          64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 20, 0.0, 100.0, VAR::kCentVZERO, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentTPC_Triggers2", "Offline triggers fired vs centrality TPC; ; centrality TPC;", kFALSE,
          64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 20, 0.0, 100.0, VAR::kCentTPC, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentSPD_Triggers2", "Offline triggers fired vs centrality SPD; ; centrality SPD;", kFALSE,
          64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 20, 0.0, 100.0, VAR::kCentSPD, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentZDC_Triggers2", "Offline triggers fired vs centrality ZDC; ; centrality ZDC;", kFALSE,
          64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 20, 0.0, 100.0, VAR::kCentZDC, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "VtxZ_Triggers2", "Offline triggers fired vs vtxZ; ; vtx Z (cm.);", kFALSE,
          64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 200, -20.0, 20.0, VAR::kVtxZ, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      continue;
    }

    // Track histograms
    if(classStr.Contains("Tracks")) {
      histos->AddHistClass(classStr.Data());
      for(Int_t ih=0; ih<6; ++ih) {
        histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO",ih+1), Form("<cos%d #varphi> vs (CentVZERO); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
            20, 0.0, 100.0, VAR::kCentVZERO, 500, -1., 1., VAR::kCosNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO",ih+1), Form("<sin%d #varphi> vs (CentVZERO); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
            20, 0.0, 100.0, VAR::kCentVZERO, 500, -1.0, 1.0, VAR::kSinNPhi+ih);
      }
    }


    // Track histograms
    if(classStr.Contains("TrackQA")) {
      histos->AddHistClass(classStr.Data());

      histos->AddHistogram(classStr.Data(), "Pt", "p_{T} distribution; p_{T} (GeV/c^{2});", kFALSE,
          1000, 0.0, 50.0, VAR::kPt);
      histos->AddHistogram(classStr.Data(), "Eta", "#eta illumination; #eta;", kFALSE,
          1000, -1.5, 1.5, VAR::kEta);
      histos->AddHistogram(classStr.Data(), "Phi", "#varphi illumination; #varphi;", kFALSE,
          1000, 0.0, 6.3, VAR::kPhi);
      histos->AddHistogram(classStr.Data(), "DCAxy", "DCAxy; DCAxy (cm.)", kFALSE,
          1000, -10.0, 10.0, VAR::kDcaXY);
      histos->AddHistogram(classStr.Data(), "DCAz", "DCAz; DCAz (cm.)", kFALSE,
          1000, -10.0, 10.0, VAR::kDcaZ);
      histos->AddHistogram(classStr.Data(), "TPCncls", "TPCncls; TPCncls", kFALSE,
          160, 0.0, 160.0, VAR::kTPCncls);

      // run dependence
      histos->AddHistogram(classStr.Data(), "Pt_Run", "<p_{T}> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 1000, 0.0, 50.0, VAR::kPt);
      histos->AddHistogram(classStr.Data(), "Eta_Run", "<#eta> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 1000, -1.5, 1.5, VAR::kEta);
      histos->AddHistogram(classStr.Data(), "Phi_Run", "<#varphi> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 1000, 0.0, 6.3, VAR::kPhi);
      histos->AddHistogram(classStr.Data(), "DCAxy_Run", "<DCAxy> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 1000, -10.0, 10.0, VAR::kDcaXY);
      histos->AddHistogram(classStr.Data(), "DCAz_Run", "<DCAz> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 1000, -10.0, 10.0, VAR::kDcaZ);

      // correlations between parameters
      histos->AddHistogram(classStr.Data(), "Eta_Pt_prof", "<p_{T}> vs #eta; #eta; p_{T} (GeV/c);", kTRUE,
          300, -1.5, +1.5, VAR::kEta, 100, 0.0, 10.0, VAR::kPt);
      histos->AddHistogram(classStr.Data(), "Phi_Pt", "p_{T} vs #varphi; #varphi (rad.); p_{T} (GeV/c)", kFALSE,
          300, -0.01, 6.3, VAR::kPhi, 100, 0.0, 2.2, VAR::kPt);
      histos->AddHistogram(classStr.Data(), "Phi_Pt_prof", "<p_{T}> vs #varphi; #varphi (rad.); p_{T} (GeV/c)", kTRUE,
          300, 0.0, 6.3, VAR::kPhi, 100, 0.0, 10.0, VAR::kPt);
      histos->AddHistogram(classStr.Data(), "Eta_Phi", "#varphi vs #eta; #eta; #varphi (rad.);", kFALSE,
          200, -1.0, +1.0, VAR::kEta, 100, 0.0, 6.3, VAR::kPhi);
      histos->AddHistogram(classStr.Data(), "TPCncls_Eta_Phi_prof", "<TPC ncls> vs #varphi vs #eta; #eta; #varphi (rad.);TPC ncls", kTRUE,
          200, -1.0, +1.0, VAR::kEta, 100, 0.0, 6.3, VAR::kPhi, 10, 0.0, 200., VAR::kTPCncls);
      histos->AddHistogram(classStr.Data(), "DCAxy_Eta_Phi_prof", "<DCAxy> vs #varphi vs #eta; #eta; #varphi (rad.);DCAxy (cm)", kTRUE,
          200, -1.0, +1.0, VAR::kEta, 100, 0.0, 6.3, VAR::kPhi, 10, 0.0, 200., VAR::kDcaXY);
      histos->AddHistogram(classStr.Data(), "DCAz_Eta_Phi_prof", "<DCAz> vs #varphi vs #eta; #eta; #varphi (rad.);DCAz (cm)", kTRUE,
          200, -1.0, +1.0, VAR::kEta, 100, 0.0, 6.3, VAR::kPhi, 10, 0.0, 200., VAR::kDcaZ);
      histos->AddHistogram(classStr.Data(), "Pt_DCAxy", "DCAxy vs p_{T}; p_{T} (GeV/c); DCA_{xy} (cm)", kFALSE,
          100, 0.0, 10.0, VAR::kPt, 500, -2.0, 2.0, VAR::kDcaXY);
      histos->AddHistogram(classStr.Data(), "Pt_DCAz", "DCAz vs p_{T}; p_{T} (GeV/c); DCA_{z} (cm)", kFALSE,
          100, 0.0, 10.0, VAR::kPt, 500, -2.0, 2.0, VAR::kDcaZ);
      histos->AddHistogram(classStr.Data(), "Eta_DCAxy", "DCAxy vs #eta; #eta; DCA_{xy} (cm)", kFALSE,
          100, -1.0, 1.0, VAR::kEta, 500, -2.0, 2.0, VAR::kDcaXY);
      histos->AddHistogram(classStr.Data(), "Eta_DCAz", "DCAz vs #eta; #eta; DCA_{z} (cm)", kFALSE,
          100, -1.0, 1.0, VAR::kEta, 500, -2.0, 2.0, VAR::kDcaZ);

      for(Int_t ih=0; ih<6; ++ih) {
        //histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO",ih+1), Form("<cos%d #varphi> vs (CentVZERO); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
        //                   20, 0.0, 100.0, VAR::kCentVZERO, 500, -1., 1., VAR::kCosNPhi+ih);
        //histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO",ih+1), Form("<sin%d #varphi> vs (CentVZERO); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
        //                   20, 0.0, 100.0, VAR::kCentVZERO, 500, -1.0, 1.0, VAR::kSinNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_Pt_Eta",ih+1), Form("<cos%d #varphi> vs (#eta,p_{T}); #eta; p_{T} (GeV/c); <cos%d #varphi>", ih+1, ih+1), kTRUE,
            20, -1.0, 1.0, VAR::kEta, 30, 0.0, 3.0, VAR::kPt, 500, -1.0, 1.0, VAR::kCosNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_Pt_Eta",ih+1), Form("<sin%d #varphi> vs (#eta,p_{T}); #eta; p_{T} (GeV/c); <sin%d #varphi>", ih+1, ih+1), kTRUE,
            20, -1.0, 1.0, VAR::kEta, 30, 0.0, 3.0, VAR::kPt, 500, -1.0, 1.0, VAR::kSinNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO_VtxZ",ih+1), Form("<cos%d #varphi> vs (CentVZERO,VtxZ); Z (cm.); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
            30, -15.0, 15.0, VAR::kVtxZ, 20, 0.0, 100.0, VAR::kCentVZERO, 500, -1., 1., VAR::kCosNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO_VtxZ",ih+1), Form("<sin%d #varphi> vs (CentVZERO,VtxZ); Z (cm.); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
            30, -15.0, 15.0, VAR::kVtxZ, 20, 0.0, 100.0, VAR::kCentVZERO, 500, -1.0, 1.0, VAR::kSinNPhi+ih);
      }
    }

    // Tracklet histograms
    if(classStr.Contains("TrackletQA")) {
      histos->AddHistClass(classStr.Data());

      histos->AddHistogram(classStr.Data(), "Eta", "#eta illumination; #eta;", kFALSE,
          1000, -3.0, 3.0, VAR::kSPDtrackletEta);
      histos->AddHistogram(classStr.Data(), "Phi", "#varphi illumination; #varphi;", kFALSE,
          300, -0.01, 6.3, VAR::kSPDtrackletPhi);
      histos->AddHistogram(classStr.Data(), "Eta_Phi", "#varphi vs #eta; #eta; #varphi (rad.);", kFALSE,
          200, -3.0, +3.0, VAR::kSPDtrackletEta, 100, 0.0, 6.3, VAR::kSPDtrackletPhi);
    }

  }

  cout << " done" << endl;
}


