AliAnalysisTask *AddTaskEMCALCorrectionForEHTask(UInt_t kPhysSel = AliVEvent::kAny)
{
//get the current analysis manager
AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
if (!mgr) {
    Error("AddTaskHFEElecHadronCorrlPbPb", "No analysis manager found.");
    return 0;
}

if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHFEElecHadronCorrlPbPb", "This task requires an input event handler");
    return NULL;
}
TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
/*  if (type=="AOD"){
 ::Error("AddTaskHFEElecHadronCorrlPbPb", "The tasks exits because AODs are in input");
 return NULL;
 }
 */
Bool_t MCthere=kTRUE;
AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
if(!mcH){
    MCthere=kFALSE;
}

    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEMCALTender.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskClusterizerFast.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusterMaker.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalAodTrackFilter.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C");
    
    TString sTracksName("AODFilterTracks");
    TString sClusName("EmcCaloClusters");
    
    // Setup task
    AliTaskCDBconnect *taskCDB = AddTaskCDBconnect();
    taskCDB->SetFallBackToRaw(kTRUE);
    
    // Tender Supplies
    const char *cPass        = 0;
    Bool_t   bDistBC         = kFALSE; //switch for recalculation cluster position from bad channel
    Bool_t   bRecalibClus    = kFALSE;
    Bool_t   bRecalcClusPos  = kFALSE;
    Bool_t   bNonLinearCorr  = kFALSE;
    Bool_t   bRemExoticCell  = kFALSE;
    Bool_t   bRemExoticClus  = kFALSE;
    Bool_t   bFidRegion      = kFALSE;
    Bool_t   bCalibEnergy    = kTRUE;
    Bool_t   bCalibTime      = kTRUE;
    Bool_t   bRemBC          = kTRUE;
    UInt_t   iNonLinFunct    = AliEMCALRecoUtils::kNoCorrection;
    Bool_t   bReclusterize   = kFALSE;
    Float_t  fSeedThresh     = 0.1;      // 100 MeV
    Float_t  fCellThresh     = 0.05;     // 50 MeV
    UInt_t   iClusterizer    = AliEMCALRecParam::kClusterizerv2;
    Bool_t   bTrackMatch     = kFALSE;
    Bool_t   bUpdateCellOnly = kTRUE;
    Float_t  fEMCtimeMin     = -50e-6;
    Float_t  fEMCtimeMax     =  50e-6;
    Float_t  fEMCtimeCut     =  1e6;
    Bool_t remExoticCell  = kFALSE;
    AliAnalysisTaskSE *pTenderTask = AddTaskEMCALTender(bDistBC, bRecalibClus, bRecalcClusPos, bNonLinearCorr, bRemExoticCell, bRemExoticClus,
                                                        bFidRegion, bCalibEnergy, bCalibTime, bRemBC, iNonLinFunct, bReclusterize, fSeedThresh,
                                                        fCellThresh, iClusterizer, bTrackMatch, bUpdateCellOnly, fEMCtimeMin, fEMCtimeMax, fEMCtimeCut, cPass);
    pTenderTask->SelectCollisionCandidates(kPhysSel);

    //Clusterizer task
    AliAnalysisTaskEMCALClusterizeFast *pClusterizerTask = AddTaskClusterizerFast("ClusterizerFast", "", "", iClusterizer,
                                                                                  fCellThresh, fSeedThresh, fEMCtimeMin, fEMCtimeMax, fEMCtimeCut,
                                                                                  remExoticCell, bDistBC, AliAnalysisTaskEMCALClusterizeFast::kFEEData);
    pClusterizerTask->SelectCollisionCandidates(kPhysSel);
    
    //Cluster maker task
    UInt_t nonLinFunct = AliEMCALRecoUtils::kBeamTestCorrected;
    Bool_t remExoticClus  = kTRUE;
    AliEmcalClusterMaker *pClusterMakerTask = AddTaskEmcalClusterMaker(nonLinFunct, remExoticClus, "usedefault", sClusName, 0., kTRUE);
    pClusterMakerTask->GetClusterContainer(0)->SetClusPtCut(0.);
    pClusterMakerTask->GetClusterContainer(0)->SetClusECut(0.);
    pClusterMakerTask->SelectCollisionCandidates(kPhysSel);
    
    //tracks maker for AOD
    char* periodstr = "16:1:includeNoITS=kTRUE doProp=kTRUE doAttemptProp=kFALSE isMC=kFALSE";
    TString period(periodstr);
    Bool_t   doAODTrackProp     = kTRUE;
    Double_t trackeff           = 1.0;
    Double_t edist = 440;
    AliEmcalAodTrackFilterTask *pHybTask = AddTaskEmcalAodTrackFilter(sTracksName,"tracks",period);
    if (doAODTrackProp) {
        pHybTask->SetDist(edist);
        pHybTask->SetAttemptPropMatch(kTRUE);
    }
    pHybTask->SelectCollisionCandidates(AliVEvent::kAny);
    pHybTask->SetTrackEfficiency(trackeff);
    
    // Cluster-track matcher task
    Double_t maxMatchR = 0.1;
    Bool_t attachEmcalPart = kFALSE;
    Bool_t updateClusters = kTRUE;
    Bool_t updateTracks = kTRUE;
    Bool_t doHistos = kTRUE;
    Int_t  nCentBins = 4;
    AliEmcalClusTrackMatcherTask *pMatcherTask = AddTaskEmcalClusTrackMatcher(sTracksName, sClusName, maxMatchR, attachEmcalPart, updateClusters, updateTracks, doHistos);
    pMatcherTask->SelectCollisionCandidates(kPhysSel);
    pMatcherTask->SetNCentBins(nCentBins);
    pMatcherTask->GetClusterContainer(0)->SetClusECut(0.0);
    pMatcherTask->GetClusterContainer(0)->SetClusPtCut(0.);
    pMatcherTask->GetParticleContainer(0)->SetParticlePtCut(0.0);

    return pMatcherTask;
}
