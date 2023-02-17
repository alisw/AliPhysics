//_____________________________________________________________________
AliAnalysisTask *AddTaskJCDijetTask(TString taskName,
                                    Bool_t isMC,
                                    TString sJCatalyst        = "JCatalystTask",
                                    TString sJCatalystDetMC   = "JCatalystDetMCTask",
                                    UInt_t flags              = 0,
                                    TString centBins          = "0.0 5.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0",
                                    TString sDijetMBins       = "0, 20, 40, 45, 55, 65, 75, 85, 100, 120, 150, 250, 400, 500, 100000",
                                    double jetCone            = 0.4,
                                    double ktjetCone          = 0.4,
                                    int ktScheme              = 1,
                                    int antiktScheme          = 1,
                                    Bool_t usePionMass        = false,
                                    Bool_t useDeltaPhiBGSubtr = true,
                                    double particleEtaCut     = 0.8,
                                    double particlePtCut      = 0.15,
                                    double leadingJetCut      = 20.0,
                                    double subleadingJetCut   = 20.0,
                                    double minJetPt           = 10.0,
                                    double constituentCut     = 5.0,
                                    double deltaPhiCut        = 2.0,
                                    double matchingR          = 0.2,
                                    double trackingIneff      = 0.0,    //If any negative number, use pt dependend tracking inefficiency.
                                    TString sAnchorPeriodTracking = "", //Only needed if trackingIneff<0.0
                                    AliJCDijetAna::jetClasses lUnfJetClassTrue = AliJCDijetAna::iAcc,
                                    AliJCDijetAna::jetClasses lUnfJetClassDet = AliJCDijetAna::iAcc,
                                    Bool_t useCoveredAreaRho  = false){

    // Load Custom Configuration and parameters
    // override values with parameters

    // flags can manipulate event selection:
    // 0: no additional events rejected.
    // AliAnalysisTask::DIJET_VERTEX13PA: use IsVertexSelected2013pA
    // AliAnalysisTask::DIJET_PILEUPSPD: use IsPileupFromSPD(3,0.6,3,2,5)
    // AliAnalysisTask::DIJET_UTILSPILEUPSPD: use IsPileUpSPD(InputEvent())
    // Combinations of these can be used by giving argument for example:
    // AliAnalysisTask::DIJET_VERTEX13PA|AliAnalysisTask::DIJET_PILEUPSPD

    // jet recombination schemes can be set with ktScheme argument:
    // E_scheme     = 0
    // pt_scheme    = 1
    // pt2_scheme   = 2
    // Et_scheme    = 3
    // Et2_scheme   = 4
    // BIpt_scheme  = 5
    // BIpt2_scheme = 6

    return AliJCDijetTask::AddTaskJCDijetTask(taskName             ,
                                              isMC                 ,
                                              sJCatalyst           ,
                                              sJCatalystDetMC      ,
                                              flags                ,
                                              centBins             ,
                                              sDijetMBins          ,
                                              jetCone              ,
                                              ktjetCone            ,
                                              ktScheme             ,
                                              antiktScheme         ,
                                              usePionMass          ,
                                              useDeltaPhiBGSubtr   ,
                                              particleEtaCut       ,
                                              particlePtCut        ,
                                              leadingJetCut        ,
                                              subleadingJetCut     ,
                                              minJetPt             ,
                                              constituentCut       ,
                                              deltaPhiCut          ,
                                              matchingR            ,
                                              trackingIneff        ,
                                              sAnchorPeriodTracking,
                                              lUnfJetClassTrue     ,
                                              lUnfJetClassDet      ,
                                              useCoveredAreaRho    );
}

