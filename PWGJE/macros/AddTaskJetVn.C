//==============================================================================
// Add task macro for jet vn task
// Will Witt 2018
//==============================================================================

/*AliAnalysisTaskJetVn* AddTaskJetVn(
    const char* nTracks                        = "usedefault",
    const char* nClusters                      = "usedefault",
    const char *nrho                           = "Rho",
    const AliJetContainer::EJetAlgo_t jetAlgo  = AliJetContainer::antikt_algorithm,
    const AliJetContainer::EJetType_t jetType  = AliJetContainer::kChargedJet, // might need to change this to kChargedJet
    const Double_t minTrPt                     = 0.15,
    const Double_t minClPt                     = 0.30,
    const Double_t ghostArea                   = 0.005,
    const AliJetContainer::ERecoScheme_t reco  = AliJetContainer::pt_scheme,
    const char* tag                            = "Jet",
    const Double_t minJetPt                    = 0.,
    Double_t   jetRadius                       = 0.2,
    Double_t   jetptcut                        = 1,
    Double_t   jetareacut                      = 0.557, // Area > 0.6*pi*R^2
    const Bool_t lockTask                      = kTRUE,
    const Bool_t bFillGhosts                   = kFALSE,
    //adding new variables from Redmer's jet v3 code
    //const char* type              = "TPC",
    //Int_t      leadhadtype        = 0,
    //const char *taskname           = "AliAnalysisTaskJetVn",
    //UInt_t     runMode            = AliAnalysisTaskJetV3::kLocal, // was kGrid
    //Bool_t     fillQA             = kTRUE,
    //TString    fitOpts            = "WLQI",
    //UInt_t     fitType            = AliAnalysisTaskJetV3::kCombined,
    //Double_t   trackptcut         = .15,
    //Bool_t     LHC10h             = kFALSE,
    Bool_t     LHC15o             = kTRUE,
    //Bool_t     addEPweights       = kFALSE,
    //Bool_t     baseClassHistos    = kTRUE,
    Float_t    minEta             = -.7,
    Float_t    maxEta             = .7,
    //UInt_t     acceptance         = AliEmcalJet::kUser,
    AliJetContainer::EJetType_t jetType = AliJetContainer::kChargedJet,
    AliJetContainer::ERecoScheme_t rscheme = AliJetContainer::pt_scheme)
)
{
    //==========================================================================
    // Get the pointer to the existing analysis manager via the static access method.
    //==========================================================================

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AddTaskJetVn", "No analysis manager to connect to.");
        return NULL;
    }

    //==========================================================================
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==========================================================================

    AliVEventHandler* handler = mgr->GetInputEventHandler();
    if (!handler) {
        ::Error("AddTaskJetVn", "This task requires an input event handler");
        return NULL;
    }

    //==========================================================================
    // Data types, names, and containers
    //==========================================================================

    enum EDataType_t {kUnknown, kESD, kAOD, kMCgen};

    EDataType_t dataType = kUnknown;
    if (handler->InheritsFrom("AliESDInputHandler")) dataType = kESD;
    else if (handler->InheritsFrom("AliAODInputHandler")) dataType = kAOD;
    else if (handler->InheritsFrom("AliMCGenHandler")) dataType = kMCgen;

    TString trackName(nTracks);
    TString clusName(nClusters);

    if (trackName == "usedefault") {
        if (dataType == kESD) trackName = "Tracks";
        else if (dataType == kAOD) trackName = "tracks";
        else if (dataType == kMCgen) trackName = "mcparticles";
    }
    if (clusName == "usedefault") {
        if (dataType == kESD) clusName = "CaloClusters";
        //else if (dataType == kAOD) clusName = "caloClusters";
        else if (dataType == kAOD) clusName = "EmcCaloClusters";
        else clusName = "";
    }

    AliParticleContainer* partCont = 0;
    if (trackName == "mcparticles") {
        AliMCParticleContainer* mcpartCont = new AliMCParticleContainer(trackName);
        partCont = mcpartCont;
    }
    else if (trackName == "tracks" || trackName == "Tracks") {
        AliTrackContainer* trackCont = new AliTrackContainer(trackName);
        partCont = trackCont;
    }
    else if (!trackName.IsNull()) {
        partCont = new AliParticleContainer(trackName);
    }
    if (partCont) partCont->SetParticlePtCut(minTrPt);

    AliClusterContainer* clusCont = 0;
    if (!clusName.IsNull()) {
        clusCont = new AliClusterContainer(clusName);
        clusCont->SetClusECut(0.);
        clusCont->SetClusPtCut(0.);
        clusCont->SetClusHadCorrEnergyCut(minClPt);
        clusCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
    }

    switch (jetType) {
    case AliJetContainer::kChargedJet:
        if (partCont) partCont->SetCharge(AliParticleContainer::kCharged);
        break;
    case AliJetContainer::kNeutralJet:
        if (partCont) partCont->SetCharge(AliParticleContainer::kNeutral);
        break;
    default:
    break;
    }

    TString name = AliJetContainer::GenerateJetName(jetType, jetAlgo, reco, jetRadius, partCont, clusCont, tag);

    Printf("Jet task name: %s", name.Data());

    AliAnalysisTaskJetVn* mgrTask = static_cast<AliAnalysisTaskJetVn *>(mgr->GetTask(name.Data()));
    if (mgrTask) return mgrTask;

    //==========================================================================
    // Init the task and do settings
    //==========================================================================

    AliAnalysisTaskJetVn* jetTask = new AliAnalysisTaskJetVn(name); // do I need more arguments here?

    jetTask->SetJetType(jetType);
    jetTask->SetJetAlgo(jetAlgo);
    jetTask->SetRecombScheme(reco);
    jetTask->SetRadius(jetRadius);
    if (partCont) jetTask->AdoptParticleContainer(partCont);
    if (clusCont) jetTask->AdoptClusterContainer(clusCont);
    jetTask->SetJetsName(tag);
    jetTask->SetMinJetPt(minJetPt);
    jetTask->SetGhostArea(ghostArea);
    jetTask->SetVzRange(-10., 10.);

    if (bFillGhosts) jetTask->SetFillGhost();
    if (lockTask) jetTask->SetLocked();

    if (LHC10h) jetTask->SetCollisionType(AliAnalysisTaskJetVn::kPbPb10h);
    else if (LHC15o) jetTask->SetCollisionType(AliAnalysisTaskJetVn::kPbPb15o);

    // a standard centrality binning is necessary (which corresponds to the binning of the ep weights)
    // so we overwrite possible existing binnings
    //Double_t cw[] = {0., 2., 4., 6., 8., 10., 30., 50., 90.};
    //jetTask->SetCentralityClasses(new TArrayD(sizeof(cw)/sizeof(cw[0]), cw));

    //==========================================================================
    // Final settings, pass to manager and set the containers
    //==========================================================================

    mgr->AddTask(jetTask);

    // Create containers for input/output
    AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
    mgr->ConnectInput(jetTask, 0, cinput);

    TObjArray* cnt = mgr->GetContainers();

    return jetTask;
}*/

AliAnalysisTaskJetVn* AddTaskJetVn(
  const char *nTracks            = "usedefault",
  const char *nClusters          = "usedefault",
  const char* nCells             = "usedefault",
  const char *suffix             = ""
)
{
  return AliAnalysisTaskJetVn::AddTaskJetVn(
      nTracks,
      nClusters,
      nCells,
      suffix);
}
