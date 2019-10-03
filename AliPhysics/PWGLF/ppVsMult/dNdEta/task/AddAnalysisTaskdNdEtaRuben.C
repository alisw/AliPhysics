const Double_t centBinsMultV0M[] = {0., 0.01, 0.1, 1., 5., 10., 15., 20., 30., 40., 50., 70., 100};
const Double_t centBinsMultRef[] = {0., 1., 4., 7., 10., 15., 20., 25., 30., 40., 50., 60., 70., 100., 200.};
const Double_t centBinsMB[] = {0., 1.};

UInt_t trigSel = AliVEvent::kMB;

//__________________________________________________________________

AliAnalysisTaskdNdEtaRuben *
AddAnalysisTaskdNdEtaRuben(const Char_t *outfilename = "AnalysisResults.root",
                           const Char_t *listname = "clist",
                           //
                           Float_t etaMin     =-1,          // min eta range to fill in histos
                           Float_t etaMax     = 1,          // max eta range to fill in histos
                           Float_t zMin       = -7,         // process events with Z vertex min
                           Float_t zMax       =  7,         //                     max positions
                           const char* useCentVar = "V0A",          // centrality variable to use
                           //
                           Float_t cutSigNStd  = 1.5,        // cut on weighed distance used to extract signal
                           Float_t cutSigDPhiS = -1,        // cut on dPhi-phiBent used to extract signal (if negative -> dphi*sqrt(cutSigNStd)
                           Bool_t  useMC  = kFALSE,          // fill MC info
                           //
                           Bool_t doRec  = kFALSE,          // fill data histos from new reco
                           Bool_t doInj  = kFALSE,          // create Inj. bg
                           Bool_t doRot  = kFALSE,          // create Rot. bg
                           //
                           // specific parameters for reconstruction
                           float  phiRot      = 3.14159e+00, // angle for bg. generation with rotation
                           float  injScale    = 1.,//0.7,    // inject injScale*Ncl(Lr1/Lr2) hits
                           Bool_t scaleDTheta = kTRUE,       // scale dTheta by 1/sin^2(theta) in trackleting
                           float  nStdDev     = 25.,         // number of st.dev. for tracklet cut to keep
                           float  dphi        = 0.06,        // dphi window (sigma of tracklet cut)
                           float  dtht        = 0.025,       // dtheta .... (if negative, abs will be used with additional cut on |dthetaX|, apart from w.distance
                           float  phishift    = 0.0045,      // bending shift
                           Bool_t remOvl      = kTRUE,
                           float  ovlPhiCut   = 0.005,
                           float  ovlZetaCut  = 0.05,
                           Bool_t checkReconstructables = kFALSE, //kTRUE, // fill histos for reconstructable (needs useMC and doRec)
//
)
{
    
    if (cutSigDPhiS<0) cutSigDPhiS = TMath::Sqrt(cutSigNStd)*dphi;
    
    /* init analysis name */
    TString analysisName = "dNdEtaRuben";
    
    /* check analysis manager */
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Error("", "cannot get analysis manager");
        return NULL;
    }
    
    /* check input event handler */
    if (!mgr->GetInputEventHandler()) {
        Error("", "cannot get input event handler");
        return NULL;
    }
    
    /* get common input data container */
    AliAnalysisDataContainer *inputc = mgr->GetCommonInputContainer();
    if (!inputc) {
        Error("", "cannot get common input container");
        return NULL;
    }
    
    /* create output data container */
    AliAnalysisDataContainer *outputc1 = mgr->CreateContainer(listname, TList::Class(), AliAnalysisManager::kOutputContainer, outfilename);

    if (!outputc1) {
        Error("", "cannot create output container \"clist\"");
        return NULL;
    }
    
    /*  create task and connect input/output */
    AliAnalysisTaskdNdEtaRuben *task = new AliAnalysisTaskdNdEtaRuben("AliAnalysisTaskdNdEtaRuben");
    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, inputc);
    mgr->ConnectOutput(task, 1, outputc1);
    
    //__________________________________________________________________
    //__________________________________________________________________
    //__________________________________________________________________
    
    /* configure task */
    
    //__________________________________________________________________
    //__________________________________________________________________
    //__________________________________________________________________
    
    //
    task->SetUseCentralityVar(useCentVar);
    
    Double_t *centBins = NULL;
    Int_t nCentBins = 0;
    switch (useCentVar) {
    case "MB": // MB
      centBins = centBinsMB;
      nCentBins = sizeof(centBinsMB) / 8 - 1;
      break;
    case "V0M": // V0M
    case "V0av": // V0av
      centBins = centBinsMultV0M;
    nCentBins = sizeof(centBinsMultV0M) / 8 - 1;
    break;
    case "Ref": // Ref
      centBins = centBinsMultRef;
      nCentBins = sizeof(centBinsMultRef) / 8 - 1;
      break;
    }
    
    task->SetCentPercentiles(centBins, nCentBins);
    task->SetTriggerSelection(trigSel);
    //
    task->SetDoNormalReco(doRec);
    task->SetDoInjection(doInj);
    task->SetDoRotation(doRot);
    //
    task->SetUseMC(useMC);
    task->SetCheckReconstructables(checkReconstructables);
    //
    task->SetEtaMin(etaMin);
    task->SetEtaMax(etaMax);
    task->SetZVertexMin(zMin);
    task->SetZVertexMax(zMax);
    //
    task->SetDPhiSCut(cutSigDPhiS);
    task->SetNStdCut(cutSigNStd);
    //
    task->SetScaleDThetaBySin2T(scaleDTheta);
    task->SetNStdDev(nStdDev);
    task->SetPhiWindow(dphi);
    task->SetThetaWindow(dtht);
    task->SetPhiShift(phishift);
    task->SetPhiOverlapCut(ovlPhiCut);
    task->SetZetaOverlapCut(ovlZetaCut);
    task->SetPhiRot(phiRot);
    task->SetInjScale(injScale);
    task->SetRemoveOverlaps(remOvl);
    //
    //    task->Dump();
    return task;
}

