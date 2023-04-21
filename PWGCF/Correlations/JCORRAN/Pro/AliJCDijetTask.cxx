/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
//______________________________________________________________________________
// Analysis task for providing various dijet informations
// author: O. Saarimaki, D.J. Kim(dong.jo.kim@cern.ch)
// ALICE Group University of Jyvaskyla
// Finland
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////
#include <AliAnalysisManager.h>
#include "AliJCDijetTask.h" 
//#include "AliProdInfo.h"

//______________________________________________________________________________
AliJCDijetTask::AliJCDijetTask() :
    AliAnalysisTaskSE("JCDijetBaseTask"),
    fJCatalystTask(NULL),
    fJCatalystDetMCTask(NULL),
    fJCatalystTaskName("JCatalystTask"),
    fJCatalystDetMCTaskName("JCatalystDetMCTask"),
    fcentralityBins(0),
    fjetCone(0),
    fktJetCone(0),
    fktScheme(0),
    fantiktScheme(0),
    fusePionMass(0),
    fuseDeltaPhiBGSubtr(0),
    fIsMC(kTRUE),
    fparticleEtaCut(0),
    fleadingJetCut(0),
    fsubleadingJetCut(0),
    fMinJetPt(0),
    fconstituentCut(0),
    fdeltaPhiCut(0),
    fmatchingR(0),
    fhistos(NULL),
    fhistosDetMC(NULL),
    fTrackEfficiencyHistogram(NULL),
    fana(NULL),
    fanaMC(NULL),
    fCBin(-1),
    fCBinDetMC(-1),
    fOutput(NULL),
    flags(0),
    fUtils(nullptr)
{
}

//______________________________________________________________________________
AliJCDijetTask::AliJCDijetTask(const char *name, TString inputformat):
    AliAnalysisTaskSE(name), 
    fJCatalystTask(NULL),
    fJCatalystDetMCTask(NULL),
    fJCatalystTaskName("JCatalystTask"),
    fJCatalystDetMCTaskName("JCatalystDetMCTask"),
    fcentralityBins(0),
    fjetCone(0),
    fktJetCone(0),
    fktScheme(0),
    fantiktScheme(0),
    fusePionMass(0),
    fuseDeltaPhiBGSubtr(0),
    fIsMC(kTRUE),
    fparticleEtaCut(0),
    fleadingJetCut(0),
    fsubleadingJetCut(0),
    fMinJetPt(0),
    fconstituentCut(0),
    fdeltaPhiCut(0),
    fmatchingR(0),
    fhistos(NULL),
    fhistosDetMC(NULL),
    fTrackEfficiencyHistogram(NULL),
    fana(NULL),
    fanaMC(NULL),
    fCBin(-1),
    fCBinDetMC(-1),
    fOutput(NULL),
    flags(0),
    fUtils(nullptr)
{
    // Constructor
    AliInfo(Form("---- AliJCDijetTask Constructor: %s ----",name));
    DefineOutput (1, TDirectory::Class());
}

//____________________________________________________________________________
AliJCDijetTask::AliJCDijetTask(const AliJCDijetTask& ap) :
    AliAnalysisTaskSE(ap.GetName()), 
    fJCatalystTask(ap.fJCatalystTask),
    fJCatalystDetMCTask(ap.fJCatalystDetMCTask),
    fJCatalystTaskName(ap.fJCatalystTaskName),
    fJCatalystDetMCTaskName(ap.fJCatalystDetMCTaskName),
    fcentralityBins(ap.fcentralityBins),
    fjetCone(ap.fjetCone),
    fktJetCone(ap.fktJetCone),
    fktScheme(ap.fktScheme),
    fantiktScheme(ap.fantiktScheme),
    fusePionMass(ap.fusePionMass),
    fuseDeltaPhiBGSubtr(ap.fuseDeltaPhiBGSubtr),
    fIsMC(ap.fIsMC),
    fparticleEtaCut(ap.fparticleEtaCut),
    fleadingJetCut(ap.fleadingJetCut),
    fsubleadingJetCut(ap.fsubleadingJetCut),
    fMinJetPt(ap.fMinJetPt),
    fconstituentCut(ap.fconstituentCut),
    fdeltaPhiCut(ap.fdeltaPhiCut),
    fmatchingR(ap.fmatchingR),
    fhistos(ap.fhistos),
    fhistosDetMC(ap.fhistosDetMC),
    fTrackEfficiencyHistogram(ap.fTrackEfficiencyHistogram),
    fana(ap.fana),
    fanaMC(ap.fanaMC),
    fCBin(ap.fCBin),
    fCBinDetMC(ap.fCBinDetMC),
    fOutput( ap.fOutput ),
    flags(ap.flags),
    fUtils(ap.fUtils)
{ 

    AliInfo("----DEBUG AliJCDijetTask COPY ----");

}

//_____________________________________________________________________________
AliJCDijetTask& AliJCDijetTask::operator = (const AliJCDijetTask& ap)
{
    // assignment operator

    AliInfo("----DEBUG AliJCDijetTask operator= ----");
    this->~AliJCDijetTask();
    new(this) AliJCDijetTask(ap);
    return *this;
}

//______________________________________________________________________________
AliJCDijetTask::~AliJCDijetTask()
{
    // destructor 
    delete fOutput;
    delete fhistos;
    delete fhistosDetMC;
    delete fTrackEfficiencyHistogram;
    delete fana;
    delete fanaMC;
    delete fUtils;

}

//________________________________________________________________________

void AliJCDijetTask::UserCreateOutputObjects()
{  
    //=== create the jcorran outputs objects
    if(fDebug > 1) printf("AliJCDijetTask::UserCreateOutPutObjects() \n");
    //=== Get AnalysisManager
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();

    fJCatalystTask = (AliJCatalystTask*)(man->GetTask( fJCatalystTaskName ));
    if(fIsMC) fJCatalystDetMCTask = (AliJCatalystTask*)(man->GetTask( fJCatalystDetMCTaskName ));

    OpenFile(1);
    fOutput = gDirectory;
    fOutput->cd();

    fEventCuts.AddQAplotsToList(fJCatalystTask->GetCataList());

    fana = new AliJCDijetAna();
    fhistos = new AliJCDijetHistos();
    fhistos->SetName("jcdijet");
    fhistos->SetCentralityBinsHistos(fcentralityBins);
    fhistos->SetDijetMBinsHistos(fsDijetMBins);
    fhistos->SetNJetClasses(fana->jetClassesSize);
    fhistos->CreateEventTrackHistos();
    fhistos->fHMG->Print();

    if(fIsMC) {
        fanaMC = new AliJCDijetAna();
        fhistosDetMC = new AliJCDijetHistos();
        fhistosDetMC->SetName("jcdijetDetMC");
        fhistosDetMC->SetCentralityBinsHistos(fcentralityBins);
        fhistosDetMC->SetDijetMBinsHistos(fsDijetMBins);
        fhistosDetMC->SetNJetClasses(fanaMC->jetClassesSize);
        fhistosDetMC->CreateEventTrackHistos();
        fhistosDetMC->fHMG->Print();
    }

    fUtils = new AliAnalysisUtils();
    fUtils->SetMaxVtxZ(10);

    TString sktScheme;
    TString santiktScheme;
    switch (fktScheme) {
        case 0:  sktScheme = "E_scheme";
                 break;
        case 1:  sktScheme = "pt_scheme";
                 break;
        case 2:  sktScheme = "pt2_scheme";
                 break;
        case 3:  sktScheme = "Et_scheme";
                 break;
        case 4:  sktScheme = "Et2_scheme";
                 break;
        case 5:  sktScheme = "BIpt_scheme";
                 break;
        case 6:  sktScheme = "BIpt2_scheme";
                 break;
        default: sktScheme = "Unknown, check macro arguments!";
                 break;
    }
    switch (fantiktScheme) {
        case 0:  santiktScheme = "E_scheme";
                 break;
        case 1:  santiktScheme = "pt_scheme";
                 break;
        case 2:  santiktScheme = "pt2_scheme";
                 break;
        case 3:  santiktScheme = "Et_scheme";
                 break;
        case 4:  santiktScheme = "Et2_scheme";
                 break;
        case 5:  santiktScheme = "BIpt_scheme";
                 break;
        case 6:  santiktScheme = "BIpt2_scheme";
                 break;
        default: santiktScheme = "Unknown, check macro arguments!";
                 break;
    }

    cout << endl;
    cout << "===========SETTINGS===========" << endl;
    cout << "MC:                         " << IsMC() << endl;
    cout << "Centrality bins:            ";
    for(unsigned i=0; i< fcentralityBins.size(); i++) cout << fcentralityBins.at(i) << " ";
    cout << endl;
    cout << "Dijet M bins:               " << fsDijetMBins.Data() << endl;
    cout << "Jet cone size:              " << fjetCone << endl;
    cout << "kt-jet cone size:           " << fktJetCone << endl;
    cout << "Using kt-jet scheme:        " << sktScheme.Data() << endl;
    cout << "Using antikt-jet scheme:    " << santiktScheme.Data() << endl;
    cout << "Using pion mass:            " << fusePionMass << endl;
    cout << "Using DeltaPhi in BG subtr: " << fuseDeltaPhiBGSubtr << endl;
    cout << "Particle eta cut:           " << fparticleEtaCut << endl;
    cout << "Particle pt min cut:        " << fparticlePtCut << endl;
    cout << "Particle pt max cut:        " << 100.0 << endl;
    cout << "Dijet leading jet cut:      " << fleadingJetCut << endl;
    cout << "Dijet subleading jet cut:   " << fsubleadingJetCut << endl;
    cout << "Jet min pt cut:             " << fMinJetPt << endl;
    cout << "Jet leading const. cut:     " << fconstituentCut << endl;
    cout << "Dijet DeltaPhi cut:         pi/" << fdeltaPhiCut << endl;
    cout << "Matching R for MC:          " << fmatchingR << endl;
    cout << "Tracking ineff for DetMC:   " << ftrackingIneff << endl;
    cout << "Unfolding with true MC set: " << iUnfJetClassTrue << endl;
    cout << "Unfolding with det  MC set: " << iUnfJetClassDet << endl;
    cout << "Use C-rho estimation:       " << bUseCrho << endl;
    cout << "Event selection flag:       " << flags << endl;
    cout << endl;

    if(fusePionMass && (fktScheme!=0 || fantiktScheme!=0)) {
        cout << "Warning: Using pion mass for jets but not using E_scheme!" << endl;
        cout << endl;
    }

    if(ftrackingIneff<0.0) {
        fYAMLConfig.Reinitialize();
        //if ftrackingIneff set as any negative number, use pt dependend ineff.
        SetArtificialTrackingEfficiencyFromYAML();
    }


#if !defined(__CINT__) && !defined(__MAKECINT__)
    //Note: for true MC we do not want to cut particles kinematically.
    fana->SetSettings(fDebug,
                      fparticleEtaCut,
                      fIsMC?0.0:fparticlePtCut,
                      fIsMC?DBL_MAX:100.0,
                      fjetCone,
                      fktJetCone,
                      fktScheme,
                      fantiktScheme,
                      fusePionMass,
                      fuseDeltaPhiBGSubtr,
                      fconstituentCut,
                      fleadingJetCut,
                      fsubleadingJetCut,
                      fMinJetPt,
                      fdeltaPhiCut,
                      fmatchingR,
                      0.0, //Tracking ineff only for det level.
                      nullptr, //Not needed if ^ not used
                      bUseCrho,
                      fIsMC); //Is this true mc

    if(fIsMC) {
        fanaMC->SetSettings(fDebug,
                            fparticleEtaCut,
                            fparticlePtCut,
                            100.0,
                            fjetCone,
                            fktJetCone,
                            fktScheme,
                            fantiktScheme,
                            fusePionMass,
                            fuseDeltaPhiBGSubtr,
                            fconstituentCut,
                            fleadingJetCut,
                            fsubleadingJetCut,
                            fMinJetPt,
                            fdeltaPhiCut,
                            fmatchingR,
                            ftrackingIneff,
                            fTrackEfficiencyHistogram,
                            bUseCrho,
                            false); //Is this true mc
    }

    // Save information about the settings used.
    // Done after SetSettings
    fana->InitHistos(fhistos, fIsMC, fcentralityBins.size(), iUnfJetClassTrue, iUnfJetClassDet);
    if(fIsMC) fanaMC->InitHistos(fhistosDetMC, fIsMC, fcentralityBins.size(), iUnfJetClassTrue, iUnfJetClassDet);

#endif

    // Load Custom Configuration and parameters
    // override values with parameters
    PostData(1, fOutput);

}

//______________________________________________________________________________
void AliJCDijetTask::UserExec(Option_t* /*option*/) 
{
    //cout << "======================== BEGIN EVENT ========================" << endl;
    // Processing of one event
    bGoodMCEvent = true;
    bGoodEvent = true;
    TClonesArray *fInputListDetMC;
    TClonesArray *fInputList;
    fhistos->fh_eventSel->Fill("events wo/ cuts",1.0);
    if(fDebug > 5) cout << "------- AliJCDijetTask Exec-------"<<endl;
    if(!((Entry()-1)%1000))  AliInfo(Form(" Processing event # %lld",  Entry())); 
    if( fJCatalystTask->GetJCatalystEntry() != fEntry) bGoodEvent=false;

    if(bGoodEvent) {
        fCBin = AliJCDijetHistos::GetCentralityClass(fJCatalystTask->GetCentrality());
        fhistos->fh_centrality->Fill(fJCatalystTask->GetCentrality());
        if(fCBin == -1) bGoodEvent=false;
    }


    if(fIsMC) {
        fhistosDetMC->fh_eventSel->Fill("events wo/ cuts",1.0);
        if(fJCatalystDetMCTask->GetJCatalystEntry() != fEntry) bGoodMCEvent=false;
        else fhistosDetMC->fh_eventSel->Fill("catalyst entry ok",1.0);

        if(bGoodMCEvent && (flags & DIJET_CATALYST)) {
            if( !fJCatalystDetMCTask->GetIsGoodEvent() ) bGoodMCEvent=false;
            else fhistosDetMC->fh_eventSel->Fill("catalyst ok",1.0);
        }

        if(bGoodMCEvent && (flags & DIJET_ALIEVENTCUT)) {
            if( !fEventCuts.AcceptEvent(InputEvent()) ) bGoodMCEvent=false;
            else fhistosDetMC->fh_eventSel->Fill("alieventcut ok",1.0);
        }

        if(bGoodMCEvent && (flags & DIJET_VERTEX13PA)) {
            if(!fUtils->IsVertexSelected2013pA(InputEvent())) bGoodMCEvent=false;
            else fhistosDetMC->fh_eventSel->Fill("vertex2013pA ok",1.0);
        }

        if(bGoodMCEvent && (flags & DIJET_PILEUPSPD)) {
            if(InputEvent()->IsPileupFromSPD(3,0.6,3,2,5)) bGoodMCEvent=false;
            else fhistosDetMC->fh_eventSel->Fill("pileupSPD ok",1.0);
        }

        if(bGoodMCEvent && (flags & DIJET_UTILSPILEUPSPD)) {
            if(fUtils->IsPileUpSPD(InputEvent())) bGoodMCEvent=false;
            else fhistosDetMC->fh_eventSel->Fill("utils pileupSPD ok",1.0);
        }

        if(bGoodMCEvent) {

            fCBinDetMC = AliJCDijetHistos::GetCentralityClass(fJCatalystTask->GetCentrality());
            fhistosDetMC->fh_centrality->Fill(fJCatalystTask->GetCentrality());
            if(fCBinDetMC == -1) bGoodMCEvent=false;
        }
        if(bGoodMCEvent) {

            TList *genHeaders = 0x0;
            AliGenEventHeader* gh = 0;
            AliMCEvent *mcEvent = MCEvent();
            if(mcEvent) genHeaders = mcEvent->GetCocktailList();
            if(genHeaders){
                for(Int_t i = 0; i<genHeaders->GetEntries(); i++){
                    gh = (AliGenEventHeader*)genHeaders->At(i);
                    AliGenPythiaEventHeader* pythiaGenHeader= dynamic_cast<AliGenPythiaEventHeader*>(gh); //identify pythia header
                    //AliGenPythiaEventHeader *pythiaGenHeader = AliAnalysisHelperJetTasks::GetPythiaEventHeader(mcEvent);
                    if(pythiaGenHeader) {
                        fptHardBin = pythiaGenHeader->GetPtHard();
                        fPythiaSigma = pythiaGenHeader->GetXsection();
                        fPythiaTrial = pythiaGenHeader->Trials();
                    }
                    else {
                        fptHardBin = 0.0;
                        fPythiaSigma = 0.0;
                        fPythiaTrial = 0.0;
                    }
                    //cout << "genHeader " << i << ": " << pythiaGenHeader << ", fptHardBin: " << fptHardBin << endl;
                }
            }
            //cout << "fptHardBin: " << fptHardBin << endl;
            //fana->SetPythiaInfo(fptHardBin, fPythiaSigma, fPythiaTrial);
            fanaMC->SetPythiaInfo(fptHardBin, fPythiaSigma, fPythiaTrial);

            fhistosDetMC->fh_eventSel->Fill("events",1.0);
            fhistosDetMC->fh_events[fCBin]->Fill("events",1.0);
            fhistosDetMC->fh_zvtx->Fill(fJCatalystTask->GetZVertex());

            fInputListDetMC = (TClonesArray*)fJCatalystDetMCTask->GetInputList();

        }


#if !defined(__CINT__) && !defined(__MAKECINT__)
        if(bGoodMCEvent) {
            //cout << "Next det level calculations:" << endl;
            fDetMCFlag = fanaMC->CalculateJets(fInputListDetMC, fhistosDetMC, fCBinDetMC, 1.0);
            //If fDetMCFlag=-1 then we want to discard whole event
            //as that means there is jet-pt > 4*pt_hard
            if(fDetMCFlag != 0) { 
                bGoodMCEvent=false;
                bGoodEvent=false;
            }
        }
        if(bGoodMCEvent) {
            fanaMC->FillJetsDijets(fhistosDetMC, fCBinDetMC, 1.0);
        }
#endif
    }

    fhistos->fh_eventSel->Fill("catalyst entry ok",1.0);
    if(bGoodEvent && (flags & DIJET_CATALYST)) {
        if( !fJCatalystTask->GetIsGoodEvent() ) bGoodEvent=false;
        else fhistos->fh_eventSel->Fill("catalyst ok",1.0);
    }

    if(bGoodEvent && (flags & DIJET_ALIEVENTCUT)) {
        if( !fEventCuts.AcceptEvent(InputEvent()) ) bGoodEvent=false;
        else fhistos->fh_eventSel->Fill("alieventcut ok",1.0);
    }

    if(bGoodEvent && (flags & DIJET_VERTEX13PA)) {
        if(!fUtils->IsVertexSelected2013pA(InputEvent())) bGoodEvent=false;
        else fhistos->fh_eventSel->Fill("vertex2013pA ok",1.0);
    }

    if(bGoodEvent && (flags & DIJET_PILEUPSPD)) {
        if(InputEvent()->IsPileupFromSPD(3,0.6,3,2,5)) bGoodEvent=false;
        else fhistos->fh_eventSel->Fill("pileupSPD ok",1.0);
    }

    if(bGoodEvent && (flags & DIJET_UTILSPILEUPSPD)) {
        if(fUtils->IsPileUpSPD(InputEvent())) bGoodEvent=false;
        else fhistos->fh_eventSel->Fill("utils pileupSPD ok",1.0);
    }

    if(bGoodEvent) {

        fhistos->fh_eventSel->Fill("events",1.0);
        fhistos->fh_events[fCBin]->Fill("events",1.0);
        fhistos->fh_zvtx->Fill(fJCatalystTask->GetZVertex());

        fInputList = (TClonesArray*)fJCatalystTask->GetInputList();
    }

#if !defined(__CINT__) && !defined(__MAKECINT__)
    //cout << "Next true level calculations:" << endl;
    if(bGoodEvent) {
        fana->CalculateJets(fInputList, fhistos, fCBin, 1.0);
        fana->FillJetsDijets(fhistos, fCBin, 1.0);
    }
#endif

    // Here response matrix calculation.
    if(fIsMC && bGoodEvent && bGoodMCEvent) {
#if !defined(__CINT__) && !defined(__MAKECINT__)
        fana->CalculateResponse(fanaMC,fhistosDetMC,AliJCDijetAna::iAcc,AliJCDijetAna::iAcc);
        fana->CalculateResponse(fanaMC,fhistosDetMC,AliJCDijetAna::iBGSubtr,AliJCDijetAna::iBGSubtr);
        fana->CalculateResponse(fanaMC,fhistosDetMC,AliJCDijetAna::iBGSubtrConstCut,AliJCDijetAna::iBGSubtrConstCut);
        fana->CalculateResponse(fanaMC,fhistosDetMC,AliJCDijetAna::iConstCut,AliJCDijetAna::iConstCut);
        fana->CalculateResponse(fanaMC,fhistosDetMC,AliJCDijetAna::iBGSubtrCutsRaw,AliJCDijetAna::iBGSubtrCutsRaw);
        fana->CalculateResponse(fanaMC,fhistosDetMC,AliJCDijetAna::iBGSubtrIndEta,AliJCDijetAna::iBGSubtrIndEta);
        fana->CalculateResponse(fanaMC,fhistosDetMC,AliJCDijetAna::iBGSubtrConstCutIndEta,AliJCDijetAna::iBGSubtrConstCutIndEta);
        fana->CalculateResponse(fanaMC,fhistosDetMC,AliJCDijetAna::iBGSubtrCommonEta,AliJCDijetAna::iBGSubtrCommonEta);
        fana->CalculateResponse(fanaMC,fhistosDetMC,AliJCDijetAna::iBGSubtrConstCutCommonEta,AliJCDijetAna::iBGSubtrConstCutCommonEta);
        //we can run custom configuration with an argument
        if( !((iUnfJetClassTrue==AliJCDijetAna::iAcc                      && iUnfJetClassDet==AliJCDijetAna::iAcc)
        ||    (iUnfJetClassTrue==AliJCDijetAna::iBGSubtr                  && iUnfJetClassDet==AliJCDijetAna::iBGSubtr)
        ||    (iUnfJetClassTrue==AliJCDijetAna::iBGSubtrConstCut          && iUnfJetClassDet==AliJCDijetAna::iBGSubtrConstCut)
        ||    (iUnfJetClassTrue==AliJCDijetAna::iConstCut                 && iUnfJetClassDet==AliJCDijetAna::iConstCut)
        ||    (iUnfJetClassTrue==AliJCDijetAna::iBGSubtrCutsRaw           && iUnfJetClassDet==AliJCDijetAna::iBGSubtrCutsRaw) 
        ||    (iUnfJetClassTrue==AliJCDijetAna::iBGSubtrIndEta            && iUnfJetClassDet==AliJCDijetAna::iBGSubtrIndEta) 
        ||    (iUnfJetClassTrue==AliJCDijetAna::iBGSubtrConstCutIndEta    && iUnfJetClassDet==AliJCDijetAna::iBGSubtrConstCutIndEta) 
        ||    (iUnfJetClassTrue==AliJCDijetAna::iBGSubtrCommonEta         && iUnfJetClassDet==AliJCDijetAna::iBGSubtrCommonEta) 
        ||    (iUnfJetClassTrue==AliJCDijetAna::iBGSubtrConstCutCommonEta && iUnfJetClassDet==AliJCDijetAna::iBGSubtrConstCutCommonEta) 
             )) {
            fana->CalculateResponse(fanaMC,fhistosDetMC,iUnfJetClassTrue,iUnfJetClassDet);
        }
#endif
    }


}

//______________________________________________________________________________
void AliJCDijetTask::Init()
{
    // Intialisation of parameters
    AliInfo("Doing initialization") ; 

}

//______________________________________________________________________________
void AliJCDijetTask::Terminate(Option_t *)
{
    // Processing when the event loop is ended
    cout<<"AliJCDijetTask Analysis DONE !!"<<endl; 
}

/**
 * Stream and initialise tracking efficiency yaml file
 * Original implementation in AliEmcalJetTask
 */
void AliJCDijetTask::AddArtificialTrackingEfficiencyConfig(TString sAnchorPeriod) {
    fsAnchorPeriod=sAnchorPeriod;
    std::string path = "$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/TrackEfficiencyConfiguration.yaml";
    Printf("Get pT-dependent Tracking efficiency from %s", path.c_str());
    int addedConfig = fYAMLConfig.AddConfiguration(path, "yamlConfig");
    if (addedConfig < 0) {
        AliFatal(TString::Format("YAML Configuration in set path %s not found!",path.c_str()).Data());
    }
    fYAMLConfig.Initialize();
}

/**
 * Set the pt-dependent tracking efficiency from the loaded YAML file
 * Original implementation in AliEmcalJetTask
 */
void AliJCDijetTask::SetArtificialTrackingEfficiencyFromYAML() {
  std::vector <Double_t> ptBinning;
  std::vector <Double_t> trackingUncertainty;
  bool res = fYAMLConfig.GetProperty("ptBinning", ptBinning, false);
  Int_t nPtBins = ptBinning.size()-1;
  double* aptBinning = ptBinning.data();

  //TODO: could implement automatic anchor period search.
  //auto userInfo = fInputHandler->GetUserInfo();
  //AliProdInfo prodInfo(userInfo);
  //std::string period = prodInfo.GetAnchorProduction().Data();
  std::string period = fsAnchorPeriod.Data();
  std::string cent = "0_100"; //Here I have implemented MB only. Search AliEmcalJetTask for cent binned implementation.
  //AliInfoStream() << "anchor production = " << prodInfo.GetAnchorProduction()<< "\n";

  res = fYAMLConfig.GetProperty({period,cent},trackingUncertainty, false);
  if(res) {
      fTrackEfficiencyHistogram = new TH1D("fTrackEfficiencyHistogram","h",nPtBins,aptBinning);
      for(Int_t i = 0; i < nPtBins; i++) {
          fTrackEfficiencyHistogram->SetBinContent(i+1, trackingUncertainty.at(i));
          AliDebug(2,TString::Format("pT %f - %f \t track uncertainty: %f", ptBinning.at(i), ptBinning.at(i+1), trackingUncertainty.at(i)).Data());
      }
  }
  else {
      fTrackEfficiencyHistogram = nullptr;
      AliFatal("not able to find any pt-dependent uncertainties for the anchored period %s of the MC that you are running over");
  }
}

//_____________________________________________________________________
AliAnalysisTask *AliJCDijetTask::AddTaskJCDijetTask(TString taskName,
                                    Bool_t isMC,
                                    TString sJCatalyst        ,
                                    TString sJCatalystDetMC   ,
                                    UInt_t flags              ,
                                    TString centBins          ,
                                    TString sDijetMBins       ,
                                    double jetCone            ,
                                    double ktjetCone          ,
                                    int ktScheme              ,
                                    int antiktScheme          ,
                                    Bool_t usePionMass        ,
                                    Bool_t useDeltaPhiBGSubtr ,
                                    double particleEtaCut     ,
                                    double particlePtCut      ,
                                    double leadingJetCut      ,
                                    double subleadingJetCut   ,
                                    double minJetPt           ,
                                    double constituentCut     ,
                                    double deltaPhiCut        ,
                                    double matchingR          ,
                                    double trackingIneff      ,
                                    TString sAnchorPeriodForTracking ,
                                    AliJCDijetAna::jetClasses lUnfJetClassTrue ,
                                    AliJCDijetAna::jetClasses lUnfJetClassDet ,
                                    Bool_t useCoveredAreaRho)
{
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

    cout<<"AddTaskJCDijetTask::flags = "<<flags<<endl;

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

    std::stringstream ss( centBins.Data() );
    double binBorder;
    vector<double> vecCentBins;
    ss >> binBorder;
    while (!ss.fail()) {
        vecCentBins.push_back(binBorder);
        ss >> binBorder;
    }

    if (vecCentBins.size() < 2) {
        ::Error("AddTaskJCDijetTask", "Centrality bins are not properly set. At least two bin borders are needed. Terminating.");
        return NULL;
    }

    for (int ivec=0; ivec < vecCentBins.size()-1; ivec++) {
        if(vecCentBins.at(ivec+1) <= vecCentBins.at(ivec)) {
            ::Error("AddTaskJCDijetTask", "Centrality bins are not properly set. Terminating.");
            return NULL;
        }
    }

    if (jetCone > 0.8 || jetCone < 0.0 || ktjetCone > 0.8 || ktjetCone < 0.0) {
        ::Error("AddTaskJCDijetTask", "Jet cones are set to be too small or too big. Terminating.");
        return NULL;
    }

    if (ktScheme < 0 || ktScheme > 6) {
        ::Error("AddTaskJCDijetTask", "Invalid ktScheme set. Please choose a setting from 0 till 6. Terminating.");
        return NULL;
    }
    if (antiktScheme < 0 || antiktScheme > 6) {
        ::Error("AddTaskJCDijetTask", "Invalid antiktScheme set. Please choose a setting from 0 till 6. Terminating.");
        return NULL;
    }
    cout << "MC: " << isMC << endl;

    //==== Set up the dijet task ====
    AliJCDijetTask *dijetTask = new AliJCDijetTask(taskName.Data(),"AOD");
    dijetTask->SetDebugLevel(5);
    dijetTask->SetJCatalystTaskName(sJCatalyst.Data());
    dijetTask->SetJCatalystTaskNameDetMC(sJCatalystDetMC.Data());
    dijetTask->SetCentralityBins(vecCentBins);
    dijetTask->SetDijetMBins(sDijetMBins);
    dijetTask->SetJetConeSize(jetCone, ktjetCone);
    dijetTask->SetBGSubtrSettings(ktScheme, antiktScheme, usePionMass, useDeltaPhiBGSubtr, useCoveredAreaRho);
    dijetTask->SetUnfoldingJetSets(lUnfJetClassTrue, lUnfJetClassDet);
    dijetTask->SetIsMC(isMC);
    dijetTask->SetCuts(particleEtaCut, particlePtCut, leadingJetCut, subleadingJetCut, constituentCut, deltaPhiCut, matchingR, trackingIneff, minJetPt);
    dijetTask->AddFlags(flags);
    if(trackingIneff<0.0) dijetTask->AddArtificialTrackingEfficiencyConfig(sAnchorPeriodForTracking);
    cout << dijetTask->GetName() << endl;


    mgr->AddTask((AliAnalysisTask*) dijetTask);

    // Create containers for input/output
    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();


    // Connect input/output
    mgr->ConnectInput(dijetTask, 0, cinput);
    AliAnalysisDataContainer *jHist = mgr->CreateContainer(Form("%scontainer",dijetTask->GetName()),  TDirectory::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s",AliAnalysisManager::GetCommonFileName(), dijetTask->GetName()));
    mgr->ConnectOutput(dijetTask, 1, jHist );

    return dijetTask;
}

