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
    fana(NULL),
    fanaMC(NULL),
    fCBin(-1),
    fCBinDetMC(-1),
    fOutput(NULL),
    flags(0),
    fUtils(nullptr)
{
    // Constructor
    AliInfo("---- AliJCDijetTask Constructor ----");
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
    delete fana;
    delete fanaMC;
    delete fUtils;

}

//________________________________________________________________________

void AliJCDijetTask::UserCreateOutputObjects()
{  
    //=== create the jcorran outputs objects
    if(fDebug > 1) printf("AliJCDijetTask::UserCreateOutPutData() \n");
    //=== Get AnalysisManager
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();

    fJCatalystTask = (AliJCatalystTask*)(man->GetTask( fJCatalystTaskName ));
    if(fIsMC) fJCatalystDetMCTask = (AliJCatalystTask*)(man->GetTask( fJCatalystDetMCTaskName ));

    OpenFile(1);
    fOutput = gDirectory;
    fOutput->cd();

    fhistos = new AliJCDijetHistos();
    fhistos->SetName("jcdijet");
    fhistos->SetCentralityBinsHistos(fcentralityBins);
    fhistos->CreateEventTrackHistos();
    fhistos->fHMG->Print();
    fana = new AliJCDijetAna();

    if(fIsMC) {
        fhistosDetMC = new AliJCDijetHistos();
        fhistosDetMC->SetName("jcdijetDetMC");
        fhistosDetMC->SetCentralityBinsHistos(fcentralityBins);
        fhistosDetMC->CreateEventTrackHistos();
        fhistosDetMC->fHMG->Print();
        fanaMC = new AliJCDijetAna();
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
    cout << "Jet cone size:              " << fjetCone << endl;
    cout << "kt-jet cone size:           " << fktJetCone << endl;
    cout << "Using kt-jet scheme:        " << sktScheme.Data() << endl;
    cout << "Using antikt-jet scheme:    " << santiktScheme.Data() << endl;
    cout << "Using pion mass:            " << fusePionMass << endl;
    cout << "Using DeltaPhi in BG subtr: " << fuseDeltaPhiBGSubtr << endl;
    cout << "Particle eta cut:           " << fparticleEtaCut << endl;
    cout << "Particle pt cut:            " << fparticlePtCut << endl;
    cout << "Dijet leading jet cut:      " << fleadingJetCut << endl;
    cout << "Dijet subleading jet cut:   " << fsubleadingJetCut << endl;
    cout << "Jet min pt cut:             " << fMinJetPt << endl;
    cout << "Jet leading const. cut:     " << fconstituentCut << endl;
    cout << "Dijet DeltaPhi cut:         pi/" << fdeltaPhiCut << endl;
    cout << "Matching R for MC:          " << fmatchingR << endl;
    cout << "Tracking ineff for DetMC:   " << ftrackingIneff << endl;
    cout << endl;

    if(fusePionMass && (fktScheme!=0 || fantiktScheme!=0)) {
        cout << "Warning: Using pion mass for jets but not using E_scheme!" << endl;
        cout << endl;
    }

#if !defined(__CINT__) && !defined(__MAKECINT__)
    fana->SetSettings(fDebug,
                      fparticleEtaCut,
                      fparticlePtCut,
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
                      0.0); //Tracking ineff only for det level.

    if(fIsMC) {
        fanaMC->SetSettings(fDebug,
                            fparticleEtaCut,
                            fparticlePtCut,
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
                            ftrackingIneff);
    }

    // Save information about the settings used.
    // Done after SetSettings
    fana->InitHistos(fhistos, fIsMC, fcentralityBins.size());
    if(fIsMC) fanaMC->InitHistos(fhistosDetMC, fIsMC, fcentralityBins.size());

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
    fhistos->fh_eventSel->Fill("events wo/ cuts",1.0);
    if(fDebug > 5) cout << "------- AliJCDijetTask Exec-------"<<endl;
    if(!((Entry()-1)%1000))  AliInfo(Form(" Processing event # %lld",  Entry())); 
    if( fJCatalystTask->GetJCatalystEntry() != fEntry) return;
    fhistos->fh_eventSel->Fill("catalyst entry ok",1.0);
    if( !fJCatalystTask->GetIsGoodEvent() ) return;
    fhistos->fh_eventSel->Fill("catalyst ok",1.0);

    if(flags & DIJET_VERTEX13PA) {
        if(!fUtils->IsVertexSelected2013pA(InputEvent())) return;
        fhistos->fh_eventSel->Fill("vertex2013pA ok",1.0);
    }

    if(flags & DIJET_PILEUPSPD) {
        if(InputEvent()->IsPileupFromSPD(3,0.6,3,2,5)) return;
        fhistos->fh_eventSel->Fill("pileupSPD ok",1.0);
    }

    if(flags & DIJET_UTILSPILEUPSPD) {
        if(fUtils->IsPileUpSPD(InputEvent())) return;
        fhistos->fh_eventSel->Fill("utils pileupSPD ok",1.0);
    }

    fCBin = AliJCDijetHistos::GetCentralityClass(fJCatalystTask->GetCentrality());
    fhistos->fh_centrality->Fill(fJCatalystTask->GetCentrality());
    if(fCBin == -1) return;

    fhistos->fh_eventSel->Fill("events",1.0);
    fhistos->fh_events[fCBin]->Fill("events",1.0);
    fhistos->fh_zvtx->Fill(fJCatalystTask->GetZVertex());
    
    TClonesArray *fInputList = (TClonesArray*)fJCatalystTask->GetInputList();

#if !defined(__CINT__) && !defined(__MAKECINT__)
    //cout << "Next true level calculations:" << endl;
    fana->CalculateJets(fInputList, fhistos, fCBin);
    fana->FillJetsDijets(fhistos, fCBin);
#endif

    if(fIsMC) {
        fhistosDetMC->fh_eventSel->Fill("events wo/ cuts",1.0);
        if(fJCatalystDetMCTask->GetJCatalystEntry() != fEntry) return;
        fhistosDetMC->fh_eventSel->Fill("catalyst entry ok",1.0);
        if( !fJCatalystDetMCTask->GetIsGoodEvent() ) return;
        fhistosDetMC->fh_eventSel->Fill("catalyst ok",1.0);

        if(flags & DIJET_VERTEX13PA) {
            if(!fUtils->IsVertexSelected2013pA(InputEvent())) return;
            fhistosDetMC->fh_eventSel->Fill("vertex2013pA ok",1.0);
        }

        if(flags & DIJET_PILEUPSPD) {
            if(InputEvent()->IsPileupFromSPD(3,0.6,3,2,5)) return;
            fhistosDetMC->fh_eventSel->Fill("pileupSPD ok",1.0);
        }

        if(flags & DIJET_UTILSPILEUPSPD) {
            if(fUtils->IsPileUpSPD(InputEvent())) return;
            fhistosDetMC->fh_eventSel->Fill("utils pileupSPD ok",1.0);
        }

        fCBinDetMC = AliJCDijetHistos::GetCentralityClass(fJCatalystTask->GetCentrality());
        fhistosDetMC->fh_centrality->Fill(fJCatalystTask->GetCentrality());
        if(fCBinDetMC == -1) return;

        fhistosDetMC->fh_eventSel->Fill("events",1.0);
        fhistosDetMC->fh_events[fCBin]->Fill("events",1.0);
        fhistosDetMC->fh_zvtx->Fill(fJCatalystTask->GetZVertex());

        TClonesArray *fInputListDetMC = (TClonesArray*)fJCatalystDetMCTask->GetInputList();

#if !defined(__CINT__) && !defined(__MAKECINT__)
        //cout << "Next det level calculations:" << endl;
        fanaMC->CalculateJets(fInputListDetMC, fhistosDetMC, fCBinDetMC);
        fanaMC->FillJetsDijets(fhistosDetMC, fCBinDetMC);

        // Here response matrix calculation.
        fana->CalculateResponse(fanaMC,fhistosDetMC);
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

