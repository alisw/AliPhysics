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
    fusePionMassInktjets(0),
    fuseDeltaPhiBGSubtr(0),
    fIsMC(kTRUE),
    fparticleEtaCut(0),
    fleadingJetCut(0),
    fsubleadingJetCut(0),
    fconstituentCut(0),
    fdeltaPhiCut(0),
    fhistos(NULL),
    fhistosDetMC(NULL),
    fana(NULL),
    fanaMC(NULL),
    fCBin(-1),
    fCBinDetMC(-1),
    fOutput(NULL)
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
    fusePionMassInktjets(0),
    fuseDeltaPhiBGSubtr(0),
    fIsMC(kTRUE),
    fparticleEtaCut(0),
    fleadingJetCut(0),
    fsubleadingJetCut(0),
    fconstituentCut(0),
    fdeltaPhiCut(0),
    fhistos(NULL),
    fhistosDetMC(NULL),
    fana(NULL),
    fanaMC(NULL),
    fCBin(-1),
    fCBinDetMC(-1),
    fOutput(NULL)
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
    fusePionMassInktjets(ap.fusePionMassInktjets),
    fuseDeltaPhiBGSubtr(ap.fuseDeltaPhiBGSubtr),
    fIsMC(ap.fIsMC),
    fparticleEtaCut(ap.fparticleEtaCut),
    fleadingJetCut(ap.fleadingJetCut),
    fsubleadingJetCut(ap.fsubleadingJetCut),
    fconstituentCut(ap.fconstituentCut),
    fdeltaPhiCut(ap.fdeltaPhiCut),
    fhistos(ap.fhistos),
    fhistosDetMC(ap.fhistosDetMC),
    fana(ap.fana),
    fanaMC(ap.fanaMC),
    fCBin(ap.fCBin),
    fCBinDetMC(ap.fCBinDetMC),
    fOutput( ap.fOutput )
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


    TString sktScheme;
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

    cout << endl;
    cout << "===========SETTINGS===========" << endl;
    cout << "MC:                         " << IsMC() << endl;
    cout << "Centrality bins:            ";
    for(unsigned i=0; i< fcentralityBins.size(); i++) cout << fcentralityBins.at(i) << " ";
    cout << endl;
    cout << "Jet cone size:              " << fjetCone << endl;
    cout << "kt-jet cone size:           " << fktJetCone << endl;
    cout << "Using kt-jet scheme:        " << sktScheme.Data() << endl;
    cout << "Using pion mass in kt-jets: " << fusePionMassInktjets << endl;
    cout << "Using DeltaPhi in BG subtr: " << fuseDeltaPhiBGSubtr << endl;
    cout << "Particle eta cut:           " << fparticleEtaCut << endl;
    cout << "Particle pt cut:            " << fparticlePtCut << endl;
    cout << "Dijet leading jet cut:      " << fleadingJetCut << endl;
    cout << "Dijet subleading jet cut:   " << fsubleadingJetCut << endl;
    cout << "Jet leading const. cut:     " << fconstituentCut << endl;
    cout << "Dijet DeltaPhi cut :        pi/" << fdeltaPhiCut << endl;
    cout << endl;

    if(fusePionMassInktjets && fktScheme!=0) {
        cout << "Warning: Using pion mass for kt-jets but not using E_scheme!" << endl;
        cout << endl;
    }

    // Save information about the settings used.
    InitHistos(fhistos);
    if(fIsMC) InitHistos(fhistosDetMC);

#if !defined(__CINT__) && !defined(__MAKECINT__)
    fana->SetSettings(fDebug,
                      fparticleEtaCut,
                      fparticlePtCut,
                      fjetCone,
                      fktJetCone,
                      fktScheme,
                      fusePionMassInktjets,
                      fuseDeltaPhiBGSubtr,
                      fconstituentCut,
                      fleadingJetCut,
                      fsubleadingJetCut,
                      fdeltaPhiCut);

    if(fIsMC) {
        fanaMC->SetSettings(fDebug,
                            fparticleEtaCut,
                            fparticlePtCut,
                            fjetCone,
                            fktJetCone,
                            fktScheme,
                            fusePionMassInktjets,
                            fuseDeltaPhiBGSubtr,
                            fconstituentCut,
                            fleadingJetCut,
                            fsubleadingJetCut,
                            fdeltaPhiCut);
    }
#endif

    // Load Custom Configuration and parameters
    // override values with parameters
    PostData(1, fOutput);

}

//______________________________________________________________________________
void AliJCDijetTask::UserExec(Option_t* /*option*/) 
{
    // Processing of one event
    if(fDebug > 5) cout << "------- AliJCDijetTask Exec-------"<<endl;
    if(!((Entry()-1)%1000))  AliInfo(Form(" Processing event # %lld",  Entry())); 
    if( fJCatalystTask->GetJCatalystEntry() != fEntry) return;
    fCBin = AliJCDijetHistos::GetCentralityClass(fJCatalystTask->GetCentrality());
    if(fCBin == -1) return;

    fhistos->fh_events[fCBin]->Fill("events",1.0);
    fhistos->fh_centrality->Fill(fJCatalystTask->GetCentrality());
    fhistos->fh_zvtx->Fill(fJCatalystTask->GetZVertex());

    TClonesArray *fInputList = (TClonesArray*)fJCatalystTask->GetInputList();

#if !defined(__CINT__) && !defined(__MAKECINT__)
    fana->CalculateJetsDijets(fInputList, fhistos, fCBin);
#endif

    if(fIsMC) {
        if(fJCatalystDetMCTask->GetJCatalystEntry() != fEntry) return;
        fCBinDetMC = AliJCDijetHistos::GetCentralityClass(fJCatalystTask->GetCentrality());
        if(fCBinDetMC == -1) return;

        fhistosDetMC->fh_events[fCBin]->Fill("events",1.0);
        fhistosDetMC->fh_centrality->Fill(fJCatalystTask->GetCentrality());
        fhistosDetMC->fh_zvtx->Fill(fJCatalystTask->GetZVertex());

        TClonesArray *fInputListDetMC = (TClonesArray*)fJCatalystDetMCTask->GetInputList();

#if !defined(__CINT__) && !defined(__MAKECINT__)
        fanaMC->CalculateJetsDijets(fInputListDetMC, fhistosDetMC, fCBinDetMC);
#endif


        // Here response matrix calculation.
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

void AliJCDijetTask::InitHistos(AliJCDijetHistos *histos) {
    histos->fh_info->Fill("Count", 1.0);
    histos->fh_info->Fill("MC", IsMC());
    for(unsigned i=0; i< fcentralityBins.size(); i++) histos->fh_info->Fill(Form("Cent bin border %02d",i), fcentralityBins.at(i));
    histos->fh_info->Fill("Jet cone", fjetCone);
    histos->fh_info->Fill("kt-jet cone", fktJetCone);
    histos->fh_info->Fill("Scheme", fktScheme);
    histos->fh_info->Fill("Use pion mass", fusePionMassInktjets);
    histos->fh_info->Fill("Use DeltaPhi BG Subtr", fuseDeltaPhiBGSubtr);
    histos->fh_info->Fill("Particle eta cut", fparticleEtaCut);
    histos->fh_info->Fill("Particle pt cut", fparticlePtCut);
    histos->fh_info->Fill("Leading jet cut", fleadingJetCut);
    histos->fh_info->Fill("Subleading jet cut", fsubleadingJetCut);
    histos->fh_info->Fill("Const. cut", fconstituentCut);
    histos->fh_info->Fill("Delta phi cut pi/",fdeltaPhiCut);

    // Initialize fh_events so that the bin order is correct
    for (unsigned iBin=0; iBin < fcentralityBins.size()-1; iBin++) {
        histos->fh_events[iBin]->Fill("events",0.0);
        histos->fh_events[iBin]->Fill("particles",0.0);
        histos->fh_events[iBin]->Fill("acc. particles",0.0);
        histos->fh_events[iBin]->Fill("no rho calc. events",0.0);
        histos->fh_events[iBin]->Fill("rho calc. events",0.0);
        histos->fh_events[iBin]->Fill("jets",0.0);
        histos->fh_events[iBin]->Fill("acc. jets",0.0);
        histos->fh_events[iBin]->Fill("const. cut jets",0.0);
        histos->fh_events[iBin]->Fill("bg. subtr. jets",0.0);
        histos->fh_events[iBin]->Fill("bg. subtr. const. cut jets",0.0);
        histos->fh_events[iBin]->Fill("kt-jets",0.0);
        histos->fh_events[iBin]->Fill("acc. kt-jets",0.0);
        histos->fh_events[iBin]->Fill("leading jet drop",0.0);
        histos->fh_events[iBin]->Fill("subleading jet drop",0.0);
        histos->fh_events[iBin]->Fill("raw dijets",0.0);
        histos->fh_events[iBin]->Fill("raw dijets leading cut",0.0);
        histos->fh_events[iBin]->Fill("raw acc. dijets",0.0);
        histos->fh_events[iBin]->Fill("raw deltaphi cut dijets",0.0);
        histos->fh_events[iBin]->Fill("bg. subtr. dijets",0.0);
        histos->fh_events[iBin]->Fill("bg. subtr. dijets leading cut",0.0);
        histos->fh_events[iBin]->Fill("bg. subtr. acc. dijets",0.0);
        histos->fh_events[iBin]->Fill("bg. subtr. deltaphi cut dijets",0.0);
        histos->fh_events[iBin]->Fill("bg. subtr. const. cut dijets",0.0);
        histos->fh_events[iBin]->Fill("bg. subtr. const. cut dijets leading cut",0.0);
        histos->fh_events[iBin]->Fill("bg. subtr. const. cut acc. dijets",0.0);
        histos->fh_events[iBin]->Fill("bg. subtr. const. cut deltaphi cut dijets",0.0);
        histos->fh_events[iBin]->Fill("const. cut dijets",0.0);
        histos->fh_events[iBin]->Fill("const. cut dijets leading cut",0.0);
        histos->fh_events[iBin]->Fill("const. cut acc. dijets",0.0);
        histos->fh_events[iBin]->Fill("const. cut deltaphi cut dijets",0.0);
        histos->fh_events[iBin]->Fill("kt dijets",0.0);
        histos->fh_events[iBin]->Fill("kt dijets leading cut",0.0);
        histos->fh_events[iBin]->Fill("kt acc. dijets",0.0);
        histos->fh_events[iBin]->Fill("kt deltaphi cut dijets",0.0);
    }
}
