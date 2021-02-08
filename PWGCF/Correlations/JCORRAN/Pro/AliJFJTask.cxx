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
#include "AliJFJTask.h" 

//______________________________________________________________________________
AliJFJTask::AliJFJTask() :
    AliAnalysisTaskSE("JFJTask"),
    fJCatalystTask(NULL),
    fJCatalystTaskName("JCatalystTask"),
    fcentralityBins(0),
    fjetCone(0),
    fktJetCone(0),
    fktScheme(0),
    fantiktScheme(0),
    fusePionMass(0),
    fuseDeltaPhiBGSubtr(0),
    fparticleEtaCut(0),
    fleadingJetCut(0),
    fsubleadingJetCut(0),
    fMinJetPt(0),
    fconstituentCut(0),
    fdeltaPhiCut(0),
    fmatchingR(0),
    fhistos(NULL),
    fana(NULL),
    fCBin(-1),
    fOutput(NULL),
    flags(0)
{
}

//______________________________________________________________________________
AliJFJTask::AliJFJTask(const char *name, TString inputformat):
    AliAnalysisTaskSE(name), 
    fJCatalystTask(NULL),
    fJCatalystTaskName("JCatalystTask"),
    fcentralityBins(0),
    fjetCone(0),
    fktJetCone(0),
    fktScheme(0),
    fantiktScheme(0),
    fusePionMass(0),
    fuseDeltaPhiBGSubtr(0),
    fparticleEtaCut(0),
    fleadingJetCut(0),
    fsubleadingJetCut(0),
    fMinJetPt(0),
    fconstituentCut(0),
    fdeltaPhiCut(0),
    fmatchingR(0),
    fhistos(NULL),
    fana(NULL),
    fCBin(-1),
    fOutput(NULL)
{
    // Constructor
    AliInfo("---- AliJFJTask Constructor ----");
    DefineOutput (1, TDirectory::Class());
}

//____________________________________________________________________________
AliJFJTask::AliJFJTask(const AliJFJTask& ap) :
    AliAnalysisTaskSE(ap.GetName()), 
    fJCatalystTask(ap.fJCatalystTask),
    fJCatalystTaskName(ap.fJCatalystTaskName),
    fcentralityBins(ap.fcentralityBins),
    fjetCone(ap.fjetCone),
    fktJetCone(ap.fktJetCone),
    fktScheme(ap.fktScheme),
    fantiktScheme(ap.fantiktScheme),
    fusePionMass(ap.fusePionMass),
    fuseDeltaPhiBGSubtr(ap.fuseDeltaPhiBGSubtr),
    fparticleEtaCut(ap.fparticleEtaCut),
    fleadingJetCut(ap.fleadingJetCut),
    fsubleadingJetCut(ap.fsubleadingJetCut),
    fMinJetPt(ap.fMinJetPt),
    fconstituentCut(ap.fconstituentCut),
    fdeltaPhiCut(ap.fdeltaPhiCut),
    fmatchingR(ap.fmatchingR),
    fhistos(ap.fhistos),
    fana(ap.fana),
    fCBin(ap.fCBin),
    fOutput( ap.fOutput )
{ 

    AliInfo("----DEBUG AliJFJTask COPY ----");

}

//_____________________________________________________________________________
AliJFJTask& AliJFJTask::operator = (const AliJFJTask& ap)
{
    // assignment operator

    AliInfo("----DEBUG AliJFJTask operator= ----");
    this->~AliJFJTask();
    new(this) AliJFJTask(ap);
    return *this;
}

//______________________________________________________________________________
AliJFJTask::~AliJFJTask()
{
    // destructor 
    delete fOutput;
    delete fhistos;
    delete fana;
}

//________________________________________________________________________

void AliJFJTask::UserCreateOutputObjects()
{  
    //=== create the jcorran outputs objects
    if(fDebug > 1) printf("AliJFJTask::UserCreateOutPutData() \n");
    //=== Get AnalysisManager
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();

    fJCatalystTask = (AliJCatalystTask*)(man->GetTask( fJCatalystTaskName ));
    if(fDebug > 1) cout << "fJCatalystTaskName = " << fJCatalystTaskName << endl;

    OpenFile(1);
    fOutput = gDirectory;
    fOutput->cd();

    fhistos = new AliJCDijetHistos();
    fhistos->SetName("jcdijet");
    fhistos->SetCentralityBinsHistos(fcentralityBins);
    fhistos->CreateEventTrackHistos();
    fhistos->fHMG->Print();
    fana = new AliJCDijetAna();

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
#endif
    // Load Custom Configuration and parameters
    // override values with parameters
    PostData(1, fOutput);

}

//______________________________________________________________________________
void AliJFJTask::UserExec(Option_t* /*option*/) 
{
    //cout << "======================== BEGIN EVENT ========================" << endl;
    // Processing of one event
    fhistos->fh_eventSel->Fill("events wo/ cuts",1.0);
    if(fDebug > 5) cout << "------- AliJFJTask Exec-------"<<endl;
    if(!((Entry()-1)%1000))  AliInfo(Form(" Processing event # %lld",  Entry())); 
    if( fJCatalystTask->GetJCatalystEntry() != fEntry) return;
    fhistos->fh_eventSel->Fill("catalyst entry ok",1.0);
    //if( !fJCatalystTask->GetIsGoodEvent() ) return;
    fhistos->fh_eventSel->Fill("catalyst ok",1.0);

    fCBin = AliJCDijetHistos::GetCentralityClass(fJCatalystTask->GetCentrality());
    fhistos->fh_centrality->Fill(fJCatalystTask->GetCentrality());
    if(fDebug ==4) cout << "fCbin = "<< fCBin << endl;

    if(fCBin == -1) return;

    fhistos->fh_eventSel->Fill("events",1.0);
    fhistos->fh_events[fCBin]->Fill("events",1.0);
    fhistos->fh_zvtx->Fill(fJCatalystTask->GetZVertex());
    
    TClonesArray *fInputList = (TClonesArray*)fJCatalystTask->GetInputList();

#if !defined(__CINT__) && !defined(__MAKECINT__)
    if(fDebug ==4) cout << "Next true level calculations:" << endl;
    fana->CalculateJets(fInputList, fhistos, fCBin);
    fana->FillJetsDijets(fhistos, fCBin);
#endif
}

//______________________________________________________________________________
void AliJFJTask::Init()
{
    // Intialisation of parameters
    AliInfo("Doing initialization") ; 

}

//______________________________________________________________________________
void AliJFJTask::Terminate(Option_t *)
{
    // Processing when the event loop is ended
    cout<<"AliJFJTask Analysis DONE !!"<<endl; 
}
