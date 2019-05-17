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
	fJCatalystTaskName("JCatalystTask"),
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
	fCBin(-1),
	fOutput(NULL)
{
}

//______________________________________________________________________________
AliJCDijetTask::AliJCDijetTask(const char *name, TString inputformat):
	AliAnalysisTaskSE(name), 
	fJCatalystTask(NULL),
	fJCatalystTaskName("JCatalystTask"),
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
	fCBin(-1),
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
	fJCatalystTaskName(ap.fJCatalystTaskName),
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
	fCBin(ap.fCBin),
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

}

//________________________________________________________________________

void AliJCDijetTask::UserCreateOutputObjects()
{  
	//=== create the jcorran outputs objects
	if(fDebug > 1) printf("AliJCDijetTask::UserCreateOutPutData() \n");
	//=== Get AnalysisManager
	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();

	fJCatalystTask = (AliJCatalystTask*)(man->GetTask( fJCatalystTaskName ));

	OpenFile(1);
	fOutput = gDirectory;
	fOutput->cd();


	fhistos = new AliJCDijetHistos();
	fhistos->SetCentralityBinsHistos(fcentralityBins);
	fhistos->CreateEventTrackHistos();
    
	fhistos->fHMG->Print();

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
    fhistos->fh_info->Fill("Count", 1.0);
    fhistos->fh_info->Fill("MC", IsMC());
    for(unsigned i=0; i< fcentralityBins.size(); i++) fhistos->fh_info->Fill(Form("Cent bin border %02d",i), fcentralityBins.at(i));
    fhistos->fh_info->Fill("Jet cone", fjetCone);
    fhistos->fh_info->Fill("kt-jet cone", fktJetCone);
    fhistos->fh_info->Fill("Scheme", fktScheme);
    fhistos->fh_info->Fill("Use pion mass", fusePionMassInktjets);
    fhistos->fh_info->Fill("Use DeltaPhi BG Subtr", fuseDeltaPhiBGSubtr);
    fhistos->fh_info->Fill("Particle eta cut", fparticleEtaCut);
    fhistos->fh_info->Fill("Particle pt cut", fparticlePtCut);
    fhistos->fh_info->Fill("Leading jet cut", fleadingJetCut);
    fhistos->fh_info->Fill("Subleading jet cut", fsubleadingJetCut);
    fhistos->fh_info->Fill("Const. cut", fconstituentCut);
    fhistos->fh_info->Fill("Delta phi cut pi/",fdeltaPhiCut);

    // Initialize fh_events so that the bin order is correct
    for (unsigned iBin=0; iBin < fcentralityBins.size()-1; iBin++) {
        fhistos->fh_events[iBin]->Fill("events",0.0);
        fhistos->fh_events[iBin]->Fill("particles",0.0);
        fhistos->fh_events[iBin]->Fill("acc. particles",0.0);
        fhistos->fh_events[iBin]->Fill("no rho calc. events",0.0);
        fhistos->fh_events[iBin]->Fill("rho calc. events",0.0);
        fhistos->fh_events[iBin]->Fill("jets",0.0);
        fhistos->fh_events[iBin]->Fill("acc. jets",0.0);
        fhistos->fh_events[iBin]->Fill("const. cut jets",0.0);
        fhistos->fh_events[iBin]->Fill("bg. subtr. jets",0.0);
        fhistos->fh_events[iBin]->Fill("bg. subtr. const. cut jets",0.0);
        fhistos->fh_events[iBin]->Fill("kt-jets",0.0);
        fhistos->fh_events[iBin]->Fill("acc. kt-jets",0.0);
        fhistos->fh_events[iBin]->Fill("leading jet drop",0.0);
        fhistos->fh_events[iBin]->Fill("subleading jet drop",0.0);
        fhistos->fh_events[iBin]->Fill("raw dijets",0.0);
        fhistos->fh_events[iBin]->Fill("raw dijets leading cut",0.0);
        fhistos->fh_events[iBin]->Fill("raw acc. dijets",0.0);
        fhistos->fh_events[iBin]->Fill("raw deltaphi cut dijets",0.0);
        fhistos->fh_events[iBin]->Fill("bg. subtr. dijets",0.0);
        fhistos->fh_events[iBin]->Fill("bg. subtr. dijets leading cut",0.0);
        fhistos->fh_events[iBin]->Fill("bg. subtr. acc. dijets",0.0);
        fhistos->fh_events[iBin]->Fill("bg. subtr. deltaphi cut dijets",0.0);
        fhistos->fh_events[iBin]->Fill("bg. subtr. const. cut dijets",0.0);
        fhistos->fh_events[iBin]->Fill("bg. subtr. const. cut dijets leading cut",0.0);
        fhistos->fh_events[iBin]->Fill("bg. subtr. const. cut acc. dijets",0.0);
        fhistos->fh_events[iBin]->Fill("bg. subtr. const. cut deltaphi cut dijets",0.0);
        fhistos->fh_events[iBin]->Fill("const. cut dijets",0.0);
        fhistos->fh_events[iBin]->Fill("const. cut dijets leading cut",0.0);
        fhistos->fh_events[iBin]->Fill("const. cut acc. dijets",0.0);
        fhistos->fh_events[iBin]->Fill("const. cut deltaphi cut dijets",0.0);
        fhistos->fh_events[iBin]->Fill("kt dijets",0.0);
        fhistos->fh_events[iBin]->Fill("kt dijets leading cut",0.0);
        fhistos->fh_events[iBin]->Fill("kt acc. dijets",0.0);
        fhistos->fh_events[iBin]->Fill("kt deltaphi cut dijets",0.0);
    }

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
	if( fJCatalystTask->GetJCatalystEntry() != fEntry ) return;
	fCBin = AliJCDijetHistos::GetCentralityClass(fJCatalystTask->GetCentrality());
	if(fCBin == -1) return;

    fhistos->fh_events[fCBin]->Fill("events",1.0);
    fhistos->fh_centrality->Fill(fJCatalystTask->GetCentrality());
    fhistos->fh_zvtx->Fill(fJCatalystTask->GetZVertex());
	
    TClonesArray *fInputList = (TClonesArray*)fJCatalystTask->GetInputList();

    AliJCDijetAna::CalculateJetsDijets(fInputList,
                        fDebug,
                        fCBin,
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
                        fdeltaPhiCut,
                        fhistos);
	
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

