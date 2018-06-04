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

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Analysis task for providing various flow informations
// author: O. Saarimaki, D.J. Kim(dong.jo.kim@cern.ch)
// ALICE Group University of Jyvaskyla
// Finland
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////
#include <AliAnalysisManager.h>
#include <AliJBaseTrack.h>
#include "AliJCDijetTask.h" 

// Fastjet includes
#include <fastjet/config.h>
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
// These will be added later
//#include <fastjet/ClusterSequence.hh>
//#include <fastjet/ClusterSequenceArea.hh>
//#include <fastjet/AreaDefinition.hh>
//#include <fastjet/Selector.hh>
//#include <fastjet/FunctionOfPseudoJet.hh>
//#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
//#include <fastjet/tools/BackgroundEstimatorBase.hh>
//#include <fastjet/tools/Subtractor.hh>

//______________________________________________________________________________
AliJCDijetTask::AliJCDijetTask() :   
	AliAnalysisTaskSE("JCDijetBaseTask"),
	fJCatalystTask(NULL),
	fJCatalystTaskName("JCatalystTask"),
	fIsMC(kTRUE),
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
	fIsMC(kTRUE),
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
	fIsMC(ap.fIsMC),
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
	fhistos->CreateEventTrackHistos();

	fhistos->fHMG->Print();


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
	if(fCBin == -1)
		return;
	if(fIsMC) {
		TClonesArray *fInputList = (TClonesArray*)fJCatalystTask->GetInputList();
		CalculateJetsDijets(fInputList);
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

//______________________________________________________________________________
void AliJCDijetTask::CalculateJetsDijets(TClonesArray *inList) {

	double const partMinPtCut = 0.15;// atlas 0.5 cms/alice 0.15
	double const coneR = 0.4; // atlas 0.6, cms 0.7 alice 0.4
	double const etaMaxCutForPart = 0.8;
	double const etaMaxCutForJet = etaMaxCutForPart-coneR;
	double const MinJetPt = 10.0; // Min Jet Pt cut to disregard low pt jets
	double const MinDijetJetPt = 20.0; // Min Jet Pt cut to disregard low pt jets

    double phi, eta, pt, mjj, ptpair, dPhi, dPhi2;
	vector<fastjet::PseudoJet> chparticles;

	//--------------------------------------------------------
	//         B e g i n    e v e n t    l o o p.
	//--------------------------------------------------------
	int noTracks = inList->GetEntries();

    chparticles.clear();
    for (int itrack = 0; itrack < noTracks; ++itrack) {//loop over all the particles in the event
        // Building input particle list for the jet reconstruction
		AliJBaseTrack *trk = (AliJBaseTrack*)inList->At(itrack);
		pt = trk->Pt();
		eta = trk->Eta();
        if (pt>partMinPtCut && TMath::Abs(eta) < etaMaxCutForPart){
            phi = trk->Phi();
            fhistos->fh_eta[fCBin]->Fill(eta);
            fhistos->fh_phi[fCBin]->Fill(phi);
            fhistos->fh_pt[fCBin]->Fill(pt);
            chparticles.push_back(fastjet::PseudoJet(trk->Px(), trk->Py(), trk->Pz(), trk->E()));
        }
    }

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Run the clustering, Reconstruct jets
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, coneR, fastjet::pt_scheme); //Other option: fastjet::E_scheme
    fastjet::ClusterSequence cs(chparticles, jet_def);
    vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets(MinJetPt)); // APPLY Min pt cut for jet

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Loop over jets and fill various histos 
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    for (unsigned ijet = 0; ijet < jets.size(); ijet++) {
        eta = jets[ijet].eta();
        // jet eta cut
        if(TMath::Abs(eta) < etaMaxCutForJet) {
            pt = jets[ijet].pt();
            phi = jets[ijet].phi();
            fhistos->fh_jetEta[fCBin]->Fill(eta);  
            fhistos->fh_jetPhi[fCBin]->Fill(phi - TMath::Pi()); //Pseudojet.phi range 0-2pi
            fhistos->fh_jetPt[fCBin]->Fill(pt);
        }
    }//end of the jet loop

    // Dijet 
    if(jets.size()>1) {
        if(jets[1].pt()>MinDijetJetPt) {
            fastjet::PseudoJet dijet = jets[0] + jets[1];
            mjj = dijet.m();
            ptpair = dijet.pt();
            fhistos->fh_DijetInvM[fCBin]->Fill(mjj);
            fhistos->fh_DijetPtPair[fCBin]->Fill(ptpair);
            dPhi = jets[1].delta_phi_to(jets[0]);
            dPhi2  = dPhi<0 ? dPhi+TMath::TwoPi() : dPhi;
            fhistos->fh_DijetDeltaPhi[fCBin]->Fill(dPhi2);

            // If subleading jet is on the opposite hemisphere compared to leading jet.
            if(TMath::Abs(dPhi2 - TMath::Pi()) < TMath::Pi()/2.0) {
                fhistos->fh_DijetInvMDeltaPhiCut[fCBin]->Fill(mjj); 
            }
        }
    }
}
