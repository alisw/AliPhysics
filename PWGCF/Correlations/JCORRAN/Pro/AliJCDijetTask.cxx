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
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/AreaDefinition.hh>
#include <fastjet/Selector.hh>
//#include <fastjet/FunctionOfPseudoJet.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
//#include <fastjet/tools/BackgroundEstimatorBase.hh>
#include <fastjet/tools/Subtractor.hh>

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
	
    TClonesArray *fInputList = (TClonesArray*)fJCatalystTask->GetInputList();
    CalculateJetsDijets(fInputList);
	
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
    double const ghost_maxrap = etaMaxCutForPart;
    unsigned int const repeat = 1; // default
    double const ghost_area   = 0.005; // ALICE=0.005 // default=0.01

    double phi, eta, pt, rho, rhom, area, mjj, ptpair, dPhi, dPhi2;
	vector<fastjet::PseudoJet> chparticles;
	vector<fastjet::PseudoJet> jets;
	vector<fastjet::PseudoJet> jets_bge;
	vector<fastjet::PseudoJet> jets_bgSubtracted;
	fastjet::PseudoJet jetAreaVector;
	fastjet::PseudoJet jet_bgSubtracted;
    fastjet::PseudoJet dijet;

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
            fhistos->fh_etaPhi[fCBin]->Fill(eta,phi);
            fhistos->fh_pt[fCBin]->Fill(pt);
            chparticles.push_back(fastjet::PseudoJet(trk->Px(), trk->Py(), trk->Pz(), trk->E()));
        }
    }

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Run the clustering, Reconstruct jets
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, coneR, fastjet::pt_scheme); //Other option: fastjet::E_scheme
    fastjet::JetDefinition jet_def_bge(fastjet::kt_algorithm, coneR, fastjet::pt_scheme); //Other option: fastjet::E_scheme

    fastjet::GhostedAreaSpec const area_spec(ghost_maxrap, repeat, ghost_area);
    fastjet::AreaDefinition const area_def(fastjet::active_area, area_spec);
    fastjet::AreaDefinition const area_def_bge(fastjet::active_area_explicit_ghosts, area_spec);

    // Selector selects first all jets inside rapidity acceptance and then all but two hardest jets.
    fastjet::Selector const selectorAllButTwo = (!fastjet::SelectorNHardest(2));
    fastjet::Selector const selectorEta = fastjet::SelectorAbsEtaMax(ghost_maxrap - coneR);
    fastjet::Selector const selectorBoth = selectorAllButTwo * selectorEta; // Here right selector is applied first, then the left one.
    fastjet::JetMedianBackgroundEstimator bge(selectorEta, jet_def_bge, area_def_bge);

    fastjet::ClusterSequenceArea cs(chparticles, jet_def, area_def);
    fastjet::ClusterSequenceArea cs_bge(chparticles, jet_def_bge, area_def_bge);
    
    jets = fastjet::sorted_by_pt(cs.inclusive_jets(MinJetPt)); // APPLY Min pt cut for jet
    jets_bge = fastjet::sorted_by_pt(cs_bge.inclusive_jets(MinJetPt)); // APPLY Min pt cut for jet
    //bge.set_particles(chparticles);
    //cout << "Testing: jets.size() = " << jets.size() << ", jets_bge.size() = " << jets_bge.size() << endl;
    //cout << "         selectorBoth(jets_bge).size() = " << selectorBoth(jets_bge).size() << ", selectorEta(jets_bge).size() = " << selectorEta(jets_bge).size() << endl;
    if( selectorBoth(jets_bge).size() < 1 ) {
        if( selectorEta(jets_bge).size() < 1) {
            rho  = 0.0;
            rhom = 0.0;
        } else {
            bge.set_jets(selectorEta(jets_bge));
            rho  = bge.rho()<0   ? 0.0 : bge.rho();
            rhom = bge.rho_m()<0 ? 0.0 : bge.rho_m();
        }
    } else { 
        bge.set_jets(selectorBoth(jets_bge));
        rho  = bge.rho()<0   ? 0.0 : bge.rho();
        rhom = bge.rho_m()<0 ? 0.0 : bge.rho_m();
    }
    

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Loop over jets and fill various histos 
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fhistos->fh_rho[fCBin]->Fill(rho);
    fhistos->fh_rhom[fCBin]->Fill(rhom);
    //std::cout << "Testing: Rho_M = " << rhom << ", has_rho_m() = " << bge.has_rho_m() << std::endl;

    for (unsigned ijet = 0; ijet < jets.size(); ijet++) {
        eta = jets[ijet].eta();
        // jet eta cut
        if(TMath::Abs(eta) < etaMaxCutForJet) {
            pt = jets[ijet].pt();
            phi = jets[ijet].phi();
            area = jets[ijet].area();
            jetAreaVector = jets[ijet].area_4vector();
            fhistos->fh_jetEta[fCBin]->Fill(eta);  
            fhistos->fh_jetPhi[fCBin]->Fill(phi - TMath::Pi()); //Pseudojet.phi range 0-2pi
            fhistos->fh_jetEtaPhi[fCBin]->Fill(eta,phi - TMath::Pi());
            fhistos->fh_jetPt[fCBin]->Fill(pt);
            fhistos->fh_jetArea[fCBin]->Fill(area);
            fhistos->fh_jetAreaRho[fCBin]->Fill(area*rho);
            
            jet_bgSubtracted = fastjet::PseudoJet(jets[ijet].px() -        rho * jetAreaVector.px(),
                                                  jets[ijet].py() -        rho * jetAreaVector.py(),
                                                  jets[ijet].pz() - (rho+rhom) * jetAreaVector.pz(),
                                                  jets[ijet].E()  - (rho+rhom) * jetAreaVector.E());

            // Check eta acceptance also for bg subtracted jets.
            eta = jet_bgSubtracted.eta();
            if(TMath::Abs(eta) < etaMaxCutForJet) {
                pt = jet_bgSubtracted.pt();
                phi = jet_bgSubtracted.phi();
                fhistos->fh_corrJetPt[fCBin]->Fill(pt);
                fhistos->fh_corrJetEta[fCBin]->Fill(eta);
                fhistos->fh_corrJetPhi[fCBin]->Fill(phi - TMath::Pi());

                jets_bgSubtracted.push_back(jet_bgSubtracted);
            }
        }
    }//end of the jet loop

    // Dijet calculations 
    if(jets.size()>1) {
        if(jets[1].pt()>MinDijetJetPt) {
            dijet = jets[0] + jets[1];
            mjj = dijet.m();
            ptpair = dijet.pt();
            fhistos->fh_dijetInvM[fCBin]->Fill(mjj);
            fhistos->fh_dijetPtPair[fCBin]->Fill(ptpair);
            dPhi = jets[1].delta_phi_to(jets[0]);
            dPhi2  = dPhi<0 ? dPhi+TMath::TwoPi() : dPhi;
            fhistos->fh_dijetDeltaPhi[fCBin]->Fill(dPhi2);

            // If subleading jet is on the opposite hemisphere compared to leading jet.
            if(TMath::Abs(dPhi2 - TMath::Pi()) < TMath::Pi()/2.0) {
                fhistos->fh_dijetInvMDeltaPhiCut[fCBin]->Fill(mjj); 
            }
        }
    }

    // Background subtracted dijets
    if(jets_bgSubtracted.size()>1) {
        if(jets_bgSubtracted[1].pt()>MinDijetJetPt) {
            dijet = jets_bgSubtracted[0] + jets_bgSubtracted[1];
            mjj = dijet.m();
            ptpair = dijet.pt();
            fhistos->fh_corrDijetInvM[fCBin]->Fill(mjj);
            fhistos->fh_corrDijetPtPair[fCBin]->Fill(ptpair);
            dPhi = jets_bgSubtracted[1].delta_phi_to(jets_bgSubtracted[0]);
            dPhi2  = dPhi<0 ? dPhi+TMath::TwoPi() : dPhi;
            fhistos->fh_corrDijetDeltaPhi[fCBin]->Fill(dPhi2);

            // If subleading jet is on the opposite hemisphere compared to leading jet.
            if(TMath::Abs(dPhi2 - TMath::Pi()) < TMath::Pi()/2.0) {
                fhistos->fh_corrDijetInvMDeltaPhiCut[fCBin]->Fill(mjj); 
            }
        }
    }
}
