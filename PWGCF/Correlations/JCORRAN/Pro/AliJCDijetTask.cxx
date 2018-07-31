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
    fhistos->fh_info->Fill("MC", IsMC());
    for(unsigned i=0; i< fcentralityBins.size(); i++) fhistos->fh_info->Fill(Form("Cent bin border %02d",i), fcentralityBins.at(i));
    fhistos->fh_info->Fill("Jet cone", fjetCone);
    fhistos->fh_info->Fill("kt-jet cone", fktJetCone);
    fhistos->fh_info->Fill("Scheme", fktScheme);
    fhistos->fh_info->Fill("Use pion mass", fusePionMassInktjets);
    fhistos->fh_info->Fill("Particle eta cut", fparticleEtaCut);
    fhistos->fh_info->Fill("Particle pt cut", fparticlePtCut);
    fhistos->fh_info->Fill("Leading jet cut", fleadingJetCut);
    fhistos->fh_info->Fill("Subleading jet cut", fsubleadingJetCut);
    fhistos->fh_info->Fill("Const. cut", fconstituentCut);
    fhistos->fh_info->Fill("Delta phi cut pi/",fdeltaPhiCut);

    // Initialize fh_events so that the bin order is correct
    for (int iBin=0; iBin < fcentralityBins.size(); iBin++) {
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
    CalculateJetsDijets(fInputList,
                        fDebug,
                        fCBin,
                        fparticleEtaCut,
                        fparticlePtCut,
                        fjetCone,
                        fktJetCone,
                        fktScheme,
                        fusePionMassInktjets,
                        fconstituentCut,
                        fleadingJetCut,
                        fsubleadingJetCut,
                        fdeltaPhiCut);
	
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
void AliJCDijetTask::CalculateJetsDijets(TClonesArray *inList,
                                         int    lDebug,
                                         int    lCBin,
                                         double lParticleEtaCut,
                                         double lParticlePtCut,
                                         double lJetCone,
                                         double lktJetCone,
                                         int    lktScheme,
                                         bool   lusePionMassInkt,
                                         double lConstituentCut,
                                         double lLeadingJetCut,
                                         double lSubleadingJetCut,
                                         double lDeltaPhiCut) {

	double const etaMaxCutForJet = lParticleEtaCut-lJetCone;
	double const MinJetPt = 10.0; // Min Jet Pt cut to disregard low pt jets
    double const ghost_maxrap = lParticleEtaCut;
    unsigned int const repeat = 1; // default
    double const ghost_area   = 0.005; // ALICE=0.005 // default=0.01
    double const pionmass = AliPID::ParticleMass(AliPID::kPion);
    enum jetClasses {iRaw, iBGSubtr, iBGSubtrConstCut, iConstCut, iktJets, jetClassesSize};

    TString sDijetTypes[jetClassesSize] = {"raw", "bg. subtr.", "bg. subtr. const. cut", "const. cut", "kt"};

    double phi, eta, pt, pt2, rho, rhom, area, mjj, ptpair, dPhi, dPhi2;
    bool leadingTrackOverThreshold = false;
	vector<fastjet::PseudoJet> chparticles;
	vector<fastjet::PseudoJet> ktchparticles;
	vector<fastjet::PseudoJet> jets[jetClassesSize];
	vector<fastjet::PseudoJet> constituents;
    fastjet::RecombinationScheme ktScheme;
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
        fhistos->fh_events[lCBin]->Fill("particles",1.0);
        if (pt>lParticlePtCut && TMath::Abs(eta) < lParticleEtaCut){
            fhistos->fh_events[lCBin]->Fill("acc. particles",1.0);
            phi = trk->Phi();
            fhistos->fh_eta[lCBin]->Fill(eta);
            fhistos->fh_phi[lCBin]->Fill(phi);
            fhistos->fh_etaPhi[lCBin]->Fill(eta,phi);
            fhistos->fh_pt[lCBin]->Fill(pt);
            chparticles.push_back(fastjet::PseudoJet(trk->Px(), trk->Py(), trk->Pz(), trk->E()));
            if(lusePionMassInkt) ktchparticles.push_back(fastjet::PseudoJet(trk->Px(), trk->Py(), trk->Pz(), TMath::Sqrt(trk->Px()*trk->Px() + trk->Py()*trk->Py() + trk->Pz()*trk->Pz() + pionmass*pionmass)));            
            else ktchparticles.push_back(fastjet::PseudoJet(trk->Px(), trk->Py(), trk->Pz(), trk->E()));
        }
    }

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Run the clustering, Reconstruct jets
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    switch (lktScheme) {
        case 0:  ktScheme = fastjet::E_scheme;
                 break;
        case 1:  ktScheme = fastjet::pt_scheme;
                 break;
        case 2:  ktScheme = fastjet::pt2_scheme;
                 break;
        case 3:  ktScheme = fastjet::Et_scheme;
                 break;
        case 4:  ktScheme = fastjet::Et2_scheme;
                 break;
        case 5:  ktScheme = fastjet::BIpt_scheme;
                 break;
        case 6:  ktScheme = fastjet::BIpt2_scheme;
                 break;
        default: ::Error("AliJCDijetTask","Unknown recombination scheme!");
    }
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, lJetCone, fastjet::pt_scheme); //Other option: fastjet::E_scheme
    fastjet::JetDefinition jet_def_bge(fastjet::kt_algorithm, lktJetCone, ktScheme);

    fastjet::GhostedAreaSpec const area_spec(ghost_maxrap, repeat, ghost_area);
    fastjet::AreaDefinition const area_def(fastjet::active_area, area_spec);
    fastjet::AreaDefinition const area_def_bge(fastjet::active_area_explicit_ghosts, area_spec);

    // Selector selects first all jets inside rapidity acceptance and then all but two hardest jets.
    fastjet::Selector const selectorAllButTwo = (!fastjet::SelectorNHardest(2));
    fastjet::Selector const selectorEta = fastjet::SelectorAbsEtaMax(ghost_maxrap - lJetCone);
    fastjet::Selector const selectorBoth = selectorAllButTwo * selectorEta; // Here right selector is applied first, then the left one.
    fastjet::JetMedianBackgroundEstimator bge(selectorEta, jet_def_bge, area_def_bge);

    fastjet::ClusterSequenceArea cs(chparticles, jet_def, area_def);
    fastjet::ClusterSequenceArea cs_bge(ktchparticles, jet_def_bge, area_def_bge);
    
    jets[iRaw]    = fastjet::sorted_by_pt(cs.inclusive_jets(MinJetPt)); // APPLY Min pt cut for jet
    jets[iktJets] = fastjet::sorted_by_pt(cs_bge.inclusive_jets(0.0)); // APPLY Min pt cut for jet

    if( selectorBoth(jets[iktJets]).size() < 1 ) {
        fhistos->fh_events[lCBin]->Fill("no rho calc. events",1.0);
        rho  = 0.0;
        rhom = 0.0;
    } else { 
        fhistos->fh_events[lCBin]->Fill("rho calc. events",1.0);
        bge.set_jets(selectorBoth(jets[iktJets]));
        rho  = bge.rho()<0   ? 0.0 : bge.rho();
        rhom = bge.rho_m()<0 ? 0.0 : bge.rho_m();
    }
    

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Loop over jets and fill various histos 
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fhistos->fh_rho[lCBin]->Fill(rho);
    fhistos->fh_rhom[lCBin]->Fill(rhom);
    if(lDebug > 9) std::cout << "Testing: Rho_M = " << rhom << ", has_rho_m() = " << bge.has_rho_m() << std::endl;

    // anti-kt jets:
    for (unsigned ijet = 0; ijet < jets[iRaw].size(); ijet++) {
        eta = jets[iRaw][ijet].eta();
        fhistos->fh_events[lCBin]->Fill("jets",1.0);
        // anti-kt-jet eta cut
        if(TMath::Abs(eta) < etaMaxCutForJet) {
            fhistos->fh_events[lCBin]->Fill("acc. jets",1.0);
            pt = jets[iRaw][ijet].pt();
            phi = jets[iRaw][ijet].phi();
            area = jets[iRaw][ijet].area();
            jetAreaVector = jets[iRaw][ijet].area_4vector();
            fhistos->fh_jetEta[lCBin][iRaw]->Fill(eta);  
            fhistos->fh_jetPhi[lCBin][iRaw]->Fill(phi - TMath::Pi()); //Pseudojet.phi range 0-2pi
            fhistos->fh_jetEtaPhi[lCBin][iRaw]->Fill(eta,phi - TMath::Pi());
            fhistos->fh_jetPt[lCBin][iRaw]->Fill(pt);
            fhistos->fh_jetArea[lCBin][iRaw]->Fill(area);
            fhistos->fh_jetAreaRho[lCBin][iRaw]->Fill(area*rho);

            leadingTrackOverThreshold=false;
            if(lDebug > 9) cout << "Jet i=" << ijet << ", jet pt=" << pt << endl;
            for(unsigned iconst=0;iconst<jets[iRaw][ijet].constituents().size(); iconst++) {
                if(lDebug > 9) cout << "Constituent i=" << iconst << ", constituent pt=" << jets[iRaw][ijet].constituents()[iconst].pt() << endl;
                if(jets[iRaw][ijet].constituents()[iconst].pt() > lConstituentCut) { // Jet leading constituent cut.
                    leadingTrackOverThreshold=true;
                    break;
                }
            }
            
            jet_bgSubtracted = fastjet::PseudoJet(jets[iRaw][ijet].px() -        rho * jetAreaVector.px(),
                                                  jets[iRaw][ijet].py() -        rho * jetAreaVector.py(),
                                                  jets[iRaw][ijet].pz() - (rho+rhom) * jetAreaVector.pz(),
                                                  jets[iRaw][ijet].E()  - (rho+rhom) * jetAreaVector.E());

            if(leadingTrackOverThreshold) {
                fhistos->fh_events[lCBin]->Fill("const. cut jets",1.0);
                fhistos->fh_jetEta[lCBin][iConstCut]->Fill(eta);
                fhistos->fh_jetPhi[lCBin][iConstCut]->Fill(phi - TMath::Pi());
                fhistos->fh_jetEtaPhi[lCBin][iConstCut]->Fill(eta,phi - TMath::Pi());
                fhistos->fh_jetPt[lCBin][iConstCut]->Fill(pt);
                fhistos->fh_jetArea[lCBin][iConstCut]->Fill(area);
                fhistos->fh_jetAreaRho[lCBin][iConstCut]->Fill(area*rho);

                jets[iConstCut].push_back(jet_bgSubtracted);
            }


            // Check eta acceptance also for bg subtracted jets.
            eta = jet_bgSubtracted.eta();
            if(TMath::Abs(eta) < etaMaxCutForJet) {
                fhistos->fh_events[lCBin]->Fill("bg. subtr. jets",1.0);
                pt2 = jet_bgSubtracted.pt();
                phi = jet_bgSubtracted.phi();
                if(ijet==0 && pt>lLeadingJetCut && pt2<=lLeadingJetCut)       fhistos->fh_events[lCBin]->Fill("leading jet drop",1.0);
                if(ijet==1 && pt>lSubleadingJetCut && pt2<=lSubleadingJetCut) fhistos->fh_events[lCBin]->Fill("subleading jet drop",1.0);
                fhistos->fh_jetEta[lCBin][iBGSubtr]->Fill(eta);
                fhistos->fh_jetPhi[lCBin][iBGSubtr]->Fill(phi - TMath::Pi());
                fhistos->fh_jetEtaPhi[lCBin][iBGSubtr]->Fill(eta,phi - TMath::Pi());
                fhistos->fh_jetPt[lCBin][iBGSubtr]->Fill(pt2);
                fhistos->fh_jetArea[lCBin][iBGSubtr]->Fill(area); // Assuming bg subtracted jet has the same area.
                fhistos->fh_jetAreaRho[lCBin][iBGSubtr]->Fill(area*rho);

                jets[iBGSubtr].push_back(jet_bgSubtracted);

                if(leadingTrackOverThreshold) {
                    fhistos->fh_events[lCBin]->Fill("bg. subtr. const. cut jets",1.0);
                    fhistos->fh_jetEta[lCBin][iBGSubtrConstCut]->Fill(eta);
                    fhistos->fh_jetPhi[lCBin][iBGSubtrConstCut]->Fill(phi - TMath::Pi());
                    fhistos->fh_jetEtaPhi[lCBin][iBGSubtrConstCut]->Fill(eta,phi - TMath::Pi());
                    fhistos->fh_jetPt[lCBin][iBGSubtrConstCut]->Fill(pt2);
                    fhistos->fh_jetArea[lCBin][iBGSubtrConstCut]->Fill(area);
                    fhistos->fh_jetAreaRho[lCBin][iBGSubtrConstCut]->Fill(area*rho);

                    jets[iBGSubtrConstCut].push_back(jet_bgSubtracted);
                }
            }
        }
    }//end of the anti-kt-jet loop

    for (unsigned ijet = 0; ijet < jets[iktJets].size(); ijet++) {
        eta = jets[iktJets][ijet].eta();
        fhistos->fh_events[lCBin]->Fill("kt-jets",1.0);
        // kt-jet eta cut
        if(TMath::Abs(eta) < etaMaxCutForJet) {
            fhistos->fh_events[lCBin]->Fill("acc. kt-jets",1.0);
            pt = jets[iktJets][ijet].pt();
            phi = jets[iktJets][ijet].phi();
            area = jets[iktJets][ijet].area();
            fhistos->fh_jetEta[lCBin][iktJets]->Fill(eta);  
            fhistos->fh_jetPhi[lCBin][iktJets]->Fill(phi - TMath::Pi()); //Pseudojet.phi range 0-2pi
            fhistos->fh_jetEtaPhi[lCBin][iktJets]->Fill(eta,phi - TMath::Pi());
            fhistos->fh_jetPt[lCBin][iktJets]->Fill(pt);
            fhistos->fh_jetArea[lCBin][iktJets]->Fill(area);
            fhistos->fh_jetAreaRho[lCBin][iktJets]->Fill(area*rho);
        }
    } //end of the kt-jet loop



    // Dijet calculations 
    for(int idijet=0; idijet < jetClassesSize; idijet++) {
        if(jets[idijet].size()>1) {
            jets[idijet] = fastjet::sorted_by_pt(jets[idijet]); // Sort in case of bg subtr messed up the order.
            fhistos->fh_events[lCBin]->Fill(Form("%s dijets",sDijetTypes[idijet].Data()),1.0);
            if(jets[idijet][0].pt()>lLeadingJetCut) {
                fhistos->fh_events[lCBin]->Fill(Form("%s dijets leading cut",sDijetTypes[idijet].Data()),1.0);
                if(jets[idijet][1].pt()>lSubleadingJetCut) {
                    fhistos->fh_events[lCBin]->Fill(Form("%s acc. dijets",sDijetTypes[idijet].Data()),1.0);
                    dijet = jets[idijet][0] + jets[idijet][1];
                    mjj = dijet.m();
                    ptpair = dijet.pt();
                    fhistos->fh_dijetInvM[lCBin][idijet]->Fill(mjj);
                    fhistos->fh_dijetPtPair[lCBin][idijet]->Fill(ptpair);
                    dPhi = jets[idijet][1].delta_phi_to(jets[idijet][0]);
                    dPhi2  = dPhi<0 ? dPhi+TMath::TwoPi() : dPhi;
                    fhistos->fh_dijetDeltaPhi[lCBin][idijet]->Fill(dPhi2);

                    // If subleading jet is on the opposite hemisphere compared to leading jet.
                    if(TMath::Abs(dPhi2 - TMath::Pi()) < TMath::Pi()/lDeltaPhiCut) {
                        fhistos->fh_events[lCBin]->Fill(Form("%s deltaphi cut dijets",sDijetTypes[idijet].Data()),1.0);
                        fhistos->fh_dijetInvMDeltaPhiCut[lCBin][idijet]->Fill(mjj); 
                        fhistos->fh_dijetPtPairDeltaPhiCut[lCBin][idijet]->Fill(ptpair); 
                    }
                }
            }
        }
    }
}
