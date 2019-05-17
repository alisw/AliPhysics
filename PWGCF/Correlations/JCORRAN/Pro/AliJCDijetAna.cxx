/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
//______________________________________________________________________________
// Analysis task for providing various dijet informations
// author: O. Saarimaki, D.J. Kim(dong.jo.kim@cern.ch)
// ALICE Group University of Jyvaskyla
// Finland
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////

// used in local and grid execution

#include <TH1D.h>
#include "AliJCDijetAna.h"

ClassImp(AliJCDijetAna)

    //----------------------------------------------------------------------------------------------------------------------------
    AliJCDijetAna::AliJCDijetAna() :
        TObject()
{
    // constructor
}

//----------------------------------------------------------------------------------------------------------------------------
AliJCDijetAna::~AliJCDijetAna(){
    // destructor
}

//----------------------------------------------------------------------------------------------------------------------------
AliJCDijetAna::AliJCDijetAna(const AliJCDijetAna& obj) :
    TObject()
{
    // copy constructor
}

//----------------------------------------------------------------------------------------------------------------------------
AliJCDijetAna& AliJCDijetAna::operator=(const AliJCDijetAna& obj){
    //----------------------------------------------------------------------------------------------------------------------------
    // equal sign operator
    return *this;
}

//______________________________________________________________________________
void AliJCDijetAna::CalculateJetsDijets(TClonesArray *inList,
                                        int    lDebug,
                                        int    lCBin,
                                        double lParticleEtaCut,
                                        double lParticlePtCut,
                                        double lJetCone,
                                        double lktJetCone,
                                        int    lktScheme,
                                        bool   lusePionMassInkt,
                                        bool   luseDeltaPhiBGSubtr,
                                        double lConstituentCut,
                                        double lLeadingJetCut,
                                        double lSubleadingJetCut,
                                        double lDeltaPhiCut,
                                        AliJCDijetHistos *fhistos){

    double const etaMaxCutForJet = lParticleEtaCut-lJetCone;
    double const MinJetPt = 10.0; // Min Jet Pt cut to disregard low pt jets
    double const ghost_maxrap = lParticleEtaCut;
    unsigned int const repeat = 1; // default
    double const ghost_area   = 0.005; // ALICE=0.005 // default=0.01
    double const pionmass = 0.139570;//AliPID::ParticleMass(AliPID::kPion);
    enum jetClasses {iRaw, iBGSubtr, iBGSubtrConstCut, iConstCut, iktJets, jetClassesSize};

    TString sDijetTypes[jetClassesSize] = {"raw", "bg. subtr.", "bg. subtr. const. cut", "const. cut", "kt"};

    double phi, eta, pt, pt2, rho, rhom, area, mjj, ptpair, dPhi, dPhi2;
    bool leadingTrackOverThreshold = false;
    vector<fastjet::PseudoJet> chparticles;
    vector<fastjet::PseudoJet> ktchparticles;
    vector<fastjet::PseudoJet> jets[jetClassesSize];
    vector<fastjet::PseudoJet> rhoEstJets;
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
    if(chparticles.size()==0) return; // We are not intereted in empty events.

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
        default: ktScheme = fastjet::external_scheme;
                 ::Error("AliJCDijetAna","Unknown recombination scheme!");
    }
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, lJetCone, fastjet::pt_scheme); //Other option: fastjet::E_scheme
    fastjet::JetDefinition jet_def_bge(fastjet::kt_algorithm, lktJetCone, ktScheme);

    fastjet::GhostedAreaSpec const area_spec(ghost_maxrap, repeat, ghost_area);
    fastjet::AreaDefinition const area_def(fastjet::active_area, area_spec);
    fastjet::AreaDefinition const area_def_bge(fastjet::active_area_explicit_ghosts, area_spec);

    // Selector selects first all jets inside rapidity acceptance and then all but two hardest jets.
    fastjet::Selector const selectorAllButTwo = (!fastjet::SelectorNHardest(2));
    fastjet::Selector const selectorEta = fastjet::SelectorAbsEtaMax(ghost_maxrap - lktJetCone);
    fastjet::Selector const selectorBoth = selectorAllButTwo * selectorEta; // Here right selector is applied first, then the left one.
    fastjet::JetMedianBackgroundEstimator bge(selectorEta, jet_def_bge, area_def_bge);

    fastjet::ClusterSequenceArea cs(chparticles, jet_def, area_def);
    fastjet::ClusterSequenceArea cs_bge(ktchparticles, jet_def_bge, area_def_bge);

    jets[iRaw]    = fastjet::sorted_by_pt(cs.inclusive_jets(MinJetPt)); // APPLY Min pt cut for jet
    jets[iktJets] = fastjet::sorted_by_pt(cs_bge.inclusive_jets(0.0)); // APPLY Min pt cut for jet

    if(luseDeltaPhiBGSubtr) {
        bool removed = false;
        for (unsigned iktJet = 1; iktJet < jets[iktJets].size(); iktJet++) { // First jet is already skipped here.
            if (!removed && TMath::Abs(jets[iktJets][iktJet].delta_phi_to(jets[iktJets][0])) > TMath::Pi()/lDeltaPhiCut) {
                removed = true;
                continue;
            }
            rhoEstJets.push_back(jets[iktJets][iktJet]); 
        }
    } else {
        rhoEstJets = selectorBoth(jets[iktJets]);
    }

    if( rhoEstJets.size() < 1 ) {
        fhistos->fh_events[lCBin]->Fill("no rho calc. events",1.0);
        rho  = 0.0;
        rhom = 0.0;
    } else { 
        fhistos->fh_events[lCBin]->Fill("rho calc. events",1.0);
        bge.set_jets(rhoEstJets);
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

                jets[iConstCut].push_back(jets[iRaw][ijet]);
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
