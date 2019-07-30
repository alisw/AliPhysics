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
    fDebug(0),
    fParticleEtaCut(0),
    fParticlePtCut(0),
    fusePionMassInkt(0),
    fUseDeltaPhiBGSubtr(0),
    fConstituentCut(0),
    fLeadingJetCut(0),
    fSubleadingJetCut(0),
    fDeltaPhiCut(0),
    etaMaxCutForJet(0),
    MinJetPt(0),
    pionmass(0)
#if !defined(__CINT__) && !defined(__MAKECINT__)
   ,chparticles(0),
    ktchparticles(0),
    jets(),
    rhoEstJets(),
    constituents(),
    ktScheme(),
    jetAreaVector(),
    jet_bgSubtracted(),
    dijet(),
    jet_def(),
    jet_def_bge(),
    area_spec(),
    area_def(),
    area_def_bge(),
    selectorAllButTwo(),
    selectorEta(),
    selectorBoth(),
    bge()
#endif
{
    // constructor
}

//----------------------------------------------------------------------------------------------------------------------------
AliJCDijetAna::~AliJCDijetAna(){
    // destructor
}

//----------------------------------------------------------------------------------------------------------------------------
AliJCDijetAna::AliJCDijetAna(const AliJCDijetAna& obj) :
    fDebug(obj.fDebug),
    fParticleEtaCut(obj.fParticleEtaCut),
    fParticlePtCut(obj.fParticlePtCut),
    fusePionMassInkt(obj.fusePionMassInkt),
    fUseDeltaPhiBGSubtr(obj.fUseDeltaPhiBGSubtr),
    fConstituentCut(obj.fConstituentCut),
    fLeadingJetCut(obj.fLeadingJetCut),
    fSubleadingJetCut(obj.fSubleadingJetCut),
    fDeltaPhiCut(obj.fDeltaPhiCut),
    etaMaxCutForJet(obj.etaMaxCutForJet),
    MinJetPt(obj.MinJetPt),
    pionmass(obj.pionmass)
#if !defined(__CINT__) && !defined(__MAKECINT__)
   ,chparticles(obj.chparticles),
    ktchparticles(obj.ktchparticles),
    jets(obj.jets),
    rhoEstJets(obj.rhoEstJets),
    constituents(obj.constituents),
    ktScheme(obj.ktScheme),
    jetAreaVector(obj.jetAreaVector),
    jet_bgSubtracted(obj.jet_bgSubtracted),
    dijet(obj.dijet),
    jet_def(obj.jet_def),
    jet_def_bge(obj.jet_def_bge),
    area_spec(obj.area_spec),
    area_def(obj.area_def),
    area_def_bge(obj.area_def_bge),
    selectorAllButTwo(obj.selectorAllButTwo),
    selectorEta(obj.selectorEta),
    selectorBoth(obj.selectorBoth),
    bge(obj.bge)
#endif
{
    // copy constructor
}

//----------------------------------------------------------------------------------------------------------------------------
AliJCDijetAna& AliJCDijetAna::operator=(const AliJCDijetAna& obj){
    //----------------------------------------------------------------------------------------------------------------------------
    // equal sign operator
    return *this;
}

#if !defined(__CINT__) && !defined(__MAKECINT__)
void AliJCDijetAna::SetSettings(int    lDebug,
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
                                double lDeltaPhiCut){
    fDebug = lDebug;
    fParticleEtaCut = lParticleEtaCut;
    fParticlePtCut = lParticlePtCut;
    fusePionMassInkt = lusePionMassInkt;
    fUseDeltaPhiBGSubtr = luseDeltaPhiBGSubtr;
    fConstituentCut = lConstituentCut;
    fLeadingJetCut = lLeadingJetCut;
    fSubleadingJetCut = lSubleadingJetCut;
    fDeltaPhiCut = lDeltaPhiCut;

    etaMaxCutForJet = lParticleEtaCut-lJetCone;
    MinJetPt = 10.0; // Min Jet Pt cut to disregard low pt jets
    double const ghost_maxrap = lParticleEtaCut;
    unsigned int const repeat = 1; // default
    double const ghost_area   = 0.005; // ALICE=0.005 // default=0.01
    pionmass = 0.139570;//AliPID::ParticleMass(AliPID::kPion);

    //Initialize jet lists
    vector<fastjet::PseudoJet> emptyJetVector;
    for(int i=0; i<jetClassesSize; i++) jets.push_back(emptyJetVector);

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

    //Jet definitions and other fastjet settings
    jet_def = fastjet::JetDefinition(fastjet::antikt_algorithm, lJetCone, fastjet::pt_scheme); //Other option: fastjet::E_scheme
    jet_def_bge = fastjet::JetDefinition(fastjet::kt_algorithm, lktJetCone, ktScheme);

    area_spec = fastjet::GhostedAreaSpec(ghost_maxrap, repeat, ghost_area);
    area_def = fastjet::AreaDefinition(fastjet::active_area, area_spec);
    area_def_bge = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, area_spec);

    // Selector selects first all jets inside rapidity acceptance and then all but two hardest jets.
    selectorAllButTwo  = fastjet::Selector(!fastjet::SelectorNHardest(2));
    selectorEta  = fastjet::SelectorAbsEtaMax(ghost_maxrap - lktJetCone);
    selectorBoth  = selectorAllButTwo * selectorEta; // Here right selector is applied first, then the left one.
    bge = fastjet::JetMedianBackgroundEstimator(selectorEta, jet_def_bge, area_def_bge);


}


//______________________________________________________________________________
void AliJCDijetAna::CalculateJetsDijets(TClonesArray *inList, AliJCDijetHistos *fhistos, int lCBin){

    chparticles.clear();
    ktchparticles.clear();
    for(int i=0; i<jetClassesSize; i++) jets[i].clear();
    rhoEstJets.clear();
    constituents.clear();

    TString sDijetTypes[jetClassesSize] = {"raw", "bg. subtr.", "bg. subtr. const. cut", "const. cut", "kt"};
    phi=0, eta=0, pt=0, pt2=0, rho=0, rhom=0, area=0, mjj=0, ptpair=0, dPhi=0, dPhi2=0;
    leadingTrackOverThreshold = false;

    //--------------------------------------------------------
    //         B e g i n    e v e n t    l o o p.
    //--------------------------------------------------------
    if (inList==0) { cout << "No list!" << endl; return;}
    noTracks = inList->GetEntriesFast();

    for (utrack = 0; utrack < noTracks; ++utrack) {//loop over all the particles in the event
        // Building input particle list for the jet reconstruction
        AliJBaseTrack *trk = (AliJBaseTrack*)inList->At(utrack);
        pt = trk->Pt();
        eta = trk->Eta();
        fhistos->fh_events[lCBin]->Fill("particles",1.0);
        if (pt>fParticlePtCut && TMath::Abs(eta) < fParticleEtaCut){
            fhistos->fh_events[lCBin]->Fill("acc. particles",1.0);
            phi = trk->Phi();
            fhistos->fh_eta[lCBin]->Fill(eta);
            fhistos->fh_phi[lCBin]->Fill(phi);
            fhistos->fh_etaPhi[lCBin]->Fill(eta,phi);
            fhistos->fh_pt[lCBin]->Fill(pt);
            chparticles.push_back(fastjet::PseudoJet(trk->Px(), trk->Py(), trk->Pz(), trk->E()));
            if(fusePionMassInkt) ktchparticles.push_back(fastjet::PseudoJet(trk->Px(), trk->Py(), trk->Pz(), TMath::Sqrt(trk->Px()*trk->Px() + trk->Py()*trk->Py() + trk->Pz()*trk->Pz() + pionmass*pionmass)));            
            else ktchparticles.push_back(fastjet::PseudoJet(trk->Px(), trk->Py(), trk->Pz(), trk->E()));
        }
    }
    if(chparticles.size()==0) return; // We are not intereted in empty events.

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Run the clustering, Reconstruct jets
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fastjet::ClusterSequenceArea cs(chparticles, jet_def, area_def);
    fastjet::ClusterSequenceArea cs_bge(ktchparticles, jet_def_bge, area_def_bge);

    jets.at(iRaw)    = fastjet::sorted_by_pt(cs.inclusive_jets(MinJetPt)); // APPLY Min pt cut for jet
    jets.at(iktJets) = fastjet::sorted_by_pt(cs_bge.inclusive_jets(0.0)); // APPLY Min pt cut for jet

    if(fUseDeltaPhiBGSubtr) {
        removed = false;
        for (uktjet = 1; uktjet < jets[iktJets].size(); uktjet++) { // First jet is already skipped here.
            if (!removed && TMath::Abs(jets[iktJets][uktjet].delta_phi_to(jets[iktJets][0])) > TMath::Pi()/fDeltaPhiCut) {
                removed = true;
                continue;
            }
            rhoEstJets.push_back(jets[iktJets][uktjet]); 
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
    if(fDebug > 9) std::cout << "Testing: Rho_M = " << rhom << ", has_rho_m() = " << bge.has_rho_m() << std::endl;

    // anti-kt jets:
    for (ujet = 0; ujet < jets[iRaw].size(); ujet++) {
        eta = jets[iRaw][ujet].eta();
        fhistos->fh_events[lCBin]->Fill("jets",1.0);
        // anti-kt-jet eta cut
        if(TMath::Abs(eta) < etaMaxCutForJet) {
            fhistos->fh_events[lCBin]->Fill("acc. jets",1.0);
            pt = jets[iRaw][ujet].pt();
            phi = jets[iRaw][ujet].phi();
            area = jets[iRaw][ujet].area();
            jetAreaVector = jets[iRaw][ujet].area_4vector();
            fhistos->fh_jetEta[lCBin][iRaw]->Fill(eta);  
            fhistos->fh_jetPhi[lCBin][iRaw]->Fill(phi - TMath::Pi()); //Pseudojet.phi range 0-2pi
            fhistos->fh_jetEtaPhi[lCBin][iRaw]->Fill(eta,phi - TMath::Pi());
            fhistos->fh_jetPt[lCBin][iRaw]->Fill(pt);
            fhistos->fh_jetArea[lCBin][iRaw]->Fill(area);
            fhistos->fh_jetAreaRho[lCBin][iRaw]->Fill(area*rho);

            leadingTrackOverThreshold=false;
            if(fDebug > 9) cout << "Jet i=" << ujet << ", jet pt=" << pt << endl;
            for(uconst=0;uconst<jets[iRaw][ujet].constituents().size(); uconst++) {
                if(fDebug > 9) cout << "Constituent i=" << uconst << ", constituent pt=" << jets[iRaw][ujet].constituents()[uconst].pt() << endl;
                if(jets[iRaw][ujet].constituents()[uconst].pt() > fConstituentCut) { // Jet leading constituent cut.
                    leadingTrackOverThreshold=true;
                    break;
                }
            }

            jet_bgSubtracted = fastjet::PseudoJet(jets[iRaw][ujet].px() -        rho * jetAreaVector.px(),
                    jets[iRaw][ujet].py() -        rho * jetAreaVector.py(),
                    jets[iRaw][ujet].pz() - (rho+rhom) * jetAreaVector.pz(),
                    jets[iRaw][ujet].E()  - (rho+rhom) * jetAreaVector.E());

            if(leadingTrackOverThreshold) {
                fhistos->fh_events[lCBin]->Fill("const. cut jets",1.0);
                fhistos->fh_jetEta[lCBin][iConstCut]->Fill(eta);
                fhistos->fh_jetPhi[lCBin][iConstCut]->Fill(phi - TMath::Pi());
                fhistos->fh_jetEtaPhi[lCBin][iConstCut]->Fill(eta,phi - TMath::Pi());
                fhistos->fh_jetPt[lCBin][iConstCut]->Fill(pt);
                fhistos->fh_jetArea[lCBin][iConstCut]->Fill(area);
                fhistos->fh_jetAreaRho[lCBin][iConstCut]->Fill(area*rho);

                jets[iConstCut].push_back(jets[iRaw][ujet]);
            }


            // Check eta acceptance also for bg subtracted jets.
            eta = jet_bgSubtracted.eta();
            if(TMath::Abs(eta) < etaMaxCutForJet) {
                fhistos->fh_events[lCBin]->Fill("bg. subtr. jets",1.0);
                pt2 = jet_bgSubtracted.pt();
                phi = jet_bgSubtracted.phi();
                if(ujet==0 && pt>fLeadingJetCut && pt2<=fLeadingJetCut)       fhistos->fh_events[lCBin]->Fill("leading jet drop",1.0);
                if(ujet==1 && pt>fSubleadingJetCut && pt2<=fSubleadingJetCut) fhistos->fh_events[lCBin]->Fill("subleading jet drop",1.0);
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

    for (uktjet = 0; uktjet < jets[iktJets].size(); uktjet++) {
        eta = jets[iktJets][uktjet].eta();
        fhistos->fh_events[lCBin]->Fill("kt-jets",1.0);
        // kt-jet eta cut
        if(TMath::Abs(eta) < etaMaxCutForJet) {
            fhistos->fh_events[lCBin]->Fill("acc. kt-jets",1.0);
            pt = jets[iktJets][uktjet].pt();
            phi = jets[iktJets][uktjet].phi();
            area = jets[iktJets][uktjet].area();
            fhistos->fh_jetEta[lCBin][iktJets]->Fill(eta);  
            fhistos->fh_jetPhi[lCBin][iktJets]->Fill(phi - TMath::Pi()); //Pseudojet.phi range 0-2pi
            fhistos->fh_jetEtaPhi[lCBin][iktJets]->Fill(eta,phi - TMath::Pi());
            fhistos->fh_jetPt[lCBin][iktJets]->Fill(pt);
            fhistos->fh_jetArea[lCBin][iktJets]->Fill(area);
            fhistos->fh_jetAreaRho[lCBin][iktJets]->Fill(area*rho);
        }
    } //end of the kt-jet loop



    // Dijet calculations 
    for(udijet=0; udijet < jetClassesSize; udijet++) {
        if(jets[udijet].size()>1) {
            jets[udijet] = fastjet::sorted_by_pt(jets[udijet]); // Sort in case of bg subtr messed up the order.
            fhistos->fh_events[lCBin]->Fill(Form("%s dijets",sDijetTypes[udijet].Data()),1.0);
            if(jets[udijet][0].pt()>fLeadingJetCut) {
                fhistos->fh_events[lCBin]->Fill(Form("%s dijets leading cut",sDijetTypes[udijet].Data()),1.0);
                if(jets[udijet][1].pt()>fSubleadingJetCut) {
                    fhistos->fh_events[lCBin]->Fill(Form("%s acc. dijets",sDijetTypes[udijet].Data()),1.0);
                    dijet = jets[udijet][0] + jets[udijet][1];
                    mjj = dijet.m();
                    ptpair = dijet.pt();
                    fhistos->fh_dijetInvM[lCBin][udijet]->Fill(mjj);
                    fhistos->fh_dijetPtPair[lCBin][udijet]->Fill(ptpair);
                    dPhi = jets[udijet][1].delta_phi_to(jets[udijet][0]);
                    dPhi2  = dPhi<0 ? dPhi+TMath::TwoPi() : dPhi;
                    fhistos->fh_dijetDeltaPhi[lCBin][udijet]->Fill(dPhi2);

                    // If subleading jet is on the opposite hemisphere compared to leading jet.
                    if(TMath::Abs(dPhi2 - TMath::Pi()) < TMath::Pi()/fDeltaPhiCut) {
                        fhistos->fh_events[lCBin]->Fill(Form("%s deltaphi cut dijets",sDijetTypes[udijet].Data()),1.0);
                        fhistos->fh_dijetInvMDeltaPhiCut[lCBin][udijet]->Fill(mjj); 
                        fhistos->fh_dijetPtPairDeltaPhiCut[lCBin][udijet]->Fill(ptpair); 
                    }
                }
            }
        }
    }
}
#endif
