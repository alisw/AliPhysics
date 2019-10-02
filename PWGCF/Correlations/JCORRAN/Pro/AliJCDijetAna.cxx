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
    fusePionMass(0),
    fUseDeltaPhiBGSubtr(0),
    fConstituentCut(0),
    fLeadingJetCut(0),
    fSubleadingJetCut(0),
    fDeltaPhiCut(0),
    etaMaxCutForJet(0),
    etaMaxCutForKtJet(0),
    MinJetPt(0),
    pionmass(0),
    matchingR(0)
#if !defined(__CINT__) && !defined(__MAKECINT__)
   ,chparticles(0),
    ktchparticles(0),
    jets(),
    rawJets(),
    rawKtJets(),
    rhoEstJets(),
    constituents(),
    dijets(),
    ktScheme(),
    antiktScheme(),
    jetAreaVector(),
    jet_bgSubtracted(),
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
    fusePionMass(obj.fusePionMass),
    fUseDeltaPhiBGSubtr(obj.fUseDeltaPhiBGSubtr),
    fConstituentCut(obj.fConstituentCut),
    fLeadingJetCut(obj.fLeadingJetCut),
    fSubleadingJetCut(obj.fSubleadingJetCut),
    fDeltaPhiCut(obj.fDeltaPhiCut),
    etaMaxCutForJet(obj.etaMaxCutForJet),
    etaMaxCutForKtJet(obj.etaMaxCutForKtJet),
    MinJetPt(obj.MinJetPt),
    pionmass(obj.pionmass),
    matchingR(obj.matchingR)
#if !defined(__CINT__) && !defined(__MAKECINT__)
   ,chparticles(obj.chparticles),
    ktchparticles(obj.ktchparticles),
    jets(obj.jets),
    rawJets(obj.rawJets),
    rawKtJets(obj.rawKtJets),
    rhoEstJets(obj.rhoEstJets),
    constituents(obj.constituents),
    dijets(obj.dijets),
    ktScheme(obj.ktScheme),
    antiktScheme(obj.antiktScheme),
    jetAreaVector(obj.jetAreaVector),
    jet_bgSubtracted(obj.jet_bgSubtracted),
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
                                int    lantiktScheme,
                                bool   lusePionMass,
                                bool   luseDeltaPhiBGSubtr,
                                double lConstituentCut,
                                double lLeadingJetCut,
                                double lSubleadingJetCut,
                                double lDeltaPhiCut,
                                double lmatchingR){
    fDebug = lDebug;
    fParticleEtaCut = lParticleEtaCut;
    fParticlePtCut = lParticlePtCut;
    fusePionMass = lusePionMass;
    fUseDeltaPhiBGSubtr = luseDeltaPhiBGSubtr;
    fConstituentCut = lConstituentCut;
    fLeadingJetCut = lLeadingJetCut;
    fSubleadingJetCut = lSubleadingJetCut;
    fDeltaPhiCut = lDeltaPhiCut;

    etaMaxCutForJet = lParticleEtaCut-lJetCone;
    etaMaxCutForKtJet = lParticleEtaCut-lktJetCone;
    MinJetPt = 10.0; // Min Jet Pt cut to disregard low pt jets
    double const ghost_maxrap = lParticleEtaCut;
    unsigned int const repeat = 1; // default
    double const ghost_area   = 0.005; // ALICE=0.005 // default=0.01
    fJetCone = lJetCone;
    fktJetCone = lktJetCone;
    pionmass = 0.139570;//AliPID::ParticleMass(AliPID::kPion);
    matchingR = lmatchingR;

    //Initialize jet lists
    fastjet::PseudoJet emptyJet;
    vector<fastjet::PseudoJet> emptyJetVector;
    vector<vector<fastjet::PseudoJet>> emptyDijetVector;
    // Fill list for jets and two jets for dijets and
    // dijets with deltaPhi cut.
    for(int i=0; i<jetClassesSize; i++) {
        jets.push_back(emptyJetVector);
        emptyJetVector.push_back(emptyJet);
        emptyJetVector.push_back(emptyJet);
        emptyDijetVector.push_back(emptyJetVector);
        emptyDijetVector.push_back(emptyJetVector);
        dijets.push_back(emptyDijetVector);
    }
    bHasDijet = false;
    bHasDeltaPhiDijet = false;
    bHasDeltaPhiSubLeadJet = false;

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
    switch (lantiktScheme) {
        case 0:  antiktScheme = fastjet::E_scheme;
                 break;
        case 1:  antiktScheme = fastjet::pt_scheme;
                 break;
        case 2:  antiktScheme = fastjet::pt2_scheme;
                 break;
        case 3:  antiktScheme = fastjet::Et_scheme;
                 break;
        case 4:  antiktScheme = fastjet::Et2_scheme;
                 break;
        case 5:  antiktScheme = fastjet::BIpt_scheme;
                 break;
        case 6:  antiktScheme = fastjet::BIpt2_scheme;
                 break;
        default: antiktScheme = fastjet::external_scheme;
                 ::Error("AliJCDijetAna","Unknown recombination scheme!");
    }

    //Jet definitions and other fastjet settings
    jet_def = fastjet::JetDefinition(fastjet::antikt_algorithm, fJetCone, antiktScheme); //Other option: fastjet::E_scheme
    jet_def_bge = fastjet::JetDefinition(fastjet::kt_algorithm, fktJetCone, ktScheme);

    area_spec = fastjet::GhostedAreaSpec(ghost_maxrap, repeat, ghost_area);
    area_def = fastjet::AreaDefinition(fastjet::active_area, area_spec);
    area_def_bge = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, area_spec);

    // Selector selects first all jets inside rapidity acceptance and then all but two hardest jets.
    selectorAllButTwo  = fastjet::Selector(!fastjet::SelectorNHardest(2));
    selectorEta  = fastjet::SelectorAbsEtaMax(ghost_maxrap - fktJetCone);
    selectorBoth  = selectorAllButTwo * selectorEta; // Here right selector is applied first, then the left one.
    bge = fastjet::JetMedianBackgroundEstimator(selectorEta, jet_def_bge, area_def_bge);


}

//______________________________________________________________________________
void AliJCDijetAna::CalculateJets(TClonesArray *inList, AliJCDijetHistos *fhistos, int lCBin){

    ResetObjects();

    //--------------------------------------------------------
    //         B e g i n    e v e n t    l o o p.
    //--------------------------------------------------------
    if (inList==0) { cout << "No list!" << endl; return;}
    noTracks = inList->GetEntriesFast();
    //cout << "Number of tracks: " << noTracks << endl;

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
            if(fusePionMass) {
                chparticles.push_back(fastjet::PseudoJet(trk->Px(), trk->Py(), trk->Pz(), TMath::Sqrt(trk->Px()*trk->Px() + trk->Py()*trk->Py() + trk->Pz()*trk->Pz() + pionmass*pionmass)));
                ktchparticles.push_back(fastjet::PseudoJet(trk->Px(), trk->Py(), trk->Pz(), TMath::Sqrt(trk->Px()*trk->Px() + trk->Py()*trk->Py() + trk->Pz()*trk->Pz() + pionmass*pionmass)));
                }
            else {
                chparticles.push_back(fastjet::PseudoJet(trk->Px(), trk->Py(), trk->Pz(), trk->E()));
                ktchparticles.push_back(fastjet::PseudoJet(trk->Px(), trk->Py(), trk->Pz(), trk->E()));
            }
        }
    }
    if(chparticles.size()==0) return; // We are not intereted in empty events.

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Run the clustering, Reconstruct jets
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    cs.reset(new fastjet::ClusterSequenceArea(chparticles, jet_def, area_def));
    cs_bge.reset(new fastjet::ClusterSequenceArea(ktchparticles, jet_def_bge, area_def_bge));

    rawJets   = fastjet::sorted_by_pt(cs->inclusive_jets(MinJetPt)); // APPLY Min pt cut for jet
    rawKtJets = fastjet::sorted_by_pt(cs_bge->inclusive_jets(0.0)); // APPLY Min pt cut for jet

    // Here one can choose to calculate background from kt jets which has left out the dijet with
    // delta phi cut used.
    if(fUseDeltaPhiBGSubtr) {
        removed = false;
        for (uktjet = 1; uktjet < rawKtJets.size(); uktjet++) { // First jet is already skipped here.
            if (!removed
             && rawKtJets[uktjet].eta() < etaMaxCutForKtJet
             && CheckDeltaPhi(rawKtJets[0], rawKtJets[uktjet], TMath::Pi()/fDeltaPhiCut)) {
                removed = true;
                continue;
            }
            rhoEstJets.push_back(rawKtJets[uktjet]); 
        }
    } else {
        rhoEstJets = selectorBoth(rawKtJets);
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
    bEvtHasAreaInfo = true;
}

void AliJCDijetAna::SetJets(vector<fastjet::PseudoJet> jetsOutside) {
    ResetObjects();
    rawJets = jetsOutside;
    bEvtHasAreaInfo = false;
}

void AliJCDijetAna::FillJetsDijets(AliJCDijetHistos *fhistos, int lCBin) {
    TString sDijetTypes[jetClassesSize] = {"raw", "bg. subtr.", "bg. subtr. const. cut", "const. cut", "kt"};
    int iAccJetCounter = 0;
    // anti-kt jets:
    for (ujet = 0; ujet < rawJets.size(); ujet++) {
        eta = rawJets[ujet].eta();
        fhistos->fh_events[lCBin]->Fill("jets",1.0);
        // anti-kt-jet eta cut
        if(TMath::Abs(eta) < etaMaxCutForJet) {
            iAccJetCounter++;
            fhistos->fh_events[lCBin]->Fill("acc. jets",1.0);
            jets[iAcc].push_back(rawJets[ujet]);
            pt = rawJets[ujet].pt();
            phi = rawJets[ujet].phi();
            if(bEvtHasAreaInfo) area = rawJets[ujet].area();
            if(bEvtHasAreaInfo) jetAreaVector = rawJets[ujet].area_4vector();
            fhistos->fh_jetEta[lCBin][iAcc]->Fill(eta);  
            fhistos->fh_jetPhi[lCBin][iAcc]->Fill(phi - TMath::Pi()); //Pseudojet.phi range 0-2pi
            fhistos->fh_jetEtaPhi[lCBin][iAcc]->Fill(eta,phi - TMath::Pi());
            fhistos->fh_jetPt[lCBin][iAcc]->Fill(pt);
            fhistos->fh_jetPt_ALICE[lCBin][iAcc]->Fill(pt);
            if(bEvtHasAreaInfo) fhistos->fh_jetArea[lCBin][iAcc]->Fill(area);
            if(bEvtHasAreaInfo) fhistos->fh_jetAreaRho[lCBin][iAcc]->Fill(area*rho);

            leadingTrackOverThreshold=false;
            if(fDebug > 9) cout << "Jet i=" << ujet << ", jet pt=" << pt << endl;
            for(uconst=0;uconst<rawJets[ujet].constituents().size(); uconst++) {
                if(fDebug > 9) cout << "Constituent i=" << uconst << ", constituent pt=" << rawJets[ujet].constituents()[uconst].pt() << endl;
                if(rawJets[ujet].constituents()[uconst].pt() > fConstituentCut) { // Jet leading constituent cut.
                    leadingTrackOverThreshold=true;
                    break;
                }
            }

            // Fill histos for jets with at least one constituent pt over fConstituentCut
            if(leadingTrackOverThreshold) {
                fhistos->fh_events[lCBin]->Fill("const. cut jets",1.0);
                fhistos->fh_jetEta[lCBin][iConstCut]->Fill(eta);
                fhistos->fh_jetPhi[lCBin][iConstCut]->Fill(phi - TMath::Pi());
                fhistos->fh_jetEtaPhi[lCBin][iConstCut]->Fill(eta,phi - TMath::Pi());
                fhistos->fh_jetPt[lCBin][iConstCut]->Fill(pt);
                fhistos->fh_jetPt_ALICE[lCBin][iConstCut]->Fill(pt);
                if(bEvtHasAreaInfo) fhistos->fh_jetArea[lCBin][iConstCut]->Fill(area);
                if(bEvtHasAreaInfo) fhistos->fh_jetAreaRho[lCBin][iConstCut]->Fill(area*rho);

                jets[iConstCut].push_back(rawJets[ujet]);
            }

            if(bEvtHasAreaInfo) {
                jet_bgSubtracted = fastjet::PseudoJet(rawJets[ujet].px() -        rho * jetAreaVector.px(),
                                                      rawJets[ujet].py() -        rho * jetAreaVector.py(),
                                                      rawJets[ujet].pz() - (rho+rhom) * jetAreaVector.pz(),
                                                      rawJets[ujet].E()  - (rho+rhom) * jetAreaVector.E());
                fhistos->fh_jetBGSubtrDeltaR[lCBin]->Fill(DeltaR(rawJets[ujet], jet_bgSubtracted));

                // Check eta acceptance also for bg subtracted jets.
                eta = jet_bgSubtracted.eta();
                if(TMath::Abs(eta) < etaMaxCutForJet) {
                    fhistos->fh_events[lCBin]->Fill("bg. subtr. jets",1.0);
                    pt2 = jet_bgSubtracted.pt();
                    phi = jet_bgSubtracted.phi();
                    if(iAccJetCounter==0 && pt>fLeadingJetCut && pt2<=fLeadingJetCut)       fhistos->fh_events[lCBin]->Fill("leading jet drop",1.0);
                    if(iAccJetCounter==1 && pt>fSubleadingJetCut && pt2<=fSubleadingJetCut) fhistos->fh_events[lCBin]->Fill("subleading jet drop",1.0);
                    fhistos->fh_jetEta[lCBin][iBGSubtr]->Fill(eta);
                    fhistos->fh_jetPhi[lCBin][iBGSubtr]->Fill(phi - TMath::Pi());
                    fhistos->fh_jetEtaPhi[lCBin][iBGSubtr]->Fill(eta,phi - TMath::Pi());
                    fhistos->fh_jetPt[lCBin][iBGSubtr]->Fill(pt2);
                    fhistos->fh_jetPt_ALICE[lCBin][iBGSubtr]->Fill(pt2);
                    if(bEvtHasAreaInfo) fhistos->fh_jetArea[lCBin][iBGSubtr]->Fill(area); // Assuming bg subtracted jet has the same area.
                    if(bEvtHasAreaInfo) fhistos->fh_jetAreaRho[lCBin][iBGSubtr]->Fill(area*rho);

                    jets[iBGSubtr].push_back(jet_bgSubtracted);

                    if(leadingTrackOverThreshold) {
                        fhistos->fh_events[lCBin]->Fill("bg. subtr. const. cut jets",1.0);
                        fhistos->fh_jetEta[lCBin][iBGSubtrConstCut]->Fill(eta);
                        fhistos->fh_jetPhi[lCBin][iBGSubtrConstCut]->Fill(phi - TMath::Pi());
                        fhistos->fh_jetEtaPhi[lCBin][iBGSubtrConstCut]->Fill(eta,phi - TMath::Pi());
                        fhistos->fh_jetPt[lCBin][iBGSubtrConstCut]->Fill(pt2);
                        fhistos->fh_jetPt_ALICE[lCBin][iBGSubtrConstCut]->Fill(pt2);
                        if(bEvtHasAreaInfo) fhistos->fh_jetArea[lCBin][iBGSubtrConstCut]->Fill(area);
                        if(bEvtHasAreaInfo) fhistos->fh_jetAreaRho[lCBin][iBGSubtrConstCut]->Fill(area*rho);

                        jets[iBGSubtrConstCut].push_back(jet_bgSubtracted);
                    }
                }
            }
        }
    }//end of the anti-kt-jet loop

    for (uktjet = 0; uktjet < rawKtJets.size(); uktjet++) {
        eta = rawKtJets[uktjet].eta();
        fhistos->fh_events[lCBin]->Fill("kt-jets",1.0);
        // kt-jet eta cut
        if(TMath::Abs(eta) < etaMaxCutForKtJet) {
            fhistos->fh_events[lCBin]->Fill("acc. kt-jets",1.0);
            jets[iktJets].push_back(rawKtJets[uktjet]);
            pt = rawKtJets[uktjet].pt();
            phi = rawKtJets[uktjet].phi();
            if(bEvtHasAreaInfo) area = rawKtJets[uktjet].area();
            fhistos->fh_jetEta[lCBin][iktJets]->Fill(eta);  
            fhistos->fh_jetPhi[lCBin][iktJets]->Fill(phi - TMath::Pi()); //Pseudojet.phi range 0-2pi
            fhistos->fh_jetEtaPhi[lCBin][iktJets]->Fill(eta,phi - TMath::Pi());
            fhistos->fh_jetPt[lCBin][iktJets]->Fill(pt);
            fhistos->fh_jetPt_ALICE[lCBin][iktJets]->Fill(pt);
            if(bEvtHasAreaInfo) fhistos->fh_jetArea[lCBin][iktJets]->Fill(area);
            if(bEvtHasAreaInfo) fhistos->fh_jetAreaRho[lCBin][iktJets]->Fill(area*rho);
        }
    } //end of the kt-jet loop

    // To measure how well the matching performs:
    // Search the smallest deltaR between jets in a single event.
    for(udijet=0; udijet < jetClassesSize; udijet++) {
        deltaRMin = 999.0;
        if(jets[udijet].size()<2) continue;
        for(ujet = 0; ujet < jets[udijet].size(); ujet++) {
            for(ujet2 = ujet+1; ujet2 < jets[udijet].size(); ujet2++) {
                deltaR = DeltaR(jets[udijet][ujet], jets[udijet][ujet2]);
                    if(deltaR < deltaRMin) {
                        deltaRMin=deltaR;
                    }
            }
        }
        fhistos->fh_jetDeltaRMin[lCBin][udijet]->Fill(deltaRMin);
    }


    // Dijet calculations 
    for(udijet=0; udijet < jetClassesSize; udijet++) {
        if(jets[udijet].size()>1) {
            jets[udijet] = fastjet::sorted_by_pt(jets[udijet]); // Sort in case of bg subtr messed up the order.
            fhistos->fh_events[lCBin]->Fill(Form("%s dijets",sDijetTypes[udijet].Data()),1.0);

            // Start by forming dijets. First without deltaPhiCut and then with the cut.
            
            // No deltaPhi cut for these jets.
            dijets[udijet][0][0] = jets[udijet][0];
            if(dijets[udijet][0][0].pt() < fLeadingJetCut) continue;
            fhistos->fh_events[lCBin]->Fill(Form("%s dijets leading cut",sDijetTypes[udijet].Data()),1.0);
            dijets[udijet][0][1] = jets[udijet][1];

            // Here we check deltaPhi cut
            dijets[udijet][1][0] = jets[udijet][0];
            bHasDeltaPhiSubLeadJet = false;
            for (ujet = 1; ujet < jets[udijet].size(); ujet++) {
                if(CheckDeltaPhi(dijets[udijet][1][0], jets[udijet][ujet], TMath::Pi()/fDeltaPhiCut)) {
                        dijets[udijet][1][1] = jets[udijet][ujet];
                        bHasDeltaPhiSubLeadJet = true;
                        break; // list is pt-ordered. The first jet to have deltaPhi check ok, is the pair.
                }
            }

            // Analysis for dijet without deltaPhi cut.
            if(dijets[udijet][0][1].pt()>fSubleadingJetCut) {
                fhistos->fh_events[lCBin]->Fill(Form("%s acc. dijets",sDijetTypes[udijet].Data()),1.0);
                if(udijet==iAcc) bHasDijet = true;
                dijet = dijets[udijet][0][0] + dijets[udijet][0][1];
                mjj = dijet.m();
                ptpair = dijet.pt();
                fhistos->fh_dijetInvM[lCBin][udijet]->Fill(mjj);
                fhistos->fh_dijetPtPair[lCBin][udijet]->Fill(ptpair);
                dPhi = GetDeltaPhi(dijets[udijet][0][0], dijets[udijet][0][1]);
                fhistos->fh_dijetDeltaPhi[lCBin][udijet]->Fill(dPhi);
            }

            // Analysis for dijet with deltaPhi cut.
            if(bHasDeltaPhiSubLeadJet && dijets[udijet][1][1].pt()>fSubleadingJetCut) {
                fhistos->fh_events[lCBin]->Fill(Form("%s deltaphi cut dijets",sDijetTypes[udijet].Data()),1.0);
                if(udijet==iAcc) bHasDeltaPhiDijet = true;
                dijet = dijets[udijet][1][0] + dijets[udijet][1][1];
                mjj = dijet.m();
                ptpair = dijet.pt();
                fhistos->fh_dijetInvMDeltaPhiCut[lCBin][udijet]->Fill(mjj);
                fhistos->fh_dijetPtPairDeltaPhiCut[lCBin][udijet]->Fill(ptpair);
                dPhi = GetDeltaPhi(dijets[udijet][1][0], dijets[udijet][1][1]);
                fhistos->fh_dijetDeltaPhiWithCut[lCBin][udijet]->Fill(dPhi);
            }
        }
    }
}

// Response matrices are calculated in this function.
void AliJCDijetAna::CalculateResponse(AliJCDijetAna *anaDetMC, AliJCDijetHistos *fhistos) {
    
    vector<vector<fastjet::PseudoJet>> jetsDetMC = anaDetMC->GetJets();

    unsigned Njets = jets[iAcc].size();
    unsigned NjetsDetMC = jetsDetMC[iAcc].size();
    double maxpt=0;
    double minR=0;
    double deltaRMatch=0;
    double ptTrue, ptDetMC;
    unsigned maxptIndex;
    bool bfound;
    std::vector<bool> bJetMatch(NjetsDetMC, false);
    bool bLeadingMatch    = false;
    bool bSubleadingMatch = false;
    bool bSubleadingMatchDeltaPhi = false;
    // for raw anti-kt jets:
    for (ujet = 0; ujet < Njets; ujet++) { //True MC jets
        maxpt=0;
        maxptIndex=-1;
        bfound=false;
        deltaR=0;
        deltaRMatch=0;
        minR=999.0;
        for (ujetDetMC = 0; ujetDetMC < NjetsDetMC; ujetDetMC++) { //Det MC jets
            deltaR = DeltaR(jets[iAcc][ujet], jetsDetMC[iAcc][ujetDetMC]);
            if(deltaR<minR) minR=deltaR;
            if(deltaR < matchingR && jetsDetMC[iAcc][ujetDetMC].pt() > maxpt) {
                maxpt = jetsDetMC[iAcc][ujetDetMC].pt();
                maxptIndex = ujetDetMC;
                deltaRMatch = deltaR;
                bfound = true;
                if(ujet==0 && ujetDetMC==0) bLeadingMatch   = true;
                if(ujet==1 && ujetDetMC==1) bSubleadingMatch= true;
                //cout << "found, detPt vs truePt: " << maxpt << " <> " << jets[iAcc][ujet].pt() << ", index: " << maxptIndex << endl;
            }
        }
        fhistos->fh_jetResponseDeltaRClosest->Fill(minR);
        if(bfound) {
            ptTrue = jets[iAcc][ujet].pt();
            ptDetMC = jetsDetMC[iAcc][maxptIndex].pt();
            fhistos->fh_jetResponse->Fill(ptDetMC, ptTrue);
            fhistos->fh_jetResponse_ALICE->Fill(ptDetMC, ptTrue);
            fhistos->fh_jetResponseDeltaR->Fill(deltaRMatch);
            fhistos->fh_jetResponseDeltaPt->Fill((ptTrue-ptDetMC)/ptTrue);
            fhistos->fh_responseInfo->Fill("True jet has pair",1.0);
            bJetMatch[maxptIndex] = true;
        } else {
            fhistos->fh_responseInfo->Fill("True jet has no pair",1.0);
        }
    }
    for (ujetDetMC = 0; ujetDetMC < NjetsDetMC; ujetDetMC++) { //Det MC jets
        if(!bJetMatch[ujetDetMC]) {
            fhistos->fh_responseInfo->Fill("Det jet has no pair",1.0);
        }
    }

    vector<vector<vector<fastjet::PseudoJet>>> dijetsDetMC = anaDetMC->GetDijets();
    fastjet::PseudoJet dijetDetMC;

    //Dijet response without deltaphi cut.
    if(bHasDijet) {
        dijet = dijets[iAcc][0][0] + dijets[iAcc][0][1];
        if(anaDetMC->HasDijet()) {
            dijetDetMC = dijetsDetMC[iAcc][0][0] + dijetsDetMC[iAcc][0][1];
            if(bLeadingMatch && bSubleadingMatch) {
                fhistos->fh_dijetResponse->Fill(dijetDetMC.m(), dijet.m());
                fhistos->fh_responseInfo->Fill("Dijet match",1.0);
            } else {
                fhistos->fh_responseInfo->Fill("Dijet not match",1.0);
            }
        } else {
            fhistos->fh_responseInfo->Fill("Dijet det not found",1.0);
        }
    } else {
        if(anaDetMC->HasDijet()) {
            //dijetDetMC = dijetsDetMC[iAcc][0][0] + dijetsDetMC[iAcc][0][1];
            fhistos->fh_responseInfo->Fill("Dijet true not found",1.0);
        }
    }

    // DeltaPhi cut dijet response.
    if(bHasDeltaPhiDijet) {
        dijet = dijets[iAcc][1][0] + dijets[iAcc][1][1];
        if(anaDetMC->HasDeltaPhiDijet()) {
            dijetDetMC = dijetsDetMC[iAcc][1][0] + dijetsDetMC[iAcc][1][1];
            // Check subleading jet match.
            if(DeltaR(dijets[iAcc][1][1], dijetsDetMC[iAcc][1][1]) < matchingR) {
                bSubleadingMatchDeltaPhi = true;
            }

            if(bLeadingMatch && bSubleadingMatchDeltaPhi) {
                fhistos->fh_dijetResponseDeltaPhiCut->Fill(dijetDetMC.m(), dijet.m());
                fhistos->fh_responseInfo->Fill("Dijet DPhi match",1.0);
            } else {
                fhistos->fh_responseInfo->Fill("Dijet DPhi not match",1.0);
            }
        } else {
            fhistos->fh_responseInfo->Fill("Dijet DPhi det not found",1.0);
        }
    } else {
        if(anaDetMC->HasDeltaPhiDijet()) {
            //dijetDetMC = dijetsDetMC[iAcc][1][0] + dijetsDetMC[iAcc][1][1];
            fhistos->fh_responseInfo->Fill("Dijet DPhi true not found",1.0);
        }
    }
    return;
}

void AliJCDijetAna::ResetObjects() {
    chparticles.clear();
    ktchparticles.clear();
    for(int i=0; i<jetClassesSize; i++) jets[i].clear();
    rawJets.clear();
    rawKtJets.clear();
    rhoEstJets.clear();
    constituents.clear();

    phi=0, eta=0, pt=0, pt2=0, rho=0, rhom=0, area=0, mjj=0, ptpair=0, dPhi=0;
    leadingTrackOverThreshold = false;
    bHasDijet = false;
    bHasDeltaPhiDijet = false;
    bHasDeltaPhiSubLeadJet = false;
    bEvtHasAreaInfo = false;
}

double AliJCDijetAna::DeltaR(fastjet::PseudoJet jet1, fastjet::PseudoJet jet2) {
    double Deta = jet1.eta() - jet2.eta();
    double Dphi = TMath::Abs(jet1.phi() - jet2.phi());
    // Make sure that Dphi is in 0-pi range.
    Dphi = Dphi>TMath::Pi() ? 2*TMath::Pi()-Dphi : Dphi;

    return TMath::Sqrt(Deta*Deta + Dphi*Dphi);
}

bool AliJCDijetAna::CheckDeltaPhi(fastjet::PseudoJet leadingJet, fastjet::PseudoJet subleadingJet, double deltaPhiCut) {
    double DeltaPhi;
    DeltaPhi = GetDeltaPhi(leadingJet, subleadingJet);
    if(TMath::Abs(DeltaPhi - TMath::Pi()) < deltaPhiCut) {
        return true;
    }
    return false;
}

double AliJCDijetAna::GetDeltaPhi(fastjet::PseudoJet jet1, fastjet::PseudoJet jet2) {
    double DeltaPhi, DeltaPhi2;
    DeltaPhi  = jet2.delta_phi_to(jet1);
    DeltaPhi2 = DeltaPhi<0 ? DeltaPhi+TMath::TwoPi() : DeltaPhi;
    return DeltaPhi2;
}

// This should be done after SetSettings.
void AliJCDijetAna::InitHistos(AliJCDijetHistos *histos, bool bIsMC, int nCentBins) {
    histos->fh_info->Fill("Count", 1.0);
    histos->fh_info->Fill("MC", bIsMC);
    for(int i=0; i< nCentBins; i++) histos->fh_info->Fill(Form("Cent bin border %02d",i), nCentBins);
    histos->fh_info->Fill("Jet cone", fJetCone);
    histos->fh_info->Fill("kt-jet cone", fktJetCone);
    histos->fh_info->Fill("kt scheme", ktScheme);
    histos->fh_info->Fill("antikt scheme", antiktScheme);
    histos->fh_info->Fill("Use pion mass", fusePionMass);
    histos->fh_info->Fill("Use DeltaPhi BG Subtr", fUseDeltaPhiBGSubtr);
    histos->fh_info->Fill("Particle eta cut", fParticleEtaCut);
    histos->fh_info->Fill("Particle pt cut", fParticlePtCut);
    histos->fh_info->Fill("Leading jet cut", fLeadingJetCut);
    histos->fh_info->Fill("Subleading jet cut", fSubleadingJetCut);
    histos->fh_info->Fill("Const. cut", fConstituentCut);
    histos->fh_info->Fill("Delta phi cut pi/",fDeltaPhiCut);
    histos->fh_info->Fill("Matching R for MC",matchingR);

    // Initialize histos->fh_events so that the bin order is correct
    for (int iBin=0; iBin < nCentBins-1; iBin++) {
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
        if(bIsMC) {
            histos->fh_responseInfo->Fill("True jet has pair",0.0);
            histos->fh_responseInfo->Fill("True jet has no pair",0.0);
            histos->fh_responseInfo->Fill("Det jet has no pair",0.0);
            histos->fh_responseInfo->Fill("Dijet match",0.0);
            histos->fh_responseInfo->Fill("Dijet not match",0.0);
            histos->fh_responseInfo->Fill("Dijet det not found",0.0);
            histos->fh_responseInfo->Fill("Dijet true not found",0.0);
            histos->fh_responseInfo->Fill("Dijet DPhi match",0.0);
            histos->fh_responseInfo->Fill("Dijet DPhi not match",0.0);
            histos->fh_responseInfo->Fill("Dijet DPhi det not found",0.0);
            histos->fh_responseInfo->Fill("Dijet DPhi true not found",0.0);
        }
    }
}
#endif
