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
                                double lMinJetPt,
                                double lDeltaPhiCut,
                                double lmatchingR,
                                double ltrackingIneff){
    fDebug = lDebug;
    fParticleEtaCut = lParticleEtaCut;
    fParticlePtCut = lParticlePtCut;
    fusePionMass = lusePionMass;
    fUseDeltaPhiBGSubtr = luseDeltaPhiBGSubtr;
    fConstituentCut = lConstituentCut;
    fLeadingJetCut = lLeadingJetCut;
    fSubleadingJetCut = lSubleadingJetCut;
    fDeltaPhiCut = lDeltaPhiCut;
    ftrackingIneff = ltrackingIneff;

    etaMaxCutForJet = lParticleEtaCut-lJetCone;
    etaMaxCutForKtJet = lParticleEtaCut-lktJetCone;
    MinJetPt = lMinJetPt; // Min Jet Pt cut to disregard low pt jets
    double const ghost_maxrap = lParticleEtaCut;
    unsigned int const repeat = 1; // default
    double const ghost_area   = 0.005; // ALICE=0.005 // default=0.01
    fJetCone = lJetCone;
    fktJetCone = lktJetCone;
    pionmass = 0.139570;//AliPID::ParticleMass(AliPID::kPion);
    matchingR = lmatchingR;
    fDeltaPt = -9999.0;
    fptHardBin = 0; // Default value

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
    randomGenerator = new TRandom3();

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
    jet_def = fastjet::JetDefinition(fastjet::antikt_algorithm, fJetCone, antiktScheme);
    jet_def_bge = fastjet::JetDefinition(fastjet::kt_algorithm, fktJetCone, ktScheme);

    area_spec = fastjet::GhostedAreaSpec(ghost_maxrap, repeat, ghost_area);
    area_def = fastjet::AreaDefinition(fastjet::active_area, area_spec);
    area_def_bge = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, area_spec);

    // Selector selects first all jets inside rapidity acceptance and then all but two hardest jets.
    selectorAllButTwo  = fastjet::Selector(!fastjet::SelectorNHardest(2));
    selectorEta  = fastjet::SelectorAbsEtaMax(ghost_maxrap - fktJetCone);
    selectorBoth  = selectorAllButTwo * selectorEta; // Here right selector is applied first, then the left one.
    bge = fastjet::JetMedianBackgroundEstimator(selectorEta, jet_def_bge, area_def_bge);


    return;
}

//______________________________________________________________________________
int AliJCDijetAna::CalculateJets(TClonesArray *inList, AliJCDijetHistos *fhistos, int lCBin){

    ResetObjects();

    //--------------------------------------------------------
    //         B e g i n    e v e n t    l o o p.
    //--------------------------------------------------------
    if (inList==0) { cout << "No list!" << endl; return 0;}
    noTracks = inList->GetEntriesFast();
    //cout << "Number of tracks: " << noTracks << endl;
    
    //Initialize random cone variables for delta-pt study
    fDeltaPt = -9999.0;
    randConePhi = randomGenerator->Uniform(-TMath::Pi(),TMath::Pi());
    randConeEta = randomGenerator->Uniform(-etaMaxCutForJet,etaMaxCutForJet);
    randConePt = 0.0;


    //cout << "Number of tracks in Ana code: " << noTracks << endl;
    for (utrack = 0; utrack < noTracks; ++utrack) {//loop over all the particles in the event
        // Building input particle list for the jet reconstruction
        AliJBaseTrack *trk = (AliJBaseTrack*)inList->At(utrack);
        pt = trk->Pt();
        eta = trk->Eta();
        if (pt>fParticlePtCut && TMath::Abs(eta) < fParticleEtaCut){
            if(ftrackingIneff>0.0 && randomGenerator->Uniform(0.0,1.0) < ftrackingIneff) continue;
            phi = trk->Phi();
            if(DeltaR(randConeEta, eta, randConePhi, phi) < fJetCone) randConePt += pt;
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

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Run the clustering, Reconstruct jets
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    cs.reset(new fastjet::ClusterSequenceArea(chparticles, jet_def, area_def));
    cs_bge.reset(new fastjet::ClusterSequenceArea(ktchparticles, jet_def_bge, area_def_bge));

    rawJets   = fastjet::sorted_by_pt(cs->inclusive_jets(MinJetPt)); // APPLY Min pt cut for jet

    // For MC runs: If we find jets with over 4 times pt_hard bin, reject the event.
    if( fptHardBin!=0 && rawJets.size()>0 ) {
        fhistos->fh_ptHard[lCBin]->Fill(fptHardBin);
        fhistos->fh_maxJetptOverPtHard[lCBin]->Fill(rawJets.at(0).pt()/fptHardBin);
        if( rawJets.at(0).pt() > fptHardBin*4 ) {
            fhistos->fh_events[lCBin]->Fill("pt_hard bin cuts",1.0);
            return -1;
        }
    }
    fhistos->fh_randConeEtaPhi[lCBin]->Fill(randConeEta,randConePhi);

    rawKtJets = fastjet::sorted_by_pt(cs_bge->inclusive_jets(0.0)); // APPLY Min pt cut for jet

    fhistos->fh_events[lCBin]->Fill("particles",noTracks);
    for (utrack = 0; utrack < chparticles.size(); utrack++) {
        fhistos->fh_events[lCBin]->Fill("acc. particles",1.0);
        pt = chparticles.at(utrack).pt();
        eta = chparticles.at(utrack).eta();
        phi = chparticles.at(utrack).phi();
        fhistos->fh_eta[lCBin]->Fill(eta);
        fhistos->fh_phi[lCBin]->Fill(phi);
        fhistos->fh_etaPhi[lCBin]->Fill(eta,phi);
        fhistos->fh_pt[lCBin]->Fill(pt);
        if(eta>0.0) fhistos->fh_ptPosEta[lCBin]->Fill(pt);
        else        fhistos->fh_ptNegEta[lCBin]->Fill(pt);
    }
    fhistos->fh_nch[lCBin]->Fill(chparticles.size());
    if(chparticles.size()==0) return 0; // We are not intereted in empty events.

    // Here one can choose to calculate background from kt jets which has left out the dijet with
    // delta phi cut used.
    if(fUseDeltaPhiBGSubtr) {
        removed = false;
        for (uktjet = 1; uktjet < rawKtJets.size(); uktjet++) { // First jet is already skipped here.
            if (!removed
             && rawKtJets.at(uktjet).eta() < etaMaxCutForKtJet
             && CheckDeltaPhi(rawKtJets.at(0), rawKtJets.at(uktjet), TMath::Pi()/fDeltaPhiCut)) {
                removed = true;
                continue;
            }
            rhoEstJets.push_back(rawKtJets.at(uktjet)); 
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


    fhistos->fh_rho[lCBin]->Fill(rho);
    fhistos->fh_rhoLin[lCBin]->Fill(rho);
    fhistos->fh_rhom[lCBin]->Fill(rhom);
    fhistos->fh_rhomLin[lCBin]->Fill(rhom);
    if(fDebug > 9) std::cout << "Testing: Rho_M = " << rhom << ", has_rho_m() = " << bge.has_rho_m() << std::endl;

    // This is the delta-phi distribution: delta-pt = sum_{tracks in rand cone}( pt ) - rho*pi*R^2
    fDeltaPt = randConePt - rho*TMath::Pi()*fJetCone*fJetCone;
    fhistos->fh_deltaPt[lCBin]->Fill(fDeltaPt);

    bEvtHasAreaInfo = true;
    return 0;
}

// If user wants to add jets from outside this function can be used.
void AliJCDijetAna::SetJets(vector<fastjet::PseudoJet> jetsOutside) {
    ResetObjects();
    rawJets = fastjet::sorted_by_pt(jetsOutside);
    bEvtHasAreaInfo = false;
    return;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Loop over jets and fill various histos 
// Note: Jets in rawJets are assumed to be sorted by pt, as is done
//       in the previous functions.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliJCDijetAna::FillJetsDijets(AliJCDijetHistos *fhistos, int lCBin) {
    TString sDijetTypes[jetClassesSize] = {"raw", "bg. subtr.", "bg. subtr. const. cut", "const. cut", "kt"};
    int iAccJetCounter = 0;
    bool bHasHighPtJet = false;

    // anti-kt jets:
    for (ujet = 0; ujet < rawJets.size(); ujet++) {
        eta = rawJets.at(ujet).eta();
        fhistos->fh_events[lCBin]->Fill("jets",1.0);
        // anti-kt-jet eta cut
        if(TMath::Abs(eta) < etaMaxCutForJet) {
            iAccJetCounter++;
            fhistos->fh_events[lCBin]->Fill("acc. jets",1.0);
            jets.at(iAcc).push_back(rawJets.at(ujet));
            pt = rawJets.at(ujet).pt();
            phi = rawJets.at(ujet).phi();
            if(iAccJetCounter==1) bHasHighPtJet = pt > fLeadingJetCut;
            if(bEvtHasAreaInfo) area = rawJets.at(ujet).area();
            if(bEvtHasAreaInfo) jetAreaVector = rawJets.at(ujet).area_4vector();
            fhistos->fh_jetEta[lCBin][iAcc]->Fill(eta);  
            fhistos->fh_jetPhi[lCBin][iAcc]->Fill(phi - TMath::Pi()); //Pseudojet.phi range 0-2pi
            fhistos->fh_jetEtaPhi[lCBin][iAcc]->Fill(eta,phi - TMath::Pi());
            fhistos->fh_jetPt[lCBin][iAcc]->Fill(pt);
            fhistos->fh_jetPt_ALICE[lCBin][iAcc]->Fill(pt);
            if(bEvtHasAreaInfo) fhistos->fh_jetArea[lCBin][iAcc]->Fill(area);
            if(bEvtHasAreaInfo) fhistos->fh_jetAreaRho[lCBin][iAcc]->Fill(area*rho);

            leadingTrackOverThreshold=false;
            if(fDebug > 9) cout << "Jet i=" << ujet << ", jet pt=" << pt << endl;
            for(uconst=0;uconst<rawJets.at(ujet).constituents().size(); uconst++) {
                if(fDebug > 9) cout << "Constituent i=" << uconst << ", constituent pt=" << rawJets.at(ujet).constituents().at(uconst).pt() << endl;
                if(rawJets.at(ujet).constituents().at(uconst).pt() > fConstituentCut) { // Jet leading constituent cut.
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

                jets.at(iConstCut).push_back(rawJets.at(ujet));
            }

            if(bEvtHasAreaInfo) {
                jet_bgSubtracted = fastjet::PseudoJet(rawJets.at(ujet).px() -        rho * jetAreaVector.px(),
                                                      rawJets.at(ujet).py() -        rho * jetAreaVector.py(),
                                                      rawJets.at(ujet).pz() - (rho+rhom) * jetAreaVector.pz(),
                                                      rawJets.at(ujet).E()  - (rho+rhom) * jetAreaVector.E());
                fhistos->fh_jetBGSubtrDeltaR[lCBin]->Fill(DeltaR(rawJets.at(ujet), jet_bgSubtracted));

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
                    fhistos->fh_jetPtTransBGSub[lCBin][iBGSubtr]->Fill(pt - rho*area);
                    if(bEvtHasAreaInfo) fhistos->fh_jetArea[lCBin][iBGSubtr]->Fill(area); // Assuming bg subtracted jet has the same area.
                    if(bEvtHasAreaInfo) fhistos->fh_jetAreaRho[lCBin][iBGSubtr]->Fill(area*rho);

                    jets.at(iBGSubtr).push_back(jet_bgSubtracted);

                    if(leadingTrackOverThreshold) {
                        fhistos->fh_events[lCBin]->Fill("bg. subtr. const. cut jets",1.0);
                        fhistos->fh_jetEta[lCBin][iBGSubtrConstCut]->Fill(eta);
                        fhistos->fh_jetPhi[lCBin][iBGSubtrConstCut]->Fill(phi - TMath::Pi());
                        fhistos->fh_jetEtaPhi[lCBin][iBGSubtrConstCut]->Fill(eta,phi - TMath::Pi());
                        fhistos->fh_jetPt[lCBin][iBGSubtrConstCut]->Fill(pt2);
                        fhistos->fh_jetPt_ALICE[lCBin][iBGSubtrConstCut]->Fill(pt2);
                        fhistos->fh_jetPtTransBGSub[lCBin][iBGSubtrConstCut]->Fill(pt - rho*area);
                        if(bEvtHasAreaInfo) fhistos->fh_jetArea[lCBin][iBGSubtrConstCut]->Fill(area);
                        if(bEvtHasAreaInfo) fhistos->fh_jetAreaRho[lCBin][iBGSubtrConstCut]->Fill(area*rho);

                        jets.at(iBGSubtrConstCut).push_back(jet_bgSubtracted);
                    }
                }
            }
        }
    }//end of the anti-kt-jet loop
    //============================================================================================
    //cout << "Jets: " << jets[iAcc].size() << endl;
    //============================================================================================
    fhistos->fh_jetN[lCBin]->Fill(iAccJetCounter);
    if(bHasHighPtJet) {
        fhistos->fh_rhoHighPt[lCBin]->Fill(rho);
        fhistos->fh_rhoLinHighPt[lCBin]->Fill(rho);
        fhistos->fh_rhomHighPt[lCBin]->Fill(rhom);
        fhistos->fh_rhomLinHighPt[lCBin]->Fill(rhom);
    }

    for (uktjet = 0; uktjet < rawKtJets.size(); uktjet++) {
        eta = rawKtJets.at(uktjet).eta();
        fhistos->fh_events[lCBin]->Fill("kt-jets",1.0);
        // kt-jet eta cut
        if(TMath::Abs(eta) < etaMaxCutForKtJet) {
            fhistos->fh_events[lCBin]->Fill("acc. kt-jets",1.0);
            jets.at(iktJets).push_back(rawKtJets.at(uktjet));
            pt = rawKtJets.at(uktjet).pt();
            phi = rawKtJets.at(uktjet).phi();
            if(bEvtHasAreaInfo) area = rawKtJets.at(uktjet).area();
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
        if(jets.at(udijet).size()<2) continue;
        for(ujet = 0; ujet < jets.at(udijet).size(); ujet++) {
            for(ujet2 = ujet+1; ujet2 < jets.at(udijet).size(); ujet2++) {
                deltaR = DeltaR(jets.at(udijet).at(ujet), jets.at(udijet).at(ujet2));
                    if(deltaR < deltaRMin) {
                        deltaRMin=deltaR;
                    }
            }
        }
        fhistos->fh_jetDeltaRMin[lCBin][udijet]->Fill(deltaRMin);
    }

    /*
    for (ujet = 0; ujet < jets.at(iAcc).size(); ujet++) {
        cout << "jet(E, pt, phi, eta) = " <<
            "(" << jets.at(iAcc).at(ujet).E()   << ", "
            << jets.at(iAcc).at(ujet).pt()  << ", "
            << jets.at(iAcc).at(ujet).phi() << ", "
            << jets.at(iAcc).at(ujet).eta() << ")" << endl;
    }
    */


    // Dijet calculations 
    for(udijet=0; udijet < jetClassesSize; udijet++) {
        if(jets.at(udijet).size()>1) {
            jets.at(udijet) = fastjet::sorted_by_pt(jets.at(udijet)); // Sort in case of bg subtr messed up the order.
            fhistos->fh_events[lCBin]->Fill(Form("%s dijets",sDijetTypes[udijet].Data()),1.0);

            // Start by forming dijets. First without deltaPhiCut and then with the cut.
            
            // No deltaPhi cut for these jets.
            dijets.at(udijet).at(0).at(0) = jets.at(udijet).at(0);
            if(dijets.at(udijet).at(0).at(0).pt() < fLeadingJetCut) continue;
            fhistos->fh_events[lCBin]->Fill(Form("%s dijets leading cut",sDijetTypes[udijet].Data()),1.0);
            dijets.at(udijet).at(0).at(1) = jets.at(udijet).at(1);

            // Here we check deltaPhi cut
            dijets.at(udijet).at(1).at(0) = jets.at(udijet).at(0);
            bHasDeltaPhiSubLeadJet = false;
            for (ujet = 1; ujet < jets.at(udijet).size(); ujet++) {
                if(CheckDeltaPhi(dijets.at(udijet).at(1).at(0), jets.at(udijet).at(ujet), TMath::Pi()/fDeltaPhiCut)) {
                        dijets.at(udijet).at(1).at(1) = jets.at(udijet).at(ujet);
                        bHasDeltaPhiSubLeadJet = true;
                        break; // list is pt-ordered. The first jet to have deltaPhi check ok, is the pair.
                }
            }

            // Analysis for dijet without deltaPhi cut.
            if(dijets.at(udijet).at(0).at(1).pt()>fSubleadingJetCut) {
                fhistos->fh_events[lCBin]->Fill(Form("%s acc. dijets",sDijetTypes[udijet].Data()),1.0);
                if(udijet==iAcc) bHasDijet = true;
                dijet = dijets.at(udijet).at(0).at(0) + dijets.at(udijet).at(0).at(1);
                mjj = dijet.m();
                ptpair = dijet.pt();
                fhistos->fh_dijetInvM[lCBin][udijet]->Fill(mjj);
                fhistos->fh_dijetInvMTrunc[lCBin][udijet]->Fill(mjj);
                fhistos->fh_dijetPtPair[lCBin][udijet]->Fill(ptpair);
                dPhi = GetDeltaPhi(dijets.at(udijet).at(0).at(0), dijets.at(udijet).at(0).at(1));
                fhistos->fh_dijetDeltaPhi[lCBin][udijet]->Fill(dPhi);
            }

            // Analysis for dijet with deltaPhi cut.
            if(bHasDeltaPhiSubLeadJet && dijets.at(udijet).at(1).at(1).pt()>fSubleadingJetCut) {
                fhistos->fh_events[lCBin]->Fill(Form("%s deltaphi cut dijets",sDijetTypes[udijet].Data()),1.0);
                if(udijet==iAcc) bHasDeltaPhiDijet = true;
                dijet = dijets.at(udijet).at(1).at(0) + dijets.at(udijet).at(1).at(1);
                mjj = dijet.m();
                ptpair = dijet.pt();
                fhistos->fh_dijetInvMDeltaPhiCut[lCBin][udijet]->Fill(mjj);
                fhistos->fh_dijetInvMDeltaPhiCutTrunc[lCBin][udijet]->Fill(mjj);
                fhistos->fh_dijetPtPairDeltaPhiCut[lCBin][udijet]->Fill(ptpair);
                dPhi = GetDeltaPhi(dijets.at(udijet).at(1).at(0), dijets.at(udijet).at(1).at(1));
                fhistos->fh_dijetDeltaPhiWithCut[lCBin][udijet]->Fill(dPhi);
            }
        }
    }
    return;
}

// Response matrices are calculated in this function.
void AliJCDijetAna::CalculateResponse(AliJCDijetAna *anaDetMC, AliJCDijetHistos *fhistos) {
    if(fDebug>8) cout << "===== BEGIN RESPONSE CALC =====" << endl;
    
    vector<vector<fastjet::PseudoJet>> jetsDetMC = anaDetMC->GetJets();

    unsigned Njets = jets.at(iAcc).size();
    unsigned NjetsDetMC = jetsDetMC.at(iAcc).size();
    //cout << "Response true jets: " << Njets << endl;
    //cout << "Response det jets size:     " << NjetsDetMC << endl;
    double maxpt=0;
    double minR=0;
    double deltaRMatch=0;
    double ptTrue, ptDetMC;
    double deltaRLL, deltaRLS, deltaRSL, deltaRSS;
    unsigned maxptIndex;
    bool bfound;
    std::vector<bool> bTrueJetMatch(Njets, false);
    std::vector<bool> bDetJetMatch(NjetsDetMC, false);
    bool bLeadingMatch    = false;
    bool bSubleadingMatch = false;
    bool bSubleadingMatchDeltaPhi = false;

    double dBinCenter = 0;
    for(int iby = 1 ; iby <= fhistos->fh_deltaPtResponseEvery->GetNbinsY(); iby++){
        dBinCenter = fhistos->fh_deltaPtResponseEvery->GetYaxis()->GetBinCenter(iby);
        fhistos->fh_deltaPtResponseEvery->Fill(dBinCenter+fDeltaPt,dBinCenter);
    }
    for(int iby = 1 ; iby <= fhistos->fh_deltaPtResponseEvery_ALICE->GetNbinsX(); iby++){
        dBinCenter = fhistos->fh_deltaPtResponseEvery_ALICE->GetYaxis()->GetBinCenter(iby);
        fhistos->fh_deltaPtResponseEvery_ALICE->Fill(dBinCenter+fDeltaPt,dBinCenter);
    }
    for (ujet = 0; ujet < Njets; ujet++) { //True MC jets
        maxpt=0;
        maxptIndex=-1;
        bfound=false;
        deltaR=0;
        deltaRMatch=0;
        minR=999.0;
        ptTrue = jets.at(iAcc).at(ujet).pt();
        fhistos->fh_deltaPtResponse->Fill(ptTrue+fDeltaPt,ptTrue);
        fhistos->fh_deltaPtResponse_ALICE->Fill(ptTrue+fDeltaPt,ptTrue);
        for (ujetDetMC = 0; ujetDetMC < NjetsDetMC; ujetDetMC++) { //Det MC jets
            deltaR = DeltaR(jets.at(iAcc).at(ujet), jetsDetMC.at(iAcc).at(ujetDetMC));
            if(deltaR<minR) minR=deltaR;
            if(deltaR < matchingR && jetsDetMC.at(iAcc).at(ujetDetMC).pt() > maxpt) {
                maxpt = jetsDetMC.at(iAcc).at(ujetDetMC).pt();
                maxptIndex = ujetDetMC;
                deltaRMatch = deltaR;
                bfound = true;
                //cout << "found, detPt vs truePt: " << maxpt << " <> " << jets[iAcc][ujet].pt() << ", index: " << maxptIndex << endl;
            }
        }
        fhistos->fh_jetResponseDeltaRClosest->Fill(minR);
        if(bfound) {
            ptDetMC = jetsDetMC.at(iAcc).at(maxptIndex).pt();
            fhistos->fh_jetResponse->Fill(ptDetMC, ptTrue);
            fhistos->fh_jetResponse_ALICE->Fill(ptDetMC, ptTrue);
            fhistos->fh_jetResponseDeltaR->Fill(deltaRMatch);
            fhistos->fh_jetResponseDeltaPt->Fill((ptTrue-ptDetMC)/ptTrue);
            fhistos->fh_responseInfo->Fill("True jet has pair",1.0);
            bTrueJetMatch.at(ujet) = true;
            bDetJetMatch.at(maxptIndex) = true;
        } else {
            fhistos->fh_responseInfo->Fill("True jet has no pair",1.0);
        }
    }
    for (ujetDetMC = 0; ujetDetMC < NjetsDetMC; ujetDetMC++) { //Det MC jets
        if(!bDetJetMatch.at(ujetDetMC)) {
            fhistos->fh_responseInfo->Fill("Det jet has no pair",1.0);
        }
    }

    if(fDebug>8) {
        cout << "True jets: " << Njets << ", det jets: " << NjetsDetMC << endl;
        cout << "True jets:" << endl;
        for (ujet = 0; ujet < Njets; ujet++) { //True MC jets
            cout << "jet(E, pt, phi, eta) = " <<
                "(" << jets.at(iAcc).at(ujet).E()   << ", "
                << jets.at(iAcc).at(ujet).pt()  << ", "
                << jets.at(iAcc).at(ujet).phi() << ", "
                << jets.at(iAcc).at(ujet).eta() << ") "
                << (bTrueJetMatch.at(ujet) ? "match" : "no match") << endl;
        }
        cout << "Det jets:" <<  endl;
        for (ujetDetMC = 0; ujetDetMC < NjetsDetMC; ujetDetMC++) { //Det MC jets
            cout << "jet(E, pt, phi, eta) = " <<
                "(" << jetsDetMC.at(iAcc).at(ujetDetMC).E()   << ", "
                << jetsDetMC.at(iAcc).at(ujetDetMC).pt()  << ", "
                << jetsDetMC.at(iAcc).at(ujetDetMC).phi() << ", "
                << jetsDetMC.at(iAcc).at(ujetDetMC).eta() << ") "
                << (bDetJetMatch.at(ujetDetMC) ? "match" : "no match") << endl;
        }
    }

    vector<vector<vector<fastjet::PseudoJet>>> dijetsDetMC = anaDetMC->GetDijets();
    fastjet::PseudoJet dijetDetMC;

    //Dijet response without deltaphi cut.
    if(bHasDijet) {
        if(fDebug>8) {
            cout << "Normal dijet:" << endl;
            cout << "True dijet" << endl;
            cout << "leading jet(E, pt, phi, eta) =    " <<
                "(" << dijets.at(iAcc).at(0).at(0).E()   << ", "
                << dijets.at(iAcc).at(0).at(0).pt()  << ", "
                << dijets.at(iAcc).at(0).at(0).phi() << ", "
                << dijets.at(iAcc).at(0).at(0).eta() << ")" << endl;
            cout << "subleading jet(E, pt, phi, eta) = " <<
                "(" << dijets.at(iAcc).at(0).at(1).E()   << ", "
                << dijets.at(iAcc).at(0).at(1).pt()  << ", "
                << dijets.at(iAcc).at(0).at(1).phi() << ", "
                << dijets.at(iAcc).at(0).at(1).eta() << ")" << endl;
        }
        dijet = dijets.at(iAcc).at(0).at(0) + dijets.at(iAcc).at(0).at(1);
        if(anaDetMC->HasDijet()) {

            // Check that leading and subleading jets match.
            // It is also ok if leading and subleading jets change places in det level.
            deltaRLL = DeltaR(dijets.at(iAcc).at(0).at(0),dijetsDetMC.at(iAcc).at(0).at(0)); //leading, leading
            deltaRLS = DeltaR(dijets.at(iAcc).at(0).at(0),dijetsDetMC.at(iAcc).at(0).at(1)); //leading, subleading
            deltaRSS = DeltaR(dijets.at(iAcc).at(0).at(1),dijetsDetMC.at(iAcc).at(0).at(1)); //subleading, subleading
            deltaRSL = DeltaR(dijets.at(iAcc).at(0).at(1),dijetsDetMC.at(iAcc).at(0).at(0)); //subleading, leading

            if(deltaRLL < matchingR || deltaRLS < matchingR) bLeadingMatch    = true;
            if(deltaRSS < matchingR || deltaRSL < matchingR) bSubleadingMatch = true;

            if(fDebug>8) {
                cout << "Det dijet" << endl;
                cout << "leading jet(E, pt, phi, eta) =    " <<
                    "(" << dijetsDetMC.at(iAcc).at(0).at(0).E()   << ", "
                    << dijetsDetMC.at(iAcc).at(0).at(0).pt()  << ", "
                    << dijetsDetMC.at(iAcc).at(0).at(0).phi() << ", "
                    << dijetsDetMC.at(iAcc).at(0).at(0).eta() << ")" << endl;
                cout << "subleading jet(E, pt, phi, eta) = " <<
                    "(" << dijetsDetMC.at(iAcc).at(0).at(1).E()   << ", "
                    << dijetsDetMC.at(iAcc).at(0).at(1).pt()  << ", "
                    << dijetsDetMC.at(iAcc).at(0).at(1).phi() << ", "
                    << dijetsDetMC.at(iAcc).at(0).at(1).eta() << ")" << endl;
            }
            dijetDetMC = dijetsDetMC.at(iAcc).at(0).at(0) + dijetsDetMC.at(iAcc).at(0).at(1);
            if(bLeadingMatch && bSubleadingMatch) {
                fhistos->fh_dijetResponse->Fill(dijetDetMC.m(), dijet.m());
                fhistos->fh_dijetResponseTrunc->Fill(dijetDetMC.m(), dijet.m());
                fhistos->fh_responseInfo->Fill("Dijet match",1.0);
            } else {
                fhistos->fh_responseInfo->Fill("Dijet not match",1.0);
            }
        } else {
            fhistos->fh_responseInfo->Fill("Dijet det not found",1.0);
        }
    } else {
        if(anaDetMC->HasDijet()) {
            //dijetDetMC = dijetsDetMC.at(iAcc).at(0).at(0) + dijetsDetMC.at(iAcc).at(0).at(1);
            fhistos->fh_responseInfo->Fill("Dijet true not found",1.0);
        }
    }
    if(fDebug>8) cout << ((anaDetMC->HasDijet() && bLeadingMatch && bSubleadingMatch) ? "Dijet match" : "Dijet no match") << endl;

    // DeltaPhi cut dijet response.
    if(bHasDeltaPhiDijet) {
        if(fDebug>8) {
            cout << "Delta Phi dijet:" << endl;
            cout << "True dijet" << endl;
            cout << "leading jet(E, pt, phi, eta) =    " <<
                "(" << dijets.at(iAcc).at(1).at(0).E()   << ", "
                << dijets.at(iAcc).at(1).at(0).pt()  << ", "
                << dijets.at(iAcc).at(1).at(0).phi() << ", "
                << dijets.at(iAcc).at(1).at(0).eta() << ")" << endl;
            cout << "subleading jet(E, pt, phi, eta) = " <<
                "(" << dijets.at(iAcc).at(1).at(1).E()   << ", "
                << dijets.at(iAcc).at(1).at(1).pt()  << ", "
                << dijets.at(iAcc).at(1).at(1).phi() << ", "
                << dijets.at(iAcc).at(1).at(1).eta() << ")" << endl;
        }
        dijet = dijets.at(iAcc).at(1).at(0) + dijets.at(iAcc).at(1).at(1);
        if(anaDetMC->HasDeltaPhiDijet()) {
        
            // Leading jets are the same as without delta phi cut.
            deltaRSS = DeltaR(dijets.at(iAcc).at(1).at(1),dijetsDetMC.at(iAcc).at(1).at(1)); //subleading, subleading
            deltaRSL = DeltaR(dijets.at(iAcc).at(1).at(1),dijetsDetMC.at(iAcc).at(1).at(0)); //subleading, leading
            if(deltaRSS < matchingR || deltaRSL < matchingR) bSubleadingMatchDeltaPhi = true;

            if(fDebug>8) {
                cout << "Det dijet" << endl;
                cout << "leading jet(E, pt, phi, eta) =    " <<
                    "(" << dijetsDetMC.at(iAcc).at(1).at(0).E()   << ", "
                    << dijetsDetMC.at(iAcc).at(1).at(0).pt()  << ", "
                    << dijetsDetMC.at(iAcc).at(1).at(0).phi() << ", "
                    << dijetsDetMC.at(iAcc).at(1).at(0).eta() << ")" << endl;
                cout << "subleading jet(E, pt, phi, eta) = " <<
                    "(" << dijetsDetMC.at(iAcc).at(1).at(1).E()   << ", "
                    << dijetsDetMC.at(iAcc).at(1).at(1).pt()  << ", "
                    << dijetsDetMC.at(iAcc).at(1).at(1).phi() << ", "
                    << dijetsDetMC.at(iAcc).at(1).at(1).eta() << ")" << endl;
            }
            dijetDetMC = dijetsDetMC.at(iAcc).at(1).at(0) + dijetsDetMC.at(iAcc).at(1).at(1);
            // Check subleading jet match.

            if(bLeadingMatch && bSubleadingMatchDeltaPhi) {
                fhistos->fh_dijetResponseDeltaPhiCut->Fill(dijetDetMC.m(), dijet.m());
                fhistos->fh_dijetResponseDeltaPhiCutTrunc->Fill(dijetDetMC.m(), dijet.m());
                fhistos->fh_responseInfo->Fill("Dijet DPhi match",1.0);
            } else {
                fhistos->fh_responseInfo->Fill("Dijet DPhi not match",1.0);
            }
        } else {
            fhistos->fh_responseInfo->Fill("Dijet DPhi det not found",1.0);
        }
    } else {
        if(anaDetMC->HasDeltaPhiDijet()) {
            //dijetDetMC = dijetsDetMC.at(iAcc).at(1).at(0) + dijetsDetMC.at(iAcc).at(1).at(1);
            fhistos->fh_responseInfo->Fill("Dijet DPhi true not found",1.0);
        }
    }
    if(fDebug>8) cout << ((anaDetMC->HasDeltaPhiDijet() && bLeadingMatch && bSubleadingMatchDeltaPhi) ? "Delta Phi Dijet match" : "Delta Phi Dijet no match") << endl;
    return;
}

void AliJCDijetAna::ResetObjects() {
    chparticles.clear();
    ktchparticles.clear();
    for(int i=0; i<jetClassesSize; i++) jets.at(i).clear();
    rawJets.clear();
    rawKtJets.clear();
    rhoEstJets.clear();

    phi=0, eta=0, pt=0, pt2=0, rho=0, rhom=0, area=0, mjj=0, ptpair=0, dPhi=0;
    leadingTrackOverThreshold = false;
    bHasDijet = false;
    bHasDeltaPhiDijet = false;
    bHasDeltaPhiSubLeadJet = false;
    bEvtHasAreaInfo = false;
    return;
}

double AliJCDijetAna::DeltaR(fastjet::PseudoJet jet1, fastjet::PseudoJet jet2) {
    return DeltaR(jet1.eta(), jet2.eta(), jet1.phi(), jet2.phi());
}

double AliJCDijetAna::DeltaR(double eta1, double eta2, double phi1, double phi2) {
    double Deta = eta1 - eta2;
    double Dphi = TMath::Abs(phi1 - phi2);
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
    histos->fh_info->Fill("tracking ineff",ftrackingIneff);

    histos->fh_eventSel->Fill("events wo/ cuts",0.0);
    histos->fh_eventSel->Fill("catalyst entry ok",0.0);
    histos->fh_eventSel->Fill("catalyst ok",0.0);
    histos->fh_eventSel->Fill("vertex2013pA ok",0.0);
    histos->fh_eventSel->Fill("pileupSPD ok",0.0);
    histos->fh_eventSel->Fill("utils pileupSPD ok",0.0);
    histos->fh_eventSel->Fill("events",0.0);

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
        histos->fh_events[iBin]->Fill("pt_hard bin cuts",0.0);
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
    return;
}
#endif
