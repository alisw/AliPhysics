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
    matchingR(0),
    ftrackingIneffHisto(NULL)
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
    selectorNoGhosts(),
    selectorNoGhostsAllButTwo(),
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
    matchingR(obj.matchingR),
    ftrackingIneffHisto(obj.ftrackingIneffHisto)
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
    selectorNoGhosts(obj.selectorNoGhosts),
    selectorNoGhostsAllButTwo(obj.selectorNoGhostsAllButTwo),
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
                                double lParticlePtCutMax,
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
                                double ltrackingIneff,
                                TH1D*  ltrackingIneffHistogram,
                                bool   luseCrho,
                                bool   lThisIsTrueMC){
    fDebug = lDebug;
    fParticleEtaCut = lParticleEtaCut;
    fParticlePtCut = lParticlePtCut;
    fParticlePtCutMax = lParticlePtCutMax;
    fusePionMass = lusePionMass;
    fUseDeltaPhiBGSubtr = luseDeltaPhiBGSubtr;
    fConstituentCut = lConstituentCut;
    fLeadingJetCut = lLeadingJetCut;
    fSubleadingJetCut = lSubleadingJetCut;
    fDeltaPhiCut = lDeltaPhiCut;
    ftrackingIneff = ltrackingIneff;
    ftrackingIneffHisto = ltrackingIneffHistogram;
    if(ftrackingIneff<0.0 && ftrackingIneffHisto==0) {
        ::Error("AliJCDijetAna","Trying to use pt dependend tracking inefficiency, but missing the histogram!");
    }
    bUseCrho = luseCrho;
    bThisIsTrueMC = lThisIsTrueMC;

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
    fPythiaSigma = 0; // Default value
    fPythiaTrial = 0; // Default value

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
        bHasDijet.push_back(false);
        bHasDeltaPhiDijet.push_back(false);
    }
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
    selectorNoGhosts  = !fastjet::SelectorIsPureGhost();
    selectorNoGhostsAllButTwo  = selectorAllButTwo * selectorNoGhosts; // Here right selector is applied first, then the left one.
    bge = fastjet::JetMedianBackgroundEstimator(selectorEta, jet_def_bge, area_def_bge);

    areaCut=0.6*TMath::Pi()*0.4*0.4;

    return;
}

//______________________________________________________________________________
int AliJCDijetAna::CalculateJets(TClonesArray *inList, AliJCDijetHistos *fhistos, int lCBin, double hisWeight){

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
            if(ftrackingIneff!=0.0) {
                if(ftrackingIneff<0.0) {
                    //Note: Histograms have been defined as portion of particles recovered whereas
                    //ftrackingIneff is defined as the portion of tracks lost.
                    ftrackingIneffTemp = 1.0-ftrackingIneffHisto->GetBinContent(ftrackingIneffHisto->FindBin(pt));
                } else {
                    ftrackingIneffTemp = ftrackingIneff;
                }
                if(randomGenerator->Uniform(0.0,1.0) < ftrackingIneffTemp) {
                    fhistos->fh_events[lCBin]->Fill("trackingIneff rejection",1.0);
                    continue;
                }
                fhistos->fh_events[lCBin]->Fill("trackingIneff accept",1.0);
            }
            phi = trk->Phi() > TMath::Pi() ? trk->Phi()-2*TMath::Pi() : trk->Phi();
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

    tempJets   = fastjet::sorted_by_pt(cs->inclusive_jets(MinJetPt)); // APPLY Min pt cut for jet

    if(bThisIsTrueMC) {
        rawJets = tempJets;
    } else {
        for (ujet = 0; ujet < tempJets.size(); ujet++) {
            if(tempJets.at(ujet).area()>areaCut) {  //Check first that jet area is large enough
                fhistos->fh_events[lCBin]->Fill("area cut ok",1.0);
                ucount=0;
                for(uconst=0;uconst<tempJets.at(ujet).constituents().size(); uconst++) {
                    if(tempJets.at(ujet).constituents().at(uconst).pt() > fParticlePtCutMax ) { // Check if jet has 100 GeV constituents
                        ucount++;
                    }
                }
                if(ucount==0) {
                    fhistos->fh_events[lCBin]->Fill("100GeV const cut ok",1.0);
                    rawJets.push_back(tempJets.at(ujet));
                } else {
                    fhistos->fh_events[lCBin]->Fill("100GeV const cut fails",1.0);
                }
            } else { //if area cut fails
                fhistos->fh_events[lCBin]->Fill("area cut fails",1.0);
            }
        }
    }

    // For MC runs: If we find jets with over 4 times pt_hard bin, reject the event.
    if( fptHardBin!=0 && rawJets.size()>0 ) {
        fhistos->fh_ptHard[lCBin]->Fill(fptHardBin,hisWeight);
        fhistos->fh_maxJetptOverPtHard[lCBin]->Fill(rawJets.at(0).pt()/fptHardBin,hisWeight);
        if( rawJets.at(0).pt() > fptHardBin*4 ) {
            fhistos->fh_events[lCBin]->Fill("pt_hard bin cuts",1.0);
            return -1;
        }
    }
    if(fPythiaSigma!=0) fhistos->fh_pythiaSigma[lCBin]->Fill(fPythiaSigma,hisWeight);
    if(fPythiaTrial!=0) fhistos->fh_pythiaTrial[lCBin]->Fill(fPythiaTrial,hisWeight);
    fhistos->fh_hisWeight[lCBin]->Fill(hisWeight);
    fhistos->fh_randConeEtaPhi[lCBin]->Fill(randConeEta,randConePhi,hisWeight);

    rawKtJets = fastjet::sorted_by_pt(selectorEta(cs_bge->inclusive_jets(0.0))); // APPLY Min pt cut for jet

    fhistos->fh_events[lCBin]->Fill("particles",noTracks);
    for (utrack = 0; utrack < chparticles.size(); utrack++) {
        fhistos->fh_events[lCBin]->Fill("acc. particles",1.0);
        pt = chparticles.at(utrack).pt();
        eta = chparticles.at(utrack).eta();
        phi = chparticles.at(utrack).phi();
        fhistos->fh_eta[lCBin]->Fill(eta,hisWeight);
        fhistos->fh_phi[lCBin]->Fill(phi-TMath::Pi(),hisWeight);
        fhistos->fh_etaPhi[lCBin]->Fill(eta,phi-TMath::Pi(),hisWeight);
        fhistos->fh_pt[lCBin]->Fill(pt,hisWeight);
        if(eta>0.0) fhistos->fh_ptPosEta[lCBin]->Fill(pt,hisWeight);
        else        fhistos->fh_ptNegEta[lCBin]->Fill(pt,hisWeight);
    }
    fhistos->fh_nch[lCBin]->Fill(chparticles.size(),hisWeight);
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
            //If using C factor do not include ghost-only jets.
            if(bUseCrho) {
                if(!rawKtJets[uktjet].is_pure_ghost()) {
                    rhoEstJets.push_back(rawKtJets.at(uktjet)); 
                }
            } else {
                rhoEstJets.push_back(rawKtJets.at(uktjet)); 
            }
        }
    } else {
        if(bUseCrho) {
            //If using C factor do not include ghost-only jets.
            rhoEstJets = selectorNoGhostsAllButTwo(rawKtJets);
        } else {
            rhoEstJets = selectorAllButTwo(rawKtJets);
        }
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
        //This will turn on the "new CMS" method, where empty areas are taken
        //into account specifically.
        if(bUseCrho) {
            double Covered=0.0;
            double Total=0.0;
            //For C factor, all kt jets are included.
            for (uktjet = 0; uktjet < rawKtJets.size(); uktjet++) { // First jet is already skipped here.
                if(rawKtJets[uktjet].is_pure_ghost()) {
                    Total+=rawKtJets[uktjet].area();
                } else {
                    Total+=rawKtJets[uktjet].area();
                    Covered+=rawKtJets[uktjet].area();
                }
            }
            Covered/=Total;
            rho=rho*Covered;
            fhistos->fh_coveredRatio[lCBin]->Fill(Covered,hisWeight);
        }
    }


    fhistos->fh_rho[lCBin]->Fill(rho,hisWeight);
    fhistos->fh_rhoLin[lCBin]->Fill(rho,hisWeight);
    fhistos->fh_rhom[lCBin]->Fill(rhom,hisWeight);
    fhistos->fh_rhomLin[lCBin]->Fill(rhom,hisWeight);
    if(fDebug > 9) std::cout << "Testing: Rho_M = " << rhom << ", has_rho_m() = " << bge.has_rho_m() << std::endl;

    // This is the delta-phi distribution: delta-pt = sum_{tracks in rand cone}( pt ) - rho*pi*R^2
    fDeltaPt = randConePt - rho*TMath::Pi()*fJetCone*fJetCone;
    fhistos->fh_deltaPt[lCBin]->Fill(fDeltaPt,hisWeight);

    return 0;
}

// If user wants to add jets from outside this function can be used.
void AliJCDijetAna::SetJets(vector<fastjet::PseudoJet> jetsOutside) {
    ResetObjects();
    rawJets = fastjet::sorted_by_pt(jetsOutside);
    return;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Loop over jets and fill various histos 
// Note: Jets in rawJets are assumed to be sorted by pt, as is done
//       in the previous functions.
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliJCDijetAna::FillJetsDijets(AliJCDijetHistos *fhistos, int lCBin, double hisWeight) {
    int iAccJetCounter = 0;
    bool bHasHighPtJet = false;

    vector<unsigned> emptyRawJetIndexVector;
    vector<vector<unsigned>> uRawJetIndex;
    for(int i=0; i<jetClassesSize; i++) {
        uRawJetIndex.push_back(emptyRawJetIndexVector);
    }

    // anti-kt jets:
    for (ujet = 0; ujet < rawJets.size(); ujet++) {
        eta = rawJets.at(ujet).eta();
        fhistos->fh_events[lCBin]->Fill("jets",1.0);
        // anti-kt-jet eta cut
        pt = rawJets.at(ujet).pt();
        phi = rawJets.at(ujet).phi();
        jetAreaVector = rawJets.at(ujet).area_4vector();
        jet_bgSubtracted = fastjet::PseudoJet(rawJets.at(ujet).px() -        rho * jetAreaVector.px(),
                                              rawJets.at(ujet).py() -        rho * jetAreaVector.py(),
                                              rawJets.at(ujet).pz() - (rho+rhom) * jetAreaVector.pz(),
                                              rawJets.at(ujet).E()  - (rho+rhom) * jetAreaVector.E());
        if(TMath::Abs(eta) < etaMaxCutForJet) {
            iAccJetCounter++;
            fhistos->fh_events[lCBin]->Fill("acc. jets",1.0);
            if(iAccJetCounter==1) bHasHighPtJet = pt > fLeadingJetCut;
            area = rawJets.at(ujet).area();
            fhistos->fh_jetEta[lCBin][iAcc]->Fill(eta,hisWeight);  
            fhistos->fh_jetPhi[lCBin][iAcc]->Fill(phi - TMath::Pi(),hisWeight); //Pseudojet.phi range 0-2pi
            fhistos->fh_jetEtaPhi[lCBin][iAcc]->Fill(eta,phi - TMath::Pi(),hisWeight);
            fhistos->fh_jetPt[lCBin][iAcc]->Fill(pt,hisWeight);
            fhistos->fh_jetPt_ALICE[lCBin][iAcc]->Fill(pt,hisWeight);
            fhistos->fh_jetPtTransBGSub[lCBin][iBGSubtr]->Fill(pt - rho*area,hisWeight); //Note, this has iBGSubtr
            fhistos->fh_jetArea[lCBin][iAcc]->Fill(area,hisWeight);
            fhistos->fh_jetAreaRho[lCBin][iAcc]->Fill(area*rho,hisWeight);

            jets.at(iAcc).push_back(rawJets.at(ujet));

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
                fhistos->fh_jetEta[lCBin][iConstCut]->Fill(eta,hisWeight);
                fhistos->fh_jetPhi[lCBin][iConstCut]->Fill(phi - TMath::Pi(),hisWeight);
                fhistos->fh_jetEtaPhi[lCBin][iConstCut]->Fill(eta,phi - TMath::Pi(),hisWeight);
                fhistos->fh_jetPt[lCBin][iConstCut]->Fill(pt,hisWeight);
                fhistos->fh_jetPt_ALICE[lCBin][iConstCut]->Fill(pt,hisWeight);
                fhistos->fh_jetPtTransBGSub[lCBin][iBGSubtrConstCut]->Fill(pt - rho*area,hisWeight); //Note, this has iBGSubtrConstCut
                fhistos->fh_jetArea[lCBin][iConstCut]->Fill(area,hisWeight);
                fhistos->fh_jetAreaRho[lCBin][iConstCut]->Fill(area*rho,hisWeight);

                jets.at(iConstCut).push_back(rawJets.at(ujet));
            }

            fhistos->fh_events[lCBin]->Fill("bg. subtr. jets (common eta)",1.0);
            pt2 = jet_bgSubtracted.pt();
            phi = jet_bgSubtracted.phi();
            eta = jet_bgSubtracted.eta();
            if(iAccJetCounter==1 && pt>fLeadingJetCut && pt2<=fLeadingJetCut)       fhistos->fh_events[lCBin]->Fill("leading jet drop",1.0);
            if(iAccJetCounter==2 && pt>fSubleadingJetCut && pt2<=fSubleadingJetCut) fhistos->fh_events[lCBin]->Fill("subleading jet drop",1.0);
            fhistos->fh_jetEta[lCBin][iBGSubtrCommonEta]->Fill(eta,hisWeight);
            fhistos->fh_jetPhi[lCBin][iBGSubtrCommonEta]->Fill(phi - TMath::Pi(),hisWeight);
            fhistos->fh_jetEtaPhi[lCBin][iBGSubtrCommonEta]->Fill(eta,phi - TMath::Pi(),hisWeight);
            fhistos->fh_jetPt[lCBin][iBGSubtrCommonEta]->Fill(pt2,hisWeight);
            fhistos->fh_jetPt_ALICE[lCBin][iBGSubtrCommonEta]->Fill(pt2,hisWeight);
            fhistos->fh_jetArea[lCBin][iBGSubtrCommonEta]->Fill(area,hisWeight); // Assuming bg subtracted jet has the same area.
            fhistos->fh_jetAreaRho[lCBin][iBGSubtrCommonEta]->Fill(area*rho,hisWeight);

            uRawJetIndex.at(iBGSubtrCommonEta).push_back(ujet);
            jets.at(iBGSubtrCommonEta).push_back(jet_bgSubtracted);

            if(leadingTrackOverThreshold) {
                fhistos->fh_events[lCBin]->Fill("bg. subtr. const. cut jets (common eta)",1.0);
                fhistos->fh_jetEta[lCBin][iBGSubtrConstCutCommonEta]->Fill(eta,hisWeight);
                fhistos->fh_jetPhi[lCBin][iBGSubtrConstCutCommonEta]->Fill(phi - TMath::Pi(),hisWeight);
                fhistos->fh_jetEtaPhi[lCBin][iBGSubtrConstCutCommonEta]->Fill(eta,phi - TMath::Pi(),hisWeight);
                fhistos->fh_jetPt[lCBin][iBGSubtrConstCutCommonEta]->Fill(pt2,hisWeight);
                fhistos->fh_jetPt_ALICE[lCBin][iBGSubtrConstCutCommonEta]->Fill(pt2,hisWeight);
                fhistos->fh_jetArea[lCBin][iBGSubtrConstCutCommonEta]->Fill(area,hisWeight);
                fhistos->fh_jetAreaRho[lCBin][iBGSubtrConstCutCommonEta]->Fill(area*rho,hisWeight);

                uRawJetIndex.at(iBGSubtrConstCutCommonEta).push_back(ujet);
                jets.at(iBGSubtrConstCutCommonEta).push_back(jet_bgSubtracted);
            }

            if(TMath::Abs(eta) < etaMaxCutForJet) {
                fhistos->fh_jetBGSubtrDeltaR[lCBin]->Fill(DeltaR(rawJets.at(ujet), jet_bgSubtracted),hisWeight);
                fhistos->fh_events[lCBin]->Fill("bg. subtr. jets",1.0);
                if(iAccJetCounter==1 && pt>fLeadingJetCut && pt2<=fLeadingJetCut)       fhistos->fh_events[lCBin]->Fill("leading jet drop",1.0);
                if(iAccJetCounter==2 && pt>fSubleadingJetCut && pt2<=fSubleadingJetCut) fhistos->fh_events[lCBin]->Fill("subleading jet drop",1.0);
                fhistos->fh_jetEta[lCBin][iBGSubtr]->Fill(eta,hisWeight);
                fhistos->fh_jetPhi[lCBin][iBGSubtr]->Fill(phi - TMath::Pi(),hisWeight);
                fhistos->fh_jetEtaPhi[lCBin][iBGSubtr]->Fill(eta,phi - TMath::Pi(),hisWeight);
                fhistos->fh_jetPt[lCBin][iBGSubtr]->Fill(pt2,hisWeight);
                fhistos->fh_jetPt_ALICE[lCBin][iBGSubtr]->Fill(pt2,hisWeight);
                fhistos->fh_jetArea[lCBin][iBGSubtr]->Fill(area,hisWeight); // Assuming bg subtracted jet has the same area.
                fhistos->fh_jetAreaRho[lCBin][iBGSubtr]->Fill(area*rho,hisWeight);

                uRawJetIndex.at(iBGSubtr).push_back(ujet);
                jets.at(iBGSubtr).push_back(jet_bgSubtracted);
                uRawJetIndex.at(iBGSubtrCutsRaw).push_back(ujet);
                jets.at(iBGSubtrCutsRaw).push_back(jet_bgSubtracted);

                if(leadingTrackOverThreshold) {
                    fhistos->fh_events[lCBin]->Fill("bg. subtr. const. cut jets",1.0);
                    fhistos->fh_jetEta[lCBin][iBGSubtrConstCut]->Fill(eta,hisWeight);
                    fhistos->fh_jetPhi[lCBin][iBGSubtrConstCut]->Fill(phi - TMath::Pi(),hisWeight);
                    fhistos->fh_jetEtaPhi[lCBin][iBGSubtrConstCut]->Fill(eta,phi - TMath::Pi(),hisWeight);
                    fhistos->fh_jetPt[lCBin][iBGSubtrConstCut]->Fill(pt2,hisWeight);
                    fhistos->fh_jetPt_ALICE[lCBin][iBGSubtrConstCut]->Fill(pt2,hisWeight);
                    fhistos->fh_jetArea[lCBin][iBGSubtrConstCut]->Fill(area,hisWeight);
                    fhistos->fh_jetAreaRho[lCBin][iBGSubtrConstCut]->Fill(area*rho,hisWeight);

                    uRawJetIndex.at(iBGSubtrConstCut).push_back(ujet);
                    jets.at(iBGSubtrConstCut).push_back(jet_bgSubtracted);
                    uRawJetIndex.at(iBGSubtrConstCutCutsRaw).push_back(ujet);
                    jets.at(iBGSubtrConstCutCutsRaw).push_back(jet_bgSubtracted);
                }
            }
        }

        // Check eta acceptance also for bg subtracted jets.
        eta = jet_bgSubtracted.eta();
        if(TMath::Abs(eta) < etaMaxCutForJet) {
            fhistos->fh_events[lCBin]->Fill("bg. subtr. jets (bg eta)",1.0);
            pt2 = jet_bgSubtracted.pt();
            phi = jet_bgSubtracted.phi();
            fhistos->fh_jetEta[lCBin][iBGSubtrIndEta]->Fill(eta,hisWeight);
            fhistos->fh_jetPhi[lCBin][iBGSubtrIndEta]->Fill(phi - TMath::Pi(),hisWeight);
            fhistos->fh_jetEtaPhi[lCBin][iBGSubtrIndEta]->Fill(eta,phi - TMath::Pi(),hisWeight);
            fhistos->fh_jetPt[lCBin][iBGSubtrIndEta]->Fill(pt2,hisWeight);
            fhistos->fh_jetPt_ALICE[lCBin][iBGSubtrIndEta]->Fill(pt2,hisWeight);
            fhistos->fh_jetArea[lCBin][iBGSubtrIndEta]->Fill(area,hisWeight); // Assuming bg subtracted jet has the same area.
            fhistos->fh_jetAreaRho[lCBin][iBGSubtrIndEta]->Fill(area*rho,hisWeight);

            uRawJetIndex.at(iBGSubtrIndEta).push_back(ujet);
            jets.at(iBGSubtrIndEta).push_back(jet_bgSubtracted);

            if(leadingTrackOverThreshold) {
                fhistos->fh_events[lCBin]->Fill("bg. subtr. const. cut jets (bg eta)",1.0);
                fhistos->fh_jetEta[lCBin][iBGSubtrConstCutIndEta]->Fill(eta,hisWeight);
                fhistos->fh_jetPhi[lCBin][iBGSubtrConstCutIndEta]->Fill(phi - TMath::Pi(),hisWeight);
                fhistos->fh_jetEtaPhi[lCBin][iBGSubtrConstCutIndEta]->Fill(eta,phi - TMath::Pi(),hisWeight);
                fhistos->fh_jetPt[lCBin][iBGSubtrConstCutIndEta]->Fill(pt2,hisWeight);
                fhistos->fh_jetPt_ALICE[lCBin][iBGSubtrConstCutIndEta]->Fill(pt2,hisWeight);
                fhistos->fh_jetArea[lCBin][iBGSubtrConstCutIndEta]->Fill(area,hisWeight);
                fhistos->fh_jetAreaRho[lCBin][iBGSubtrConstCutIndEta]->Fill(area*rho,hisWeight);

                uRawJetIndex.at(iBGSubtrConstCutIndEta).push_back(ujet);
                jets.at(iBGSubtrConstCutIndEta).push_back(jet_bgSubtracted);
            }
        }
    }//end of the anti-kt-jet loop
    //============================================================================================
    //cout << "Jets: " << jets[iAcc].size() << endl;
    //============================================================================================
    fhistos->fh_jetN[lCBin]->Fill(iAccJetCounter,hisWeight);
    if(bHasHighPtJet) {
        fhistos->fh_rhoHighPt[lCBin]->Fill(rho,hisWeight);
        fhistos->fh_rhoLinHighPt[lCBin]->Fill(rho,hisWeight);
        fhistos->fh_rhomHighPt[lCBin]->Fill(rhom,hisWeight);
        fhistos->fh_rhomLinHighPt[lCBin]->Fill(rhom,hisWeight);
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
            area = rawKtJets.at(uktjet).area();
            fhistos->fh_jetEta[lCBin][iktJets]->Fill(eta,hisWeight);  
            fhistos->fh_jetPhi[lCBin][iktJets]->Fill(phi - TMath::Pi(),hisWeight); //Pseudojet.phi range 0-2pi
            fhistos->fh_jetEtaPhi[lCBin][iktJets]->Fill(eta,phi - TMath::Pi(),hisWeight);
            fhistos->fh_jetPt[lCBin][iktJets]->Fill(pt,hisWeight);
            fhistos->fh_jetPt_ALICE[lCBin][iktJets]->Fill(pt,hisWeight);
            fhistos->fh_jetArea[lCBin][iktJets]->Fill(area,hisWeight);
            fhistos->fh_jetAreaRho[lCBin][iktJets]->Fill(area*rho,hisWeight);
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
        fhistos->fh_jetDeltaRMin[lCBin][udijet]->Fill(deltaRMin,hisWeight);
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
    double jetCutCompare;
    double maxPt;
    double secondMaxPt;
    int iDijetBin;
    std::vector<unsigned> uLeading(jetClassesSize,10000);
    std::vector<unsigned> uSubleading(jetClassesSize,10000);
    std::vector<unsigned> uSubleadingDeltaPhi(jetClassesSize,10000);
    for(udijet=0; udijet < jetClassesSize; udijet++) {
        if(jets.at(udijet).size()>1) {
            // After modifications: we don't want to mess up the order.
            // Leading and subleading will be figured out.
            //jets.at(udijet) = fastjet::sorted_by_pt(jets.at(udijet)); // Sort in case of bg subtr messed up the order.
            fhistos->fh_events[lCBin]->Fill(Form("%s dijets",sDijetTypes[udijet].Data()),1.0);

            // Start by forming dijets. First without deltaPhiCut and then with the cut.
            maxPt=0;
            secondMaxPt=0;
            for (int ijet = jets.at(udijet).size()-1; ijet >= 0 ; ijet--) {
                pt = jets.at(udijet).at(ijet).pt();
                if(pt>maxPt) {
                    secondMaxPt=maxPt;
                    maxPt=pt;
                    uSubleading.at(udijet)=uLeading.at(udijet); //Previous highest is now subleading
                    uLeading.at(udijet)=ijet;                //New highest becomes leading
                } else if(pt>secondMaxPt) {  //If new jet is not highest, it still might higher than second highest.
                    secondMaxPt=pt;
                    uSubleading.at(udijet)=ijet;
                }
            }
            
            // No deltaPhi cut for these jets.
            if(udijet==iBGSubtrCutsRaw || udijet==iBGSubtrConstCutCutsRaw) { //For these sets, the cut is made with raw jets.
                if(rawJets.size()>uLeading.at(udijet)) {
                    //Ensure we still check the same jet which was selected as leading in 'udijet' set
                    jetCutCompare = rawJets.at(uRawJetIndex.at(udijet).at(uLeading.at(udijet))).pt();  
                } else continue;
            } else {
                jetCutCompare = jets.at(udijet).at(uLeading.at(udijet)).pt();
            }
            if(jetCutCompare < fLeadingJetCut) continue;
            //if(udijet==iAcc)cout << "maxpt: " << maxPt << ", secondMax: " << secondMaxPt << endl;
            fhistos->fh_events[lCBin]->Fill(Form("%s dijets leading cut",sDijetTypes[udijet].Data()),1.0);
            dijets.at(udijet).at(0).at(0) = jets.at(udijet).at(uLeading.at(udijet));
            dijets.at(udijet).at(0).at(1) = jets.at(udijet).at(uSubleading.at(udijet));

            // Here we check deltaPhi cut
            // The leading is always the same as without deltaPhi cut
            dijets.at(udijet).at(1).at(0) = jets.at(udijet).at(uLeading.at(udijet));
            bHasDeltaPhiSubLeadJet = false;
            maxPt=0;
            for (ujet = 0; ujet < jets.at(udijet).size(); ujet++) {
                if(ujet==uLeading.at(udijet)) continue; //The known leading jet will be skipped
                if(CheckDeltaPhi(dijets.at(udijet).at(1).at(0), jets.at(udijet).at(ujet), TMath::Pi()/fDeltaPhiCut)) {
                    pt = jets.at(udijet).at(ujet).pt();
                    if(pt>maxPt) {
                        maxPt=pt;
                        uSubleadingDeltaPhi.at(udijet) = ujet;
                    }
                }
            }
            if(uSubleadingDeltaPhi.at(udijet)!=10000) {
                dijets.at(udijet).at(1).at(1) = jets.at(udijet).at(uSubleadingDeltaPhi.at(udijet));
                bHasDeltaPhiSubLeadJet = true;
                //if(udijet==iAcc)cout << "secondMax(deltaPhi): " << maxPt << endl;
            }

            // Analysis for dijet without deltaPhi cut.
            if(udijet==iBGSubtrCutsRaw || udijet==iBGSubtrConstCutCutsRaw) {
                if(rawJets.size()>uSubleading.at(udijet)) {
                    jetCutCompare = rawJets.at(uRawJetIndex.at(udijet).at(uSubleading.at(udijet))).pt();
                } else continue;
            } else {
                jetCutCompare = dijets.at(udijet).at(0).at(1).pt();
            }
            if(jetCutCompare>fSubleadingJetCut) {
                fhistos->fh_events[lCBin]->Fill(Form("%s acc. dijets",sDijetTypes[udijet].Data()),1.0);
                bHasDijet.at(udijet) = true;
                dijet = dijets.at(udijet).at(0).at(0) + dijets.at(udijet).at(0).at(1);
                fhistos->fh_jetPtLeadSublead_ALICE[lCBin][udijet]->Fill(dijets.at(udijet).at(0).at(0).pt(),hisWeight);
                fhistos->fh_jetPtLeadSublead_ALICE[lCBin][udijet]->Fill(dijets.at(udijet).at(0).at(1).pt(),hisWeight);
                mjj = dijet.m();
                ptpair = dijet.pt();
                iDijetBin = fhistos->GetDijetMClass(mjj);
                fhistos->fh_jetPtLeadSubleadMBin_ALICE[lCBin][udijet][iDijetBin]->Fill(dijets.at(udijet).at(0).at(0).pt(),hisWeight);
                fhistos->fh_jetPtLeadSubleadMBin_ALICE[lCBin][udijet][iDijetBin]->Fill(dijets.at(udijet).at(0).at(1).pt(),hisWeight);
                fhistos->fh_dijetInvM[lCBin][udijet]->Fill(mjj,hisWeight);
                fhistos->fh_dijetInvMLin[lCBin][udijet]->Fill(mjj,hisWeight);
                fhistos->fh_dijetInvMTrunc[lCBin][udijet]->Fill(mjj,hisWeight);
                fhistos->fh_dijetInvMTrunc2[lCBin][udijet]->Fill(mjj,hisWeight);
                fhistos->fh_dijetPtPair[lCBin][udijet]->Fill(ptpair,hisWeight);
                dPhi = GetDeltaPhi(dijets.at(udijet).at(0).at(0), dijets.at(udijet).at(0).at(1));
                fhistos->fh_dijetDeltaPhi[lCBin][udijet]->Fill(dPhi,hisWeight);
                fhistos->fh_dijetCosDeltaPhi[lCBin][udijet]->Fill(TMath::Cos(dPhi),hisWeight);
                fhistos->fh_dijetDeltaEta[lCBin][udijet]->Fill(dijets.at(udijet).at(0).at(0).eta()-dijets.at(udijet).at(0).at(1).eta(),hisWeight);
                fhistos->fh_dijetCoshDeltaEta[lCBin][udijet]->Fill(TMath::CosH(dijets.at(udijet).at(0).at(0).eta()-dijets.at(udijet).at(0).at(1).eta()),hisWeight);
                fhistos->fh_dijetSqrt2pt12[lCBin][udijet]->Fill(TMath::Sqrt(2*dijets.at(udijet).at(0).at(0).pt()*dijets.at(udijet).at(0).at(1).pt()),hisWeight);
                fhistos->fh_dijetSqrtGeometry[lCBin][udijet]->Fill(TMath::Sqrt(TMath::CosH(dijets.at(udijet).at(0).at(0).eta()-dijets.at(udijet).at(0).at(1).eta())-TMath::Cos(dPhi)),hisWeight);
            }

            // Analysis for dijet with deltaPhi cut.
            if(bHasDeltaPhiSubLeadJet) {
                if(udijet==iBGSubtrCutsRaw || udijet==iBGSubtrConstCutCutsRaw) {
                    if(rawJets.size()>uSubleadingDeltaPhi.at(udijet)) {
                        jetCutCompare = rawJets.at(uRawJetIndex.at(udijet).at(uSubleadingDeltaPhi.at(udijet))).pt();
                    } else continue;
                } else {
                    jetCutCompare = dijets.at(udijet).at(1).at(1).pt();
                }
                if(jetCutCompare>fSubleadingJetCut) {
                    //if(udijet==iAcc)cout << "DeltaPhi cut successfull" << endl;
                    fhistos->fh_events[lCBin]->Fill(Form("%s deltaphi cut dijets",sDijetTypes[udijet].Data()),1.0);
                    bHasDeltaPhiDijet.at(udijet) = true;
                    dijet = dijets.at(udijet).at(1).at(0) + dijets.at(udijet).at(1).at(1);
                    fhistos->fh_jetPtLeadSubleadDeltaPhi_ALICE[lCBin][udijet]->Fill(dijets.at(udijet).at(1).at(0).pt(),hisWeight);
                    fhistos->fh_jetPtLeadSubleadDeltaPhi_ALICE[lCBin][udijet]->Fill(dijets.at(udijet).at(1).at(1).pt(),hisWeight);
                    mjj = dijet.m();
                    ptpair = dijet.pt();
                    iDijetBin = fhistos->GetDijetMClass(mjj);
                    fhistos->fh_jetPtLeadSubleadDeltaPhiMBin_ALICE[lCBin][udijet][iDijetBin]->Fill(dijets.at(udijet).at(0).at(0).pt(),hisWeight);
                    fhistos->fh_jetPtLeadSubleadDeltaPhiMBin_ALICE[lCBin][udijet][iDijetBin]->Fill(dijets.at(udijet).at(0).at(1).pt(),hisWeight);
                    fhistos->fh_dijetInvMDeltaPhiCut[lCBin][udijet]->Fill(mjj,hisWeight);
                    fhistos->fh_dijetInvMDeltaPhiCutLin[lCBin][udijet]->Fill(mjj,hisWeight);
                    fhistos->fh_dijetInvMDeltaPhiCutTrunc[lCBin][udijet]->Fill(mjj,hisWeight);
                    fhistos->fh_dijetInvMDeltaPhiCutTrunc2[lCBin][udijet]->Fill(mjj,hisWeight);
                    fhistos->fh_dijetPtPairDeltaPhiCut[lCBin][udijet]->Fill(ptpair,hisWeight);
                    dPhi = GetDeltaPhi(dijets.at(udijet).at(1).at(0), dijets.at(udijet).at(1).at(1));
                    fhistos->fh_dijetDeltaPhiWithCut[lCBin][udijet]->Fill(dPhi,hisWeight);
                    fhistos->fh_dijetCosDeltaPhiWithCut[lCBin][udijet]->Fill(TMath::Cos(dPhi),hisWeight);
                    fhistos->fh_dijetDeltaEtaWithCut[lCBin][udijet]->Fill(dijets.at(udijet).at(1).at(0).eta()-dijets.at(udijet).at(1).at(1).eta(),hisWeight);
                    fhistos->fh_dijetCoshDeltaEtaWithCut[lCBin][udijet]->Fill(TMath::CosH(dijets.at(udijet).at(1).at(0).eta()-dijets.at(udijet).at(1).at(1).eta()),hisWeight);
                    fhistos->fh_dijetSqrt2pt12WithCut[lCBin][udijet]->Fill(TMath::Sqrt(2*dijets.at(udijet).at(1).at(0).pt()*dijets.at(udijet).at(1).at(1).pt()),hisWeight);
                    fhistos->fh_dijetSqrtGeometryWithCut[lCBin][udijet]->Fill(TMath::Sqrt(TMath::CosH(dijets.at(udijet).at(1).at(0).eta()-dijets.at(udijet).at(1).at(1).eta())-TMath::Cos(dPhi)),hisWeight);
                }
            }
        }
    }


    if(bHasDeltaPhiDijet.at(iBGSubtr)) {
        CalculateDeltaM(iBGSubtr, 
                        uRawJetIndex.at(iBGSubtr).at(uLeading.at(iBGSubtr)),
                        uRawJetIndex.at(iBGSubtr).at(uSubleadingDeltaPhi.at(iBGSubtr)),
                        lCBin, fhistos, hisWeight);
    }
    if(bHasDeltaPhiDijet.at(iBGSubtrConstCut)) {
        CalculateDeltaM(iBGSubtrConstCut, 
                        uRawJetIndex.at(iBGSubtrConstCut).at(uLeading.at(iBGSubtrConstCut)),
                        uRawJetIndex.at(iBGSubtrConstCut).at(uSubleadingDeltaPhi.at(iBGSubtrConstCut)),
                        lCBin, fhistos, hisWeight);
    }
    if(bHasDeltaPhiDijet.at(iBGSubtrCutsRaw)) {
        CalculateDeltaM(iBGSubtrCutsRaw, 
                        uRawJetIndex.at(iBGSubtrCutsRaw).at(uLeading.at(iBGSubtrCutsRaw)),
                        uRawJetIndex.at(iBGSubtrCutsRaw).at(uSubleadingDeltaPhi.at(iBGSubtrCutsRaw)),
                        lCBin, fhistos, hisWeight);
    }
    if(bHasDeltaPhiDijet.at(iBGSubtrConstCutCutsRaw)) {
        CalculateDeltaM(iBGSubtrConstCutCutsRaw, 
                        uRawJetIndex.at(iBGSubtrConstCutCutsRaw).at(uLeading.at(iBGSubtrConstCutCutsRaw)),
                        uRawJetIndex.at(iBGSubtrConstCutCutsRaw).at(uSubleadingDeltaPhi.at(iBGSubtrConstCutCutsRaw)),
                        lCBin, fhistos, hisWeight);
    }
    if(bHasDeltaPhiDijet.at(iBGSubtrCommonEta)) {
        CalculateDeltaM(iBGSubtrCommonEta, 
                        uRawJetIndex.at(iBGSubtrCommonEta).at(uLeading.at(iBGSubtrCommonEta)),
                        uRawJetIndex.at(iBGSubtrCommonEta).at(uSubleadingDeltaPhi.at(iBGSubtrCommonEta)),
                        lCBin, fhistos, hisWeight);
    }
    if(bHasDeltaPhiDijet.at(iBGSubtrConstCutCommonEta)) {
        CalculateDeltaM(iBGSubtrConstCutCommonEta, 
                        uRawJetIndex.at(iBGSubtrConstCutCommonEta).at(uLeading.at(iBGSubtrConstCutCommonEta)),
                        uRawJetIndex.at(iBGSubtrConstCutCommonEta).at(uSubleadingDeltaPhi.at(iBGSubtrConstCutCommonEta)),
                        lCBin, fhistos, hisWeight);
    }
    if(bHasDeltaPhiDijet.at(iBGSubtrIndEta)) {
        CalculateDeltaM(iBGSubtrIndEta, 
                        uRawJetIndex.at(iBGSubtrIndEta).at(uLeading.at(iBGSubtrIndEta)),
                        uRawJetIndex.at(iBGSubtrIndEta).at(uSubleadingDeltaPhi.at(iBGSubtrIndEta)),
                        lCBin, fhistos, hisWeight);
    }
    if(bHasDeltaPhiDijet.at(iBGSubtrConstCutIndEta)) {
        CalculateDeltaM(iBGSubtrConstCutIndEta, 
                        uRawJetIndex.at(iBGSubtrConstCutIndEta).at(uLeading.at(iBGSubtrConstCutIndEta)),
                        uRawJetIndex.at(iBGSubtrConstCutIndEta).at(uSubleadingDeltaPhi.at(iBGSubtrConstCutIndEta)),
                        lCBin, fhistos, hisWeight);
    }
    return;
}

void AliJCDijetAna::CalculateDeltaM(int iJetSet, unsigned uLead, unsigned uSublead, int lcentBin, AliJCDijetHistos *fhistos, double hisWeight) {
    fastjet::PseudoJet tempLeading = rawJets.at(uLead);
    fastjet::PseudoJet tempSubLead = rawJets.at(uSublead);

    double firstJetPhi = tempLeading.phi() - TMath::Pi(); //-pi to pi
    double secndJetPhi = tempSubLead.phi() - TMath::Pi(); //-pi to pi
    double firstJetEta = tempLeading.eta();
    double secndJetEta = tempSubLead.eta();
    double firstConePhi     = firstJetPhi-(TMath::Pi()/2.0) < -TMath::Pi() ? firstJetPhi+(3.0*TMath::Pi()/2.0) : firstJetPhi-(TMath::Pi()/2.0);
    double secondConePhiAlt = secndJetPhi-(TMath::Pi()/2.0) < -TMath::Pi() ? secndJetPhi+(3.0*TMath::Pi()/2.0) : secndJetPhi-(TMath::Pi()/2.0);
    double firstConeEta = firstJetEta;
    double secondConeEta = secndJetEta;
    double distanceBtwJet2Cone1 = DeltaR(firstConeEta, secndJetEta, firstConePhi, secndJetPhi);
    double distanceBtwJet1Cone2Alt = DeltaR(secondConeEta, firstJetEta, secondConePhiAlt, firstJetPhi);
    bool bConeNearJet = distanceBtwJet2Cone1 < 2*fJetCone;
    bool bConeNearJetAlt = (distanceBtwJet2Cone1 < 2*fJetCone) || (distanceBtwJet1Cone2Alt < 2*fJetCone);
    bool bConesOverlapAlt = DeltaR(firstConeEta, secondConeEta, firstConePhi, secondConePhiAlt) < 0.8;
    fastjet::PseudoJet holderJet;
    fastjet::PseudoJet firstConeP;
    fastjet::PseudoJet secondConePAlt;
    pt = 0;
    double pt2Alt = 0;
    for (utrack = 0; utrack < chparticles.size(); utrack++) {
        phi = chparticles.at(utrack).phi() - TMath::Pi();
        eta = chparticles.at(utrack).eta();
        if(DeltaR(firstConeEta, eta, firstConePhi, phi) < fJetCone) {
            holderJet = fastjet::PseudoJet(chparticles.at(utrack).px(), chparticles.at(utrack).py(), chparticles.at(utrack).pz(), chparticles.at(utrack).E());
            firstConeP += holderJet;
            pt += holderJet.pt();

        }
        if(DeltaR(secondConeEta, eta, secondConePhiAlt, phi) < fJetCone) {
            holderJet = fastjet::PseudoJet(chparticles.at(utrack).px(), chparticles.at(utrack).py(), chparticles.at(utrack).pz(), chparticles.at(utrack).E());
            secondConePAlt += holderJet;
            pt2Alt += holderJet.pt();
        }
    }

    fastjet::PseudoJet doubleConeAlt = firstConeP + secondConePAlt;

    double localRho1 = pt/(TMath::Pi()*0.4*0.4);
    double localRho2Alt = pt2Alt/(TMath::Pi()*0.4*0.4);

    fastjet::PseudoJet bgsubtrJeFifth1 = fastjet::PseudoJet(tempLeading.px() - localRho1 * tempLeading.area_4vector().px(),
            tempLeading.py() - localRho1 * tempLeading.area_4vector().py(),
            tempLeading.pz() - localRho1 * tempLeading.area_4vector().pz(),
            tempLeading.E()  - localRho1 * tempLeading.area_4vector().E());
    fastjet::PseudoJet bgsubtrJeFifth2Alt = fastjet::PseudoJet(tempSubLead.px() - localRho2Alt * tempSubLead.area_4vector().px(),
            tempSubLead.py() - localRho2Alt * tempSubLead.area_4vector().py(),
            tempSubLead.pz() - localRho2Alt * tempSubLead.area_4vector().pz(),
            tempSubLead.E()  - localRho2Alt * tempSubLead.area_4vector().E());
    fastjet::PseudoJet doubleDeltaCone_fifthAlt = bgsubtrJeFifth1 + bgsubtrJeFifth2Alt;

    //If the local background is so large that m()<0.0, it is not meaningful any more to look
    //at the fluctuations, so these events will be removed.
    if(doubleDeltaCone_fifthAlt.m()<0.0) {
        fhistos->fh_events[lcentBin]->Fill("M-prime had negative mass",1.0);
        return;
    }
    fhistos->fh_jet2Cone1Dist->Fill(distanceBtwJet2Cone1,hisWeight);
    fhistos->fh_jet1Cone2AltDist->Fill(distanceBtwJet1Cone2Alt,hisWeight);
    fhistos->fh_doubleConeMAlt->Fill(doubleConeAlt.m(),hisWeight);
    fhistos->fh_localRho1->Fill(localRho1,hisWeight);
    fhistos->fh_localRho2Alt->Fill(localRho2Alt,hisWeight);
    fhistos->fh_deltaRho1->Fill(localRho1-rho,hisWeight); //localrho is usually bigger
    fhistos->fh_deltaRho2Alt->Fill(localRho2Alt-rho,hisWeight);
    fhistos->fh_deltaLocalRhoAlt->Fill(localRho1-localRho2Alt,hisWeight);
    if(bgsubtrJeFifth1.pt()<fLeadingJetCut)    fhistos->fh_events[lcentBin]->Fill("localRho dropped leading jet",1.0);
    if(bgsubtrJeFifth2Alt.pt()<fSubleadingJetCut) fhistos->fh_events[lcentBin]->Fill("localRho dropped subleading jet (Alt)",1.0);

    dijet = dijets.at(iJetSet).at(1).at(0) + dijets.at(iJetSet).at(1).at(1);
    mjj=dijet.m();

    double fDeltaMAlt=mjj-doubleDeltaCone_fifthAlt.m();
    fhistos->fh_dijetdeltaM5Alt[iJetSet]->Fill(fDeltaMAlt,hisWeight);
    fhistos->fh_dijetdeltaMScaled[iJetSet]->Fill(fDeltaMAlt/mjj,hisWeight);
    int iDijetBin2 = fhistos->GetDijetMClass(mjj);
    fhistos->fh_dijetdeltaM5Binned[iJetSet][iDijetBin2]->Fill(fDeltaMAlt,hisWeight);
    fhistos->fh_dijetdeltaMScaledBinned[iJetSet][iDijetBin2]->Fill(fDeltaMAlt/mjj,hisWeight);
    if(bConeNearJetAlt) fhistos->fh_dijetdeltaM5NearConeAlt[iJetSet]->Fill(fDeltaMAlt,hisWeight);
    if(bConesOverlapAlt) fhistos->fh_events[lcentBin]->Fill("cones overlap alt", 1.0);
    fhistos->fh_dijetMLocalRhoAlt[iJetSet]->Fill(doubleDeltaCone_fifthAlt.m(),hisWeight);
    fhistos->fh_deltaMResponse[iJetSet]->Fill(mjj+fDeltaMAlt, mjj,hisWeight);

    return;
}

// Response matrices are calculated in this function.
void AliJCDijetAna::CalculateResponse(AliJCDijetAna *anaDetMC, AliJCDijetHistos *fhistos, int iJetSetPart, int iJetSetDet, double hisWeight) {
    if(fDebug>8) cout << "===== BEGIN RESPONSE CALC with sets " << iJetSetPart << " and " << iJetSetDet << " =====" << endl;
    
    vector<vector<fastjet::PseudoJet>> jetsDetMC = anaDetMC->GetJets();

    unsigned Njets = jets.at(iJetSetPart).size();
    unsigned NjetsDetMC = jetsDetMC.at(iJetSetDet).size();
    //cout << "Response true jets: " << Njets << endl;
    //cout << "Response det jets size:     " << NjetsDetMC << endl;
    double minRGlobal=0; //Not restricted by matching
    double deltaRMatch=0; //Only matched jets
    double ptTrue, ptDetMC, ptTrueScalar, ptDetMCScalar;
    double deltaRLL, deltaRLS, deltaRSL, deltaRSS;
    unsigned chosenMatchIndex;
    bool bfound;
    std::vector<bool> bTrueJetMatch(Njets, false);
    std::vector<bool> bDetJetMatch(NjetsDetMC, false);
    bool bLeadingMatch    = false;
    bool bSubleadingMatch = false;
    bool bSubleadingMatchDeltaPhi = false;

    //Single jet pt matching by closest in eta-phi plane
    for (ujet = 0; ujet < Njets; ujet++) { //True MC jets
        chosenMatchIndex=-1;
        bfound=false;
        deltaR=0;
        deltaRMatch=999.0;
        minRGlobal=999.0;
        ptTrue = jets.at(iJetSetPart).at(ujet).pt();
        fhistos->fh_deltaPtResponse[iJetSetPart]->Fill(ptTrue+fDeltaPt,ptTrue,hisWeight);
        fhistos->fh_deltaPtResponse_ALICE[iJetSetPart]->Fill(ptTrue+fDeltaPt,ptTrue,hisWeight);
        for (ujetDetMC = 0; ujetDetMC < NjetsDetMC; ujetDetMC++) { //Det MC jets
            if(bDetJetMatch.at(ujetDetMC)) continue; //We want to match jets only once.
            deltaR = DeltaR(jets.at(iJetSetPart).at(ujet), jetsDetMC.at(iJetSetDet).at(ujetDetMC));
            if(deltaR<minRGlobal) minRGlobal=deltaR;
            if(deltaR < matchingR && deltaR<deltaRMatch) {
                chosenMatchIndex = ujetDetMC;
                deltaRMatch = deltaR;
                bfound = true;
            }
        }
        fhistos->fh_jetResponseDeltaRClosest[iJetSetPart][iJetSetDet]->Fill(minRGlobal,hisWeight);
        if(bfound) {
            ptDetMC = jetsDetMC.at(iJetSetDet).at(chosenMatchIndex).pt();
            fhistos->fh_jetResponse[iJetSetPart][iJetSetDet]->Fill(ptDetMC, ptTrue,hisWeight);
            fhistos->fh_jetResponse_ALICE[iJetSetPart][iJetSetDet]->Fill(ptDetMC, ptTrue,hisWeight);
            fhistos->fh_jetResponseDeltaR[iJetSetPart][iJetSetDet]->Fill(deltaRMatch,hisWeight);
            fhistos->fh_jetResponseDeltaPt[iJetSetPart][iJetSetDet]->Fill((ptTrue-ptDetMC)/ptTrue,hisWeight);
            fhistos->fh_responseInfo[iJetSetPart][iJetSetDet]->Fill("True jet has pair",1.0);
            //For comparisons, we want to save also the scalar BG removal response:
            if(iJetSetDet==iAcc && iJetSetPart==iAcc) {
                area   = jets.at(iAcc).at(ujet).area();
                ptTrueScalar = jets.at(iAcc).at(ujet).pt() - rho*area;
                area   = jetsDetMC.at(iAcc).at(chosenMatchIndex).area();
                ptDetMCScalar = jetsDetMC.at(iAcc).at(chosenMatchIndex).pt() - anaDetMC->GetRho()*area;
                fhistos->fh_jetResponse_ALICEScalar->Fill(ptDetMCScalar, ptTrueScalar,hisWeight);
            }
            bTrueJetMatch.at(ujet) = true;
            bDetJetMatch.at(chosenMatchIndex) = true;
        } else {
            fhistos->fh_responseInfo[iJetSetPart][iJetSetDet]->Fill("True jet has no pair",1.0);
        }
    }
    for (ujetDetMC = 0; ujetDetMC < NjetsDetMC; ujetDetMC++) { //Det MC jets
        if(!bDetJetMatch.at(ujetDetMC)) {
            fhistos->fh_responseInfo[iJetSetPart][iJetSetDet]->Fill("Det jet has no pair",1.0);
        }
    }

    if(fDebug>8) {
        cout << "True jets: " << Njets << ", det jets: " << NjetsDetMC << endl;
        cout << "True jets:" << endl;
        for (ujet = 0; ujet < Njets; ujet++) { //True MC jets
            cout << "jet(E, pt, phi, eta) = " <<
                "(" << jets.at(iJetSetPart).at(ujet).E()   << ", "
                << jets.at(iJetSetPart).at(ujet).pt()  << ", "
                << jets.at(iJetSetPart).at(ujet).phi() << ", "
                << jets.at(iJetSetPart).at(ujet).eta() << ") "
                << (bTrueJetMatch.at(ujet) ? "match" : "no match") << endl;
        }
        cout << "Det jets:" <<  endl;
        for (ujetDetMC = 0; ujetDetMC < NjetsDetMC; ujetDetMC++) { //Det MC jets
            cout << "jet(E, pt, phi, eta) = " <<
                "(" << jetsDetMC.at(iJetSetDet).at(ujetDetMC).E()   << ", "
                << jetsDetMC.at(iJetSetDet).at(ujetDetMC).pt()  << ", "
                << jetsDetMC.at(iJetSetDet).at(ujetDetMC).phi() << ", "
                << jetsDetMC.at(iJetSetDet).at(ujetDetMC).eta() << ") "
                << (bDetJetMatch.at(ujetDetMC) ? "match" : "no match") << endl;
        }
    }

    vector<vector<vector<fastjet::PseudoJet>>> dijetsDetMC = anaDetMC->GetDijets();
    fastjet::PseudoJet dijetDetMC;

    //Dijet response without deltaphi cut.
    if(bHasDijet.at(iJetSetPart)) {
        if(fDebug>8) {
            cout << "Normal dijet:" << endl;
            cout << "True dijet" << endl;
            cout << "leading jet(E, pt, phi, eta) =    " <<
                "(" << dijets.at(iJetSetPart).at(0).at(0).E()   << ", "
                << dijets.at(iJetSetPart).at(0).at(0).pt()  << ", "
                << dijets.at(iJetSetPart).at(0).at(0).phi() << ", "
                << dijets.at(iJetSetPart).at(0).at(0).eta() << ")" << endl;
            cout << "subleading jet(E, pt, phi, eta) = " <<
                "(" << dijets.at(iJetSetPart).at(0).at(1).E()   << ", "
                << dijets.at(iJetSetPart).at(0).at(1).pt()  << ", "
                << dijets.at(iJetSetPart).at(0).at(1).phi() << ", "
                << dijets.at(iJetSetPart).at(0).at(1).eta() << ")" << endl;
        }
        dijet = dijets.at(iJetSetPart).at(0).at(0) + dijets.at(iJetSetPart).at(0).at(1);
        if(anaDetMC->HasDijet(iJetSetDet)) {

            // Check that leading and subleading jets match.
            // It is also ok if leading and subleading jets change places in det level.
            deltaRLL = DeltaR(dijets.at(iJetSetPart).at(0).at(0),dijetsDetMC.at(iJetSetDet).at(0).at(0)); //leading, leading
            deltaRLS = DeltaR(dijets.at(iJetSetPart).at(0).at(0),dijetsDetMC.at(iJetSetDet).at(0).at(1)); //leading, subleading
            deltaRSS = DeltaR(dijets.at(iJetSetPart).at(0).at(1),dijetsDetMC.at(iJetSetDet).at(0).at(1)); //subleading, subleading
            deltaRSL = DeltaR(dijets.at(iJetSetPart).at(0).at(1),dijetsDetMC.at(iJetSetDet).at(0).at(0)); //subleading, leading

            if(deltaRLL < matchingR || deltaRLS < matchingR) bLeadingMatch    = true;
            if(deltaRSS < matchingR || deltaRSL < matchingR) bSubleadingMatch = true;

            if(fDebug>8) {
                cout << "Det dijet" << endl;
                cout << "leading jet(E, pt, phi, eta) =    " <<
                    "(" << dijetsDetMC.at(iJetSetDet).at(0).at(0).E()   << ", "
                    << dijetsDetMC.at(iJetSetDet).at(0).at(0).pt()  << ", "
                    << dijetsDetMC.at(iJetSetDet).at(0).at(0).phi() << ", "
                    << dijetsDetMC.at(iJetSetDet).at(0).at(0).eta() << ")" << endl;
                cout << "subleading jet(E, pt, phi, eta) = " <<
                    "(" << dijetsDetMC.at(iJetSetDet).at(0).at(1).E()   << ", "
                    << dijetsDetMC.at(iJetSetDet).at(0).at(1).pt()  << ", "
                    << dijetsDetMC.at(iJetSetDet).at(0).at(1).phi() << ", "
                    << dijetsDetMC.at(iJetSetDet).at(0).at(1).eta() << ")" << endl;
            }
            dijetDetMC = dijetsDetMC.at(iJetSetDet).at(0).at(0) + dijetsDetMC.at(iJetSetDet).at(0).at(1);
            fhistos->fh_dijetResponseLinNoMatching[iJetSetPart][iJetSetDet]->Fill(dijetDetMC.m(), dijet.m(),hisWeight);
            if(bLeadingMatch && bSubleadingMatch) {
                fhistos->fh_dijetResponse[iJetSetPart][iJetSetDet]->Fill(dijetDetMC.m(), dijet.m(),hisWeight);
                fhistos->fh_dijetResponseLin[iJetSetPart][iJetSetDet]->Fill(dijetDetMC.m(), dijet.m(),hisWeight);
                fhistos->fh_dijetResponseTrunc[iJetSetPart][iJetSetDet]->Fill(dijetDetMC.m(), dijet.m(),hisWeight);
                fhistos->fh_dijetResponseTrunc2[iJetSetPart][iJetSetDet]->Fill(dijetDetMC.m(), dijet.m(),hisWeight);
                fhistos->fh_responseInfo[iJetSetPart][iJetSetDet]->Fill("Dijet match",1.0);
            } else {
                fhistos->fh_responseInfo[iJetSetPart][iJetSetDet]->Fill("Dijet not match",1.0);
            }
        } else {
            fhistos->fh_responseInfo[iJetSetPart][iJetSetDet]->Fill("Dijet det not found",1.0);
        }
    } else {
        if(anaDetMC->HasDijet(iJetSetDet)) {
            //dijetDetMC = dijetsDetMC.at(iJetSetDet).at(0).at(0) + dijetsDetMC.at(iJetSetDet).at(0).at(1);
            fhistos->fh_responseInfo[iJetSetPart][iJetSetDet]->Fill("Dijet true not found",1.0);
        }
    }
    if(fDebug>8) cout << ((anaDetMC->HasDijet(iJetSetDet) && bLeadingMatch && bSubleadingMatch) ? "Dijet match" : "Dijet no match") << endl;

    // DeltaPhi cut dijet response.
    if(bHasDeltaPhiDijet.at(iJetSetPart)) {
        if(fDebug>8) {
            cout << "Delta Phi dijet:" << endl;
            cout << "True dijet" << endl;
            cout << "leading jet(E, pt, phi, eta) =    " <<
                "(" << dijets.at(iJetSetPart).at(1).at(0).E()   << ", "
                << dijets.at(iJetSetPart).at(1).at(0).pt()  << ", "
                << dijets.at(iJetSetPart).at(1).at(0).phi() << ", "
                << dijets.at(iJetSetPart).at(1).at(0).eta() << ")" << endl;
            cout << "subleading jet(E, pt, phi, eta) = " <<
                "(" << dijets.at(iJetSetPart).at(1).at(1).E()   << ", "
                << dijets.at(iJetSetPart).at(1).at(1).pt()  << ", "
                << dijets.at(iJetSetPart).at(1).at(1).phi() << ", "
                << dijets.at(iJetSetPart).at(1).at(1).eta() << ")" << endl;
        }
        dijet = dijets.at(iJetSetPart).at(1).at(0) + dijets.at(iJetSetPart).at(1).at(1);
        if(anaDetMC->HasDeltaPhiDijet(iJetSetDet)) {
        
            // Leading jets are the same as without delta phi cut.
            deltaRSS = DeltaR(dijets.at(iJetSetPart).at(1).at(1),dijetsDetMC.at(iJetSetDet).at(1).at(1)); //subleading, subleading
            deltaRSL = DeltaR(dijets.at(iJetSetPart).at(1).at(1),dijetsDetMC.at(iJetSetDet).at(1).at(0)); //subleading, leading
            if(deltaRSS < matchingR || deltaRSL < matchingR) bSubleadingMatchDeltaPhi = true;

            if(fDebug>8) {
                cout << "Det dijet" << endl;
                cout << "leading jet(E, pt, phi, eta) =    " <<
                    "(" << dijetsDetMC.at(iJetSetDet).at(1).at(0).E()   << ", "
                    << dijetsDetMC.at(iJetSetDet).at(1).at(0).pt()  << ", "
                    << dijetsDetMC.at(iJetSetDet).at(1).at(0).phi() << ", "
                    << dijetsDetMC.at(iJetSetDet).at(1).at(0).eta() << ")" << endl;
                cout << "subleading jet(E, pt, phi, eta) = " <<
                    "(" << dijetsDetMC.at(iJetSetDet).at(1).at(1).E()   << ", "
                    << dijetsDetMC.at(iJetSetDet).at(1).at(1).pt()  << ", "
                    << dijetsDetMC.at(iJetSetDet).at(1).at(1).phi() << ", "
                    << dijetsDetMC.at(iJetSetDet).at(1).at(1).eta() << ")" << endl;
            }
            dijetDetMC = dijetsDetMC.at(iJetSetDet).at(1).at(0) + dijetsDetMC.at(iJetSetDet).at(1).at(1);
            // Check subleading jet match.

            fhistos->fh_dijetResponseDeltaPhiCutLinNoMatching[iJetSetPart][iJetSetDet]->Fill(dijetDetMC.m(), dijet.m(),hisWeight);
            if(bLeadingMatch && bSubleadingMatchDeltaPhi) {
                fhistos->fh_dijetResponseDeltaPhiCut[iJetSetPart][iJetSetDet]->Fill(dijetDetMC.m(), dijet.m(),hisWeight);
                fhistos->fh_dijetResponseDeltaPhiCutLin[iJetSetPart][iJetSetDet]->Fill(dijetDetMC.m(), dijet.m(),hisWeight);
                fhistos->fh_dijetResponseDeltaPhiCutTrunc[iJetSetPart][iJetSetDet]->Fill(dijetDetMC.m(), dijet.m(),hisWeight);
                fhistos->fh_dijetResponseDeltaPhiCutTrunc2[iJetSetPart][iJetSetDet]->Fill(dijetDetMC.m(), dijet.m(),hisWeight);
                fhistos->fh_responseInfo[iJetSetPart][iJetSetDet]->Fill("Dijet DPhi match",1.0);
            } else {
                fhistos->fh_responseInfo[iJetSetPart][iJetSetDet]->Fill("Dijet DPhi not match",1.0);
            }
        } else {
            fhistos->fh_responseInfo[iJetSetPart][iJetSetDet]->Fill("Dijet DPhi det not found",1.0);
        }
    } else {
        if(anaDetMC->HasDeltaPhiDijet(iJetSetDet)) {
            //dijetDetMC = dijetsDetMC.at(iJetSetDet).at(1).at(0) + dijetsDetMC.at(iJetSetDet).at(1).at(1);
            fhistos->fh_responseInfo[iJetSetPart][iJetSetDet]->Fill("Dijet DPhi true not found",1.0);
        }
    }
    if(fDebug>8) cout << ((anaDetMC->HasDeltaPhiDijet(iJetSetDet) && bLeadingMatch && bSubleadingMatchDeltaPhi) ? "Delta Phi Dijet match" : "Delta Phi Dijet no match") << endl;
    return;
}

void AliJCDijetAna::ResetObjects() {
    chparticles.clear();
    ktchparticles.clear();
    for(int i=0; i<jetClassesSize; i++) {
        jets.at(i).clear();
        bHasDijet.at(i) = false;
        bHasDeltaPhiDijet.at(i) = false;
    }
    rawJets.clear();
    tempJets.clear();
    rawKtJets.clear();
    rhoEstJets.clear();

    phi=0, eta=0, pt=0, pt2=0, rho=0, rhom=0, area=0, mjj=0, ptpair=0, dPhi=0;
    leadingTrackOverThreshold = false;
    bHasDeltaPhiSubLeadJet = false;
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
void AliJCDijetAna::InitHistos(AliJCDijetHistos *histos, bool bIsMC, int nCentBins, int iJetClassTrue, int iJetClassDet) {
    cout << "Initing histograms for AliJCDijetAna" << endl;
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
    histos->fh_eventSel->Fill("alieventcut ok",0.0);
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
        histos->fh_events[iBin]->Fill("bg. subtr. jets (common eta)",0.0);
        histos->fh_events[iBin]->Fill("bg. subtr. const. cut jets (common eta)",0.0);
        histos->fh_events[iBin]->Fill("bg. subtr. jets (bg eta)",0.0);
        histos->fh_events[iBin]->Fill("bg. subtr. const. cut jets (bg eta)",0.0);
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
        histos->fh_events[iBin]->Fill("bg. subtr. cuts raw dijets",0.0);
        histos->fh_events[iBin]->Fill("bg. subtr. cuts raw dijets leading cut",0.0);
        histos->fh_events[iBin]->Fill("bg. subtr. cuts raw acc. dijets",0.0);
        histos->fh_events[iBin]->Fill("bg. subtr. cuts raw deltaphi cut dijets",0.0);
        histos->fh_events[iBin]->Fill("bg. subtr. const. cut cuts raw dijets",0.0);
        histos->fh_events[iBin]->Fill("bg. subtr. const. cut cuts raw dijets leading cut",0.0);
        histos->fh_events[iBin]->Fill("bg. subtr. const. cut cuts raw acc. dijets",0.0);
        histos->fh_events[iBin]->Fill("bg. subtr. const. cut cuts raw deltaphi cut dijets",0.0);
        histos->fh_events[iBin]->Fill("kt dijets",0.0);
        histos->fh_events[iBin]->Fill("kt dijets leading cut",0.0);
        histos->fh_events[iBin]->Fill("kt acc. dijets",0.0);
        histos->fh_events[iBin]->Fill("kt deltaphi cut dijets",0.0);
        histos->fh_events[iBin]->Fill("pt_hard bin cuts",0.0);
        histos->fh_events[iBin]->Fill("area cut ok",0.0);
        histos->fh_events[iBin]->Fill("area cut fails",0.0);
        histos->fh_events[iBin]->Fill("100GeV const cut ok",0.0);
        histos->fh_events[iBin]->Fill("100GeV const cut fails",0.0);
        histos->fh_events[iBin]->Fill("localRho dropped leading jet",0.0);
        histos->fh_events[iBin]->Fill("localRho dropped subleading jet (Alt)",0.0);
        histos->fh_events[iBin]->Fill("M-prime had negative mass",0.0);
        histos->fh_events[iBin]->Fill("localRho increased leading jet",0.0);
        histos->fh_events[iBin]->Fill("localRho increased subleading jet",0.0);
        histos->fh_events[iBin]->Fill("cones overlap",0.0);
        histos->fh_events[iBin]->Fill("cones overlap alt",0.0);
        histos->fh_events[iBin]->Fill("trackingIneff rejection",0.0);
        histos->fh_events[iBin]->Fill("trackingIneff accept",0.0);
    }
    if(bIsMC) {
        for(int iBin=0; iBin<jetClassesSize-1; iBin++) {
            histos->fh_responseInfo[iBin][iBin]->Fill("True jet has pair",0.0);
            histos->fh_responseInfo[iBin][iBin]->Fill("True jet has no pair",0.0);
            histos->fh_responseInfo[iBin][iBin]->Fill("Det jet has no pair",0.0);
            histos->fh_responseInfo[iBin][iBin]->Fill("Dijet match",0.0);
            histos->fh_responseInfo[iBin][iBin]->Fill("Dijet not match",0.0);
            histos->fh_responseInfo[iBin][iBin]->Fill("Dijet det not found",0.0);
            histos->fh_responseInfo[iBin][iBin]->Fill("Dijet true not found",0.0);
            histos->fh_responseInfo[iBin][iBin]->Fill("Dijet DPhi match",0.0);
            histos->fh_responseInfo[iBin][iBin]->Fill("Dijet DPhi not match",0.0);
            histos->fh_responseInfo[iBin][iBin]->Fill("Dijet DPhi det not found",0.0);
            histos->fh_responseInfo[iBin][iBin]->Fill("Dijet DPhi true not found",0.0);
        }
        if(iJetClassDet!=iJetClassTrue) {
            histos->fh_responseInfo[iJetClassTrue][iJetClassDet]->Fill("True jet has pair",0.0);
            histos->fh_responseInfo[iJetClassTrue][iJetClassDet]->Fill("True jet has no pair",0.0);
            histos->fh_responseInfo[iJetClassTrue][iJetClassDet]->Fill("Det jet has no pair",0.0);
            histos->fh_responseInfo[iJetClassTrue][iJetClassDet]->Fill("Dijet match",0.0);
            histos->fh_responseInfo[iJetClassTrue][iJetClassDet]->Fill("Dijet not match",0.0);
            histos->fh_responseInfo[iJetClassTrue][iJetClassDet]->Fill("Dijet det not found",0.0);
            histos->fh_responseInfo[iJetClassTrue][iJetClassDet]->Fill("Dijet true not found",0.0);
            histos->fh_responseInfo[iJetClassTrue][iJetClassDet]->Fill("Dijet DPhi match",0.0);
            histos->fh_responseInfo[iJetClassTrue][iJetClassDet]->Fill("Dijet DPhi not match",0.0);
            histos->fh_responseInfo[iJetClassTrue][iJetClassDet]->Fill("Dijet DPhi det not found",0.0);
            histos->fh_responseInfo[iJetClassTrue][iJetClassDet]->Fill("Dijet DPhi true not found",0.0);
        }
    }
    return;
}
#endif
