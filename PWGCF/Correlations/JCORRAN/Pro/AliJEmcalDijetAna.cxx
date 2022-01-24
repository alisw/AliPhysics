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
#include "AliJEmcalDijetAna.h"

ClassImp(AliJEmcalDijetAna)

//----------------------------------------------------------------------------------------------------------------------------
AliJEmcalDijetAna::AliJEmcalDijetAna() :
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
AliJEmcalDijetAna::~AliJEmcalDijetAna(){
    // destructor
}

//----------------------------------------------------------------------------------------------------------------------------
AliJEmcalDijetAna::AliJEmcalDijetAna(const AliJEmcalDijetAna& obj) :
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
AliJEmcalDijetAna& AliJEmcalDijetAna::operator=(const AliJEmcalDijetAna& obj){
    //----------------------------------------------------------------------------------------------------------------------------
    // equal sign operator
    return *this;
}

#if !defined(__CINT__) && !defined(__MAKECINT__)
void AliJEmcalDijetAna::SetSettings(int    lDebug,
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
                 ::Error("AliJEmcalDijetAna","Unknown recombination scheme!");
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
                 ::Error("AliJEmcalDijetAna","Unknown recombination scheme!");
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
    
    //OpenFile(1);
    //fOutput = gDirectory;
    //fOutput->cd();

    AllocateTrackHistograms();

    TIter next(fhistos.GetListOfHistograms());
    TObject* obj = 0;
    while ((obj = next())) {
        fOutput->Add(obj);
    }

    return;
}

//______________________________________________________________________________
int AliJEmcalDijetAna::CalculateJets(TClonesArray *inList, int lCBin){

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

    rawJets   = fastjet::sorted_by_pt(cs->inclusive_jets(MinJetPt)); // APPLY Min pt cut for jet

    // For MC runs: If we find jets with over 4 times pt_hard bin, reject the event.
    if( fptHardBin!=0 && rawJets.size()>0 ) {
        fHistname = TString::Format("%s/h_ptHard_%d", fGroupname.Data(), lCBin);
        fhistos.FillTH1(fHistname, fptHardBin);
        //fhistos->fh_ptHard[lCBin]->Fill(fptHardBin);
        fHistname = TString::Format("%s/h_maxJetptOverPtHard_%d", fGroupname.Data(), lCBin);
        fhistos.FillTH1(fHistname, rawJets.at(0).pt()/fptHardBin);
        if( rawJets.at(0).pt() > fptHardBin*4 ) {
            fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), lCBin);
            fhistos.FillTH1(fHistname, "pt_hard bin cuts",1.0);
            return -1;
        }
    }
    fHistname = TString::Format("%s/h_randConeEtaPhi_%d", fGroupname.Data(), lCBin);
    fhistos.FillTH2(fHistname, randConeEta,randConePhi);

    rawKtJets = fastjet::sorted_by_pt(cs_bge->inclusive_jets(0.0)); // APPLY Min pt cut for jet

    fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), lCBin);
    fhistos.FillTH1(fHistname, "particles",noTracks);
    for (utrack = 0; utrack < chparticles.size(); utrack++) {
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), lCBin);
        fhistos.FillTH1(fHistname, "acc. particles",1.0);
        pt = chparticles.at(utrack).pt();
        eta = chparticles.at(utrack).eta();
        phi = chparticles.at(utrack).phi();
        fHistname = TString::Format("%s/h_eta_%d", fGroupname.Data(), lCBin);
        fhistos.FillTH1(fHistname, eta);
        fHistname = TString::Format("%s/h_phi_%d", fGroupname.Data(), lCBin);
        fhistos.FillTH1(fHistname, phi);
        fHistname = TString::Format("%s/h_etaPhi_%d", fGroupname.Data(), lCBin);
        fhistos.FillTH2(fHistname, eta,phi);
        fHistname = TString::Format("%s/h_pt_%d", fGroupname.Data(), lCBin);
        fhistos.FillTH1(fHistname, pt);
        if(eta>0.0) {
            fHistname = TString::Format("%s/h_ptPosEta_%d", fGroupname.Data(), lCBin);
            fhistos.FillTH1(fHistname, pt);
        } else {
            fHistname = TString::Format("%s/h_ptNegEta_%d", fGroupname.Data(), lCBin);
            fhistos.FillTH1(fHistname, pt);
        }
    }
    fHistname = TString::Format("%s/h_nch_%d", fGroupname.Data(), lCBin);
    fhistos.FillTH1(fHistname, chparticles.size());
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
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), lCBin);
        fhistos.FillTH1(fHistname, "no rho calc. events",1.0);
        rho  = 0.0;
        rhom = 0.0;
    } else { 
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), lCBin);
        fhistos.FillTH1(fHistname, "rho calc. events",1.0);
        bge.set_jets(rhoEstJets);
        rho  = bge.rho()<0   ? 0.0 : bge.rho();
        rhom = bge.rho_m()<0 ? 0.0 : bge.rho_m();
    }


    fHistname = TString::Format("%s/h_rho_%d", fGroupname.Data(), lCBin);
    fhistos.FillTH1(fHistname, rho);
    fHistname = TString::Format("%s/h_rhoLin_%d", fGroupname.Data(), lCBin);
    fhistos.FillTH1(fHistname, rho);
    fHistname = TString::Format("%s/h_rhom_%d", fGroupname.Data(), lCBin);
    fhistos.FillTH1(fHistname, rhom);
    fHistname = TString::Format("%s/h_rhomLin_%d", fGroupname.Data(), lCBin);
    fhistos.FillTH1(fHistname, rhom);
    if(fDebug > 9) std::cout << "Testing: Rho_M = " << rhom << ", has_rho_m() = " << bge.has_rho_m() << std::endl;

    // This is the delta-phi distribution: delta-pt = sum_{tracks in rand cone}( pt ) - rho*pi*R^2
    fDeltaPt = randConePt - rho*TMath::Pi()*fJetCone*fJetCone;
    fHistname = TString::Format("%s/h_deltaPt_%d", fGroupname.Data(), lCBin);
    fhistos.FillTH1(fHistname, fDeltaPt);

    bEvtHasAreaInfo = true;
    return 0;
}

// If user wants to add jets from outside this function can be used.
void AliJEmcalDijetAna::SetJets(vector<fastjet::PseudoJet> jetsOutside) {
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
void AliJEmcalDijetAna::FillJetsDijets(int lCBin) {
    TString sDijetTypes[jetClassesSize] = {"raw", "bg. subtr.", "bg. subtr. const. cut", "const. cut", "kt"};
    int iAccJetCounter = 0;
    bool bHasHighPtJet = false;

    // anti-kt jets:
    for (ujet = 0; ujet < rawJets.size(); ujet++) {
        eta = rawJets.at(ujet).eta();
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), lCBin);
        fhistos.FillTH1(fHistname, "jets",1.0);
        // anti-kt-jet eta cut
        if(TMath::Abs(eta) < etaMaxCutForJet) {
            iAccJetCounter++;
            fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), lCBin);
            fhistos.FillTH1(fHistname, "acc. jets",1.0);
            jets.at(iAcc).push_back(rawJets.at(ujet));
            pt = rawJets.at(ujet).pt();
            phi = rawJets.at(ujet).phi();
            if(iAccJetCounter==1) bHasHighPtJet = pt > fLeadingJetCut;
            if(bEvtHasAreaInfo) area = rawJets.at(ujet).area();
            if(bEvtHasAreaInfo) jetAreaVector = rawJets.at(ujet).area_4vector();
            fHistname = TString::Format("%s/h_jetEta_%d_jetC%d",fGroupname.Data(), lCBin, iAcc);
            fhistos.FillTH1(fHistname, eta);  
            fHistname = TString::Format("%s/h_jetPhi_%d_jetC%d",fGroupname.Data(), lCBin, iAcc);
            fhistos.FillTH1(fHistname, phi - TMath::Pi()); //Pseudojet.phi range 0-2pi
            fHistname = TString::Format("%s/h_jetEtaPhi_%d_jetC%d",fGroupname.Data(), lCBin, iAcc);
            fhistos.FillTH2(fHistname, eta,phi - TMath::Pi());
            fHistname = TString::Format("%s/h_jetPt_%d_jetC%d",fGroupname.Data(), lCBin, iAcc);
            fhistos.FillTH1(fHistname, pt);
            fHistname = TString::Format("%s/h_jetPt_ALICE_%d_jetC%d",fGroupname.Data(), lCBin, iAcc);
            fhistos.FillTH1(fHistname, pt);
            if(bEvtHasAreaInfo) {
                fHistname = TString::Format("%s/h_jetArea_%d_jetC%d",fGroupname.Data(), lCBin, iAcc);
                fhistos.FillTH1(fHistname, area);
                fHistname = TString::Format("%s/h_jetAreaRho_%d_jetC%d",fGroupname.Data(), lCBin, iAcc);
                fhistos.FillTH1(fHistname, area*rho);
            }

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
                fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), lCBin);
                fhistos.FillTH1(fHistname, "const. cut jets",1.0);
                fHistname = TString::Format("%s/h_jetEta_%d_jetC%d",fGroupname.Data(), lCBin, iConstCut);
                fhistos.FillTH1(fHistname, eta);
                fHistname = TString::Format("%s/h_jetPhi_%d_jetC%d",fGroupname.Data(), lCBin, iConstCut);
                fhistos.FillTH1(fHistname, phi - TMath::Pi());
                fHistname = TString::Format("%s/h_jetEtaPhi_%d_jetC%d",fGroupname.Data(), lCBin, iConstCut);
                fhistos.FillTH2(fHistname, eta,phi - TMath::Pi());
                fHistname = TString::Format("%s/h_jetPt_%d_jetC%d",fGroupname.Data(), lCBin, iConstCut);
                fhistos.FillTH1(fHistname, pt);
                fHistname = TString::Format("%s/h_jetPt_ALICE_%d_jetC%d",fGroupname.Data(), lCBin, iConstCut);
                fhistos.FillTH1(fHistname, pt);
                if(bEvtHasAreaInfo) {
                    fHistname = TString::Format("%s/h_jetArea_%d_jetC%d",fGroupname.Data(), lCBin, iConstCut);
                    fhistos.FillTH1(fHistname, area);
                    fHistname = TString::Format("%s/h_jetAreaRho_%d_jetC%d",fGroupname.Data(), lCBin, iConstCut);
                    fhistos.FillTH1(fHistname, area*rho);
                }

                jets.at(iConstCut).push_back(rawJets.at(ujet));
            }

            if(bEvtHasAreaInfo) {
                jet_bgSubtracted = fastjet::PseudoJet(rawJets.at(ujet).px() -        rho * jetAreaVector.px(),
                                                      rawJets.at(ujet).py() -        rho * jetAreaVector.py(),
                                                      rawJets.at(ujet).pz() - (rho+rhom) * jetAreaVector.pz(),
                                                      rawJets.at(ujet).E()  - (rho+rhom) * jetAreaVector.E());
                fHistname = TString::Format("%s/h_jetBGSubtrDeltaR_%d", fGroupname.Data(), lCBin);
                fhistos.FillTH1(fHistname, DeltaR(rawJets.at(ujet),jet_bgSubtracted));

                // Check eta acceptance also for bg subtracted jets.
                eta = jet_bgSubtracted.eta();
                if(TMath::Abs(eta) < etaMaxCutForJet) {
                    fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), lCBin);
                    fhistos.FillTH1(fHistname, "bg. subtr. jets",1.0);
                    pt2 = jet_bgSubtracted.pt();
                    phi = jet_bgSubtracted.phi();
                    if(iAccJetCounter==0 && pt>fLeadingJetCut && pt2<=fLeadingJetCut) {
                        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), lCBin);
                        fhistos.FillTH1(fHistname, "leading jet drop",1.0);
                    }
                    if(iAccJetCounter==1 && pt>fSubleadingJetCut && pt2<=fSubleadingJetCut) {
                        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), lCBin);
                        fhistos.FillTH1(fHistname, "subleading jet drop",1.0);
                    }
                    fHistname = TString::Format("%s/h_jetEta_%d_jetC%d",fGroupname.Data(), lCBin, iBGSubtr);
                    fhistos.FillTH1(fHistname, eta);
                    fHistname = TString::Format("%s/h_jetPhi_%d_jetC%d",fGroupname.Data(), lCBin, iBGSubtr);
                    fhistos.FillTH1(fHistname, phi - TMath::Pi());
                    fHistname = TString::Format("%s/h_jetEtaPhi_%d_jetC%d",fGroupname.Data(), lCBin, iBGSubtr);
                    fhistos.FillTH2(fHistname, eta,phi - TMath::Pi());
                    fHistname = TString::Format("%s/h_jetPt_%d_jetC%d",fGroupname.Data(), lCBin, iBGSubtr);
                    fhistos.FillTH1(fHistname, pt2);
                    fHistname = TString::Format("%s/h_jetPt_ALICE_%d_jetC%d",fGroupname.Data(), lCBin, iBGSubtr);
                    fhistos.FillTH1(fHistname, pt2);
                    fHistname = TString::Format("%s/h_jetPtTransBGSub_%d_jetC%d",fGroupname.Data(), lCBin, iBGSubtr);
                    fhistos.FillTH1(fHistname, pt - rho*area);
                    if(bEvtHasAreaInfo) {
                        fHistname = TString::Format("%s/h_jetArea_%d_jetC%d",fGroupname.Data(), lCBin, iBGSubtr);
                        fhistos.FillTH1(fHistname, area); // Assuming bg subtracted jet has the same area.
                        fHistname = TString::Format("%s/h_jetAreaRho_%d_jetC%d",fGroupname.Data(), lCBin, iBGSubtr);
                        fhistos.FillTH1(fHistname, area*rho);
                    }

                    jets.at(iBGSubtr).push_back(jet_bgSubtracted);

                    if(leadingTrackOverThreshold) {
                        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), lCBin);
                        fhistos.FillTH1(fHistname, "bg. subtr. const. cut jets",1.0);
                        fHistname = TString::Format("%s/h_jetEta_%d_jetC%d",fGroupname.Data(), lCBin, iBGSubtrConstCut);
                        fhistos.FillTH1(fHistname, eta);
                        fHistname = TString::Format("%s/h_jetPhi_%d_jetC%d",fGroupname.Data(), lCBin, iBGSubtrConstCut);
                        fhistos.FillTH1(fHistname, phi - TMath::Pi());
                        fHistname = TString::Format("%s/h_jetEtaPhi_%d_jetC%d",fGroupname.Data(), lCBin, iBGSubtrConstCut);
                        fhistos.FillTH2(fHistname, eta,phi - TMath::Pi());
                        fHistname = TString::Format("%s/h_jetPt_%d_jetC%d",fGroupname.Data(), lCBin, iBGSubtrConstCut);
                        fhistos.FillTH1(fHistname, pt2);
                        fHistname = TString::Format("%s/h_jetPt_ALICE_%d_jetC%d",fGroupname.Data(), lCBin, iBGSubtrConstCut);
                        fhistos.FillTH1(fHistname, pt2);
                        fHistname = TString::Format("%s/h_jetPtTransBGSub_%d_jetC%d",fGroupname.Data(), lCBin, iBGSubtrConstCut);
                        fhistos.FillTH1(fHistname, pt - rho*area);
                        if(bEvtHasAreaInfo) {
                            fHistname = TString::Format("%s/h_jetArea_%d_jetC%d",fGroupname.Data(), lCBin, iBGSubtrConstCut);
                            fhistos.FillTH1(fHistname, area);
                            fHistname = TString::Format("%s/h_jetAreaRho_%d_jetC%d",fGroupname.Data(), lCBin, iBGSubtrConstCut);
                            fhistos.FillTH1(fHistname, area*rho);
                        }

                        jets.at(iBGSubtrConstCut).push_back(jet_bgSubtracted);
                    }
                }
            }
        }
    }//end of the anti-kt-jet loop
    //============================================================================================
    //cout << "Jets: " << jets[iAcc].size() << endl;
    //============================================================================================
    fHistname = TString::Format("%s/h_jetN_%d", fGroupname.Data(), lCBin);
    fhistos.FillTH1(fHistname, iAccJetCounter);
    if(bHasHighPtJet) {
        fHistname = TString::Format("%s/h_rhoHighPt_%d", fGroupname.Data(), lCBin);
        fhistos.FillTH1(fHistname, rho);
        fHistname = TString::Format("%s/h_rhoLinHighPt_%d", fGroupname.Data(), lCBin);
        fhistos.FillTH1(fHistname, rho);
        fHistname = TString::Format("%s/h_rhomHighPt_%d", fGroupname.Data(), lCBin);
        fhistos.FillTH1(fHistname, rhom);
        fHistname = TString::Format("%s/h_rhomLinHighPt_%d", fGroupname.Data(), lCBin);
        fhistos.FillTH1(fHistname, rhom);
    }

    for (uktjet = 0; uktjet < rawKtJets.size(); uktjet++) {
        eta = rawKtJets.at(uktjet).eta();
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), lCBin);
        fhistos.FillTH1(fHistname, "kt-jets",1.0);
        // kt-jet eta cut
        if(TMath::Abs(eta) < etaMaxCutForKtJet) {
            fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), lCBin);
            fhistos.FillTH1(fHistname, "acc. kt-jets",1.0);
            jets.at(iktJets).push_back(rawKtJets.at(uktjet));
            pt = rawKtJets.at(uktjet).pt();
            phi = rawKtJets.at(uktjet).phi();
            if(bEvtHasAreaInfo) area = rawKtJets.at(uktjet).area();
            fHistname = TString::Format("%s/h_jetEta_%d_jetC%d",fGroupname.Data(), lCBin, iktJets);
            fhistos.FillTH1(fHistname, eta);  
            fHistname = TString::Format("%s/h_jetPhi_%d_jetC%d",fGroupname.Data(), lCBin, iktJets);
            fhistos.FillTH1(fHistname, phi - TMath::Pi()); //Pseudojet.phi range 0-2pi
            fHistname = TString::Format("%s/h_jetEtaPhi_%d_jetC%d",fGroupname.Data(), lCBin, iktJets);
            fhistos.FillTH2(fHistname, eta,phi - TMath::Pi());
            fHistname = TString::Format("%s/h_jetPt_%d_jetC%d",fGroupname.Data(), lCBin, iktJets);
            fhistos.FillTH1(fHistname, pt);
            fHistname = TString::Format("%s/h_jetPt_ALICE_%d_jetC%d",fGroupname.Data(), lCBin, iktJets);
            fhistos.FillTH1(fHistname, pt);
            if(bEvtHasAreaInfo) {
                fHistname = TString::Format("%s/h_jetArea_%d_jetC%d",fGroupname.Data(), lCBin, iktJets);
                fhistos.FillTH1(fHistname, area);
                fHistname = TString::Format("%s/h_jetAreaRho_%d_jetC%d",fGroupname.Data(), lCBin, iktJets);
                fhistos.FillTH1(fHistname, area*rho);
            }
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
        fHistname = TString::Format("%s/h_jetDeltaRMin_%d_jetC%d",fGroupname.Data(), lCBin, udijet);
        fhistos.FillTH1(fHistname, deltaRMin);
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
            fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), lCBin);
            fhistos.FillTH1(fHistname, Form("%s dijets",sDijetTypes[udijet].Data()),1.0);

            // Start by forming dijets. First without deltaPhiCut and then with the cut.
            
            // No deltaPhi cut for these jets.
            dijets.at(udijet).at(0).at(0) = jets.at(udijet).at(0);
            if(dijets.at(udijet).at(0).at(0).pt() < fLeadingJetCut) continue;
            fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), lCBin);
            fhistos.FillTH1(fHistname, Form("%s dijets leading cut",sDijetTypes[udijet].Data()),1.0);
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
                fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), lCBin);
                fhistos.FillTH1(fHistname, Form("%s acc. dijets",sDijetTypes[udijet].Data()),1.0);
                if(udijet==iAcc) bHasDijet = true;
                dijet = dijets.at(udijet).at(0).at(0) + dijets.at(udijet).at(0).at(1);
                mjj = dijet.m();
                ptpair = dijet.pt();
                fHistname = TString::Format("%s/h_dijetInvM_%d_jetC%d",fGroupname.Data(), lCBin, udijet);
                fhistos.FillTH1(fHistname, mjj);
                fHistname = TString::Format("%s/h_dijetInvMTrunc_%d_jetC%d",fGroupname.Data(), lCBin, udijet);
                fhistos.FillTH1(fHistname, mjj);
                fHistname = TString::Format("%s/h_dijetInvMTrunc2_%d_jetC%d",fGroupname.Data(), lCBin, udijet);
                fhistos.FillTH1(fHistname, mjj);
                fHistname = TString::Format("%s/h_dijetPtPair_%d_jetC%d",fGroupname.Data(), lCBin, udijet);
                fhistos.FillTH1(fHistname, ptpair);
                dPhi = GetDeltaPhi(dijets.at(udijet).at(0).at(0), dijets.at(udijet).at(0).at(1));
                fHistname = TString::Format("%s/h_dijetDeltaPhi_%d_jetC%d",fGroupname.Data(), lCBin, udijet);
                fhistos.FillTH1(fHistname, dPhi);
            }

            // Analysis for dijet with deltaPhi cut.
            if(bHasDeltaPhiSubLeadJet && dijets.at(udijet).at(1).at(1).pt()>fSubleadingJetCut) {
                fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), lCBin);
                fhistos.FillTH1(fHistname, Form("%s deltaphi cut dijets",sDijetTypes[udijet].Data()),1.0);
                if(udijet==iAcc) bHasDeltaPhiDijet = true;
                dijet = dijets.at(udijet).at(1).at(0) + dijets.at(udijet).at(1).at(1);
                mjj = dijet.m();
                ptpair = dijet.pt();
                fHistname = TString::Format("%s/h_dijetInvMDeltaPhiCut_%d_jetC%d",fGroupname.Data(), lCBin, udijet);
                fhistos.FillTH1(fHistname, mjj);
                fHistname = TString::Format("%s/h_dijetInvMDeltaPhiCutTrunc_%d_jetC%d",fGroupname.Data(), lCBin, udijet);
                fhistos.FillTH1(fHistname, mjj);
                fHistname = TString::Format("%s/h_dijetInvMDeltaPhiCutTrunc2_%d_jetC%d",fGroupname.Data(), lCBin, udijet);
                fhistos.FillTH1(fHistname, mjj);
                fHistname = TString::Format("%s/h_dijetPtPairDeltaPhiCut_%d_jetC%d",fGroupname.Data(), lCBin, udijet);
                fhistos.FillTH1(fHistname, ptpair);
                dPhi = GetDeltaPhi(dijets.at(udijet).at(1).at(0), dijets.at(udijet).at(1).at(1));
                fHistname = TString::Format("%s/h_dijetDeltaPhiWithCut_%d_jetC%d",fGroupname.Data(), lCBin, udijet);
                fhistos.FillTH1(fHistname, dPhi);
            }
        }
    }

    // deltaM calculations here:
    if(bHasDeltaPhiDijet) {
        double firstJetPhi = dijets.at(iAcc).at(1).at(0).phi() - TMath::Pi(); //-pi to pi
        double firstConePhi = firstJetPhi-TMath::Pi()/2.0 < -TMath::Pi() ? firstJetPhi+3.0*TMath::Pi()/2.0 : firstJetPhi-TMath::Pi()/2.0;
        double coneDeltaPhi = GetDeltaPhi(dijets.at(iAcc).at(1).at(0), dijets.at(iAcc).at(1).at(1)); // 0-2pi
        double secondConePhi = firstConePhi+coneDeltaPhi > TMath::Pi() ? firstConePhi+coneDeltaPhi-2*TMath::Pi() : firstConePhi+coneDeltaPhi;
        double firstConeEta = dijets.at(iAcc).at(1).at(0).eta();
        double secondConeEta = dijets.at(iAcc).at(1).at(1).eta();
        fastjet::PseudoJet holderJet;
        fastjet::PseudoJet firstConeP;
        fastjet::PseudoJet secondConeP;
        fastjet::PseudoJet firstConeDeltaP;
        fastjet::PseudoJet secondConeDeltaP;
        pt = 0;
        pt2 = 0;
        for (utrack = 0; utrack < chparticles.size(); utrack++) {
            phi = chparticles.at(utrack).phi();
            eta = chparticles.at(utrack).eta();
            if(DeltaR(firstConeEta, eta, firstConePhi, phi) < fJetCone) {
                holderJet = fastjet::PseudoJet(chparticles.at(utrack).px(), chparticles.at(utrack).py(), chparticles.at(utrack).pz(), chparticles.at(utrack).E());
                firstConeP += holderJet;
                pt += holderJet.pt();
                    
            }
            if(DeltaR(secondConeEta, eta, secondConePhi, phi) < fJetCone) {
                holderJet = fastjet::PseudoJet(chparticles.at(utrack).px(), chparticles.at(utrack).py(), chparticles.at(utrack).pz(), chparticles.at(utrack).E());
                secondConeP += holderJet;
                pt2 += holderJet.pt();
            }
        }
        //Turn the cones back to jet direction
        firstConeP.reset_PtYPhiM(firstConeP.perp(),firstConeP.rap(),dijets.at(iAcc).at(1).at(0).phi(),firstConeP.m());
        secondConeP.reset_PtYPhiM(secondConeP.perp(),secondConeP.rap(),dijets.at(iAcc).at(1).at(1).phi(),secondConeP.m());

        jetAreaVector = dijets.at(iAcc).at(1).at(0).area_4vector();
        firstConeDeltaP = firstConeP - fastjet::PseudoJet(rho*jetAreaVector.px(),rho*jetAreaVector.py(),(rho+rhom)*jetAreaVector.pz(),(rho+rhom)*jetAreaVector.E());
        jetAreaVector = dijets.at(iAcc).at(1).at(1).area_4vector();
        secondConeDeltaP = secondConeP - fastjet::PseudoJet(rho*jetAreaVector.px(),rho*jetAreaVector.py(),(rho+rhom)*jetAreaVector.pz(),(rho+rhom)*jetAreaVector.E());
        fastjet::PseudoJet doubleDeltaCone = firstConeDeltaP + secondConeDeltaP;

        fHistname = TString::Format("%s/h_dijetdeltaM1", fGroupname.Data());
        fhistos.FillTH1(fHistname, doubleDeltaCone.m());

        // The second method
        dijet = dijets.at(iAcc).at(1).at(0) + dijets.at(iAcc).at(1).at(1);
        double doubSquared = doubleDeltaCone.E()*doubleDeltaCone.E() - doubleDeltaCone.px()*doubleDeltaCone.px() - doubleDeltaCone.py()*doubleDeltaCone.py() - doubleDeltaCone.pz()*doubleDeltaCone.pz();
        double dijetTimesDeltaCones = dijet.E()*doubleDeltaCone.E() - dijet.px()*doubleDeltaCone.px() - dijet.py()*doubleDeltaCone.py() - dijet.pz()*doubleDeltaCone.pz();

        double mass2;
        if(doubSquared + 2*dijetTimesDeltaCones > 0.0) mass2 = TMath::Sqrt(doubSquared + 2*dijetTimesDeltaCones);
        else mass2 = -TMath::Sqrt(-(doubSquared + 2*dijetTimesDeltaCones));
        fHistname = TString::Format("%s/h_dijetdeltaM2", fGroupname.Data());
        fhistos.FillTH1(fHistname, mass2);

        // The third method
        fastjet::PseudoJet doubleCone = firstConeP + secondConeP;
        fHistname = TString::Format("%s/h_doubleConeM", fGroupname.Data());
        fhistos.FillTH1(fHistname, doubleCone.m());
        double conesTimesDeltaCones = doubleCone.E()*doubleDeltaCone.E() - doubleCone.px()*doubleDeltaCone.px() - doubleCone.py()*doubleDeltaCone.py() - doubleCone.pz()*doubleDeltaCone.pz();

        double mass3;
        if(doubSquared + 2*conesTimesDeltaCones > 0.0) mass3 = TMath::Sqrt(doubSquared + 2*conesTimesDeltaCones);
        else mass3 = -TMath::Sqrt(-(doubSquared + 2*conesTimesDeltaCones));
        fHistname = TString::Format("%s/h_dijetdeltaM3", fGroupname.Data());
        fhistos.FillTH1(fHistname, mass3);

        // The fourth method
        double localRho1 = pt/(TMath::Pi()*0.4*0.4);
        double localRho2 = pt2/(TMath::Pi()*0.4*0.4);
        double deltaRho1 = rho-localRho1;
        double deltaRho2 = rho-localRho2;
        fastjet::PseudoJet deltaPmu1_fourth = deltaRho1*dijets.at(iAcc).at(1).at(0).area_4vector();
        fastjet::PseudoJet deltaPmu2_fourth = deltaRho2*dijets.at(iAcc).at(1).at(1).area_4vector();
        fastjet::PseudoJet doubleDeltaCone_fourth = deltaPmu1_fourth + deltaPmu2_fourth;

        fHistname = TString::Format("%s/h_dijetdeltaM4", fGroupname.Data());
        fhistos.FillTH1(fHistname, doubleDeltaCone_fourth.m());

        // THe fifth method
        deltaPmu1_fourth = localRho1*dijets.at(iAcc).at(1).at(0).area_4vector();
        deltaPmu2_fourth = localRho2*dijets.at(iAcc).at(1).at(1).area_4vector();
        fastjet::PseudoJet bgsubtrJeFifth1 = fastjet::PseudoJet(dijets.at(iAcc).at(1).at(0).px() - localRho1 * dijets.at(iAcc).at(1).at(0).area_4vector().px(),
                                                                dijets.at(iAcc).at(1).at(0).py() - localRho1 * dijets.at(iAcc).at(1).at(0).area_4vector().py(),
                                                                dijets.at(iAcc).at(1).at(0).pz() - localRho1 * dijets.at(iAcc).at(1).at(0).area_4vector().pz(),
                                                                dijets.at(iAcc).at(1).at(0).E()  - localRho1 * dijets.at(iAcc).at(1).at(0).area_4vector().E());
        fastjet::PseudoJet bgsubtrJeFifth2 = fastjet::PseudoJet(dijets.at(iAcc).at(1).at(1).px() - localRho2 * dijets.at(iAcc).at(1).at(1).area_4vector().px(),
                                                                dijets.at(iAcc).at(1).at(1).py() - localRho2 * dijets.at(iAcc).at(1).at(1).area_4vector().py(),
                                                                dijets.at(iAcc).at(1).at(1).pz() - localRho2 * dijets.at(iAcc).at(1).at(1).area_4vector().pz(),
                                                                dijets.at(iAcc).at(1).at(1).E()  - localRho2 * dijets.at(iAcc).at(1).at(1).area_4vector().E());
        doubleDeltaCone_fourth = bgsubtrJeFifth1 + bgsubtrJeFifth2;
        //doubleDeltaCone_fourth =   dijets.at(iAcc).at(1).at(0) - deltaPmu1_fourth
        //                         + dijets.at(iAcc).at(1).at(1) - deltaPmu2_fourth;

        dijet = dijets.at(iBGSubtr).at(1).at(0) + dijets.at(iBGSubtr).at(1).at(1);

        fDeltaM=dijet.m()-doubleDeltaCone_fourth.m();
        fHistname = TString::Format("%s/h_dijetdeltaM5", fGroupname.Data());
        fhistos.FillTH1(fHistname, fDeltaM);
        fHistname = TString::Format("%s/h_deltaMResponse", fGroupname.Data());
        fhistos.FillTH2(fHistname, dijet.m()+fDeltaM, dijet.m());

        double dBinCenter = 0;
        fHistname = TString::Format("%s/h_deltaMResponseEvery", fGroupname.Data());
        for(int iby = 1 ; iby <= ((TH1F*)fhistos.FindObject(fHistname))->GetNbinsX(); iby++){
            dBinCenter = ((TH1F*)fhistos.FindObject(fHistname))->GetYaxis()->GetBinCenter(iby);
            fhistos.FillTH2(fHistname, dBinCenter+fDeltaPt,dBinCenter);
        }
    }
    return;
}

// Response matrices are calculated in this function.
void AliJEmcalDijetAna::CalculateResponse(AliJEmcalDijetAna *anaDetMC, int iJetSetPart, int iJetSetDet) {
    if(fDebug>8) cout << "===== BEGIN RESPONSE CALC =====" << endl;
    
    vector<vector<fastjet::PseudoJet>> jetsDetMC = anaDetMC->GetJets();

    unsigned Njets = jets.at(iJetSetPart).size();
    unsigned NjetsDetMC = jetsDetMC.at(iJetSetDet).size();
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
    fHistname = TString::Format("%s/h_deltaPtResponseEvery", fGroupname.Data());
    for(int iby = 1 ; iby <= ((TH1F*)fhistos.FindObject(fHistname))->GetNbinsX(); iby++){
        dBinCenter = ((TH1F*)fhistos.FindObject(fHistname))->GetYaxis()->GetBinCenter(iby);
        fhistos.FillTH1(fHistname, dBinCenter+fDeltaPt,dBinCenter);
    }
    fHistname = TString::Format("%s/h_deltaPtResponseEvery_ALICE", fGroupname.Data());
    for(int iby = 1 ; iby <= ((TH1F*)fhistos.FindObject(fHistname))->GetNbinsX(); iby++){
        dBinCenter = ((TH1F*)fhistos.FindObject(fHistname))->GetYaxis()->GetBinCenter(iby);
        fhistos.FillTH2(fHistname, dBinCenter+fDeltaPt,dBinCenter);
    }
    for (ujet = 0; ujet < Njets; ujet++) { //True MC jets
        maxpt=0;
        maxptIndex=-1;
        bfound=false;
        deltaR=0;
        deltaRMatch=0;
        minR=999.0;
        ptTrue = jets.at(iJetSetPart).at(ujet).pt();
        fHistname = TString::Format("%s/h_deltaPtResponse", fGroupname.Data());
        fhistos.FillTH2(fHistname, ptTrue+fDeltaPt,ptTrue);
        fHistname = TString::Format("%s/h_deltaPtResponse_ALICE", fGroupname.Data());
        fhistos.FillTH2(fHistname, ptTrue+fDeltaPt,ptTrue);
        for (ujetDetMC = 0; ujetDetMC < NjetsDetMC; ujetDetMC++) { //Det MC jets
            deltaR = DeltaR(jets.at(iJetSetPart).at(ujet), jetsDetMC.at(iJetSetDet).at(ujetDetMC));
            if(deltaR<minR) minR=deltaR;
            if(deltaR < matchingR && jetsDetMC.at(iJetSetDet).at(ujetDetMC).pt() > maxpt) {
                maxpt = jetsDetMC.at(iJetSetDet).at(ujetDetMC).pt();
                maxptIndex = ujetDetMC;
                deltaRMatch = deltaR;
                bfound = true;
                //cout << "found, detPt vs truePt: " << maxpt << " <> " << jets[iJetSetPart][ujet].pt() << ", index: " << maxptIndex << endl;
            }
        }
        fHistname = TString::Format("%s/h_jetResponseDeltaRClosest", fGroupname.Data());
        fhistos.FillTH1(fHistname, minR);
        if(bfound) {
            ptDetMC = jetsDetMC.at(iJetSetDet).at(maxptIndex).pt();
            fHistname = TString::Format("%s/h_jetResponse", fGroupname.Data());
            fhistos.FillTH2(fHistname, ptDetMC, ptTrue);
            fHistname = TString::Format("%s/h_jetResponse_ALICE", fGroupname.Data());
            fhistos.FillTH2(fHistname, ptDetMC, ptTrue);
            fHistname = TString::Format("%s/h_jetResponseDeltaR", fGroupname.Data());
            fhistos.FillTH1(fHistname, deltaRMatch);
            fHistname = TString::Format("%s/h_jetResponseDeltaPt", fGroupname.Data());
            fhistos.FillTH1(fHistname, (ptTrue-ptDetMC)/ptTrue);
            fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
            fhistos.FillTH1(fHistname, "True jet has pair",1.0);
            bTrueJetMatch.at(ujet) = true;
            bDetJetMatch.at(maxptIndex) = true;
        } else {
            fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
            fhistos.FillTH1(fHistname, "True jet has no pair",1.0);
        }
    }
    for (ujetDetMC = 0; ujetDetMC < NjetsDetMC; ujetDetMC++) { //Det MC jets
        if(!bDetJetMatch.at(ujetDetMC)) {
            fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
            fhistos.FillTH1(fHistname, "Det jet has no pair",1.0);
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
    if(bHasDijet) {
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
        if(anaDetMC->HasDijet()) {

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
            if(bLeadingMatch && bSubleadingMatch) {
                fHistname = TString::Format("%s/h_dijetResponse", fGroupname.Data());
                fhistos.FillTH2(fHistname, dijetDetMC.m(), dijet.m());
                fHistname = TString::Format("%s/h_dijetResponseTrunc", fGroupname.Data());
                fhistos.FillTH2(fHistname, dijetDetMC.m(), dijet.m());
                fHistname = TString::Format("%s/h_dijetResponseTrunc2", fGroupname.Data());
                fhistos.FillTH2(fHistname, dijetDetMC.m(), dijet.m());
                fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
                fhistos.FillTH1(fHistname, "Dijet match",1.0);
            } else {
                fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
                fhistos.FillTH1(fHistname, "Dijet not match",1.0);
            }
        } else {
            fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
            fhistos.FillTH1(fHistname, "Dijet det not found",1.0);
        }
    } else {
        if(anaDetMC->HasDijet()) {
            //dijetDetMC = dijetsDetMC.at(iJetSetDet).at(0).at(0) + dijetsDetMC.at(iJetSetDet).at(0).at(1);
            fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
            fhistos.FillTH1(fHistname, "Dijet true not found",1.0);
        }
    }
    if(fDebug>8) cout << ((anaDetMC->HasDijet() && bLeadingMatch && bSubleadingMatch) ? "Dijet match" : "Dijet no match") << endl;

    // DeltaPhi cut dijet response.
    if(bHasDeltaPhiDijet) {
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
        if(anaDetMC->HasDeltaPhiDijet()) {
        
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

            if(bLeadingMatch && bSubleadingMatchDeltaPhi) {
                fHistname = TString::Format("%s/h_dijetResponseDeltaPhiCut", fGroupname.Data());
                fhistos.FillTH2(fHistname, dijetDetMC.m(), dijet.m());
                fHistname = TString::Format("%s/h_dijetResponseDeltaPhiCutTrunc", fGroupname.Data());
                fhistos.FillTH2(fHistname, dijetDetMC.m(), dijet.m());
                fHistname = TString::Format("%s/h_dijetResponseDeltaPhiCutTrunc2", fGroupname.Data());
                fhistos.FillTH2(fHistname, dijetDetMC.m(), dijet.m());
                fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
                fhistos.FillTH1(fHistname, "Dijet DPhi match",1.0);
            } else {
                fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
                fhistos.FillTH1(fHistname, "Dijet DPhi not match",1.0);
            }
        } else {
            fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
            fhistos.FillTH1(fHistname, "Dijet DPhi det not found",1.0);
        }
    } else {
        if(anaDetMC->HasDeltaPhiDijet()) {
            //dijetDetMC = dijetsDetMC.at(iJetSetDet).at(1).at(0) + dijetsDetMC.at(iJetSetDet).at(1).at(1);
            fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
            fhistos.FillTH1(fHistname, "Dijet DPhi true not found",1.0);
        }
    }
    if(fDebug>8) cout << ((anaDetMC->HasDeltaPhiDijet() && bLeadingMatch && bSubleadingMatchDeltaPhi) ? "Delta Phi Dijet match" : "Delta Phi Dijet no match") << endl;
    return;
}

void AliJEmcalDijetAna::ResetObjects() {
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

double AliJEmcalDijetAna::DeltaR(fastjet::PseudoJet jet1, fastjet::PseudoJet jet2) {
    return DeltaR(jet1.eta(), jet2.eta(), jet1.phi(), jet2.phi());
}

double AliJEmcalDijetAna::DeltaR(double eta1, double eta2, double phi1, double phi2) {
    double Deta = eta1 - eta2;
    double Dphi = TMath::Abs(phi1 - phi2);
    // Make sure that Dphi is in 0-pi range.
    Dphi = Dphi>TMath::Pi() ? 2*TMath::Pi()-Dphi : Dphi;

    return TMath::Sqrt(Deta*Deta + Dphi*Dphi);
}


bool AliJEmcalDijetAna::CheckDeltaPhi(fastjet::PseudoJet leadingJet, fastjet::PseudoJet subleadingJet, double deltaPhiCut) {
    double DeltaPhi;
    DeltaPhi = GetDeltaPhi(leadingJet, subleadingJet);
    if(TMath::Abs(DeltaPhi - TMath::Pi()) < deltaPhiCut) {
        return true;
    }
    return false;
}

double AliJEmcalDijetAna::GetDeltaPhi(fastjet::PseudoJet jet1, fastjet::PseudoJet jet2) {
    double DeltaPhi, DeltaPhi2;
    DeltaPhi  = jet2.delta_phi_to(jet1);
    DeltaPhi2 = DeltaPhi<0 ? DeltaPhi+TMath::TwoPi() : DeltaPhi;
    return DeltaPhi2;
}

// This should be done after SetSettings.
void AliJEmcalDijetAna::InitHistos(bool bIsMC) {
    cout << "Initing histograms for AliJEmcalDijetAna" << endl;
    fHistname = TString::Format("%s/h_info", fGroupname.Data());
    fhistos.FillTH1(fHistname, "Count", 1.0);
    fHistname = TString::Format("%s/h_info", fGroupname.Data());
    fhistos.FillTH1(fHistname, "MC", bIsMC);
    for(int i=0; i< fNcentBins; i++) {
        fHistname = TString::Format("%s/h_info", fGroupname.Data());
        fhistos.FillTH1(fHistname, Form("Cent bin border %02d",i), fNcentBins);
    }
    fHistname = TString::Format("%s/h_info", fGroupname.Data());
    fhistos.FillTH1(fHistname, "Jet cone", fJetCone);
    fHistname = TString::Format("%s/h_info", fGroupname.Data());
    fhistos.FillTH1(fHistname, "kt-jet cone", fktJetCone);
    fHistname = TString::Format("%s/h_info", fGroupname.Data());
    fhistos.FillTH1(fHistname, "kt scheme", ktScheme);
    fHistname = TString::Format("%s/h_info", fGroupname.Data());
    fhistos.FillTH1(fHistname, "antikt scheme", antiktScheme);
    fHistname = TString::Format("%s/h_info", fGroupname.Data());
    fhistos.FillTH1(fHistname, "Use pion mass", fusePionMass);
    fHistname = TString::Format("%s/h_info", fGroupname.Data());
    fhistos.FillTH1(fHistname, "Use DeltaPhi BG Subtr", fUseDeltaPhiBGSubtr);
    fHistname = TString::Format("%s/h_info", fGroupname.Data());
    fhistos.FillTH1(fHistname, "Particle eta cut", fParticleEtaCut);
    fHistname = TString::Format("%s/h_info", fGroupname.Data());
    fhistos.FillTH1(fHistname, "Particle pt cut", fParticlePtCut);
    fHistname = TString::Format("%s/h_info", fGroupname.Data());
    fhistos.FillTH1(fHistname, "Leading jet cut", fLeadingJetCut);
    fHistname = TString::Format("%s/h_info", fGroupname.Data());
    fhistos.FillTH1(fHistname, "Subleading jet cut", fSubleadingJetCut);
    fHistname = TString::Format("%s/h_info", fGroupname.Data());
    fhistos.FillTH1(fHistname, "Const. cut", fConstituentCut);
    fHistname = TString::Format("%s/h_info", fGroupname.Data());
    fhistos.FillTH1(fHistname, "Delta phi cut pi/",fDeltaPhiCut);
    fHistname = TString::Format("%s/h_info", fGroupname.Data());
    fhistos.FillTH1(fHistname, "Matching R for MC",matchingR);
    fHistname = TString::Format("%s/h_info", fGroupname.Data());
    fhistos.FillTH1(fHistname, "tracking ineff",ftrackingIneff);

    fHistname = TString::Format("%s/h_eventSel", fGroupname.Data());
    fhistos.FillTH1(fHistname, "events wo/ cuts",0.0);
    fHistname = TString::Format("%s/h_eventSel", fGroupname.Data());
    fhistos.FillTH1(fHistname, "catalyst entry ok",0.0);
    fHistname = TString::Format("%s/h_eventSel", fGroupname.Data());
    fhistos.FillTH1(fHistname, "catalyst ok",0.0);
    fHistname = TString::Format("%s/h_eventSel", fGroupname.Data());
    fhistos.FillTH1(fHistname, "vertex2013pA ok",0.0);
    fHistname = TString::Format("%s/h_eventSel", fGroupname.Data());
    fhistos.FillTH1(fHistname, "pileupSPD ok",0.0);
    fHistname = TString::Format("%s/h_eventSel", fGroupname.Data());
    fhistos.FillTH1(fHistname, "utils pileupSPD ok",0.0);
    fHistname = TString::Format("%s/h_eventSel", fGroupname.Data());
    fhistos.FillTH1(fHistname, "events",0.0);

    // Initialize fhistos.FillTH1(fh_events so that the bin order is correct
    for (int iBin=0; iBin < fNcentBins-1; iBin++) {
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "events",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "particles",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "acc. particles",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "no rho calc. events",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "rho calc. events",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "jets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "acc. jets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "const. cut jets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "bg. subtr. jets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "bg. subtr. const. cut jets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "kt-jets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "acc. kt-jets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "leading jet drop",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "subleading jet drop",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "raw dijets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "raw dijets leading cut",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "raw acc. dijets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "raw deltaphi cut dijets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "bg. subtr. dijets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "bg. subtr. dijets leading cut",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "bg. subtr. acc. dijets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "bg. subtr. deltaphi cut dijets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "bg. subtr. const. cut dijets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "bg. subtr. const. cut dijets leading cut",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "bg. subtr. const. cut acc. dijets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "bg. subtr. const. cut deltaphi cut dijets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "const. cut dijets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "const. cut dijets leading cut",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "const. cut acc. dijets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "const. cut deltaphi cut dijets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "kt dijets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "kt dijets leading cut",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "kt acc. dijets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "kt deltaphi cut dijets",0.0);
        fHistname = TString::Format("%s/h_events_%d", fGroupname.Data(), iBin);
        fhistos.FillTH1(fHistname, "pt_hard bin cuts",0.0);
        if(bIsMC) {
            fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
            fhistos.FillTH1(fHistname, "True jet has pair",0.0);
            fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
            fhistos.FillTH1(fHistname, "True jet has no pair",0.0);
            fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
            fhistos.FillTH1(fHistname, "Det jet has no pair",0.0);
            fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
            fhistos.FillTH1(fHistname, "Dijet match",0.0);
            fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
            fhistos.FillTH1(fHistname, "Dijet not match",0.0);
            fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
            fhistos.FillTH1(fHistname, "Dijet det not found",0.0);
            fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
            fhistos.FillTH1(fHistname, "Dijet true not found",0.0);
            fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
            fhistos.FillTH1(fHistname, "Dijet DPhi match",0.0);
            fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
            fhistos.FillTH1(fHistname, "Dijet DPhi not match",0.0);
            fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
            fhistos.FillTH1(fHistname, "Dijet DPhi det not found",0.0);
            fHistname = TString::Format("%s/h_responseInfo", fGroupname.Data());
            fhistos.FillTH1(fHistname, "Dijet DPhi true not found",0.0);
        }
    }
    return;
}

/*
 * This function allocates the histograms for basic tracking QA.
 * A set of histograms (pT, eta, phi, difference between kinematic properties
 * at the vertex and at the EMCal surface, number of tracks) is allocated
 * per each particle container and per each centrality bin.
 */
void AliJEmcalDijetAna::AllocateTrackHistograms()
{
    cout << "allo1" << endl;
    TString histname;
    TString histtitle;
    TString groupname;
    //AliParticleContainer* partCont = 0;
    TIter next(&fParticleCollArray);

    int NBINSJet=150;
    double LogBinsXJet[NBINSJet+1], LimLJet=0.1, LimHJet=500;
    double logBWJet = (log(LimHJet)-log(LimLJet))/NBINSJet;
    for(int ijetBin=0;ijetBin<=NBINSJet;ijetBin++) LogBinsXJet[ijetBin]=LimLJet*exp(ijetBin*logBWJet);

    int NBINSDijet=170;
    double logBinsXDijet[NBINSDijet+1], LimLDijet=0.1, LimHDijet=1000;
    double logBWDijet = (log(LimHDijet)-log(LimLDijet))/NBINSDijet;
    for(int iDijet=0;iDijet<=NBINSDijet;iDijet++) logBinsXDijet[iDijet]=LimLDijet*exp(iDijet*logBWDijet);

    int NBINSAlice=300;
    double NBINLowAlice=10;
    double NBINHighAlice=310;

    /*
       hBinning = new TH1F("hBinning","hBinning", NBINSDijet, logBinsXDijet);
       hBinningALICE = new TH1F("hBinningALICE","hBinningALICE", NBINSAlice, NBINLowAlice, NBINHighAlice);
       fana->SetBinning(hBinning);
       fana->SetBinningALICE(hBinningALICE);
       */
    cout << "allo2" << endl;

    int fJetBin=5;
    //while ((partCont = static_cast<AliParticleContainer*>(next()))) {
        groupname = "dijet";//partCont->GetName();
        // Protect against creating the histograms twice
        cout << "allo21" << endl;
        if (fhistos.FindObject(groupname)) {
            AliWarning(TString::Format("%s: Found groupname %s in hist manager. The track containers will be filled into the same histograms.", GetName(), groupname.Data()));
            return;
        }
        
        cout << "allo22 group " << groupname.Data() << endl;
        fhistos.CreateHistoGroup(groupname);
        cout << "allo3" << endl;

        histname = TString::Format("%s/h_eventSel", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH1(histname, histtitle, 10, 0.0, 10.0 );

        histname = TString::Format("%s/h_info", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH1(histname, histtitle, 40, 0.0, 40.0 );

        histname = TString::Format("%s/h_centrality", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH1(histname, histtitle, 100, 0.0, 100.0 );

        histname = TString::Format("%s/h_zvtx", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH1(histname, histtitle, 40, -20.0, 20.0 );

        for (Int_t cent = 0; cent < fNcentBins; cent++) {
            histname = TString::Format("%s/histTrackPt_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;#it{p}_{T,track} (GeV/#it{c});counts", histname.Data());
            fhistos.CreateTH1(histname, histtitle, 200 / 2, 1, 200 / 2);

            /*
               histname = TString::Format("%s/histTrackPhi_%d", groupname.Data(), cent);
               histtitle = TString::Format("%s;#it{#phi}_{track};counts", histname.Data());
               fhistos.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

               histname = TString::Format("%s/histTrackEta_%d", groupname.Data(), cent);
               histtitle = TString::Format("%s;#it{#eta}_{track};counts", histname.Data());
               fhistos.CreateTH1(histname, histtitle, fNbins / 6, -1, 1);

               if (TClass(partCont->GetClassName()).InheritsFrom("AliVTrack")) {
               histname = TString::Format("%s/fHistDeltaEtaPt_%d", groupname.Data(), cent);
               histtitle = TString::Format("%s;#it{p}_{T,track}^{vertex} (GeV/#it{c});#it{#eta}_{track}^{vertex} - #it{#eta}_{track}^{EMCal};counts", histname.Data());
               fhistos.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, 50, -0.5, 0.5);

               histname = TString::Format("%s/fHistDeltaPhiPt_%d", groupname.Data(), cent);
               histtitle = TString::Format("%s;#it{p}_{T,track}^{vertex} (GeV/#it{c});#it{#phi}_{track}^{vertex} - #it{#phi}_{track}^{EMCal};counts", histname.Data());
               fhistos.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, 200, -2, 2);

               histname = TString::Format("%s/fHistDeltaPtvsPt_%d", groupname.Data(), cent);
               histtitle = TString::Format("%s;#it{p}_{T,track}^{vertex} (GeV/#it{c});#it{p}_{T,track}^{vertex} - #it{p}_{T,track}^{EMCal} (GeV/#it{c});counts", histname.Data());
               fhistos.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, fNbins / 2, -fMaxBinPt/2, fMaxBinPt/2);

               histname = TString::Format("%s/fHistEoverPvsP_%d", groupname.Data(), cent);
               histtitle = TString::Format("%s;#it{P}_{track} (GeV/#it{c});#it{E}_{cluster} / #it{P}_{track} #it{c};counts", histname.Data());
               fhistos.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, fNbins / 2, 0, 4);
               }

               histname = TString::Format("%s/histNTracks_%d", groupname.Data(), cent);
               histtitle = TString::Format("%s;number of tracks;events", histname.Data());
               if (fForceBeamType != kpp) {
               fhistos.CreateTH1(histname, histtitle, 500, 0, 5000);
               }
               else {
               fhistos.CreateTH1(histname, histtitle, 200, 0, 200);
               }
               */
            histname = TString::Format("%s/h_events_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH1(histname, histtitle, 40, 0.0, 40.0 );

            histname = TString::Format("%s/h_nch_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH1(histname, histtitle, 101, -0.5, 100.5 );

            // ============= CHARGED PARTICLE HISTOS ============= 
            histname = TString::Format("%s/h_pt_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH1(histname, histtitle, NBINSJet, LogBinsXJet );

            histname = TString::Format("%s/h_ptPosEta_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH1(histname, histtitle,NBINSJet, LogBinsXJet );

            histname = TString::Format("%s/h_ptNegEta_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH1(histname, histtitle,NBINSJet, LogBinsXJet );

            histname = TString::Format("%s/h_eta_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH1(histname, histtitle, 100, -1.0, 1.0 );

            histname = TString::Format("%s/h_phi_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH1(histname, histtitle, 100, -TMath::Pi(), TMath::Pi());

            histname = TString::Format("%s/h_etaPhi_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH2(histname, histtitle, 100, -1.0, 1.0, 100, -TMath::Pi(), TMath::Pi());

            // ============= JET HISTOS ============= 
            histname = TString::Format("%s/h_rho_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH1(histname, histtitle, NBINSJet, LogBinsXJet);

            histname = TString::Format("%s/h_rhoHighPt_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH1(histname, histtitle, NBINSJet, LogBinsXJet);

            histname = TString::Format("%s/h_rhom_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH1(histname, histtitle, NBINSJet, LogBinsXJet);

            histname = TString::Format("%s/h_rhomHighPt_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH1(histname, histtitle, NBINSJet, LogBinsXJet);

            histname = TString::Format("%s/h_rhoLin_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH1(histname, histtitle, 501, -0.1, 100.1);

            histname = TString::Format("%s/h_rhoLinHighPt_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH1(histname, histtitle, 501, -0.1, 100.1);

            histname = TString::Format("%s/h_rhomLin_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH1(histname, histtitle, 501, -0.1, 100.1);

            histname = TString::Format("%s/h_rhomLinHighPt_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH1(histname, histtitle, 501, -0.1, 100.1);

            histname = TString::Format("%s/h_jetN_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH1(histname, histtitle, 51, -0.5, 50.5 );

            histname = TString::Format("%s/h_randConeEtaPhi_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH2(histname, histtitle, 100, -1.0, 1.0, 100, -TMath::Pi(), TMath::Pi());

            for (Int_t jetbin = 0; jetbin < fJetBin; jetbin++) {
                histname = TString::Format("%s/h_jetPt_%d_jetC%d", groupname.Data(), cent, jetbin);
                histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
                fhistos.CreateTH1(histname, histtitle,NBINSJet, LogBinsXJet );

                histname = TString::Format("%s/h_jetPt_ALICE_%d_jetC%d", groupname.Data(), cent, jetbin);
                histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
                fhistos.CreateTH1(histname, histtitle, NBINSAlice, NBINLowAlice, NBINHighAlice );

                histname = TString::Format("%s/h_jetPtTransBGSub_%d_jetC%d", groupname.Data(), cent, jetbin);
                histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
                fhistos.CreateTH1(histname, histtitle, NBINSAlice, NBINLowAlice, NBINHighAlice );

                histname = TString::Format("%s/h_jetEta_%d_jetC%d", groupname.Data(), cent, jetbin);
                histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
                fhistos.CreateTH1(histname, histtitle, 100, -1.0, 1.0);

                histname = TString::Format("%s/h_jetPhi_%d_jetC%d", groupname.Data(), cent, jetbin);
                histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
                fhistos.CreateTH1(histname, histtitle, 100, -TMath::Pi(), TMath::Pi());

                histname = TString::Format("%s/h_jetEtaPhi_%d_jetC%d", groupname.Data(), cent, jetbin);
                histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
                fhistos.CreateTH2(histname, histtitle, 100, -1.0, 1.0, 100, -TMath::Pi(), TMath::Pi());

                histname = TString::Format("%s/h_jetArea_%d_jetC%d", groupname.Data(), cent, jetbin);
                histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
                fhistos.CreateTH1(histname, histtitle, NBINSJet, LogBinsXJet );

                histname = TString::Format("%s/h_jetAreaRho_%d_jetC%d", groupname.Data(), cent, jetbin);
                histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
                fhistos.CreateTH1(histname, histtitle, NBINSJet, LogBinsXJet );
            }

            histname = TString::Format("%s/h_deltaPt_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH1(histname, histtitle, 321, -20.5, 300.5);

            histname = TString::Format("%s/h_maxJetptOverPtHard_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH1(histname, histtitle, 100, 0, 20);

            histname = TString::Format("%s/h_ptHard_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH1(histname, histtitle,NBINSJet, LogBinsXJet );

            // ============= DIJET HISTOS ============= 
            for (Int_t jetbin = 0; jetbin < fJetBin; jetbin++) {
                histname = TString::Format("%s/h_dijetInvM_%d_jetC%d", groupname.Data(), cent, jetbin);
                histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
                fhistos.CreateTH1(histname, histtitle, NBINSDijet, logBinsXDijet);

                histname = TString::Format("%s/h_dijetInvMTrunc_%d_jetC%d", groupname.Data(), cent, jetbin);
                histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
                fhistos.CreateTH1(histname, histtitle, 50, 30, 280);

                histname = TString::Format("%s/h_dijetInvMTrunc2_%d_jetC%d", groupname.Data(), cent, jetbin);
                histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
                fhistos.CreateTH1(histname, histtitle, 100, 30, 530);

                histname = TString::Format("%s/h_dijetPtPair_%d_jetC%d", groupname.Data(), cent, jetbin);
                histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
                fhistos.CreateTH1(histname, histtitle, NBINSDijet, logBinsXDijet );

                histname = TString::Format("%s/h_dijetDeltaPhi_%d_jetC%d", groupname.Data(), cent, jetbin);
                histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
                fhistos.CreateTH1(histname, histtitle, 100, 0, 10);

                histname = TString::Format("%s/h_dijetPtPairDeltaPhiCut_%d_jetC%d", groupname.Data(), cent, jetbin);
                histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
                fhistos.CreateTH1(histname, histtitle, NBINSDijet, logBinsXDijet );

                histname = TString::Format("%s/h_dijetInvMDeltaPhiCut_%d_jetC%d", groupname.Data(), cent, jetbin);
                histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
                fhistos.CreateTH1(histname, histtitle, NBINSDijet, logBinsXDijet);

                histname = TString::Format("%s/h_dijetInvMDeltaPhiCutTrunc_%d_jetC%d", groupname.Data(), cent, jetbin);
                histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
                fhistos.CreateTH1(histname, histtitle, 50, 30, 280);

                histname = TString::Format("%s/h_dijetInvMDeltaPhiCutTrunc2_%d_jetC%d", groupname.Data(), cent, jetbin);
                histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
                fhistos.CreateTH1(histname, histtitle, 100, 30, 530);

                histname = TString::Format("%s/h_dijetDeltaPhiWithCut_%d_jetC%d", groupname.Data(), cent, jetbin);
                histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
                fhistos.CreateTH1(histname, histtitle, 100, 0, 10);
            }
        }


        // ============ Response histograms ===========
        histname = TString::Format("%s/h_responseInfo", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH1(histname, histtitle, 40, 0.0, 40.0 );

        histname = TString::Format("%s/h_jetResponseDeltaR", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH1(histname, histtitle, 100, 0.0, 1.0);

        histname = TString::Format("%s/h_jetResponseDeltaRClosest", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH1(histname, histtitle, 100, 0.0, 1.0);

        histname = TString::Format("%s/h_jetResponseDeltaPt", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH1(histname, histtitle, 200, -2.0, 1.0);

        for (Int_t cent = 0; cent < fNcentBins; cent++) {
            for (Int_t jetbin = 0; jetbin < fJetBin; jetbin++) {
                histname = TString::Format("%s/h_jetDeltaRMin_%d_jetC%d", groupname.Data(), cent, jetbin);
                histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
                fhistos.CreateTH1(histname, histtitle, 400, 0.0, 4.0);
            }

            histname = TString::Format("%s/h_jetBGSubtrDeltaR_%d", groupname.Data(), cent);
            histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
            fhistos.CreateTH1(histname, histtitle, 400, 0.0, 4.0);
        }

        histname = TString::Format("%s/h_jetResponse", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH2(histname, histtitle, NBINSJet, LogBinsXJet, NBINSJet, LogBinsXJet );

        histname = TString::Format("%s/h_jetResponse_ALICE", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH2(histname, histtitle, NBINSAlice, NBINLowAlice, NBINHighAlice, NBINSAlice, NBINLowAlice, NBINHighAlice );

        histname = TString::Format("%s/h_deltaPtResponse", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH2(histname, histtitle, NBINSJet, LogBinsXJet, NBINSJet, LogBinsXJet );

        histname = TString::Format("%s/h_deltaPtResponse_ALICE", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH2(histname, histtitle, NBINSAlice, NBINLowAlice, NBINHighAlice, NBINSAlice, NBINLowAlice, NBINHighAlice );

        histname = TString::Format("%s/h_deltaPtResponseEvery", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH2(histname, histtitle, NBINSJet, LogBinsXJet, NBINSJet, LogBinsXJet );

        histname = TString::Format("%s/h_deltaPtResponseEvery_ALICE", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH2(histname, histtitle, NBINSAlice, NBINLowAlice, NBINHighAlice, NBINSAlice, NBINLowAlice, NBINHighAlice );

        histname = TString::Format("%s/h_dijetResponse", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH2(histname, histtitle, NBINSDijet, logBinsXDijet, NBINSDijet, logBinsXDijet );

        histname = TString::Format("%s/h_doubleConeM", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH1(histname, histtitle, 500, 0, 500);

        histname = TString::Format("%s/h_dijetdeltaM1", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH1(histname, histtitle, 500, -250, 250);

        histname = TString::Format("%s/h_dijetdeltaM2", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH1(histname, histtitle, 500, -250, 250);

        histname = TString::Format("%s/h_dijetdeltaM3", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH1(histname, histtitle, 500, -250, 250);

        histname = TString::Format("%s/h_dijetdeltaM4", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH1(histname, histtitle, 500, -250, 250);

        histname = TString::Format("%s/h_dijetdeltaM5", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH1(histname, histtitle, 500, -250, 250);

        histname = TString::Format("%s/h_deltaMResponse", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH2(histname, histtitle, NBINSDijet, logBinsXDijet, NBINSDijet, logBinsXDijet );

        histname = TString::Format("%s/h_deltaMResponseEvery", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH2(histname, histtitle, NBINSDijet, logBinsXDijet, NBINSDijet, logBinsXDijet );

        histname = TString::Format("%s/h_dijetResponseTrunc", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH2(histname, histtitle, 50, 30, 280, 50, 30, 280);

        histname = TString::Format("%s/h_dijetResponseTrunc2", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH2(histname, histtitle, 100, 30, 530, 100, 30, 530);

        histname = TString::Format("%s/h_dijetResponseDeltaPhiCut", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH2(histname, histtitle, NBINSDijet, logBinsXDijet, NBINSDijet, logBinsXDijet );

        histname = TString::Format("%s/h_dijetResponseDeltaPhiCutTrunc", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH2(histname, histtitle, 50, 30, 280, 50, 30, 280);

        histname = TString::Format("%s/h_dijetResponseDeltaPhiCutTrunc2", groupname.Data());
        histtitle = TString::Format("%s;x_axis;y_axis", histname.Data());
        fhistos.CreateTH2(histname, histtitle, 100, 30, 530, 100, 30, 530);

    //}
}




#endif
