//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for MC analysis
//  - MC output
//  implementation file
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________

#include "AliAnalysisEtMonteCarlo.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDtrack.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TParticle.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "TList.h"
#include "AliESDCaloCluster.h"
#include "AliLog.h"
#include <iostream>
#include <AliCentrality.h>

using namespace std;

ClassImp(AliAnalysisEtMonteCarlo);


// ctor
AliAnalysisEtMonteCarlo::AliAnalysisEtMonteCarlo():AliAnalysisEt()
        ,fImpactParameter(0)
        ,fNcoll(0)
        ,fNpart(0)

        ,fHistDecayVertexNonRemovedCharged(0)
        ,fHistDecayVertexRemovedCharged(0)
        ,fHistDecayVertexNonRemovedNeutral(0)
        ,fHistDecayVertexRemovedNeutral(0)

        ,fHistRemovedOrNot(0)
        ,fHistEtNonRemovedProtons(0)
        ,fHistEtNonRemovedAntiProtons(0)
        ,fHistEtNonRemovedPiPlus(0)
        ,fHistEtNonRemovedPiMinus(0)
        ,fHistEtNonRemovedKaonPlus(0)
        ,fHistEtNonRemovedKaonMinus(0)
        ,fHistEtNonRemovedK0s(0)
        ,fHistEtNonRemovedLambdas(0)
        ,fHistEtNonRemovedElectrons(0)
        ,fHistEtNonRemovedPositrons(0)
        ,fHistEtNonRemovedMuPlus(0)
        ,fHistEtNonRemovedMuMinus(0)
        ,fHistEtNonRemovedNeutrons(0)
        ,fHistEtNonRemovedAntiNeutrons(0)
        ,fHistEtNonRemovedGammas(0)
        ,fHistEtNonRemovedGammasFromPi0(0)
        ,fHistEtRemovedGammas(0)
        ,fHistEtRemovedNeutrons(0)
        ,fHistEtRemovedAntiNeutrons(0)
        ,fHistMultNonRemovedProtons(0)
        ,fHistMultNonRemovedAntiProtons(0)
        ,fHistMultNonRemovedPiPlus(0)
        ,fHistMultNonRemovedPiMinus(0)
        ,fHistMultNonRemovedKaonPlus(0)
        ,fHistMultNonRemovedKaonMinus(0)
        ,fHistMultNonRemovedK0s(0)
        ,fHistMultNonRemovedLambdas(0)
        ,fHistMultNonRemovedElectrons(0)
        ,fHistMultNonRemovedPositrons(0)
        ,fHistMultNonRemovedMuPlus(0)
        ,fHistMultNonRemovedMuMinus(0)
        ,fHistMultNonRemovedNeutrons(0)
        ,fHistMultNonRemovedAntiNeutrons(0)
        ,fHistMultNonRemovedGammas(0)
        ,fHistMultRemovedGammas(0)
        ,fHistMultRemovedNeutrons(0)
        ,fHistMultRemovedAntiNeutrons(0)
        ,fHistTrackMultvsNonRemovedCharged(0)
        ,fHistTrackMultvsNonRemovedNeutral(0)
        ,fHistTrackMultvsRemovedGamma(0)
        ,fHistClusterMultvsNonRemovedCharged(0)
        ,fHistClusterMultvsNonRemovedNeutral(0)
        ,fHistClusterMultvsRemovedGamma(0)
        ,fHistMultvsNonRemovedChargedE(0)
        ,fHistMultvsNonRemovedNeutralE(0)
        ,fHistMultvsRemovedGammaE(0)
        ,fEtNonRemovedProtons(0)
        ,fEtNonRemovedAntiProtons(0)
        ,fEtNonRemovedPiPlus(0)
        ,fEtNonRemovedPiMinus(0)
        ,fEtNonRemovedKaonPlus(0)
        ,fEtNonRemovedKaonMinus(0)
        ,fEtNonRemovedK0s(0)
        ,fEtNonRemovedLambdas(0)
        ,fEtNonRemovedElectrons(0)
        ,fEtNonRemovedPositrons(0)
        ,fEtNonRemovedMuMinus(0)
        ,fEtNonRemovedMuPlus(0)
        ,fEtNonRemovedGammas(0)
        ,fEtNonRemovedGammasFromPi0(0)
        ,fEtNonRemovedNeutrons(0)
        ,fEtNonRemovedAntiNeutrons(0)
        ,fEtRemovedGammas(0)
        ,fEtRemovedNeutrons(0)
        ,fEtRemovedAntiNeutrons(0)
        ,fMultNonRemovedProtons(0)
        ,fMultNonRemovedAntiProtons(0)
        ,fMultNonRemovedPiPlus(0)
        ,fMultNonRemovedPiMinus(0)
        ,fMultNonRemovedKaonPlus(0)
        ,fMultNonRemovedKaonMinus(0)
        ,fMultNonRemovedK0s(0)
        ,fMultNonRemovedLambdas(0)
        ,fMultNonRemovedElectrons(0)
        ,fMultNonRemovedPositrons(0)
        ,fMultNonRemovedMuMinus(0)
        ,fMultNonRemovedMuPlus(0)
        ,fMultNonRemovedGammas(0)
        ,fMultNonRemovedNeutrons(0)
        ,fMultNonRemovedAntiNeutrons(0)
        ,fMultRemovedGammas(0)
        ,fMultRemovedNeutrons(0)
        ,fMultRemovedAntiNeutrons(0)
        ,fTrackMultInAcc(0)
        ,fHistDxDzNonRemovedCharged(0)
        ,fHistDxDzRemovedCharged(0)
        ,fHistDxDzNonRemovedNeutral(0)
        ,fHistDxDzRemovedNeutral(0)
	,fHistPiPlusMult(0)
	,fHistPiMinusMult(0)
	,fHistPiZeroMult(0)
	,fHistPiPlusMultAcc(0)
	,fHistPiMinusMultAcc(0)
	,fHistPiZeroMultAcc(0)
	,fPiPlusMult(0)
	,fPiMinusMult(0)
	,fPiZeroMult(0)
	,fPiPlusMultAcc(0)
	,fPiMinusMultAcc(0)
	,fPiZeroMultAcc(0)
        ,fNeutralRemoved(0)
        ,fChargedRemoved(0)
        ,fChargedNotRemoved(0)
        ,fNeutralNotRemoved(0)
        ,fEnergyNeutralRemoved(0)
        ,fEnergyChargedRemoved(0)
        ,fEnergyChargedNotRemoved(0)
        ,fEnergyNeutralNotRemoved(0)

{
    fTrackDistanceCut = 7.0;
}

// dtor
AliAnalysisEtMonteCarlo::~AliAnalysisEtMonteCarlo()
{//Destructor
	delete fHistDecayVertexNonRemovedCharged; // Decay vertex for non-removed charged particles
	delete fHistDecayVertexRemovedCharged; // Decay vertex for non-removed charged particles
	delete fHistDecayVertexNonRemovedNeutral; // Decay vertex for non-removed charged particles
	delete fHistDecayVertexRemovedNeutral; // Decay vertex for non-removed charged particles
	
	delete fHistRemovedOrNot; // If charged/neutral particles were removed or not
	
	delete fHistEtNonRemovedProtons; // enter comment here
	delete fHistEtNonRemovedAntiProtons; // enter comment here
	delete fHistEtNonRemovedPiPlus; // enter comment here
	delete fHistEtNonRemovedPiMinus; // enter comment here
	delete fHistEtNonRemovedKaonPlus; // enter comment here
	delete fHistEtNonRemovedKaonMinus; // enter comment here
	delete fHistEtNonRemovedK0s; // enter comment here
	delete fHistEtNonRemovedLambdas; // enter comment here
	delete fHistEtNonRemovedElectrons; // enter comment here
	delete fHistEtNonRemovedPositrons; // enter comment here
	delete fHistEtNonRemovedMuPlus; // enter comment here
	delete fHistEtNonRemovedMuMinus; // enter comment here
	delete fHistEtNonRemovedNeutrons; // enter comment here
	delete fHistEtNonRemovedAntiNeutrons; // enter comment here
	delete fHistEtNonRemovedGammas; // enter comment here
	delete fHistEtNonRemovedGammasFromPi0; // enter comment here
	
	delete fHistEtRemovedGammas; // enter comment here
	delete fHistEtRemovedNeutrons; // enter comment here
	delete fHistEtRemovedAntiNeutrons; // enter comment here
	
	
	delete fHistMultNonRemovedProtons; // enter comment here 
	delete fHistMultNonRemovedAntiProtons; // enter comment here 
	delete fHistMultNonRemovedPiPlus; // enter comment here 
	delete fHistMultNonRemovedPiMinus; // enter comment here 
	delete fHistMultNonRemovedKaonPlus; // enter comment here 
	delete fHistMultNonRemovedKaonMinus; // enter comment here 
	delete fHistMultNonRemovedK0s; // enter comment here 
	delete fHistMultNonRemovedLambdas; // enter comment here 
	delete fHistMultNonRemovedElectrons; // enter comment here 
	delete fHistMultNonRemovedPositrons; // enter comment here 
	delete fHistMultNonRemovedMuPlus; // enter comment here
	delete fHistMultNonRemovedMuMinus; // enter comment here
	delete fHistMultNonRemovedNeutrons; // enter comment here
	delete fHistMultNonRemovedAntiNeutrons; // enter comment here
	delete fHistMultNonRemovedGammas; // enter comment here
	
	delete fHistMultRemovedGammas; // enter comment here
	delete fHistMultRemovedNeutrons; // enter comment here
	delete fHistMultRemovedAntiNeutrons; // enter comment here
	
	delete fHistTrackMultvsNonRemovedCharged; // enter comment here
	delete fHistTrackMultvsNonRemovedNeutral; // enter comment here
	delete fHistTrackMultvsRemovedGamma; // enter comment here
	
	delete fHistClusterMultvsNonRemovedCharged; // enter comment here
	delete fHistClusterMultvsNonRemovedNeutral; // enter comment here
	delete fHistClusterMultvsRemovedGamma; // enter comment here
	
	delete fHistMultvsNonRemovedChargedE; // enter comment here
	delete fHistMultvsNonRemovedNeutralE; // enter comment here
	delete fHistMultvsRemovedGammaE; // enter comment here
	delete fHistDxDzNonRemovedCharged; // enter comment here
	delete fHistDxDzRemovedCharged; // enter comment here
	delete fHistDxDzNonRemovedNeutral; // enter comment here
	delete fHistDxDzRemovedNeutral; // enter comment here
	
	delete fHistPiPlusMult; // enter comment here
	delete fHistPiMinusMult; // enter comment here
	delete fHistPiZeroMult; // enter comment here

	delete fHistPiPlusMultAcc; // enter comment here
	delete fHistPiMinusMultAcc; // enter comment here
	delete fHistPiZeroMultAcc; // enter comment here
}

Int_t AliAnalysisEtMonteCarlo::AnalyseEvent(AliVEvent* ev)
{ // analyse MC event
    ResetEventValues();

    fPiPlusMult = 0;
    fPiMinusMult = 0;
    fPiZeroMult = 0;
    if (fCentrality)
    {
        fCentClass = fCentrality->GetCentralityClass10(fCentralityMethod);

    }
    fSparseTracks[0] = 0;
    fSparseTracks[1] = 0;
    fSparseTracks[2] = 0;
    fSparseTracks[3] = 0;
    fSparseTracks[4] = 0;
    fSparseTracks[5] = 0;
    fSparseTracks[6] = fCentClass;


    // Get us an mc event
    if (!ev) {
        AliFatal("ERROR: Event does not exist");
        return 0;
    }
    AliMCEvent *event = dynamic_cast<AliMCEvent*>(ev);
    if (!event) {
        AliFatal("ERROR: MC Event does not exist");
        return 0;
    }

    Double_t protonMass =fgProtonMass;

    // Hijing header
    AliGenEventHeader* genHeader = event->GenEventHeader();
    if (!genHeader) {
        Printf("ERROR: Event generation header does not exist");
        return 0;
    }
    AliGenHijingEventHeader* hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(genHeader);
    if (hijingGenHeader) {
        fImpactParameter = hijingGenHeader->ImpactParameter();
        fNcoll = hijingGenHeader->HardScatters(); // or should this be some combination of NN() NNw() NwN() NwNw() ?
        fNpart = hijingGenHeader->ProjectileParticipants() + hijingGenHeader->TargetParticipants();
        /*
          printf("Hijing: ImpactParameter %g ReactionPlaneAngle %g \n",
          hijingGenHeader->ImpactParameter(), hijingGenHeader->ReactionPlaneAngle());
          printf("HardScatters %d ProjecileParticipants %d TargetParticipants %d\n",
          hijingGenHeader->HardScatters(), hijingGenHeader->ProjectileParticipants(), hijingGenHeader->TargetParticipants());
          printf("ProjSpectatorsn %d ProjSpectatorsp %d TargSpectatorsn %d TargSpectatorsp %d\n",
          hijingGenHeader->ProjSpectatorsn(), hijingGenHeader->ProjSpectatorsp(), hijingGenHeader->TargSpectatorsn(), hijingGenHeader->TargSpectatorsp());
          printf("NN %d NNw %d NwN %d, NwNw %d\n",
          hijingGenHeader->NN(), hijingGenHeader->NNw(), hijingGenHeader->NwN(), hijingGenHeader->NwNw());
        */
    }

    /* // placeholder if we want to get some Pythia info later
       AliGenPythiaEventHeader* pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(genHeader);
       if (pythiaGenHeader) { // not Hijing; try with Pythia
       printf("Pythia: ProcessType %d  GetPtHard %g \n",
       pythiaGenHeader->ProcessType(), pythiaGenHeader->GetPtHard());
       }
    */

    // Let's play with the stack!
    AliStack *stack = event->Stack();

    Int_t nPrim = stack->GetNtrack();

    Int_t partCount = 0;
    for (Int_t iPart = 0; iPart < nPrim; iPart++)
    {

        TParticle *part = stack->Particle(iPart);
//         TParticle *partMom = 0;
//         TParticle *partMomLast = 0;

        if (!part)
        {
            Printf("ERROR: Could not get particle %d", iPart);
            continue;
        }

//         Int_t iPartMom = part->GetMother(0);
//         Int_t iPartLastMom = part->GetMother(1);

        TParticlePDG *pdg = part->GetPDG(0);
        //TParticlePDG *pdgMom = 0;
//         TParticlePDG *pdgMomLast = 0;

        if (!pdg)
        {
            //Printf("ERROR: Could not get particle PDG %d", iPart);
            continue;
        }

//         if (iPartMom>0)
//         {
//             partMom = stack->Particle(iPartMom);
//             //pdgMom = partMom->GetPDG(0);
//         }

//         if (iPartLastMom>0)
//         {
//             partMomLast = stack->Particle(iPartLastMom);
//             pdgMomLast = partMomLast->GetPDG(0);
//         }


        Double_t particleMassPart = 0; //The mass part in the Et calculation for this particle

        // Check if it is a primary particle
        //if (!stack->IsPhysicalPrimary(iPart)) continue;

        //printf("MC: iPart %03d eta %4.3f phi %4.3f code %d charge %g \n", iPart, part->Eta(), part->Phi(), pdg->PdgCode(), pdg->Charge()); // tmp/debug printout

        // Check for reasonable (for now neutral and singly charged) charge on the particle
        //TODO:Maybe not only singly charged?
        if ((stack->IsPhysicalPrimary(iPart))&&( TMath::Abs(TMath::Abs(pdg->Charge()) - fCuts->GetMonteCarloSingleChargedParticle())<1e-3 || TMath::Abs(TMath::Abs(pdg->Charge()) - fCuts->GetMonteCarloNeutralParticle())<1e-3))
        {

            fMultiplicity++;

            if (TMath::Abs(part->Eta()) < fCuts->GetCommonEtaCut())
            {
                // calculate E_T
                if (
                    TMath::Abs(pdg->PdgCode()) == fgProtonCode ||
                    TMath::Abs(pdg->PdgCode()) == fgNeutronCode ||
                    TMath::Abs(pdg->PdgCode()) == fgLambdaCode ||
                    TMath::Abs(pdg->PdgCode()) == fgXiCode ||
                    TMath::Abs(pdg->PdgCode()) == fgXi0Code ||
                    TMath::Abs(pdg->PdgCode()) == fgOmegaCode
                )
                {
                    if (pdg->PdgCode() > 0) {
                        particleMassPart = - protonMass;
                    }
                    if (pdg->PdgCode() < 0) {
                        particleMassPart = protonMass;
                    }
                }
                Double_t et = part->Energy() * TMath::Sin(part->Theta()) + particleMassPart;

                fSparseTracks[0] = pdg->PdgCode();
                fSparseTracks[1] = pdg->Charge()*3;
                fSparseTracks[2] = pdg->Mass();
                fSparseTracks[3] = et;
                fSparseTracks[4] = part->Pt();
                fSparseTracks[5] = part->Eta();
                fSparseHistTracks->Fill(fSparseTracks);

                // Fill up total E_T counters for each particle species
                if (pdg->PdgCode() == fgProtonCode || pdg->PdgCode() == fgAntiProtonCode)
                {
                    fProtonEt += et;
                }
                if (pdg->PdgCode() == fgPiPlusCode || pdg->PdgCode() == fgPiMinusCode)
                {
                    fPionEt += et;
                    if (pdg->PdgCode() == fgPiPlusCode)
                    {
                        fPiPlusMult++;
                    }
                    else
                    {
                        fPiMinusMult++;
                    }
                }
                if (pdg->PdgCode() == fgGammaCode)
                {
                    fPiZeroMult++;
                }
                if (pdg->PdgCode() == fgKPlusCode || pdg->PdgCode() == fgKMinusCode)
                {
                    fChargedKaonEt += et;
                }
                if (pdg->PdgCode() == fgMuPlusCode || pdg->PdgCode() == fgMuMinusCode)
                {
                    fMuonEt += et;
                }
                if (pdg->PdgCode() == fgEPlusCode || pdg->PdgCode() == fgEMinusCode)
                {
                    fElectronEt += et;
                }

                // some neutrals also
                if (pdg->PdgCode() == fgNeutronCode)
                {
                    fNeutronEt += et;
                }
                if (pdg->PdgCode() == fgAntiNeutronCode)
                {
                    fAntiNeutronEt += et;
                }
                if (pdg->PdgCode() == fgGammaCode)
                {
                    fGammaEt += et;
                }

                // Neutral particles
                if (TMath::Abs(pdg->Charge() - fCuts->GetMonteCarloNeutralParticle()) <1e-3 )
                {

                    if (et > fCuts->GetCommonClusterEnergyCut()) fTotNeutralEt += et;

                    // inside EMCal acceptance
                    if (TMath::Abs(part->Eta()) < fEtaCutAcc && part->Phi() < fPhiCutAccMax && part->Phi() > fPhiCutAccMin)
                    {
                        fNeutralMultiplicity++;
                        //std::cout << pdg->PdgCode() << ", " << iPart << ", " <<  part->GetMother(0) << ", " << stack->IsPhysicalPrimary(iPart) << ", " << part->GetNDaughters() << ", " << part->Energy() << ", " << part->Eta() << ", " << part->Phi() << std::endl;

                        //if(part->GetDaughter(0) > 0) std::cout << stack->Particle(part->GetDaughter(0))->GetPdgCode() << " " << stack->Particle(part->GetDaughter(1))->GetPdgCode() << std::endl;
                        if (et > fCuts->GetCommonClusterEnergyCut()) fTotNeutralEtAcc += et;
                        if (et > fCuts->GetCommonClusterEnergyCut()) fTotEtAcc += et;
                        if (part->Energy() > 0.05) partCount++;
                    }
                }
                //Charged particles
                else if (TMath::Abs( pdg->Charge() - fCuts->GetMonteCarloNeutralParticle())>1e-3 )
                {
                    fChargedMultiplicity++;
                    fTotChargedEt += et;

                    // inside EMCal acceptance
                    if (TMath::Abs(part->Eta()) < fEtaCutAcc && part->Phi() < fPhiCutAccMax && part->Phi() > fPhiCutAccMin)
                    {
                        fTotChargedEtAcc += et;
                        fTotEtAcc += et;

                        if (pdg->PdgCode() == fgProtonCode || pdg->PdgCode() == fgAntiProtonCode)
                        {
                            fProtonEtAcc += et;
                        }
                        if (pdg->PdgCode() == fgPiPlusCode || pdg->PdgCode() == fgPiMinusCode)
                        {
                            fPionEtAcc += et;
                        }
                        if (pdg->PdgCode() == fgKPlusCode || pdg->PdgCode() == fgKMinusCode)
                        {
                            fChargedKaonEtAcc += et;
                        }
                        if (pdg->PdgCode() == fgMuPlusCode || pdg->PdgCode() == fgMuMinusCode)
                        {
                            fMuonEtAcc += et;
                        }
                        
                        if (pdg->PdgCode() == fgEPlusCode || pdg->PdgCode() == fgEMinusCode)
                        {
                          fElectronEtAcc += et;
                        }
                    } // inside EMCal acceptance

                    //	  if (TrackHitsCalorimeter(part, event->GetMagneticField()))
                    if (TrackHitsCalorimeter(part)) // magnetic field info not filled?
                    {
                        if (pdg->Charge() > 0) fHistPhivsPtPos->Fill(part->Phi(),part->Pt());
                        else fHistPhivsPtNeg->Fill(part->Phi(), part->Pt());
                    }
                }
            }

        }
    }

    fTotEt = fTotChargedEt + fTotNeutralEt;
    fTotEtAcc = fTotChargedEtAcc + fTotNeutralEtAcc;


    //std::cout << "Event done! # of particles: " << partCount << std::endl;
    return 0;
}
//Int_t AliAnalysisEtMonteCarlo::AnalyseEvent(AliMCEvent* mcEvent,AliESDEvent* realEvent)
Int_t AliAnalysisEtMonteCarlo::AnalyseEvent(AliVEvent* ev,AliVEvent* ev2)
{ // analyse MC and real event info
    //if(!mcEvent || !realEvent){
    if (!ev || !ev2) {
        AliFatal("ERROR: Event does not exist");
        return 0;
    }

    fSparseClusters[0] = 0;
    fSparseClusters[1] = 0;
    fSparseClusters[2] = 0;
    fSparseClusters[3] = 0;
    fSparseClusters[4] = 0;
    fSparseClusters[5] = 0;
    fSparseClusters[6] = 0;
    fSparseClusters[7] = 0;
    fSparseClusters[8] = 0;
    fSparseClusters[9] = fCentClass;
    fSparseClusters[10] = 0;

    //AnalyseEvent(mcEvent);
    AnalyseEvent(ev);
    AliMCEvent *mcEvent = dynamic_cast<AliMCEvent*>(ev);
    AliESDEvent *realEvent = dynamic_cast<AliESDEvent*>(ev2);
    if (!mcEvent || !realEvent) {
        AliFatal("ERROR: mcEvent or realEvent does not exist");
        return 0;
    }

    AliStack *stack = mcEvent->Stack();

    // get all detector clusters
    TRefArray* caloClusters = new TRefArray();
    if (fDetector == fCuts->GetDetectorEmcal()) realEvent->GetEMCALClusters( caloClusters );
    else if (fDetector == fCuts->GetDetectorPhos()) realEvent->GetPHOSClusters( caloClusters );
    else {
        AliFatal("Detector ID has not been specified");
        return -1;
    }

    // Track loop to fill a pT spectrum
    for (Int_t iTracks = 0; iTracks < realEvent->GetNumberOfTracks(); iTracks++)
    {
        AliESDtrack* track = realEvent->GetTrack(iTracks);
        if (!track)
        {
            Printf("ERROR: Could not receive track %d", iTracks);
            continue;
        }

        if ( track->Phi() < fCuts->GetGeometryPhosPhiAccMaxCut() * TMath::Pi()/180. && track->Phi() > fCuts->GetGeometryPhosPhiAccMinCut()  * TMath::Pi()/180. && TMath::Abs(track->Eta()) < fCuts->GetGeometryPhosEtaAccCut())
        {
            fTrackMultInAcc++;
        }
    }

    Int_t nCluster = caloClusters->GetEntries();

    // loop the clusters
    for (int iCluster = 0; iCluster < nCluster; iCluster++ )
    {
        AliESDCaloCluster* caloCluster = ( AliESDCaloCluster* )caloClusters->At( iCluster );
        //Float_t caloE = caloCluster->E();

        UInt_t iPart = (UInt_t)TMath::Abs(caloCluster->GetLabel());
        TParticle *part  = stack->Particle(iPart);
        TParticle *partMom = 0;

        if (!part)
        {
            Printf("No MC particle %d", iCluster);
            continue;
        }

        // compare MC and Rec energies for all particles
        //fHistAllERecEMC->Fill(part->Energy(),caloE);

        TParticlePDG *pdg = part->GetPDG(0);
        TParticlePDG *pdgMom = 0;

        if (!pdg)
        {
            Printf("ERROR: Could not get particle PDG %d", iPart);
            continue;
        }

        Int_t iPartMom = part->GetMother(0);

        if (iPartMom>0)
        {
            partMom = stack->Particle(iPartMom);
            pdgMom = partMom->GetPDG(0);
        }

        // Check if it is a primary particle
        //if (!stack->IsPhysicalPrimary(iPart)) continue;
        if (!stack->IsPhysicalPrimary(iPart)) // check whether particle is primary. we keep secondary electron and gamma for testing.
        {
            if (!((pdg->PdgCode() == fgEPlusCode) || (pdg->PdgCode() == fgEMinusCode) || (pdg->PdgCode() == fgGammaCode)))
                continue;
        } // end of primary particle check

        if (TMath::Abs(TMath::Abs(pdg->Charge()) - fCuts->GetMonteCarloSingleChargedParticle())<1e-3 && TMath::Abs(TMath::Abs(pdg->Charge()) - fCuts->GetMonteCarloNeutralParticle())<1e-3) continue;

        Double_t clEt = CalculateTransverseEnergy(caloCluster);
        if (caloCluster->E() < fCuts->GetCommonClusterEnergyCut()) continue;
        Float_t pos[3];

        caloCluster->GetPosition(pos);
        TVector3 cp(pos);

        fSparseClusters[0] = pdg->PdgCode();
        fSparseClusters[1] = pdg->Charge()/3;
        fSparseClusters[2] = pdg->Mass();
        fSparseClusters[3] = clEt;
        fSparseClusters[4] = caloCluster->E();
        fSparseClusters[5] = cp.Eta();
        fSparseClusters[6] = part->Energy() * TMath::Sin(part->Theta());
        fSparseClusters[7] = part->Pt();
        fSparseClusters[8] = part->Eta();
        fSparseClusters[9] = fCentClass;
        fSparseClusters[10] = 0;
        fSparseHistClusters->Fill(fSparseClusters);

        //if(caloCluster->GetEmcCpvDistance() < fTrackDistanceCut)

        if (TMath::Abs(caloCluster->GetTrackDx()) < fTrackDxCut && TMath::Abs(caloCluster->GetTrackDz()) < fTrackDzCut)
        {

            if (pdg->Charge() != 0)
            {
                //std::cout << "Removing charged: " << pdg->PdgCode() << ", with energy: " << caloCluster->E() << ", dist matched: " << caloCluster->GetTrackDx() << " " << caloCluster->GetTrackDz() << std::endl;
                fChargedRemoved++;
                fEnergyChargedRemoved += caloCluster->E();
                fHistRemovedOrNot->Fill(0.0, fCentClass);
                //std::cout << part->Vx() << " " << part->Vy() << " " << part->Vz() << " " << std::endl;
                //fHistDecayVertexRemovedCharged->Fill(part->Vy(), part->Vx(), part->Vz());
                fHistDxDzRemovedCharged->Fill(caloCluster->GetTrackDx(), caloCluster->GetTrackDz());


            }
            else
            {
                //std::cout << "Removing neutral: " << pdg->PdgCode() << ", with energy: " << caloCluster->E() << ", dist matched: " << caloCluster->GetTrackDx() << " " << caloCluster->GetTrackDz() << std::endl;
                //std::cout << "Mother is: " << stack->Particle(part->GetMother(0))->GetPdgCode() << std::endl;
                fNeutralRemoved++;
                fEnergyNeutralRemoved += caloCluster->E();
                fHistRemovedOrNot->Fill(1.0, fCentClass);
                //std::cout << part->Vx() << " " << part->Vy() << " " << part->Vz() << " " << std::endl;
                //fHistDecayVertexRemovedNeutral->Fill(part->Vy(), part->Vx(), part->Vz());
                fHistDxDzRemovedNeutral->Fill(caloCluster->GetTrackDx(), caloCluster->GetTrackDz());

                if (pdg->PdgCode() == fgGammaCode)
                {
                    fEtRemovedGammas += clEt;
                    fMultRemovedGammas++;
                }
                else if (pdg->PdgCode() == fgNeutronCode)
                {
                    fEtRemovedNeutrons += clEt;
                    fMultRemovedNeutrons++;
                }
                else if (pdg->PdgCode() == fgAntiNeutronCode)
                {
                    fEtRemovedAntiNeutrons += clEt;
                    fMultRemovedAntiNeutrons++;
                }
                //else std::cout << "Hmm, what is this (neutral: " << pdg->PdgCode() << std::endl;
            }
        }
        else
        {

            if (pdg->Charge() != 0)
            {
                //std::cout << "Not removing charged: " << pdg->PdgCode() << ", with energy: " << caloCluster->E() << ", dist matched: " << caloCluster->GetTrackDx() << " " << caloCluster->GetTrackDz() << std::endl;
                //std::cout << "Mother is: " << stack->Particle(part->GetMother(0))->GetPdgCode() << std::endl;
                fChargedNotRemoved++;
                fEnergyChargedNotRemoved += caloCluster->E();
                fHistRemovedOrNot->Fill(2.0, fCentClass);
                //std::cout << fHistDecayVertexNonRemovedCharged << std::endl;
                //std::cout << part->Vx() << " " << part->Vy() << " " << part->Vz() << " " << std::endl;
                //fHistDecayVertexNonRemovedCharged->Fill(part->Vy(), part->Vx(), part->Vz());
                fHistDxDzNonRemovedCharged->Fill(caloCluster->GetTrackDx(), caloCluster->GetTrackDz());
                if (pdg->PdgCode() == fgProtonCode)
                {
                    //std::cout << clEt << std::endl;
                    fEtNonRemovedProtons += clEt;
                    fMultNonRemovedProtons++;
                }
                else if (pdg->PdgCode() == fgAntiProtonCode)
                {
                    //std::cout << clEt << std::endl;
                    fEtNonRemovedAntiProtons += clEt;
                    fMultNonRemovedAntiProtons++;
                }
                else if (pdg->PdgCode() == fgPiPlusCode)
                {
                    //std::cout << "PI+" <<  clEt << std::endl;
                    fEtNonRemovedPiPlus += clEt;
                    fMultNonRemovedPiPlus++;
                }
                else if (pdg->PdgCode() == fgPiMinusCode)
                {
                    // std::cout << "PI-"  << clEt << std::endl;
                    fEtNonRemovedPiMinus += clEt;
                    fMultNonRemovedPiMinus++;
                }
                else if (pdg->PdgCode() == fgKPlusCode)
                {
                    //std::cout << clEt << std::endl;
                    fEtNonRemovedKaonPlus += clEt;
                    fMultNonRemovedKaonPlus++;
                }
                else if (pdg->PdgCode() == fgKMinusCode)
                {
                    //std::cout << clEt << std::endl;
                    fEtNonRemovedKaonMinus += clEt;
                    fMultNonRemovedKaonMinus++;
                }
                else if (pdg->PdgCode() == fgEPlusCode)
                {
                    //std::cout << clEt << std::endl;
                    if (TMath::Sqrt(part->Vx()*part->Vx() + part->Vy()*part->Vy()) < 430)
                    {
                        fEtNonRemovedPositrons += clEt;
                        fMultNonRemovedPositrons++;
                    }
                }
                else if (pdg->PdgCode() == fgEMinusCode)
                {
                    //std::cout << clEt << std::endl;
                    if (TMath::Sqrt(part->Vx()*part->Vx() + part->Vy()*part->Vy()) < 430)
                    {
                        fEtNonRemovedElectrons += clEt;
                        fMultNonRemovedElectrons++;
                    }
                }
                else if (pdg->PdgCode() == fgMuPlusCode)
                {
                    std::cout << clEt << std::endl;
                    fEtNonRemovedMuPlus += clEt;
                    fMultNonRemovedMuPlus++;
                }
                else if (pdg->PdgCode() == fgMuMinusCode)
                {
                    std::cout << clEt << std::endl;
                    fEtNonRemovedMuMinus += clEt;
                    fMultNonRemovedMuMinus++;
                }

            }
            else
            {
                //std::cout << "Accepted: " << pdg->PdgCode() << ", with energy: " << caloCluster->E() << ", dist matched: " << caloCluster->GetTrackDx() << " " << caloCluster->GetTrackDz() << std::endl;
                //std::cout << "Not removing charged: " << pdg->PdgCode() << ", with energy: " << caloCluster->E() << ", dist matched: " << caloCluster->GetTrackDx() << " " << caloCluster->GetTrackDz() << std::endl;
                //std::cout << "Mother is: " << stack->Particle(part->GetMother(0))->GetPdgCode() << stack->Particle(part->GetMother(0))->GetPDG()->GetName() << ", E: " << part->Energy() << std::endl;
                fNeutralNotRemoved++;
                fEnergyNeutralNotRemoved += caloCluster->E();
                fHistRemovedOrNot->Fill(3.0, fCentClass);
                //std::cout << part->Vx() << " " << part->Vy() << " " << part->Vz() << " " << std::endl;
                //fHistDecayVertexNonRemovedNeutral->Fill(part->Vy(), part->Vx(), part->Vz());
                fHistDxDzNonRemovedNeutral->Fill(caloCluster->GetTrackDx(), caloCluster->GetTrackDz());
                if (pdg->PdgCode() == fgGammaCode)
                {
                    fEtNonRemovedGammas += clEt;
                    fMultNonRemovedGammas++;
                    if (pdgMom)
                    {
                        if (TMath::Abs(pdgMom->PdgCode()) == fgPi0Code || TMath::Abs(pdgMom->PdgCode()) == fgEtaCode || TMath::Abs(pdgMom->PdgCode()) == 331)
                        {
//			std::cout << "Mother of gamma: " << pdgMom->PdgCode() << " " << pdgMom->GetName() <<  ", E: " << part->Energy() << std::endl;
                            fEtNonRemovedGammasFromPi0 += clEt;
                        }
                    }
                }
                else if (pdg->PdgCode() == fgNeutronCode)
                {
                    fEtNonRemovedNeutrons += clEt;
                    fMultNonRemovedNeutrons++;
                }
                else if (pdg->PdgCode() == fgAntiNeutronCode)
                {
                    fEtNonRemovedAntiNeutrons += clEt;
                    fMultNonRemovedAntiNeutrons++;
                }
                else if (pdg->PdgCode() == fgK0LCode || pdg->PdgCode() == fgK0SCode)
                {
                    fEtNonRemovedK0s += clEt;
                    fMultNonRemovedK0s++;

                }
                else if (TMath::Abs(pdg->PdgCode()) == fgLambdaCode)
                {
                    fEtNonRemovedLambdas += clEt;
                    fMultNonRemovedLambdas++;
                }
                else std::cout << "Hmm, what is this (neutral not removed): " << pdg->PdgCode() << " " << pdg->GetName() << ", ET: " << clEt << std::endl;
            }
        }
    } // end of loop over clusters

    delete caloClusters;
    //std::cout << "Distance cut: " << fTrackDistanceCut << std::endl;
    //std::cout << "Number of removed neutrals: " << fNeutralRemoved << std::endl;
    //std::cout << "Number of removed charged: " << fChargedRemoved << std::endl;
    //std::cout << "Number of non-removed charged: " << fChargedNotRemoved << std::endl;
    //std::cout << "Number of non-removed neutral: " << fNeutralNotRemoved << std::endl;

//  std::cout << "Energy of removed neutrals: " << fEnergyNeutralRemoved << std::endl;
//  std::cout << "Energy of removed charged: " << fEnergyChargedRemoved << std::endl;
//  std::cout << "Energy of non-removed charged: " << fEnergyChargedNotRemoved << std::endl;
//  std::cout << "Energy of non-removed neutral: " << fEnergyNeutralNotRemoved << std::endl;

    FillHistograms();
    return 0;
}

void AliAnalysisEtMonteCarlo::Init()
{ // init
    AliAnalysisEt::Init();
}

void AliAnalysisEtMonteCarlo::ResetEventValues()
{ // reset event values
    AliAnalysisEt::ResetEventValues();

    // collision geometry defaults for p+p:
    fImpactParameter = 0;
    fNcoll = 1;
    fNpart = 2;

    fEtNonRemovedProtons = 0;
    fEtNonRemovedAntiProtons = 0;
    fEtNonRemovedPiPlus = 0;
    fEtNonRemovedPiMinus = 0;
    fEtNonRemovedKaonPlus = 0;
    fEtNonRemovedKaonMinus = 0;
    fEtNonRemovedK0s = 0;
    fEtNonRemovedLambdas = 0;
    fEtNonRemovedElectrons = 0;
    fEtNonRemovedPositrons = 0;
    fEtNonRemovedMuPlus = 0;
    fEtNonRemovedMuMinus = 0;
    fEtNonRemovedNeutrons = 0;
    fEtNonRemovedAntiNeutrons = 0;
    fEtNonRemovedGammas = 0;
    fEtNonRemovedGammasFromPi0 = 0;

    fEtRemovedGammas = 0;
    fEtRemovedNeutrons = 0;
    fEtRemovedAntiNeutrons = 0;

    fMultNonRemovedProtons = 0;
    fMultNonRemovedAntiProtons = 0;
    fMultNonRemovedPiPlus = 0;
    fMultNonRemovedPiMinus = 0;
    fMultNonRemovedKaonPlus = 0;
    fMultNonRemovedKaonMinus = 0;
    fMultNonRemovedK0s = 0;
    fMultNonRemovedLambdas = 0;
    fMultNonRemovedElectrons = 0;
    fMultNonRemovedPositrons = 0;
    fMultNonRemovedMuPlus = 0;
    fMultNonRemovedMuMinus = 0;
    fMultNonRemovedNeutrons = 0;
    fMultNonRemovedAntiNeutrons = 0;
    fMultNonRemovedGammas = 0;

    fMultRemovedGammas = 0;
    fMultRemovedNeutrons = 0;
    fMultRemovedAntiNeutrons = 0;

    fTrackMultInAcc = 0;

}

void AliAnalysisEtMonteCarlo::CreateHistograms()
{ // histogram related additions
    AliAnalysisEt::CreateHistograms();
    if (fTree) {
        fTree->Branch("fImpactParameter",&fImpactParameter,"fImpactParameter/D");
        fTree->Branch("fNcoll",&fNcoll,"fNcoll/I");
        fTree->Branch("fNpart",&fNpart,"fNpart/I");
    }
	
    //fHistDecayVertexNonRemovedCharged = new TH3F("fHistDecayVertexNonRemovedCharged","fHistDecayVertexNonRemovedCharged", 500, -470, 30, 500, -300, 300, 40, -20, 20);
    //fHistDecayVertexRemovedCharged = new TH3F("fHistDecayVertexRemovedCharged","fHistDecayVertexRemovedCharged", 500, -470, 30, 500, -300, 300, 40, -20, 20);
    //fHistDecayVertexNonRemovedNeutral = new TH3F("fHistDecayVertexNonRemovedNeutral","fHistDecayVertexNonRemovedNeutral", 500, -470, 30, 500, -300, 300, 40, -20, 20);
    //fHistDecayVertexRemovedNeutral = new TH3F("fHistDecayVertexRemovedNeutral","fHistDecayVertexRemovedNeutral", 500, -470, 30, 500, -300, 300, 40, -20, 20);

    fHistRemovedOrNot = new TH2F("fHistRemovedOrNot", "fHistRemovedOrNot", 4, -0.5, 3.5, 10, -0.5, 9.5);

    fHistEtNonRemovedProtons = new TH2F("fHistEtNonRemovedProtons", "fHistEtNonRemovedProtons", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedAntiProtons = new TH2F("fHistEtNonRemovedAntiProtons", "fHistEtNonRemovedAntiProtons", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedPiPlus = new TH2F("fHistEtNonRemovedPiPlus", "fHistEtNonRemovedPiPlus", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedPiMinus = new TH2F("fHistEtNonRemovedPiMinus", "fHistEtNonRemovedPiMinus", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedKaonPlus = new TH2F("fHistEtNonRemovedKaonPlus", "fHistEtNonRemovedKaonPlus", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedKaonMinus = new TH2F("fHistEtNonRemovedKaonMinus", "fHistEtNonRemovedKaonMinus", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedK0s = new TH2F("fHistEtNonRemovedK0s", "fHistEtNonRemovedK0s", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedLambdas = new TH2F("fHistEtNonRemovedLambdas", "fHistEtNonRemovedLambdas", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedElectrons = new TH2F("fHistEtNonRemovedElectrons", "fHistEtNonRemovedElectrons", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedPositrons = new TH2F("fHistEtNonRemovedPositrons", "fHistEtNonRemovedPositrons", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedMuPlus = new TH2F("fHistEtNonRemovedMuPlus", "fHistEtNonRemovedMuPlus", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedMuMinus = new TH2F("fHistEtNonRemovedMuMinus", "fHistEtNonRemovedMuMinus", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedNeutrons = new TH2F("fHistEtNonRemovedNeutrons", "fHistEtNonRemovedNeutrons", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedAntiNeutrons = new TH2F("fHistEtNonRemovedAntiNeutrons", "fHistEtNonRemovedAntiNeutrons", 1500, 0, 30, 10, -0.5, 9.5);

    fHistEtNonRemovedGammas = new  TH2F("fHistEtNonRemovedGammas", "fHistEtNonRemovedGammas", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedGammasFromPi0 = new  TH2F("fHistEtNonRemovedGammasFromPi0", "fHistEtNonRemovedGammasFromPi0", 1500, 0, 30, 10, -0.5, 9.5);

    fHistEtRemovedGammas = new  TH2F("fHistEtRemovedGammas", "fHistEtRemovedGammas", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtRemovedNeutrons = new  TH2F("fHistEtRemovedNeutrons", "fHistEtRemovedNeutrons", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtRemovedAntiNeutrons = new  TH2F("fHistEtRemovedAntiNeutrons", "fHistEtRemovedAntiNeutrons", 1500, 0, 30, 10, -0.5, 9.5);

    fHistMultNonRemovedProtons = new TH2F("fHistMultNonRemovedProtons", "fHistMultNonRemovedProtons", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedAntiProtons = new TH2F("fHistMultNonRemovedAntiProtons", "fHistMultNonRemovedAntiProtons", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedPiPlus = new TH2F("fHistMultNonRemovedPiPlus", "fHistMultNonRemovedPiPlus", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedPiMinus = new TH2F("fHistMultNonRemovedPiMinus", "fHistMultNonRemovedPiMinus", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedKaonPlus = new TH2F("fHistMultNonRemovedKaonPlus", "fHistMultNonRemovedKaonPlus", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedKaonMinus = new TH2F("fHistMultNonRemovedKaonMinus", "fHistMultNonRemovedKaonMinus", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedK0s = new TH2F("fHistMultNonRemovedK0s", "fHistMultNonRemovedK0s", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedLambdas = new TH2F("fHistMultNonRemovedLambdas", "fHistMultNonRemovedLambdas", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedElectrons = new TH2F("fHistMultNonRemovedElectrons", "fHistMultNonRemovedElectrons", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedPositrons = new TH2F("fHistMultNonRemovedPositrons", "fHistMultNonRemovedPositrons", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedMuPlus = new TH2F("fHistMultNonRemovedMuPlus", "fHistMultNonRemovedMuPlus", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedMuMinus = new TH2F("fHistMultNonRemovedMuMinus", "fHistMultNonRemovedMuMinus", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedNeutrons = new TH2F("fHistMultNonRemovedNeutrons", "fHistMultNonRemovedNeutrons", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedAntiNeutrons = new TH2F("fHistMultNonRemovedAntiNeutrons", "fHistMultNonRemovedAntiNeutrons", 100, -0.5, 99.5, 10, -0.5, 9.5);

    fHistMultNonRemovedGammas = new  TH2F("fHistMultNonRemovedGammas", "fHistMultNonRemovedGammas", 100, -0.5, 99.5, 10, -0.5, 9.5);

    fHistMultRemovedGammas = new  TH2F("fHistMultRemovedGammas", "fHistMultRemovedGammas", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultRemovedNeutrons = new  TH2F("fHistMultRemovedNeutrons", "fHistMultRemovedNeutrons", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultRemovedAntiNeutrons = new  TH2F("fHistMultRemovedAntiNeutrons", "fHistMultRemovedAntiNeutrons", 100, -0.5, 99.5, 10, -0.5, 9.5);

    fHistTrackMultvsNonRemovedCharged = new TH2F("fHistTrackMultvsNonRemovedCharged", "fHistTrackMultvsNonRemovedCharged", 1000, -0.5, 999.5, 100, -0.5, 99.5);
    fHistTrackMultvsNonRemovedNeutral = new TH2F("fHistTrackMultvsNonRemovedNeutral", "fHistTrackMultvsNonRemovedNeutral", 1000, -0.5, 999.5, 100, -0.5, 99.5);
    fHistTrackMultvsRemovedGamma = new TH2F("fHistTrackMultvsRemovedGamma", "fHistTrackMultvsRemovedGamma", 1000, -0.5, 999.5, 100, -0.5, 99.5);

    fHistClusterMultvsNonRemovedCharged = new TH2F("fHistClusterMultvsNonRemovedCharged", "fHistClusterMultvsNonRemovedCharged", 1000, -0.5, 999.5, 100, -0.5, 99.5);
    fHistClusterMultvsNonRemovedNeutral = new TH2F("fHistClusterMultvsNonRemovedNeutral", "fHistClusterMultvsNonRemovedNeutral", 1000, -0.5, 999.5, 100, -0.5, 99.5);
    fHistClusterMultvsRemovedGamma = new TH2F("fHistClusterMultvsRemovedGamma", "fHistClusterMultvsRemovedGamma", 1000, -0.5, 999.5, 100, -0.5, 99.5);

    fHistDxDzNonRemovedCharged = new TH2F("fHistDxDzNonRemovedCharged", "fHistDxDzNonRemovedCharged", 800, -200, 200, 800, -200, 200);
    fHistDxDzRemovedCharged = new TH2F("fHistDxDzRemovedCharged", "fHistDxDzRemovedCharged", 800, -200, 200, 800, -200, 200);
    fHistDxDzNonRemovedNeutral = new TH2F("fHistDxDzNonRemovedNeutral", "fHistDxDzNonRemovedNeutral", 800, -200, 200, 800, -200, 200);
    fHistDxDzRemovedNeutral = new TH2F("fHistDxDzRemovedNeutral", "fHistDxDzRemovedNeutral", 800, -200, 200, 800, -200, 200);

    fHistPiPlusMult = new TH1F("fHistPiPlusMult", "fHistPiPlusMult", 2000, -0.5, 1999.5);
    fHistPiMinusMult = new TH1F("fHistPiMinusMult", "fHistPiMinusMult", 2000, -0.5, 1999.5);
    fHistPiZeroMult = new TH1F("fHistPiZeroMult", "fHistPiZeroMult", 2000, -0.5, 1999.5);

    fHistPiPlusMultAcc = new TH1F("fHistPiPlusMultAcc", "fHistPiPlusMultAcc", 2000, -0.5, 1999.5);
    fHistPiMinusMultAcc = new TH1F("fHistPiMinusMultAcc", "fHistPiMinusMultAcc", 2000, -0.5, 1999.5);
    fHistPiZeroMultAcc = new TH1F("fHistPiZeroMultAcc", "fHistPiZeroMultAcc", 2000, -0.5, 1999.5);

}

void AliAnalysisEtMonteCarlo::FillOutputList(TList *list)
{//fill the output list
    AliAnalysisEt::FillOutputList(list);

    list->Add(fHistRemovedOrNot);

    list->Add(fHistEtNonRemovedProtons);
    list->Add(fHistEtNonRemovedAntiProtons);
    list->Add(fHistEtNonRemovedPiPlus);
    list->Add(fHistEtNonRemovedPiMinus);
    list->Add(fHistEtNonRemovedKaonPlus);
    list->Add(fHistEtNonRemovedKaonMinus);
    list->Add(fHistEtNonRemovedK0s);
    list->Add(fHistEtNonRemovedLambdas);
    list->Add(fHistEtNonRemovedElectrons);
    list->Add(fHistEtNonRemovedPositrons);
    list->Add(fHistEtNonRemovedMuPlus);
    list->Add(fHistEtNonRemovedMuMinus);
    list->Add(fHistEtNonRemovedNeutrons);
    list->Add(fHistEtNonRemovedAntiNeutrons);
    list->Add(fHistEtNonRemovedGammas);
    list->Add(fHistEtNonRemovedGammasFromPi0);

    list->Add(fHistEtRemovedGammas);
    list->Add(fHistEtRemovedNeutrons);
    list->Add(fHistEtRemovedAntiNeutrons);


    list->Add(fHistMultNonRemovedProtons);
    list->Add(fHistMultNonRemovedAntiProtons);
    list->Add(fHistMultNonRemovedPiPlus);
    list->Add(fHistMultNonRemovedPiMinus);
    list->Add(fHistMultNonRemovedKaonPlus);
    list->Add(fHistMultNonRemovedKaonMinus);
    list->Add(fHistMultNonRemovedK0s);
    list->Add(fHistMultNonRemovedLambdas);
    list->Add(fHistMultNonRemovedElectrons);
    list->Add(fHistMultNonRemovedPositrons);
    list->Add(fHistMultNonRemovedMuPlus);
    list->Add(fHistMultNonRemovedMuMinus);
    list->Add(fHistMultNonRemovedNeutrons);
    list->Add(fHistMultNonRemovedAntiNeutrons);
    list->Add(fHistMultNonRemovedGammas);

    list->Add(fHistMultRemovedGammas);
    list->Add(fHistMultRemovedNeutrons);
    list->Add(fHistMultRemovedAntiNeutrons);

    list->Add(fHistTrackMultvsNonRemovedCharged);
    list->Add(fHistTrackMultvsNonRemovedNeutral);
    list->Add(fHistTrackMultvsRemovedGamma);

    list->Add(fHistClusterMultvsNonRemovedCharged);
    list->Add(fHistClusterMultvsNonRemovedNeutral);
    list->Add(fHistClusterMultvsRemovedGamma);

    //list->Add(fHistDecayVertexNonRemovedCharged);
    //list->Add(fHistDecayVertexNonRemovedNeutral);
    //list->Add(fHistDecayVertexRemovedCharged);
    //list->Add(fHistDecayVertexRemovedNeutral);

    list->Add(fHistDxDzNonRemovedCharged);
    list->Add(fHistDxDzRemovedCharged);
    list->Add(fHistDxDzNonRemovedNeutral);
    list->Add(fHistDxDzRemovedNeutral);

    list->Add(fHistPiPlusMult);
    list->Add(fHistPiMinusMult);
    list->Add(fHistPiZeroMult);
    list->Add(fHistPiPlusMultAcc);
    list->Add(fHistPiMinusMultAcc);
    list->Add(fHistPiZeroMultAcc);

}


bool AliAnalysisEtMonteCarlo::TrackHitsCalorimeter(TParticle* part, Double_t magField)
{
    //  printf(" TrackHitsCalorimeter - magField %f\n", magField);
    AliESDtrack *esdTrack = new AliESDtrack(part);
    // Printf("MC Propagating track: eta: %f, phi: %f, pt: %f", esdTrack->Eta(), esdTrack->Phi(), esdTrack->Pt());

    Bool_t prop = esdTrack->PropagateTo(fDetectorRadius, magField);

    // if(prop) Printf("Track propagated, eta: %f, phi: %f, pt: %f", esdTrack->Eta(), esdTrack->Phi(), esdTrack->Pt());

    bool status = prop &&
                  TMath::Abs(esdTrack->Eta()) < fEtaCutAcc &&
                  esdTrack->Phi() > fPhiCutAccMin*TMath::Pi()/180. &&
                  esdTrack->Phi() < fPhiCutAccMax*TMath::Pi()/180.;
    delete esdTrack;

    return status;
}

void AliAnalysisEtMonteCarlo::FillHistograms()
{ // let base class fill its histograms, and us fill the local ones
    AliAnalysisEt::FillHistograms();
    //std::cout << fEtNonRemovedPiPlus << " " << fCentClass << std::endl;

    fHistEtNonRemovedProtons->Fill(fEtNonRemovedProtons, fCentClass);
    fHistEtNonRemovedAntiProtons->Fill(fEtNonRemovedAntiProtons, fCentClass);
    fHistEtNonRemovedKaonPlus->Fill(fEtNonRemovedKaonPlus, fCentClass);
    fHistEtNonRemovedKaonMinus->Fill(fEtNonRemovedKaonMinus, fCentClass);
    fHistEtNonRemovedK0s->Fill(fEtNonRemovedK0s, fCentClass);
    fHistEtNonRemovedLambdas->Fill(fEtNonRemovedLambdas, fCentClass);
    fHistEtNonRemovedPiPlus->Fill(fEtNonRemovedPiPlus, fCentClass);
    fHistEtNonRemovedPiMinus->Fill(fEtNonRemovedPiMinus, fCentClass);
    fHistEtNonRemovedElectrons->Fill(fEtNonRemovedElectrons, fCentClass);
    fHistEtNonRemovedPositrons->Fill(fEtNonRemovedPositrons, fCentClass);
    fHistEtNonRemovedMuPlus->Fill(fEtNonRemovedMuPlus, fCentClass);
    fHistEtNonRemovedMuMinus->Fill(fEtNonRemovedMuMinus, fCentClass);
    fHistEtNonRemovedNeutrons->Fill(fEtNonRemovedNeutrons, fCentClass);
    fHistEtNonRemovedAntiNeutrons->Fill(fEtNonRemovedAntiNeutrons, fCentClass);
    fHistEtNonRemovedGammas->Fill(fEtNonRemovedGammas, fCentClass);
    fHistEtNonRemovedGammasFromPi0->Fill(fEtNonRemovedGammasFromPi0, fCentClass);

    fHistEtRemovedGammas->Fill(fEtRemovedGammas, fCentClass);
    fHistEtRemovedNeutrons->Fill(fEtRemovedNeutrons, fCentClass);
    fHistEtRemovedAntiNeutrons->Fill(fEtRemovedAntiNeutrons, fCentClass);

    fHistMultNonRemovedProtons->Fill(fMultNonRemovedProtons, fCentClass);
    fHistMultNonRemovedAntiProtons->Fill(fMultNonRemovedAntiProtons, fCentClass);
    fHistMultNonRemovedKaonPlus->Fill(fMultNonRemovedKaonPlus, fCentClass);
    fHistMultNonRemovedKaonMinus->Fill(fMultNonRemovedKaonMinus, fCentClass);
    fHistMultNonRemovedK0s->Fill(fMultNonRemovedK0s, fCentClass);
    fHistMultNonRemovedLambdas->Fill(fMultNonRemovedLambdas, fCentClass);
    fHistMultNonRemovedPiPlus->Fill(fMultNonRemovedPiPlus, fCentClass);
    fHistMultNonRemovedPiMinus->Fill(fMultNonRemovedPiMinus, fCentClass);
    fHistMultNonRemovedElectrons->Fill(fMultNonRemovedElectrons, fCentClass);
    fHistMultNonRemovedPositrons->Fill(fMultNonRemovedPositrons, fCentClass);
    fHistMultNonRemovedMuPlus->Fill(fMultNonRemovedMuPlus, fCentClass);
    fHistMultNonRemovedMuMinus->Fill(fMultNonRemovedMuMinus, fCentClass);
    fHistMultNonRemovedNeutrons->Fill(fMultNonRemovedNeutrons, fCentClass);
    fHistMultNonRemovedAntiNeutrons->Fill(fMultNonRemovedAntiNeutrons, fCentClass);
    fHistMultNonRemovedGammas->Fill(fMultNonRemovedGammas, fCentClass);

    fHistMultRemovedGammas->Fill(fMultRemovedGammas, fCentClass);
    fHistMultRemovedNeutrons->Fill(fMultRemovedNeutrons, fCentClass);
    fHistMultRemovedAntiNeutrons->Fill(fMultRemovedAntiNeutrons, fCentClass);

    fHistTrackMultvsNonRemovedCharged->Fill(fTrackMultInAcc,
                                            fMultNonRemovedAntiProtons+fMultNonRemovedElectrons+fMultNonRemovedKaonMinus+fMultNonRemovedKaonPlus
                                            +fMultNonRemovedMuMinus+fMultNonRemovedMuPlus+fMultNonRemovedPiMinus+fMultNonRemovedPiPlus+fMultNonRemovedPositrons
                                            +fMultNonRemovedProtons);

    fHistTrackMultvsNonRemovedNeutral->Fill(fTrackMultInAcc,
                                            fMultNonRemovedNeutrons+fMultNonRemovedAntiNeutrons+fMultNonRemovedK0s+fMultNonRemovedLambdas);

    fHistTrackMultvsRemovedGamma->Fill(fTrackMultInAcc,
                                       fMultRemovedGammas);

    fHistClusterMultvsNonRemovedCharged->Fill(fNeutralMultiplicity,
            fMultNonRemovedAntiProtons+fMultNonRemovedElectrons+fMultNonRemovedKaonMinus
            +fMultNonRemovedKaonPlus+fMultNonRemovedMuMinus+fMultNonRemovedMuPlus
            +fMultNonRemovedPiMinus+fMultNonRemovedPiPlus+fMultNonRemovedPositrons+fMultNonRemovedProtons);

    fHistClusterMultvsNonRemovedNeutral->Fill(fTrackMultInAcc,
            fMultNonRemovedNeutrons+fMultNonRemovedAntiNeutrons+fMultNonRemovedK0s+fMultNonRemovedLambdas);

    fHistClusterMultvsRemovedGamma->Fill(fTrackMultInAcc,
                                         fMultRemovedGammas);

}
















