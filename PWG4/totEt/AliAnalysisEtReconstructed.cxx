//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for ESD analysis
//  - reconstruction output
//  implementation file
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________

#include "AliAnalysisEtReconstructed.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDtrack.h"
#include "AliEMCALTrack.h"
#include "AliESDCaloCluster.h"
#include "TVector3.h"
#include "TGeoGlobalMagField.h"
#include "AliMagF.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliVParticle.h"
#include "TDatabasePDG.h"
#include "TList.h"
#include "AliESDpid.h"
#include <iostream>
#include "TH2F.h"
#include "AliAnalysisHadEtCorrections.h"
#include "AliLog.h"
#include "AliCentrality.h"




using namespace std;

ClassImp(AliAnalysisEtReconstructed);


AliAnalysisEtReconstructed::AliAnalysisEtReconstructed() :
        AliAnalysisEt()
        ,fCorrections(0)
        ,fPidCut(0)
        ,fClusterType(0)
        ,fHistChargedPionEnergyDeposit(0)
        ,fHistProtonEnergyDeposit(0)
        ,fHistAntiProtonEnergyDeposit(0)
        ,fHistChargedKaonEnergyDeposit(0)
        ,fHistMuonEnergyDeposit(0)
	,fHistRemovedEnergy(0)
        ,fGeomCorrection(1.0)
        ,fEMinCorrection(1.0)
{

}

AliAnalysisEtReconstructed::~AliAnalysisEtReconstructed()
{//destructor
    delete fCorrections;
    delete fHistChargedPionEnergyDeposit; /** Energy deposited in calorimeter by charged pions */    
    delete fHistProtonEnergyDeposit; /** Energy deposited in calorimeter by protons */    
    delete fHistAntiProtonEnergyDeposit; /** Energy deposited in calorimeter by anti-protons */    
    delete fHistChargedKaonEnergyDeposit; /** Energy deposited in calorimeter by charged kaons */    
    delete fHistMuonEnergyDeposit; /** Energy deposited in calorimeter by muons */

    delete fHistRemovedEnergy; // removed energy

}

Int_t AliAnalysisEtReconstructed::AnalyseEvent(AliVEvent* ev)
{
    // analyse ESD event
    ResetEventValues();

    if (!ev) {
        AliFatal("ERROR: Event does not exist");
        return 0;
    }

    AliESDEvent *event = dynamic_cast<AliESDEvent*>(ev);
    if (!event) {
        AliFatal("ERROR: ESD Event does not exist");
        return 0;
    }

    Int_t cent = -1;
    if (fCentrality)
    {
        cent = fCentrality->GetCentralityClass10("V0M");
	fCentClass = fCentrality->GetCentralityClass10("V0M");
    }

    Double_t protonMass = fgProtonMass;

    //for PID
    AliESDpid pID;
    pID.MakePID(event);
    TObjArray* list = fEsdtrackCutsTPC->GetAcceptedTracks(event);

    Int_t nGoodTracks = list->GetEntries();
    // printf("nGoodTracks %d nCaloClusters %d\n", nGoodTracks, event->GetNumberOfCaloClusters());

    for (Int_t iTrack = 0; iTrack < nGoodTracks; iTrack++)
    {
        AliESDtrack *track = dynamic_cast<AliESDtrack*> (list->At(iTrack));
        if (!track)
        {
            AliError(Form("ERROR: Could not get track %d", iTrack));
            continue;
        }

        fMultiplicity++;


        Float_t nSigmaPion,nSigmaProton,nSigmaKaon,nSigmaElectron;
        nSigmaPion = TMath::Abs(pID.NumberOfSigmasTPC(track,AliPID::kPion));
        nSigmaProton = TMath::Abs(pID.NumberOfSigmasTPC(track,AliPID::kProton));
        nSigmaKaon = TMath::Abs(pID.NumberOfSigmasTPC(track,AliPID::kKaon));
        nSigmaElectron = TMath::Abs(pID.NumberOfSigmasTPC(track,AliPID::kElectron));
        /*
        bool isPion = (nSigmaPion<3.0 && nSigmaProton>2.0 && nSigmaKaon>2.0);
        bool isElectron = (nSigmaElectron<2.0 && nSigmaPion>4.0 && nSigmaProton>3.0 && nSigmaKaon>3.0);
        bool isKaon = (nSigmaPion>3.0 && nSigmaProton>2.0 && nSigmaKaon<2.0);
        bool isProton = (nSigmaPion>3.0 && nSigmaProton<2.0 && nSigmaKaon>2.0);
        */

        Int_t nItsClusters = dynamic_cast<AliESDtrack*>(track)->GetNcls(0);
        Int_t nTPCClusters = dynamic_cast<AliESDtrack*>(track)->GetNcls(1);

        Float_t massPart = 0;

        const Double_t *pidWeights = track->PID();
        Int_t maxpid = -1;
        Double_t maxpidweight = 0;

        if (pidWeights)
        {
            for (Int_t p =0; p < AliPID::kSPECIES; p++)
            {
                if (pidWeights[p] > maxpidweight)
                {
                    maxpidweight = pidWeights[p];
                    maxpid = p;
                }
            }
            if (maxpid == AliPID::kProton)
            {
                //by definition of ET
                massPart = -protonMass*track->Charge();
            }

        }

        Double_t et = track->E() * TMath::Sin(track->Theta()) + massPart;

        fSparseTracks[0] = maxpid;
        fSparseTracks[1] = track->Charge();
        fSparseTracks[2] = track->M();
        fSparseTracks[3] = et;
        fSparseTracks[4] = track->Pt();
        fSparseTracks[5] = track->Eta();
        fSparseTracks[6] = cent;
        fSparseHistTracks->Fill(fSparseTracks);
        //printf("Rec track: iTrack %03d eta %4.3f phi %4.3f nITSCl %d nTPCCl %d\n", iTrack, track->Eta(), track->Phi(), nItsClusters, nTPCClusters); // tmp/debug printout

        if (TMath::Abs(track->Eta()) < fCuts->GetCommonEtaCut() && CheckGoodVertex(track) && nItsClusters > fCuts->GetReconstructedNItsClustersCut() && nTPCClusters > fCuts->GetReconstructedNTpcClustersCut() )
        {
            fTotChargedEt +=  et;
            fChargedMultiplicity++;
            if (maxpid != -1)
            {
                if (maxpid == AliPID::kProton)
                {
                    fProtonEt += et;
                }
                if (maxpid == AliPID::kPion)
                {
                    fPionEt += et;
                }
                if (maxpid == AliPID::kKaon)
                {
                    fChargedKaonEt += et;
                }
                if (maxpid == AliPID::kMuon)
                {
                    fMuonEt += et;
                }
                if (maxpid == AliPID::kElectron)
                {
                    fElectronEt += et;
                }
            }

            if (TMath::Abs(track->Eta()) < fEtaCutAcc && track->Phi() < fPhiCutAccMax && track->Phi() > fPhiCutAccMin)
            {
                fTotChargedEtAcc += track->E()*TMath::Sin(track->Theta()) + massPart;
                if (maxpid != -1)
                {
                    if (maxpid == AliPID::kProton)
                    {
                        fProtonEtAcc += et;
                    }
                    if (maxpid == AliPID::kPion)
                    {
                        fPionEtAcc += et;
                    }
                    if (maxpid == AliPID::kKaon)
                    {
                        fChargedKaonEtAcc += et;
                    }
                    if (maxpid == AliPID::kMuon)
                    {
                        fMuonEtAcc += et;
                    }
                    if (maxpid == AliPID::kElectron)
                    {
                        fElectronEtAcc += et;
                    }
                }

            }
        }

        if (TrackHitsCalorimeter(track, event->GetMagneticField()))
        {
            Double_t phi = track->Phi();
            Double_t pt = track->Pt();
            // printf("Rec track hit: iTrack %03d phi %4.3f pt %4.3f\n", iTrack, phi, pt); // tmp/debug printout
            if (track->Charge() > 0) fHistPhivsPtPos->Fill(phi, pt);
            else fHistPhivsPtNeg->Fill(phi, pt);
        }
    }

    for (Int_t iCluster = 0; iCluster < event->GetNumberOfCaloClusters(); iCluster++)
    {
        AliESDCaloCluster* cluster = event->GetCaloCluster(iCluster);
        if (!cluster)
        {
            AliError(Form("ERROR: Could not get cluster %d", iCluster));
            continue;
        }
        if (cluster->GetType() != fClusterType) continue;

        //if(cluster->GetTracksMatched() > 0)
//	printf("Rec Cluster: iCluster %03d E %4.3f type %.qd NCells %d, nmatched: %d, distance to closest: %f\n", iCluster, cluster->E(), (int)(cluster->GetType()), cluster->GetNCells(), cluster->GetNTracksMatched(), cluster->GetEmcCpvDistance()); // tmp/debug printout


        if (cluster->E() < fClusterEnergyCut) continue;

        Float_t pos[3];

        cluster->GetPosition(pos);
        TVector3 cp(pos);

        Double_t distance = cluster->GetEmcCpvDistance();
        Int_t trackMatchedIndex = cluster->GetTrackMatchedIndex();
        if ( cluster->IsEMCAL() ) {
            distance = CalcTrackClusterDistance(pos, &trackMatchedIndex, event);
        }

        fSparseClusters[0] = 0;
        fSparseClusters[1] = 0;
        fSparseClusters[2] = 0;
        fSparseClusters[6] = 0;
        fSparseClusters[7] = 0;
        fSparseClusters[8] = 0;
        fSparseClusters[9] = cent;
        fSparseClusters[10] = 0;


        if (cluster->GetNTracksMatched() > 0 && trackMatchedIndex > -1)
        {
            AliVTrack *tmptrack = event->GetTrack(trackMatchedIndex);
            if (!tmptrack)
            {
                AliError("Error: track does not exist");
                return -1;
            }
            const Double_t *pidWeights = tmptrack->PID();

            Double_t maxpidweight = 0;
            Int_t maxpid = 0;
            Double_t massPart = 0;
            if (pidWeights)
            {
                for (Int_t p =0; p < AliPID::kSPECIES; p++)
                {
                    if (pidWeights[p] > maxpidweight)
                    {
                        maxpidweight = pidWeights[p];
                        maxpid = p;
                    }
                }
                if (maxpid == AliPID::kProton)
                {
                    //by definition of ET
                    massPart = -protonMass*tmptrack->Charge();
                }
            }
            fSparseClusters[0] = maxpid;
            fSparseClusters[1] = tmptrack->Charge();
            fSparseClusters[2] = tmptrack->M();
            fSparseClusters[6] = tmptrack->E() * TMath::Sin(tmptrack->Theta()) + massPart;;
            fSparseClusters[7] = tmptrack->Pt();
            fSparseClusters[8] = tmptrack->Eta();
        }

        fSparseClusters[10] = distance;

        fHistTMDeltaR->Fill(distance);
        fHistTMDxDz->Fill(cluster->GetTrackDx(), cluster->GetTrackDz());

//        Float_t clusteret = cluster->E() * TMath::Sin(cp.Theta());

        Bool_t matched = false;

        if (cluster->IsEMCAL()) matched = distance < fTrackDistanceCut;
        else matched = (TMath::Abs(cluster->GetTrackDx()) < fTrackDxCut && TMath::Abs(cluster->GetTrackDz()) < fTrackDzCut);

        if (matched)
        {
            if (cluster->GetNTracksMatched() > 0 && trackMatchedIndex>=0)
            {
                AliVTrack *track = event->GetTrack(trackMatchedIndex);
                if (!track) {
                    AliError("Error: track does not exist");
                }
                else {
                    const Double_t *pidWeights = track->PID();

                    Double_t maxpidweight = 0;
                    Int_t maxpid = 0;

                    if (pidWeights)
                    {
                        for (Int_t p =0; p < AliPID::kSPECIES; p++)
                        {
                            if (pidWeights[p] > maxpidweight)
                            {
                                maxpidweight = pidWeights[p];
                                maxpid = p;
                            }
                        }
                        if (fCuts->GetHistMakeTreeDeposit() && fTreeDeposit)
                        {
                            fEnergyDeposited = cluster->E();
                            fEnergyTPC = track->E();
                            fCharge = track->Charge();
                            fParticlePid = maxpid;
                            fPidProb = maxpidweight;
                            AliESDtrack *esdTrack = dynamic_cast<AliESDtrack*>(track);
                            if (!esdTrack) {
                                AliError("Error: track does not exist");
                            }
                            else {
                                if (esdTrack) fTrackPassedCut = fEsdtrackCutsTPC->AcceptTrack(esdTrack);
                                fTreeDeposit->Fill();
                            }
                        }

                        if (maxpidweight > fPidCut)
                        {
                            Float_t dist = TMath::Sqrt(pos[0]*pos[0] + pos[1]*pos[1]);

                            Float_t theta = TMath::ATan(pos[2]/dist)+TMath::Pi()/2;

                            Float_t et = cluster->E() * TMath::Sin(theta);
                            if (maxpid == AliPID::kProton)
                            {

                                if (track->Charge() == 1)
                                {
                                    fBaryonEt += et;
                                    fHistProtonEnergyDeposit->Fill(cluster->E(), track->E());
                                }
                                else if (track->Charge() == -1)
                                {
                                    fAntiBaryonEt += et;
                                    fHistAntiProtonEnergyDeposit->Fill(cluster->E(), track->E());
                                }
                            }
                            else if (maxpid == AliPID::kPion)
                            {
                                fMesonEt += et;
                                fHistChargedPionEnergyDeposit->Fill(cluster->E(), track->E());
                            }
                            else if (maxpid == AliPID::kKaon)
                            {
                                fMesonEt += et;
                                fHistChargedKaonEnergyDeposit->Fill(cluster->E(), track->E());
                            }
                            else if (maxpid == AliPID::kMuon)
                            {
                                fHistMuonEnergyDeposit->Fill(cluster->E(), track->E());
                            }
                        }
                    }
                }
            }
            //continue;
        } // distance
        else
        {
            fSparseClusters[0] = AliPID::kPhoton;
            fSparseClusters[1] = 0;

            if (cluster->E() >  fSingleCellEnergyCut && cluster->GetNCells() == fCuts->GetCommonSingleCell()) continue;
            if (cluster->E() < fClusterEnergyCut) continue;
            cluster->GetPosition(pos);

            // TODO: replace with TVector3, too lazy now...

            float dist = TMath::Sqrt(pos[0]*pos[0] + pos[1]*pos[1]);

            float theta = TMath::ATan(pos[2]/dist)+TMath::Pi()/2;
            // float eta = TMath::Log(TMath::Abs( TMath::Tan( 0.5 * theta ) ) );
            fTotNeutralEt += cluster->E() * TMath::Sin(theta);
            fNeutralMultiplicity++;
        }

        cluster->GetPosition(pos);

        // TODO: replace with TVector3

        float dist = TMath::Sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
        float theta = TMath::ATan(pos[2]/dist)+TMath::Pi()/2;
        float eta = TMath::Log(TMath::Abs( TMath::Tan( 0.5 * theta ) ) );
        fSparseClusters[3] = cluster->E() * TMath::Sin(theta);
        fSparseClusters[4] = cluster->E();
        fSparseClusters[5] = eta;

        fSparseHistClusters->Fill(fSparseClusters);

        fMultiplicity++;
    }
    Double_t removedEnergy = GetChargedContribution(fNeutralMultiplicity) + GetNeutralContribution(fNeutralMultiplicity) - GetGammaContribution(fNeutralMultiplicity);
    fHistRemovedEnergy->Fill(removedEnergy);
//  std::cout << "fTotNeutralEt: " << fTotNeutralEt << ", Contribution from non-removed charged: " << GetChargedContribution(fNeutralMultiplicity) << ", neutral: " << GetNeutralContribution(fNeutralMultiplicity) << ", gammas: " << GetGammaContribution(fNeutralMultiplicity) << ", multiplicity: " << fNeutralMultiplicity<< std::endl;
    fTotNeutralEt = fGeomCorrection * fEMinCorrection * fTotNeutralEt - removedEnergy;
    fTotNeutralEtAcc = fTotNeutralEt/fGeomCorrection;
    fTotEt = fTotChargedEt + fTotNeutralEt;
    fTotEtAcc = fTotChargedEtAcc + fTotNeutralEtAcc;
    fSparseEt[0] = fTotEt;
    fSparseEt[1] = fTotNeutralEt;
    fSparseEt[2] = fTotChargedEtAcc;
    fSparseEt[3] = fMultiplicity;
    fSparseEt[4] = fNeutralMultiplicity;
    fSparseEt[5] = fChargedMultiplicity;
    fSparseEt[6] = cent;
    // Fill the histograms...
    FillHistograms();

    return 0;
}

bool AliAnalysisEtReconstructed::CheckGoodVertex(AliVParticle* track)
{ // check vertex

    Float_t bxy = 999.;
    Float_t bz = 999.;
    if (!track) {
        AliError("ERROR: no track");
        return kFALSE;
    }
    AliESDtrack *esdTrack = dynamic_cast<AliESDtrack*>(track);
    if (!esdTrack) {
        AliError("ERROR: no track");
        return kFALSE;
    }
    esdTrack->GetImpactParametersTPC(bxy,bz);


    bool status = (TMath::Abs(track->Xv()) < fCuts->GetReconstructedVertexXCut()) &&
                  (TMath::Abs(track->Yv()) < fCuts->GetReconstructedVertexYCut()) &&
                  (TMath::Abs(track->Zv()) < fCuts->GetReconstructedVertexZCut()) &&
                  (TMath::Abs(bxy) < fCuts->GetReconstructedIPxyCut()) &&
                  (TMath::Abs(bz) < fCuts->GetReconstructedIPzCut());

    return status;
}

void AliAnalysisEtReconstructed::Init()
{ // Init
    AliAnalysisEt::Init();
    fPidCut = fCuts->GetReconstructedPidCut();
    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 1., 1., AliMagF::k5kG));
    if (!fCorrections) {
        cout<<"Warning!  You have not set corrections.  Your code will crash.  You have to set the corrections."<<endl;
    }
}

bool AliAnalysisEtReconstructed::TrackHitsCalorimeter(AliVParticle* track, Double_t magField)
{ // propagate track to detector radius

    if (!track) {
        cout<<"Warning: track empty"<<endl;
        return kFALSE;
    }
    AliESDtrack *esdTrack= dynamic_cast<AliESDtrack*>(track);
    if (!esdTrack) {
        AliError("ERROR: no ESD track");
        return kFALSE;
    }
    // Printf("Propagating track: eta: %f, phi: %f, pt: %f", esdTrack->Eta(), esdTrack->Phi(), esdTrack->Pt());

    Bool_t prop = esdTrack->PropagateTo(fDetectorRadius, magField);

    // if (prop) Printf("Track propagated, eta: %f, phi: %f, pt: %f", esdTrack->Eta(), esdTrack->Phi(), esdTrack->Pt());
    return prop &&
           TMath::Abs(esdTrack->Eta()) < fEtaCutAcc &&
           esdTrack->Phi() > fPhiCutAccMin*TMath::Pi()/180. &&
           esdTrack->Phi() < fPhiCutAccMax*TMath::Pi()/180.;
}

void AliAnalysisEtReconstructed::FillOutputList(TList* list)
{ // add some extra histograms to the ones from base class
    AliAnalysisEt::FillOutputList(list);

    list->Add(fHistChargedPionEnergyDeposit);
    list->Add(fHistProtonEnergyDeposit);
    list->Add(fHistAntiProtonEnergyDeposit);
    list->Add(fHistChargedKaonEnergyDeposit);
    list->Add(fHistMuonEnergyDeposit);

    list->Add(fHistRemovedEnergy);
}

void AliAnalysisEtReconstructed::CreateHistograms()
{ // add some extra histograms to the ones from base class
    AliAnalysisEt::CreateHistograms();

    Int_t nbinsEt = 1000;
    Double_t minEt = 0;
    Double_t maxEt = 10;

    // possibly change histogram limits
    if (fCuts) {
        nbinsEt = fCuts->GetHistNbinsParticleEt();
        minEt = fCuts->GetHistMinParticleEt();
        maxEt = fCuts->GetHistMaxParticleEt();
    }

    TString histname;
    histname = "fHistChargedPionEnergyDeposit" + fHistogramNameSuffix;
    fHistChargedPionEnergyDeposit = new TH2F(histname.Data(), "Energy deposited by #pi^{+/-}", nbinsEt, minEt, maxEt, nbinsEt, minEt, maxEt);
    fHistChargedPionEnergyDeposit->SetXTitle("Energy deposited in calorimeter");
    fHistChargedPionEnergyDeposit->SetYTitle("Energy of track");

    histname = "fHistProtonEnergyDeposit" + fHistogramNameSuffix;
    fHistProtonEnergyDeposit = new TH2F(histname.Data(), "Energy deposited by protons", nbinsEt, minEt, maxEt, nbinsEt, minEt, maxEt);
    fHistProtonEnergyDeposit->SetXTitle("Energy deposited in calorimeter");
    fHistProtonEnergyDeposit->SetYTitle("Energy of track");

    histname = "fHistAntiProtonEnergyDeposit" + fHistogramNameSuffix;
    fHistAntiProtonEnergyDeposit = new TH2F(histname.Data(), "Energy deposited by anti-protons", nbinsEt, minEt, maxEt, nbinsEt, minEt, maxEt);
    fHistAntiProtonEnergyDeposit->SetXTitle("Energy deposited in calorimeter");
    fHistAntiProtonEnergyDeposit->SetYTitle("Energy of track");

    histname = "fHistChargedKaonEnergyDeposit" + fHistogramNameSuffix;
    fHistChargedKaonEnergyDeposit = new TH2F(histname.Data(), "Energy deposited by K^{+/-}", nbinsEt, minEt, maxEt, nbinsEt, minEt, maxEt);
    fHistChargedKaonEnergyDeposit->SetXTitle("Energy deposited in calorimeter");
    fHistChargedKaonEnergyDeposit->SetYTitle("Energy of track");

    histname = "fHistMuonEnergyDeposit" + fHistogramNameSuffix;
    fHistMuonEnergyDeposit = new TH2F(histname.Data(), "Energy deposited by #mu^{+/-}", nbinsEt, minEt, maxEt, nbinsEt, minEt, maxEt);
    fHistMuonEnergyDeposit->SetXTitle("Energy deposited in calorimeter");
    fHistMuonEnergyDeposit->SetYTitle("Energy of track");

    histname = "fHistRemovedEnergy" + fHistogramNameSuffix;
    fHistRemovedEnergy = new TH1F(histname.Data(), histname.Data(), 1000, 0, 20);
    //fHistMuonEnergyDeposit->SetXTitle("Energy deposited in calorimeter");
    //fHistMuonEnergyDeposit->SetYTitle("Energy of track");



}

Double_t
AliAnalysisEtReconstructed::CalcTrackClusterDistance(const Float_t clsPos[3],
        Int_t *trkMatchId,
        const AliESDEvent *event)
{ // calculate distance between cluster and closest track

    Double_t trkPos[3] = {0,0,0};

    Int_t bestTrkMatchId = -1;
    Double_t distance = 9999; // init to a big number

    Double_t dist = 0;
    Double_t distX = 0, distY = 0, distZ = 0;

    for (Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++) {
        AliESDtrack *track = event->GetTrack(iTrack);
        if (!track) {
            AliError(Form("ERROR: Could not get track %d", iTrack));
            continue;
        }

        // check for approx. eta and phi range before we propagate..
        // TBD

        AliEMCALTrack *emctrack = new AliEMCALTrack(*track);
        if (!emctrack->PropagateToGlobal(clsPos[0],clsPos[1],clsPos[2],0.,0.) ) {
            continue;
        }
        emctrack->GetXYZ(trkPos);
        delete emctrack;

        distX = clsPos[0]-trkPos[0];
        distY = clsPos[1]-trkPos[1];
        distZ = clsPos[2]-trkPos[2];
        dist = TMath::Sqrt(distX*distX + distY*distY + distZ*distZ);

        if (dist < distance) {
            distance = dist;
            bestTrkMatchId = iTrack;
        }
    } // iTrack

    // printf("CalcTrackClusterDistance: bestTrkMatch %d origTrkMatch %d distance %f\n", bestTrkMatchId, *trkMatchId, distance);
    *trkMatchId = bestTrkMatchId;
    return distance;
}

