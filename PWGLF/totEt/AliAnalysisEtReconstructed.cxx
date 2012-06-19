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
#include "TH2I.h"
#include "TH1I.h"
#include "TFile.h"
#include "AliAnalysisHadEtCorrections.h"
#include "AliAnalysisEtSelector.h"
#include "AliLog.h"
#include "AliCentrality.h"
#include "AliPHOSGeoUtils.h"
#include "AliPHOSGeometry.h"


using namespace std;

ClassImp(AliAnalysisEtReconstructed);


AliAnalysisEtReconstructed::AliAnalysisEtReconstructed() :
        AliAnalysisEt()
        ,fCorrections(0)
        ,fPidCut(0)
        ,fHistChargedPionEnergyDeposit(0)
        ,fHistProtonEnergyDeposit(0)
        ,fHistAntiProtonEnergyDeposit(0)
        ,fHistChargedKaonEnergyDeposit(0)
        ,fHistMuonEnergyDeposit(0)
        ,fHistRemovedEnergy(0)
        ,fGeomCorrection(1.0)
        ,fEMinCorrection(1.0/0.89)
	,fRecEffCorrection(1.0)
	,fClusterPosition(0)
	,fHistChargedEnergyRemoved(0)
	,fHistNeutralEnergyRemoved(0)
	,fHistGammaEnergyAdded(0)
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
    delete fClusterPosition;
    delete fHistChargedEnergyRemoved;
    delete fHistNeutralEnergyRemoved;
    delete fHistGammaEnergyAdded;

}

Int_t AliAnalysisEtReconstructed::AnalyseEvent(AliVEvent* ev)
{

    //AliAnalysisEt::AnalyseEvent(ev);
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
    fSelector->SetEvent(event);
    
    Int_t cent = -1;
    if (fCentrality)
    {
        cent = fCentrality->GetCentralityClass10("V0M");
        fCentClass = fCentrality->GetCentralityClass10("V0M");
    }

    for (Int_t iCluster = 0; iCluster < event->GetNumberOfCaloClusters(); iCluster++)
    {
        AliESDCaloCluster* cluster = event->GetCaloCluster(iCluster);
        if (!cluster)
        {
            AliError(Form("ERROR: Could not get cluster %d", iCluster));
            continue;
        }
        int x = 0;
	fCutFlow->Fill(x++);
	if(cluster->IsEMCAL()) continue;
	fCutFlow->Fill(x++);
	if(!fSelector->CutMinEnergy(*cluster)) continue;
	fCutFlow->Fill(x++);
        if (!fSelector->CutDistanceToBadChannel(*cluster)) continue;
	fCutFlow->Fill(x++);

        Float_t pos[3];

        cluster->GetPosition(pos);
        TVector3 cp(pos);
	
	Double_t distance = cluster->GetEmcCpvDistance();
        Int_t trackMatchedIndex = cluster->GetTrackMatchedIndex();
        if ( cluster->IsEMCAL() ) {
            distance = CalcTrackClusterDistance(pos, &trackMatchedIndex, event);
        }

        Bool_t matched = false;

	
	matched = !fSelector->CutTrackMatching(*cluster);

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
                        if (fCuts->GetHistMakeTreeDeposit() && fDepositTree)
                        {
                            fEnergyDeposited = cluster->E();
                            fMomentumTPC = track->P();
                            fCharge = track->Charge();
                            fParticlePid = maxpid;
                            fPidProb = maxpidweight;
                            AliESDtrack *esdTrack = dynamic_cast<AliESDtrack*>(track);
                            if (!esdTrack) {
                                AliError("Error: track does not exist");
                            }
                            else {
                                if (esdTrack) fTrackPassedCut = fEsdtrackCutsTPC->AcceptTrack(esdTrack);
                                fDepositTree->Fill();
                            }
                        }

                        if (maxpidweight > fPidCut)
                        {
                            //Float_t dist = TMath::Sqrt(pos[0]*pos[0] + pos[1]*pos[1]);

                            //Float_t theta = TMath::ATan(pos[2]/dist)+TMath::Pi()/2;

                            //Float_t et = cluster->E() * TMath::Sin(theta);
                            if (maxpid == AliPID::kProton)
                            {

                                if (track->Charge() == 1)
                                {
                                    fHistProtonEnergyDeposit->Fill(cluster->E(), track->E());
                                }
                                else if (track->Charge() == -1)
                                {
                                    fHistAntiProtonEnergyDeposit->Fill(cluster->E(), track->E());
                                }
                            }
                            else if (maxpid == AliPID::kPion)
                            {
                                fHistChargedPionEnergyDeposit->Fill(cluster->E(), track->E());
                            }
                            else if (maxpid == AliPID::kKaon)
                            {
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
	  fCutFlow->Fill(x++);
	  //std::cout << x++ << std::endl;

            //if (cluster->E() >  fSingleCellEnergyCut && cluster->GetNCells() == fCuts->GetCommonSingleCell()) continue;
            //if (cluster->E() < fClusterEnergyCut) continue;
            cluster->GetPosition(pos);
	    
	    TVector3 p2(pos);
	    
	    fClusterPosition->Fill(p2.Phi(), p2.PseudoRapidity());

	    fTotNeutralEt += CalculateTransverseEnergy(cluster);
            fNeutralMultiplicity++;
        }
        fMultiplicity++;
    }
    
    fChargedEnergyRemoved = GetChargedContribution(fNeutralMultiplicity);
    fNeutralEnergyRemoved = GetNeutralContribution(fNeutralMultiplicity);
    fHistChargedEnergyRemoved->Fill(fChargedEnergyRemoved, fNeutralMultiplicity);
    fHistNeutralEnergyRemoved->Fill(fNeutralEnergyRemoved, fNeutralMultiplicity);
    
    fGammaEnergyAdded = GetGammaContribution(fNeutralMultiplicity);
    fHistGammaEnergyAdded->Fill(fGammaEnergyAdded, fNeutralMultiplicity);

    Double_t removedEnergy = GetChargedContribution(fNeutralMultiplicity) + GetNeutralContribution(fNeutralMultiplicity) - GetGammaContribution(fNeutralMultiplicity) + GetSecondaryContribution(fNeutralMultiplicity);
    fHistRemovedEnergy->Fill(removedEnergy);
    
    fTotNeutralEt = fGeomCorrection * fEMinCorrection * (fTotNeutralEt - removedEnergy);
    fTotEt = fTotChargedEt + fTotNeutralEt;
// Fill the histograms...0
    FillHistograms();
   // std::cout << "fTotNeutralEt: " << fTotNeutralEt << ", Contribution from non-removed charged: " << GetChargedContribution(fNeutralMultiplicity) << ", neutral: " << GetNeutralContribution(fNeutralMultiplicity) << ", gammas: " << GetGammaContribution(fNeutralMultiplicity) << ", multiplicity: " << fNeutralMultiplicity<< std::endl;
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
    return prop && fSelector->CutGeometricalAcceptance(*esdTrack);
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
    list->Add(fClusterPosition);
    
    list->Add(fHistChargedEnergyRemoved);
    list->Add(fHistNeutralEnergyRemoved);
    list->Add(fHistGammaEnergyAdded);
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

    histname = "fClusterPosition" + fHistogramNameSuffix;
    fClusterPosition = new TH2D(histname.Data(), "Position of accepted neutral clusters",1000, -2.0, -.5, 1000, -.13 , 0.13);
    fClusterPosition->SetXTitle("Energy deposited in calorimeter");
    fClusterPosition->SetYTitle("Energy of track");


    histname = "fHistChargedEnergyRemoved" + fHistogramNameSuffix;
    fHistChargedEnergyRemoved = new TH2D(histname.Data(), histname.Data(), 1000, .0, 30, 100, -0.5 , 99.5);

    histname = "fHistNeutralEnergyRemoved" + fHistogramNameSuffix;
    fHistNeutralEnergyRemoved = new TH2D(histname.Data(), histname.Data(), 1000, .0, 30, 100, -0.5 , 99.5);

    histname = "fHistGammaEnergyAdded" + fHistogramNameSuffix;
    fHistGammaEnergyAdded = new TH2D(histname.Data(), histname.Data(), 1000, .0, 30, 100, -0.5 , 99.5);


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
