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
    ,fGeoUtils(0)
    ,fBadMapM2(0)
    ,fBadMapM3(0)
    ,fBadMapM4(0)
    ,fClusterPosition(0)
    ,fHistChargedEnergyRemoved(0)
    ,fHistNeutralEnergyRemoved(0)
    ,fHistGammaEnergyAdded(0)
{

}

AliAnalysisEtReconstructed::~AliAnalysisEtReconstructed()
{   //destructor
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
    fSelector->SetEvent(event);
    if (!fMatrixInitialized)
    {
        for (Int_t mod=0; mod<5; mod++) {
            if (!event->GetPHOSMatrix(mod)) continue;
            fGeoUtils->SetMisalMatrix(event->GetPHOSMatrix(mod),mod) ;
//	    std::cout << event->GetPHOSMatrix(mod) << std::endl;
            Printf("PHOS geo matrix %p for module # %d is set\n", event->GetPHOSMatrix(mod), mod);
        }
        fMatrixInitialized = kTRUE;
    }

    Int_t cent = -1;
    if (fCentrality)
    {
        cent = fCentrality->GetCentralityClass10("V0M");
        fCentClass = fCentrality->GetCentralityClass10("V0M");
    }

    //Double_t protonMass = fgProtonMass;

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

        Double_t r = 0;

        matched = !fSelector->CutTrackMatching(*cluster, r);

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
            fCutFlow->Fill(x++);
            //std::cout << x++ << std::endl;
            fSparseClusters[0] = AliPID::kPhoton;
            fSparseClusters[1] = 0;

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

    Double_t removedEnergy = GetChargedContribution(fNeutralMultiplicity) + GetNeutralContribution(fNeutralMultiplicity) - GetGammaContribution(fNeutralMultiplicity);
    fHistRemovedEnergy->Fill(removedEnergy);

    fTotNeutralEt = fGeomCorrection * fEMinCorrection * (fTotNeutralEt - removedEnergy);
    fTotNeutralEtAcc = fTotNeutralEt;
    fTotEt = fTotChargedEt + fTotNeutralEt;
    fTotEtAcc = fTotChargedEtAcc + fTotNeutralEtAcc;
    if(fMakeSparse) {
        fSparseEt[0] = fTotEt;
        fSparseEt[1] = fTotNeutralEt;
        fSparseEt[2] = fTotChargedEtAcc;
        fSparseEt[3] = fMultiplicity;
        fSparseEt[4] = fNeutralMultiplicity;
        fSparseEt[5] = fChargedMultiplicity;
        fSparseEt[6] = cent;
    }
    // Fill the histograms...
    FillHistograms();

    return 0;
}

bool AliAnalysisEtReconstructed::CheckGoodVertex(AliVParticle* track)
{   // check vertex

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
{   // Init
    AliAnalysisEt::Init();
    fPidCut = fCuts->GetReconstructedPidCut();
    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 1., 1., AliMagF::k5kG));
    if (!fCorrections) {
        cout<<"Warning!  You have not set corrections.  Your code will crash.  You have to set the corrections."<<endl;
    }
    //fGeoUtils = new AliPHOSGeoUtils("PHOS", "noCPV");
    fGeoUtils = AliPHOSGeometry::GetInstance("IHEP");
    // ifstream f("badchannels.txt", ios::in);
    TFile *f = TFile::Open("badchannels.root", "READ");

    fBadMapM2 = (TH2I*)f->Get("bad_channels_m2");
    fBadMapM3 = (TH2I*)f->Get("bad_channels_m3");
    fBadMapM4 = (TH2I*)f->Get("bad_channels_m4");
}

bool AliAnalysisEtReconstructed::TrackHitsCalorimeter(AliVParticle* track, Double_t magField)
{   // propagate track to detector radius

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
{   // add some extra histograms to the ones from base class
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
{   // add some extra histograms to the ones from base class
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
{   // calculate distance between cluster and closest track

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


Bool_t AliAnalysisEtReconstructed::TooCloseToBadChannel(const AliESDCaloCluster &cluster) const
{

    Float_t gPos[3];
    cluster.GetPosition(gPos);
    Int_t relId[4];
    TVector3 glVec(gPos);
    fGeoUtils->GlobalPos2RelId(glVec, relId);

    TVector3 locVec;
    fGeoUtils->Global2Local(locVec, glVec, relId[0]);
//    std::cout << fGeoUtils << std::endl;
    //std::cout << relId[0] << " " << cluster.IsPHOS() <<  std::endl;
    //std::cout << locVec[0] << " " << " " << locVec[1] << " " << locVec[2] << std::endl;
    for (Int_t x = 0; x < fBadMapM2->GetNbinsX(); x++)
    {
        for (Int_t z = 0; z < fBadMapM2->GetNbinsY(); z++)
        {
            if (relId[0] == 3)
            {
                if (fBadMapM2->GetBinContent(x+1, z+1) != 0)
                {
                    Int_t tmpRel[4];
                    tmpRel[0] = 3;
                    tmpRel[1] = 0;
                    tmpRel[2] = x+1;
                    tmpRel[3] = z+1;

                    Float_t tmpX;
                    Float_t tmpZ;
                    fGeoUtils->RelPosInModule(tmpRel, tmpX, tmpZ);

                    Float_t distance = TMath::Sqrt((tmpX-locVec[0])*(tmpX-locVec[0]) + (tmpZ - locVec[2])*(tmpZ-locVec[2]));
                    //Float_t distance = TMath::Sqrt((x-relId[3])*(x-relId[3]) + (z - relId[2])*(z-relId[2]));

                    if (distance < fCuts->GetPhosBadDistanceCut())
                    {
//		      std::cout << "Module 2, position: " << locVec[0] << ", " << locVec[2] << ", distance to bad channel: " << distance << ", number of cells: " << cluster.GetNCells() <<  std::endl;
                        return kTRUE;
                    }
                }
            }
            if (relId[0] == 2)
            {
                if (fBadMapM3->GetBinContent(x+1, z+1) != 0)
                {
                    Int_t tmpRel[4];
                    tmpRel[0] = 2;
                    tmpRel[1] = 0;
                    tmpRel[2] = x+1;
                    tmpRel[3] = z+1;

                    Float_t tmpX;
                    Float_t tmpZ;
                    fGeoUtils->RelPosInModule(tmpRel, tmpX, tmpZ);

                    Float_t distance = TMath::Sqrt((tmpX-locVec[0])*(tmpX-locVec[0]) + (tmpZ - locVec[2])*(tmpZ-locVec[2]));

//                    Float_t distance = TMath::Sqrt((x-locVec[0])*(x-locVec[0]) + (z - locVec[2])*(z-locVec[2]));
                    //Float_t distance = TMath::Sqrt((x-relId[3])*(x-relId[3]) + (z - relId[2])*(z-relId[2]));
                    if (distance < fCuts->GetPhosBadDistanceCut())
                    {
//		      std::cout << "Module 3, position: " << locVec[0] << ", " << locVec[2] << ", distance to bad channel: " << distance << ", number of cells: " << cluster.GetNCells() <<  std::endl;
                        return kTRUE;
                    }
                }
            }
            if (relId[0] == 1)
            {
                if (fBadMapM4->GetBinContent(x+1, z+1) != 0)
                {
                    Int_t tmpRel[4];
                    tmpRel[0] = 1;
                    tmpRel[1] = 0;
                    tmpRel[2] = x+1;
                    tmpRel[3] = z+1;

                    Float_t tmpX;
                    Float_t tmpZ;
                    fGeoUtils->RelPosInModule(tmpRel, tmpX, tmpZ);

                    Float_t distance = TMath::Sqrt((tmpX-locVec[0])*(tmpX-locVec[0]) + (tmpZ - locVec[2])*(tmpZ-locVec[2]));

//                    Float_t distance = TMath::Sqrt((x-locVec[0])*(x-locVec[0]) + (z - locVec[2])*(z-locVec[2]));
                    //Float_t distance = TMath::Sqrt((x-relId[3])*(x-relId[3]) + (z - relId[2])*(z-relId[2]));
                    if (distance < fCuts->GetPhosBadDistanceCut())
                    {
//			std::cout << "Module 4, position: " << locVec[0] << ", " << locVec[2] << ", distance to bad channel: " << distance << ", number of cells: " << cluster.GetNCells() <<  std::endl;
                        return kTRUE;
                    }
                }
            }

        }
    }

    return kFALSE;


}




/*
Bool_t AliAnalysisEtReconstructed::TooCloseToBadChannel(const AliESDCaloCluster &cluster) const
{

    Float_t gPos[3];

    cluster.GetPosition(gPos);
    Int_t relId[4];
    TVector3 glVec(gPos);
    fGeoUtils->GlobalPos2RelId(glVec, relId);
    TVector3 locVec;
    fGeoUtils->Global2Local(locVec, glVec, relId[0]);

    std::vector<Int_t>::const_iterator badIt;

    for (Int_t x = 0; x < fBadMapM2->GetNbinsX(); x++)
    {
        for (Int_t z = 0; z < fBadMapM2->GetNbinsY(); z++)
        {
            if (relId[0] == 3)
            {
                if (fBadMapM2->GetBinContent(x+1, z+1) != 0)
                {

                    Float_t distance = TMath::Sqrt((x-locVec[0])*(x-locVec[0]) + (z - locVec[2])*(z-locVec[2]));
                    if (distance < fCuts->GetPhosBadDistanceCut())
                    {
                        return kTRUE;
                    }
                }
            }
            if (relId[0] == 2)
            {
                if (fBadMapM3->GetBinContent(x+1, z+1) != 0)
                {

                    Float_t distance = TMath::Sqrt((x-locVec[0])*(x-locVec[0]) + (z - locVec[2])*(z-locVec[2]));
                    if (distance < fCuts->GetPhosBadDistanceCut())
                    {
                        return kTRUE;
                    }
                }
            }
            if (relId[0] == 1)
            {
                if (fBadMapM4->GetBinContent(x+1, z+1) != 0)
                {

                    Float_t distance = TMath::Sqrt((x-locVec[0])*(x-locVec[0]) + (z - locVec[2])*(z-locVec[2]));
                    if (distance < fCuts->GetPhosBadDistanceCut())
                    {
                        return kTRUE;
                    }
                }
            }
        }
    }

    return kFALSE;
}
*/
