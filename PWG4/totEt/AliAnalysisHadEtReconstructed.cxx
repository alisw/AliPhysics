
//Create by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
#include "AliAnalysisHadEtReconstructed.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDtrack.h"
#include "AliESDCaloCluster.h"
#include "TVector3.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliVParticle.h"
#include <iostream>
#include "TH2F.h"

AliAnalysisHadEtReconstructed::AliAnalysisHadEtReconstructed() :
        AliAnalysisHadEt()
        ,fNTpcClustersCut(EtReconstructedCuts::kNTpcClustersCut)
        ,fNItsClustersCut(EtReconstructedCuts::knItsClustersCut)
        ,fTrackDistanceCut(0)
{

}

Int_t AliAnalysisHadEtReconstructed::AnalyseEvent(AliVEvent* ev)
{
    ResetEventValues();
    AliESDEvent *event = dynamic_cast<AliESDEvent*>(ev);

    for (Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++)
    {
        AliVParticle *track = event->GetTrack(iTrack);
        if (!track)
        {
            Printf("ERROR: Could not get track %d", iTrack);
            continue;
        }

        fMultiplicity++;

        const Double_t *pidWeights = track->PID();
        if (pidWeights)
        {
            Int_t maxpid = -1;
            Float_t maxpidweight = 0;
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
                //     massPart = -0.938*track->Charge();
            }

        }

    }

    fTotNeutralEtAcc = fTotNeutralEt;
    fTotEt = fTotChargedEt + fTotNeutralEt;
    fTotEtAcc = fTotChargedEtAcc + fTotNeutralEtAcc;

    std::cout << fTotChargedEtAcc << std::endl;
    // Fill the histograms...
    FillHistograms();

    return 1;
}

bool AliAnalysisHadEtReconstructed::CheckGoodVertex(AliVParticle* track)
{

    Float_t bxy = 999.;
    Float_t bz = 999.;
    dynamic_cast<AliESDtrack*>(track)->GetImpactParametersTPC(bxy,bz);

    return TMath::Abs(track->Xv()) < fVertexXCut && TMath::Abs(track->Yv()) < fVertexYCut && TMath::Abs(track->Zv()) < fVertexZCut && TMath::Abs(bxy) < fIPxyCut && TMath::Abs(bz) < fIPzCut;;

}

void AliAnalysisHadEtReconstructed::Init()
{

    AliAnalysisHadEt::Init();

    fVertexXCut = EtReconstructedCuts::kVertexXCut;
    fVertexYCut = EtReconstructedCuts::kVertexYCut;
    fVertexZCut = EtReconstructedCuts::kVertexZCut;
    fIPxyCut = EtReconstructedCuts::kIPxyCut;
    fIPzCut = EtReconstructedCuts::kIPzCut;

}


