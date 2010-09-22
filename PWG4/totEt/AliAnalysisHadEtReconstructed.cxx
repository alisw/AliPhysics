//_________________________________________________________________________
//  Utility Class for transverse energy studies, charged hadrons
//  Base class for ESD analysis
//  - reconstruction output
// implementation file
//
//Created by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
//_________________________________________________________________________
#include "AliAnalysisHadEtReconstructed.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDtrack.h"
#include "AliESDCaloCluster.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliVParticle.h"
#include <iostream>

using namespace std;

ClassImp(AliAnalysisHadEtReconstructed);


AliAnalysisHadEtReconstructed::AliAnalysisHadEtReconstructed() :
        AliAnalysisHadEt()
{

}

AliAnalysisHadEtReconstructed::~AliAnalysisHadEtReconstructed() 
{
}

Int_t AliAnalysisHadEtReconstructed::AnalyseEvent(AliVEvent* ev)
{ // analyse ESD event
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
{ // check vertex

    Float_t bxy = 999.;
    Float_t bz = 999.;
    dynamic_cast<AliESDtrack*>(track)->GetImpactParametersTPC(bxy,bz);

    bool status = (TMath::Abs(track->Xv()) < fCuts->GetReconstructedVertexXCut()) && 
      (TMath::Abs(track->Yv()) < fCuts->GetReconstructedVertexYCut()) && 
      (TMath::Abs(track->Zv()) < fCuts->GetReconstructedVertexZCut()) && 
      (TMath::Abs(bxy) < fCuts->GetReconstructedIPxyCut()) && 
      (TMath::Abs(bz) < fCuts->GetReconstructedIPzCut()); 

    return status;
}

void AliAnalysisHadEtReconstructed::Init()
{ // Init
    AliAnalysisHadEt::Init();
}


