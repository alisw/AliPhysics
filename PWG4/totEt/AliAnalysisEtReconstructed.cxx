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
#include "AliESDCaloCluster.h"
#include "TVector3.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliVParticle.h"
#include <iostream>
#include "TH2F.h"

AliAnalysisEtReconstructed::AliAnalysisEtReconstructed() :
        AliAnalysisEt()
        ,fTrackDistanceCut(0)
        ,fClusterType(0)
{

}

AliAnalysisEtReconstructed::~AliAnalysisEtReconstructed() 
{
}

Int_t AliAnalysisEtReconstructed::AnalyseEvent(AliVEvent* ev)
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

        Int_t nItsClusters = dynamic_cast<AliESDtrack*>(track)->GetNcls(0);
        Int_t nTPCClusters = dynamic_cast<AliESDtrack*>(track)->GetNcls(1);

        Float_t massPart = 0;

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

        Double_t et = track->E() * TMath::Sin(track->Theta()) + massPart;
	// printf("Rec track: iTrack %03d eta %4.3f phi %4.3f nITSCl %d nTPCCl %d\n", iTrack, track->Eta(), track->Phi(), nItsClusters, nTPCClusters); // tmp/debug printout

        if (TMath::Abs(track->Eta()) < fCuts->GetCommonEtaCut() && CheckGoodVertex(track) && nItsClusters > fCuts->GetReconstructedNItsClustersCut() && nTPCClusters > fCuts->GetReconstructedNTpcClustersCut() )
        {
	    fTotChargedEt +=  et;
            fChargedMultiplicity++;

            if (TMath::Abs(track->Eta()) < fEtaCutAcc && track->Phi() < fPhiCutAccMax && track->Phi() > fPhiCutAccMin)
            {
                fTotChargedEtAcc += track->E()*TMath::Sin(track->Theta()) + massPart;
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
            Printf("ERROR: Could not get cluster %d", iCluster);
            continue;
        }

	// printf("Rec Cluster: iCluster %03d E %4.3f type %d NCells %d\n", iCluster, cluster->E(), (int)(cluster->GetType()), cluster->GetNCells()); // tmp/debug printout
        if (cluster->GetType() != fClusterType) continue;

        if (cluster->E() < fClusterEnergyCut) continue;
        Float_t pos[3];
        TVector3 cp(pos);
        cluster->GetPosition(pos);
        //if (pos[0] < -(32.0*2.2)) continue; //Ensure that modules 0 and 1 are not used
        // if(cp.Phi() < 260.*TMath::Pi()/180.) continue;
        fHistTMDeltaR->Fill(cluster->GetEmcCpvDistance());
        if (cluster->GetEmcCpvDistance() < fTrackDistanceCut)
        {
            continue;
            //AliVParticle *matchedTrack = event->GetTrack(cluster->GetTrackMatched());
// 	    if(CheckGoodVertex(matchedTrack))
// 	    {
// 	       totChargedEnergy +=  matchedTrack->E();;
// 	       totChargedEt += matchedTrack->E()*TMath::Sin(matchedTrack);
// 	    }
        }

        if (cluster->E() >  fSingleCellEnergyCut && cluster->GetNCells() == fCuts->GetCommonSingleCell()) continue;

        cluster->GetPosition(pos);
      
	// TODO: replace with TVector3, too lazy now...

        float dist = TMath::Sqrt(pos[0]*pos[0] + pos[1]*pos[1]);

        float theta = TMath::ATan(pos[2]/dist)+TMath::Pi()/2;
        // float eta = TMath::Log(TMath::Abs( TMath::Tan( 0.5 * theta ) ) );
        fTotNeutralEt += cluster->E() * TMath::Sin(theta);
	fNeutralMultiplicity++;

        fMultiplicity++;
    }

    fTotNeutralEtAcc = fTotNeutralEt;

    fTotEt = fTotChargedEt + fTotNeutralEt;
    fTotEtAcc = fTotChargedEtAcc + fTotNeutralEtAcc;

    std::cout << fTotChargedEtAcc << std::endl;
    // Fill the histograms...
    FillHistograms();

    return 0;
}

bool AliAnalysisEtReconstructed::CheckGoodVertex(AliVParticle* track)
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

void AliAnalysisEtReconstructed::Init()
{ // Init
    AliAnalysisEt::Init();
}

bool AliAnalysisEtReconstructed::TrackHitsCalorimeter(AliVParticle* track, Double_t magField)
{ // propagate track to detector radius

   AliESDtrack *esdTrack = dynamic_cast<AliESDtrack*>(track);
   // Printf("Propagating track: eta: %f, phi: %f, pt: %f", esdTrack->Eta(), esdTrack->Phi(), esdTrack->Pt());

    Bool_t prop = esdTrack->PropagateTo(fDetectorRadius, magField);

    // if (prop) Printf("Track propagated, eta: %f, phi: %f, pt: %f", esdTrack->Eta(), esdTrack->Phi(), esdTrack->Pt());
    return prop && 
		   TMath::Abs(esdTrack->Eta()) < fEtaCutAcc && 
		   esdTrack->Phi() > fPhiCutAccMin*TMath::Pi()/180. && 
		   esdTrack->Phi() < fPhiCutAccMax*TMath::Pi()/180.;
}

