#include <vector>
#include <iostream>

#include "TMath.h"

#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenPythiaEventHeader.h"

#include "AliEventClassifierSpherocity.h"
#include "AliIsPi0PhysicalPrimary.h"

using namespace std;

ClassImp(AliEventClassifierSpherocity)

AliEventClassifierSpherocity::AliEventClassifierSpherocity(const char* name, const char* title,
					     TList *taskOutputList)
  : AliEventClassifierBase(name, title, taskOutputList)
{
  fExpectedMinValue = 0;
  fExpectedMaxValue = 1;
}

Bool_t AliEventClassifierSpherocity::TrackPassesSelection(AliMCParticle* track, AliStack *stack, Int_t iTrack) {
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTrack);
      return false;
    }
    // discard unphysical particles from some generators
    if (track->Pt() == 0 || track->E() <= 0)
      return false;

    // Only calculate for primaries (Aliroot definition excluding Pi0)
    if (!stack->IsPhysicalPrimary(iTrack)) return false;

    // Restrict to |eta| < 0.8
    if (TMath::Abs(track->Eta()) > 0.8) return false;

    return true;
}

void AliEventClassifierSpherocity::CalculateClassifierValue(AliMCEvent *event, AliStack *stack) {
  // This implementation is adapted from PWGLF/SPECTRA/Spherocity/AliTransverseEventShape.cxx
  fClassifierValue = 0.0;

  Float_t sumRatioSquare = 0;
  Float_t minimalSumRatioSquare = 2;
  
  // Step size in phi unit vector used to find m spherocity
  Float_t phiStepSize = 0.1;

  // Computing total pt
  Float_t sumapt = 0;
  Int_t ntracks = event->GetNumberOfTracks();
  for (Int_t iTrack = 0; iTrack < ntracks; iTrack++) {
    AliMCParticle *track = static_cast<AliMCParticle*>(event->GetTrack(iTrack));
    if (!TrackPassesSelection(track, stack, iTrack)) continue;
    sumapt += track->Pt();
  }
  
  // Getting thrust
  for(Int_t i = 0; i < 360/(phiStepSize); ++i){
    Float_t numerator = 0;
    Float_t phiparam  = 0;
    Float_t nx = 0;
    Float_t ny = 0;
    phiparam=((TMath::Pi()) * i * phiStepSize) / 180; // parametrization of the angle
    nx = TMath::Cos(phiparam);            // x component of an unitary vector n
    ny = TMath::Sin(phiparam);            // y component of an unitary vector n
    for(Int_t iTrack = 0; iTrack < ntracks; ++iTrack){
      AliMCParticle *track = static_cast<AliMCParticle*>(event->GetTrack(iTrack));
      if (!TrackPassesSelection(track, stack, iTrack)) continue;
      
      Float_t pxA = track->Pt() * TMath::Cos(track->Phi());
      Float_t pyA = track->Pt() * TMath::Sin(track->Phi());
      //product between p projection in XY plane and the unitary vector
      numerator += TMath::Abs( ny * pxA - nx * pyA );
    }
    
    sumRatioSquare = TMath::Power((numerator / sumapt), 2);
    if(sumRatioSquare < minimalSumRatioSquare)  //minimization of pFull
      {
	minimalSumRatioSquare = sumRatioSquare;
      }
  }

  // Compute the final spherocity:
  fClassifierValue = (minimalSumRatioSquare * TMath::Pi() * TMath::Pi()) / 4.0;
}
