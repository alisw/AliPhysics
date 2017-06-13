#include <vector>
#include <iostream>

#include "TMath.h"

#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenPythiaEventHeader.h"

#include "AliEventClassifierSphericity.h"
#include "AliIsPi0PhysicalPrimary.h"

using namespace std;

ClassImp(AliEventClassifierSphericity)

AliEventClassifierSphericity::AliEventClassifierSphericity(const char* name, const char* title,
							   TList *taskOutputList)
  : AliEventClassifierBase(name, title, taskOutputList)
{
  fExpectedMinValue = 0;
  fExpectedMaxValue = 1;
}

void AliEventClassifierSphericity::CalculateClassifierValue(AliMCEvent *event, AliStack *stack) {
  // This implementation is adapted from PWGLF/SPECTRA/Spherocity/AliTransverseEventShape.cxx
  fClassifierValue = -1.0;

  Float_t sphericity = -1.0;
  Float_t s00=0;
  Float_t s01=0;
  Float_t s11=0;
  Float_t totalpt=0;

  Int_t ntracks = event->GetNumberOfTracks();
  for (Int_t iTrack = 0; iTrack < ntracks; iTrack++) {
    AliMCParticle *track = static_cast<AliMCParticle*>(event->GetTrack(iTrack));
    if (!track)
      continue;
    // Only calculate for primaries (Aliroot definition excluding Pi0)
    if (!stack->IsPhysicalPrimary(iTrack))
      continue;
    // discard unphysical particles from some generators
    if (track->Pt() == 0 || track->E() <= 0)
      continue;
    
    Float_t px = track->Pt() * TMath::Cos( track->Phi() );
    Float_t py = track->Pt() * TMath::Sin( track->Phi() );
    s00 += (px * px) / track->Pt();
    s01 += (py * px) / track->Pt();
    s11 += (py * py) / track->Pt();
    totalpt += track->Pt();
  }
  // did we have valid tracks or did we never reach the bottom of the for loop?
  if (!(totalpt > 0)) {
    fClassifierValue = -1;
    return;
  }

  Double_t S00=s00/totalpt;
  Double_t S01=s01/totalpt;
  Double_t S11=s11/totalpt;

  Float_t lambda1=((S00+S11)+TMath::Sqrt((S00+S11)*(S00+S11)-4*(S00*S11-S01*S01)))/2;
  Float_t lambda2=((S00+S11)-TMath::Sqrt((S00+S11)*(S00+S11)-4*(S00*S11-S01*S01)))/2;
  if((lambda2==0)&&(lambda1==0))
    sphericity=0;
  if(lambda1+lambda2!=0)
    sphericity=2*TMath::Min( lambda1,lambda2 )/( lambda1+lambda2 );
  
  // Compute the final sphericity:
  fClassifierValue = sphericity;
}
