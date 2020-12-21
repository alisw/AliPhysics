#include <vector>
#include <iostream>

#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliEventClassifierMult.h"


using namespace std;

ClassImp(AliEventClassifierMult)

AliEventClassifierMult::AliEventClassifierMult(const char* name, const char* title,
						 //std::vector<std::vector<Float_t> > regions,
						 Float_t regions[],
						 Int_t lengthRegions,
						 Bool_t regionsAreInclusive,
						 Bool_t countCharged,
						 TList *taskOutputList,
                                                 Int_t collisionSystem)
  : AliEventClassifierBase(name, title, taskOutputList, collisionSystem),
  fRegionsAreInclusive(regionsAreInclusive),
  fCountCharged(countCharged)
{
  if(fCollisionSystem == 0) {
      fExpectedMinValue = 0;
      fExpectedMaxValue = 250;
  } else if (fCollisionSystem == 1) {
      fExpectedMinValue = 0;
      fExpectedMaxValue = 250;
  } else {
      fExpectedMinValue = 0;
      fExpectedMaxValue = 1500;
  }
  std::vector<Float_t> v(2);
  for (Int_t i=0; i<lengthRegions; i+=2){
    v[0] = regions[i];
    v[1] = regions[i+1];
    fRegions.push_back(v);
  }
}

void AliEventClassifierMult::CalculateClassifierValue(AliMCEvent *event, AliStack *stack) {
  fClassifierValue = 0.0;
  for (Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++) {
    AliMCParticle *track = static_cast<AliMCParticle*>(event->GetTrack(iTrack));
    // load track
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTrack);
      continue;
    }
    // discard unphysical particles from some generators
    if (track->Pt() == 0 || track->E() <= 0)
      continue;

    // Only calculate for primaries (Aliroot definition excluding Pi0)
    if (!stack->IsPhysicalPrimary(iTrack)) continue;

    // do we count charged or neutral?
    if (track->Charge() == 0 && fCountCharged) continue;

    // does this track fall into any of the defined regions?
    Double_t eta = track->Eta();
    Bool_t trackIsInRegion = false;
    for(Int_t i = 0; i != fRegions.size(); i++) {
      if(eta >= fRegions[i][0] && eta <=fRegions[i][1]) {
	trackIsInRegion = true;
      }
    }
    // Are we counting tracks inside or outside of the region?
    // increment counter accordingly
    if (trackIsInRegion && fRegionsAreInclusive) fClassifierValue += 1.0;
    else if (!trackIsInRegion && !fRegionsAreInclusive) fClassifierValue += 1.0;   
  }
}
