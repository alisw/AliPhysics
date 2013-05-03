/**************************************************************************
 * Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

///////////////////////////////////////////////////////////////////////////////
//
// Track matching between TRD online tracks and ESD tracks.
//
// Author: Felix Rettig <rettig@compeng.uni-frankfurt.de>
//
///////////////////////////////////////////////////////////////////////////////

#include <TH1.h>
#include <AliESDEvent.h>
#include <AliExternalTrackParam.h>
#include "AliESDtrack.h"
#include "AliESDTrdTrack.h"
#include <AliGeomManager.h>
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"
#include "AliTRDonlineTrackMatching.h"

const Float_t AliTRDonlineTrackMatching::fgkSaveInnerRadius = 290.5;
const Float_t AliTRDonlineTrackMatching::fgkSaveOuterRadius = 364.5;

Float_t AliTRDonlineTrackMatching::fEsdTrackCutMinTPCrows = 0.;
Float_t AliTRDonlineTrackMatching::fEsdTrackCutMinRatioRowsFindableClusters = 0.;
Float_t AliTRDonlineTrackMatching::fEsdTrackCutMaxChi2TPCclusters = 0.;
Float_t AliTRDonlineTrackMatching::fEsdTrackCutMaxChi2ITSclusters = 0.;
Float_t AliTRDonlineTrackMatching::fEsdTrackCutMaxDCAtoVertexXY = 0.;
Float_t AliTRDonlineTrackMatching::fEsdTrackCutMaxDCAtoVertexZ = 0.;
UShort_t AliTRDonlineTrackMatching::fEsdTrackCutsITSlayerMask = 0;  // similar to 2011 default cut: 0x3
Float_t AliTRDonlineTrackMatching::fEsdTrackVCutsChi2TPCconstrainedVsGlobal = 0.;
Float_t	AliTRDonlineTrackMatching::fEsdTrackCutPtDCAOfs = 0.;
Float_t	AliTRDonlineTrackMatching::fEsdTrackCutPtDCACoeff = 0.;
Bool_t AliTRDonlineTrackMatching::fEsdTrackCutMinimal = kFALSE;
Bool_t AliTRDonlineTrackMatching::fEsdTrackCutRequireTPCrefit = kTRUE;
Bool_t AliTRDonlineTrackMatching::fEsdTrackCutRequireITSrefit = kFALSE;
Bool_t AliTRDonlineTrackMatching::fEsdTrackCutPrim = kFALSE;

AliTRDonlineTrackMatching::AliTRDonlineTrackMatching() :
  TObject(),
  fTRDgeo(NULL),
  fMinMatchRating(0.25),
  fHistMatchRating(NULL)
{
  // default ctor
  SetEsdTrackDefaultCuts("minimal");
}

AliTRDonlineTrackMatching::AliTRDonlineTrackMatching(const AliTRDonlineTrackMatching &c) :
  TObject(c),
  fTRDgeo(c.fTRDgeo),
  fMinMatchRating(c.fMinMatchRating),
  fHistMatchRating(c.fHistMatchRating)
{
  // copy ctor
}

AliTRDonlineTrackMatching::~AliTRDonlineTrackMatching() {

  // dtor

  delete fTRDgeo;
  fTRDgeo = NULL;
}

Short_t AliTRDonlineTrackMatching::EstimateSector(const Double_t globalCoords[3]) {

  // estimates sector by phi angle in x-y plane

  if ((TMath::Abs(globalCoords[0]) > 600) || (TMath::Abs(globalCoords[0]) > 600) || (TMath::Sqrt(globalCoords[0]*globalCoords[0] + globalCoords[1]*globalCoords[1]) < 0.01)){
    //printf("GGG %.3f/%.3f\n", globalCoords[0], globalCoords[1]);
    return -1;
  } else {
    Double_t ang = TMath::ATan2(globalCoords[1], globalCoords[0]);
    if (ang > 0){
#ifdef TRD_TM_DEBUG
      printf("	  es: %.2f/%.2f	 -> phi: %.2fdeg -> Sec %02d  (A)\n",
	     globalCoords[0], globalCoords[1], TMath::ATan2(globalCoords[1], globalCoords[0])*180./TMath::Pi(),
	     TMath::FloorNint(ang/(20./180.*TMath::Pi())));
#endif
      return TMath::FloorNint(ang/(20./180.*TMath::Pi()));
    } else {
#ifdef TRD_TM_DEBUG
      printf("	  es: %.2f/%.2f	 -> phi: %.2fdeg -> Sec %02d  (B)\n",
	     globalCoords[0], globalCoords[1], TMath::ATan2(globalCoords[1], globalCoords[0])*180./TMath::Pi(),
	     17 - TMath::FloorNint(TMath::Abs(ang)/(20./180.*TMath::Pi())));
#endif
      return 17 - TMath::FloorNint(TMath::Abs(ang)/(20./180.*TMath::Pi()));
    }

  }
}

Short_t AliTRDonlineTrackMatching::EstimateLayer(Double_t radius) {

  // estimates layer by radial distance (for virtual stack at phi = 0)

  const Float_t rBoundaries[7] = {290.80, 302.20, 315.06, 327.55, 340.3, 352.80, 364.15}; // radial border lines centered between anode plane and successing radiator
  const Short_t rLayers[7] = {-1, 0, 1, 2, 3, 4, 5};
  for (UShort_t i = 0; i < 7; ++i){
    if (radius < rBoundaries[i])
      return rLayers[i];
  }
  return -2; // radius larger than outmost layer
}

Short_t AliTRDonlineTrackMatching::EstimateLocalStack(const Double_t globalCoords[3]) {

  // determines stack within sector by z position

  Double_t absZ = TMath::Abs(globalCoords[2]);
  Short_t signZ = (globalCoords[2] > 0.) ? 1 : -1;
  Double_t r = TMath::Sqrt(globalCoords[0]*globalCoords[0] + globalCoords[1]*globalCoords[1]);
  Short_t layer = EstimateLayer(r);

#ifdef TRD_TM_DEBUG
  printf("EstimateLocalStack A	r: %.2f	  x: %.2f/%.2f/%.2f  -> layer: %i    absZ = %.2f\n",
	 r, globalCoords[0], globalCoords[1], globalCoords[2], layer, absZ);
#endif

  if (layer < 0)
    return -1;

  Double_t innerStackHalfLength = AliTRDgeometry::GetChamberLength(0, 2) / 2.;  // same for all layers
  if (absZ < innerStackHalfLength)
    return 2;

  Double_t outerStackLength = AliTRDgeometry::GetChamberLength(layer, 1);

  absZ -= innerStackHalfLength;

#ifdef TRD_TM_DEBUG
  printf("EstimateLocalStack B	r: %.2f	  x: %.2f/%.2f/%.2f  -> layer: %i    absZ = %.2f    il: %.2f   ol: %.2f\n",
	 r, globalCoords[0], globalCoords[1], globalCoords[2], layer, absZ, 2.*innerStackHalfLength, outerStackLength);
#endif

  if (absZ > 2.05*outerStackLength)
    return (signZ > 0) ? -2 : -1; // outside supermodule in z direction

  if (absZ < outerStackLength)
    return (signZ > 0) ? 1 : 3;
  else
    return (signZ > 0) ? 0 : 4;

}

Short_t AliTRDonlineTrackMatching::EstimateStack(const Double_t globalCoords[3]) {

  // returns the closest TRD stack to a 3D position in global coordinates

  Short_t sec = EstimateSector(globalCoords);
  Short_t st = EstimateLocalStack(globalCoords);
#ifdef TRD_TM_DEBUG
  printf("EstimateStack sec %d  st %d\n", sec, st);
#endif
  if ((sec < 0) || (st < 0))
    return -1;
  else
    return 5*sec + st;
}

Bool_t AliTRDonlineTrackMatching::StackToTrack(const AliExternalTrackParam *track, Short_t &stack, UShort_t &layersWithTracklet, Double_t magFieldinKiloGauss){

  // returns stack to track param

  stack = -1;
  layersWithTracklet = 0;

  UInt_t stackHits[fgkTrdStacks];
  Double_t x[3];
  memset(stackHits, 0, fgkTrdStacks*sizeof(UInt_t));

#ifdef TRD_TM_DEBUG
  printf("STACK-TO-TRACK\n");
#endif

  Double_t r = fgkSaveInnerRadius;
  while (r < fgkSaveOuterRadius){
    track->GetXYZAt(r, magFieldinKiloGauss, x);
    stack = EstimateStack(x);
    if (stack >= 0){
      stackHits[stack]++;
      if (stackHits[stack] > 16) // experimental
	break;
#ifdef TRD_TM_DEBUG
      printf(" r=%.3fcm  %.2f/%.2f  -  %d hits for stack %d  S%02d-%d   (mag=%.1f)\n",
	     r, x[0], x[1], stackHits[stack], stack, stack/5, stack%5, magFieldinKiloGauss);
#endif
    }
    r += 1.;
  }

  // find stack with most hits
  UInt_t bestHits = 0;
  for (UShort_t iStack = 0; iStack < fgkTrdStacks; ++iStack){
    if (stackHits[iStack] == 0)
      continue;
#ifdef TRD_TM_DEBUG
    printf("  finally %d hits in stack S%02d-%d\n", stackHits[iStack], iStack/5, iStack%5);
#endif
    if (stackHits[iStack] > bestHits){
      bestHits = stackHits[iStack];
      stack = iStack;
    }
  }

  if (stack >= 0){
#ifdef TRD_TM_DEBUG
    printf("best stack: S%02d-%d\n", TrdLsiSec(stack), TrdLsiSi(stack));
#endif
    return kTRUE;
  }

  return kFALSE;
}

Bool_t AliTRDonlineTrackMatching::StackToTrack(const AliESDtrack* track, Short_t &stack, UShort_t &layersWithTracklet, Double_t magFieldinKiloGauss){

  // returns stack to ESD track

  if (track->GetOuterParam())
    return StackToTrack(track->GetOuterParam(), stack, layersWithTracklet, magFieldinKiloGauss);
  else if (track->GetInnerParam())
    return StackToTrack(track->GetInnerParam(), stack, layersWithTracklet, magFieldinKiloGauss);
  else
    return StackToTrack(track, stack, layersWithTracklet, magFieldinKiloGauss);
}

Bool_t AliTRDonlineTrackMatching::AcceptTrack(const AliESDtrack* esdTrack, const AliESDEvent* esdEvent){

  // returns result ESD track cuts

  if (!esdTrack)
    return kFALSE;

  UInt_t status = esdTrack->GetStatus();

  if (fEsdTrackCutMinimal){
    return ((status & AliESDtrack::kTPCout) > 0);
  }

  // require TPC fit
  if ((fEsdTrackCutRequireTPCrefit) && (!(status & AliESDtrack::kTPCrefit)))
    return kFALSE;

  // require ITS re-fit
  if ((fEsdTrackCutRequireITSrefit) && (!(status & AliESDtrack::kITSrefit)))
    return kFALSE;

  // TPC requirements
  Float_t nCrossedRowsTPC = esdTrack->GetTPCCrossedRows();
  Float_t ratioCrossedRowsOverFindableClustersTPC =
    (esdTrack->GetTPCNclsF() > 0) ? (nCrossedRowsTPC / esdTrack->GetTPCNclsF()) : 1.0;
  Float_t chi2PerClusterTPC =
    (esdTrack->GetTPCclusters(0) > 0) ? (esdTrack->GetTPCchi2()/Float_t(esdTrack->GetTPCclusters(0))) : 100.;

  if (
      (nCrossedRowsTPC < fEsdTrackCutMinTPCrows) ||
      (ratioCrossedRowsOverFindableClustersTPC < fEsdTrackCutMinRatioRowsFindableClusters) ||
      (chi2PerClusterTPC > fEsdTrackCutMaxChi2TPCclusters)
      )
    return kFALSE;

  // ITS requirements
  Float_t chi2PerClusterITS = (esdTrack->GetITSclusters(0) > 0) ? esdTrack->GetITSchi2()/Float_t(esdTrack->GetITSclusters(0)) : 1000.;
  UShort_t clustersInAnyITSlayer = kFALSE;
  for (UShort_t layer = 0; layer < 6; ++layer)
    clustersInAnyITSlayer += (esdTrack->HasPointOnITSLayer(layer) & ((fEsdTrackCutsITSlayerMask >> layer) & 1));

  if ((fEsdTrackCutsITSlayerMask != 0) &&
      ((clustersInAnyITSlayer == 0) || (chi2PerClusterITS >= fEsdTrackCutMaxChi2ITSclusters))
      )
    return kFALSE;

  // geometric requirements
  Float_t impactPos[2], impactCov[3];
  esdTrack->GetImpactParameters(impactPos, impactCov);

  if (TMath::Abs(impactPos[0]) > fEsdTrackCutMaxDCAtoVertexXY)
    return kFALSE;

  if (TMath::Abs(impactPos[1]) > fEsdTrackCutMaxDCAtoVertexZ)
    return kFALSE;

  if (fEsdTrackCutPrim){
    // additional requirements for primary tracks

    const AliESDVertex* vertex = esdEvent->GetPrimaryVertexTracks();
    if ((!vertex) || (!vertex->GetStatus()))
      vertex = esdEvent->GetPrimaryVertexSPD();

    Float_t chi2TPCConstrainedVsGlobal =
      (vertex->GetStatus()) ? esdTrack->GetChi2TPCConstrainedVsGlobal(vertex) : (fEsdTrackVCutsChi2TPCconstrainedVsGlobal + 10.);

    if (chi2TPCConstrainedVsGlobal > fEsdTrackVCutsChi2TPCconstrainedVsGlobal)
      return kFALSE;

    Float_t cutDCAToVertexXYPtDep =
      fEsdTrackCutPtDCAOfs + fEsdTrackCutPtDCACoeff/((TMath::Abs(esdTrack->Pt()) > 0.0001) ? esdTrack->Pt() : 0.0001);

    if (TMath::Abs(impactPos[0]) >= cutDCAToVertexXYPtDep)
      return kFALSE;

  }

  return kTRUE;
}

Bool_t AliTRDonlineTrackMatching::ProcessEvent(AliESDEvent *esdEvent) {

  // performs track matching for all TRD online tracks of the ESD event

  UInt_t numTrdTracks = esdEvent->GetNumberOfTrdTracks();
  if (numTrdTracks <= 0)
    return kTRUE;

  if (!AliGeomManager::GetGeometry()){
    AliError("Geometry not available! Skipping TRD track matching.");
    return kFALSE;
  }

  if (!fTRDgeo){
    fTRDgeo = new AliTRDgeometry();
  }

  //
  // ESD track selection and sorting by TRD stack
  //

  UInt_t esdTracksByStack[fgkTrdStacks][fgkMaxEsdTracksPerStack];
  UInt_t esdTrackNumByStack[fgkTrdStacks];
  memset(esdTrackNumByStack, 0, fgkTrdStacks*sizeof(UInt_t));

  UInt_t numEsdTracks = esdEvent->GetNumberOfTracks();
#ifdef TRD_TM_DEBUG
  UInt_t numEsdTracksAccepted = 0;
#endif
  Short_t stack;
  UShort_t layers;
  AliESDtrack* esdTrack;

  for (UInt_t iEsdTrack = 0; iEsdTrack < numEsdTracks; ++iEsdTrack){
    esdTrack = esdEvent->GetTrack(iEsdTrack);

    if (!esdTrack){
      AliError("invalid ESD track!");
      continue;
    }

    // track filter here
    if (!AcceptTrack(esdTrack, esdEvent))
      continue;
#ifdef TRD_TM_DEBUG
    else
      numEsdTracksAccepted++;
#endif

    // assign ESD track to TRD stack
    if (StackToTrack(esdTrack, stack, layers, esdEvent->GetMagneticField())){

      if (stack < 0){
#ifdef TRD_TM_DEBUG
	printf("#TRACKMATCHING - invalid stack for ESD track\n");
#endif
	continue;
      }

      // register track in relevant stacks
      Int_t stacksForReg[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
      stacksForReg[0] = stack; // stack hit
      stacksForReg[1] = (stack + 5) % 90; // same stack in next supermodule
      stacksForReg[2] = (stack - 5); // same stack in previous supermodule
      if (stacksForReg[2] < 0)
	stacksForReg[2] += 90;

      switch(TrdLsiSi(stack)){
      case 0:
	// stack 0
	stacksForReg[3] = stack + 1; // next stack in same supermodule
	stacksForReg[4] =  stacksForReg[1] + 1; // next stack in next supermodule
	stacksForReg[5] =  stacksForReg[2] + 1; // next stack in previous supermodule
	break;
      case 1:
      case 2:
      case 3:
	stacksForReg[3] = stack + 1; // next stack in same supermodule
	stacksForReg[4] =  stacksForReg[1] + 1; // next stack in next supermodule
	stacksForReg[5] =  stacksForReg[2] + 1; // next stack in previous supermodule
	stacksForReg[6] = stack - 1; // previous stack in same supermodule
	stacksForReg[7] =  stacksForReg[1] - 1; // previous stack in next supermodule
	stacksForReg[8] =  stacksForReg[2] - 1; // previous stack in previous supermodule
	break;
      case 4:
	stacksForReg[3] = stack - 1; // previous stack in same supermodule
	stacksForReg[4] =  stacksForReg[1] - 1; // previous stack in next supermodule
	stacksForReg[5] =  stacksForReg[2] - 1; // previous stack in previous supermodule
	break;
      default:
	break;
      }

#ifdef TRD_TM_DEBUG
      printf("#TRACKMATCHING - assigned ESD track %d to following TRD stacks:", iEsdTrack);
#endif

      // register for stacks
      for (UShort_t iReg = 0; iReg < 9; ++iReg){
	if (stacksForReg[iReg] < 0)
	  break;

	if (stacksForReg[iReg] >= 90){
	  AliError(Form("invalid stack for registration: %i", stacksForReg[iReg]));
	  continue;
	}

	if (esdTrackNumByStack[stacksForReg[iReg]] < fgkMaxEsdTracksPerStack - 1)
	  esdTracksByStack[stacksForReg[iReg]][esdTrackNumByStack[stacksForReg[iReg]]++] = iEsdTrack;
#ifdef TRD_TM_DEBUG
	else
	  printf("#TRACKMATCHING - maximum number (%d) of ESD tracks per stack reached for S%02d-%d (%d tracks total). Skipping track!\n",
		 fgkMaxEsdTracksPerStack, TrdLsiSec(stacksForReg[iReg]), TrdLsiSi(stacksForReg[iReg]), numEsdTracks);
	printf(" S%02d-%d", TrdLsiSec(stacksForReg[iReg]), TrdLsiSi(stacksForReg[iReg]));
#endif
      }
#ifdef TRD_TM_DEBUG
      printf(" (ESD-ASSIGN)\n");
#endif

//      if (esdTrackNumByStack[stack] >= fgkMaxEsdTracksPerStack){
//#ifdef TRD_TM_DEBUG
// 	printf("#TRACKMATCHING - maximum number (%d) of ESD tracks per stack reached for S%02d-%d (%d tracks total). Skipping track!\n",
// 	       fgkMaxEsdTracksPerStack, TrdLsiSec(stack), TrdLsiSi(stack), numEsdTracks);
//#endif
// 	continue;
//      }
//
//      esdTracksByStack[stack][esdTrackNumByStack[stack]++] = iEsdTrack;
//#ifdef TRD_TM_DEBUG
//      printf("#TRACKMATCHING - assigned ESD track %d to TRD stack S%02d-%d\n",
// 	     iEsdTrack, TrdLsiSec(stack), TrdLsiSi(stack));
//#endif
    }

  } // loop over esd tracks

#ifdef TRD_TM_DEBUG
  printf("#TRACKMATCHING - %d ESD tracks accepted, %d rejected\n",
	 numEsdTracksAccepted, numEsdTracks - numEsdTracksAccepted);
#endif

  //
  // search matching ESD track for each TRD online track
  //
  AliESDTrdTrack* trdTrack;
  Double_t trdPt;
  AliESDtrack* matchCandidate;
  AliESDtrack* matchTrack;
  Int_t matchEsdTrackIndexInStack;
  Double_t matchRating;
  Int_t matchCandidateCount;
  Double_t distY, distZ;

  for (UInt_t iTrdTrack = 0; iTrdTrack < numTrdTracks; ++iTrdTrack){

    trdTrack = esdEvent->GetTrdTrack(iTrdTrack);
    stack = TrdSecSiLsi(trdTrack->GetSector(), trdTrack->GetStack());
    trdPt = (esdEvent->GetMagneticField() > 0.) ? (-1.*trdTrack->Pt()) : trdTrack->Pt();
    matchTrack = NULL;
    matchEsdTrackIndexInStack = -1;
    matchRating = 0.;
    matchCandidateCount = 0;

#ifdef TRD_TM_DEBUG
    printf("#TRACKMATCHING - trying to match TRD online track %d in S%02d-%d\n",
	   iTrdTrack, trdTrack->GetSector(), trdTrack->GetStack());
#endif

    // loop over all esd tracks in the same stack and check distance
    for (UInt_t iEsdTrack = 0; iEsdTrack < esdTrackNumByStack[stack]; ++iEsdTrack){
      matchCandidate = esdEvent->GetTrack(esdTracksByStack[stack][iEsdTrack]);

      if (EstimateTrackDistance(matchCandidate, trdTrack, esdEvent->GetMagneticField(), &distY, &distZ) == 0){
	Double_t rating = RateTrackMatch(distY, distZ, matchCandidate->GetSignedPt(), trdPt);
#ifdef TRD_TM_DEBUG
	printf("#TRACKMATCHING  S%02d-%d  trd %d - esd %d   dy: %.3f    dz: %.3f   r: %.3f    pt e: %.2f  t: %.2f   match: %d\n",
	       trdTrack->GetSector(), trdTrack->GetStack(), iTrdTrack, iEsdTrack,
	       distY, distZ, rating, matchCandidate->GetSignedPt(), trdPt,
	       (rating >= fMinMatchRating) ? 1 : 0);
#endif
	if (rating > 0.){
	  // possibly matching pair found
	  matchCandidateCount++;
	  if ((matchTrack == NULL) || (rating > matchRating)){
	    // new best match
	    matchTrack = matchCandidate;
	    matchEsdTrackIndexInStack = iEsdTrack;
	    matchRating = rating;
	  }
	}

      } else {
	// estimation of distance failed
#ifdef TRD_TM_DEBUG
	printf("TRACKMATCHING  S%02d-%d  trd %d - esd %d   failed\n",
	       trdTrack->GetSector(), trdTrack->GetStack(), iTrdTrack, iEsdTrack);
#endif
      }
    } // loop over esd tracks in same stack

    if (fHistMatchRating){
      fHistMatchRating->Fill(matchRating);
    }

    if ((matchTrack) && (matchRating >= fMinMatchRating)){
#ifdef TRD_TM_DEBUG
      printf("#TRACKMATCHING  S%02d-%d  trd %d - esd %d   match!    pt:  %.2f  %.2f\n",
	     trdTrack->GetSector(), trdTrack->GetStack(), iTrdTrack, matchEsdTrackIndexInStack,
	     trdPt, matchTrack->GetSignedPt());
#endif
      trdTrack->SetTrackMatchReference(matchTrack);
    } else
      trdTrack->SetTrackMatchReference(NULL);

  } // loop over TRD online tracks

  return kTRUE;
}

Bool_t AliTRDonlineTrackMatching::TrackPlaneIntersect(AliExternalTrackParam *trk, Double_t pnt[3], Double_t norm[3], Double_t mag){

  // calculates the intersection point of a track param and a plane defined by point pnt and normal vector norm

  UInt_t its = 0;
  Double_t r = 290.;
  Double_t step = 10.;
  Int_t flag = 0;
  Double_t dist = 0, dist_prev = 0;
  Double_t x[3] = {0., 0., 0.};

  dist = (x[0] - pnt[0]) * norm[0] + (x[1] - pnt[1]) *norm[1] + (x[2] - pnt[2]) * norm[2];

  while(TMath::Abs(dist) > 0.1) {

    trk->GetXYZAt(r, mag, x);

    if ((x[0] * x[0] + x[1] * x[1]) < 100.)  // extrapolation to radius failed
      return kFALSE;

    //distance between current track position and plane
    dist_prev = TMath::Abs(dist);
    dist = (x[0] - pnt[0]) * norm[0] + (x[1] - pnt[1]) * norm[1];
    if ((flag) && (TMath::Abs(dist) > dist_prev)){
      step /= -2.;
    }
    flag=1;
    r += step;
    its++;
    if ((r > 380.) || (r < 100.) || (its > 100) || (TMath::Abs(step) < 0.00001)){
      break;
    }
  }
  for (Int_t i=0; i<3; i++)
    pnt[i] = x[i];

  return kTRUE;
}

Int_t AliTRDonlineTrackMatching::EstimateTrackDistance(AliESDtrack *esd_track, AliESDTrdTrack* gtu_track, Double_t mag, Double_t *ydist, Double_t *zdist){

  // returns an estimate for the spatial distance between TPC offline track and GTU online track

  if ((!esd_track) || (!gtu_track))
    return -3;

  // AssertTRDGeometry();
  if (!fTRDgeo)
    fTRDgeo = new AliTRDgeometry();

  Float_t diff_y = 0;
  Float_t diff_z = 0;
  Int_t nLayers = 0;
  Double_t xtrkl[3];
  Double_t ptrkl[3];
  Double_t ptrkl2[3];
  UInt_t trklDet;
  UShort_t trklLayer;
  UInt_t stack_gtu;
  UShort_t stackInSector;

  for (UShort_t iLayer = 0; iLayer < 6; iLayer++){
    AliESDTrdTracklet* trkl = gtu_track->GetTracklet(iLayer);
    if (trkl){
      trklDet = trkl->GetDetector();
      trklLayer = TrdDetLyr(trklDet);
      stack_gtu = TrdDetLsi(trklDet);
      stackInSector = TrdDetSi(trklDet);

      // local coordinates of the outer end point of the tracklet
      xtrkl[0] = AliTRDgeometry::AnodePos();
      xtrkl[1] = trkl->GetLocalY();

      if(stackInSector == 2){ // corrected version by Felix Muecke
	xtrkl[2] = fTRDgeo->GetPadPlane(trklLayer, stackInSector)->GetRowPos(trkl->GetBinZ()) -
	  (fTRDgeo->GetPadPlane(trklLayer, stackInSector)->GetRowSize(trkl->GetBinZ()))/2. -
	  fTRDgeo->GetPadPlane(trklLayer, stackInSector)->GetRowPos(6);
      } else {
	xtrkl[2] = fTRDgeo->GetPadPlane(trklLayer, stackInSector)->GetRowPos(trkl->GetBinZ()) -
	  (fTRDgeo->GetPadPlane(trklLayer, stackInSector)->GetRowSize(trkl->GetBinZ()))/2. -
	  fTRDgeo->GetPadPlane(trklLayer, stackInSector)->GetRowPos(8);
      }

      // old draft version
      // xtrkl[2] = fTRDgeo->GetPadPlane(trklLayer, (trklDet/6) % 5)->GetRowPos(trkl->GetBinZ()) -
      //  	fTRDgeo->GetPadPlane(trklLayer, (trklDet/6) % 5)->GetRowSize(trkl->GetBinZ()) -
      //  	fTRDgeo->GetPadPlane(trklLayer, (trklDet/6) % 5)->GetRowPos(8);

      // transform to global coordinates
      TGeoHMatrix *matrix = fTRDgeo->GetClusterMatrix(trklDet);
      if (!matrix){
	if ((stack_gtu != 13*5+2) && (stack_gtu != 14*5+2) && (stack_gtu != 15*5+2))
	  AliError(Form("invalid TRD cluster matrix in EstimateTrackDistance for detector %i", trklDet));
	return -5;
      }
      matrix->LocalToMaster(xtrkl, ptrkl);
      fTRDgeo->RotateBack(gtu_track->GetSector() * 30, ptrkl, ptrkl2);  // ptrkl2 now contains the global position of the outer end point of the tracklet

      // calculate parameterization of plane representing the tracklets layer
      Double_t layer_zero_local[3] = {0., 0.,  0.};
      Double_t layer_zero_global[3], layer_zero_global2[3];

      matrix->LocalToMaster(layer_zero_local, layer_zero_global);
      fTRDgeo->RotateBack(trklDet, layer_zero_global, layer_zero_global2); // layer_zero_global2 points to chamber origin in global coords

      Double_t layer_ref_local[3] = {AliTRDgeometry::AnodePos(), 0.,  0.};
      Double_t layer_ref_global[3], layer_ref_global2[3];

      matrix->LocalToMaster(layer_ref_local, layer_ref_global);
      fTRDgeo->RotateBack(trklDet, layer_ref_global, layer_ref_global2); // layer_ref_global2 points to center anode pos within plane in global coords

      Double_t n0[3] = {layer_ref_global2[0]-layer_zero_global2[0],
			layer_ref_global2[1]-layer_zero_global2[1],
			layer_ref_global2[2]-layer_zero_global2[2]};

      Double_t n_len = TMath::Sqrt(n0[0]*n0[0] + n0[1]*n0[1] + n0[2]*n0[2]);
      if (n_len == 0.){ // This should never happen
	AliError("divison by zero in estimate_track_distance!");
	n_len = 1.;
      }
      Double_t n[3] = {n0[0]/n_len, n0[1]/n_len, n0[2]/n_len}; // normal vector of plane

      AliExternalTrackParam *outerTPC = new AliExternalTrackParam(*(esd_track->GetOuterParam()));
      Bool_t isects = TrackPlaneIntersect(outerTPC, layer_ref_global2, n, mag); // find intersection point between track and TRD layer
      delete outerTPC;
      outerTPC = NULL;

      if (isects == kFALSE){ // extrapolation fails, because track never reaches the TRD radius
	return -1;
      }

      Double_t m[2] = {ptrkl2[0] - layer_ref_global2[0], ptrkl2[1] - layer_ref_global2[1]};
      Double_t len_m = TMath::Sqrt(m[0]*m[0] + m[1]*m[1]);
      diff_y += len_m;
      diff_z += TMath::Abs(ptrkl2[2] - layer_ref_global2[2]);
      nLayers++;
    }
  }

  if (nLayers > 0){
    *ydist = diff_y / nLayers;
    *zdist = diff_z / nLayers;
    return 0;
  }
  else
    return -4;
}

Double_t AliTRDonlineTrackMatching::PtDiffRel(Double_t refPt, Double_t gtuPt){

  // return relative pt difference

  if (TMath::Abs(refPt) > 0.000001){
    return (gtuPt - refPt) / refPt;
  } else
    return 0.;
}


Double_t AliTRDonlineTrackMatching::RateTrackMatch(Double_t distY, Double_t distZ, Double_t rpt, Double_t gpt){

  // returns a match rating derived from Y and Z distance as well as pt difference

  // maximum limits for spatial distance
  if ((distY > 5.) || (distZ > 20.))
    return 0.;

  // same pt sign required
  if ((rpt * gpt) < 0.)
    return 0.;

  Double_t rating_distY = -0.1 * distY + 1.;
  Double_t rating_distZ = -0.025 * distZ + 1.;
  Double_t rating_ptDiff = 1. - TMath::Abs(PtDiffRel(rpt, gpt));

  if (rating_ptDiff <  0.)
    rating_ptDiff = 0.2;

  Double_t total = rating_distY * rating_distZ * rating_ptDiff;

#ifdef TRD_TM_DEBUG
  if (total > 1.){
    printf("<ERROR> track match rating exceeds limit of 1.0: %.3f", total);
  }
#endif

  return total;
}


void AliTRDonlineTrackMatching::SetEsdTrackDefaultCuts(const char* cutIdent) {

  if (strcmp(cutIdent, "strict") == 0){

#ifdef TRD_TM_DEBUG
    printf("AliTRDonlineTrackMatching -- default track cuts selected");
#endif

    fEsdTrackCutMinimal = kFALSE;
    fEsdTrackCutPrim = kFALSE;

    fEsdTrackCutMinTPCrows = 70;
    fEsdTrackCutRequireTPCrefit = kTRUE;
    fEsdTrackCutMinRatioRowsFindableClusters = 0.8;
    fEsdTrackCutMaxChi2TPCclusters = 4.;
    fEsdTrackVCutsChi2TPCconstrainedVsGlobal = 36.;

    fEsdTrackCutRequireITSrefit = kFALSE;
    fEsdTrackCutMaxChi2ITSclusters = 36.;

    fEsdTrackCutMaxDCAtoVertexXY = 1000.;
    fEsdTrackCutMaxDCAtoVertexZ = 2.;
    fEsdTrackCutsITSlayerMask = 0x0;

    fEsdTrackCutPtDCAOfs = 0.0105;
    fEsdTrackCutPtDCACoeff = 0.0350;
  } else if (strcmp(cutIdent, "minimal") == 0){

#ifdef TRD_TM_DEBUG
    printf("AliTRDonlineTrackMatching -- minimal track cuts selected\n");
#endif

    fEsdTrackCutMinimal = kFALSE;
    fEsdTrackCutPrim = kFALSE;

    fEsdTrackCutMinTPCrows = 70;
    fEsdTrackCutRequireTPCrefit = kTRUE;
    fEsdTrackCutMinRatioRowsFindableClusters = 0.;
    fEsdTrackCutMaxChi2TPCclusters = 100.;
    fEsdTrackVCutsChi2TPCconstrainedVsGlobal = 1000.;

    fEsdTrackCutRequireITSrefit = kFALSE;
    fEsdTrackCutMaxChi2ITSclusters = 0.;

    fEsdTrackCutMaxDCAtoVertexXY = 1000.;
    fEsdTrackCutMaxDCAtoVertexZ = 1000.;
    fEsdTrackCutsITSlayerMask = 0x0;
  } else
    AliErrorClass("invalid cut set");

}
