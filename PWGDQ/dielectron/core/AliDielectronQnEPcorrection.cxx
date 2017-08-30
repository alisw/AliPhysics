/**
* @Author: Pascal Dillenseger <pascaldillenseger>
* @Date:   2016-11-11, 11:39:13
* @Email:  pdillens@cern.ch
 * @Last modified by:   pascaldillenseger
 * @Last modified time: 2017-08-07, 10:18:48
*/



/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////////////////
//          Class provides functions used for the eventplane
//          estimation with the 2016 est. Qn framework.
//          Can be used to configure and enable the auto-correlation
//          removal for the TPC eventplane.
// Authors:
//
//   Pascal Dillenseger <pdillens@cern.ch>
//
/////////////////////////////////////////////////////////////////////////


#include "AliDielectronQnEPcorrection.h"


ClassImp(AliDielectronQnEPcorrection)


AliDielectronQnEPcorrection::AliDielectronQnEPcorrection() :
  TNamed(),
  fPairCutActive{0},
  fPairCutValues{0.}
{
  //
  // Default costructor
  //
  fPairCutValues[kPairMass][kLow] = 2.92;
  fPairCutValues[kPairMass][kHigh] = 3.16;
}

//_______________________________________________________________________
AliDielectronQnEPcorrection::~AliDielectronQnEPcorrection()
{
  //
  // Destructor
  //
}

//_______________________________________________________________________
Bool_t AliDielectronQnEPcorrection::IsSelected(const AliDielectronPair *pair){
  // Apply pair cuts to check if the pair could create auto-correlation
  //  and therefore needs a eventplane Correction
  Bool_t activeCut = kFALSE;
  if(fPairCutActive[kPairMass]){
    if(pair->M() < fPairCutValues[kPairMass][kLow] || pair->M() > fPairCutValues[kPairMass][kHigh]) return kFALSE;
    else activeCut = kTRUE;
  }
  if(fPairCutActive[kPairPt]){
    if(pair->Pt() < fPairCutValues[kPairPt][kLow] || pair->Pt() > fPairCutValues[kPairPt][kHigh]) return kFALSE;
    else activeCut = kTRUE;
  }
  return activeCut;
}

//_______________________________________________________________________
Bool_t AliDielectronQnEPcorrection::IsSelected(const AliVTrack *track){
  if(track->IsA() == AliAODTrack::Class()) return IsSelectedAODtrack((AliAODTrack*) track);
  if(track->IsA() == AliESDtrack::Class()) return IsSelectedESDtrack((AliESDtrack*) track);
  AliFatal("No known aliroot track construct passed");
  return kFALSE;
}

//_______________________________________________________________________
Bool_t AliDielectronQnEPcorrection::IsSelectedESDtrack(const AliESDtrack *track){
  AliFatal("Code is not implemented for ESDs!!");
  return kFALSE;
}

//_______________________________________________________________________
Bool_t AliDielectronQnEPcorrection::IsSelectedAODtrack(const AliAODTrack *track){
  // Currently only the QnFramework track cuts are used, could be expanded if neccessary
  if(!(track->TestFilterBit(256) || track->TestFilterBit(512))) return kFALSE;
  if((track->Eta() < -.8) || (track->Eta() > .8)) return kFALSE;
  if((track->Pt() < .2) || (track->Pt() > 5.)) return kFALSE;
  return kTRUE;
}

//_______________________________________________________________________
Double_t AliDielectronQnEPcorrection::GetACcorrectedQnTPCEventplane(const AliDielectronPair* pair, TList* qnlist){
  //
  // Function to remove auto-correlations from the eventplanes estimated with the QnFramework est. 2016
  //  Checks the track array for particle pairs in the Jpsi mass window and removes their QnVector contributions
  //
  Double_t tpcEPangle = -999.;
  Int_t nRemovedTracks = 0;
  Float_t cQX=0., cQY=0.;

  // POI (particle of interest) rejection
  if(!pair) return tpcEPangle;
  // if(!AliDielectronQnEPcorrection::IsSelected(pair)) return tpcEPangle; // happens already in the VarManager
  const AliAODTrack *track1 = (AliAODTrack*) pair->GetFirstDaughterP();
  const AliAODTrack *track2 = (AliAODTrack*) pair->GetSecondDaughterP();

  Bool_t acceptedTrack1 = AliDielectronQnEPcorrection::IsSelected(track1);
  Bool_t acceptedTrack2 = AliDielectronQnEPcorrection::IsSelected(track2);

  if(!acceptedTrack1 && !acceptedTrack2) return tpcEPangle;
  //sum the contribution to the qVector of the tracks passing the EventPlanePOIPreFilter
  if(acceptedTrack1){
    cQX += TMath::Cos(2*track1->Phi());
    cQY += TMath::Sin(2*track1->Phi());
    nRemovedTracks++;
  }
  if(acceptedTrack2){
    cQX += TMath::Cos(2*track2->Phi());
    cQY += TMath::Sin(2*track2->Phi());
    nRemovedTracks++;
  }

  TVector2 qvec;
  AliQnCorrectionsQnVector *qnVectorPlain = AliDielectronQnEPcorrection::GetQnVectorFromList(qnlist,"TPC","plain","plain");
  AliQnCorrectionsQnVector *qnVectorRec = AliDielectronQnEPcorrection::GetQnVectorFromList(qnlist,"TPC","rec","rec");
  AliQnCorrectionsQnVector *qnVectorTwist = AliDielectronQnEPcorrection::GetQnVectorFromList(qnlist,"TPC","twist","twist");
  if(qnVectorPlain != NULL){
    Float_t qnPlainX = qnVectorPlain->Qx(2);
    Float_t qnPlainY = qnVectorPlain->Qy(2);
    if(qnVectorTwist != NULL){
      Float_t qnTwistX = qnVectorTwist->Qx(2);
      Float_t qnTwistY = qnVectorTwist->Qy(2);
      Float_t qnCorrX = qnPlainX - qnTwistX; // Correction from QnFramework on X value
      Float_t qnCorrY = qnPlainY - qnTwistY; // Correction from QnFramework on Y value
      Int_t M = qnVectorPlain->GetN(); // Number of used tracks
      Float_t qnWOautoCorrX = (qnPlainX*M - cQX)/(M - nRemovedTracks) -qnCorrX;
      Float_t qnWOautoCorrY = (qnPlainY*M - cQY)/(M - nRemovedTracks) -qnCorrY;
      qvec.Set(qnWOautoCorrX,qnWOautoCorrY);
      tpcEPangle = TVector2::Phi_mpi_pi(qvec.Phi())/2;
      return tpcEPangle;
    }
    if(qnVectorRec != NULL){
      Float_t qnRecX = qnVectorRec->Qx(2);
      Float_t qnRecY = qnVectorRec->Qy(2);
      Float_t qnCorrX = qnPlainX - qnRecX; // Correction from QnFramework on X value
      Float_t qnCorrY = qnPlainY - qnRecY; // Correction from QnFramework on Y value
      Int_t M = qnVectorPlain->GetN(); // Number of used tracks
      Float_t qnWOautoCorrX = (qnPlainX*M - cQX)/(M - nRemovedTracks) -qnCorrX;
      Float_t qnWOautoCorrY = (qnPlainY*M - cQY)/(M - nRemovedTracks) -qnCorrY;
      qvec.Set(qnWOautoCorrX,qnWOautoCorrY);
      tpcEPangle = TVector2::Phi_mpi_pi(qvec.Phi())/2;
      return tpcEPangle;
    }
  }
  return -tpcEPangle;
}




//________________________________________________________________
AliQnCorrectionsQnVector* AliDielectronQnEPcorrection::GetQnVectorFromList(
    TList *list,
    const char *subdetector,
    const char *expectedstep,
    const char *altstep){

  AliQnCorrectionsQnVector *theQnVector = NULL;

  TList *pQvecList = dynamic_cast<TList*> (list->FindObject(subdetector));
  if (pQvecList != NULL) {
    /* the detector is present */
    if (TString(expectedstep).EqualTo("latest"))
      theQnVector = (AliQnCorrectionsQnVector*) pQvecList->First();
    else
      theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(expectedstep);

    if (theQnVector == NULL) {
      /* the Qn vector for the expected step was not there */
      if (TString(altstep).EqualTo("latest"))
        theQnVector = (AliQnCorrectionsQnVector*) pQvecList->First();
      else
        theQnVector = (AliQnCorrectionsQnVector*) pQvecList->FindObject(altstep);
    }
  }
  if (theQnVector != NULL) {
    /* check the Qn vector quality */
    if (!(theQnVector->IsGoodQuality()) || !(theQnVector->GetN() != 0))
      /* not good quality, discarded */
      theQnVector = NULL;
  }
  return theQnVector;
}
