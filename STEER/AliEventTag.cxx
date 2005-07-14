/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//-----------------------------------------------------------------
//           Implementation of the EventTag class
//   This is the class to deal with the tags in the event level
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

#include "AliEventTag.h"

ClassImp(AliEventTag)

//______________________________________________________________________________
AliEventTag::AliEventTag()
{
  fAliceEventId = 0;
  fGUID = 0;

  fNumberOfParticipants = -10;
  fImpactParameter = -10.0;

  fPrimaryVertexX = -100.0;
  fPrimaryVertexY = -100.0;
  fPrimaryVertexZ = -100.0;

  fPrimaryVertexFlag = -1;
  fPrimaryVertexZError = -100.0;

  fTriggerInfo = -10;

  fZDCNeutronEnergy = -10.0;
  fZDCProtonEnergy = -10.0;
  fZDCEMEnergy = -10.0;

  fT0VertexZ = -10.0;

  fNumberOfTracks = -10;
  fNumberOfPositiveTracks = -10;
  fNumberOfNegativeTracks = -10;
  fNumberOfNeutralTracks = -10;
  
  fNumberOfV0s = -10;
  fNumberOfCascades = -10;
  fNumberOfKinks = -10;
  
  fNumberOfPMDTracks = -10;
  fNumberOfPHOSTracks = -10;
  fNumberOfEMCALTracks = -10;
  fNumberOfFMDTracks = -10;
  
  fNumberOfJetCandidates = -10;
  fNumberOfHardPhotonsCandidates = -10;

  fMaxJetEnergy = -100.0;
  fMaxNeutralEnergy = -100.0; 
  fNumberOfChargedAbovePtRange = -10;
  fNumberOfMuonsAbovePtRange = -10;
  fNumberOfElectronsAbovePtRange = -10;

  fNumberOfElectrons = -10;
  fNumberOfMuons = -10;
  fNumberOfPions = -10;
  fNumberOfKaons = -10;
  fNumberOfProtons = -10;
  fNumberOfLambdas = -10;

  fNumberOfPhotons = -10;
  fNumberOfPi0s = -10;
  fNumberOfNeutrons = -10;
  fNumberOfKaon0s = -10;

  fNumberOfJPsiCandidates = -10;
  fNumberOfPsiPrimeCandidates = -10;
  fNumberOfUpsilonCandidates = -10;
  fNumberOfUpsilonPrimeCandidates = -10;
  fNumberOfUpsilonDoublePrimeCandidates = -10;
  fNumberOfCharmParticleCandidates = -10;
  fNumberOfBeautyParticleCandidates = -10;
  
  fTotalP = -10.0;
  fMeanPt = -10.0;
  fMaxPt = -10.0;

  fTotalNeutralP = -10.0;
  fMeanNeutralPt = -10.0;
  fMaxNeutralPt = -10.0;

  fEventPlaneAngle = -10.0;
  fHBTRadii = -10.0;

}


//______________________________________________________________________________
AliEventTag::AliEventTag(AliEventTag *EvTag)
{
  // EventTag copy constructor
  SetEventId(EvTag->GetEventId());
  SetGUID(EvTag->GetGUID());
  
  SetNumOfParticipants(EvTag->GetNumOfParticipants());
  SetImpactParameter(EvTag->GetImpactParameter());
  
  SetVertexX(EvTag->GetVertexX());
  SetVertexY(EvTag->GetVertexY());
  SetVertexZ(EvTag->GetVertexZ());

  SetVertexFlag(EvTag->GetVertexFlag());
  SetVertexZError(EvTag->GetVertexZError());

  SetTrigger(EvTag->GetTrigger());
  
  SetZDCNeutronEnergy(EvTag->GetZDCNeutronEnergy());
  SetZDCProtonEnergy(EvTag->GetZDCProtonEnergy());
  SetZDCEMEnergy(EvTag->GetZDCEMEnergy());
  
  SetT0VertexZ(EvTag->GetT0VertexZ());
  
  SetNumOfTracks(EvTag->GetNumOfTracks());
  SetNumOfPosTracks(EvTag->GetNumOfPosTracks());
  SetNumOfNegTracks(EvTag->GetNumOfNegTracks());
  SetNumOfNeutrTracks(EvTag->GetNumOfNeutrTracks());
  
  SetNumOfV0s(EvTag->GetNumOfV0s());
  SetNumOfCascades(EvTag->GetNumOfCascades());
  SetNumOfKinks(EvTag->GetNumOfKinks());
  
  SetNumOfPMDTracks(EvTag->GetNumOfPMDTracks());
  SetNumOfPHOSTracks(EvTag->GetNumOfPHOSTracks());
  SetNumOfEMCALTracks(EvTag->GetNumOfEMCALTracks());
  SetNumOfFMDTracks(EvTag->GetNumOfFMDTracks());
  
  SetNumOfJetCandidates(EvTag->GetNumOfJetCandidates());
  SetNumOfHardPhotonsCandidates(EvTag->GetNumOfHardPhotonsCandidates());

  SetMaxJetEnergy(EvTag->GetMaxJetEnergy());
  SetMaxNeutralEnergy(EvTag->GetMaxNeutralEnergy());
  SetNumOfChargedAbovePtRange(EvTag->GetNumOfChargedAbovePtRange());
  SetNumOfMuonsAbovePtRange(EvTag->GetNumOfMuonsAbovePtRange());
  SetNumOfElectronsAbovePtRange(EvTag->GetNumOfElectronsAbovePtRange());



  SetNumOfJPsiCandidates(EvTag->GetNumOfJPsiCandidates());
  SetNumOfPsiPrimeCandidates(EvTag->GetNumOfPsiPrimeCandidates());
  SetNumOfUpsilonCandidates(EvTag->GetNumOfUpsilonCandidates());
  SetNumOfUpsilonPrimeCandidates(EvTag->GetNumOfUpsilonPrimeCandidates());
  SetNumOfUpsilonDoublePrimeCandidates(EvTag->GetNumOfUpsilonDoublePrimeCandidates());
  SetNumOfCharmCandidates(EvTag->GetNumOfCharmCandidates());
  SetNumOfBeautyCandidates(EvTag->GetNumOfBeautyCandidates());
  
  SetNumOfElectrons(EvTag->GetNumOfElectrons());
  SetNumOfMuons(EvTag->GetNumOfMuons());
  SetNumOfPions(EvTag->GetNumOfPions());
  SetNumOfKaons(EvTag->GetNumOfKaons());
  SetNumOfProtons(EvTag->GetNumOfProtons());
  SetNumOfLambdas(EvTag->GetNumOfLambdas());
 

  SetNumOfPhotons(EvTag->GetNumOfPhotons());
  SetNumOfPi0s(EvTag->GetNumOfPi0s());
  SetNumOfNeutrons(EvTag->GetNumOfNeutrons());
  SetNumOfKaon0s(EvTag->GetNumOfKaon0s());
  
  SetTotalMomentum(EvTag->GetTotalMomentum());
  SetMeanPt(EvTag->GetMeanPt());
  SetMaxPt(EvTag->GetMaxPt());

  SetNeutralTotalMomentum(EvTag->GetNeutralTotalMomentum());
  SetNeutralMeanPt(EvTag->GetNeutralMeanPt());
  SetNeutralMaxPt(EvTag->GetNeutralMaxPt());
  
  SetEventPlaneAngle(EvTag->GetEventPlaneAngle());
  SetHBTRadii(EvTag->GetHBTRadii());
}
//______________________________________________________________________________
AliEventTag::~AliEventTag()
{
}
