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
  AliEventTag::AliEventTag() : 
    TObject(),
    fPeriodNumber(0),
    fOrbitNumber(0),
    fBunchCrossNumber(0),
    fFiredTriggerClasses(),
    fEventType(0),
    fGUID(0),
    fPath(0),
    fsize(0),
    fmd5(0),
    fturl(0),
    fNumberOfParticipants(-10),
    fNumberOfParticipants2(-10),
    fImpactParameter(-10.0),
    fPrimaryVertexFlag(-1),
    fPrimaryVertexX(-100.0),
    fPrimaryVertexY(-100.0),
    fPrimaryVertexZ(-100.0),
    fPrimaryVertexZError(-100.0),
    fTriggerMask(0),
    fTriggerCluster(0),
    fZDCNeutron1Energy(-10.0),
    fZDCProton1Energy(-10.0),
    fZDCNeutron2Energy(-10.0),
    fZDCProton2Energy(-10.0),
    fT0VertexZ(-10.0),
    fNumberOfTracks(-10),
    fNumberOfPositiveTracks(-10),
    fNumberOfNegativeTracks(-10),
    fNumberOfNeutralTracks(-10),  
    fNumberOfV0s(-10),
    fNumberOfCascades(-10),
    fNumberOfKinks(-10),
    fNumberOfPMDTracks(-10),
    fNumberOfFMDTracks(-10),
    fNumberOfPHOSClusters(-10),
    fNumberOfEMCALClusters(-10),
    fNumberOfJetCandidates(-10),
    fMaxJetEnergy(-100.0),
    fNumberOfHardPhotonsCandidates(-10),
    fMaxNeutralEnergy(-100.0),
    fNumberOfChargedAbove1GeV(-10),
    fNumberOfChargedAbove3GeV(-10),
    fNumberOfChargedAbove10GeV(-10),
    fNumberOfMuonsAbove1GeV(-10),
    fNumberOfMuonsAbove3GeV(-10),
    fNumberOfMuonsAbove10GeV(-10),
    fNumberOfElectronsAbove1GeV(-10),
    fNumberOfElectronsAbove3GeV(-10),
    fNumberOfElectronsAbove10GeV(-10),
    fNumberOfElectrons(-10),
    fNumberOfFWMuons(-10),
    fNumberOfMuons(-10),
    fNumberOfPions(-10),
    fNumberOfKaons(-10),
    fNumberOfProtons(-10),
    fNumberOfLambdas(-10),
    fNumberOfPhotons(-10),
    fNumberOfPi0s(-10),
    fNumberOfNeutrons(-10),
    fNumberOfKaon0s(-10),
    fTotalP(-10.0),
    fMeanPt(-10.0),
    fMaxPt(-10.0),
    fEtaMaxPt(-13.0),
    fPhiMaxPt(+13.0),
    fTotalNeutralP(-10.0),
    fMeanNeutralPt(-10.0),
    fMaxNeutralPt(-10.0),
    fEventPlaneAngle(-10.0),
    fHBTRadii(-10.0),
    fNumberOfFiredChipsLayer1(0),
    fNumberOfFiredChipsLayer2(0),
    fNumberOfSPDTracklets(0),
    fMTotV0A(0),
    fMTotV0C(0),
    fNbPMV0A(0),
    fNbPMV0C(0)
{
  // AliEventTag default constructor
  for(Int_t i=0; i<2; i++)     fZDCEMEnergy[i] = -10.0;
}


//___________________________________________________________________________
AliEventTag::AliEventTag(const AliEventTag & evTag) :
  TObject(evTag),
  fPeriodNumber(evTag.fPeriodNumber),
  fOrbitNumber(evTag.fOrbitNumber),
  fBunchCrossNumber(evTag.fBunchCrossNumber),
  fFiredTriggerClasses(evTag.fFiredTriggerClasses),
  fEventType(evTag.fEventType),
  fGUID(evTag.fGUID),
  fPath(evTag.fPath),
  fsize(evTag.fsize),
  fmd5(evTag.fmd5),
  fturl(evTag.fturl),
  fNumberOfParticipants(evTag.fNumberOfParticipants),
  fNumberOfParticipants2(evTag.fNumberOfParticipants2),
  fImpactParameter(evTag.fImpactParameter),
  fPrimaryVertexFlag(evTag.fPrimaryVertexFlag),
  fPrimaryVertexX(evTag.fPrimaryVertexX),
  fPrimaryVertexY(evTag.fPrimaryVertexY),
  fPrimaryVertexZ(evTag.fPrimaryVertexZ),
  fPrimaryVertexZError(evTag.fPrimaryVertexZError),
  fTriggerMask(evTag.fTriggerMask),
  fTriggerCluster(evTag.fTriggerCluster),
  fZDCNeutron1Energy(evTag.fZDCNeutron1Energy),
  fZDCProton1Energy(evTag.fZDCProton1Energy),
  fZDCNeutron2Energy(evTag.fZDCNeutron2Energy),
  fZDCProton2Energy(evTag.fZDCProton2Energy),
  fT0VertexZ(evTag.fT0VertexZ),
  fNumberOfTracks(evTag.fNumberOfTracks),
  fNumberOfPositiveTracks(evTag.fNumberOfPositiveTracks),
  fNumberOfNegativeTracks(evTag.fNumberOfNegativeTracks),
  fNumberOfNeutralTracks(evTag.fNumberOfNeutralTracks),  
  fNumberOfV0s(evTag.fNumberOfV0s),
  fNumberOfCascades(evTag.fNumberOfCascades),
  fNumberOfKinks(evTag.fNumberOfKinks),
  fNumberOfPMDTracks(evTag.fNumberOfPMDTracks),
  fNumberOfFMDTracks(evTag.fNumberOfFMDTracks),
  fNumberOfPHOSClusters(evTag.fNumberOfPHOSClusters),
  fNumberOfEMCALClusters(evTag.fNumberOfEMCALClusters),
  fNumberOfJetCandidates(evTag.fNumberOfJetCandidates),
  fMaxJetEnergy(evTag.fMaxJetEnergy),
  fNumberOfHardPhotonsCandidates(evTag.fNumberOfHardPhotonsCandidates),
  fMaxNeutralEnergy(evTag.fMaxNeutralEnergy),
  fNumberOfChargedAbove1GeV(evTag.fNumberOfChargedAbove1GeV),
  fNumberOfChargedAbove3GeV(evTag.fNumberOfChargedAbove3GeV),
  fNumberOfChargedAbove10GeV(evTag.fNumberOfChargedAbove10GeV),
  fNumberOfMuonsAbove1GeV(evTag.fNumberOfMuonsAbove1GeV),
  fNumberOfMuonsAbove3GeV(evTag.fNumberOfMuonsAbove3GeV),
  fNumberOfMuonsAbove10GeV(evTag.fNumberOfMuonsAbove10GeV),
  fNumberOfElectronsAbove1GeV(evTag.fNumberOfElectronsAbove1GeV),
  fNumberOfElectronsAbove3GeV(evTag.fNumberOfElectronsAbove3GeV),
  fNumberOfElectronsAbove10GeV(evTag.fNumberOfElectronsAbove10GeV),
  fNumberOfElectrons(evTag.fNumberOfElectrons),
  fNumberOfFWMuons(evTag.fNumberOfFWMuons),
  fNumberOfMuons(evTag.fNumberOfMuons),
  fNumberOfPions(evTag.fNumberOfPions),
  fNumberOfKaons(evTag.fNumberOfKaons),
  fNumberOfProtons(evTag.fNumberOfProtons),
  fNumberOfLambdas(evTag.fNumberOfLambdas),
  fNumberOfPhotons(evTag.fNumberOfPhotons),
  fNumberOfPi0s(evTag.fNumberOfPi0s),
  fNumberOfNeutrons(evTag.fNumberOfNeutrons),
  fNumberOfKaon0s(evTag.fNumberOfKaon0s),
  fTotalP(evTag.fTotalP),
  fMeanPt(evTag.fMeanPt),
  fMaxPt(evTag.fMaxPt),
  fEtaMaxPt(evTag.fEtaMaxPt),
  fPhiMaxPt(evTag.fPhiMaxPt),
  fTotalNeutralP(evTag.fTotalNeutralP),
  fMeanNeutralPt(evTag.fMeanNeutralPt),
  fMaxNeutralPt(evTag.fMaxNeutralPt),
  fEventPlaneAngle(evTag.fEventPlaneAngle),
  fHBTRadii(evTag.fHBTRadii),
  fNumberOfFiredChipsLayer1(evTag.fNumberOfFiredChipsLayer1),
  fNumberOfFiredChipsLayer2(evTag.fNumberOfFiredChipsLayer2),
  fNumberOfSPDTracklets(evTag.fNumberOfSPDTracklets),
  fMTotV0A(evTag.fMTotV0A),
  fMTotV0C(evTag.fMTotV0C),
  fNbPMV0A(evTag.fNbPMV0A),
  fNbPMV0C(evTag.fNbPMV0C)
 {
  // EventTag copy constructor
  for(Int_t i=0; i<2; i++)     fZDCEMEnergy[i] = evTag.fZDCEMEnergy[i];
}

//___________________________________________________________________________
AliEventTag & AliEventTag::operator=(const AliEventTag &evTag) {
  // EventTag assignment operator
  if (this != &evTag) {
    TObject::operator=(evTag);
    
    SetPeriodNumber(evTag.GetPeriodNumber());
    SetOrbitNumber(evTag.GetOrbitNumber());
    SetBunchCrossNumber(evTag.GetBunchCrossNumber());
    SetFiredTriggerClasses(evTag.GetFiredTriggerClasses());
    SetEventType(evTag.GetEventType());
    SetGUID(evTag.GetGUID());
    SetPath(evTag.GetPath());
    SetMD5(evTag.GetMD5());
    SetTURL(evTag.GetTURL());
    SetSize(evTag.GetSize());
    SetNumOfParticipants(evTag.GetNumOfParticipants());
    SetImpactParameter(evTag.GetImpactParameter());
    SetVertexX(evTag.GetVertexX());
    SetVertexY(evTag.GetVertexY());
    SetVertexZ(evTag.GetVertexZ());
    SetVertexFlag(evTag.GetVertexFlag());
    SetVertexZError(evTag.GetVertexZError());
    SetTriggerMask(evTag.GetTriggerMask());
    SetTriggerCluster(evTag.GetTriggerCluster());
    SetZDCNeutron1Energy(evTag.GetZDCNeutron1Energy());
    SetZDCProton1Energy(evTag.GetZDCProton1Energy());
    SetZDCNeutron2Energy(evTag.GetZDCNeutron2Energy());
    SetZDCProton2Energy(evTag.GetZDCProton2Energy());
    SetZDCEMEnergy(evTag.GetZDCEMEnergy(0),evTag.GetZDCEMEnergy(1));
    SetT0VertexZ(evTag.GetT0VertexZ());
    SetNumOfTracks(evTag.GetNumOfTracks());
    SetNumOfPosTracks(evTag.GetNumOfPosTracks());
    SetNumOfNegTracks(evTag.GetNumOfNegTracks());
    SetNumOfNeutrTracks(evTag.GetNumOfNeutrTracks());
    SetNumOfV0s(evTag.GetNumOfV0s());
    SetNumOfCascades(evTag.GetNumOfCascades());
    SetNumOfKinks(evTag.GetNumOfKinks());
    SetNumOfPMDTracks(evTag.GetNumOfPMDTracks());
    SetNumOfFMDTracks(evTag.GetNumOfFMDTracks());
    SetNumOfPHOSClusters(evTag.GetNumOfPHOSClusters());
    SetNumOfEMCALClusters(evTag.GetNumOfEMCALClusters());
    SetNumOfJetCandidates(evTag.GetNumOfJetCandidates());
    SetNumOfHardPhotonsCandidates(evTag.GetNumOfHardPhotonsCandidates());
    SetMaxJetEnergy(evTag.GetMaxJetEnergy());
    SetMaxNeutralEnergy(evTag.GetMaxNeutralEnergy());
    SetNumOfChargedAbove1GeV(evTag.GetNumOfChargedAbove1GeV());
    SetNumOfChargedAbove3GeV(evTag.GetNumOfChargedAbove3GeV());
    SetNumOfChargedAbove10GeV(evTag.GetNumOfChargedAbove10GeV());
    SetNumOfMuonsAbove1GeV(evTag.GetNumOfMuonsAbove1GeV());
    SetNumOfMuonsAbove3GeV(evTag.GetNumOfMuonsAbove3GeV());
    SetNumOfMuonsAbove10GeV(evTag.GetNumOfMuonsAbove10GeV());
    SetNumOfElectronsAbove1GeV(evTag.GetNumOfElectronsAbove1GeV());
    SetNumOfElectronsAbove3GeV(evTag.GetNumOfElectronsAbove3GeV());
    SetNumOfElectronsAbove10GeV(evTag.GetNumOfElectronsAbove10GeV());
    SetNumOfElectrons(evTag.GetNumOfElectrons());
    SetNumOfFWMuons(evTag.GetNumOfFWMuons());
    SetNumOfMuons(evTag.GetNumOfMuons());
    SetNumOfPions(evTag.GetNumOfPions());
    SetNumOfKaons(evTag.GetNumOfKaons());
    SetNumOfProtons(evTag.GetNumOfProtons());
    SetNumOfLambdas(evTag.GetNumOfLambdas());
    SetNumOfPhotons(evTag.GetNumOfPhotons());
    SetNumOfPi0s(evTag.GetNumOfPi0s());
    SetNumOfNeutrons(evTag.GetNumOfNeutrons());
    SetNumOfKaon0s(evTag.GetNumOfKaon0s());
    SetTotalMomentum(evTag.GetTotalMomentum());
    SetMeanPt(evTag.GetMeanPt());
    SetMaxPt(evTag.GetMaxPt());
    SetEtaMaxPt(evTag.GetEtaMaxPt());
    SetPhiMaxPt(evTag.GetPhiMaxPt());
    SetNeutralTotalMomentum(evTag.GetNeutralTotalMomentum());
    SetNeutralMeanPt(evTag.GetNeutralMeanPt());
    SetNeutralMaxPt(evTag.GetNeutralMaxPt());
    SetEventPlaneAngle(evTag.GetEventPlaneAngle());
    SetHBTRadii(evTag.GetHBTRadii());
    SetNumberOfFiredChipsLayer1(evTag.GetNumberOfFiredChipsLayer1());
    SetNumberOfFiredChipsLayer2(evTag.GetNumberOfFiredChipsLayer2());
    SetNumberOfSPDTracklets(evTag.GetNumberOfSPDTracklets());
    SetMTotV0A(evTag.GetMTotV0A());
    SetMTotV0C(evTag.GetMTotV0C());
    SetNbPMV0A(evTag.GetNbPMV0A());
    SetNbPMV0C(evTag.GetNbPMV0C());
  }
  return *this;
}

//___________________________________________________________________________
AliEventTag::~AliEventTag() {
  // AliEventTag destructor
}
