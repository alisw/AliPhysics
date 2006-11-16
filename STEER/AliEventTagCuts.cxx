/**************************************************************************
 * Author: Panos Christakoglou.                                           *
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
//           AliEventTagCuts class
//   This is the class to deal with the event tag level cuts
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

class AliLog;
class AliESD;

#include "AliEventTag.h"
#include "AliEventTagCuts.h"

ClassImp(AliEventTagCuts)


//___________________________________________________________________________
AliEventTagCuts::AliEventTagCuts() :
  TObject(),



  
  fVxMin(-1000.0), fVxMax(1000.0), 
  fVxFlag(kFALSE),
  fVyMin(-1000.0), fVyMax(1000.0),  
  fVyFlag(kFALSE),
  fVzMin(-1000.0), fVzMax(1000.0),
  fVzFlag(kFALSE),
  fParticipantsMin(-1), fParticipantMax(10000),
  fParticipantsFlag(kFALSE),
  fImpactParamMin(-1.0), fImpactParamMax(1000.0),
  fImpactParamFlag(kFALSE),
  fPrimaryVertexFlag(1),
  fPVFlag(kFALSE),

  fPrimaryVertexZErrorMin(-10000.), fPrimaryVertexZErrorMax(10000.),
  fPVzErrorFlag(kFALSE),
  fTriggerMask(0),
  fTriggerMaskFlag(kFALSE),
  fTriggerCluster(0),
  fTriggerClusterFlag(kFALSE),
 
  fZDCNeutron1EnergyMin(-1.0), fZDCNeutron1EnergyMax(100000.0),
  fZDCNeutron1EnergyFlag(kFALSE),
  fZDCProton1EnergyMin(-1.0), fZDCProton1EnergyMax(100000.0),
  fZDCProton1EnergyFlag(kFALSE),
  fZDCNeutron2EnergyMin(-1.0), fZDCNeutron2EnergyMax(100000.0),
  fZDCNeutron2EnergyFlag(kFALSE),
  fZDCProton2EnergyMin(-1.0), fZDCProton2EnergyMax(100000.0),
  fZDCProton2EnergyFlag(kFALSE),
  fZDCEMEnergyMin(-1.0), fZDCEMEnergyMax(100000.0),
  fZDCEMEnergyFlag(kFALSE),
  fT0VertexZMin(-10000.0), fT0VertexZMax(10000.0),  
  fT0VertexZFlag(kFALSE),
  fMultMin(0), fMultMax(100000),  
  fMultFlag(kFALSE),
  fMultPosMin(-1), fMultPosMax(100000),
  fMultPosFlag(kFALSE),
  fMultNegMin(-1), fMultNegMax(100000),
  fMultNegFlag(kFALSE),
  fMultNeutrMin(-1), fMultNeutrMax(100000),
  fMultNeutrFlag(kFALSE),
  fV0sMin(-1), fV0sMax(1000000),
  fV0sFlag(kFALSE),
  fCascadesMin(-1), fCascadesMax(100000),
  fCascadesFlag(kFALSE),
  fkinksMin(-1), fkinksMax(1000000),
  fkinksFlag(kFALSE),

  fPMDTracksMin(-1), fPMDTracksMax(100000),
  fPMDTracksFlag(kFALSE),
  fFMDTracksMin(-1), fFMDTracksMax(100000),
  fFMDTracksFlag(kFALSE),
  fPHOSClustersMin(-1), fPHOSClustersMax(100000),
  fPHOSClustersFlag(kFALSE),
  fEMCALClustersMin(-1), fEMCALClustersMax(100000),
  fEMCALClustersFlag(kFALSE),
  fJetCandidatesMin(-1), fJetCandidatesMax(100000),
  fJetCandidatesFlag(kFALSE),

  fMaxJetEnergy(-1.0), 
  fMaxJetEnergyFlag(kFALSE),
  fNHardPhotonsCandidatesMin(-1), fNHardPhotonsCandidatesMax(100000),
  fNHardPhotonsCandidatesFlag(kFALSE),
  fMaxNeutralEnergy(-1.0), 
  fMaxNeutralFlag(kFALSE),
  fChargedAbove1GeVMin(-1), fChargedAbove1GeVMax(100000),
  fChargedAbove1GeVFlag(kFALSE),
  fChargedAbove3GeVMin(-1), fChargedAbove3GeVMax(100000),
  fChargedAbove3GeVFlag(kFALSE),
  fChargedAbove10GeVMin(-1), fChargedAbove10GeVMax(100000),
  fChargedAbove10GeVFlag(kFALSE),
  fMuonsAbove1GeVMin(-1), fMuonsAbove1GeVMax(100000),
  fMuonsAbove1GeVFlag(kFALSE),
  fMuonsAbove3GeVMin(-1), fMuonsAbove3GeVMax(100000),
  fMuonsAbove3GeVFlag(kFALSE),
  fMuonsAbove10GeVMin(-1), fMuonsAbove10GeVMax(100000), 
  fMuonsAbove10GeVFlag(kFALSE),
  fElectronsAbove1GeVMin(-1), fElectronsAbove1GeVMax(100000),
  fElectronsAbove1GeVFlag(kFALSE),
  fElectronsAbove3GeVMin(-1), fElectronsAbove3GeVMax(100000),
  fElectronsAbove3GeVFlag(kFALSE),
  fElectronsAbove10GeVMin(-1), fElectronsAbove10GeVMax(100000),
  fElectronsAbove10GeVFlag(kFALSE),
  fElectronsMin(-1), fElectronsMax(100000),
  fElectronsFlag(kFALSE),
  fMuonsMin(-1), fMuonsMax(100000),
  fMuonsFlag(kFALSE),
  fPionsMin(-1), fPionsMax(100000),
  fPionsFlag(kFALSE),
  fKaonsMin(-1), fKaonsMax(100000),
  fKaonsFlag(kFALSE),
  fProtonsMin(-1), fProtonsMax(100000),
  fProtonsFlag(kFALSE),
  fLambdasMin(-1), fLambdasMax(100000),
  fLambdasFlag(kFALSE),
  fPhotonsMin(-1), fPhotonsMax(100000),
  fPhotonFlag(kFALSE),
  fPi0sMin(-1), fPi0sMax(100000), 
  fPi0sFlag(kFALSE),
  fNeutronsMin(-1), fNeutronsMax(100000), 
  fNeutronsFlag(kFALSE),
  fKaon0sMin(-1), fKaon0sMax(100000), 
  fKaon0sFlag(kFALSE),
  fTotalPMin(-1.0), fTotalPMax(1000000.0),
  fTotalPFlag(kFALSE),
  fMeanPtMin(-1.0), fMeanPtMax(100000.0),
  fMeanPtFlag(kFALSE),
  fMaxPt(-1.0),
  fMaxPtFlag(kFALSE),
  fTotalNeutralPMin(-1.0), fTotalNeutralPMax(1000000.0),   
  fTotalNeutralPFlag(kFALSE),
  fMeanNeutralPtMin(-1.0), fMeanNeutralPtMax(1000000.0), 
  fMeanNeutralPtFlag(kFALSE),
  fMaxNeutralPt(-1.0), 
  fMaxNeutralPtFlag(kFALSE),
  fEventPlaneAngleMin(-10000000.0), fEventPlaneAngleMax(10000000.0), 
  fEventPlaneAngleFlag(kFALSE),
  fHBTRadiiMin(-1.0), fHBTRadiiMax(100000.0), 
  fHBTRadiiFlag(kFALSE)
{
  //Default constructor which calls the Reset method.
  Reset();
}

//___________________________________________________________________________
AliEventTagCuts::~AliEventTagCuts() {  
  //Defaut destructor.
}

//___________________________________________________________________________
void AliEventTagCuts::Reset() {
  //Sets dummy values to every private member.
  fVxFlag = kFALSE;
  fVyFlag = kFALSE;
  fVzFlag = kFALSE;
  fParticipantsFlag = kFALSE;
  fImpactParamFlag = kFALSE;
  fPVFlag = kFALSE;

  fPVzErrorFlag = kFALSE;
  fTriggerMaskFlag = kFALSE;
  fTriggerClusterFlag = kFALSE;

  fZDCNeutron1EnergyFlag = kFALSE;
  fZDCProton1EnergyFlag = kFALSE;
  fZDCNeutron2EnergyFlag = kFALSE;
  fZDCProton2EnergyFlag = kFALSE;
  fZDCEMEnergyFlag = kFALSE;
  fT0VertexZFlag = kFALSE;
  fMultFlag = kFALSE;
  fMultPosFlag = kFALSE;
  fMultNegFlag = kFALSE;
  fMultNeutrFlag = kFALSE;
  fV0sFlag = kFALSE;
  fCascadesFlag = kFALSE;
  fkinksFlag = kFALSE;

  fPMDTracksFlag = kFALSE;
  fFMDTracksFlag = kFALSE;
  fPHOSClustersFlag = kFALSE;
  fEMCALClustersFlag = kFALSE;
  fJetCandidatesFlag = kFALSE;

  fMaxJetEnergyFlag = kFALSE;
  fNHardPhotonsCandidatesFlag = kFALSE;
  fMaxNeutralFlag = kFALSE;
  fChargedAbove1GeVFlag = kFALSE;
  fChargedAbove3GeVFlag = kFALSE;
  fChargedAbove10GeVFlag = kFALSE;
  fMuonsAbove1GeVFlag = kFALSE;
  fMuonsAbove3GeVFlag = kFALSE;
  fMuonsAbove10GeVFlag = kFALSE;
  fElectronsAbove1GeVFlag = kFALSE;
  fElectronsAbove3GeVFlag = kFALSE;
  fElectronsAbove10GeVFlag = kFALSE;
  fElectronsFlag = kFALSE;
  fMuonsFlag = kFALSE;
  fPionsFlag = kFALSE;
  fKaonsFlag = kFALSE;
  fProtonsFlag = kFALSE;
  fLambdasFlag = kFALSE;
  fPhotonFlag = kFALSE;
  fPi0sFlag = kFALSE;
  fNeutronsFlag = kFALSE;
  fKaon0sFlag = kFALSE;
  fTotalPFlag = kFALSE;
  fMeanPtFlag = kFALSE;
  fMaxPtFlag = kFALSE;
  fTotalNeutralPFlag = kFALSE;
  fMeanNeutralPtFlag = kFALSE;
  fMaxNeutralPtFlag = kFALSE;
  fEventPlaneAngleFlag = kFALSE;
  fHBTRadiiFlag = kFALSE;
  
  fVxMin = -1000.0; fVxMax = 1000.0; 
  fVyMin = -1000.0; fVyMax = 1000.0;  
  fVzMin = -1000.0; fVzMax = 1000.0;
  fParticipantsMin = -1; fParticipantMax = 10000;
  fImpactParamMin = -1.0; fImpactParamMax = 1000.0;
  fPrimaryVertexFlag = 1;

  fPrimaryVertexZErrorMin = -10000.; fPrimaryVertexZErrorMax = 10000.;
  fTriggerMask = 0;
  fTriggerCluster = 0;
 
  fZDCNeutron1EnergyMin = -1.0; fZDCNeutron1EnergyMax = 100000.0;
  fZDCProton1EnergyMin = -1.0; fZDCProton1EnergyMax = 100000.0;
  fZDCNeutron2EnergyMin = -1.0; fZDCNeutron2EnergyMax = 100000.0;
  fZDCProton2EnergyMin = -1.0; fZDCProton2EnergyMax = 100000.0;
  fZDCEMEnergyMin = -1.0; fZDCEMEnergyMax = 100000.0;
  fT0VertexZMin = -10000.0; fT0VertexZMax = 10000.0;  
  fMultMin = 0; fMultMax = 100000;  
  fMultPosMin = -1; fMultPosMax = 100000;
  fMultNegMin = -1; fMultNegMax = 100000;
  fMultNeutrMin = -1; fMultNeutrMax = 100000;
  fV0sMin = -1; fV0sMax = 1000000;
  fCascadesMin = -1; fCascadesMax = 100000;
  fkinksMin = -1; fkinksMax = 1000000;

  fPMDTracksMin = -1, fPMDTracksMax = 100000;
  fFMDTracksMin = -1, fFMDTracksMax = 100000;
  fPHOSClustersMin = -1, fPHOSClustersMax = 100000;
  fEMCALClustersMin = -1, fEMCALClustersMax = 100000;
  fJetCandidatesMin = -1, fJetCandidatesMax = 100000;

  fMaxJetEnergy = -1.0; 
  fNHardPhotonsCandidatesMin = -1; fNHardPhotonsCandidatesMax = 100000;
  fMaxNeutralEnergy = -1.0; 
  fChargedAbove1GeVMin = -1; fChargedAbove1GeVMax = 100000;
  fChargedAbove3GeVMin = -1; fChargedAbove3GeVMax = 100000;
  fChargedAbove10GeVMin = -1; fChargedAbove10GeVMax = 100000;
  fMuonsAbove1GeVMin = -1; fMuonsAbove1GeVMax = 100000;
  fMuonsAbove3GeVMin = -1; fMuonsAbove3GeVMax = 100000;
  fMuonsAbove10GeVMin = -1; fMuonsAbove10GeVMax = 100000; 
  fElectronsAbove1GeVMin = -1; fElectronsAbove1GeVMax = 100000;
  fElectronsAbove3GeVMin = -1; fElectronsAbove3GeVMax = 100000;
  fElectronsAbove10GeVMin = -1; fElectronsAbove10GeVMax = 100000;
  fElectronsMin = -1; fElectronsMax = 100000;
  fMuonsMin = -1; fMuonsMax = 100000;
  fPionsMin = -1; fPionsMax = 100000;
  fKaonsMin = -1; fKaonsMax = 100000;
  fProtonsMin = -1; fProtonsMax = 100000;
  fLambdasMin = -1; fLambdasMax = 100000;
  fPhotonsMin = -1; fPhotonsMax = 100000;
  fPi0sMin = -1; fPi0sMax = 100000; 
  fNeutronsMin = -1; fNeutronsMax = 100000; 
  fKaon0sMin = -1; fKaon0sMax = 100000; 
  fTotalPMin = -1.0; fTotalPMax = 1000000.0;
  fMeanPtMin = -1.0; fMeanPtMax = 100000.0;
  fMaxPt = -1.0; fTotalNeutralPMin = -1.0;
  fTotalNeutralPMax = 1000000.0;   
  fMeanNeutralPtMin = -1.0; fMeanNeutralPtMax = 1000000.0; 
  fMaxNeutralPt = -1.0; 
  fEventPlaneAngleMin = -10000000.0; fEventPlaneAngleMax = 10000000.0; 
  fHBTRadiiMin = -1.0; fHBTRadiiMax = 100000.0; 
}

//___________________________________________________________________________
void AliEventTagCuts::SetPrimaryVertexXRange(Float_t low, Float_t high) {
  //Sets the primary vertex x range 
  //and the corresponding flag to kTRUE if the cut is used.
  fVxMin = low;
  fVxMax = high; 
  fVxFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetPrimaryVertexYRange(Float_t low, Float_t high) {
  //Sets the primary vertex y range 
  //and the corresponding flag to kTRUE if the cut is used.
  fVyMin = low;
  fVyMax = high; 
  fVyFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetPrimaryVertexZRange(Float_t low, Float_t high) {
  //Sets the primary vertex z range 
  //and the corresponding flag to kTRUE if the cut is used.
  fVzMin = low;
  fVzMax = high; 
  fVzFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetPrimaryVertexZErrorRange(Float_t low, Float_t high) {
  //Sets the primary vertex z error range 
  //and the corresponding flag to kTRUE if the cut is used.
  fPrimaryVertexZErrorMin = low;
  fPrimaryVertexZErrorMax = high; 
  fPVzErrorFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetTriggerMask(ULong64_t trmask) {
  //Sets the trigger mask 
  //and the corresponding flag to kTRUE if the cut is used.
  fTriggerMask = trmask;
  fTriggerMaskFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetTriggerCluster(UChar_t trcluster) {
  //Sets the trigger cluster 
  //and the corresponding flag to kTRUE if the cut is used.
  fTriggerCluster = trcluster;
  fTriggerClusterFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetMultiplicityRange(Int_t low, Int_t high) {
  //Sets the primary multiplicity range 
  //and the corresponding flag to kTRUE if the cut is used.
  fMultMin = low;
  fMultMax = high;
  fMultFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNParticipantsRange(Int_t low, Int_t high) {
  //Sets the number of participants range 
  //and the corresponding flag to kTRUE if the cut is used.
  fParticipantsMin = low;
  fParticipantMax = high;
  fParticipantsFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetImpactParamRange(Float_t low, Float_t high) {
  //Sets the impact parameter range 
  //and the corresponding flag to kTRUE if the cut is used.
  fImpactParamMin = low;
  fImpactParamMax = high;
  fImpactParamFlag = kTRUE;
}
 

//___________________________________________________________________________
void AliEventTagCuts::SetPrimaryVertexFlag(Int_t flag) {
  //Sets the primary vertex flag cut 
  //and the corresponding flag to kTRUE if the cut is used.
  fPrimaryVertexFlag = flag;
  fPVFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetZDCNeutron1Range(Float_t low, Float_t high) {
  //Sets the ZDC's neutron energy range 
  //and the corresponding flag to kTRUE if the cut is used.
  fZDCNeutron1EnergyMin = low;
  fZDCNeutron1EnergyMax = high;
  fZDCNeutron1EnergyFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetZDCProton1Range(Float_t low, Float_t high) {
  //Sets the ZDC's proton energy range 
  //and the corresponding flag to kTRUE if the cut is used.
  fZDCProton1EnergyMin = low;
  fZDCProton1EnergyMax = high;
  fZDCProton1EnergyFlag = kTRUE;
}
//___________________________________________________________________________
void AliEventTagCuts::SetZDCNeutron2Range(Float_t low, Float_t high) {
  //Sets the ZDC's neutron energy range 
  //and the corresponding flag to kTRUE if the cut is used.
  fZDCNeutron2EnergyMin = low;
  fZDCNeutron2EnergyMax = high;
  fZDCNeutron2EnergyFlag = kTRUE;
}
//___________________________________________________________________________
void AliEventTagCuts::SetZDCProton2Range(Float_t low, Float_t high) {
  //Sets the ZDC's proton energy range 
  //and the corresponding flag to kTRUE if the cut is used.
  fZDCProton2EnergyMin = low;
  fZDCProton2EnergyMax = high;
  fZDCProton2EnergyFlag = kTRUE;
}
//___________________________________________________________________________
void AliEventTagCuts::SetZDCEMRange(Float_t low, Float_t high) {
  //Sets the ZDC's e/m energy range 
  //and the corresponding flag to kTRUE if the cut is used.
  fZDCEMEnergyMin = low;
  fZDCEMEnergyMax = high;
  fZDCEMEnergyFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetT0VertexZRange(Float_t low, Float_t high) {
  //Sets the T0's Vz range 
  //and the corresponding flag to kTRUE if the cut is used.
  fT0VertexZMin = low;
  fT0VertexZMax = high;
  fT0VertexZFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetPosMultiplicityRange(Int_t low, Int_t high) {
  //Sets the positive multiplicity range 
  //and the corresponding flag to kTRUE if the cut is used.
  fMultPosMin = low;
  fMultPosMax = high;
  fMultPosFlag = kTRUE;
}


//___________________________________________________________________________
void AliEventTagCuts::SetNegMultiplicityRange(Int_t low, Int_t high) {
  //Sets the negative multiplicity range 
  //and the corresponding flag to kTRUE if the cut is used.
  fMultNegMin = low;
  fMultNegMax = high;
  fMultNegFlag = kTRUE;
}


//___________________________________________________________________________
void AliEventTagCuts::SetNeutrMultiplicityRange(Int_t low, Int_t high) {
  //Sets the neutral particle multiplicity range 
  //and the corresponding flag to kTRUE if the cut is used.
  fMultNeutrMin = low;
  fMultNeutrMax = high;
  fMultNeutrFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNV0sRange(Int_t low, Int_t high) {
  //Sets the v0s multiplicity range 
  //and the corresponding flag to kTRUE if the cut is used.
  fV0sMin = low;
  fV0sMax = high;
  fV0sFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNCascadesRange(Int_t low, Int_t high) {
  //Sets the cascades multiplicity range 
  //and the corresponding flag to kTRUE if the cut is used.
  fCascadesMin = low;
  fCascadesMax = high;
  fCascadesFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNKinksRange(Int_t low, Int_t high) {
  //Sets the kinks multiplicity range 
  //and the corresponding flag to kTRUE if the cut is used.
  fkinksMin = low;
  fkinksMax = high;
  fkinksFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNPMDTracksRange(Int_t low, Int_t high) {
  //Sets the number of PMD tracks range 
  //and the corresponding flag to kTRUE if the cut is used.
  fPMDTracksMin = low;
  fPMDTracksMax = high;
  fPMDTracksFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNFMDTracksRange(Int_t low, Int_t high) {
  //Sets the number of FMD tracks range 
  //and the corresponding flag to kTRUE if the cut is used.
  fFMDTracksMin = low;
  fFMDTracksMax = high;
  fFMDTracksFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNPHOSClustersRange(Int_t low, Int_t high) {
  //Sets the number of PHOS clusters range 
  //and the corresponding flag to kTRUE if the cut is used.
  fPHOSClustersMin = low;
  fPHOSClustersMax = high;
  fPHOSClustersFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNEMCALClustersRange(Int_t low, Int_t high) {
  //Sets the number of EMCAL clusters range 
  //and the corresponding flag to kTRUE if the cut is used.
  fEMCALClustersMin = low;
  fEMCALClustersMax = high;
  fEMCALClustersFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNJetCandidatesRange(Int_t low, Int_t high) {
  //Sets the number of jet candidates range 
  //and the corresponding flag to kTRUE if the cut is used.
  fJetCandidatesMin = low;
  fJetCandidatesMax = high;
  fJetCandidatesFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetTopJetEnergyMin(Float_t low) {
  //Sets the lower limit of the maximum jet energy
  //and the corresponding flag to kTRUE if the cut is used.
  fMaxJetEnergy = low; 
  fMaxJetEnergyFlag = kTRUE;
}
//___________________________________________________________________________
void AliEventTagCuts::SetTopNeutralEnergyMin(Float_t low) {
  //Sets the lower limit of the maximum neutral jet energy
  //and the corresponding flag to kTRUE if the cut is used.
  fMaxNeutralEnergy = low; 
  fMaxNeutralFlag = kTRUE;
}
//___________________________________________________________________________
void AliEventTagCuts::SetNHardPhotonsRange(Int_t low, Int_t high) {
  //Sets the hard photons multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fNHardPhotonsCandidatesMin = low;
  fNHardPhotonsCandidatesMax = high;
  fNHardPhotonsCandidatesFlag = kTRUE;
} 

//___________________________________________________________________________
void AliEventTagCuts::SetNChargedAbove1GeVRange(Int_t low, Int_t high) {
  //Sets the number of charged above 1GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fChargedAbove1GeVMin = low;
  fChargedAbove1GeVMax = high;
  fChargedAbove1GeVFlag = kTRUE;
}

//___________________________________________________________________________
 void AliEventTagCuts::SetNChargedAbove3GeVRange(Int_t low, Int_t high) {
  //Sets the number of charged above 3GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fChargedAbove3GeVMin = low;
  fChargedAbove3GeVMax = high;
  fChargedAbove3GeVFlag = kTRUE;
}


//___________________________________________________________________________
void AliEventTagCuts::SetNChargedAbove10GeVRange(Int_t low, Int_t high) {
  //Sets the number of charged above 10GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fChargedAbove10GeVMin = low;
  fChargedAbove10GeVMax = high;
  fChargedAbove10GeVFlag = kTRUE;
}


//___________________________________________________________________________
void AliEventTagCuts::SetNMuonsAbove1GeVRange(Int_t low, Int_t high) {
  //Sets the number of muons above 1GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fMuonsAbove1GeVMin = low;
  fMuonsAbove1GeVMax = high;
  fMuonsAbove1GeVFlag = kTRUE;
}


//___________________________________________________________________________
void AliEventTagCuts::SetNMuonsAbove3GeVRange(Int_t low, Int_t high) {
  //Sets the number of muons above 3GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fMuonsAbove3GeVMin = low;
  fMuonsAbove3GeVMax = high;
  fMuonsAbove3GeVFlag = kTRUE;
} 

//___________________________________________________________________________
void AliEventTagCuts::SetNMuonsAbove10GeVRange(Int_t low, Int_t high) {
  //Sets the number of muons above 10GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fMuonsAbove10GeVMin = low;
  fMuonsAbove10GeVMax = high; 
  fMuonsAbove10GeVFlag = kTRUE;
}


//___________________________________________________________________________
void AliEventTagCuts::SetNElectronsAbove1GeVRange(Int_t low, Int_t high) {
  //Sets the number of electrons above 1GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fElectronsAbove1GeVMin = low;
  fElectronsAbove1GeVMax = high;
  fElectronsAbove1GeVFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNElectronsAbove3GeVRange(Int_t low, Int_t high) {
  //Sets the number of electrons above 3GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fElectronsAbove3GeVMin = low;
  fElectronsAbove3GeVMax = high;
  fElectronsAbove3GeVFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNElectronsAbove10GeVRange(Int_t low, Int_t high) {  
  //Sets the number of electrons above 10GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fElectronsAbove10GeVMin = low;
  fElectronsAbove10GeVMax = high;
  fElectronsAbove10GeVFlag = kTRUE;
}
//___________________________________________________________________________
void AliEventTagCuts::SetNElectronRange(Int_t low, Int_t high) {
  //Sets the electron multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fElectronsMin = low;
  fElectronsMax = high;
  fElectronsFlag = kTRUE;
}
//___________________________________________________________________________
void AliEventTagCuts::SetNMuonRange(Int_t low, Int_t high) {
  //Sets the muon multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fMuonsMin = low;
  fMuonsMax = high;
  fMuonsFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNPionRange(Int_t low, Int_t high) {
  //Sets the pion multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fPionsMin = low;
  fPionsMax = high;
  fPionsFlag = kTRUE;
} 

//___________________________________________________________________________
void AliEventTagCuts::SetNKaonRange(Int_t low, Int_t high) {
  //Sets the kaon multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fKaonsMin = low;
  fKaonsMax = high;
  fKaonsFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNProtonRange(Int_t low, Int_t high) {
  //Sets the proton multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fProtonsMin = low;
  fProtonsMax = high;
  fProtonsFlag = kTRUE;
} 

//___________________________________________________________________________
void AliEventTagCuts::SetNLambdaRange(Int_t low, Int_t high) {
  //Sets the lambda multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fLambdasMin = low;
  fLambdasMax = high;
  fLambdasFlag = kTRUE;
} 
//___________________________________________________________________________
void AliEventTagCuts::SetNPhotonRange(Int_t low, Int_t high) {
  //Sets the photon multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fPhotonsMin = low;
  fPhotonsMax = high;
  fPhotonFlag = kTRUE;
} 
//___________________________________________________________________________
void AliEventTagCuts::SetNPi0Range(Int_t low, Int_t high) {
  //Sets the pi0 multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fPi0sMin = low;
  fPi0sMax = high; 
  fPi0sFlag = kTRUE;
}  

//___________________________________________________________________________
void AliEventTagCuts::SetNNeutronRange(Int_t low, Int_t high) {
  //Sets the neutron multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fNeutronsMin = low;
  fNeutronsMax = high; 
  fNeutronsFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNKaon0Range(Int_t low, Int_t high) {  
  //Sets the K0s multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fKaon0sMin = low;
  fKaon0sMax = high; 
  fKaon0sFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetTotalPRange(Float_t low, Float_t high) {
  //Sets the total momentum range
  //and the corresponding flag to kTRUE if the cut is used.
  fTotalPMin = low;
  fTotalPMax = high;
  fTotalPFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetMeanPtRange(Float_t low, Float_t high) {
  //Sets the mean Pt range
  //and the corresponding flag to kTRUE if the cut is used.
  fMeanPtMin = low;
  fMeanPtMax = high;
  fMeanPtFlag = kTRUE;
}  

//___________________________________________________________________________
void AliEventTagCuts::SetTopPtMin(Float_t low) {
  //Sets the lower limit of the max Pt value
  //and the corresponding flag to kTRUE if the cut is used.
  fMaxPt = low; 
  fMaxPtFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetTotalNeutralPRange(Float_t low, Float_t high) {  
  //Sets the total momentum of neutral particles range
  //and the corresponding flag to kTRUE if the cut is used.
  fTotalNeutralPMin =low ;
  fTotalNeutralPMax = high;  
  fTotalNeutralPFlag = kTRUE;
}
//___________________________________________________________________________
void AliEventTagCuts::SetMeanNeutralPtPRange(Float_t low, Float_t high) {  
  //Sets the mean Pt of neutral particles range
  //and the corresponding flag to kTRUE if the cut is used.
  fMeanNeutralPtMin = low;
  fMeanNeutralPtMax = high; 
  fMeanNeutralPtFlag = kTRUE;
} 
//___________________________________________________________________________
void AliEventTagCuts::SetTopNeutralPtMin(Float_t low) {  
  //Sets the lower limit of the maximum Pt of neutral particles
  //and the corresponding flag to kTRUE if the cut is used.
  fMaxNeutralPt = low; 
  fMaxNeutralPtFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetEvPlaneAngleRange(Float_t low, Float_t high) {
  //Sets the event plane range
  //and the corresponding flag to kTRUE if the cut is used.
  fEventPlaneAngleMin = low;
  fEventPlaneAngleMax = high; 
  fEventPlaneAngleFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetHBTRadiiRange(Float_t low, Float_t high) {
  //Sets the HBT radii range
  //and the corresponding flag to kTRUE if the cut is used.
  fHBTRadiiMin = low;
  fHBTRadiiMax = high; 
  fHBTRadiiFlag = kTRUE;
}

//___________________________________________________________________________
Bool_t AliEventTagCuts::IsAccepted(AliEventTag *EvTag) const {
  //Returns true if the event is accepted otherwise false.
  if(fVzFlag)
    if((EvTag->GetVertexZ() < fVzMin) || (EvTag->GetVertexZ() > fVzMax))
      return kFALSE;
  
  if(fVyFlag)
    if((EvTag->GetVertexY() < fVyMin) || (EvTag->GetVertexY() > fVyMax))
      return kFALSE;
  
  if(fVxFlag)
    if((EvTag->GetVertexX() < fVxMin) || (EvTag->GetVertexX() > fVxMax))
      return kFALSE;
  
  if(fParticipantsFlag)
    if((EvTag->GetNumOfParticipants() < fParticipantsMin) || (EvTag->GetNumOfParticipants() > fParticipantMax))
      return kFALSE; 
  
  if(fImpactParamFlag)
    if((EvTag->GetImpactParameter() < fImpactParamMin) || (EvTag->GetImpactParameter() > fImpactParamMax))
      return kFALSE; 
  
  if(fPVFlag)
    if((EvTag->GetVertexFlag() != fPrimaryVertexFlag))
      return kFALSE; 
  
  if(fPVzErrorFlag)
    if((EvTag->GetVertexZError() < fPrimaryVertexZErrorMin) || (EvTag->GetVertexZError() > fPrimaryVertexZErrorMax))
      return kFALSE; 
  if(fTriggerMaskFlag)
    if((EvTag->GetTriggerMask() != fTriggerMask))
      return kFALSE; 
  if(fTriggerClusterFlag)
    if((EvTag->GetTriggerMask() != fTriggerMask))
      return kFALSE; 

  if(fZDCNeutron1EnergyFlag)
    if((EvTag->GetZDCNeutron1Energy() < fZDCNeutron1EnergyMin) || (EvTag->GetZDCNeutron1Energy() > fZDCNeutron1EnergyMax))
      return kFALSE; 
  
  if(fZDCProton1EnergyFlag)
    if((EvTag->GetZDCProton1Energy() < fZDCProton1EnergyMin) || (EvTag->GetZDCProton1Energy() > fZDCProton1EnergyMax))
      return kFALSE; 
  
  if(fZDCNeutron2EnergyFlag)
    if((EvTag->GetZDCNeutron2Energy() < fZDCNeutron2EnergyMin) || (EvTag->GetZDCNeutron2Energy() > fZDCNeutron2EnergyMax))
      return kFALSE; 
  
  if(fZDCProton2EnergyFlag)
    if((EvTag->GetZDCProton2Energy() < fZDCProton2EnergyMin) || (EvTag->GetZDCProton2Energy() > fZDCProton2EnergyMax))
      return kFALSE; 
  
  if(fZDCEMEnergyFlag)
    if((EvTag->GetZDCEMEnergy() < fZDCEMEnergyMin) || (EvTag->GetZDCEMEnergy() > fZDCEMEnergyMax))
      return kFALSE; 
  
  if(fT0VertexZFlag)
    if((EvTag->GetT0VertexZ() < fT0VertexZMin) || (EvTag->GetT0VertexZ() > fT0VertexZMax))
      return kFALSE; 
  
  if(fMultFlag)
    if((EvTag->GetNumOfTracks() < fMultMin) || (EvTag->GetNumOfTracks() > fMultMax))
      return kFALSE; 
  if(fMultPosFlag)
    if((EvTag->GetNumOfPosTracks() < fMultPosMin) || (EvTag->GetNumOfPosTracks() > fMultPosMax))
      return kFALSE; 
  
  if(fMultNegFlag)
    if((EvTag->GetNumOfNegTracks() < fMultNegMin) || (EvTag->GetNumOfNegTracks() > fMultNegMax))
      return kFALSE; 
  
  if(fMultNeutrFlag)
    if((EvTag->GetNumOfNeutrTracks() < fMultNeutrMin) || (EvTag->GetNumOfNeutrTracks() > fMultNeutrMax))
      return kFALSE; 
  
  if(fV0sFlag)
    if((EvTag->GetNumOfV0s() < fV0sMin) || (EvTag->GetNumOfV0s() > fV0sMax))
      return kFALSE; 
  
  if(fCascadesFlag)
    if((EvTag->GetNumOfCascades() < fCascadesMin) || (EvTag->GetNumOfCascades() > fCascadesMax))
      return kFALSE; 
  
  if(fkinksFlag)
    if((EvTag->GetNumOfKinks() < fkinksMin) || (EvTag->GetNumOfKinks() > fkinksMax))
      return kFALSE; 


  if(fPMDTracksFlag)
    if((EvTag->GetNumOfPMDTracks() < fPMDTracksMin) || (EvTag->GetNumOfPMDTracks() > fPMDTracksMax))
      return kFALSE; 
  if(fFMDTracksFlag)
    if((EvTag->GetNumOfFMDTracks() < fFMDTracksMin) || (EvTag->GetNumOfFMDTracks() > fFMDTracksMax))
      return kFALSE; 
  if(fPHOSClustersFlag)
    if((EvTag->GetNumOfPHOSClusters() < fPHOSClustersMin) || (EvTag->GetNumOfPHOSClusters() > fPHOSClustersMax))
      return kFALSE; 
  if(fEMCALClustersFlag)
    if((EvTag->GetNumOfEMCALClusters() < fEMCALClustersMin) || (EvTag->GetNumOfEMCALClusters() > fEMCALClustersMax))
      return kFALSE; 
  if(fJetCandidatesFlag)
    if((EvTag->GetNumOfJetCandidates() < fJetCandidatesMin) || (EvTag->GetNumOfJetCandidates() > fJetCandidatesMax))
      return kFALSE; 


  if(fMaxJetEnergyFlag)
    if((EvTag->GetMaxJetEnergy() < fMaxJetEnergy))
      return kFALSE; 
  
  if(fNHardPhotonsCandidatesFlag)
    if((EvTag->GetNumOfHardPhotonsCandidates() < fNHardPhotonsCandidatesMin) || (EvTag->GetNumOfHardPhotonsCandidates() > fNHardPhotonsCandidatesMax))
      return kFALSE; 
  
  if(fMaxNeutralFlag)
    if((EvTag->GetMaxNeutralEnergy() < fMaxNeutralEnergy))
      return kFALSE; 
  
  if(fChargedAbove1GeVFlag)
    if((EvTag->GetNumOfChargedAbove1GeV() < fChargedAbove1GeVMin) || (EvTag->GetNumOfChargedAbove1GeV() > fChargedAbove1GeVMax))
      return kFALSE; 
  
  if(fChargedAbove3GeVFlag)
    if((EvTag->GetNumOfChargedAbove3GeV() < fChargedAbove3GeVMin) || (EvTag->GetNumOfChargedAbove3GeV() > fChargedAbove3GeVMax))
      return kFALSE; 
  
  if(fChargedAbove10GeVFlag)
    if((EvTag->GetNumOfChargedAbove10GeV() < fChargedAbove10GeVMin) || (EvTag->GetNumOfChargedAbove10GeV() > fChargedAbove10GeVMax))
      return kFALSE; 
  
  if(fMuonsAbove1GeVFlag)
    if((EvTag->GetNumOfMuonsAbove1GeV() < fMuonsAbove1GeVMin) || (EvTag->GetNumOfMuonsAbove1GeV() > fMuonsAbove1GeVMax))
      return kFALSE; 
  
  if(fMuonsAbove3GeVFlag)
    if((EvTag->GetNumOfMuonsAbove3GeV() < fMuonsAbove3GeVMin) || (EvTag->GetNumOfMuonsAbove3GeV() > fMuonsAbove3GeVMax))
      return kFALSE; 
  
  if(fMuonsAbove10GeVFlag)
    if((EvTag->GetNumOfMuonsAbove10GeV() < fMuonsAbove10GeVMin) || (EvTag->GetNumOfMuonsAbove10GeV() > fMuonsAbove10GeVMax))
      return kFALSE; 
  
  if(fElectronsAbove1GeVFlag)
    if((EvTag->GetNumOfElectronsAbove1GeV()  < fElectronsAbove1GeVMin) || (EvTag->GetNumOfElectronsAbove1GeV()  > fElectronsAbove1GeVMax))
      return kFALSE; 
  
  if(fElectronsAbove3GeVFlag)
    if((EvTag->GetNumOfElectronsAbove3GeV() < fElectronsAbove3GeVMin) || (EvTag->GetNumOfElectronsAbove3GeV() > fElectronsAbove3GeVMax))
      return kFALSE; 
  
  if(fElectronsAbove10GeVFlag)
    if((EvTag->GetNumOfElectronsAbove10GeV() < fElectronsAbove10GeVMin) || (EvTag->GetNumOfElectronsAbove10GeV() > fElectronsAbove10GeVMax))
      return kFALSE; 
  
  if(fElectronsFlag)
    if((EvTag->GetNumOfElectrons() < fElectronsMin) || (EvTag->GetNumOfElectrons() > fElectronsMax))
      return kFALSE; 
  
  if(fMuonsFlag)
    if((EvTag->GetNumOfMuons() < fMuonsMin) || (EvTag->GetNumOfMuons() > fMuonsMax))
      return kFALSE; 
  
  if(fPionsFlag)
    if((EvTag->GetNumOfPions() < fPionsMin) || (EvTag->GetNumOfPions() > fPionsMax))
      return kFALSE; 
  
  if(fKaonsFlag)
    if((EvTag->GetNumOfKaons() < fKaonsMin) || (EvTag->GetNumOfKaons() > fKaonsMax))
      return kFALSE; 
  
  if(fProtonsFlag)
    if((EvTag->GetNumOfProtons() < fProtonsMin) || (EvTag->GetNumOfProtons() > fProtonsMax))
      return kFALSE; 
  
  if(fLambdasFlag)
    if((EvTag->GetNumOfLambdas() < fLambdasMin) || (EvTag->GetNumOfLambdas() > fLambdasMax))
      return kFALSE; 
  
  if(fPhotonFlag)
    if((EvTag->GetNumOfPhotons() < fPhotonsMin) || (EvTag->GetNumOfPhotons() > fPhotonsMax))
      return kFALSE; 
  
  if(fPi0sFlag)
    if((EvTag->GetNumOfPi0s() < fPi0sMin) || (EvTag->GetNumOfPi0s() > fPi0sMax))
      return kFALSE; 
  
  if(fNeutronsFlag)
    if((EvTag->GetNumOfNeutrons() < fNeutronsMin) || (EvTag->GetNumOfNeutrons() > fNeutronsMax))
      return kFALSE; 
  
  if(fKaon0sFlag)
    if((EvTag->GetNumOfKaon0s() < fKaon0sMin) || (EvTag->GetNumOfKaon0s() > fKaon0sMax))
      return kFALSE; 
  
  if(fTotalPFlag)
    if((EvTag->GetTotalMomentum() < fTotalPMin) || (EvTag->GetTotalMomentum() > fTotalPMax))
      return kFALSE; 
  
  if(fMeanPtFlag)
    if((EvTag->GetMeanPt() < fMeanPtMin) || (EvTag->GetMeanPt() > fMeanPtMax))
      return kFALSE; 
  
  if(fMaxPtFlag)
    if((EvTag->GetMaxPt() < fMaxPt))
      return kFALSE; 
  
  if(fTotalNeutralPFlag)
    if((EvTag->GetNeutralTotalMomentum() < fTotalNeutralPMin) || (EvTag->GetNeutralTotalMomentum() > fTotalNeutralPMax))
      return kFALSE; 
  
  if(fMeanNeutralPtFlag)
    if((EvTag->GetNeutralMeanPt() < fMeanNeutralPtMin) || (EvTag->GetNeutralMeanPt() >fMeanNeutralPtMax ))
      return kFALSE; 
  
  if(fMaxNeutralPtFlag)
    if((EvTag->GetNeutralMaxPt() < fMaxNeutralPt))
      return kFALSE; 
  
  if(fEventPlaneAngleFlag)
    if((EvTag->GetEventPlaneAngle() < fEventPlaneAngleMin) || (EvTag->GetEventPlaneAngle() > fEventPlaneAngleMax))
      return kFALSE; 
  
  if(fHBTRadiiFlag)
    if((EvTag->GetHBTRadii() < fHBTRadiiMin) || (EvTag->GetHBTRadii() > fHBTRadiiMax))
      return kFALSE; 
  
  return kTRUE;
}
