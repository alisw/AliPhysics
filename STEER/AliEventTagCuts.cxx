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
  fNParticipantsMin(-1), fNParticipantsMax(10000),
  fNParticipantsFlag(kFALSE),
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
  fPosMultMin(-1), fPosMultMax(100000),
  fPosMultFlag(kFALSE),
  fNegMultMin(-1), fNegMultMax(100000),
  fNegMultFlag(kFALSE),
  fNeutrMultMin(-1), fNeutrMultMax(100000),
  fNeutrMultFlag(kFALSE),
  fNV0sMin(-1), fNV0sMax(1000000),
  fNV0sFlag(kFALSE),
  fNCascadesMin(-1), fNCascadesMax(100000),
  fNCascadesFlag(kFALSE),
  fNKinksMin(-1), fNKinksMax(1000000),
  fNKinksFlag(kFALSE),

  fNPMDTracksMin(-1), fNPMDTracksMax(100000),
  fNPMDTracksFlag(kFALSE),
  fNFMDTracksMin(-1), fNFMDTracksMax(100000),
  fNFMDTracksFlag(kFALSE),
  fNPHOSClustersMin(-1), fNPHOSClustersMax(100000),
  fNPHOSClustersFlag(kFALSE),
  fNEMCALClustersMin(-1), fNEMCALClustersMax(100000),
  fNEMCALClustersFlag(kFALSE),
  fNJetCandidatesMin(-1), fNJetCandidatesMax(100000),
  fNJetCandidatesFlag(kFALSE),

  fTopJetEnergyMin(-1.0), 
  fTopJetEnergyMinFlag(kFALSE),
  fNHardPhotonCandidatesMin(-1), fNHardPhotonCandidatesMax(100000),
  fNHardPhotonCandidatesFlag(kFALSE),
  fTopNeutralEnergyMin(-1.0), 
  fTopNeutralEnergyMinFlag(kFALSE),
  fNChargedAbove1GeVMin(-1), fNChargedAbove1GeVMax(100000),
  fNChargedAbove1GeVFlag(kFALSE),
  fNChargedAbove3GeVMin(-1), fNChargedAbove3GeVMax(100000),
  fNChargedAbove3GeVFlag(kFALSE),
  fNChargedAbove10GeVMin(-1), fNChargedAbove10GeVMax(100000),
  fNChargedAbove10GeVFlag(kFALSE),
  fNMuonsAbove1GeVMin(-1), fNMuonsAbove1GeVMax(100000),
  fNMuonsAbove1GeVFlag(kFALSE),
  fNMuonsAbove3GeVMin(-1), fNMuonsAbove3GeVMax(100000),
  fNMuonsAbove3GeVFlag(kFALSE),
  fNMuonsAbove10GeVMin(-1), fNMuonsAbove10GeVMax(100000), 
  fNMuonsAbove10GeVFlag(kFALSE),
  fNElectronsAbove1GeVMin(-1), fNElectronsAbove1GeVMax(100000),
  fNElectronsAbove1GeVFlag(kFALSE),
  fNElectronsAbove3GeVMin(-1), fNElectronsAbove3GeVMax(100000),
  fNElectronsAbove3GeVFlag(kFALSE),
  fNElectronsAbove10GeVMin(-1), fNElectronsAbove10GeVMax(100000),
  fNElectronsAbove10GeVFlag(kFALSE),
  fNElectronsMin(-1), fNElectronsMax(100000),
  fNElectronsFlag(kFALSE),
  fNMuonsMin(-1), fNMuonsMax(100000),
  fNMuonsFlag(kFALSE),
  fNPionsMin(-1), fNPionsMax(100000),
  fNPionsFlag(kFALSE),
  fNKaonsMin(-1), fNKaonsMax(100000),
  fNKaonsFlag(kFALSE),
  fNProtonsMin(-1), fNProtonsMax(100000),
  fNProtonsFlag(kFALSE),
  fNLambdasMin(-1), fNLambdasMax(100000),
  fNLambdasFlag(kFALSE),
  fNPhotonsMin(-1), fNPhotonsMax(100000),
  fNPhotonFlag(kFALSE),
  fNPi0sMin(-1), fNPi0sMax(100000), 
  fNPi0sFlag(kFALSE),
  fNNeutronsMin(-1), fNNeutronsMax(100000), 
  fNNeutronsFlag(kFALSE),
  fNKaon0sMin(-1), fNKaon0sMax(100000), 
  fNKaon0sFlag(kFALSE),
  fTotalPMin(-1.0), fTotalPMax(1000000.0),
  fTotalPFlag(kFALSE),
  fMeanPtMin(-1.0), fMeanPtMax(100000.0),
  fMeanPtFlag(kFALSE),
  fTopPtMin(-1.0),
  fTopPtMinFlag(kFALSE),
  fTotalNeutralPMin(-1.0), fTotalNeutralPMax(1000000.0),   
  fTotalNeutralPFlag(kFALSE),
  fMeanNeutralPtMin(-1.0), fMeanNeutralPtMax(1000000.0), 
  fMeanNeutralPtFlag(kFALSE),
  fTopNeutralPtMin(-1.0), 
  fTopNeutralPtMinFlag(kFALSE),
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
  fNParticipantsFlag = kFALSE;
  fImpactParamFlag = kFALSE;

  fVxFlag = kFALSE;
  fVyFlag = kFALSE;
  fVzFlag = kFALSE;
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
  fPosMultFlag = kFALSE;
  fNegMultFlag = kFALSE;
  fNeutrMultFlag = kFALSE;
  fNV0sFlag = kFALSE;
  fNCascadesFlag = kFALSE;
  fNKinksFlag = kFALSE;

  fNPMDTracksFlag = kFALSE;
  fNFMDTracksFlag = kFALSE;
  fNPHOSClustersFlag = kFALSE;
  fNEMCALClustersFlag = kFALSE;
  fNJetCandidatesFlag = kFALSE;

  fTopJetEnergyMinFlag = kFALSE;
  fNHardPhotonCandidatesFlag = kFALSE;
  fTopNeutralEnergyMinFlag = kFALSE;
  fNChargedAbove1GeVFlag = kFALSE;
  fNChargedAbove3GeVFlag = kFALSE;
  fNChargedAbove10GeVFlag = kFALSE;
  fNMuonsAbove1GeVFlag = kFALSE;
  fNMuonsAbove3GeVFlag = kFALSE;
  fNMuonsAbove10GeVFlag = kFALSE;
  fNElectronsAbove1GeVFlag = kFALSE;
  fNElectronsAbove3GeVFlag = kFALSE;
  fNElectronsAbove10GeVFlag = kFALSE;
  fNElectronsFlag = kFALSE;
  fNMuonsFlag = kFALSE;
  fNPionsFlag = kFALSE;
  fNKaonsFlag = kFALSE;
  fNProtonsFlag = kFALSE;
  fNLambdasFlag = kFALSE;
  fNPhotonFlag = kFALSE;
  fNPi0sFlag = kFALSE;
  fNNeutronsFlag = kFALSE;
  fNKaon0sFlag = kFALSE;
  fTotalPFlag = kFALSE;
  fMeanPtFlag = kFALSE;
  fTopPtMinFlag = kFALSE;
  fTotalNeutralPFlag = kFALSE;
  fMeanNeutralPtFlag = kFALSE;
  fTopNeutralPtMinFlag = kFALSE;
  fEventPlaneAngleFlag = kFALSE;
  fHBTRadiiFlag = kFALSE;
  
  fVxMin = -1000.0; fVxMax = 1000.0; 
  fVyMin = -1000.0; fVyMax = 1000.0;  
  fVzMin = -1000.0; fVzMax = 1000.0;
  fNParticipantsMin = -1; fNParticipantsMax = 10000;
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
  fPosMultMin = -1; fPosMultMax = 100000;
  fNegMultMin = -1; fNegMultMax = 100000;
  fNeutrMultMin = -1; fNeutrMultMax = 100000;
  fNV0sMin = -1; fNV0sMax = 1000000;
  fNCascadesMin = -1; fNCascadesMax = 100000;
  fNKinksMin = -1; fNKinksMax = 1000000;

  fNPMDTracksMin = -1, fNPMDTracksMax = 100000;
  fNFMDTracksMin = -1, fNFMDTracksMax = 100000;
  fNPHOSClustersMin = -1, fNPHOSClustersMax = 100000;
  fNEMCALClustersMin = -1, fNEMCALClustersMax = 100000;
  fNJetCandidatesMin = -1, fNJetCandidatesMax = 100000;

  fTopJetEnergyMin = -1.0; 
  fNHardPhotonCandidatesMin = -1; fNHardPhotonCandidatesMax = 100000;
  fTopNeutralEnergyMin = -1.0; 
  fNChargedAbove1GeVMin = -1; fNChargedAbove1GeVMax = 100000;
  fNChargedAbove3GeVMin = -1; fNChargedAbove3GeVMax = 100000;
  fNChargedAbove10GeVMin = -1; fNChargedAbove10GeVMax = 100000;
  fNMuonsAbove1GeVMin = -1; fNMuonsAbove1GeVMax = 100000;
  fNMuonsAbove3GeVMin = -1; fNMuonsAbove3GeVMax = 100000;
  fNMuonsAbove10GeVMin = -1; fNMuonsAbove10GeVMax = 100000; 
  fNElectronsAbove1GeVMin = -1; fNElectronsAbove1GeVMax = 100000;
  fNElectronsAbove3GeVMin = -1; fNElectronsAbove3GeVMax = 100000;
  fNElectronsAbove10GeVMin = -1; fNElectronsAbove10GeVMax = 100000;
  fNElectronsMin = -1; fNElectronsMax = 100000;
  fNMuonsMin = -1; fNMuonsMax = 100000;
  fNPionsMin = -1; fNPionsMax = 100000;
  fNKaonsMin = -1; fNKaonsMax = 100000;
  fNProtonsMin = -1; fNProtonsMax = 100000;
  fNLambdasMin = -1; fNLambdasMax = 100000;
  fNPhotonsMin = -1; fNPhotonsMax = 100000;
  fNPi0sMin = -1; fNPi0sMax = 100000; 
  fNNeutronsMin = -1; fNNeutronsMax = 100000; 
  fNKaon0sMin = -1; fNKaon0sMax = 100000; 
  fTotalPMin = -1.0; fTotalPMax = 1000000.0;
  fMeanPtMin = -1.0; fMeanPtMax = 100000.0;
  fTopPtMin = -1.0; fTotalNeutralPMin = -1.0;
  fTotalNeutralPMax = 1000000.0;   
  fMeanNeutralPtMin = -1.0; fMeanNeutralPtMax = 1000000.0; 
  fTopNeutralPtMin = -1.0; 
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
  fNParticipantsMin = low;
  fNParticipantsMax = high;
  fNParticipantsFlag = kTRUE;
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
  fPosMultMin = low;
  fPosMultMax = high;
  fPosMultFlag = kTRUE;
}


//___________________________________________________________________________
void AliEventTagCuts::SetNegMultiplicityRange(Int_t low, Int_t high) {
  //Sets the negative multiplicity range 
  //and the corresponding flag to kTRUE if the cut is used.
  fNegMultMin = low;
  fNegMultMax = high;
  fNegMultFlag = kTRUE;
}


//___________________________________________________________________________
void AliEventTagCuts::SetNeutrMultiplicityRange(Int_t low, Int_t high) {
  //Sets the neutral particle multiplicity range 
  //and the corresponding flag to kTRUE if the cut is used.
  fNeutrMultMin = low;
  fNeutrMultMax = high;
  fNeutrMultFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNV0sRange(Int_t low, Int_t high) {
  //Sets the v0s multiplicity range 
  //and the corresponding flag to kTRUE if the cut is used.
  fNV0sMin = low;
  fNV0sMax = high;
  fNV0sFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNCascadesRange(Int_t low, Int_t high) {
  //Sets the cascades multiplicity range 
  //and the corresponding flag to kTRUE if the cut is used.
  fNCascadesMin = low;
  fNCascadesMax = high;
  fNCascadesFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNKinksRange(Int_t low, Int_t high) {
  //Sets the kinks multiplicity range 
  //and the corresponding flag to kTRUE if the cut is used.
  fNKinksMin = low;
  fNKinksMax = high;
  fNKinksFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNPMDTracksRange(Int_t low, Int_t high) {
  //Sets the number of PMD tracks range 
  //and the corresponding flag to kTRUE if the cut is used.
  fNPMDTracksMin = low;
  fNPMDTracksMax = high;
  fNPMDTracksFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNFMDTracksRange(Int_t low, Int_t high) {
  //Sets the number of FMD tracks range 
  //and the corresponding flag to kTRUE if the cut is used.
  fNFMDTracksMin = low;
  fNFMDTracksMax = high;
  fNFMDTracksFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNPHOSClustersRange(Int_t low, Int_t high) {
  //Sets the number of PHOS clusters range 
  //and the corresponding flag to kTRUE if the cut is used.
  fNPHOSClustersMin = low;
  fNPHOSClustersMax = high;
  fNPHOSClustersFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNEMCALClustersRange(Int_t low, Int_t high) {
  //Sets the number of EMCAL clusters range 
  //and the corresponding flag to kTRUE if the cut is used.
  fNEMCALClustersMin = low;
  fNEMCALClustersMax = high;
  fNEMCALClustersFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNJetCandidatesRange(Int_t low, Int_t high) {
  //Sets the number of jet candidates range 
  //and the corresponding flag to kTRUE if the cut is used.
  fNJetCandidatesMin = low;
  fNJetCandidatesMax = high;
  fNJetCandidatesFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetTopJetEnergyMin(Float_t low) {
  //Sets the lower limit of the maximum jet energy
  //and the corresponding flag to kTRUE if the cut is used.
  fTopJetEnergyMin = low; 
  fTopJetEnergyMinFlag = kTRUE;
}
//___________________________________________________________________________
void AliEventTagCuts::SetTopNeutralEnergyMin(Float_t low) {
  //Sets the lower limit of the maximum neutral jet energy
  //and the corresponding flag to kTRUE if the cut is used.
  fTopNeutralEnergyMin = low; 
  fTopNeutralEnergyMinFlag = kTRUE;
}
//___________________________________________________________________________
void AliEventTagCuts::SetNHardPhotonsRange(Int_t low, Int_t high) {
  //Sets the hard photons multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fNHardPhotonCandidatesMin = low;
  fNHardPhotonCandidatesMax = high;
  fNHardPhotonCandidatesFlag = kTRUE;
} 

//___________________________________________________________________________
void AliEventTagCuts::SetNChargedAbove1GeVRange(Int_t low, Int_t high) {
  //Sets the number of charged above 1GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fNChargedAbove1GeVMin = low;
  fNChargedAbove1GeVMax = high;
  fNChargedAbove1GeVFlag = kTRUE;
}

//___________________________________________________________________________
 void AliEventTagCuts::SetNChargedAbove3GeVRange(Int_t low, Int_t high) {
  //Sets the number of charged above 3GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fNChargedAbove3GeVMin = low;
  fNChargedAbove3GeVMax = high;
  fNChargedAbove3GeVFlag = kTRUE;
}


//___________________________________________________________________________
void AliEventTagCuts::SetNChargedAbove10GeVRange(Int_t low, Int_t high) {
  //Sets the number of charged above 10GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fNChargedAbove10GeVMin = low;
  fNChargedAbove10GeVMax = high;
  fNChargedAbove10GeVFlag = kTRUE;
}


//___________________________________________________________________________
void AliEventTagCuts::SetNMuonsAbove1GeVRange(Int_t low, Int_t high) {
  //Sets the number of muons above 1GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fNMuonsAbove1GeVMin = low;
  fNMuonsAbove1GeVMax = high;
  fNMuonsAbove1GeVFlag = kTRUE;
}


//___________________________________________________________________________
void AliEventTagCuts::SetNMuonsAbove3GeVRange(Int_t low, Int_t high) {
  //Sets the number of muons above 3GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fNMuonsAbove3GeVMin = low;
  fNMuonsAbove3GeVMax = high;
  fNMuonsAbove3GeVFlag = kTRUE;
} 

//___________________________________________________________________________
void AliEventTagCuts::SetNMuonsAbove10GeVRange(Int_t low, Int_t high) {
  //Sets the number of muons above 10GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fNMuonsAbove10GeVMin = low;
  fNMuonsAbove10GeVMax = high; 
  fNMuonsAbove10GeVFlag = kTRUE;
}


//___________________________________________________________________________
void AliEventTagCuts::SetNElectronsAbove1GeVRange(Int_t low, Int_t high) {
  //Sets the number of electrons above 1GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fNElectronsAbove1GeVMin = low;
  fNElectronsAbove1GeVMax = high;
  fNElectronsAbove1GeVFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNElectronsAbove3GeVRange(Int_t low, Int_t high) {
  //Sets the number of electrons above 3GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fNElectronsAbove3GeVMin = low;
  fNElectronsAbove3GeVMax = high;
  fNElectronsAbove3GeVFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNElectronsAbove10GeVRange(Int_t low, Int_t high) {  
  //Sets the number of electrons above 10GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fNElectronsAbove10GeVMin = low;
  fNElectronsAbove10GeVMax = high;
  fNElectronsAbove10GeVFlag = kTRUE;
}
//___________________________________________________________________________
void AliEventTagCuts::SetNElectronRange(Int_t low, Int_t high) {
  //Sets the electron multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fNElectronsMin = low;
  fNElectronsMax = high;
  fNElectronsFlag = kTRUE;
}
//___________________________________________________________________________
void AliEventTagCuts::SetNMuonRange(Int_t low, Int_t high) {
  //Sets the muon multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fNMuonsMin = low;
  fNMuonsMax = high;
  fNMuonsFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNPionRange(Int_t low, Int_t high) {
  //Sets the pion multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fNPionsMin = low;
  fNPionsMax = high;
  fNPionsFlag = kTRUE;
} 

//___________________________________________________________________________
void AliEventTagCuts::SetNKaonRange(Int_t low, Int_t high) {
  //Sets the kaon multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fNKaonsMin = low;
  fNKaonsMax = high;
  fNKaonsFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNProtonRange(Int_t low, Int_t high) {
  //Sets the proton multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fNProtonsMin = low;
  fNProtonsMax = high;
  fNProtonsFlag = kTRUE;
} 

//___________________________________________________________________________
void AliEventTagCuts::SetNLambdaRange(Int_t low, Int_t high) {
  //Sets the lambda multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fNLambdasMin = low;
  fNLambdasMax = high;
  fNLambdasFlag = kTRUE;
} 
//___________________________________________________________________________
void AliEventTagCuts::SetNPhotonRange(Int_t low, Int_t high) {
  //Sets the photon multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fNPhotonsMin = low;
  fNPhotonsMax = high;
  fNPhotonFlag = kTRUE;
} 
//___________________________________________________________________________
void AliEventTagCuts::SetNPi0Range(Int_t low, Int_t high) {
  //Sets the pi0 multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fNPi0sMin = low;
  fNPi0sMax = high; 
  fNPi0sFlag = kTRUE;
}  

//___________________________________________________________________________
void AliEventTagCuts::SetNNeutronRange(Int_t low, Int_t high) {
  //Sets the neutron multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fNNeutronsMin = low;
  fNNeutronsMax = high; 
  fNNeutronsFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetNKaon0Range(Int_t low, Int_t high) {  
  //Sets the K0s multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fNKaon0sMin = low;
  fNKaon0sMax = high; 
  fNKaon0sFlag = kTRUE;
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
  fTopPtMin = low; 
  fTopPtMinFlag = kTRUE;
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
  fTopNeutralPtMin = low; 
  fTopNeutralPtMinFlag = kTRUE;
}

//___________________________________________________________________________
void AliEventTagCuts::SetEventPlaneAngleRange(Float_t low, Float_t high) {
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
  
  if(fNParticipantsFlag)
    if((EvTag->GetNumOfParticipants() < fNParticipantsMin) || (EvTag->GetNumOfParticipants() > fNParticipantsMax))
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
  
  if(fPosMultFlag)
    if((EvTag->GetNumOfPosTracks() < fPosMultMin) || (EvTag->GetNumOfPosTracks() > fPosMultMax))
      return kFALSE; 
  
  if(fNegMultFlag)
    if((EvTag->GetNumOfNegTracks() < fNegMultMin) || (EvTag->GetNumOfNegTracks() > fNegMultMax))
      return kFALSE; 
  
  if(fNeutrMultFlag)
    if((EvTag->GetNumOfNeutrTracks() < fNeutrMultMin) || (EvTag->GetNumOfNeutrTracks() > fNeutrMultMax))
      return kFALSE; 
  
  if(fNV0sFlag)
    if((EvTag->GetNumOfV0s() < fNV0sMin) || (EvTag->GetNumOfV0s() > fNV0sMax))
      return kFALSE; 
  
  if(fNCascadesFlag)
    if((EvTag->GetNumOfCascades() < fNCascadesMin) || (EvTag->GetNumOfCascades() > fNCascadesMax))
      return kFALSE; 
  
  if(fNKinksFlag)
    if((EvTag->GetNumOfKinks() < fNKinksMin) || (EvTag->GetNumOfKinks() > fNKinksMax))
      return kFALSE; 


  if(fNPMDTracksFlag)
    if((EvTag->GetNumOfPMDTracks() < fNPMDTracksMin) || (EvTag->GetNumOfPMDTracks() > fNPMDTracksMax))
      return kFALSE; 
  if(fNFMDTracksFlag)
    if((EvTag->GetNumOfFMDTracks() < fNFMDTracksMin) || (EvTag->GetNumOfFMDTracks() > fNFMDTracksMax))
      return kFALSE; 
  if(fNPHOSClustersFlag)
    if((EvTag->GetNumOfPHOSClusters() < fNPHOSClustersMin) || (EvTag->GetNumOfPHOSClusters() > fNPHOSClustersMax))
      return kFALSE; 
  if(fNEMCALClustersFlag)
    if((EvTag->GetNumOfEMCALClusters() < fNEMCALClustersMin) || (EvTag->GetNumOfEMCALClusters() > fNEMCALClustersMax))
      return kFALSE; 
  if(fNJetCandidatesFlag)
    if((EvTag->GetNumOfJetCandidates() < fNJetCandidatesMin) || (EvTag->GetNumOfJetCandidates() > fNJetCandidatesMax))
      return kFALSE; 


  if(fTopJetEnergyMinFlag)
    if((EvTag->GetMaxJetEnergy() < fTopJetEnergyMin))
      return kFALSE; 
  
  if(fNHardPhotonCandidatesFlag)
    if((EvTag->GetNumOfHardPhotonsCandidates() < fNHardPhotonCandidatesMin) || (EvTag->GetNumOfHardPhotonsCandidates() > fNHardPhotonCandidatesMax))
      return kFALSE; 
  
  if(fTopNeutralEnergyMinFlag)
    if((EvTag->GetMaxNeutralEnergy() < fTopNeutralEnergyMin))
      return kFALSE; 
  
  if(fNChargedAbove1GeVFlag)
    if((EvTag->GetNumOfChargedAbove1GeV() < fNChargedAbove1GeVMin) || (EvTag->GetNumOfChargedAbove1GeV() > fNChargedAbove1GeVMax))
      return kFALSE; 
  
  if(fNChargedAbove3GeVFlag)
    if((EvTag->GetNumOfChargedAbove3GeV() < fNChargedAbove3GeVMin) || (EvTag->GetNumOfChargedAbove3GeV() > fNChargedAbove3GeVMax))
      return kFALSE; 
  
  if(fNChargedAbove10GeVFlag)
    if((EvTag->GetNumOfChargedAbove10GeV() < fNChargedAbove10GeVMin) || (EvTag->GetNumOfChargedAbove10GeV() > fNChargedAbove10GeVMax))
      return kFALSE; 
  
  if(fNMuonsAbove1GeVFlag)
    if((EvTag->GetNumOfMuonsAbove1GeV() < fNMuonsAbove1GeVMin) || (EvTag->GetNumOfMuonsAbove1GeV() > fNMuonsAbove1GeVMax))
      return kFALSE; 
  
  if(fNMuonsAbove3GeVFlag)
    if((EvTag->GetNumOfMuonsAbove3GeV() < fNMuonsAbove3GeVMin) || (EvTag->GetNumOfMuonsAbove3GeV() > fNMuonsAbove3GeVMax))
      return kFALSE; 
  
  if(fNMuonsAbove10GeVFlag)
    if((EvTag->GetNumOfMuonsAbove10GeV() < fNMuonsAbove10GeVMin) || (EvTag->GetNumOfMuonsAbove10GeV() > fNMuonsAbove10GeVMax))
      return kFALSE; 
  
  if(fNElectronsAbove1GeVFlag)
    if((EvTag->GetNumOfElectronsAbove1GeV()  < fNElectronsAbove1GeVMin) || (EvTag->GetNumOfElectronsAbove1GeV()  > fNElectronsAbove1GeVMax))
      return kFALSE; 
  
  if(fNElectronsAbove3GeVFlag)
    if((EvTag->GetNumOfElectronsAbove3GeV() < fNElectronsAbove3GeVMin) || (EvTag->GetNumOfElectronsAbove3GeV() > fNElectronsAbove3GeVMax))
      return kFALSE; 
  
  if(fNElectronsAbove10GeVFlag)
    if((EvTag->GetNumOfElectronsAbove10GeV() < fNElectronsAbove10GeVMin) || (EvTag->GetNumOfElectronsAbove10GeV() > fNElectronsAbove10GeVMax))
      return kFALSE; 
  
  if(fNElectronsFlag)
    if((EvTag->GetNumOfElectrons() < fNElectronsMin) || (EvTag->GetNumOfElectrons() > fNElectronsMax))
      return kFALSE; 
  
  if(fNMuonsFlag)
    if((EvTag->GetNumOfMuons() < fNMuonsMin) || (EvTag->GetNumOfMuons() > fNMuonsMax))
      return kFALSE; 
  
  if(fNPionsFlag)
    if((EvTag->GetNumOfPions() < fNPionsMin) || (EvTag->GetNumOfPions() > fNPionsMax))
      return kFALSE; 
  
  if(fNKaonsFlag)
    if((EvTag->GetNumOfKaons() < fNKaonsMin) || (EvTag->GetNumOfKaons() > fNKaonsMax))
      return kFALSE; 
  
  if(fNProtonsFlag)
    if((EvTag->GetNumOfProtons() < fNProtonsMin) || (EvTag->GetNumOfProtons() > fNProtonsMax))
      return kFALSE; 
  
  if(fNLambdasFlag)
    if((EvTag->GetNumOfLambdas() < fNLambdasMin) || (EvTag->GetNumOfLambdas() > fNLambdasMax))
      return kFALSE; 
  
  if(fNPhotonFlag)
    if((EvTag->GetNumOfPhotons() < fNPhotonsMin) || (EvTag->GetNumOfPhotons() > fNPhotonsMax))
      return kFALSE; 
  
  if(fNPi0sFlag)
    if((EvTag->GetNumOfPi0s() < fNPi0sMin) || (EvTag->GetNumOfPi0s() > fNPi0sMax))
      return kFALSE; 
  
  if(fNNeutronsFlag)
    if((EvTag->GetNumOfNeutrons() < fNNeutronsMin) || (EvTag->GetNumOfNeutrons() > fNNeutronsMax))
      return kFALSE; 
  
  if(fNKaon0sFlag)
    if((EvTag->GetNumOfKaon0s() < fNKaon0sMin) || (EvTag->GetNumOfKaon0s() > fNKaon0sMax))
      return kFALSE; 
  
  if(fTotalPFlag)
    if((EvTag->GetTotalMomentum() < fTotalPMin) || (EvTag->GetTotalMomentum() > fTotalPMax))
      return kFALSE; 
  
  if(fMeanPtFlag)
    if((EvTag->GetMeanPt() < fMeanPtMin) || (EvTag->GetMeanPt() > fMeanPtMax))
      return kFALSE; 
  
  if(fTopPtMinFlag)
    if((EvTag->GetMaxPt() < fTopPtMin))
      return kFALSE; 
  
  if(fTotalNeutralPFlag)
    if((EvTag->GetNeutralTotalMomentum() < fTotalNeutralPMin) || (EvTag->GetNeutralTotalMomentum() > fTotalNeutralPMax))
      return kFALSE; 
  
  if(fMeanNeutralPtFlag)
    if((EvTag->GetNeutralMeanPt() < fMeanNeutralPtMin) || (EvTag->GetNeutralMeanPt() >fMeanNeutralPtMax ))
      return kFALSE; 
  
  if(fTopNeutralPtMinFlag)
    if((EvTag->GetNeutralMaxPt() < fTopNeutralPtMin))
      return kFALSE; 
  
  if(fEventPlaneAngleFlag)
    if((EvTag->GetEventPlaneAngle() < fEventPlaneAngleMin) || (EvTag->GetEventPlaneAngle() > fEventPlaneAngleMax))
      return kFALSE; 
  
  if(fHBTRadiiFlag)
    if((EvTag->GetHBTRadii() < fHBTRadiiMin) || (EvTag->GetHBTRadii() > fHBTRadiiMax))
      return kFALSE; 
  
  return kTRUE;
}
