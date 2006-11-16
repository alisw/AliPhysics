#ifndef ALIEVENTTAGCUTS_H
#define ALIEVENTTAGCUTS_H
/*  See cxx source for full Copyright notice */


/* $Id$ */

//-------------------------------------------------------------------------
//                       Class AliEventTagCuts
//   This is the class for the cuts in event tags
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------

#include <TObject.h>

class AliEventTag;

//___________________________________________________________________________
class AliEventTagCuts : public TObject {
 public:
  AliEventTagCuts();
  ~AliEventTagCuts();
  void Reset();
  
 //____________________________________________________//
  void SetNParticipantsRange(Int_t low, Int_t high);
  void SetImpactParamRange(Float_t low, Float_t high);

  void SetPrimaryVertexXRange(Float_t low, Float_t high);
  void SetPrimaryVertexYRange(Float_t low, Float_t high);
  void SetPrimaryVertexZRange(Float_t low, Float_t high);
  void SetPrimaryVertexFlag(Int_t flag);
  void SetPrimaryVertexZErrorRange(Float_t low, Float_t high);

  void SetTriggerMask(ULong64_t trmask);
  void SetTriggerCluster(UChar_t trcluster);

  void SetZDCNeutron1Range(Float_t low, Float_t high);
  void SetZDCProton1Range(Float_t low, Float_t high);
  void SetZDCEMRange(Float_t low, Float_t high);
  void SetZDCNeutron2Range(Float_t low, Float_t high);
  void SetZDCProton2Range(Float_t low, Float_t high);
  void SetT0VertexZRange(Float_t low, Float_t high);

  void SetMultiplicityRange(Int_t low, Int_t high);
  void SetPosMultiplicityRange(Int_t low, Int_t high);
  void SetNegMultiplicityRange(Int_t low, Int_t high);
  void SetNeutrMultiplicityRange(Int_t low, Int_t high);
  void SetNV0sRange(Int_t low, Int_t high);
  void SetNCascadesRange(Int_t low, Int_t high);
  void SetNKinksRange(Int_t low, Int_t high);
 
  void SetNPMDTracksRange(Int_t low, Int_t high);
  void SetNFMDTracksRange(Int_t low, Int_t high);
  void SetNPHOSClustersRange(Int_t low, Int_t high);
  void SetNEMCALClustersRange(Int_t low, Int_t high);
  void SetNJetCandidatesRange(Int_t low, Int_t high);

  void SetTopJetEnergyMin(Float_t low);
  void SetTopNeutralEnergyMin(Float_t low);
  void SetNHardPhotonsRange(Int_t low, Int_t high);
  void SetNChargedAbove1GeVRange(Int_t low, Int_t high);
  void SetNChargedAbove3GeVRange(Int_t low, Int_t high);
  void SetNChargedAbove10GeVRange(Int_t low, Int_t high);
  void SetNMuonsAbove1GeVRange(Int_t low, Int_t high);
  void SetNMuonsAbove3GeVRange(Int_t low, Int_t high);
  void SetNMuonsAbove10GeVRange(Int_t low, Int_t high);
  void SetNElectronsAbove1GeVRange(Int_t low, Int_t high);
  void SetNElectronsAbove3GeVRange(Int_t low, Int_t high);
  void SetNElectronsAbove10GeVRange(Int_t low, Int_t high);
  void SetNElectronRange(Int_t low, Int_t high);
  void SetNMuonRange(Int_t low, Int_t high);
  void SetNPionRange(Int_t low, Int_t high);
  void SetNKaonRange(Int_t low, Int_t high);
  void SetNProtonRange(Int_t low, Int_t high);
  void SetNLambdaRange(Int_t low, Int_t high);
  void SetNPhotonRange(Int_t low, Int_t high);
  void SetNPi0Range(Int_t low, Int_t high);
  void SetNNeutronRange(Int_t low, Int_t high);
  void SetNKaon0Range(Int_t low, Int_t high); 
  void SetTotalPRange(Float_t low, Float_t high);
  void SetMeanPtRange(Float_t low, Float_t high);
  void SetTopPtMin(Float_t low);
  void SetTotalNeutralPRange(Float_t low, Float_t high);
  void SetMeanNeutralPtPRange(Float_t low, Float_t high);
  void SetTopNeutralPtMin(Float_t low);
  void SetEvPlaneAngleRange(Float_t low, Float_t high);
  void SetHBTRadiiRange(Float_t low, Float_t high);
 
  Bool_t IsAccepted(AliEventTag *EvTag) const;

  //____________________________________________________//
 private:
  Float_t fVxMin, fVxMax;  //Definition of the range of the Vx
  Bool_t fVxFlag;          //Shows whether this cut is used or not
  Float_t fVyMin, fVyMax;  //Definition of the range of the Vy
  Bool_t fVyFlag;          //Shows whether this cut is used or not
  Float_t fVzMin, fVzMax;  //Definition of the range of the Vz
  Bool_t fVzFlag;          //Shows whether this cut is used or not
  Int_t fParticipantsMin, fParticipantMax; //# participants range
  Bool_t fParticipantsFlag;//Shows whether this cut is used or not
  Float_t fImpactParamMin, fImpactParamMax; //Impact parameter range
  Bool_t fImpactParamFlag; //Shows whether this cut is used or not
  Int_t fPrimaryVertexFlag; //Primary vertex flag: 0->not found, 1->found
  Bool_t fPVFlag;          //Shows whether this cut is used or not

  Float_t fPrimaryVertexZErrorMin, fPrimaryVertexZErrorMax; //Range of the primary vertex z error
  Bool_t fPVzErrorFlag;          //Shows whether this cut is used or not
  ULong64_t fTriggerMask;  //trigger mask definition
  Bool_t fTriggerMaskFlag; //Shows whether this cut is used or not
  UChar_t fTriggerCluster;  //trigger cluster definition
  Bool_t fTriggerClusterFlag; //Shows whether this cut is used or not
  
  Float_t fZDCNeutron1EnergyMin, fZDCNeutron1EnergyMax; //ZDC min,max - neutron
  Bool_t fZDCNeutron1EnergyFlag;//Shows whether this cut is used or not
  Float_t fZDCProton1EnergyMin, fZDCProton1EnergyMax; //ZDC min,max - proton
  Bool_t fZDCProton1EnergyFlag;//Shows whether this cut is used or not
  Float_t fZDCNeutron2EnergyMin, fZDCNeutron2EnergyMax; //ZDC min,max - neutron
  Bool_t fZDCNeutron2EnergyFlag;//Shows whether this cut is used or not
  Float_t fZDCProton2EnergyMin, fZDCProton2EnergyMax; //ZDC min,max - proton
  Bool_t fZDCProton2EnergyFlag;//Shows whether this cut is used or not
  Float_t fZDCEMEnergyMin, fZDCEMEnergyMax; //ZDC min,max - em
  Bool_t fZDCEMEnergyFlag;//Shows whether this cut is used or not
  Float_t fT0VertexZMin, fT0VertexZMax; //T0 min, max
  Bool_t fT0VertexZFlag;//Shows whether this cut is used or not  
  Int_t fMultMin, fMultMax;  //Definition of the range of the multiplicity
  Bool_t fMultFlag;//Shows whether this cut is used or not
  Int_t fMultPosMin, fMultPosMax; //Positive tracks multiplicity range
  Bool_t fMultPosFlag;//Shows whether this cut is used or not
  Int_t fMultNegMin, fMultNegMax; //Negative tracks multiplicity range
  Bool_t fMultNegFlag;//Shows whether this cut is used or not
  Int_t fMultNeutrMin, fMultNeutrMax; //Neutral tracks multiplicity range
  Bool_t fMultNeutrFlag;//Shows whether this cut is used or not
  Int_t fV0sMin, fV0sMax; //Range of V0s
  Bool_t fV0sFlag;//Shows whether this cut is used or not
  Int_t fCascadesMin, fCascadesMax; //Range of cascades
  Bool_t fCascadesFlag;//Shows whether this cut is used or not
  Int_t fkinksMin, fkinksMax; //Range of kinks
  Bool_t fkinksFlag;//Shows whether this cut is used or not
  
  Int_t fPMDTracksMin, fPMDTracksMax; //Range of PMD tracks
  Bool_t fPMDTracksFlag;//Shows whether this cut is used or not
  Int_t fFMDTracksMin, fFMDTracksMax; //Range of FMD tracks
  Bool_t fFMDTracksFlag;//Shows whether this cut is used or not
  Int_t fPHOSClustersMin, fPHOSClustersMax; //Range of PHOS clusters
  Bool_t fPHOSClustersFlag;//Shows whether this cut is used or not
  Int_t fEMCALClustersMin, fEMCALClustersMax; //Range of EMCAL clusters
  Bool_t fEMCALClustersFlag;//Shows whether this cut is used or not
  Int_t fJetCandidatesMin, fJetCandidatesMax; //Range of jet candidates
  Bool_t fJetCandidatesFlag;//Shows whether this cut is used or not

  Float_t fMaxJetEnergy; //max jet energy info
  Bool_t fMaxJetEnergyFlag;//Shows whether this cut is used or not
  
  Int_t fNHardPhotonsCandidatesMin, fNHardPhotonsCandidatesMax; //Hard photons candidates
  Bool_t fNHardPhotonsCandidatesFlag;//Shows whether this cut is used or not
  Float_t fMaxNeutralEnergy; //max neutral energy info
  Bool_t fMaxNeutralFlag;//Shows whether this cut is used or not  
  Int_t fChargedAbove1GeVMin, fChargedAbove1GeVMax;//Definition of the range of the number of charged above 1GeV
  Bool_t fChargedAbove1GeVFlag;//Shows whether this cut is used or not
  Int_t fChargedAbove3GeVMin, fChargedAbove3GeVMax;//Definition of the range of the number of charged above 3GeV
  Bool_t fChargedAbove3GeVFlag;//Shows whether this cut is used or not
  Int_t fChargedAbove10GeVMin, fChargedAbove10GeVMax;//Definition of the range of the number of charged above 10GeV
  Bool_t fChargedAbove10GeVFlag;//Shows whether this cut is used or not
  Int_t fMuonsAbove1GeVMin, fMuonsAbove1GeVMax;//Definition of the range of the number of muons above 1GeV
  Bool_t fMuonsAbove1GeVFlag;//Shows whether this cut is used or not
  Int_t fMuonsAbove3GeVMin, fMuonsAbove3GeVMax;//Definition of the range of the number of muons above 3GeV
  Bool_t fMuonsAbove3GeVFlag;//Shows whether this cut is used or not
  Int_t fMuonsAbove10GeVMin, fMuonsAbove10GeVMax; //Definition of the range of the number of muons above 10GeV
  Bool_t fMuonsAbove10GeVFlag;//Shows whether this cut is used or not
  Int_t fElectronsAbove1GeVMin, fElectronsAbove1GeVMax;//Definition of the range of the number of electorns above 1GeV
  Bool_t fElectronsAbove1GeVFlag;//Shows whether this cut is used or not
  Int_t fElectronsAbove3GeVMin, fElectronsAbove3GeVMax;//Definition of the range of the number of electorns above 3GeV
  Bool_t fElectronsAbove3GeVFlag;//Shows whether this cut is used or not
  Int_t fElectronsAbove10GeVMin,fElectronsAbove10GeVMax;//Definition of the range of the number of electorns above 10GeV
  Bool_t fElectronsAbove10GeVFlag;//Shows whether this cut is used or not  
  Int_t fElectronsMin, fElectronsMax; //Number of electrons range
  Bool_t fElectronsFlag;//Shows whether this cut is used or not
  Int_t fMuonsMin, fMuonsMax;  //Number of muons range
  Bool_t fMuonsFlag;//Shows whether this cut is used or not
  Int_t fPionsMin, fPionsMax; //Number of pions range
  Bool_t fPionsFlag;//Shows whether this cut is used or not
  Int_t fKaonsMin, fKaonsMax; //Number of kaons range
  Bool_t fKaonsFlag;//Shows whether this cut is used or not
  Int_t fProtonsMin, fProtonsMax; //Number of protons range
  Bool_t fProtonsFlag;//Shows whether this cut is used or not
  Int_t fLambdasMin, fLambdasMax; //Number of lambdas range
  Bool_t fLambdasFlag;//Shows whether this cut is used or not
  Int_t fPhotonsMin, fPhotonsMax; //Number of photons range
  Bool_t fPhotonFlag;//Shows whether this cut is used or not
  Int_t fPi0sMin, fPi0sMax; //Number of Pi0s range
  Bool_t fPi0sFlag;//Shows whether this cut is used or not
  Int_t fNeutronsMin, fNeutronsMax; //Number of neutrons range
  Bool_t fNeutronsFlag;//Shows whether this cut is used or not
  Int_t fKaon0sMin, fKaon0sMax; //Number of K0s range
  Bool_t fKaon0sFlag;//Shows whether this cut is used or not  
  Float_t fTotalPMin, fTotalPMax; //Range of the sum of the momentum per event
  Bool_t fTotalPFlag;//Shows whether this cut is used or not
  Float_t fMeanPtMin, fMeanPtMax; //Range of mean Pt per event
  Bool_t fMeanPtFlag;//Shows whether this cut is used or not
  Float_t fMaxPt; //Max Pt for each event
  Bool_t fMaxPtFlag;//Shows whether this cut is used or not
  Float_t fTotalNeutralPMin, fTotalNeutralPMax; //Sum of the momentum per event for neutral
  Bool_t fTotalNeutralPFlag;//Shows whether this cut is used or not
  Float_t fMeanNeutralPtMin, fMeanNeutralPtMax; //Mean Pt per event for neutral
  Bool_t fMeanNeutralPtFlag;//Shows whether this cut is used or not
  Float_t fMaxNeutralPt; //Max Pt for each event for neutral
  Bool_t fMaxNeutralPtFlag;//Shows whether this cut is used or not
  Float_t fEventPlaneAngleMin, fEventPlaneAngleMax; //event plane info
  Bool_t fEventPlaneAngleFlag;//Shows whether this cut is used or not
  Float_t fHBTRadiiMin, fHBTRadiiMax; //HBT info
  Bool_t fHBTRadiiFlag;//Shows whether this cut is used or not

  ClassDef(AliEventTagCuts, 2)
};

#endif
