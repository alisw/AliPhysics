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
  void SetEventPlaneAngleRange(Float_t low, Float_t high);
  void SetHBTRadiiRange(Float_t low, Float_t high);
 
  Bool_t IsAccepted(AliEventTag *EvTag) const;

  //____________________________________________________//
 private:
  Int_t fNParticipantsMin, fNParticipantsMax;                 // # participants range
  Bool_t fNParticipantsFlag;                                  // Shows whether this cut is used or not
  Float_t fImpactParamMin, fImpactParamMax;                   // Impact parameter range
  Bool_t fImpactParamFlag;                                    // Shows whether this cut is used or not

  Float_t fVxMin, fVxMax;                                     // Definition of the range of the Vx
  Bool_t fVxFlag;                                             // Shows whether this cut is used or not
  Float_t fVyMin, fVyMax;                                     // Definition of the range of the Vy
  Bool_t fVyFlag;                                             // Shows whether this cut is used or not
  Float_t fVzMin, fVzMax;                                     // Definition of the range of the Vz
  Bool_t fVzFlag;                                             // Shows whether this cut is used or not
  Int_t fPrimaryVertexFlag;                                   // Primary vertex flag: 0->not found, 1->found
  Bool_t fPVFlag;                                             // Shows whether this cut is used or not
  Float_t fPrimaryVertexZErrorMin, fPrimaryVertexZErrorMax;   // Range of the primary vertex z error
  Bool_t fPVzErrorFlag;                                       // Shows whether this cut is used or not

  ULong64_t fTriggerMask;                                     // trigger mask definition
  Bool_t fTriggerMaskFlag;                                    // Shows whether this cut is used or not
  UChar_t fTriggerCluster;                                    // trigger cluster definition
  Bool_t fTriggerClusterFlag;                                 // Shows whether this cut is used or not
  
  Float_t fZDCNeutron1EnergyMin, fZDCNeutron1EnergyMax;       // ZDC min,max - neutron
  Bool_t fZDCNeutron1EnergyFlag;                              // Shows whether this cut is used or not
  Float_t fZDCProton1EnergyMin, fZDCProton1EnergyMax;         // ZDC min,max - proton
  Bool_t fZDCProton1EnergyFlag;                               // Shows whether this cut is used or not
  Float_t fZDCNeutron2EnergyMin, fZDCNeutron2EnergyMax;       // ZDC min,max - neutron
  Bool_t fZDCNeutron2EnergyFlag;                              // Shows whether this cut is used or not
  Float_t fZDCProton2EnergyMin, fZDCProton2EnergyMax;         // ZDC min,max - proton
  Bool_t fZDCProton2EnergyFlag;                               // Shows whether this cut is used or not
  Float_t fZDCEMEnergyMin, fZDCEMEnergyMax;                   // ZDC min,max - em
  Bool_t fZDCEMEnergyFlag;                                    // Shows whether this cut is used or not
  Float_t fT0VertexZMin, fT0VertexZMax;                       // T0 min, max
  Bool_t fT0VertexZFlag;                                      // Shows whether this cut is used or not  

  Int_t fMultMin, fMultMax;                                   // Definition of the range of the multiplicity
  Bool_t fMultFlag;                                           // Shows whether this cut is used or not
  Int_t fPosMultMin, fPosMultMax;                             // Positive tracks multiplicity range
  Bool_t fPosMultFlag;                                        // Shows whether this cut is used or not
  Int_t fNegMultMin, fNegMultMax;                             // Negative tracks multiplicity range
  Bool_t fNegMultFlag;                                        // Shows whether this cut is used or not
  Int_t fNeutrMultMin, fNeutrMultMax;                         // Neutral tracks multiplicity range
  Bool_t fNeutrMultFlag;                                      // Shows whether this cut is used or not
  Int_t fNV0sMin, fNV0sMax;                                   // Range of # of V0s
  Bool_t fNV0sFlag;                                           // Shows whether this cut is used or not
  Int_t fNCascadesMin, fNCascadesMax;                         // Range of # of cascades
  Bool_t fNCascadesFlag;                                      // Shows whether this cut is used or not
  Int_t fNKinksMin, fNKinksMax;                               // Range of # of kinks
  Bool_t fNKinksFlag;                                         // Shows whether this cut is used or not
  
  Int_t fNPMDTracksMin, fNPMDTracksMax;                       // Range of # of PMD tracks
  Bool_t fNPMDTracksFlag;                                     // Shows whether this cut is used or not
  Int_t fNFMDTracksMin, fNFMDTracksMax;                       // Range of # of FMD tracks
  Bool_t fNFMDTracksFlag;                                     // Shows whether this cut is used or not
  Int_t fNPHOSClustersMin, fNPHOSClustersMax;                 // Range of # of PHOS clusters
  Bool_t fNPHOSClustersFlag;                                  // Shows whether this cut is used or not
  Int_t fNEMCALClustersMin, fNEMCALClustersMax;               // Range of # of EMCAL clusters
  Bool_t fNEMCALClustersFlag;                                 // Shows whether this cut is used or not
  Int_t fNJetCandidatesMin, fNJetCandidatesMax;               // Range of # of jet candidates
  Bool_t fNJetCandidatesFlag;                                 // Shows whether this cut is used or not

  Float_t fTopJetEnergyMin;                                   // top jet energy minimum value
  Bool_t fTopJetEnergyMinFlag;                                // Shows whether this cut is used or not
  Float_t fTopNeutralEnergyMin;                               // top neutral energy minimum value
  Bool_t fTopNeutralEnergyMinFlag;                            // Shows whether this cut is used or not  
  
  Int_t fNHardPhotonCandidatesMin, fNHardPhotonCandidatesMax; // # of hard photons candidates
  Bool_t fNHardPhotonCandidatesFlag;                          // Shows whether this cut is used or not
  Int_t fNChargedAbove1GeVMin, fNChargedAbove1GeVMax;         // Definition of the range of the # of charged above 1GeV
  Bool_t fNChargedAbove1GeVFlag;                              // Shows whether this cut is used or not
  Int_t fNChargedAbove3GeVMin, fNChargedAbove3GeVMax;         // Definition of the range of the # of charged above 3GeV
  Bool_t fNChargedAbove3GeVFlag;                              // Shows whether this cut is used or not
  Int_t fNChargedAbove10GeVMin, fNChargedAbove10GeVMax;       // Definition of the range of the # of charged above 10GeV
  Bool_t fNChargedAbove10GeVFlag;                             // Shows whether this cut is used or not
  Int_t fNMuonsAbove1GeVMin, fNMuonsAbove1GeVMax;             // Definition of the range of the # of muons above 1GeV
  Bool_t fNMuonsAbove1GeVFlag;                                // Shows whether this cut is used or not
  Int_t fNMuonsAbove3GeVMin, fNMuonsAbove3GeVMax;             // Definition of the range of the # of muons above 3GeV
  Bool_t fNMuonsAbove3GeVFlag;                                // Shows whether this cut is used or not
  Int_t fNMuonsAbove10GeVMin, fNMuonsAbove10GeVMax;           // Definition of the range of the # of muons above 10GeV
  Bool_t fNMuonsAbove10GeVFlag;                               // Shows whether this cut is used or not
  Int_t fNElectronsAbove1GeVMin, fNElectronsAbove1GeVMax;     // Definition of the range of the # of electorns above 1GeV
  Bool_t fNElectronsAbove1GeVFlag;                            // Shows whether this cut is used or not
  Int_t fNElectronsAbove3GeVMin, fNElectronsAbove3GeVMax;     // Definition of the range of the # of electorns above 3GeV
  Bool_t fNElectronsAbove3GeVFlag;                            // Shows whether this cut is used or not
  Int_t fNElectronsAbove10GeVMin,fNElectronsAbove10GeVMax;    // Definition of the range of the # of electorns above 10GeV
  Bool_t fNElectronsAbove10GeVFlag;                           // Shows whether this cut is used or not  
  Int_t fNElectronsMin, fNElectronsMax;                       // # of electrons range
  Bool_t fNElectronsFlag;                                     // Shows whether this cut is used or not
  Int_t fNMuonsMin, fNMuonsMax;                               // # of muons range
  Bool_t fNMuonsFlag;                                         // Shows whether this cut is used or not
  Int_t fNPionsMin, fNPionsMax;                               // # of pions range
  Bool_t fNPionsFlag;                                         // Shows whether this cut is used or not
  Int_t fNKaonsMin, fNKaonsMax;                               // # of kaons range
  Bool_t fNKaonsFlag;                                         // Shows whether this cut is used or not
  Int_t fNProtonsMin, fNProtonsMax;                           // # of protons range
  Bool_t fNProtonsFlag;                                       // Shows whether this cut is used or not
  Int_t fNLambdasMin, fNLambdasMax;                           // # of lambdas range
  Bool_t fNLambdasFlag;                                       // Shows whether this cut is used or not
  Int_t fNPhotonsMin, fNPhotonsMax;                           // # of photons range
  Bool_t fNPhotonFlag;                                        // Shows whether this cut is used or not
  Int_t fNPi0sMin, fNPi0sMax;                                 // # of Pi0s range
  Bool_t fNPi0sFlag;                                          // Shows whether this cut is used or not
  Int_t fNNeutronsMin, fNNeutronsMax;                         // # of neutrons range
  Bool_t fNNeutronsFlag;                                      // Shows whether this cut is used or not
  Int_t fNKaon0sMin, fNKaon0sMax;                             // # of K0s range
  Bool_t fNKaon0sFlag;                                        // Shows whether this cut is used or not  
  Float_t fTotalPMin, fTotalPMax;                             // Range of the sum of the momentum per event
  Bool_t fTotalPFlag;                                         // Shows whether this cut is used or not
  Float_t fMeanPtMin, fMeanPtMax;                             // Range of mean Pt per event
  Bool_t fMeanPtFlag;                                         // Shows whether this cut is used or not
  Float_t fTopPtMin;                                          // Max Pt for each event
  Bool_t fTopPtMinFlag;                                       // Shows whether this cut is used or not
  Float_t fTotalNeutralPMin, fTotalNeutralPMax;               // Sum of the momentum per event for neutral
  Bool_t fTotalNeutralPFlag;                                  // Shows whether this cut is used or not
  Float_t fMeanNeutralPtMin, fMeanNeutralPtMax;               // Mean Pt per event for neutral
  Bool_t fMeanNeutralPtFlag;                                  // Shows whether this cut is used or not
  Float_t fTopNeutralPtMin;                                   // Minimum value for highest Pt for the event for neutral
  Bool_t fTopNeutralPtMinFlag;                                // Shows whether this cut is used or not
  Float_t fEventPlaneAngleMin, fEventPlaneAngleMax;           // event plane info
  Bool_t fEventPlaneAngleFlag;                                // Shows whether this cut is used or not
  Float_t fHBTRadiiMin, fHBTRadiiMax;                         // HBT info
  Bool_t fHBTRadiiFlag;                                       // Shows whether this cut is used or not

  ClassDef(AliEventTagCuts, 2)
};

#endif
