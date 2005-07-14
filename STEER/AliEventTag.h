#ifndef ALIEVENTTAG_H
#define ALIEVENTTAG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliEventTag
//   This is the class to deal with the tags for the event level
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------

#include "TObject.h"
#include "TClonesArray.h"

//______________________________________________________________________________
class AliEventTag : public TObject
{
 private:
  Int_t    fAliceEventId;                      //The event id
  Int_t    fGUID;			       //The unique identifier of the file
  Int_t    fNumberOfParticipants;    	       //Number of participants
  Float_t  fImpactParameter;		       //The impact parameter
 
  Int_t    fPrimaryVertexFlag;		       //Primary vertex flag: 0->not found, 1->found

  Float_t  fPrimaryVertexX;		       //Primary vertex - X coordinate
  Float_t  fPrimaryVertexY;		       //Primary vertex - Y coordinate
  Float_t  fPrimaryVertexZ;		       //Primary vertex - Z coordinate

  Float_t  fPrimaryVertexZError;	       //Primary vertex - Z coordinate - error

  Int_t    fTriggerInfo;		       //Information from trigger
  Float_t  fZDCNeutronEnergy;		       //ZDC info - neutron
  Float_t  fZDCProtonEnergy;		       //ZDC info - proton
  Float_t  fZDCEMEnergy;		       //ZDC info - em
  Float_t  fT0VertexZ;			       //T0 info
  Int_t    fNumberOfTracks;		       //Multiplicity
  Int_t    fNumberOfPositiveTracks;	       //Multiplicity of positive tracks
  Int_t    fNumberOfNegativeTracks;	       //Multiplicity of negative tracks
  Int_t    fNumberOfNeutralTracks;	       //Multiplicity of neutral tracks
  Int_t    fNumberOfV0s;		       //Number of V0s
  Int_t    fNumberOfCascades;		       //Number of cascades
  Int_t    fNumberOfKinks;		       //Number of kinks
  Int_t    fNumberOfPMDTracks;		       //PMD tracks
  Int_t    fNumberOfPHOSTracks;		       //PHOS tracks
  Int_t    fNumberOfEMCALTracks;	       //EMCAL tracks
  Int_t    fNumberOfFMDTracks;		       //FMD tracks
  Int_t    fNumberOfJetCandidates;	       //Jet candidates

  Float_t  fMaxJetEnergy;                      //jet energy info

  Int_t    fNumberOfHardPhotonsCandidates;     //Hard photons candidates

  Float_t  fMaxNeutralEnergy;                   //neutral energy info
  Int_t    fNumberOfChargedAbovePtRange;
  Int_t    fNumberOfMuonsAbovePtRange;
  Int_t    fNumberOfElectronsAbovePtRange;



  Int_t    fNumberOfElectrons;		       //Number of electrons
  Int_t    fNumberOfMuons;		       //Number of muons
  Int_t    fNumberOfPions;		       //Number of pions
  Int_t    fNumberOfKaons;		       //Number of kaons
  Int_t    fNumberOfProtons;		       //Number of protons
  Int_t    fNumberOfLambdas;		       //Number of lambdas

  Int_t    fNumberOfPhotons;
  Int_t    fNumberOfPi0s;
  Int_t    fNumberOfNeutrons;
  Int_t    fNumberOfKaon0s;


  Int_t    fNumberOfJPsiCandidates;	       //JPsi candidates
  Int_t    fNumberOfPsiPrimeCandidates;	       //Psi prime candidates
  Int_t    fNumberOfUpsilonCandidates;	       //Upsilon candidates
  Int_t    fNumberOfUpsilonPrimeCandidates;    //Upsilon prime candidates
  Int_t    fNumberOfUpsilonDoublePrimeCandidates;
  Int_t    fNumberOfCharmParticleCandidates;
  Int_t    fNumberOfBeautyParticleCandidates;
 
  Float_t  fTotalP;			       //Sum of the momentum per event
  Float_t  fMeanPt;			       //Mean Pt per event
  Float_t  fMaxPt;			       //Max Pt for each event

  Float_t  fTotalNeutralP;		       //Sum of the momentum per event for neutral
  Float_t  fMeanNeutralPt;		       //Mean Pt per event for neutral
  Float_t  fMaxNeutralPt;		       //Max Pt for each event for neutral

  Float_t  fEventPlaneAngle;		       //event plane info
  Float_t  fHBTRadii;                          //HBT info

 public:
  AliEventTag();
  AliEventTag(AliEventTag *t);
  virtual ~AliEventTag();
  
  void   SetEventId(Int_t Pid) {fAliceEventId = Pid;}
  void   SetGUID(Int_t Pid) {fGUID = Pid;}

  void   SetNumOfParticipants(Int_t P) {fNumberOfParticipants = P;}
  void   SetImpactParameter(Float_t Pimpact) {fImpactParameter = Pimpact;}

  void   SetVertexX(Float_t Pvx) {fPrimaryVertexX = Pvx;}
  void   SetVertexY(Float_t Pvy) {fPrimaryVertexY = Pvy;}
  void   SetVertexZ(Float_t Pvz) {fPrimaryVertexZ = Pvz;}

  void SetVertexFlag(Int_t i) {fPrimaryVertexFlag = i;}
  void SetVertexZError(Float_t f) { fPrimaryVertexZError = f;}

  void   SetTrigger(Int_t Ptr) {fTriggerInfo = Ptr;}

  void   SetZDCNeutronEnergy(Float_t Pen) {fZDCNeutronEnergy = Pen;}
  void   SetZDCProtonEnergy(Float_t Pen) {fZDCProtonEnergy = Pen;}
  void   SetZDCEMEnergy(Float_t Pen) {fZDCEMEnergy = Pen;}

  void   SetT0VertexZ(Float_t Pvz) {fT0VertexZ = Pvz;}

  void   SetNumOfTracks(Int_t Ptr) {fNumberOfTracks = Ptr;}
  void   SetNumOfPosTracks(Int_t Ptr) {fNumberOfPositiveTracks = Ptr;}
  void   SetNumOfNegTracks(Int_t Ptr) {fNumberOfNegativeTracks = Ptr;}
  void   SetNumOfNeutrTracks(Int_t Ptr) {fNumberOfNeutralTracks = Ptr;}

  void   SetNumOfV0s(Int_t Ptr) {fNumberOfV0s = Ptr;}
  void   SetNumOfCascades(Int_t Ptr) {fNumberOfCascades = Ptr;}
  void   SetNumOfKinks(Int_t Ptr) {fNumberOfKinks = Ptr;}

  void   SetNumOfPMDTracks(Int_t Ptr) {fNumberOfPMDTracks = Ptr;}
  void   SetNumOfPHOSTracks(Int_t Ptr) {fNumberOfPHOSTracks = Ptr;}
  void   SetNumOfEMCALTracks(Int_t Ptr) {fNumberOfEMCALTracks = Ptr;}
  void   SetNumOfFMDTracks(Int_t Ptr) {fNumberOfFMDTracks = Ptr;}

  void   SetNumOfJetCandidates(Int_t Ptr) {fNumberOfJetCandidates = Ptr;}
  void   SetNumOfHardPhotonsCandidates(Int_t Ptr) {fNumberOfHardPhotonsCandidates = Ptr;}


  void   SetMaxJetEnergy(Float_t f) {fMaxJetEnergy = f;}
  void   SetMaxNeutralEnergy(Float_t f) {fMaxNeutralEnergy = f;}
  void   SetNumOfChargedAbovePtRange(Int_t i) {fNumberOfChargedAbovePtRange = i;}
  void   SetNumOfMuonsAbovePtRange(Int_t i) {fNumberOfMuonsAbovePtRange = i;}
  void   SetNumOfElectronsAbovePtRange(Int_t i) {fNumberOfElectronsAbovePtRange = i;}

  void   SetNumOfJPsiCandidates(Int_t Ptr) {fNumberOfJPsiCandidates = Ptr;}
  void   SetNumOfPsiPrimeCandidates(Int_t Ptr) {fNumberOfPsiPrimeCandidates = Ptr;}
  void   SetNumOfUpsilonCandidates(Int_t Ptr) {fNumberOfUpsilonCandidates = Ptr;}
  void   SetNumOfUpsilonPrimeCandidates(Int_t Ptr) {fNumberOfUpsilonPrimeCandidates = Ptr;}
  void   SetNumOfUpsilonDoublePrimeCandidates(Int_t Ptr) {fNumberOfUpsilonDoublePrimeCandidates = Ptr;}
  void   SetNumOfCharmCandidates(Int_t Ptr) {fNumberOfCharmParticleCandidates = Ptr;}
  void   SetNumOfBeautyCandidates(Int_t Ptr) {fNumberOfBeautyParticleCandidates = Ptr;}

  void   SetNumOfElectrons(Int_t Ptr) {fNumberOfElectrons = Ptr;}
  void   SetNumOfMuons(Int_t Ptr) {fNumberOfMuons = Ptr;}
  void   SetNumOfPions(Int_t Ptr) {fNumberOfPions = Ptr;}
  void   SetNumOfKaons(Int_t Ptr) {fNumberOfKaons = Ptr;}
  void   SetNumOfProtons(Int_t Ptr) {fNumberOfProtons = Ptr;}
  void   SetNumOfLambdas(Int_t Ptr) {fNumberOfLambdas = Ptr;}


  void   SetNumOfPhotons(Int_t Ptr) {fNumberOfPhotons = Ptr;}
  void   SetNumOfPi0s(Int_t Ptr) {fNumberOfPi0s = Ptr;}
  void   SetNumOfNeutrons(Int_t Ptr) {fNumberOfNeutrons = Ptr;}
  void   SetNumOfKaon0s(Int_t Ptr) {fNumberOfKaon0s = Ptr;}

  void   SetTotalMomentum(Float_t P) {fTotalP = P;}
  void   SetMeanPt(Float_t Pt) {fMeanPt = Pt;}
  void   SetMaxPt(Float_t Pt) {fMaxPt = Pt;}

  void SetNeutralTotalMomentum(Float_t f) {fTotalNeutralP = f;}
  void SetNeutralMeanPt(Float_t f) {fMeanNeutralPt = f;}
  void SetNeutralMaxPt(Float_t f) {fMaxNeutralPt = f;}

  void SetEventPlaneAngle(Float_t f) {fEventPlaneAngle = f;}
  void SetHBTRadii(Float_t f) {fHBTRadii = f;}



  Int_t   GetEventId() {return fAliceEventId;}
  Int_t   GetGUID() {return fGUID;}

  Int_t   GetNumOfParticipants() {return fNumberOfParticipants;}
  Float_t GetImpactParameter() {return fImpactParameter;}

  Float_t GetVertexX() {return fPrimaryVertexX;}
  Float_t GetVertexY() {return fPrimaryVertexY;}
  Float_t GetVertexZ() {return fPrimaryVertexZ;}

  Int_t GetVertexFlag() {return fPrimaryVertexFlag;}
  Float_t GetVertexZError() {return fPrimaryVertexZError;}



  Int_t   GetTrigger() {return fTriggerInfo;}

  Float_t GetZDCNeutronEnergy() {return fZDCNeutronEnergy;}
  Float_t GetZDCProtonEnergy() {return fZDCProtonEnergy;}
  Float_t GetZDCEMEnergy() {return fZDCEMEnergy;}

  Float_t GetT0VertexZ() {return fT0VertexZ;}

  Int_t   GetNumOfTracks() {return fNumberOfTracks;}
  Int_t   GetNumOfPosTracks() {return fNumberOfPositiveTracks;}
  Int_t   GetNumOfNegTracks() {return fNumberOfNegativeTracks;}
  Int_t   GetNumOfNeutrTracks() {return fNumberOfNeutralTracks;}

  Int_t   GetNumOfV0s() {return fNumberOfV0s;}
  Int_t   GetNumOfCascades() {return fNumberOfCascades;}
  Int_t   GetNumOfKinks() {return fNumberOfKinks;}

  Int_t   GetNumOfPMDTracks() {return fNumberOfPMDTracks;}
  Int_t   GetNumOfPHOSTracks() {return fNumberOfPHOSTracks;}
  Int_t   GetNumOfEMCALTracks() {return fNumberOfEMCALTracks;}
  Int_t   GetNumOfFMDTracks() {return fNumberOfFMDTracks;}

  Int_t   GetNumOfJetCandidates() {return fNumberOfJetCandidates;}
  Int_t   GetNumOfHardPhotonsCandidates() {return fNumberOfHardPhotonsCandidates;}

  Float_t GetMaxJetEnergy() {return fMaxJetEnergy;}
  Float_t GetMaxNeutralEnergy() {return fMaxNeutralEnergy;}
  Int_t   GetNumOfChargedAbovePtRange() {return fNumberOfChargedAbovePtRange;}
  Int_t   GetNumOfMuonsAbovePtRange() {return fNumberOfMuonsAbovePtRange;}
  Int_t   GetNumOfElectronsAbovePtRange() {return fNumberOfElectronsAbovePtRange;}


  Int_t   GetNumOfJPsiCandidates() {return fNumberOfJPsiCandidates;}
  Int_t   GetNumOfPsiPrimeCandidates() {return fNumberOfPsiPrimeCandidates;}
  Int_t   GetNumOfUpsilonCandidates() {return fNumberOfUpsilonCandidates;}
  Int_t   GetNumOfUpsilonPrimeCandidates() {return fNumberOfUpsilonPrimeCandidates;}
  Int_t   GetNumOfUpsilonDoublePrimeCandidates() {return fNumberOfUpsilonDoublePrimeCandidates;}
  Int_t   GetNumOfCharmCandidates() {return fNumberOfCharmParticleCandidates;}
  Int_t   GetNumOfBeautyCandidates() {return fNumberOfBeautyParticleCandidates;}

  Int_t   GetNumOfElectrons() {return fNumberOfElectrons;}
  Int_t   GetNumOfMuons() {return fNumberOfMuons;}
  Int_t   GetNumOfPions() {return fNumberOfPions;}
  Int_t   GetNumOfKaons() {return fNumberOfKaons;}
  Int_t   GetNumOfProtons() {return fNumberOfProtons;}
  Int_t   GetNumOfLambdas() {return fNumberOfLambdas;}


  Int_t   GetNumOfPhotons() {return fNumberOfPhotons;}
  Int_t   GetNumOfPi0s() {return fNumberOfPi0s;}
  Int_t   GetNumOfNeutrons() {return fNumberOfNeutrons;}
  Int_t   GetNumOfKaon0s() {return fNumberOfKaon0s;}


  Float_t GetTotalMomentum() {return fTotalP;}
  Float_t GetMeanPt() {return fMeanPt;}
  Float_t GetMaxPt() {return fMaxPt;}

  Float_t GetNeutralTotalMomentum() {return fTotalNeutralP;}
  Float_t GetNeutralMeanPt() {return fMeanNeutralPt;}
  Float_t GetNeutralMaxPt() {return fMaxNeutralPt;}

  Float_t GetEventPlaneAngle() {return fEventPlaneAngle;}
  Float_t GetHBTRadii() {return fHBTRadii;}

  ClassDef(AliEventTag,1)  //(ClassName, ClassVersion)
    };
//______________________________________________________________________________


#endif
