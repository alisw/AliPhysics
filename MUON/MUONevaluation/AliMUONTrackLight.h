#ifndef ALIMUONTRACKLIGHT_H
#define ALIMUONTRACKLIGHT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 06/09/2006

/// \ingroup evaluation
/// \class AliMUONTrackLight
/// \brief Compact information for the muon generated tracks
/// 
/// Compact information for the muon generated tracks in the MUON arm 
/// useful at the last stage of the analysis chain
/// provides a link between the reconstructed track and the generated particle 
/// stores kinematical information at gen. and rec. level and 
/// the decay history of the muon, allowing the identification of the 
/// mother process 
/// 
/// To be used together with AliMUONPairLight
///
/// \author This class was prepared by INFN Cagliari, July 2006
/// (authors: H.Woehri, A.de Falco)

// ROOT classes
#include "TLorentzVector.h"

class AliMUONTrack;
class AliESDMuonTrack;
class AliStack;
class TParticle;
class AliMUONVTrackStore;

class AliMUONTrackLight : public TObject {
 public: 
  AliMUONTrackLight(); 
  AliMUONTrackLight(AliESDMuonTrack* muonTrack); 
  AliMUONTrackLight(const AliMUONTrackLight &muonCopy);
  AliMUONTrackLight& operator=(const AliMUONTrackLight&);
  virtual ~AliMUONTrackLight(); 
  
  /// Set 4-momentum of the generated particle
  void SetPGen(TLorentzVector pgen) {fPgen = pgen;}
  /// Return 4-momentum of the generated particle  
  TLorentzVector GetPGen() const {return fPgen;}
  /// Set reconstructed 4-momentum
  void SetPRec(TLorentzVector prec) {fPrec = prec;}
  /// Return reconstructed 4-momentum 
  TLorentzVector GetPRec() const {return fPrec;}
  /// Set primary vertex position from the ITS
  void SetVertex(Double_t *xyz) {for (Int_t i=0; i<3; i++) fXYZ[i]=xyz[i];}
  /// Return primary vertex x position from the ITS 
  Double_t GetX() const { return fXYZ[0]; } 
  /// Return primary vertex y position from the ITS 
  Double_t GetY() const { return fXYZ[1]; } 
  /// Return primary vertex z position from the ITS 
  Double_t GetZ() const { return fXYZ[2]; } 
  /// Return  primary vertex position from the ITS
  Double_t* GetVertex() { return fXYZ; } 
  /// Set chi2 / ndf in the MUON track fit
  void SetChi2(Double_t chi2) {fChi2=chi2;}
  /// Return chi2 / ndf in the MUON track fit 
  Double_t GetChi2() const { return fChi2; } 
  /// Set weight assigned to the muon
  void SetWeight(Double_t w) {fWeight=w;}
  /// Return weight assigned to the muon 
  Double_t GetWeight() const { return fWeight; } 
  
  /// Set muon charge 
  void SetCharge(Int_t charge) {fCharge = charge;}
  /// Return muon charge 
  Int_t GetCharge() const {return fCharge;}
  /// Return hadronised parents and grandparents 
  Int_t GetParentPDGCode(Int_t index = 0) const { return fParentPDGCode[index]; } 
  /// Return line of Pythia output for hadronised parents & grandparents 
  Int_t GetParentPythiaLine(Int_t index = 0) const { return fParentPythiaLine[index]; } 
  /// Return pdg of the string [0], quarks/gluons [1,2], sometimes proton [3] 
  Int_t GetQuarkPDGCode(Int_t index = 0) const { return fQuarkPDGCode[index]; } 
  /// Return line of Pythia output for string [0] and quarks [1,2], sometimes proton [3] 
  Int_t GetQuarkPythiaLine(Int_t index = 0) const { return fQuarkPythiaLine[index]; } 
  /// Return line of Pythia output for string [0] and quarks [1,2], sometimes proton [3] 
  Int_t GetTrackPythiaLine() const {return fTrackPythiaLine;}
  /// Return pdg code of the rec. track (in general will be a muon) 
  Int_t GetTrackPDGCode() const {return fTrackPDGCode;}
  /// Set line of kin. stack where rec. track (in general, the muon) is stored
  void SetTrackPythiaLine(Int_t trackLine) {fTrackPythiaLine = trackLine;}
  /// Set pdg code of the rec. track (in general will be a muon)
  void SetTrackPDGCode(Int_t trackPdg) {fTrackPDGCode = trackPdg;}
  void FillFromESD(AliESDMuonTrack* muonTrack, Double_t zvert=-9999);
  void FillFromAliMUONTrack(AliMUONTrack *trackReco,Double_t zvert=-9999);
  /// Return info if is a muon 
  Bool_t IsAMuon() const { return (TMath::Abs(fTrackPDGCode)==13); }
  Bool_t IsParentPionOrKaon(Int_t idParent = 0);
  void SetPxPyPz(Double_t px, Double_t py, Double_t pz); 
  /// Set flag for trigger 
  void SetTriggered(Bool_t isTriggered) { fIsTriggered = isTriggered; } 
  /// Return flag for trigger  
  Bool_t IsTriggered() const { return fIsTriggered; } 
  /// Return acually filled no. of *fragmented* parents 
  Int_t GetNParents() const {return fNParents;}
  void FillMuonHistory(AliStack *stack, TParticle *part);
  Bool_t IsB0(Int_t intTest) const;//checks if the provided PDG code corresponds to a neutral B meson
  Bool_t IsMotherAResonance(Int_t index=0) const;
  /// Return flag for oscillation 
  Bool_t GetOscillation() const {return fOscillation;}
  virtual void PrintInfo(const Option_t* opt); //"H" muon's decay history
  //"K" muon kinematics
  //"A" all variables
  Int_t GetParentFlavour(Int_t idParent=0) const;
  Bool_t IsDiquark(Int_t pdg) const;
protected:
  static const Int_t fgkNParentsMax = 5;   ///< maximum number of parents
  TLorentzVector fPrec; ///< reconstructed 4-momentum
  Double_t fXYZ[3];     ///< primary vertex position from the ITS 
  Bool_t fIsTriggered;  ///< flag for trigger 
  Int_t fCharge;        ///< muon charge 
  Double_t fChi2;       ///< chi2 / ndf in the MUON track fit 
  Float_t fCentr;       ///< centrality 
  TLorentzVector fPgen;   ///< 4-momentum of the generated particle
  Int_t fTrackPythiaLine; ///< line of kin. stack where rec. track (in general, the muon) is stored
  Int_t fTrackPDGCode; ///< pdg code of the rec. track (in general will be a muon)
  Int_t fParentPDGCode[fgkNParentsMax]; ///< hadronised parents and grandparents 
  Int_t fParentPythiaLine[fgkNParentsMax];///< line of Pythia output for hadronised parents & grandparents
  Int_t fQuarkPDGCode[4]; ///< pdg of the string [0], quarks/gluons [1,2], sometimes proton [3] 
  Int_t fQuarkPythiaLine[4]; ///< line of Pythia output for string [0] and quarks [1,2], sometimes proton [3]
  Bool_t fOscillation; ///< flag for oscillation
  Int_t fNParents; ///< acually filled no. of *fragmented* parents
  Double_t fWeight; ///< weight assigned to the muon
  /// Set flag for oscillation
  void SetOscillation(Bool_t oscillation) { fOscillation = oscillation; }
  /// Set hadronised parents and grandparents 
  void SetParentPDGCode(Int_t index, Int_t pdg) { fParentPDGCode[index] = pdg; } 
  /// Set line of Pythia output for hadronised parents & grandparents
  void SetParentPythiaLine(Int_t index, Int_t line) { fParentPythiaLine[index] = line; } 
  /// Set pdg of the string [0], quarks/gluons [1,2], sometimes proton [3] 
  void SetQuarkPDGCode(Int_t index, Int_t pdg){ fQuarkPDGCode[index] = pdg; }
  /// Set line of Pythia output for string [0] and quarks [1,2], sometimes proton [3]
  void SetQuarkPythiaLine(Int_t index, Int_t line){ fQuarkPythiaLine[index] = line; }
  void ResetQuarkInfo();

  ClassDef(AliMUONTrackLight,1)  /// Muon Track for analysis
}; 

#endif
