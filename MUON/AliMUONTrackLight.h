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

#include <TClonesArray.h>
#include <TLorentzVector.h>

class AliMUONTrack;
class AliESDMuonTrack;
class AliRunLoader;

class TParticle;

class AliMUONTrackLight : public TObject { 
 public: 
  AliMUONTrackLight(); 
  AliMUONTrackLight(AliESDMuonTrack* muonTrack); 
  AliMUONTrackLight(const AliMUONTrackLight &muonCopy);
  virtual ~AliMUONTrackLight(){;} 
  
  void SetPGen(TLorentzVector pgen) {fPgen = pgen;}
  TLorentzVector GetPGen() const {return fPgen;}
  void SetPRec(TLorentzVector prec) {fPrec = prec;}
  TLorentzVector GetPRec() const {return fPrec;}
  void SetVertex(Double_t *xyz) {for (Int_t i=0; i<3; i++) fXYZ[i]=xyz[i];}
  Double_t* GetVertex() { return fXYZ; } 
  void SetCharge(Int_t charge) {fCharge = charge;}
  Int_t GetCharge() const {return fCharge;}
  Int_t GetParentPDGCode(Int_t index = 0) const { return fParentPDGCode[index]; } 
  Int_t GetParentPythiaLine(Int_t index = 0) const { return fParentPythiaLine[index]; } 
  Int_t GetQuarkPDGCode(Int_t index = 0) const { return fQuarkPDGCode[index]; } 
  Int_t GetQuarkPythiaLine(Int_t index = 0) const { return fQuarkPythiaLine[index]; } 
  Int_t GetTrackPythiaLine() const {return fTrackPythiaLine;}
  Int_t GetTrackPDGCode() const {return fTrackPDGCode;}
  void SetTrackPythiaLine(Int_t trackLine) {fTrackPythiaLine = trackLine;}
  void SetTrackPDGCode(Int_t trackPdg) {fTrackPDGCode = trackPdg;}
  void ComputePRecAndChargeFromESD(AliESDMuonTrack* muonTrack);
  Bool_t IsAMuon() const { return (TMath::Abs(fTrackPDGCode)==13); }
  void SetPxPyPz(Double_t px, Double_t py, Double_t pz); 
  void SetTriggered(Bool_t isTriggered) { fIsTriggered = isTriggered; } 
  Bool_t IsTriggered() const { return fIsTriggered; } 
  TParticle* FindRefTrack(AliMUONTrack* trackReco, TClonesArray* trackRefArray, AliRunLoader *runLoader); 
  Int_t TrackCheck(Bool_t *compTrack);
  Int_t GetNParents() const {return fNParents;}
  void FillMuonHistory(AliRunLoader *runLoader, TParticle *part);
  Bool_t IsB0(Int_t intTest) const;//checks if the provided PDG code corresponds to a neutral B meson
  Bool_t IsMotherAResonance(Int_t index=0) const;
  Bool_t GetOscillation() const {return fOscillation;}
  virtual void PrintInfo(Option_t* opt); //"H" muon's decay history
  //"K" muon kinematics
  //"A" all variables
  Int_t GetParentFlavour(Int_t idParent=0) const;
protected:
  static const Int_t fgkNParentsMax = 5;   /// maximum number of parents
  TLorentzVector fPrec; /// reconstructed 4-momentum
  Double_t fXYZ[3];     /// primary vertex position from the ITS 
  Bool_t fIsTriggered;  /// flag for trigger 
  Int_t fCharge;        /// muon charge 
  Float_t fCentr;       /// centrality 
  TLorentzVector fPgen;   /// 4-momentum of the generated particle
  Int_t fTrackPythiaLine; ///line of kin. stack where rec. track (in general, the muon) is stored
  Int_t fTrackPDGCode; ///pdg code of the rec. track (in general will be a muon)
  Int_t fParentPDGCode[fgkNParentsMax]; /// hadronised parents and grandparents 
  Int_t fParentPythiaLine[fgkNParentsMax];///line of Pythia output for hadronised parents & grandparents
  Int_t fQuarkPDGCode[4]; /// pdg of the string [0], quarks/gluons [1,2], sometimes proton [3] 
  Int_t fQuarkPythiaLine[4]; ///line of Pythia output for string [0] and quarks [1,2], sometimes proton [3]
  Bool_t fOscillation; /// flag for oscillation
  Int_t fNParents; ///acually filled no. of *fragmented* parents
  void SetOscillation(Bool_t oscillation) { fOscillation = oscillation; }
  void SetParentPDGCode(Int_t index, Int_t pdg) { fParentPDGCode[index] = pdg; } 
  void SetParentPythiaLine(Int_t index, Int_t line) { fParentPythiaLine[index] = line; } 
  void SetQuarkPDGCode(Int_t index, Int_t pdg){ fQuarkPDGCode[index] = pdg; }
  void SetQuarkPythiaLine(Int_t index, Int_t line){ fQuarkPythiaLine[index] = line; }
  void ResetQuarkInfo();

  ClassDef(AliMUONTrackLight,1)  /// Muon Track for analysis
}; 

#endif
