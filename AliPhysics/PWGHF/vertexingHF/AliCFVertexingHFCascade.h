#ifndef ALICFVERTEXINGHFCASCADE_H
#define ALICFVERTEXINGHFCASCADE_H

/**************************************************************************
 * Copyright(c) 1998-2011, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------------
/// \class AliCFVertexingHFCascade
/// \brief Class for HF corrections as a function of many variables and steps
/// For D* and other cascades
///
/// \author Author : A.GRELLI - a.grelli@uu.nl UTRECHT
//-----------------------------------------------------------------------


#include "AliCFVertexingHF.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"

class AliAODMCParticle;
class TClonesArray;
class AliCFVertexingHF;
class AliESDtrack;
class TDatabasePDG;
class AliPIDResponse;

class AliCFVertexingHFCascade : public AliCFVertexingHF{
 public:
		
  AliCFVertexingHFCascade();
  AliCFVertexingHFCascade(TClonesArray *mcArray, UShort_t originDselection);
	
  //virtual ~AliCFVertexingHFCascade(){};
  
  
  Bool_t GetGeneratedValuesFromMCParticle(Double_t* /*vectorMC*/);
  Bool_t GetRecoValuesFromCandidate(Double_t* /*vectorReco*/ ) const;
  Bool_t CheckMCChannelDecay()const;
  
  Bool_t SetRecoCandidateParam(AliAODRecoDecayHF *recoCand);
  //Bool_t EvaluateIfD0toKpi(AliAODMCParticle* neutralDaugh, Double_t* VectorD0)const;
  Bool_t EvaluateIfCorrectNeutrDaugh(AliAODMCParticle* neutralDaugh, Double_t* VectorD0)const;

  void SetPtAccCut(Float_t* ptAccCut);
  void SetEtaAccCut(Float_t* etaAccCut);
  void SetAccCut(Float_t* ptAccCut, Float_t* etaAccCut);
  void SetAccCut();

  Double_t GetEtaProng(Int_t iProng)const;
  Double_t GetPtProng(Int_t iProng) const;

  void SetPDGcascade(Int_t pdg)    {fPDGcascade = pdg;}
  void SetPDGbachelor(Int_t pdg)   {fPDGbachelor = pdg;}
  void SetPDGneutrDaugh(Int_t pdg)         {fPDGneutrDaugh = pdg;}
  void SetPDGneutrDaughForMC(Int_t pdg)         {fPDGneutrDaughForMC = pdg;}
  void SetPDGneutrDaughPositive(Int_t pdg) {fPDGneutrDaughPositive = pdg;}
  void SetPDGneutrDaughNegative(Int_t pdg) {fPDGneutrDaughNegative = pdg;}
  void SetPrimaryVertex(AliAODVertex* vtx) {fPrimVtx = vtx;}

  Int_t GetPDGcascade()    const {return fPDGcascade;}
  Int_t GetPDGbachelor()   const {return fPDGbachelor;}
  Int_t GetPDGneutrDaugh()         const {return fPDGneutrDaugh;}
  Int_t GetPDGneutrDaughForMC()         const {return fPDGneutrDaughForMC;}
  Int_t GetPDGneutrDaughPositive() const {return fPDGneutrDaughPositive;}
  Int_t GetPDGneutrDaughNegative() const {return fPDGneutrDaughNegative;}
  AliAODVertex* GetPrimaryVertex() const {return fPrimVtx;}

  Bool_t CheckAdditionalCuts(AliPIDResponse* pidResponse) const;

  void SetUseCutsForTMVA(Bool_t useCutsForTMVA) {fUseCutsForTMVA = useCutsForTMVA;}
  Bool_t GetUseCutsForTMVA() const {return fUseCutsForTMVA;}

  void SetCutOnMomConservation(Float_t cut) {fCutOnMomConservation = cut;}
  Bool_t GetCutOnMomConservation() const {return fCutOnMomConservation;}

 protected:
  
  
 private:	
  AliCFVertexingHFCascade(const AliCFVertexingHFCascade& c);
  AliCFVertexingHFCascade& operator= (const AliCFVertexingHFCascade& other);

  Int_t fPDGcascade;   /// pdg code of the cascade
  Int_t fPDGbachelor;  /// pdg code of the bachelor
  Int_t fPDGneutrDaugh;        /// pdg code of the V0
  Int_t fPDGneutrDaughForMC;     /// pdg code of the V0
  Int_t fPDGneutrDaughPositive;  /// pdg code of the positive daughter of the V0
  Int_t fPDGneutrDaughNegative;  /// pdg code of the negative daughter of the V0
  AliAODVertex* fPrimVtx;        /// primaryVertex
  Bool_t fUseCutsForTMVA;        /// flag to decide whether to use or not the preselection
                                 /// cuts of the TMVA when filling the CF
  Float_t fCutOnMomConservation; /// cut on momentum conservation

  /// \cond CLASSIMP    
  ClassDef(AliCFVertexingHFCascade,4); /// CF class for D* and other cascades
    /// \endcond
};

#endif
