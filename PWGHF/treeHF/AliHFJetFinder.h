#ifndef ALIHFJETFINDER_H
#define ALIHFJETFINDER_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFJetFinder
// \helper class to handle jet finding, matching and substructure
// \authors:
// N. Zardoshti, nima.zardoshti@cern.ch
/////////////////////////////////////////////////////////////

#include "TObject.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODMCParticle.h"
#include "AliHFJet.h"
#include "AliTLorentzVector.h"
#include "AliFJWrapper.h"
#include "FJ_includes.h"

class AliFJWrapper;

class AliHFJetFinder : public TObject
{
  public:

  enum JetAlgorithm{
    antikt,
    kt,
    ca
  };
  enum RecombScheme{
    e_scheme,
    pt_scheme
  };
  enum AreaType{
    active,
    passive,
    voronoi
  };
  enum Charge{
    charged,
    full,
    neutral
  };

  AliHFJetFinder();
  AliHFJetFinder(char *name);

  virtual ~AliHFJetFinder();

  void SetFJWrapper(); 
  AliHFJet GetHFJet(TClonesArray *array, AliAODRecoDecayHF *cand, Double_t invmass=0);
  AliHFJet GetHFMCJet(TClonesArray *array, AliAODMCParticle *mcpart);
  std::vector<AliHFJet> GetHFJets(TClonesArray *array, AliAODRecoDecayHF *cand, Double_t invmass);
  std::vector<AliHFJet> GetHFMCJets(TClonesArray *array, AliAODMCParticle *mcpart);
  std::vector<AliHFJet> GetJets(TClonesArray *array);
  std::vector<AliHFJet> GetMCJets(TClonesArray *array);
  void FindJets(TClonesArray *array, AliAODRecoDecayHF *cand=nullptr, Double_t invmass=0);
  void FindMCJets(TClonesArray *array, AliAODMCParticle *mcpart=nullptr);
  #if !defined(__CINT__) && !defined(__MAKECINT__)
  void SetJetVariables(AliHFJet& hfjet, const std::vector<fastjet::PseudoJet>& constituents, const fastjet::PseudoJet& jet, Int_t jetID, AliAODRecoDecayHF *cand=nullptr);
  void SetMCJetVariables(AliHFJet& hfjet, const std::vector<fastjet::PseudoJet>& constituents, const fastjet::PseudoJet& jet, Int_t jetID, AliAODMCParticle *mcpart=nullptr);
  void SetJetSubstructureVariables(AliHFJet& hfjet, const fastjet::PseudoJet& jet, const std::vector<fastjet::PseudoJet>& constituents);
  fastjet::JetFinder JetAlgorithm(Int_t jetalgo);
  fastjet::RecombinationScheme RecombinationScheme(Int_t recombscheme);
  fastjet::AreaType AreaType(Int_t area);
  #endif
  Bool_t CheckTrack(AliAODTrack *track);
  Bool_t CheckFilterBits(AliAODTrack *track);
  Bool_t CheckParticle(AliAODMCParticle *particle);
  Int_t Find_Candidate_Jet();
  Float_t RelativePhi(Float_t phi_1, Float_t phi_2);


  void SetDoJetSubstructure(Bool_t b)      {fDoJetSubstructure=b;}
  void SetMinJetPt(Float_t f)              {fMinJetPt = f;}
  void SetJetRadius(Float_t f)             {fJetRadius = f;}
  void SetJetAlgorithm(Int_t i)            {fJetAlgorithm=i;}
  void SetJetRecombScheme(Int_t i)         {fJetRecombScheme = i;}
  void SetGhostArea(Float_t f)             {fJetGhostArea = f;}
  void SetJetAreaType(Int_t i)             {fJetAreaType = i;}
  void SetMinSubJetPt(Float_t f)           {fMinSubJetPt = f;}
  void SetSubJetRadius(Float_t f)          {fSubJetRadius = f;}
  void SetSubJetAlgorithm(Int_t i)         {fSubJetAlgorithm=i;}
  void SetSubJetRecombScheme(Int_t i)      {fSubJetRecombScheme = i;}
  void SetSoftDropParams(Float_t f1, Float_t f2) {fSoftDropZCut = f1; fSoftDropBeta = f2;}
  void SetMinTrackPt(Float_t f)            {fMinTrackPt = f;}
  void SetMaxTrackPt(Float_t f)            {fMaxTrackPt = f;}
  void SetMaxTrackEta(Float_t f)           {fMaxTrackEta = f;}
  void SetMinParticlePt(Float_t f)         {fMinParticlePt = f;}
  void SetMaxParticlePt(Float_t f)         {fMaxParticlePt = f;}
  void SetMaxParticleEta(Float_t f)        {fMaxParticleEta = f;}
  void SetCharged(Int_t i)                 {fCharged = i;}
  void SetTrackingEfficiency(Double_t d)   {fTrackingEfficiency = d;}
  

  Float_t                  fMinJetPt;
  Float_t                  fJetRadius;
  Int_t                    fJetAlgorithm;
  Int_t                    fJetRecombScheme;
  Float_t                  fJetGhostArea;
  Int_t                    fJetAreaType;

  Float_t                  fMinSubJetPt;
  Float_t                  fSubJetRadius;
  Int_t                    fSubJetAlgorithm;
  Int_t                    fSubJetRecombScheme;
  Float_t                  fSoftDropZCut;
  Float_t                  fSoftDropBeta;

  Float_t                  fMinTrackPt;
  Float_t                  fMaxTrackPt;
  Float_t                  fMinTrackE;
  Float_t                  fMaxTrackE;
  Float_t                  fMaxTrackEta;
  Float_t                  fMaxTrackPhi;

  Float_t                  fMinParticlePt;
  Float_t                  fMaxParticlePt;
  Float_t                  fMaxParticleEta;
  Float_t                  fCharged;

  Double_t                 fTrackingEfficiency;
  
  Bool_t                   fDoJetSubstructure;
  AliFJWrapper            *fFastJetWrapper;
  std::vector<std::pair<Int_t,Double_t>>   fConstituentCharge;




  /// \cond CLASSIMP
  ClassDef(AliHFJetFinder,1); ///
  /// \endcond
};
#endif
