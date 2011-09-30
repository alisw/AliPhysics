#ifndef ALIFASTJETHEADERV1_H
#define ALIFASTJETHEADERV1_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//---------------------------------------------------------------------
// FastJet v2.3.4 finder algorithm interface
// Finder Header Class 
// Author: Rafael.Diaz.Valdes@cern.ch
//---------------------------------------------------------------------
 

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/JetDefinition.hh"

#include "AliJetHeader.h"

 
class AliFastJetHeaderV1 : public AliJetHeader
{
 public:
 
  AliFastJetHeaderV1();
  virtual ~AliFastJetHeaderV1() { }

  // Getters
  Double_t                     GetRparam()            const {return fRparam;}
  fastjet::JetAlgorithm        GetAlgorithm()         const {return fAlgorithm;}
  fastjet::JetAlgorithm        GetBGAlgorithm()       const {return fBGAlgorithm;} 
  fastjet::Strategy            GetStrategy()          const {return fStrategy;}
  fastjet::RecombinationScheme GetRecombScheme()      const {return fRecombScheme;}
  Double_t                     GetGhostEtaMax()       const {return fGhostEtaMax;}
  Double_t                     GetGhostArea()         const {return fGhostArea;}
  Int_t                        GetActiveAreaRepeats() const {return fActiveAreaRepeats;}
  fastjet::AreaType            GetAreaType()          const {return fAreaType;}
  Double_t                     GetPtMin()             const {return fPtMin;}
  Double_t                     GetRapMax()            const {return fRapMax;}
  Double_t                     GetRapMin()            const {return fRapMin;}
  Double_t                     GetPhiMax()            const {return fPhiMax;}
  Double_t                     GetPhiMin()            const {return fPhiMin;}
  // Added temporarily !!! To be removed if not necessary
  Float_t                      GetMinCellEt()         const {return fMinCellEt;} 
  Bool_t                       GetBGMode()            const {return fBGMode;}
  Double_t                     GetRparamBkg()            const {return fRparamBkg;}
  Bool_t                       Use4VectorArea()       const {return fUse4VectorArea;}  

  // Setters
  void SetRparam(Double_t f)                           {fRparam = f;}
  void SetAlgorithm(fastjet::JetAlgorithm f)           {fAlgorithm = f;}
  void SetBGAlgorithm(fastjet::JetAlgorithm f)         {fBGAlgorithm = f;}
  void SetStrategy(fastjet::Strategy f)                {fStrategy = f;}
  void SetRecombScheme(fastjet::RecombinationScheme f) {fRecombScheme = f;}
  void SetGhostEtaMax(Double_t f)                      {fGhostEtaMax = f;}
  void SetGhostArea(Double_t f)                        {fGhostArea = f;}
  void SetActiveAreaRepeats(Int_t f)                   {fActiveAreaRepeats =f;}
  void SetAreaType(fastjet::AreaType f)                {fAreaType = f;}
  void SetRapRange(Double_t fmin, Double_t fmax)       {fRapMin = fmin; fRapMax = fmax;}
  void SetPhiRange(Double_t fmin, Double_t fmax)       {fPhiMin = fmin; fPhiMax = fmax;}
  void SetPtMin(Double_t ptmin)                        {fPtMin = ptmin;}
  void SetBGMode(Bool_t bgmode)                        {fBGMode = bgmode;}
  void SetUse4VectorArea()                             {fUse4VectorArea = kTRUE;}
 
  void SetComment(TString com) {fComment=com;}
  void SetComment(const char* com) {AliJetHeader::SetComment(com);}
  
  void SetRparamBkg(Double_t f)                           {fRparamBkg = f;}

  // others
  void PrintParameters() const;

 protected:

  //fastjet::JetDefinition parameters
  Double_t fRparam;   // R param
  Double_t fRparamBkg;//R param for bkg calculation
  fastjet::JetAlgorithm fAlgorithm; //fastjet::kt_algorithm
  fastjet::JetAlgorithm fBGAlgorithm; //fastjet::kt_algorithm
  fastjet::Strategy fStrategy;  //= fastjet::Best;
  fastjet::RecombinationScheme fRecombScheme; // = fastjet::BIpt_scheme;
  
  //fastjet::GhostedAreaSpec parameters
  Double_t fGhostEtaMax;       // Max eta for ghosts
  Double_t fGhostArea;         // Ghost area 
  Int_t    fActiveAreaRepeats; // Active are repetitions
  
  //fastjet::AreaDefinition parameters
  fastjet::AreaType fAreaType; // area types
  
  //fastjet::ClusterSequenceArea options parameters
  Double_t fPtMin; //jets with pt > ptmin
  Float_t  fMinCellEt;       //  Min Et in one cell

  //fastjet::RangeDefinition parameters 
  Double_t fRapMax, fRapMin; // rapidity range of background sub 
  Double_t fPhiMax, fPhiMin; // phi range of background sub
  Bool_t   fBGMode;          // Do we subtract BG or not?
  Bool_t   fUse4VectorArea;  // Toggle use of 4-vector area 

  ClassDef(AliFastJetHeaderV1,3)
};
 
#endif
