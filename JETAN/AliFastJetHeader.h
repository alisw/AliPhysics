#ifndef ALIFASTJETHEADER_H
#define ALIFASTJETHEADER_H
 
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

 
class AliFastJetHeader : public AliJetHeader
{
 public:
 
  AliFastJetHeader();
  virtual ~AliFastJetHeader() { }

  // Getters
  Double_t                     GetRparam()            const {return fRparam;}
  fastjet::JetAlgorithm        GetAlgorithm()         const {return fAlgorithm;}
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
  
  // Setters
  void SetRparam(Double_t f)                           {fRparam = f;}
  void SetAlgorithm(fastjet::JetAlgorithm f)           {fAlgorithm = f;}
  void SetStrategy(fastjet::Strategy f)                {fStrategy = f;}
  void SetRecombScheme(fastjet::RecombinationScheme f) {fRecombScheme = f;}
  void SetGhostEtaMax(Double_t f)                      {fGhostEtaMax = f;}
  void SetGhostArea(Double_t f)                        {fGhostArea = f;}
  void SetActiveAreaRepeats(Int_t f)                   {fActiveAreaRepeats =f;}
  void SetAreaType(fastjet::AreaType f)                {fAreaType = f;}
  void SetRapRange(Double_t fmin, Double_t fmax)       {fRapMin = fmin; fRapMax = fmax;}
  void SetPhiRange(Double_t fmin, Double_t fmax)       {fPhiMin = fmin; fPhiMax = fmax;}
  
  void SetComment(TString com)       {fComment=com;}
  void SetComment(const char* com)   {fComment=com;}
  
  // others
  void PrintParameters() const;

 protected:

  //fastjet::JetDefinition parameters
  Double_t fRparam;
  fastjet::JetAlgorithm fAlgorithm; //fastjet::kt_algorithm
  fastjet::Strategy fStrategy;  //= fastjet::Best;
  fastjet::RecombinationScheme fRecombScheme; // = fastjet::BIpt_scheme;
  
  //fastjet::GhostedAreaSpec parameters
  Double_t fGhostEtaMax;       // max area of ghosts
  Double_t fGhostArea;         // ghost area
  Int_t    fActiveAreaRepeats; // number of repetitions of active area 
  
  //fastjet::AreaDefinition parameters
  fastjet::AreaType fAreaType; // the are type
  
  //fastjet::ClusterSequenceArea options parameters
  Double_t fPtMin; //jets with pt > ptmin
  
  //fastjet::RangeDefinition parameters 
  Double_t fRapMax, fRapMin; // rapidity range of background sub 
  Double_t fPhiMax, fPhiMin; // phi range of background sub
  
  
  ClassDef(AliFastJetHeader,2)
};
 
#endif
