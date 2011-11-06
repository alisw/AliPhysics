#ifndef ALIANACHARGEDPARTICLES_H
#define ALIANACHARGEDPARTICLES_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: $ */

//_________________________________________________________________________
// Example class on how to read AODCaloClusters, ESDCaloCells and AODTracks and how 
// fill AODs with PWG4PartCorr analysis frame
// Select the type of detector information that you want to analyze, CTS (tracking), PHOS or EMCAL
// Select the PID custer type of the calorimeters
// Set min momentum of the cluster/tracks
// Fill few histograms
//
//-- Author: Gustavo Conesa (INFN-LNF)
// Root system
class TH2F; 

// Analysis system
#include "AliAnaPartCorrBaseClass.h"
 
class AliAnaChargedParticles : public AliAnaPartCorrBaseClass {
  
 public: 
  AliAnaChargedParticles() ; // default ctor
  virtual ~AliAnaChargedParticles() { ; } //virtual dtor
 private:  
  AliAnaChargedParticles(const AliAnaChargedParticles & g) ; // cpy ctor
  AliAnaChargedParticles & operator = (const AliAnaChargedParticles & g) ;//cpy assignment

 public:

  TList * GetCreateOutputObjects();
  
  void    Init();
  
  void    InitParameters();
  
  void    Print(const Option_t * opt) const;
  
  void    MakeAnalysisFillAOD()  ;
  
  void    MakeAnalysisFillHistograms() ; 
  
  Int_t   GetPdgOfSelectedCharged()  const  { return fPdg ; }
  void    SelectChargedWithPdg( Int_t pdg ) { fPdg = pdg  ; }
    
 private:
  
  Int_t  fPdg ; //identified particle id
  
  //Histograms 
  TH1F * fhNtracks;     //! track multiplicity distribution
  TH1F * fhPt;          //! pT distribution
  TH2F * fhPhiNeg;      //! phi distribution vs pT, negative
  TH2F * fhEtaNeg;      //! eta distribution vs pT, negative
  TH2F * fhPhiPos;      //! phi distribution vs pT, positive
  TH2F * fhEtaPos;      //! eta distribution vs pT, positive
  TH2F * fhEtaPhiPos;   //! eta vs phi distribution of positive charge  
  TH2F * fhEtaPhiNeg;   //! eta vs phi distribution of negative charge  

  //MC
  TH1F * fhPtPion;      //! pT distribution
  TH2F * fhPhiPion;     //! phi distribution vs pT
  TH2F * fhEtaPion;     //! eta distribution vs pT
  
  TH1F * fhPtProton;    //! pT distribution
  TH2F * fhPhiProton;   //! phi distribution vs pT
  TH2F * fhEtaProton;   //! eta distribution vs pT
  
  TH1F * fhPtElectron;  //! pT distribution
  TH2F * fhPhiElectron; //! phi distribution vs pT
  TH2F * fhEtaElectron; //! eta distribution vs pT
  
  TH1F * fhPtKaon;      //! pT distribution
  TH2F * fhPhiKaon;     //! phi distribution vs pT
  TH2F * fhEtaKaon;     //! eta distribution vs pT
  
  TH1F * fhPtUnknown;   //! pT distribution
  TH2F * fhPhiUnknown;  //! phi distribution vs pT
  TH2F * fhEtaUnknown;  //! eta distribution vs pT
  
  ClassDef(AliAnaChargedParticles,3)

} ;


#endif //ALIANACHARGEDPARTICLES_H



