#ifndef ALICDFJETFINDER_H
#define ALICDFJETFINDER_H

/*
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.
 * See cxx source for full Copyright notice
 *
*/

/* $Id$ */

//  Definition of constants, structures and functions for jet finder

#include "AliJetFinder.h"

using namespace std ;


// structure of jet and particles container
struct varContainer
{
  Double_t pt;   // pt of particle/jet
  Double_t eta;  // eta of particle/jet
  Double_t phi;  // phi of particle/jet
  Int_t    njet; // njet is multiplicity of jet or if particles are stored , njet is index number of jet 
};

class AliCdfJetHeader;

class AliCdfJetFinder : public AliJetFinder
{
 public:
  AliCdfJetFinder();
  virtual        ~AliCdfJetFinder();

  void           CreateOutputObjects(TList * const histos);
  void           FindJets();
  void           InitData();
  void           FindCones();
  void           ComputeConesWeight();
  void           WriteJets() ;
  void           AnalizeJets();
  void           Clean();
    
 protected:
  AliCdfJetFinder ( const AliCdfJetFinder& jf );
  AliCdfJetFinder& operator = ( const AliCdfJetFinder& jf );

  TList*         fHistos;          //  List of histograms

  Bool_t         fAODwrite ;       //  write jets to AOD
  Bool_t         fAODtracksWrite;  //  write jet tracks to AOD
  Bool_t         fAnalyseJets;     //  analyse jets
	
  Int_t          fNJets;           //! counter of number of jets
  Int_t          fNPart;           //! number of particles in event
  Int_t          fNInC;            //! number of charged particles in event
  Int_t          fNInN;            //! number of neutral cells in event

  Double_t       fRadius;          // radius of jet 

  Int_t          fMinJetParticles; //  leading jet must contain AT LEAST fMinJetParticles
  Double_t       fJetPtCut;        //  leading jet must have AT LEAST fJetPtCut

  varContainer** fVectParticle;    //! container for Particles
  varContainer** fVectJet;         //! container for Jets

  Double_t*      fPtArray;         //! momentum array
  Int_t*         fIdxArray;        //! index array of sorted pts

  ClassDef(AliCdfJetFinder,3)      //  CDF jet finder

};
#endif

