#ifndef ALIPHOSANALYZE_H
#define ALIPHOSANALYZE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// Algorythm class to analyze PHOSv1 events:
// Construct histograms and displays them.
// Use the macro EditorBar.C for best access to the functionnalities
//*--
//*-- Author : Yves Schutz (SUBATECH)

// --- ROOT system ---

class TFile ;
class TH1F ;
class TH2F ;

// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSv1 ;
class AliPHOSGeometry ;
class AliPHOSLoader ;
class AliRunLoader ; 

class AliPHOSAnalyze : public TObject {

public:

  AliPHOSAnalyze() ;              // ctor
  AliPHOSAnalyze(Text_t * name) ; // ctor
  AliPHOSAnalyze(const AliPHOSAnalyze & ana) ; // cpy ctor                   
  virtual ~AliPHOSAnalyze() ;     // dtor

  void DrawRecon(Int_t Nevent= 0,Int_t Nmod = 1) ; 
  // draws positions of entering of primaries and reconstructed objects in PHOS

  void InvariantMass() ;      // Photons invariant mass distributions

  void EnergyResolution () ;  // analyzes Energy resolution ;

  void PositionResolution() ; // analyzes Position resolution ;

  void Contamination() ;      // Counts contamination of photon spectrum

  void Ls() ; //Prints PHOS-related contents of TreeS, TreeD and TreeR

  void SetEnergyCorrection(Float_t ecor){fCorrection = ecor ;} 

  AliPHOSAnalyze & operator = (const AliPHOSAnalyze & /*rvalue*/)  {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented") ;
    return *this ; 
  }
 
private:

  Float_t CorrectedEnergy(Float_t ReconstEnergy)const
    {return ReconstEnergy * fCorrection;} 
  //Converts reconstructed energy (energy of the EMCRecPoint) to the energy of primary
  //The coeficient shoud be (and was) calculated usin Erec vs. Eprim plot 
  //(see Energy Resolution function). However, if one change parameters of reconstruction 
  //or geometry, one have to recalculate coefficient!

private:

  Float_t fCorrection ;               //! Conversion coefficient between True and Reconstructed energies
  Int_t   fEvt ;                      //! the evt number being processed 
  TString ffileName ;                 //! the root file that contains the data

  AliRunLoader* fRunLoader;           //! run loader of the specified filename

ClassDef(AliPHOSAnalyze,1)  // PHOSv1 event analyzis algorithm

};

#endif // AliPHOSANALYZE_H
