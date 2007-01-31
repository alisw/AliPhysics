#ifndef ALIGENGEVSIM_H
#define ALIGENGEVSIM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////////
//
// AliGenGeVSim is a class implementing GeVSim event generator.
// 
// GeVSim is a simple Monte-Carlo event generator for testing detector and 
// algorythm performance especialy concerning flow and event-by-event studies
//
// In this event generator particles are generated from thermal distributions 
// without any dynamics and addicional constrains. Distribution parameters like
// multiplicity, particle type yields, inverse slope parameters, flow coeficients 
// and expansion velocities are expleicite defined by the user.
//
// GeVSim contains four thermal distributions the same as
// MevSim event generator developed for STAR experiment.
//
// In addition custom distributions can be used be the mean 
// either two dimensional formula (TF2), a two dimensional histogram or
// two one dimensional histograms.
//  
// Azimuthal distribution is deconvoluted from (Pt,Y) distribution
// and is described by two Fourier coefficients representing 
// Directed and Elliptic flow. 
// 
////////////////////////////////////////////////////////////////////////////////
//
// To apply flow to event ganerated by an arbitraly event generator
// refer to AliGenAfterBurnerFlow class.
//
////////////////////////////////////////////////////////////////////////////////
//
// For examples, parameters and testing macros refer to:
// http:/home.cern.ch/radomski
// 
// for more detailed description refer to ALICE NOTE
// "GeVSim Monte-Carlo Event Generator"
// S.Radosmki, P. Foka.
//  
// Author:
// Sylwester Radomski,
// GSI, March 2002
//  
// S.Radomski@gsi.de
//
////////////////////////////////////////////////////////////////////////////////
//
// Updated and revised: September 2002, S. Radomski, GSI
//
////////////////////////////////////////////////////////////////////////////////


class TFormula;
class TF1;
class TF2;
class TH1D;
class TH2D;
class TObjArray;
class AliGeVSimParticle;

#include "AliGenerator.h"


class AliGenGeVSim : public AliGenerator {

 public:
  
  AliGenGeVSim();
  AliGenGeVSim(Float_t psi, Bool_t isMultTotal = kTRUE);
  
  virtual ~AliGenGeVSim();
  
  /////////////////////////////////////////////////////////////////
  
  void AddParticleType(AliGeVSimParticle *part);
  void SetMultTotal(Bool_t isTotal = kTRUE);
  
  void Init();
  void Generate();
  
  /////////////////////////////////////////////////////////////////
  
 private:

  AliGenGeVSim(const AliGenGeVSim & ggs);
  AliGenGeVSim & operator=(const AliGenGeVSim & /*ggs*/);


  Int_t   fModel;            // Selected model (1-7)
  Float_t fPsi;              // Reaction Plane angle (0-2pi)
  Bool_t  fIsMultTotal;      // Mode od multiplicity: total, dN/dY

  TF1 *fPtFormula;           //! Pt formula for model (1)
  TF1 *fYFormula;            //! Y formula for model (1)
  TF2 *fPtYFormula[4];       //! Pt,Y formulae for model (2)-(4)
  TF1 *fPhiFormula;          //! phi formula 
  
  TFormula *fCurrentForm;    //! currently used formula
  TH1D *fHist[2];            //! two 1D histograms (fModel == 6)
  TH2D *fPtYHist;            //! two-dimensional histogram (fModel == 7)

  TObjArray *fPartTypes;     // Registered particles
  
  void InitFormula();
  void SetFormula(Int_t pdg);
  void AdjustFormula();
  void DetermineReactionPlane();
  void GetRandomPtY(Double_t &pt, Double_t &y);

  Float_t GetdNdYToTotal();
  
  Bool_t CheckPtYPhi(Float_t pt, Float_t y, Float_t phi) const;  // for histograms only
  Bool_t CheckAcceptance(Float_t p[3]);
  
  Float_t FindScaler(Int_t paramId, Int_t pdg);
  
  /////////////////////////////////////////////////////////////////


  ClassDef(AliGenGeVSim, 2)

};

#endif
