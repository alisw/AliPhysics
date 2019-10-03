//-*- Mode: C++ -*-
// $Id$

//* This file is property of and copyright by the ALICE Project        * 
//* ALICE Experiment at CERN, All rights reserved.                     *
//* See cxx source for full Copyright notice                           *

/// @file   AliDxHFEParticleSelectionD0.h
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-03-19
/// @brief  D0 selection for D-HFE correlations
///

#ifndef ALIDXHFEPARTICLESELECTIOND0_H
#define ALIDXHFEPARTICLESELECTIOND0_H

#include "AliDxHFEParticleSelection.h"

class AliVEvent;
class AliRDHFCuts;
class TList;
class AliAnalysisCuts;

/**
 * @class AliDxHFEParticleSelectionD0
 * D0 selection for D-HFE correlations, implements the specific
 * selection criteria.
 */
class AliDxHFEParticleSelectionD0 : public AliDxHFEParticleSelection {
  public:
  /// constructor
  AliDxHFEParticleSelectionD0(const char* opt="");
  /// destructor
  virtual ~AliDxHFEParticleSelectionD0();

  enum {
    kCutD0 = 0,
    kCutList=1,
    kNCutObjects
  };

  enum {
    kDstar=0,
    kCandSelTrack,
    kIsInFinducialAcceptance,
    kNegPtbin,
    kNoDaugthers,
    kSelected0,
    kSelectedD0,
    kSelectedD0bar,
    kSelectedboth,
    kNCutLabels
  };

  enum {
    kD0Prompt=0,
    kD0FeedDown=1
  };


  /// overloaded from AliDxHFEParticleSelection: init the control objects
  virtual int InitControlObjects();
  virtual THnSparse* DefineTHnSparse();
  virtual int FillParticleProperties(AliVParticle* p, Double_t* date, int dimension) const;
  virtual AliVParticle* CreateParticle(AliVParticle* track);

  //Function for daughter control objects
  //TODO: move to AliDxHFEParticleSelection to be used for several particles?
  virtual int InitControlObjectsDaughters(TString name, int daughter);

  //Overloaded from AliDxHFEParticleSelection
  virtual TObjArray* Select(TObjArray* particles, AliVEvent* pEvent);
  using AliDxHFEParticleSelection::Select;

  // Get the D0 efficiency
  virtual double GetD0Eff(AliVParticle* tr);

  /// overloaded from AliDxHFEParticleSelection: check particle
  virtual int IsSelected(AliVParticle* p, const AliVEvent *pEvent=NULL);
  /// overloaded from AliDxHFEParticleSelection: set cuts
  virtual void SetCuts(TObject* /*cuts*/, int level=0);
  void SetInvMass(Double_t mass){fD0InvMass=mass;}
  void SetPtBin(Int_t ptbin){fPtBin=ptbin;}

  //AliRDHFCutsD0toKpi GetCuts()const {return fCuts;}
  Int_t  GetFillOnlyD0D0bar() const {return fFillOnlyD0D0bar;}
  Int_t GetPtBin() const {return fPtBin;}
  Double_t GetInvMass() const {return fD0InvMass;}
  AliRDHFCuts *GetHFCuts() const {return fCuts;}
  Double_t GetEventMult() const {return fMultEv;}
  AliAnalysisCuts *GetEffCutObject() const{ return fD0Eff;}

 protected:
  /// overloaded from AliDxHFEParticleSelection: histogram particle properties
  virtual int HistogramParticleProperties(AliVParticle* p, int selected=1);

 private:

  int ParseArguments(const char* arguments);

  /// copy contructor prohibited
  AliDxHFEParticleSelectionD0(const AliDxHFEParticleSelectionD0&);
  /// assignment operator prohibited
  AliDxHFEParticleSelectionD0& operator=(const AliDxHFEParticleSelectionD0&);

  THnSparse* fD0Properties;     //  the particle properties of selected particles
  THnSparse* fD0Daughter0;      //  the particle properties of selected particles
  THnSparse* fD0Daughter1;      //  the particle properties of selected particles
  AliRDHFCuts* fCuts;           //! pointer to external cuts object 
  AliAnalysisCuts *fD0Eff;      //! Object holding D0 eff maps 
  Int_t     fFillOnlyD0D0bar;   //  flag to set what to fill (0 = both, 1 = D0 only, 2 = D0bar only)
  Double_t  fD0InvMass;         // D0InvMass
  Int_t     fPtBin;             // Pt Bin
  TList*    fHistoList;         // list of histograms
  Bool_t    fUseD0Efficiency;   // Whether or not to correct for D0 efficiency
  Double_t  fMultEv;            // multiplicity of event
  Short_t fSystem;              // Which system pp/PbPb
  Short_t fUseCentrality;       // Using centrality or not

  static const char* fgkDgTrackControlBinNames[]; //! bin labels for track control histogram
  static const char* fgkCutBinNames[];            //! bin labels for cut histogram

  // TODO: at the moment the dimensions of the different THnSparse objects are different
  // needs to be consolidated
  // TODO: one might need particle properties of all and/or at different cut stages

  ClassDef(AliDxHFEParticleSelectionD0, 4);
};

#endif
