//-*- Mode: C++ -*-
// $Id$

//* This file is property of and copyright by the ALICE Project        * 
//* ALICE Experiment at CERN, All rights reserved.                     *
//* See cxx source for full Copyright notice                           *

/// @file   AliDxHFEParticleSelection.h
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-03-19
/// @brief  Base class for particle selection
///

#ifndef ALIDXHFEPARTICLESELECTION_H
#define ALIDXHFEPARTICLESELECTION_H

#include "TNamed.h"
#include "TString.h"
class AliVEvent;
class AliVParticle;
class TObjArray;
class TH1;
class THnSparse;

/**
 * @class AliDxHFEParticleSelection
 * This is the base class for particle selections for the D0 - HFE
 * correlation studies.
 *
 * Ideas:
 * - configurable particle selection
 * - eventually histogramming of particle properties before vs. after
 *   selection
 *
 * Might be that there is already something similar, then this class
 * can be merged with some other class.
 */
class AliDxHFEParticleSelection : public TNamed {
  public:
  /// constructor
  AliDxHFEParticleSelection(const char* name=NULL, const char* opt="");
  /// destructor
  virtual ~AliDxHFEParticleSelection();

  /// set options
  void SetOption(const char* opt) { fOption = opt; }
  /// overloaded from TObject: get option
  virtual Option_t* GetOption() const { return fOption;}

  /// init the control objects
  virtual int InitControlObjects();

  /// create selection, array contains only pointers but does not own the objects
  /// object array needs to be deleted by caller
  virtual TObjArray* Select(const AliVEvent* pEvent);

  /// create selection, array contains only pointers but does not own the objects
  /// object array needs to be deleted by caller
  virtual TObjArray* Select(TObjArray* pTracks);

  /// check and add track to internal array
  int CheckAndAdd(AliVParticle* p);
  /// get selected tracks
  const TObjArray* GetSelected() const {return fSelectedTracks;}

  /// check particle if it passes the selection criteria
  virtual bool IsSelected(AliVParticle* p);

  /// inherited from TObject: cleanup
  virtual void Clear(Option_t * option ="");
  /// inherited from TObject: print info
  virtual void Print(Option_t *option="") const;
  /// inherited from TObject: safe selection criteria
  virtual void SaveAs(const char *filename="",Option_t *option="") const;
  /// inherited from TObject: draw content
  virtual void Draw(Option_t *option="");
  /// inherited from TObject: find object by name
  virtual TObject* FindObject(const char *name) const;
  /// inherited from TObject: find object by pointer
  virtual TObject* FindObject(const TObject *obj) const;

 protected:
  /// add control object to list, the base class becomes owner of the object
  int AddControlObject(TObject* pObj);

  /// histogram particle properties
  virtual int HistogramParticleProperties(AliVParticle* p, bool selected=true);

 private:
  /// copy contructor prohibited
  AliDxHFEParticleSelection(const AliDxHFEParticleSelection&);
  /// assignment operator prohibited
  AliDxHFEParticleSelection& operator=(const AliDxHFEParticleSelection&);

  TString fOption; // option
  TObjArray* fSelectedTracks; //! array of selected tracks

  // control histograms, note: only the list is saved, pointers only used for fast access
  TList* fControlObjects; // list of control objects
  TH1* fhEventControl; //! event control histogram
  TH1* fhTrackControl; //! track control histogram

  ClassDef(AliDxHFEParticleSelection, 1);
};

#endif
