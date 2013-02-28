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
class AliPIDResponse;
class TObjArray;
class TH1;
class TH2;
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

  enum {
    kEventsAll = 0,
    kEventsSel,
    kEventsWithParticle,
    kNEventPropertyLabels
  };
  
  enum {
    kTrackAll = 0,
    kTrackSel,
    kNTrackPropertyLabels
  };

  /// set options
  void SetOption(const char* opt) { fOption = opt; }
  /// overloaded from TObject: get option
  virtual Option_t* GetOption() const { return fOption;}

  /// init the control objects
  virtual int Init();
  virtual int InitControlObjects();

  /// create selection from 'Tracks' member of the event,
  /// array contains only pointers but does not own the objects
  /// object array needs to be deleted by caller
  virtual TObjArray* Select(const AliVEvent* pEvent);
  /// create selection from the array of particles,
  /// array contains only pointers but does not own the objects
  /// object array needs to be deleted by caller
  virtual TObjArray* Select(TObjArray* particles, const AliVEvent* pEvent);

  virtual void SetPIDResponse(const AliPIDResponse* /*const pidresp*/){}

  // Get the list fControlObjects. 
  const TList* GetControlObjects() const {return fControlObjects;}

  /// histogram event properties
  virtual int HistogramEventProperties(int bin);

  virtual int FillParticleProperties(AliVParticle* p, Double_t* date, int dimension) const;
  virtual AliVParticle* CreateParticle(AliVParticle* track);

  /// check and add track to internal array
  int CheckAndAdd(AliVParticle* p);

  /// set cuts object: general TObject pointer is used as argument to support
  // different types; a type cast check is implemented in the method
  virtual void SetCuts(TObject* /*cuts*/, int /*level*/=0) {}

  // TODO: check whether that is needed, should be covered by the specific
  // child implementation
  Bool_t GetUseMC() const {return fUseMC;}

  /// get selected tracks
  const TObjArray* GetSelected() const {return fSelectedTracks;}

  /// check particle if it passes the selection criteria
  virtual int IsSelected(AliVParticle* p, const AliVEvent *pEvent=NULL);

  /// inherited from TObject: cleanup
  virtual void Clear(Option_t * option ="");
  /// inherited from TObject: print info
  virtual void Print(Option_t *option="") const;
  /// inherited from TObject: safe selection criteria
  virtual void SaveAs(const char *filename="", Option_t *option="") const;
  /// inherited from TObject: draw content
  virtual void Draw(Option_t *option="");
  /// inherited from TObject: find object by name
  virtual TObject* FindObject(const char *name) const;
  /// inherited from TObject: find object by pointer
  virtual TObject* FindObject(const TObject *obj) const;

  /// set verbosity
  void SetVerbosity(int verbosity) {fVerbosity=verbosity;}

  /// get verbosity
  inline int GetVerbosity() const {return fVerbosity;}

  /// get the dimension of THn, fixed
  inline int GetDimTHnSparse() const {return fDimThn;}

  /// create 2D control histogram
  TH2* CreateControl2DHistogram(const char* name,
				const char* title,
				double* nBins,
				const char* xaxis,
				const char* yaxis) const;

  /// create control histogram
  TH1* CreateControlHistogram(const char* name,
			      const char* title,
			      int nBins,
			      const char** binLabels) const;
  
  /// create control THnSparse
  THnSparse* CreateControlTHnSparse(const char* name,
				    int thnSize,
				    int* thnBins,
				    double* thnMin,
				    double* thnMax,
				    const char** binLabels) const;

  // define and create the THnSparse object
  // initializes also the dimension to be used further
  virtual THnSparse* DefineTHnSparse();

 protected:
  /// add control object to list, the base class becomes owner of the object
  int AddControlObject(TObject* pObj);

  /// histogram particle properties
  virtual int HistogramParticleProperties(AliVParticle* p, int selected=1);

  /// set the dimension of THn and allocate filling array
  void InitTHnSparseArray(int dimension) {
    fDimThn=dimension; 
    if (fParticleProperties) delete[] fParticleProperties; fParticleProperties=NULL;
    if (dimension>0) fParticleProperties=new Double_t[dimension];
  }

  inline Double_t* ParticleProperties() const {return fParticleProperties;}

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
  bool fUseMC;         // specific implementation for MC selection
  int fVerbosity;      //! verbosity
  int fDimThn;         //  dim of thnsparse
  Double_t* fParticleProperties;  //! filling array for THnSparse

  static const char* fgkEventControlBinNames[]; //! bin labels for event control histogram
  static const char* fgkTrackControlBinNames[]; //! bin labels for track control histogram

  ClassDef(AliDxHFEParticleSelection, 3);

};

#endif
