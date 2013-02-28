//-*- Mode: C++ -*-
// $Id$

//* This file is property of and copyright by the ALICE Project        * 
//* ALICE Experiment at CERN, All rights reserved.                     *
//* See cxx source for full Copyright notice                           *

/// @file   AliDxHFECorrelation.h
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-04-25
/// @brief  Worker class for D0-HF electron correlation
///

#ifndef ALIDXHFECORRELATION_H
#define ALIDXHFECORRELATION_H

#include "TNamed.h"

class TH1;
class THnSparse;
class TObject;
class TList;
class AliHFCorrelator;
class AliVParticle;
class TObjArray;
class AliVEvent;
class AliAnalysisCuts;

class AliDxHFECorrelation : public TNamed {
 public:
  /// default constructor
  AliDxHFECorrelation(const char* name=NULL);
  /// destructor
  virtual ~AliDxHFECorrelation();

  // event control histogram
  enum {
    kEventsAll = 0, // all events
    kEventsSel,     // selected events
    kEventsTriggered,     // events with D0s
    kEventsCorrelated,     // events with correlated D0s
    kNEventControlLabels
  };

  // Enums for setting trigger particle type
  enum {
    kD=0,       
    kElectron=1
  } ;

  // init
  int Init(const char* arguments="");

  // parse argument string
  int ParseArguments(const char* arguments);

  /// fill histograms from particles
  int Fill(const TObjArray* candidatesD0, const TObjArray* candidatesElectron, const AliVEvent* pEvent);

  /// histogram event properties
  virtual int HistogramEventProperties(int bin);
  virtual THnSparse* DefineTHnSparse();
  virtual int FillParticleProperties(AliVParticle* tr, AliVParticle* as, Double_t* data, int dimension) const;

  /// create control THnSparse
  THnSparse* CreateControlTHnSparse(const char* name,
				    int thnSize,
				    int* thnBins,
				    double* thnMin,
				    double* thnMax,
				    const char** binLabels) const;

  /// overloaded from TObject: cleanup
  virtual void Clear(Option_t * option ="");
  /// overloaded from TObject: print info
  virtual void Print(Option_t *option="") const;
  /// overloaded from TObject: draw histograms
  virtual void Draw(Option_t *option="");
  /// overloaded from TObject: find object by name
  virtual TObject* FindObject(const char *name) const;
  /// overloaded from TObject: find object by pointer
  virtual TObject* FindObject(const TObject *obj) const;
  /// overloaded from TObject: save to file
  virtual void     SaveAs(const char *filename="",Option_t *option="") const; // *MENU*

  virtual void SetCuts(AliAnalysisCuts* cuts) {fCuts=cuts;}
  virtual void SetUseMC(Bool_t useMC){fUseMC=useMC;}
  //void SetUseEventMixing(Bool_t useMixing) {fUseEventMixing=useMixing;}
  //void SetSystem(Bool_t system){fSystem=system;}
  //void SetPhiRange(Double_t min, Double_t max){fMinPhi=min; fMaxPhi=max;}
  // TODO: SetEventType only needed for MC. How to avoid this?
  virtual void SetEventType(int type){fEventType=type;}

  Bool_t GetUseMC() const {return fUseMC;}
  const TList* GetControlObjects() const {return fControlObjects;}
  Double_t GetMinPhi() const {return fMinPhi;}
  Double_t GetMaxPhi() const {return fMaxPhi;}
  Double_t GetDeltaPhi() const {return fDeltaPhi;}
  Double_t GetDeltaEta() const {return fDeltaEta;}
  inline int GetDimTHnSparse() const {return fDimThn;}
  Int_t GetTriggerParticleType() const {return fTriggerParticleType;}

  void EventMixingChecks(const AliVEvent* pEvent);

  AliDxHFECorrelation& operator+=(const AliDxHFECorrelation& other);

  // Probably not needed anymore, since code was changed to THnSparse
  // but keep here in case we need it later
  enum {
    khD0pT,         // TH1F
    khD0Phi,        // TH1F
    khD0Eta,        // TH1F
    khElectronpT,   // TH1F
    khElectronPhi,  // TH1F
    khElectronEta,  // TH1F
    kNofHistograms
  };

 protected:
  /// add control object to list, the base class becomes owner of the object
  int AddControlObject(TObject* pObj);

 /// set the dimension of THn and allocate filling array
  void InitTHnSparseArray(int dimension) {
    fDimThn=dimension; 
    if (fCorrArray) delete[] fCorrArray; fCorrArray=NULL;
    if (dimension>0) fCorrArray=new Double_t[dimension];
  }

  inline Double_t* ParticleProperties() const {return fCorrArray;}

 private:
  /// copy constructor
  AliDxHFECorrelation(const AliDxHFECorrelation& other);
  /// assignment operator
  AliDxHFECorrelation& operator=(const AliDxHFECorrelation& other);

  // 2012-09-18: when running on Grid the histograms were empty. We encountered
  // messages "cannot create object of class TH1" when writing the analysis manager
  // to file for Grid analysis.
  // This class had a TH1 member marked to be saved, the object though was part of
  // a list, also a member of the class. Root has a problem with the schema info
  // in that case.
  // Solved by marking fhEventControlCorr as transient, the cause, though, is not
  // understood

  TObjArray* fHistograms;        //  the histograms - for the moment not in use. 
  TList* fControlObjects;        //  list of control objects
  THnSparse* fCorrProperties;    //  the Correlation properties of selected particles
  TH1* fhEventControlCorr;       //! event control histogram (saved via control object list)
  AliAnalysisCuts *fCuts;        //! Cuts 
  Bool_t fUseMC;                 // use MC info
  AliHFCorrelator *fCorrelator;  //! object for correlations
  Bool_t fUseEventMixing;        // Run Event Mixing analysis
  Short_t fSystem;               // Which system pp/PbPb
  Double_t fMinPhi;              // Holds min phi
  Double_t fMaxPhi;              // Holds maxa phi
  Double_t fDeltaPhi;            // Delta Phi  
  Double_t fDeltaEta;            // Delta Eta
  int fDimThn;                   // Holds dim of THnSparse
  Double_t* fCorrArray;          //! filling array for THnSparse
  Int_t fEventType;              // Event type. Only needed for MC (fix)
  Int_t fTriggerParticleType;    // Which particle to trigger on

  static const char* fgkEventControlBinNames[];

  ClassDef(AliDxHFECorrelation, 5)
};
#endif
