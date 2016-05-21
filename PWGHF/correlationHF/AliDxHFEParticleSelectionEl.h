//-*- Mode: C++ -*-
// // $Id$

//* This file is property of and copyright by the ALICE Project        * 
//* ALICE Experiment at CERN, All rights reserved.                     *
//* See cxx source for full Copyright notice                           *

/// @file   AliDxHFEParticleSelectionEl.h
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-03-19
/// @brief  Electron selection for D-HFE correlations
///

#ifndef ALIDXHFEPARTICLESELECTIONEL_H
#define ALIDXHFEPARTICLESELECTIONEL_H

#include "AliDxHFEParticleSelection.h"

// 2012-09-17: there has been a problem in the dictionary generation for par file
// compilation, so we have to include the header files indicated below
//  Generating dictionary ...
//  In file included from $ROOTSYS/include/TObject.h:230:0,
//                   from G__PWGHFcorrelationHF.h:32,
//                   from G__PWGHFcorrelationHF.cxx:17:
//  $ROOTSYS/include/TBuffer.h: In function ‘TBuffer& operator>>(TBuffer&, Tmpl*&) [with Tmpl = AliHFEcuts]’:
//  G__PWGHFcorrelationHF.cxx:1658:15:   instantiated from here
//  $ROOTSYS/include/TBuffer.h:373:47: error: invalid use of incomplete type ‘struct AliHFEcuts’
//  correlationHF/AliDxHFEParticleSelectionEl.h:25:7: error: forward declaration of ‘struct AliHFEcuts’

#include "AliHFEcuts.h" // need to directly include to avoid compilation error in the dictionary
#include "AliHFEpid.h"  // need to directly include to avoid compilation error in the dictionary
#include "AliAODTrack.h"

class AliAODTrack;
class AliPID;
class AliPIDResponse;
class AliHFEvarManager;
class AliHFEpidBase;
class AliHFEtools;
class AliVEvent;
class AliSelectNonHFE;
class AliCFManager;
class TList;
class AliESDtrackCuts;

/**
 * @class AliDxHFEParticleSelectionEl
 * Electron selection for D-HFE correlations, implements the specific
 * selection criteria.
 */
class AliDxHFEParticleSelectionEl : public AliDxHFEParticleSelection {
  public:
  /// constructor
  AliDxHFEParticleSelectionEl(const char* opt="");
  /// destructor
  virtual ~AliDxHFEParticleSelectionEl();
  
  enum {
    kCutHFE = 0,
    kCutPIDTOFTPC = 1,
    kCutPIDTOF = 2,
    kCutPIDTPC = 3,
    kCutList=3,
    kNCuts
  };

  enum {
    kNotSelected=-2,
    kNoCuts=-1,
    kRecKineITSTPC=0,
    kRecPrim,
    kHFEcutsITS,
    kHFEcutsTOF,
    kHFEcutsTPC,
    kPIDTOF,
    kPIDTPC,
    kPIDTOFTPC,
    kRejPi,
    kRejProton,
    kINVMASS,
    kSelected,
    kNCutLabels
  };

  enum {
    kNoInvMass=0,
    kInvMassSingleSelected, // based on AliSelectNonHFE, one partcle passing track cuts
                            // and 2nd particle passing weaker cut
    kInvMassTwoSelected     // pairs of particles both passing the track cuts
  };


  ///overloaded from AliDxHFEParticleSelection: Init
  virtual int Init();

  /// overloaded from AliDxHFEParticleSelection: init the control objects
  virtual int InitControlObjects();
  virtual THnSparse* DefineTHnSparse();

  /// overloaded from AliDxHFEParticleSelection: check particle
  virtual TObjArray* Select(const AliVEvent* pEvent);
  using AliDxHFEParticleSelection::Select;

  virtual int IsSelected(AliVParticle* p, const AliVEvent* pEvent);

  virtual int HistogramParticleProperties(AliVParticle* p, int selected);

  virtual int FillParticleProperties(AliVParticle* p, Double_t* date, int dimension) const;

  /// set cuts object: a type cast check is implemented in the method
  virtual void SetCuts(TObject* /*cuts*/, int /*level*/=0);
  virtual void SetFinalCutStep(int cutstep){fFinalCutStep=cutstep;}
  virtual void SetStoreLastCutStep(bool storecut){fStoreCutStepInfo=storecut;}

  virtual void SetPIDResponse(const AliPIDResponse* const pidresp){fPIDResponse=(AliPIDResponse*)(pidresp);}
  
  Int_t GetLastSurvivedCutsStep() const{return fSurvivedCutStep;}

 protected:
  Bool_t CheckTOFPIDStatus(AliAODTrack *track) const;

 private:

  /// copy contructor prohibited
  AliDxHFEParticleSelectionEl(const AliDxHFEParticleSelectionEl&);
  /// assignment operator prohibited
  AliDxHFEParticleSelectionEl& operator=(const AliDxHFEParticleSelectionEl&);
 
  /// check cut of specified step, e.g.
  bool ProcessCutStep(Int_t cutStep, AliVParticle *track);
  int ParseArguments(const char* arguments);

  virtual void InvMassFilter(TList* elList, Bool_t* selIndx);

  AliHFEpid*    fPIDTOFTPC;          //! the PID object
  AliHFEpid*    fPIDTOF;             //! the PID TOF object
  AliHFEpid*    fPIDTPC;             //! the PID TOF object
  AliHFEpid*    fPIDTPCEMCAL;             //! the PID TOF object
  THnSparse*    fElectronProperties; // the particle properties of selected particles
  TList*        fHistoList;          // list of histograms
  TList*        fCutPidList;         // list for pid and cut objects
  AliPIDResponse* fPIDResponse;      // fPIDResponse
  AliHFEcuts*   fCuts;               //! Cuts for HF electrons
  AliSelectNonHFE* fSelNHFE;         // AliSelectNonHFE, inv-mass selection
  AliESDtrackCuts* fTrackCuts;       //! Track cuts for associated track for Inv Mass
  AliCFManager* fCFM;                //! Correction Framework Manager
  Int_t         fFinalCutStep;       // Holds the final cutstep
  Double_t      fInvMassLow;         // lower inv-mass cut
  Int_t         fUseInvMassCut;      // whether to use inv mass cut
  Int_t         fSystem;             // whether running on pp(0,default), PbPb(1) or pPb(2, not used yet)
  Int_t         fTrackNum;           // Current track number
  Double_t      fImpactParamCutRadial; // cut on Impact parameter (when storing all tracks)
  Double_t      fEtaCut;             // cut on eta (when storing all tracks)
  Int_t         fSurvivedCutStep;    // Last cutstep it survived
  Bool_t        fStoreCutStepInfo;   // Whether to store which cut step is the effective one
  Bool_t        fSetFilterBit;       // Whether or not to use filter bits 
  Int_t         fBit;                // Which bit to use
  Double_t      fMaxPtCombinedPID;   // Max pt for the combined PID of TPC and TOF (default = 999, thus combined for all pt)
  Double_t      fMaxPTOFWhenPresent; // Max p for using tof when present (without tofmatching requirements). Default -1
  Bool_t        fUseEMCAL;           // Whether or not to use EMCAL PID
  Double_t      fMinM20;             // Lower cut value for Shower shape, EMCAL, M20
  Double_t      fMaxM20;             // Upper cut value for Shower shape, EMCAL, M20
  Double_t      fMinM02;             // Lower cut value for Shower shape, EMCAL, M02
  Double_t      fMaxM02;             // Upper cut value for Shower shape, EMCAL, M02
  Double_t      fDispersion;         // Shower shape dispersion cut
  Double_t      fEovPMin;            // Lower cut value for E over P, EMCAL
  Double_t      fEovPMax;            // Upper cut value for E over P, EMCAL
  Bool_t        fUseTOFonlyWhenPresent; // use tof only when it's there
  Bool_t        fStopAfterFilterBit; // Stop selection after checking filterbits
  Bool_t        fCutLS;              // Cut also electrons from LS distribution (set to false by default)

  static const char* fgkCutBinNames[]; //! bin labels for cuts histogram
  
  ClassDef(AliDxHFEParticleSelectionEl, 7); 
};

#endif
