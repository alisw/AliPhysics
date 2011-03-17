/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-----------------------------------------------------------------------
// Author : Marta Verweij - UU
//-----------------------------------------------------------------------

#ifndef ALIPWG4HIGHPTSPECTRA_H
#define ALIPWG4HIGHPTSPECTRA_H

#include "AliAnalysisTask.h"
#include "AliCFManager.h"

class TH1I;
class TH1F;
class TH1D;
class TProfile;
class TFile;
class TList;

//class AliCFManager;
class AliESDtrackCuts;
class AliESDEvent;
class AliESDVertex;
class AliMCEvent;
class AliStack;
class AliGenPythiaEventHeader;

class AliPWG4HighPtSpectra : public AliAnalysisTask {
 public:

  enum {
    kStepReconstructed          = 0,
    kStepSecondaries            = 1,
    kStepReconstructedMC        = 2,
    kStepMCAcceptance           = 3
  };

  AliPWG4HighPtSpectra();
  AliPWG4HighPtSpectra(const Char_t* name);
  // AliPWG4HighPtSpectra& operator= (const AliPWG4HighPtSpectra& c);
  //  AliPWG4HighPtSpectra(const AliPWG4HighPtSpectra& c);
  ~AliPWG4HighPtSpectra() {;};

  // ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  virtual void   LocalInit();
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual Bool_t Notify(); //Copied from AliAnalysisTaskJetSpectrum2

  Bool_t IsPbPb() {return fIsPbPb;}  //is PbPb data?
  Bool_t SelectEvent();    //decides if event is used for analysis
  Int_t CalculateCentrality(AliESDEvent *esd);

  //Setters
  void SetIsPbPb(Bool_t cs)                {fIsPbPb = cs;}
  void SetCentralityClass(int cent)        {fCentClass=cent;}
 
  // CORRECTION FRAMEWORK RELATED FUNCTIONS
  void     SetCFManagerPos(const AliCFManager* io1) {fCFManagerPos = io1;}   // global correction manager 
  const AliCFManager * GetCFManagerPos() const {return fCFManagerPos;}           // get corr manager 
  void     SetCFManagerNeg(const AliCFManager* io2) {fCFManagerNeg = io2;}   // global correction manager 
  const AliCFManager * GetCFManagerNeg() const {return fCFManagerNeg;}            // get corr manager

  //if fTrackType=0 (Global)
  //if fTrackType=1 (TPConly)
  //if fTrackType=2 (TPConly constrained)
  void SetTrackType(Int_t trackType) {fTrackType = trackType;}
  //AliESDtrackCuts setters
  void SetCuts(AliESDtrackCuts* trackCuts) {fTrackCuts = trackCuts;}

  void SetSigmaConstrainedMax(Double_t sigma) {fSigmaConstrainedMax=sigma;}

  // Data types
  Bool_t IsReadAODData()   const {return fReadAODData;}
  void   SetReadAODData(Bool_t flag=kTRUE) {fReadAODData=flag;}

  static AliGenPythiaEventHeader*  GetPythiaEventHeader(AliMCEvent *mcEvent);
  static Bool_t PythiaInfoFromFile(const char* currFile,Float_t &fXsec,Float_t &fTrials);// get the cross section and the trails either from pyxsec.root or from pysec_hists.root
  
 protected:
  Bool_t              fReadAODData ;       // flag for AOD/ESD input files
  const AliCFManager  *fCFManagerPos    ;  // pointer to the CF manager for positive charged particles
  const AliCFManager  *fCFManagerNeg    ;  // pointer to the CF manager for negative charged particles
 
  AliESDEvent *fESD;      //! ESD object
  AliMCEvent  *fMC;       //! MC event object
  AliStack    *fStack;    //! stack object

  const AliESDVertex   *fVtx;     //! vertex object

  Bool_t   fIsPbPb;               //  kTRUE if PbPb
  Int_t fCentClass;               // Select only events from predefined centrality class

  Int_t   fTrackType;     // Type of track to be used in analysis
  //AliESDtrackCuts options. Must be setted in AddTaskPWG4HighPTSpectra.C. They correspond with different steps in container.
  AliESDtrackCuts *fTrackCuts;           // trackCuts applied to global tracks

  Double_t fSigmaConstrainedMax;  // max sigma on constrained fit

 private:
  AliPWG4HighPtSpectra(const AliPWG4HighPtSpectra&);
  AliPWG4HighPtSpectra& operator=(const AliPWG4HighPtSpectra&);

  Float_t fAvgTrials;             // Average number of trials

  // Histograms
  //Number of events
  TList *fHistList;             //! List of output histograms
  TH1F  *fNEventAll;            //! Event counter
  TH1F  *fNEventSel;            //! Event counter: Selected events for analysis
  TH1F  *fNEventReject;         //! Book keeping of reason of rejecting events

  TH1F *fh1Centrality;                         //! Centrality

  TProfile*     fh1Xsec;                       //! pythia cross section and trials
  TH1F*         fh1Trials;                     //! trials which are added
  TH1F*         fh1PtHard;                     //! pt hard of the event
  TH1F*         fh1PtHardTrials;               //! pt hard of the event

  ClassDef(AliPWG4HighPtSpectra,3);
};

#endif
