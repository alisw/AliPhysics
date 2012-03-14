#ifndef TOFSPECTRAPPANALYSIS_H
#define TOFSPECTRAPOANALYSIS_H

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// This analysis extracts pT-spectra of charged kaons, protons, and pions.  //
// It is based on particles identifation via the dE/dx signal of the TPC.   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

class TH1;
class TH1F;
class TH2F;
class TH3F;
class TList;
class TObjArray;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliHeader;
class AliESDpid;
#include "AliTOFT0v1.h"
#include "AliTOFT0maker.h"
#include "AliTOFcalib.h"
#include "AliCDBManager.h"
#include <TTree.h>

class AliAnalysisFilter;
class AliCFContainer;
class TDatabasePDG;

#include "AliAnalysisTask.h"
#include "AliESDVertex.h"
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"
#include "AliTOFT0v1.h"
#include "AliTOFT0maker.h"
#include "AliTOFcalib.h"
#include "AliCDBManager.h"



#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"

class TOFSpectrappAnalysis : public AliAnalysisTaskSE {
 public:
  TOFSpectrappAnalysis(const char *name);
  TOFSpectrappAnalysis();
  virtual ~TOFSpectrappAnalysis() {}
  //
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  //
  Bool_t         SelectOnImpPar(AliESDtrack* t);
  //
  void           SetESDtrackCuts(AliESDtrackCuts * trackCuts){fESDtrackCuts = trackCuts;};
  //void           SetAlephParameters(const Double_t * parameters){for(Int_t j=0;j<5;j++) fAlephParameters[j] = parameters[j]; Initialize();};
  Int_t           Mult();
  
  //
  
 private:
  //
  //void  BinLogAxis(const THnSparse *h, Int_t axisNumber);
  
  //
  AliESDEvent *fESD;                  //! ESD object
  AliESDtrackCuts * fESDtrackCuts;    // basic cut variables
  AliESDpid       * fESDpid;          // basic TPC object for n-sigma cuts
  Bool_t        fMCtrue;              // flag if real data or MC is processed
  Double_t      fAlephParameters[5];  // Aleph Parameters for Bethe-Bloch
  Float_t spdCorr;
  Int_t multiplicity;
  Double_t ZPrimVertex;
  Int_t frun; 
  Int_t frunOld;
  Bool_t fLoadOCDB;
  Bool_t correctTExp;  
  Bool_t calibrateESD;
  Bool_t useT0TOF;
  Double_t timeResolution; 
  Bool_t tuneTOFMC;
  TTree *fTreeTrack;
  //TTree *fTreeEv;
  TH1D *hNumEv;
  AliTOFcalib *tofCalib;
  AliTOFT0maker *t0maker;
  
  Float_t fDCAXY;
  Float_t fDCAZ;
  Int_t kselimppar;
  Int_t fmatch;
  Double_t fmom;
  Double_t flength;
  Int_t fsign;
  Double_t ftimetof;
  Double_t fexptimeel;
  Double_t fexptimemu;
  Double_t fexptimepi;
  Double_t fexptimeka;
  Double_t fexptimepr;
  Int_t ftofchan;
  Double_t feta;
  Double_t fphi;
  Double_t fmomtrasv;
  Float_t t0track;
  Float_t t0trackSigma;
  Double_t sigmael;
  Double_t sigmamu;
  Double_t sigmapi;
  Double_t sigmaka;
  Double_t sigmapr;
  Double_t TPCSignal;
  Float_t TPCSigmaPI;
  Float_t TPCSigmaKA;
  Float_t TPCSigmaPR;
  Int_t TOFlabel0;
  Int_t TOFlabel1;
  Int_t TOFlabel2;

    
  Double_t r1[5];
 
  
  //
  TOFSpectrappAnalysis(const TOFSpectrappAnalysis&); 
  TOFSpectrappAnalysis & operator=(const TOFSpectrappAnalysis&); 

  ClassDef(TOFSpectrappAnalysis, 1); 
};

#endif
