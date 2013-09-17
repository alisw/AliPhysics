#ifndef ALIANALYSISANTINUCLEI_H
#define ALIANALYSISANTINUCLEI_H

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


#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"

class AliAnalysisAntiNuclei : public AliAnalysisTaskSE {
 public:
  AliAnalysisAntiNuclei(const char *name);
  AliAnalysisAntiNuclei();
  virtual ~AliAnalysisAntiNuclei() {}
  //
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  //
  void           SetESDtrackCuts(AliESDtrackCuts * trackCuts){fESDtrackCuts = trackCuts;};
  void           SetAlephParameters(const Double_t * parameters){for(Int_t j=0;j<5;j++) fAlephParameters[j] = parameters[j]; Initialize();};
  void           SetIsMCtrue(Bool_t isMCdata = kTRUE){fMCtrue = isMCdata;};
  //
  void           Initialize();
  //
  
 private:
  //
  void  BinLogAxis(const TH1 *h);
  //
  AliESDEvent *fESD;                   //! ESD object
  TList       *fListHist;              //! list for histograms
  //
  AliESDtrackCuts * fESDtrackCuts;     // basic cut variables
  AliESDtrackCuts * fESDTrackCutsMult; // cuts for the MULTIPLICITY DETERMINATION
  AliESDpid       * fESDpid;           // basic TPC object for n-sigma cuts
  Bool_t        fMCtrue;               // flag if real data or MC is processed
  Double_t      fAlephParameters[5];   // Aleph Parameters for Bethe-Bloch
  //
  //
  //
  THnSparseF * fHistRealTracks;        //! histogram with all necessary information for real tracks
  THnSparseF * fHistMCparticles;       //! histogram with all necessary information for MC particles
  //
  TH3F       * fHistPidQA;             //! histogram for the QA of the PID
  TH2F       * fHistTofQA;             //! histogram for the QA of the PID
  TH2F       * fHistMult;              //! control histogram for multiplicity
  TH1F       * fHistCentrality;        //! control histogram for centrality
  TH2F       * fHistMomCorr;           //! histogram for momentum and rapidity correction due to wrong propagation mass
  TH3F       * fHistEtaPtGen;          //! histogram for rapidity correction due to cuts on generation level
  //
  AliAnalysisAntiNuclei(const AliAnalysisAntiNuclei&); 
  AliAnalysisAntiNuclei& operator=(const AliAnalysisAntiNuclei&); 

  ClassDef(AliAnalysisAntiNuclei, 1); 
};

#endif
