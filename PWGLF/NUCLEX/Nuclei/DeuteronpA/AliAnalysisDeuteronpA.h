#ifndef ALIANALYSISDEUTERONPA_H
#define ALIANALYSISDEUTERONPA_H

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
class AliAnalysisUtils;


#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"

class AliAnalysisDeuteronpA : public AliAnalysisTaskSE {
 public:
  AliAnalysisDeuteronpA(const char *name);
  AliAnalysisDeuteronpA();
  virtual ~AliAnalysisDeuteronpA() {}
  //
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  //
  void           SetESDtrackCuts(AliESDtrackCuts * trackCuts){fESDtrackCuts = trackCuts;};
  void           SetIsMCtrue(Bool_t isMCdata = kTRUE){fMCtrue = isMCdata;};
  void           SetRapCMSpA(Bool_t isRapCMSpA = kTRUE){fRapCMSpA = isRapCMSpA;};
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
  Bool_t        fRapCMSpA;             // flag if shift to CMS_NN system for pA
  AliAnalysisUtils  *fUtils;           // For vertex cut and pileup rejection
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
  TH2F       * fHistVertex;            //! histogram to monitor the vertex position and resolution
  TH1F       * fHistVertexRes;         //! histogram for difference between MC truth and rec vertex z position
  TH2F       * fHistVertexResTracks;   //! histogram to control vertex positon vs no. of primary tracks
  //
  AliAnalysisDeuteronpA(const AliAnalysisDeuteronpA&); 
  AliAnalysisDeuteronpA& operator=(const AliAnalysisDeuteronpA&); 

  ClassDef(AliAnalysisDeuteronpA, 2); 
};

#endif
