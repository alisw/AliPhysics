#ifndef ALIANALYSISHE4_H
#define ALIANALYSISHE4_H

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// This analysis extracts pT-spectra of He4.                                //
// It is based on particles identifation via the dE/dx signal of the TPC.   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

class TH1;
class TH1F;
class TH2F;
class TH3F;
class TProfile2D;

class TList;
class TObjArray;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliHeader;
class AliESDpid;
class AliAnalysisUtils;

#include "THnSparse.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisHe4 : public AliAnalysisTaskSE {
 public:
  AliAnalysisHe4(const char *name);
  AliAnalysisHe4();
  virtual ~AliAnalysisHe4() {}
  //
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  //
  void           SetESDtrackCuts(AliESDtrackCuts * trackCuts){fESDtrackCuts = trackCuts;};
//   void           SetAlephParameters(const Double_t * parameters){for(Int_t j=0;j<5;j++) fAlephParameters[j] = parameters[j]; Initialize();};
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
  TH2F       * fHistMomCorr;           //! histogram for momentum and rapidity correction due to wrong propagation mass
  TH3F       * fHistEtaPtGen;          //! histogram for rapidity correction due to cuts on generation level
  TH2F       * fHistVertex;            //! histogram to monitor the vertex position and resolution
  TH1F       * fHistVertexRes;         //! histogram for difference between MC truth and rec vertex z position
  TH2F       * fHistVertexResTracks;   //! histogram to control vertex positon vs no. of primary tracks
  
  
  TH2F       * fHistDeDx;              //! histo for a dE/dx   
//   THnSparse  * fHistDeDx;              //! histo for a dE/dx   

  
  TH1F       * fHistHe3;              //! histo for He3 canidates
  TH1F       * fHistHe4;              //! histo for He4 canidates   
  TH1F       * fHistAntiHe3;          //! histo for anti He3 canidates
  TH1F       * fHistAntiHe4;          //! histo for anti He4 canidates   

  TProfile2D       * fHistTPCsigHe3;            //! histo to control tpc signal sigma cut for he3
  TProfile2D       * fHistTPCsigHe4;            //! histo to control tpc signal sigma cut for he3
  
  
  //
  AliAnalysisHe4(const AliAnalysisHe4&); 
  AliAnalysisHe4& operator=(const AliAnalysisHe4&); 

  ClassDef(AliAnalysisHe4, 1); 
};

#endif
